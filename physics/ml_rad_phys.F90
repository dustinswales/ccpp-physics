! ###########################################################################################
!> \file ml_rad_phys.F90
!!
!! This module includes the ONNX based machine-learning emulator for longwave radiation detailed
!! in 10.22541/essoar.168319865.58439449/v1
!!
!! There are 24 predictors for this emulator:
!!         Name                       Units
!!    1  - 'pressure_pascals'                [Pa]
!!    2  - 'temperature_kelvins'             [K]
!!    3  - 'specific_humidity_kg_kg01'       [kg/kg]
!!    4  - 'relative_humidity_unitless'      [1]
!!    5  - 'liquid_water_content_kg_m03'     [kg/m3]
!!    6  - 'ice_water_content_kg_m03'        [kg/m3]
!!    7  - 'liquid_water_path_kg_m02'        [kg/m2]
!!    8  - 'ice_water_path_kg_m02'           [kg/m2]
!!    9  - 'vapour_path_kg_m02'              [kg/m2]
!!    10 - 'upward_liquid_water_path_kg_m02' [kg/m2]
!!    11 - 'upward_liquid_ice_path_kg_m02'   [kg/m2]
!!    12 - 'upward_vapour_path_kg_m02'       [kg/m2]
!!    13 - 'liquid_effective_radius_metres'  [m]
!!    14 - 'ice_effective_radius_metres'     [m]
!!    15 - 'o3_mixing_ratio_kg_kg01'         [kg/kg]
!!    16 - 'co2_concentration_ppmv'          [ppmv]
!!    17 - 'ch4_concentration_ppmv'          [ppmv]
!!    18 - 'n2o_concentration_ppmv'          [ppmv]
!!    19 - 'height_m_agl'                    [m]
!!    20 - 'height_thickness_metres'         [m]
!!    21 - 'pressure_thickness_pascals'      [Pa]
!!    22 - 'zenith_angle_radians'            [1]
!!    23 - 'surface_temperature_kelvins'     [K]
!!    24 - 'surface_emissivity'              [1]
!!
! ###########################################################################################
module ml_rad_phys
  use machine,  only: kind_phys, kind_dbl_prec
  use funcphys, only: fpvs
  use module_radiation_gases, only: NF_VGAS, getgases, getozn
  use module_rad_ml, only: ty_rad_ml_data, ty_rad_ml_ref_data, ip2io, is2D, predictor_names
  use mersenne_twister, only: random_setseed, random_number, random_stat
  use netcdf
  use inferof
  use iso_c_binding, only : c_double, c_int, c_float, c_char, c_null_char, c_ptr
  implicit none

  type(infero_model) :: model

  ! Default effective radii, used when coupling radiation to single-moment cloud microphysics
  ! This is set by the host using, effr_in=F
  real (kind_phys), parameter :: &
       reliq_def  = 10.0, & ! Default liquid radius: 10 microns
       reice_def  = 50.0    ! Default ice radius:  50 microns

  public ml_rad_phys_init, ml_rad_phys_run

contains
! ########################################################################################
!! \section arg_table_ml_rad_phys_init
!! \htmlinclude ml_rad_phys_init.html
!!
! ########################################################################################
  subroutine ml_rad_phys_init(do_ml_rad, infero_model_path, infero_model_type, rad_ml_data, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         do_ml_rad            ! Use ML emulator for LW radiation?
    character(len=128), intent(in) :: &
         infero_model_path, & ! 
         infero_model_type    !
    type(ty_rad_ml_data), intent(inout) :: &
         rad_ml_data          ! DDT containing training data (IN:raw,OUT:sorted)
    
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg               !
    integer, intent(out) :: &
         errflg               !

    ! Locals
    character(1024) :: yaml_config
    integer :: iCol, iPred, iLay, iName, count
    logical,dimension(rad_ml_data%nCol) :: isort
    real(kind_phys), dimension(rad_ml_data%npred_scalar, rad_ml_data%nCol) :: &
         scalar_predictor
    real(kind_phys), dimension(rad_ml_data%npred_vector, rad_ml_data%nLev, rad_ml_data%nCol) :: &
         vector_predictor

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ml_rad) return

    ! ######################################################################################
    !
    ! Sort training data. Needed for normalization process for predictor matrix.
    !
    ! ######################################################################################
    ! Scalars
    do iPred= 1,rad_ml_data%npred_scalar
       isort(:) = .true.
       do iCol = 1,rad_ml_data%nCol
          scalar_predictor(iPred,iCol) = minval(rad_ml_data%scalar_predictor(iPred,:),isort)
          isort(minloc(rad_ml_data%scalar_predictor(iPred,:),isort)) = .false.
       end do
       rad_ml_data%scalar_predictor(iPred,:) = scalar_predictor(iPred,:)
    enddo

    ! Two-dimensional arrays
    do iPred = 1,rad_ml_data%npred_vector
       do iLay = 1,rad_ml_data%nLev
          isort(:) = .true.
          do iCol = 1,rad_ml_data%nCol
             vector_predictor(iPred,iLay,iCol) = minval(rad_ml_data%vector_predictor(iPred,iLay,:),isort)
             isort(minloc(rad_ml_data%vector_predictor(iPred,iLay,:),isort)) = .false.
          end do
          rad_ml_data%vector_predictor(iPred,iLay,:) = vector_predictor(iPred,iLay,:)
       enddo
    enddo

    ! ######################################################################################
    !
    ! Determine the mapping between the fields in training data and the order of expected 
    ! inputs to the emulator. Store in "ip2io"
    ! Also, store which data container the predictors are in (1D vs 2D)
    ! (The emulator expects the predictors in order. The inputs may not be in this order)
    !
    ! ######################################################################################
    count = 0
    do iPred= 1,rad_ml_data%npred_scalar
       do iName = 1,size(predictor_names)
          print*,'     ',iName,predictor_names(iName),len(trim(predictor_names(iName)))
          if (trim(predictor_names(iName)) == trim(rad_ml_data%scalar_predictor_names(iPred))) then
             ip2io(iName) = iPred
             is2D(iName)  = .false.
             count = count + 1
          endif
       enddo
    enddo

    do iPred = 1,rad_ml_data%npred_vector
       do iName = 1,size(predictor_names)
          if (trim(predictor_names(iName)) == trim(rad_ml_data%vector_predictor_names(iPred))) then
             ip2io(iName) = iPred
             is2D(iName)  = .true.
             count = count + 1
          endif
       enddo
    enddo

    ! ######################################################################################
    !
    ! Initialize Infero
    !
    ! ######################################################################################
    call infero_check(infero_initialise())

    ! YAML config string
    yaml_config = "---"//NEW_LINE('A') &
         //"  path: "//TRIM(infero_model_path)//NEW_LINE('A') &
         //"  type: "//TRIM(infero_model_type)//c_null_char

    ! Get a infero model
    call infero_check(model%initialise_from_yaml_string(yaml_config))

  end subroutine ml_rad_phys_init

! #########################################################################################
!! \section arg_table_ml_rad_phys_run
!! \htmlinclude ml_rad_phys_run.html
!!
! #########################################################################################
  subroutine ml_rad_phys_run(do_ml_rad, effr_in, debug, nCol, nLev, i_cldliq, i_cldice,   &
       i_ozone, ico2, isubc, icseed, semis, lon, lat, prsl, tgrs, prslk, prsi, cld_reliq, cld_reice,     &
       qgrs, rad_ml_data, ref_data, con_epsqs, con_eps, con_epsm1, con_rd, con_fvirt,     &
       con_g, con_pi, errmsg, errflg, &
! JUST USED FOR INLINE COMPARISION WITH ACTIVE RADIATION SCHEME
       htrlw, sfcflw, sfcfsw, topflw, topfsw)
    use module_radsw_parameters,             only: topfsw_type, sfcfsw_type
    use module_radlw_parameters,             only: topflw_type, sfcflw_type

    ! Inputs
    type(ty_rad_ml_data), intent(in) :: &
         rad_ml_data    ! DDT containing training data
    type(ty_rad_ml_ref_data), intent(in) :: &
         ref_data       ! DDT containing reference data.
    logical, intent(in) :: &
         do_ml_rad,   & ! Use ML emulator for LW radiation?
         effr_in,     & ! Provide hydrometeor radii from macrophysics? 
         debug          ! Debug mode?
    integer, intent(in) ::  &
         nCol,        & ! Number of horizontal grid points
         nLev,        & ! Number of vertical layers
         i_cldliq,    & ! Index into tracer array for cloud liquid.
         i_cldice,    & ! Index into tracer array for cloud ice.
         i_ozone,     & ! Index into tracer array for ozone concentration.
         ico2,        & ! Flag for co2 radiation scheme
         isubc          ! Flag for cloud-seeding (rng) for cloud-sampling
    integer,intent(in),dimension(:) :: &
         icseed         ! Seed for random number generation for longwave radiation
    real(kind_phys), dimension(:), intent(in) :: &
         lon,         & ! Longitude
         lat,         & ! Latitude
         semis          ! Longwave surface emissivity
    real(kind_phys), dimension(:,:), intent(in) :: &
         prsl,        & ! Pressure at model-layer centers (Pa)
         tgrs,        & ! Temperature at model-layer centers (K)
         prslk,       & ! Exner function at model layer centers (1)
         prsi,        & ! Pressure at model-interfaces (Pa)
         cld_reliq,   & ! Effective radius (m)
         cld_reice      ! Effective radius (m)
    real(kind_phys), dimension(:,:,:), intent(in) :: &
         qgrs           ! Tracer concentrations (kg/kg)
    real(kind_phys), intent(in) :: &
         con_epsqs,   & ! Physical constant: Minimum saturation mixing-ratio (kg/kg)
         con_eps,     & ! Physical constant: Epsilon (Rd/Rv)
         con_epsm1,   & ! Physical constant: Epsilon (Rd/Rv) minus one
         con_rd,      & ! Physical constant: gas-constant for dry air
         con_fvirt,   & ! Physical constant: Inverse of epsilon minus one
         con_g,       & ! Physical constant: gravitational constant
         con_pi         ! Physical constant: Pi

    ! JUST USED FOR INLINE COMPARISION WITH ACTIVE RADIATION SCHEME
    real(kind_phys), dimension(:,:), intent(in) :: &
         htrlw
    type(sfcflw_type), dimension(:), intent(in) :: sfcflw
    type(sfcfsw_type), dimension(:), intent(in) :: sfcfsw
    type(topflw_type), dimension(:), intent(in) :: topflw
    type(topfsw_type), dimension(:), intent(in) :: topfsw

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg

    ! Locals
    logical :: top_at_1
    integer :: iSFC, iTOA, stride, ipred, ilev, icase, iinf, iLay, iCol
    integer, dimension(nCol) :: ipseed
    real(kind_phys) :: es, qs, dp, tem1, tem2, pfac, ranku(nCol), rankn(nCol), rankk(1)
    real(kind_phys), dimension(nLev) :: tv, rho, tempVar
    real(kind_phys), dimension(nLev+1) :: hgtb, zo, zi
    real(kind_phys), dimension(nCol, nLev) :: o3_lay
    real(kind_phys), dimension(nCol, nLev, NF_VGAS) :: gas_vmr
    real(c_float), allocatable :: it2f(:,:,:) ! data for inference in profile, height, value order
    real(c_float), allocatable :: ot2f(:,:)   ! data from inference in profile, height order
    real(c_float), dimension(nCol, nLev, size(predictor_names)) :: predictor_matrix, upredictor_matrix
    real(c_float), dimension(nCol, nLev+2) :: target_matrix

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Are we supposed to be here?
    if (.not. do_ml_rad) return

    ! What is vertical ordering of the host?
    top_at_1 = (prsi(1,1) .lt.  prsi(1, nLev))
    if (top_at_1) then
       iSFC   = nLev
       iTOA   = 1
       stride = -1
    else
       iSFC   = 1
       iTOA   = nLev
       stride = 1
    endif

    ! Change random number seed value for each radiation invocation (isubc =1 or 2).
    if(isubc == 1) then      ! advance prescribed permutation seed
       do iCol=1,nCol
          ipseed(iCol) = 128 + iCol - 1
       enddo
    elseif (isubc == 2) then ! use input array of permutaion seeds
       do iCol=1,nCol
          ipseed(iCol) = icseed(iCol)
       enddo
    endif

    ! ######################################################################################
    ! Set up trace gas concentrations
    ! ######################################################################################

    ! Do we need to get climatological ozone? (only if not using prognostic ozone)
    if (i_ozone .le. 0) then
       call getozn (prslk, lat, nCol, nLev, top_at_1, o3_lay)
    endif

    ! Get trace-gas concentrations.
    call getgases (prsi/100., lon, lat, nCol, nLev, ico2, top_at_1, con_pi, gas_vmr)

    ! ######################################################################################
    ! Compute fields for predictor matrix and populate,
    ! ######################################################################################
    tem1 = 1._kind_phys/con_g
    tem2 = con_rd/con_g
    do iCol=1,nCol
       do iLay=1,nLev
          ! Thermodynamic temporaries
          es        = min( prsl(iCol,iLay),  fpvs( tgrs(iCol,iLay) ) )
          qs        = max( con_epsqs, con_eps * es / (prsl(iCol,iLay) + con_epsm1*es) )
          tv(iLay)  = tgrs(iCol,iLay) * (1._kind_phys + con_fvirt*qgrs(iCol,iLay,1))
          rho(iLay) = prsl(iCol,iLay)/(con_rd*tv(iLay))
          dp        = abs(prsi(iCol,iLay+1)-prsi(iCol,iLay))

          ! Layer pressure (Pa)
          predictor_matrix(iCol,iLay,1)  = prsl(iCol,iLay)

          ! Layer temperature (K)
          predictor_matrix(iCol,iLay,2)  = tgrs(iCol,iLay)

          ! Layer specific-humidity (kg/kg)
          predictor_matrix(iCol,iLay,3)  = qgrs(iCol,iLay,1)

          ! Compute layer relative-humidity (kg/kg)->(1)
          predictor_matrix(iCol,iLay,4)  = max( 0._kind_phys, min( 1._kind_phys, max(con_epsqs, qgrs(iCol,iLay,1))/qs ) )

          ! Compute layer liquid/Ice water content (kg/kg)->(kg/m3).
          predictor_matrix(iCol,iLay,5)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldliq)*rho(iLay))
          predictor_matrix(iCol,iLay,6)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldice)*rho(iLay))

          ! Compute layer liquid/ice/vapor condensate path, from mixing ratios (kg/kg)->(kg/m2).
          predictor_matrix(iCol,iLay,7)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldliq) * tem1 * dp)
          predictor_matrix(iCol,iLay,8)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldice) * tem1 * dp)
          predictor_matrix(iCol,iLay,9)  = max(0._kind_phys, qgrs(iCol,iLay,1) * tem1 * dp)

          ! Cloud effective radii (m).
          if (effr_in) then
             predictor_matrix(iCol,iLay,13) = cld_reliq(iCol,iLay)
             predictor_matrix(iCol,iLay,14) = cld_reice(iCol,iLay)
          else
             predictor_matrix(iCol,iLay,13) = reliq_def*1e-6 ! (microns)->(meters)
             predictor_matrix(iCol,iLay,14) = reice_def*1e-6 ! (microns)->(meters)
          endif

          ! Ozone mixing-ratio (kg/kg).
          if (i_ozone .gt. 0) then
             predictor_matrix(iCol,iLay,15) = qgrs(iCol,iLay,i_ozone)
          else
             predictor_matrix(iCol,iLay,15) = o3_lay(iCol,iLay)
          endif

          ! Trace gases (ppmv)
          predictor_matrix(iCol,iLay,16) = gas_vmr(iCol,iLay,1)*1e6
          predictor_matrix(iCol,iLay,17) = gas_vmr(iCol,iLay,3)*1e6
          predictor_matrix(iCol,iLay,18) = gas_vmr(iCol,iLay,2)*1e6

          ! Pressure-thickness (Pa)
          predictor_matrix(iCol,iLay,21) = dp

       enddo ! END vertical loop
       
       ! Zenith-angle (1) Why do we need this for Longwave?
       predictor_matrix(iCol,:,22) = 0._kind_phys

       ! Surface temperature (K)
       predictor_matrix(iCol,:,23) = tgrs(iCol,iSFC)

       ! Surface emissivity (1)
       predictor_matrix(iCol,:,24) = semis(iCol)

       ! Compute vertically integrated (in both directions) liquid/ice/vapor condensate path (kg/m2).
       if (top_at_1) then
          do iLay=1,nLev
             predictor_matrix(iCol,iLay,10) = sum(predictor_matrix(iCol,iLay:nLev,7))
             predictor_matrix(iCol,iLay,11) = sum(predictor_matrix(iCol,iLay:nLev,8))
             predictor_matrix(iCol,iLay,12) = sum(predictor_matrix(iCol,iLay:nLev,9))
          enddo
          do iLay=1,nLev
             predictor_matrix(iCol,iLay,7)  = sum(predictor_matrix(iCol,1:iLay,7))
             predictor_matrix(iCol,iLay,8)  = sum(predictor_matrix(iCol,1:iLay,8))
             predictor_matrix(iCol,iLay,9)  = sum(predictor_matrix(iCol,1:iLay,9))
          enddo
       else
          do iLay=1,nLev
             predictor_matrix(iCol,iLay,10) = sum(predictor_matrix(iCol,1:iLay,7))
             predictor_matrix(iCol,iLay,11) = sum(predictor_matrix(iCol,1:iLay,8))
             predictor_matrix(iCol,iLay,12) = sum(predictor_matrix(iCol,1:iLay,9))
          enddo
          do iLay=1,nLev
             predictor_matrix(iCol,iLay,7)  = sum(predictor_matrix(iCol,iLay:nLev,7))
             predictor_matrix(iCol,iLay,8)  = sum(predictor_matrix(iCol,iLay:nLev,8))
             predictor_matrix(iCol,iLay,9)  = sum(predictor_matrix(iCol,iLay:nLev,9))
          enddo
       endif

       ! Layer thickness and height above ground (m)
       if (top_at_1) then
          ! Layer thickness (m)
          do iLay=1,nLev
             predictor_matrix(iCol,iLay,20) = tem2 * abs(log(prsi(iCol,iLay+1)) - log(prsi(iCol,iLay))) * tv(iLay)
          enddo
          ! Height at layer boundaries
          hgtb(nLev+1) = 0._kind_phys
          do iLay=nLev,1,-1
             hgtb(iLay)= hgtb(iLay+1) + predictor_matrix(iCol,iLay,20)
          enddo
          ! Height at layer centers
          do iLay = nLev, 1, -1
             pfac = abs(log(prsi(iCol,iLay+1)) - log(prsl(iCol,iLay))) /  &
                    abs(log(prsi(iCol,iLay+1)) - log(prsi(iCol,iLay)))
             predictor_matrix(iCol,iLay,19) = hgtb(iLay+1) + pfac * (hgtb(iLay) - hgtb(iLay+1))
          enddo
       else ! SFC->TOA
          ! Layer thickness (m)
          do iLay=nLev,1,-1
             predictor_matrix(iCol,iLay,20) = tem2 * abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay+1))) * tv(iLay)
          enddo
          ! Height at layer boundaries
          hgtb(1) = 0._kind_phys
          do iLay=1,nLev
             hgtb(iLay+1)= hgtb(iLay) + predictor_matrix(iCol,iLay,20)
          enddo
          ! Height at layer centers
          do iLay = 1, nLev
             pfac = abs(log(prsi(iCol,iLay)) - log(prsl(iCol,iLay)  )) /  &
                    abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay+1)))
             predictor_matrix(iCol,iLay,19) = hgtb(iLay) + pfac * (hgtb(iLay+1) - hgtb(iLay))
          enddo ! END Vertical loop
       endif    ! END TOA->SFC 
    enddo       ! END column loop

    ! ######################################################################################
    ! Regrid in the vertical
    ! ######################################################################################
    !
    do iPred = 1,size(predictor_names)
       print*,'NAME:',predictor_names(iPred)
       do iCol=1,nCol
          zo(1)        = 0._kind_phys
          zo(2:nLev+1) = rad_ml_data%vector_predictor_mean(12,:)
          if (is2D(iPred)) then
             call rad_ml_data%vert_int(real(predictor_matrix(iCol,:,iPred),kind=kind_phys), hgtb, tempVar, zo)
             do iLay=1,nLev
                write(*,'(i5,4f12.4)') iLay, hgtb(iLay), predictor_matrix(iCol,iLay,iPred), zo(iLay), tempVar(iLay)
             enddo
             predictor_matrix(iCol,:,iPred) = tempVar
          endif
       enddo
    enddo

    ! ######################################################################################
    ! Normalize predictor matrix
    ! ######################################################################################
    upredictor_matrix = predictor_matrix
    do iPred = 1,size(predictor_names)
       do iLay=1,nLev
          do iCol=1,nCol
             ! Find uniform rank of predictor value wrt training data.
             if (is2D(iPred)) then
                rankk = real(minloc(abs(predictor_matrix(iCol,iLay,iPred) - rad_ml_data%vector_predictor(ip2io(iPred),iLay,:))),kind=kind_phys)
             else
                rankk = real(minloc(abs(predictor_matrix(iCol,iLay,iPred) - rad_ml_data%scalar_predictor(ip2io(iPred),:))),kind=kind_phys)
             endif
             ranku(iCol) = rankk(1)/rad_ml_data%nCol
             predictor_matrix(iCol,iLay,iPred) = sqrt(2._kind_phys)*erfinv(2.*ranku(iCol)-1)
             !predictor_matrix(iCol,iLay,iPred) = boxmuller_transform(2.*ranku(iCol)-1,ipseed(iCol),con_pi)
          enddo ! END Column loop
       enddo    ! END Vertical loop
    enddo       ! END Predictor loop

    ! ######################################################################################
    ! Begin Inference...
    ! ######################################################################################
    ! Run inferences
    call infero_check(model%infer(predictor_matrix, target_matrix))

    ! DEBUGGING INFO
    do iCol=1,nCol
       print*,'heating profile (K/day): '
       do iLay=1,nLev
          write(*,'(a17,i3,3f12.2)')'    ',iLay,prsl(iCol,iLay),target_matrix(iCol,iLay),htrlw(iCol,iLay)*3600.*24.
       enddo
       write(*,'(a20,3f12.2)') 'surface(down):    ',target_matrix(iCol,nLev+1),sfcflw(iCol)%dnfx0,sfcflw(iCol)%dnfxc
       write(*,'(a20,3f12.2)') 'toa(up):          ',target_matrix(iCol,nLev+2),topflw(iCol)%upfx0,topflw(iCol)%upfxc
    enddo

    write(*, '(a5, 25a10)')  'Layer', 'p', 'T', 'q', 'rh', 'LWC',                       &
         'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                            &
         'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                            &
         'dp', 'sza', 'Tsfc', 'sfc_emiss'
    write(*, '(a5, 25a10)')  '', '(Pa)', '(k)', '(g/kg)', '(1)', '(mg/m3)',             &
         '(mg/m3)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)',         &
         '(liq)','(ice)', '(mg/kg)','(ppmv)', '(ppbv)', '(ppmv)', '(m)', '(m)',         &
         '(Pa)', '(1)', '(K)', '(1)'
    do iCol=1,nCol
       do iLay=1,nLev
          write (*,'(i5,25f10.2)') iLay,             &
               upredictor_matrix(iCol,iLay,1),       &
               upredictor_matrix(iCol,iLay,2),       &
               1.e3*upredictor_matrix(iCol,iLay,3),  &
               upredictor_matrix(iCol,iLay,4),       &
               1.e6*upredictor_matrix(iCol,iLay,5),  &
               1.e6*upredictor_matrix(iCol,iLay,6),  &
               upredictor_matrix(iCol,iLay,7),       &
               upredictor_matrix(iCol,iLay,8),       &
               upredictor_matrix(iCol,iLay,9),       &
               upredictor_matrix(iCol,iLay,10),      &
               upredictor_matrix(iCol,iLay,11),      &
               upredictor_matrix(iCol,iLay,12),      &
               1.e6*upredictor_matrix(iCol,iLay,13), &
               1.e6*upredictor_matrix(iCol,iLay,14), &
               1.e6*upredictor_matrix(iCol,iLay,15), &
               upredictor_matrix(iCol,iLay,16),      &
               1.e3*upredictor_matrix(iCol,iLay,17), &
               upredictor_matrix(iCol,iLay,18),      &
               upredictor_matrix(iCol,iLay,19),      &
               upredictor_matrix(iCol,iLay,20),      &
               upredictor_matrix(iCol,iLay,21),      &
               upredictor_matrix(iCol,iLay,22),      &
               upredictor_matrix(iCol,iLay,23),      &
               upredictor_matrix(iCol,iLay,24)
          write (*,'(a5,25f10.2)')     '',&
               ref_data%unnorm_predictor_matrix(1,iLay,iCol),       &
               ref_data%unnorm_predictor_matrix(2,iLay,iCol),       &
               1.e3*ref_data%unnorm_predictor_matrix(3,iLay,iCol),  &
               ref_data%unnorm_predictor_matrix(4,iLay,iCol),       &
               1.e6*ref_data%unnorm_predictor_matrix(5,iLay,iCol),  &
               1.e6*ref_data%unnorm_predictor_matrix(6,iLay,iCol),  &
               ref_data%unnorm_predictor_matrix(7,iLay,iCol),       &
               ref_data%unnorm_predictor_matrix(8,iLay,iCol),       &
               ref_data%unnorm_predictor_matrix(9,iLay,iCol),       &
               ref_data%unnorm_predictor_matrix(10,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(11,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(12,iLay,iCol),      &
               1.e6*ref_data%unnorm_predictor_matrix(13,iLay,iCol), &
               1.e6*ref_data%unnorm_predictor_matrix(14,iLay,iCol), &
               1.e6*ref_data%unnorm_predictor_matrix(15,iLay,iCol), &
               ref_data%unnorm_predictor_matrix(16,iLay,iCol),      &
               1.e3*ref_data%unnorm_predictor_matrix(17,iLay,iCol), &
               ref_data%unnorm_predictor_matrix(18,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(19,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(20,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(21,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(22,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(23,iLay,iCol),      &
               ref_data%unnorm_predictor_matrix(24,iLay,iCol)
       enddo
    enddo
    write (*,'(a50)') '################################'
    write(*, '(a5, 25a10)')  'Layer', 'p', 'T', 'q', 'rh', 'LWC',                       &
         'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                            &
         'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                            &
         'dp', 'sza', 'Tsfc', 'sfc_emiss'
    do iCol=1,nCol
       do iLay=1,nLev
          write (*,'(i5,25f10.2)') iLay,                &
               predictor_matrix(iCol,iLay,1),           &
               predictor_matrix(iCol,iLay,2),           &
               predictor_matrix(iCol,iLay,3),           &
               predictor_matrix(iCol,iLay,4),           &
               predictor_matrix(iCol,iLay,5),           &
               predictor_matrix(iCol,iLay,6),           &
               predictor_matrix(iCol,iLay,7),           &
               predictor_matrix(iCol,iLay,8),           &
               predictor_matrix(iCol,iLay,9),           &
               predictor_matrix(iCol,iLay,10),          &
               predictor_matrix(iCol,iLay,11),          &
               predictor_matrix(iCol,iLay,12),          &
               predictor_matrix(iCol,iLay,13),          &
               predictor_matrix(iCol,iLay,14),          &
               predictor_matrix(iCol,iLay,15),          &
               predictor_matrix(iCol,iLay,16),          &
               predictor_matrix(iCol,iLay,17),          &
               predictor_matrix(iCol,iLay,18),          &
               predictor_matrix(iCol,iLay,19),          &
               predictor_matrix(iCol,iLay,20),          &
               predictor_matrix(iCol,iLay,21),          &
               predictor_matrix(iCol,iLay,22),          &
               predictor_matrix(iCol,iLay,23),          &
               predictor_matrix(iCol,iLay,24)
          write (*,'(i5,25f10.2)') iLay,                &
               ref_data%predictor_matrix(1,iLay,iCol+9),  &
               ref_data%predictor_matrix(2,iLay,iCol+9),  &
               ref_data%predictor_matrix(3,iLay,iCol+9),  &
               ref_data%predictor_matrix(4,iLay,iCol+9),  &
               ref_data%predictor_matrix(5,iLay,iCol+9),  &
               ref_data%predictor_matrix(6,iLay,iCol+9),  &
               ref_data%predictor_matrix(7,iLay,iCol+9),  &
               ref_data%predictor_matrix(8,iLay,iCol+9),  &
               ref_data%predictor_matrix(9,iLay,iCol+9),  &
               ref_data%predictor_matrix(10,iLay,iCol+9), &
               ref_data%predictor_matrix(11,iLay,iCol+9), &
               ref_data%predictor_matrix(12,iLay,iCol+9), &
               ref_data%predictor_matrix(13,iLay,iCol+9), &
               ref_data%predictor_matrix(14,iLay,iCol+9), &
               ref_data%predictor_matrix(15,iLay,iCol+9), &
               ref_data%predictor_matrix(16,iLay,iCol+9), &
               ref_data%predictor_matrix(17,iLay,iCol+9), &
               ref_data%predictor_matrix(18,iLay,iCol+9), &
               ref_data%predictor_matrix(19,iLay,iCol+9), &
               ref_data%predictor_matrix(20,iLay,iCol+9), &
               ref_data%predictor_matrix(21,iLay,iCol+9), &
               ref_data%predictor_matrix(22,iLay,iCol+9), &
               ref_data%predictor_matrix(23,iLay,iCol+9), &
               ref_data%predictor_matrix(24,iLay,iCol+9)
          print*,'-'
      enddo
    enddo
    !stop
    ! ######################################################################################
    ! Debug mode (Run using reference data).
    ! ######################################################################################
    if (debug) then
!       write(*, '(a5, 25a10)')  'Layer', 'p', 'T', 'q', 'rh', 'LWC',                       &
!            'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                            &
!            'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                            &
!            'dp', 'sza', 'Tsfc', 'sfc_emiss'
!       write(*, '(a5, 25a10)')  '', '(Pa)', '(k)', '(g/kg)', '(1)', '(mg/m3)',             &
!            '(mg/m3)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)',         &
!            '(liq)','(ice)', '(mg/kg)','(ppmv)', '(ppbv)', '(ppmv)', '(m)', '(m)',         &
!            '(Pa)', '(1)', '(K)', '(1)'
!       do iCol=1,nCol
!          do iLay=1,nLev
!             write (*,'(i5,25f10.2)') iLay, &
!                  ref_data%unnorm_predictor_matrix(1,iLay,iCol),       &
!                  ref_data%unnorm_predictor_matrix(2,iLay,iCol),       &
!                  1.e3*ref_data%unnorm_predictor_matrix(3,iLay,iCol),  &
!                  ref_data%unnorm_predictor_matrix(4,iLay,iCol),       &
!                  1.e6*ref_data%unnorm_predictor_matrix(5,iLay,iCol),  &
!                  1.e6*ref_data%unnorm_predictor_matrix(6,iLay,iCol),  &
!                  ref_data%unnorm_predictor_matrix(7,iLay,iCol),       &
!                  ref_data%unnorm_predictor_matrix(8,iLay,iCol),       &
!                  ref_data%unnorm_predictor_matrix(9,iLay,iCol),       &
!                  ref_data%unnorm_predictor_matrix(10,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(11,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(12,iLay,iCol),      &
!                  1.e6*ref_data%unnorm_predictor_matrix(13,iLay,iCol), &
!                  1.e6*ref_data%unnorm_predictor_matrix(14,iLay,iCol), &
!                  1.e6*ref_data%unnorm_predictor_matrix(15,iLay,iCol), &
!                  ref_data%unnorm_predictor_matrix(16,iLay,iCol),      &
!                  1.e3*ref_data%unnorm_predictor_matrix(17,iLay,iCol), &
!                  ref_data%unnorm_predictor_matrix(18,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(19,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(20,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(21,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(22,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(23,iLay,iCol),      &
!                  ref_data%unnorm_predictor_matrix(24,iLay,iCol)
!          enddo
!       enddo
!
       ! Bundle this into subroutine.
!       allocate(it2f(ref_data%nCol, ref_data%nLev, ref_data%nPred))
!       do ipred=1,ref_data%nPred
!          do ilev=1,ref_data%nLev
!             do icase=1,ref_data%nCol
!                it2f(icase, ilev, ipred) = ref_data%predictor_matrix(ipred, ilev, icase)
!             enddo
!          enddo
!       enddo
       
!       allocate(ot2f(ref_data%nCol, ref_data%nLev+ 2))
!       call infero_check(model%infer(it2f, ot2f ))
       
!       do icase = 1, ref_data%nCol
!          do ilev = 1, ref_data%nLev + 2
!             do ipred = 1, 1
!                if (ilev <= ref_data%nLev) then
!                   if (abs(ot2f(icase,ilev) - ref_data%vector_prediction(ipred,ilev,icase)) .gt. 1e-3) then
!                      write(*,*) "ERROR: output element ",icase,ilev, " (", ot2f(icase,ilev) ,") ", &
!                           "is different from expected value ", ref_data%vector_prediction(ipred,ilev,icase)
!                      stop 1
!                   end if
!                   print *, "( ",icase , ", ", ilev, ", ", ipred,") = ", ot2f(icase, ilev), ref_data%vector_prediction(ipred, ilev, icase)
!                else
!                   if (abs(ot2f(icase,ilev) - ref_data%scalar_prediction(ilev-127,icase)) .gt. 1e-3) then
!                      write(*,*) "ERROR: output element ",icase,ilev, " (", ot2f(icase,ilev) ,") ", &
!                           "is different from expected value ", ref_data%scalar_prediction(ilev-127,icase)
!                      stop 1
!                   end if
!                   print *,  "( ",icase , ", ", ilev, ", ", ipred,") = ", ot2f(icase, ilev), ref_data%scalar_prediction(ilev-127, icase)
!                end if
!             end do
!          end do
!       end do

    endif
    !stop
  end subroutine ml_rad_phys_run

! #########################################################################################
!! \section arg_table_ml_rad_phys_finalize
!! \htmlinclude ml_rad_phys_finalize.html
!!
! #########################################################################################
  subroutine ml_rad_phys_finalize(do_ml_rad, errmsg, errflg)
    ! Inputs
    logical,           intent(in) :: do_ml_rad
    ! Outputs
    character(len=*),  intent(out) :: errmsg
    integer,           intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ml_rad) return

    ! Free the model
    call infero_check(model%free())

    ! Finalize
    call infero_check(infero_finalise())

  end subroutine ml_rad_phys_finalize

  ! #########################################################################################
  ! 
  ! #########################################################################################
  function boxmuller_transform(var_u, seed, con_pi) result(var_n)
    ! IN
    real(kind_phys),intent(in) :: con_pi
    real(kind_phys),intent(in) :: var_u
    integer, intent(in) :: seed
    ! OUT
    real(kind_phys) :: var_n
    ! Locals
    real(kind_phys) :: rho, theta
    type(random_stat) :: rng_stat
    real(kind_dbl_prec) :: rng(1)
    
    call random_setseed(seed,rng_stat)
    call random_number(rng,rng_stat)
    rho   = sqrt(-2._kind_phys*log(rng(1)))
    theta = 2._kind_phys*con_pi*var_u
    var_n = rho*cos(theta)
    write(*,'(a12,4f12.4)') 'box-muller: ',rng,var_u,var_n

  end function boxmuller_transform

  function erfinv(var_in) result(var_out)
    real(kind_phys),intent(in) :: var_in
    real(kind_phys) :: var_out
    real(kind_dbl_prec),parameter,dimension(4) :: &
         a = (/ 0.886226899, -1.645349621,  0.914624893, -0.140543331/), &
         b = (/-2.118377725,  1.442710462, -0.329097515,  0.012229801/),&
         c = (/-1.970840454, -1.624906493,  3.429567803,  1.641345311/)
    real(kind_dbl_prec),parameter,dimension(2) :: &
         d = (/ 3.543889200,  1.637067800/)
    real(kind_phys) :: x, z

    if (abs(var_in) == 1._kind_phys) then
       x = var_in
    elseif(var_in < -0.7) then
       z = sqrt(-log((1._kind_phys+var_in)/2._kind_phys))
       x = -(((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0)
    else
       if (var_in > 0.7) then
          z = var_in*var_in
          x = var_in*(((a(4)*z+a(3))*z+a(2))*z+a(1))/((((b(4)*z+b(4))*z+b(2))*z+b(1))*z+1.0)
       else
          z = sqrt(-Log((1._kind_phys-var_in)/2._kind_phys));
          x = (((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0)
       endif
    endif
    var_out = x - (erf(x) - var_in) / (2.0/sqrt(3.14) * exp(-x*x))

  end function erfinv

end module ml_rad_phys
