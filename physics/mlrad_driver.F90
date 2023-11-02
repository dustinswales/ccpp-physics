! ###########################################################################################
!> \file mlrad_driver.F90
!!
!! This module includes the ONNX based machine-learning emulator for longwave and shortwave
!! radiation detailed in 10.22541/essoar.168319865.58439449/v1
!!
!! There are 24 predictors for the Longwave (LW) emulator:
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
!! There are 26 predictors for the Shortwave (SW) emulator:
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
!!    19 - 'aerosol_extinction_metres01'     [1]
!!    20 - 'height_m_agl'                    [m]
!!    21 - 'height_thickness_metres'         [m]
!!    22 - 'pressure_thickness_pascals'      [Pa]
!!    22 - 'zenith_angle_radians'            [1]
!!    23 - 'surface_temperature_kelvins'     [K]
!!    24 - 'albedo'                          [1]
!!    25 - 'aerosol_single_scattering_albedo'[1]
!!    26 - 'aerosol_asymmetry_param'         [1]
!!
! ###########################################################################################
module mlrad_driver
  use machine,  only: kind_phys, kind_dbl_prec
  use funcphys, only: fpvs
  use module_radiation_gases, only: NF_VGAS, getgases, getozn
  use module_mlrad, only: ty_mlrad_data, ip2io_lw, is2D_lw, pnames_lw, pnames_sw, ip2io_sw, is2D_sw
  use mersenne_twister, only: random_setseed, random_number, random_stat
  use inferof
  use iso_c_binding, only : c_double, c_int, c_float, c_char, c_null_char, c_ptr
  implicit none

  type(infero_model) :: model

  ! Default effective radii, used when coupling radiation to single-moment cloud microphysics
  ! This is set by the host using, effr_in=F
  real (kind_phys), parameter :: &
       reliq_def  = 10.0, & ! Default liquid radius: 10 microns
       reice_def  = 50.0    ! Default ice radius:  50 microns

  public mlrad_driver_init, mlrad_driver_run

contains
! ########################################################################################
!! \section arg_table_mlrad_driver_init
!! \htmlinclude mlrad_driver_init.html
!!
! ########################################################################################
  subroutine mlrad_driver_init(do_mlrad, infero_mpath, infero_mtype, mlrad_data,         &
       errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         do_mlrad        ! Use ML emulator for radiation?
    character(len=128), intent(in) :: &
         infero_mpath, & ! Infero model path
         infero_mtype    ! Infero model type
    type(ty_mlrad_data), intent(inout) :: &
         mlrad_data      ! DDT containing training data (IN:raw,OUT:sorted)
    
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg          ! CCPP error message
    integer, intent(out) :: &
         errflg          ! CCPP error flag

    ! Locals
    character(1024) :: yaml_config
    integer :: iCol, iPred, iLay, iName, count

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_mlrad) return

    ! ######################################################################################
    !
    ! Determine the mapping between the fields in training data and the order of expected 
    ! inputs to the emulator. Store in "ip2io"
    ! Also, store which data container the predictors are in (1D vs 2D)
    ! (The emulator expects the predictors in order. The inputs may not be in this order)
    !
    ! ######################################################################################
    ! Longwave
    count = 0
    do iPred= 1,mlrad_data%lw%npreds
       do iName = 1,size(pnames_lw)
          if (trim(pnames_lw(iName)) == trim(mlrad_data%lw%scalar_pnames(iPred))) then
             ip2io_lw(iName) = iPred
             is2D_lw(iName)  = .false.
             count = count + 1
          endif
       enddo
    enddo
    !
    do iPred = 1,mlrad_data%lw%npredv
       do iName = 1,size(pnames_lw)
          if (trim(pnames_lw(iName)) == trim(mlrad_data%lw%vector_pnames(iPred))) then
             ip2io_lw(iName) = iPred
             is2D_lw(iName)  = .true.
             count = count + 1
          endif
       enddo
    enddo

    ! Shortwave
    do iPred= 1,mlrad_data%sw%npreds
       do iName = 1,size(pnames_sw)
          if (trim(pnames_sw(iName)) == trim(mlrad_data%sw%scalar_pnames(iPred))) then
             ip2io_sw(iName) = iPred
             is2D_sw(iName)  = .false.
             count = count + 1
          endif
       enddo
    enddo
    !
    do iPred = 1,mlrad_data%sw%npredv
       do iName = 1,size(pnames_sw)
          if (trim(pnames_sw(iName)) == trim(mlrad_data%sw%vector_pnames(iPred))) then
             ip2io_sw(iName) = iPred
             is2D_sw(iName)  = .true.
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
         //"  path: "//TRIM(infero_mpath)//NEW_LINE('A') &
         //"  type: "//TRIM(infero_mtype)//c_null_char

    ! Get a infero model
    call infero_check(model%initialise_from_yaml_string(yaml_config))

  end subroutine mlrad_driver_init

! #########################################################################################
!! \section arg_table_mlrad_driver_run
!! \htmlinclude mlrad_driver_run.html
!!
! #########################################################################################
  subroutine mlrad_driver_run(do_mlrad, effr_in, debug, nCol, nLev, nDay, i_cldliq, i_cldice,  &
       i_ozone, ico2, isubc, icseed, idx, semis, lon, lat, prsl, tgrs, prslk, prsi, cld_reliq, &
       cld_reice, qgrs, mlrad_data, con_epsqs, con_eps, con_epsm1, con_rd,   &
       con_fvirt, con_g, con_pi, htrlw, sfcflw, sfcfsw, topflw, topfsw, errmsg, errflg)
    use module_radsw_parameters, only: topfsw_type, sfcfsw_type
    use module_radlw_parameters, only: topflw_type, sfcflw_type

    ! Inputs
    type(ty_mlrad_data), intent(in) :: &
         mlrad_data  ! DDT containing training data
    logical, intent(in) :: &
         do_mlrad,   & ! Use ML emulator for LW radiation?
         effr_in,     & ! Provide hydrometeor radii from macrophysics? 
         debug          ! Debug mode?
    integer, intent(in) ::  &
         nCol,        & ! Number of horizontal grid points
         nLev,        & ! Number of vertical layers
         nDay,        & ! Number of daylit columns
         i_cldliq,    & ! Index into tracer array for cloud liquid.
         i_cldice,    & ! Index into tracer array for cloud ice.
         i_ozone,     & ! Index into tracer array for ozone concentration.
         ico2,        & ! Flag for co2 radiation scheme
         isubc          ! Flag for cloud-seeding (rng) for cloud-sampling
    integer,intent(in),dimension(:) :: &
         icseed,      &  ! Seed for random number generation for longwave radiation
         idx             ! Index array for daytime points
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

    ! Outputs (fluxes to ccpp and host)
    real(kind_phys), dimension(:,:), intent(inout) :: &
         htrlw          ! Longwave heating-rate (K/s)
    type(sfcflw_type), dimension(:), intent(inout) :: &
         sfcflw         ! Longwave fluxes at surface (W/m2)
    type(sfcfsw_type), dimension(:), intent(inout) :: &
         sfcfsw         ! Shortwave fluxes at surface (W/m2)
    type(topflw_type), dimension(:), intent(inout) :: &
         topflw         ! Longwave fluxes at TOA (W/m2)
    type(topfsw_type), dimension(:), intent(inout) :: &
         topfsw         ! Shortwave fluxes at TOA (W/m2)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg         ! CCPP error message
    integer, intent(out) :: &
         errflg         ! CCPP error flag

    ! Locals
    logical :: top_at_1
    integer :: ipred, ilev, icase, iinf, iLay, iCol, iDay, ncol_pred
    integer, dimension(nCol) :: ipseed
    real(kind_phys) :: es, qs, dp, tem1, tem2, pfac, ranku(nCol), rankn(nCol), rankk(1), &
         bin(1), edge(nLev+1), bot, top
    real(kind_phys), dimension(nLev) :: tv, rho, tempVar
    real(kind_phys), dimension(nLev+1) :: hgtb, zo, zi
    real(kind_phys), dimension(nCol, nLev) :: o3_lay
    real(kind_phys), dimension(nCol, nLev, NF_VGAS) :: gas_vmr
    real(c_float), allocatable :: it2f(:,:,:) ! data for inference in profile, height, value order
    real(c_float), allocatable :: ot2f(:,:)   ! data from inference in profile, height order
    real(c_float), dimension(nCol, nLev, size(pnames_lw)) :: predictor_matrix_lw, upredictor_matrix_lw
    real(c_float), dimension(nDay, nLev, size(pnames_sw)) :: predictor_matrix_sw, upredictor_matrix_sw
    real(c_float), dimension(nCol, nLev+2) :: target_matrix_lw
    real(c_float), dimension(nDay, nLev+2) :: target_matrix_sw

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Are we supposed to be here?
    if (.not. do_mlrad) return

    ! What is vertical ordering of the host?
    top_at_1 = (prsi(1,1) .lt.  prsi(1, nLev))

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
    !
    ! LONGWAVE
    !
    ! ######################################################################################

    !
    ! Compute Longwave predictor matrix...
    !
    tem1 = 1._kind_phys/con_g
    tem2 = con_rd/con_g
    do iCol=1,nCol
       ! Pressure business (DJS: Will revisit after retraining)
       predictor_matrix_lw(iCol,:,1)  = prsi(iCol,1:nLev) - 1._kind_phys
       do iLay=1,nLev
          predictor_matrix_lw(iCol,iLay,21) = predictor_matrix_lw(iCol,iLay,1) - predictor_matrix_lw(iCol,iLay+1,1)
       enddo
       edge(1)      = predictor_matrix_lw(iCol,1,1)        - predictor_matrix_lw(iCol,1,21)/2.0
       edge(2:nLev) = predictor_matrix_lw(iCol,1:nLev-1,1) + predictor_matrix_lw(iCol,1:nLev-1,21)/2.0
       edge(nLev+1) = predictor_matrix_lw(iCol,nLev,1)     + predictor_matrix_lw(iCol,nLev,21)/2.0
       predictor_matrix_lw(iCol,:,21) = abs(edge(2:nlev+1)-edge(1:nlev))

       do iLay=1,nLev
          ! (DJS: Will revisit after retraining)
          ! Layer pressure (Pa)
          !predictor_matrix_lw(iCol,iLay,1)  = prsl(iCol,iLay)

          ! Pressure thickness
          !predictor_matrix_lw(iCol,iLay,21) = abs((prsi(iCol,iLay+1))-(prsi(iCol,iLay)))

          ! Thermodynamic temporaries
          es        = min( predictor_matrix_lw(iCol,iLay,1),  fpvs( tgrs(iCol,iLay) ) )
          qs        = max( con_epsqs, con_eps * es / (predictor_matrix_lw(iCol,iLay,1) + con_epsm1*es) )
          tv(iLay)  = tgrs(iCol,iLay) * (1._kind_phys + con_fvirt*qgrs(iCol,iLay,1))
          rho(iLay) = predictor_matrix_lw(iCol,iLay,1)/(con_rd*tv(iLay))

          ! Layer temperature (K)
          predictor_matrix_lw(iCol,iLay,2)  = tgrs(iCol,iLay)

          ! Layer specific-humidity (kg/kg)
          predictor_matrix_lw(iCol,iLay,3)  = qgrs(iCol,iLay,1)

          ! Compute layer relative-humidity (kg/kg)->(1)
          predictor_matrix_lw(iCol,iLay,4)  = max( 0._kind_phys, min( 1._kind_phys, max(con_epsqs, qgrs(iCol,iLay,1))/qs ) )

          ! Compute layer liquid/Ice water content (kg/kg)->(kg/m3).
          predictor_matrix_lw(iCol,iLay,5)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldliq)*rho(iLay))
          predictor_matrix_lw(iCol,iLay,6)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldice)*rho(iLay))

          ! Compute layer liquid/ice/vapor condensate path, from mixing ratios (kg/kg)->(kg/m2).
          predictor_matrix_lw(iCol,iLay,7)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldliq) * tem1 * predictor_matrix_lw(iCol,iLay,21))
          predictor_matrix_lw(iCol,iLay,8)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldice) * tem1 * predictor_matrix_lw(iCol,iLay,21))
          predictor_matrix_lw(iCol,iLay,9)  = max(0._kind_phys, qgrs(iCol,iLay,1) * tem1 * predictor_matrix_lw(iCol,iLay,21))

          ! Cloud effective radii (m).
          if (effr_in) then
             predictor_matrix_lw(iCol,iLay,13) = cld_reliq(iCol,iLay)
             predictor_matrix_lw(iCol,iLay,14) = cld_reice(iCol,iLay)
          else
             predictor_matrix_lw(iCol,iLay,13) = reliq_def*1e-6 ! (microns)->(meters)
             predictor_matrix_lw(iCol,iLay,14) = reice_def*1e-6 ! (microns)->(meters)
          endif

          ! Ozone mixing-ratio (kg/kg).
          if (i_ozone .gt. 0) then
             predictor_matrix_lw(iCol,iLay,15) = qgrs(iCol,iLay,i_ozone)
          else
             predictor_matrix_lw(iCol,iLay,15) = o3_lay(iCol,iLay)
          endif

          ! Trace gases (ppmv)
          predictor_matrix_lw(iCol,iLay,16) = gas_vmr(iCol,iLay,1)*1e6
          predictor_matrix_lw(iCol,iLay,17) = gas_vmr(iCol,iLay,3)*1e6
          predictor_matrix_lw(iCol,iLay,18) = gas_vmr(iCol,iLay,2)*1e6

       enddo ! END vertical loop
       
       ! Zenith-angle (1) Why do we need this for Longwave?
       predictor_matrix_lw(iCol,:,22) = 0._kind_phys

       ! Surface temperature (K)
       predictor_matrix_lw(iCol,:,23) = tgrs(iCol,1)

       ! Surface emissivity (1)
       predictor_matrix_lw(iCol,:,24) = semis(iCol)

       ! Compute vertically integrated (in both directions) liquid/ice/vapor condensate path (kg/m2).
       do iLay=1,nLev
          predictor_matrix_lw(iCol,iLay,10) = sum(predictor_matrix_lw(iCol,1:iLay,7))
          predictor_matrix_lw(iCol,iLay,11) = sum(predictor_matrix_lw(iCol,1:iLay,8))
          predictor_matrix_lw(iCol,iLay,12) = sum(predictor_matrix_lw(iCol,1:iLay,9))
       enddo
       do iLay=1,nLev
          predictor_matrix_lw(iCol,iLay,7)  = sum(predictor_matrix_lw(iCol,iLay:nLev,7))
          predictor_matrix_lw(iCol,iLay,8)  = sum(predictor_matrix_lw(iCol,iLay:nLev,8))
          predictor_matrix_lw(iCol,iLay,9)  = sum(predictor_matrix_lw(iCol,iLay:nLev,9))
       enddo

       ! Layer thickness and height above ground (m) (assumes SFC->TOA vertical ordering)
       ! Layer thickness (m)
       do iLay=nLev,1,-1
          predictor_matrix_lw(iCol,iLay,20) = tem2 * abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay+1))) * tv(iLay)
       enddo
       ! Height at layer boundaries
       hgtb(1) = 0._kind_phys
       do iLay=1,nLev
          hgtb(iLay+1)= hgtb(iLay) + predictor_matrix_lw(iCol,iLay,20)
          predictor_matrix_lw(iCol,iLay,19) = hgtb(iLay+1)
       enddo

       ! Height at layer centers
       !do iLay = 1, nLev
       !   pfac = abs(log(prsi(iCol,iLay)) - log(prsl(iCol,iLay))) /  &
       !        abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay+1)))
       !   predictor_matrix_lw(iCol,iLay,19) = hgtb(iLay) + pfac * (hgtb(iLay+1) - hgtb(iLay))
       !enddo
    enddo    ! END column loop

    !
    ! Normalize Longwave predictor matrix...
    !
    if (debug) upredictor_matrix_lw = predictor_matrix_lw
    do iPred = 1,size(pnames_lw)
       do iLay=1,nLev
          do iCol=1,nCol
             ! Find uniform rank of predictor value wrt training data.
             if (is2D_lw(iPred)) then
                rankk     = real(minloc(abs(predictor_matrix_lw(iCol,iLay,iPred) - &
                                            mlrad_data%lw%vector_breakpoint(:,iLay, ip2io_lw(iPred)))),kind=kind_phys)
                ncol_pred = count(mlrad_data%lw%vector_breakpoint(:, iLay, ip2io_lw(iPred)) .eq. &
                                  mlrad_data%lw%vector_breakpoint(:, iLay, ip2io_lw(iPred)))
             else
                rankk     = real(minloc(abs(predictor_matrix_lw(iCol,iLay,iPred) - &
                                            mlrad_data%lw%scalar_breakpoint(:, ip2io_lw(iPred)))),kind=kind_phys)
                ncol_pred = count(mlrad_data%lw%scalar_breakpoint(:, ip2io_lw(iPred)) .eq. &
                                  mlrad_data%lw%scalar_breakpoint(:, ip2io_lw(iPred)))
             endif
             rankk       = max(1.e-6_kind_phys,rankk)
             ranku(iCol) = rankk(1)/min(mlrad_data%lw%nCol,ncol_pred)
             predictor_matrix_lw(iCol,iLay,iPred) = sqrt(2._kind_phys)*erfinv(2.*ranku(iCol)-1)
          enddo ! END Column loop
       enddo    ! END Vertical loop
    enddo       ! END Predictor loop


    ! ######################################################################################
    !
    ! SHORTWAVE
    !
    ! ######################################################################################
    if (nDay > 0) then
       !
       ! Compute Shortwave predictor matrix...
       !  *NOTE* Many of the same fields are used for SW as in LW, just copy these over.
       !
       do iDay=1,nDay
          do iLay=1,nLev
             predictor_matrix_sw(iDay,iLay,1)  = predictor_matrix_lw(idx(iCol),iLay,1)
             predictor_matrix_sw(iDay,iLay,2)  = predictor_matrix_lw(idx(iCol),iLay,2)
             predictor_matrix_sw(iDay,iLay,3)  = predictor_matrix_lw(idx(iCol),iLay,3)
             predictor_matrix_sw(iDay,iLay,4)  = predictor_matrix_lw(idx(iCol),iLay,4)
             predictor_matrix_sw(iDay,iLay,5)  = predictor_matrix_lw(idx(iCol),iLay,5)
             predictor_matrix_sw(iDay,iLay,6)  = predictor_matrix_lw(idx(iCol),iLay,6)
             predictor_matrix_sw(iDay,iLay,7)  = predictor_matrix_lw(idx(iCol),iLay,7)
             predictor_matrix_sw(iDay,iLay,8)  = predictor_matrix_lw(idx(iCol),iLay,8)
             predictor_matrix_sw(iDay,iLay,9)  = predictor_matrix_lw(idx(iCol),iLay,9)
             predictor_matrix_sw(iDay,iLay,10) = predictor_matrix_lw(idx(iCol),iLay,10)
             predictor_matrix_sw(iDay,iLay,11) = predictor_matrix_lw(idx(iCol),iLay,11)
             predictor_matrix_sw(iDay,iLay,12) = predictor_matrix_lw(idx(iCol),iLay,12)
             predictor_matrix_sw(iDay,iLay,13) = predictor_matrix_lw(idx(iCol),iLay,13)
             predictor_matrix_sw(iDay,iLay,14) = predictor_matrix_lw(idx(iCol),iLay,14)
             predictor_matrix_sw(iDay,iLay,15) = predictor_matrix_lw(idx(iCol),iLay,15)
             predictor_matrix_sw(iDay,iLay,16) = predictor_matrix_lw(idx(iCol),iLay,16)
             predictor_matrix_sw(iDay,iLay,17) = predictor_matrix_lw(idx(iCol),iLay,17)
             predictor_matrix_sw(iDay,iLay,18) = predictor_matrix_lw(idx(iCol),iLay,18)
          enddo
       enddo

       !
       ! Normalize Shortwave predictor matrix.
       !
       do iPred = 1,size(pnames_sw)
          do iLay=1,nLev
             do iDay=1,nDay
                ! Find uniform rank of predictor value wrt training data.
                if (is2D_sw(iPred)) then
                   rankk     = real(minloc(abs(predictor_matrix_sw(iDay,iLay,iPred) - &
                                               mlrad_data%sw%vector_breakpoint(:,iLay, ip2io_sw(iPred)))),kind=kind_phys)
                   ncol_pred = count(mlrad_data%sw%vector_breakpoint(:, iLay, ip2io_sw(iPred)) .eq. &
                                     mlrad_data%sw%vector_breakpoint(:, iLay, ip2io_sw(iPred)))
                else
                   rankk     = real(minloc(abs(predictor_matrix_sw(iDay,iLay,iPred) - &
                                               mlrad_data%sw%scalar_breakpoint(:, ip2io_sw(iPred)))),kind=kind_phys)
                   ncol_pred = count(mlrad_data%sw%scalar_breakpoint(:, ip2io_sw(iPred)) .eq. &
                                     mlrad_data%sw%scalar_breakpoint(:, ip2io_sw(iPred)))
                endif
                rankk       = max(1.e-6_kind_phys,rankk)
                ranku(iDay) = rankk(1)/min(mlrad_data%sw%nCol,ncol_pred)
                predictor_matrix_sw(iDay,iLay,iPred) = sqrt(2._kind_phys)*erfinv(2.*ranku(iDay)-1)
             enddo ! END Column loop
          enddo    ! END Vertical loop
       enddo       ! END Predictor loop
       !
    endif          ! END Daylit columns


    ! ######################################################################################
    ! Begin Inference...
    ! ######################################################################################
    ! Run inferences
    ! Longwave (all columns)
    call infero_check(model%infer(predictor_matrix_lw, target_matrix_lw))

    ! Shortwave (daylit-columns only)
    !call infero_check(model%infer(predictor_matrix_sw, target_matrix_sw))

    do iCol=1,nCol
       print*,'heating profile (K/day): '
       write(*,'(a20,3a12)'  ) '                  ','ML(train)','ML(all)','G(all)'
       do iLay=1,nLev
          write(*,'(a17,i3,3f12.2)')'    ',iLay,predictor_matrix_lw(iCol,iLay,1),target_matrix_lw(iCol,iLay),htrlw(iCol,iLay)*3600.*24.
       enddo
       write(*,'(a20,3a12)'  ) '                  ','ML(all)','G(all)','G(clr)'
       write(*,'(a20,3f12.2)') 'surface(down):    ',target_matrix_lw(iCol,nLev+1),sfcflw(iCol)%dnfxc,sfcflw(iCol)%dnfx0
       write(*,'(a20,3f12.2)') 'toa(up):          ',target_matrix_lw(iCol,nLev+2),topflw(iCol)%upfxc,topflw(iCol)%upfx0
    enddo

    ! Copy from prediction-matrix to ccpp interstitials.
    do iCol=1,nCol
       htrlw(iCol,:)      = target_matrix_lw(iCol,iLay)/(3600.*24.) ! K/day -> K/sec
       sfcflw(iCol)%dnfxc = target_matrix_lw(iCol,nLev+1)
       topflw(iCol)%upfxc = target_matrix_lw(iCol,nLev+2)
       !htrsw(iCol,:)      = target_matrix_sw(iCol,iLay)/(3600.*24.) ! K/day -> K/sec 
       !sfcfsw(iCol)%dnfxc = target_matrix_sw(iCol,nLev+1)
       !topfsw(iCol)%upfxc = target_matrix_sw(iCol,nLev+2)
    enddo

    if (debug) then
       do iCol=1,nCol
          print*,'heating profile (K/day): '
          write(*,'(a20,3a12)'  ) '                  ','ML(train)','ML(all)','G(all)'
          do iLay=1,nLev
             write(*,'(a17,i3,3f12.2)')'    ',iLay,predictor_matrix_lw(iCol,iLay,1),target_matrix_lw(iCol,iLay),htrlw(iCol,iLay)*3600.*24.
          enddo
          write(*,'(a20,3a12)'  ) '                  ','ML(all)','G(all)','G(clr)'
          write(*,'(a20,3f12.2)') 'surface(down):    ',target_matrix_lw(iCol,nLev+1),sfcflw(iCol)%dnfxc,sfcflw(iCol)%dnfx0
          write(*,'(a20,3f12.2)') 'toa(up):          ',target_matrix_lw(iCol,nLev+2),topflw(iCol)%upfxc,topflw(iCol)%upfx0
       enddo

       write(*, '(a5, 25a10)')  'Layer', 'p', 'T', 'q', 'rh', 'LWC',                      &
            'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                           &
            'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                           &
            'dp', 'sza', 'Tsfc', 'sfc_emiss'
       write(*, '(a5, 25a10)')  '', '(Pa)', '(k)', '(g/kg)', '(1)', '(g/m3)',             &
            '(g/m3)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)',         &
            '(liq)','(ice)', '(mg/kg)','(ppmv)', '(ppmv)', '(ppmv)', '(m)', '(m)',        &
            '(Pa)', '(1)', '(K)', '(1)'
       do iCol=1,nCol
          do iLay=1,nLev
             write (*,'(i5,25f10.2)') iLay,                &
                  upredictor_matrix_lw(iCol,iLay,1),       &
                  upredictor_matrix_lw(iCol,iLay,2),       &
                  1.e3*upredictor_matrix_lw(iCol,iLay,3),  &
                  upredictor_matrix_lw(iCol,iLay,4),       &
                  1.e6*upredictor_matrix_lw(iCol,iLay,5),  &
                  1.e6*upredictor_matrix_lw(iCol,iLay,6),  &
                  upredictor_matrix_lw(iCol,iLay,7),       &
                  upredictor_matrix_lw(iCol,iLay,8),       &
                  upredictor_matrix_lw(iCol,iLay,9),       &
                  upredictor_matrix_lw(iCol,iLay,10),      &
                  upredictor_matrix_lw(iCol,iLay,11),      &
                  upredictor_matrix_lw(iCol,iLay,12),      &
                  1.e6*upredictor_matrix_lw(iCol,iLay,13), &
                  1.e6*upredictor_matrix_lw(iCol,iLay,14), &
                  1.e6*upredictor_matrix_lw(iCol,iLay,15), &
                  upredictor_matrix_lw(iCol,iLay,16),      &
                  upredictor_matrix_lw(iCol,iLay,17),      &
                  upredictor_matrix_lw(iCol,iLay,18),      &
                  upredictor_matrix_lw(iCol,iLay,19),      &
                  upredictor_matrix_lw(iCol,iLay,20),      &
                  upredictor_matrix_lw(iCol,iLay,21),      &
                  upredictor_matrix_lw(iCol,iLay,22),      &
                  upredictor_matrix_lw(iCol,iLay,23),      &
                  upredictor_matrix_lw(iCol,iLay,24)
          enddo
       enddo
       write (*,'(a50)') '################################'
       write(*, '(a5, 25a10)')  'Layer', 'p', 'T', 'q', 'rh', 'LWC',   &
            'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',        &
            'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',        &
            'dp', 'sza', 'Tsfc', 'sfc_emiss'
       do iCol=1,nCol
          do iLay=1,nLev
             write (*,'(i5,25f10.2)') iLay,                &
                  predictor_matrix_lw(iCol,iLay,1),        &
                  predictor_matrix_lw(iCol,iLay,2),        &
                  predictor_matrix_lw(iCol,iLay,3),        &
                  predictor_matrix_lw(iCol,iLay,4),        &
                  predictor_matrix_lw(iCol,iLay,5),        &
                  predictor_matrix_lw(iCol,iLay,6),        &
                  predictor_matrix_lw(iCol,iLay,7),        &
                  predictor_matrix_lw(iCol,iLay,8),        &
                  predictor_matrix_lw(iCol,iLay,9),        &
                  predictor_matrix_lw(iCol,iLay,10),       &
                  predictor_matrix_lw(iCol,iLay,11),       &
                  predictor_matrix_lw(iCol,iLay,12),       &
                  predictor_matrix_lw(iCol,iLay,13),       &
                  predictor_matrix_lw(iCol,iLay,14),       &
                  predictor_matrix_lw(iCol,iLay,15),       &
                  predictor_matrix_lw(iCol,iLay,16),       &
                  predictor_matrix_lw(iCol,iLay,17),       &
                  predictor_matrix_lw(iCol,iLay,18),       &
                  predictor_matrix_lw(iCol,iLay,19),       &
                  predictor_matrix_lw(iCol,iLay,20),       &
                  predictor_matrix_lw(iCol,iLay,21),       &
                  predictor_matrix_lw(iCol,iLay,22),       &
                  predictor_matrix_lw(iCol,iLay,23),       &
                  predictor_matrix_lw(iCol,iLay,24)
          enddo
       enddo
    endif

  end subroutine mlrad_driver_run

! #########################################################################################
!! \section arg_table_mlrad_driver_finalize
!! \htmlinclude mlrad_driver_finalize.html
!!
! #########################################################################################
  subroutine mlrad_driver_finalize(do_mlrad, errmsg, errflg)
    ! Inputs
    logical,           intent(in) :: do_mlrad
    ! Outputs
    character(len=*),  intent(out) :: errmsg
    integer,           intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_mlrad) return

    ! Free the model
    call infero_check(model%free())

    ! Finalize
    call infero_check(infero_finalise())

    !close(98)

  end subroutine mlrad_driver_finalize

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

end module mlrad_driver
