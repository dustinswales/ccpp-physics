! ###########################################################################################
!> \file ml_rad_phys.F90
!!
!! This module includes the ONNX based machine-learning emulator for longwave radiation detailed
!! in 10.22541/essoar.168319865.58439449/v1
!!
!! There are 24 predictors for this emulator:
!!         Name                       Units
!!    1  - 'pressure'                 [Pa]
!!    2  - 'temperature'              [K]
!!    3  - 'specific_humidity'        [kg/kg]
!!    4  - 'relative_humidity'        [1]
!!    5  - 'liquid_water_content'     [kg/m3]
!!    6  - 'ice_water_content'        [kg/m3]
!!    7  - 'liquid_water_path'        [kg/m2]
!!    8  - 'ice_water_path'           [kg/m2]
!!    9  - 'vapour_path'              [kg/m2]
!!    10 - 'upward_liquid_water_path' [kg/m2]
!!    11 - 'upward_ice_water_path'    [kg/m2]
!!    12 - 'upward_vapour_path'       [kg/m2]
!!    13 - 'liquid_effective_radius'  [m]
!!    14 - 'ice_effective_radius'     [m]
!!    15 - 'o3_mixing_ratio'          [kg/kg]
!!    16 - 'co2_concentration'        [ppmv]
!!    17 - 'ch4_concentration'        [ppmv]
!!    18 - 'n2o_concentration'        [ppmv]
!!    19 - 'height_agl'               [m]
!!    20 - 'height_thickness'         [m]
!!    21 - 'pressure_thickness'       [Pa]
!!    22 - 'zenith_angle'             [1]
!!    23 - 'surface_temperature'      [K]
!!    24 - 'surface_emissivity'       [1]
!!
! ###########################################################################################
module ml_rad_phys
  use machine,  only: kind_phys
  use funcphys, only: fpvs
  use module_radiation_gases, only: NF_VGAS, getgases, getozn
  use netcdf
  use inferof
  use iso_c_binding, only : c_double, c_int, c_float, c_char, c_null_char, c_ptr
  implicit none

  type(infero_model) :: model

  real (kind_phys), parameter :: &
       reliq_def  = 10.0 ,  & ! Default liq radius to 10 micron (used when effr_in=F)
       reice_def  = 50.0,   & ! Default ice radius to 50 micron (used when effr_in=F)
       rerain_def = 1000.0, & ! Default rain radius to 1000 micron (used when effr_in=F)
       resnow_def = 250.0     ! Default snow radius to 250 micron (used when effr_in=F)

  public ml_rad_phys_init, ml_rad_phys_run

contains
! ########################################################################################
!! \section arg_table_ml_rad_phys_init
!! \htmlinclude ml_rad_phys_init.html
!!
! #########################################################################################
  subroutine ml_rad_phys_init(do_ml_rad, infero_model_path, infero_model_type, errmsg, errflg)
    ! Inputs
    logical,            intent(in) :: do_ml_rad
    character(len=128), intent(in) :: infero_model_path
    character(len=128), intent(in) :: infero_model_type

    ! Outputs
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Locals
    character(1024) :: yaml_config

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ml_rad) return

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
  subroutine ml_rad_phys_run(do_ml_rad, effr_in, debug, nCol, nLev, i_cldliq, i_cldice, i_ozone, ico2, semis, lon, lat,&
       prsl, tgrs, prslk, prsi, cld_reliq, cld_reice, qgrs, predictor_matrix, vector_prediction_matrix,       &
       scalar_prediction_matrix, con_epsqs, con_eps, con_epsm1, con_rd, con_fvirt, con_g, con_pi, errmsg, errflg)
    ! Inputs
    logical, intent(in) :: &
         do_ml_rad,         & ! Use ML emulator for LW radiation?
         effr_in,           & ! Provide hydrometeor radii from macrophysics? 
         debug                ! Debug mode?
    integer, intent(in) ::  &
         nCol,              & ! Number of horizontal grid points
         nLev,              & ! Number of vertical layers
         i_cldliq,          & ! Index into tracer array for cloud liquid.
         i_cldice,          & ! Index into tracer array for cloud ice.
         i_ozone,           & ! Index into tracer array for ozone concentration.
         ico2                 ! Flag for co2 radiation scheme
    real(kind_phys), dimension(:), intent(in) :: &
         lon,               & ! Longitude
         lat,               & ! Latitude
         semis                ! Longwave surface emissivity
    real(kind_phys), dimension(:,:), intent(in) :: &
         prsl,              & ! Pressure at model-layer centers (Pa)
         tgrs,              & ! Temperature at model-layer centers (K)
         prslk,             & ! Exner function at model layer centers (1)
         prsi,              & ! Pressure at model-interfaces (Pa)
         cld_reliq,         & ! Effective radius (m)
         cld_reice            ! Effective radius (m)
    real(kind_phys), dimension(:,:,:), intent(in) :: &
         qgrs                 ! Tracer concentrations (kg/kg)
    real(kind_phys), intent(in) :: &
         con_epsqs,         & ! Physical constant: Minimum saturation mixing-ratio (kg/kg)
         con_eps,           & ! Physical constant: Epsilon (Rd/Rv)
         con_epsm1,         & ! Physical constant: Epsilon (Rd/Rv) minus one
         con_rd,            & ! Physical constant: gas-constant for dry air
         con_fvirt,         & ! Physical constant: Inverse of epsilon minus one
         con_g,             & ! Physical constant: gravitational constant
         con_pi               ! Physical constant: Pi

    real(kind_phys), intent(in), dimension(:,:,:) :: &
         predictor_matrix, &
         vector_prediction_matrix
    real(kind_phys),intent(in), dimension(:,:) :: &
         scalar_prediction_matrix

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg

    ! Locals
    integer :: npred_ml, nlev_ml, ncase_ml, ipred, ilev, icase, iinf, iLay, iCol
    real(c_float), allocatable :: it2f(:,:,:) ! data for inference in profile, height, value order
    real(c_float), allocatable :: ot2f(:,:)   ! data from inference in profile, height order

    logical :: top_at_1
    integer :: iSFC, iTOA, stride
    real(kind_phys), dimension(nCol, nLev, 24) :: predictor_array
    real(kind_phys) :: es, qs, dp, tem1, tem2, pfac
    real(kind_phys), allocatable :: o3_lay(:,:)
    real(kind_phys), dimension(nLev) :: tv, rho
    real(kind_phys), dimension(nLev+1) :: hgtb
    !real(kind_phys), dimension(nCol, nLev) :: 
    real(kind_phys), dimension(nCol, nLev, NF_VGAS) :: gas_vmr

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ml_rad) return

    ! Dimensions
    npred_ml = size(predictor_matrix(:,0,0))
    nlev_ml  = size(predictor_matrix(0,:,0))
    ncase_ml = size(predictor_matrix(0,0,:))

    ! What is vertical ordering?
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

    ! Do we need to get climatological ozone? (only if not using prognostic ozone)
    if (i_ozone .le. 0) then
       allocate(o3_lay(nCol,nLev))
       call getozn (prslk, lat, nCol, nLev, top_at_1, o3_lay)
    endif

    ! Get trace-gas concentrations.
    call getgases (prsi/100., lon, lat, nCol, nLev, ico2, top_at_1, con_pi, gas_vmr)

    ! ######################################################################################
    ! Compute fields for predictor matrix and populate predictor matrix.
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
          dp        = abs(prsl(iCol,iLay+1)-prsl(iCol,iLay))

          ! Layer pressure (Pa)
          predictor_array(iCol,iLay,1)  = prsl(iCol,iLay)

          ! Layer temperature (K)
          predictor_array(iCol,iLay,2)  = tgrs(iCol,iLay)

          ! Layer specific-humidity (kg/kg)
          predictor_array(iCol,iLay,3)  = qgrs(iCol,iLay,1)

          ! Compute layer relative-humidity (kg/kg)->(1)
          predictor_array(iCol,iLay,4)  = max( 0._kind_phys, min( 1._kind_phys, max(con_epsqs, qgrs(iCol,iLay,1))/qs ) )

          ! Compute layer liquid/Ice water content (kg/kg)->(kg/m3).
          predictor_array(iCol,iLay,5)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldliq)*rho(iLay))
          predictor_array(iCol,iLay,6)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldice)*rho(iLay))

          ! Compute layer liquid/ice/vapor condensate path, from mixing ratios (kg/kg)->(kg/m2).
          predictor_array(iCol,iLay,7)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldliq) * tem1 * dp)
          predictor_array(iCol,iLay,8)  = max(0._kind_phys, qgrs(iCol,iLay,i_cldice) * tem1 * dp)
          predictor_array(iCol,iLay,9)  = max(0._kind_phys, qgrs(iCol,iLay,1) * tem1 * dp)

          ! Cloud effective radii (m).
          if (effr_in) then
             predictor_array(iCol,iLay,13) = cld_reliq(iCol,iLay)
             predictor_array(iCol,iLay,14) = cld_reice(iCol,iLay)
          else
             predictor_array(iCol,iLay,13) = reliq_def*1e-6 ! (microns)->(meters)
             predictor_array(iCol,iLay,14) = reice_def*1e-6 ! (microns)->(meters)
          endif

          ! Ozone mixing-ratio (kg/kg).
          if (i_ozone .gt. 0) then
             predictor_array(iCol,iLay,15) = qgrs(iCol,iLay,i_ozone)
          else
             predictor_array(iCol,iLay,15) = o3_lay(iCol,iLay)
          endif

          ! Trace gases (ppmv)
          predictor_array(iCol,iLay,16) = gas_vmr(iCol,iLay,1)
          predictor_array(iCol,iLay,17) = gas_vmr(iCol,iLay,3)
          predictor_array(iCol,iLay,18) = gas_vmr(iCol,iLay,2)

          ! Pressure-thickness (Pa)
          predictor_array(iCol,iLay,21) = dp

          ! Zenith-angle (1) Why do we need this for Longwave?
          predictor_array(iCol,:,22) = 0._kind_phys

          ! Surface temperature (K)
          predictor_array(iCol,:,23) = tgrs(iCol,iSFC)

          ! Surface emissivity (1)
          predictor_array(iCol,:,24) = semis(iCol)
       enddo ! END vertical loop

       ! Compute vertically integrated liquid/ice/vapor condensate path (kg/m2).
       if (top_at_1) then
          do iLay=1,nLev
             predictor_array(iCol,iLay,10) = sum(predictor_array(iCol,1:iLay,7))
             predictor_array(iCol,iLay,11) = sum(predictor_array(iCol,1:iLay,8))
             predictor_array(iCol,iLay,12) = sum(predictor_array(iCol,1:iLay,9))
          enddo
       else
          do iLay=1,nLev
             predictor_array(iCol,iLay,10) = sum(predictor_array(iCol,iLay:nLev,7))
             predictor_array(iCol,iLay,11) = sum(predictor_array(iCol,iLay:nLev,8))
             predictor_array(iCol,iLay,12) = sum(predictor_array(iCol,iLay:nLev,9))
          enddo
       endif

       ! Layer thickness and height above ground (m)
       if (top_at_1) then
          ! Layer thickness (m)
          do iLay=1,nLev
             predictor_array(iCol,iLay,20) = tem2 * abs(log(prsi(iCol,iLay+1)) - log(prsi(iCol,iLay))) * tv(iLay)
          enddo
          ! Height at layer boundaries
          hgtb(nLev+1) = 0._kind_phys
          do iLay=nLev,1,-1
             hgtb(iLay)= hgtb(iLay+1) + predictor_array(iCol,iLay,20)
          enddo
          ! Height at layer centers
          do iLay = nLev, 1, -1
             pfac = abs(log(prsi(iCol,iLay+1)) - log(prsi(iCol,iLay))) /  &
                    abs(log(prsi(iCol,iLay+1)) - log(prsi(iCol,iLay)))
             predictor_array(iCol,iLay,19) = hgtb(iLay+1) + pfac * (hgtb(iLay) - hgtb(iLay+1))
          enddo
       else
          ! Layer thickness (m)
          do iLay=nLev,1,-1
             predictor_array(iCol,iLay,20) = tem2 * abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay+1))) * tv(iLay)
          enddo
          ! Height at layer boundaries
          hgtb(1) = 0._kind_phys
          do iLay=1,nLev
             hgtb(iLay+1)= hgtb(iLay) + predictor_array(iCol,iLay,20)
          enddo
          ! Height at layer centers
          do iLay = 1, nLev
             pfac = abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay)  )) /  &
                    abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay+1)))
             predictor_array(iCol,iLay,19) = hgtb(iLay) + pfac * (hgtb(iLay+1) - hgtb(iLay))
          enddo
       endif
    enddo ! END column loop

    ! ######################################################################################
    ! Debugging information
    ! ######################################################################################
    if (debug) then
       write(*, '(26a12)') 'Layer', 'p(Pa)', 'T(k)', 'q(kg/kg)', 'rh(1)', 'LWC(mg/m3)',    &
            'IWC(mg/m3)', 'LWP(g/m2)', 'IWP(g/m2)', 'WVP(g/m2)', 'iLWP(g/m2)',             &
            'iIWP(g/m2)', 'iWVP(g/m2)', 'Reff(liq)','Reff(ice)', 'o3(mg/kg)','co2(ppmv)',  &
            'ch4(ppmv)', 'n2o(ppmv)', 'z(m)', 'dz(m)', 'dp(Pa)', 'sza', 'Tsfc(K)',         &
            'sfc_emiss(1)', 'rho(kg/m3)'
       do iCol=1,nCol
          do iLay=1,nLev
             write (*,'(i12,26f12.4)') iLay,          &
                  predictor_array(iCol,iLay,1),       &
                  predictor_array(iCol,iLay,2),       &
                  predictor_array(iCol,iLay,3),       &
                  predictor_array(iCol,iLay,4),       &
                  1.e6*predictor_array(iCol,iLay,5),  &
                  1.e6*predictor_array(iCol,iLay,6),  &
                  predictor_array(iCol,iLay,7),       &
                  predictor_array(iCol,iLay,8),       &
                  predictor_array(iCol,iLay,9),       &
                  predictor_array(iCol,iLay,10),      &
                  predictor_array(iCol,iLay,11),      &
                  predictor_array(iCol,iLay,12),      &
                  1.e6*predictor_array(iCol,iLay,13), &
                  1.e6*predictor_array(iCol,iLay,14), &
                  1.e6*predictor_array(iCol,iLay,15), &
                  1.e6*predictor_array(iCol,iLay,16), &
                  1.e6*predictor_array(iCol,iLay,17), &
                  1.e6*predictor_array(iCol,iLay,18), &
                  predictor_array(iCol,iLay,19),      &
                  predictor_array(iCol,iLay,20),      &
                  predictor_array(iCol,iLay,21),      &
                  predictor_array(iCol,iLay,22),      &
                  predictor_array(iCol,iLay,23),      &
                  predictor_array(iCol,iLay,24),      &
                  rho(iLay)
          enddo
       enddo
    endif

    ! ######################################################################################
    ! Begin Inference...
    ! ######################################################################################

    ! Reverse order of input predictor matrix for inference.
    allocate(it2f(ncase_ml, nlev_ml, npred_ml))
    do ipred=1,npred_ml
       do ilev=1,nlev_ml
          do icase=1,ncase_ml
             it2f(icase, ilev, ipred) = predictor_matrix(ipred, ilev, icase)
          enddo
       enddo
    enddo

    ! Run inferences
    allocate(ot2f(ncase_ml, nlev_ml + 2))
    call infero_check(model%infer(it2f, ot2f ))

    ! ######################################################################################
    ! Check the data
    ! ######################################################################################
    if (debug) then
       do icase = 1, ncase_ml
          do ilev = 1, nlev_ml + 2
             do ipred = 1, 1
                if (ilev <= nlev_ml) then
                   if (abs(ot2f(icase,ilev) - vector_prediction_matrix(ipred,ilev,icase)) .gt. 1e-3) then
                      write(*,*) "ERROR: output element ",icase,ilev, " (", ot2f(icase,ilev) ,") ", &
                           "is different from expected value ", vector_prediction_matrix(ipred,ilev,icase)
                      stop 1
                   end if
                   print *, "( ",icase , ", ", ilev, ", ", ipred,") = ", ot2f(icase, ilev), vector_prediction_matrix(ipred, ilev, icase)
                else
                   if (abs(ot2f(icase,ilev) - scalar_prediction_matrix(ilev-127,icase)) .gt. 1e-3) then
                      write(*,*) "ERROR: output element ",icase,ilev, " (", ot2f(icase,ilev) ,") ", &
                           "is different from expected value ", scalar_prediction_matrix(ilev-127,icase)
                      stop 1
                   end if
                   print *,  "( ",icase , ", ", ilev, ", ", ipred,") = ", ot2f(icase, ilev), scalar_prediction_matrix(ilev-127, icase)
                end if
             end do
          end do
       end do
    endif

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
  
end module ml_rad_phys
