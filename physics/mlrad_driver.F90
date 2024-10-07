! ###########################################################################################
!> \file mlrad_river.F90
!!
!! This module includes the ONNX based machine-learning emulator for longwave and shortwave
!! radiation detailed in 10.22541/essoar.168319865.58439449/v1
!!
!! There are 24 predictors for the Longwave (LW) emulator:
!!         Name                                Units
!!    1  - 'zenith_angle_radians'              [1]
!!    2  - 'surface_temperature_kelvins'       [K]
!!    3  - 'surface_emissivity'                [1]  
!!    4  - 'pressure_pascals'                  [Pa]
!!    5  - 'temperature_kelvins'               [K]
!!    6  - 'specific_humidity_kg_kg01'         [kg/kg]
!!    7  - 'relative_humidity_unitless'        [1]
!!    8  - 'liquid_water_content_kg_m03'       [kg/m3]
!!    9  - 'ice_water_content_kg_m03'          [kg/m3]
!!    10 - 'downward_liquid_water_path_kg_m02' [kg/m2]
!!    11 - 'downward_ice_water_path_kg_m02'    [kg/m2]
!!    12 - 'downward_vapour_path_kg_m02'       [kg/m2]
!!    13 - 'upward_liquid_water_path_kg_m02'   [kg/m2]
!!    14 - 'upward_liquid_ice_path_kg_m02'     [kg/m2]
!!    15 - 'upward_vapour_path_kg_m02'         [kg/m2]
!!    16 - 'liquid_effective_radius_metres'    [m]
!!    17 - 'ice_effective_radius_metres'       [m]
!!    18 - 'o3_mixing_ratio_kg_kg01'           [kg/kg]
!!    19 - 'co2_concentration_ppmv'            [ppmv]
!!    20 - 'ch4_concentration_ppmv'            [ppmv]
!!    21 - 'n2o_concentration_ppmv'            [ppmv]
!!    22 - 'height_m_agl'                      [m]
!!    23 - 'height_thickness_metres'           [m]
!!    24 - 'pressure_thickness_pascals'        [Pa]
!!
!! There are 26 predictors for the Shortwave (SW) emulator:
!!         Name                                Units
!!    1  - 'zenith_angle_radians'              [1]
!!    2  - 'albedo'                            [1]
!!    3  - 'aerosol_single_scattering_albedo'  [1]
!!    4  - 'aerosol_asymmetry_param'           [1] 
!!    5  - 'pressure_pascals'                  [Pa]
!!    6  - 'temperature_kelvins'               [K]
!!    7  - 'specific_humidity_kg_kg01'         [kg/kg]
!!    8  - 'relative_humidity_unitless'        [1]
!!    9  - 'liquid_water_content_kg_m03'       [kg/m3]
!!    10 - 'ice_water_content_kg_m03'          [kg/m3]
!!    11 - 'downward_liquid_water_path_kg_m02' [kg/m2]
!!    12 - 'downward_ice_water_path_kg_m02'    [kg/m2]
!!    13 - 'downward_vapour_path_kg_m02'       [kg/m2]
!!    14 - 'upward_liquid_water_path_kg_m02'   [kg/m2]
!!    15 - 'upward_liquid_ice_path_kg_m02'     [kg/m2]
!!    16 - 'upward_vapour_path_kg_m02'         [kg/m2]
!!    17 - 'liquid_effective_radius_metres'    [m]
!!    18 - 'ice_effective_radius_metres'       [m]
!!    19 - 'o3_mixing_ratio_kg_kg01'           [kg/kg]
!!    20 - 'co2_concentration_ppmv'            [ppmv]
!!    21 - 'ch4_concentration_ppmv'            [ppmv]
!!    22 - 'n2o_concentration_ppmv'            [ppmv]
!!    23 - 'aerosol_extinction_metres01'       [1]
!!    24 - 'height_m_agl'                      [m]
!!    25 - 'height_thickness_metres'           [m]
!!    26 - 'pressure_thickness_pascals'        [Pa]
!!
! ###########################################################################################
module mlrad_driver
  use machine,                   only: kind_phys, kind_dbl_prec
  use funcphys,                  only: fpvs
  use module_radiation_gases,    only: NF_VGAS, getgases
  use module_radiation_aerosols, only: setaer
  use module_mlrad,              only: ty_mlrad_data, ip2io_lw, is2D_lw, pnames_lw,         &
                                       pnames_sw, ip2io_sw, is2D_sw,ty_rad_ml_ref_data,     &
                                       ilw_sza, ilw_sfct, ilw_emiss, ilw_p, ilw_t, ilw_q,   &
                                       ilw_rh, ilw_lwc, ilw_iwc, ilw_dlwp, ilw_diwp,        &
                                       ilw_dwvp, ilw_ulwp, ilw_uiwp, ilw_uwvp, ilw_reliq,   &
                                       ilw_reice, ilw_o3mr, ilw_co2, ilw_ch4, ilw_n2o,      &
                                       ilw_z, ilw_dz, ilw_dp, isw_sza, isw_alb, isw_aerssa, &
                                       isw_aerasy, isw_p, isw_t, isw_q, isw_rh, isw_lwc,    &
                                       isw_iwc, isw_dlwp, isw_diwp, isw_dwvp, isw_ulwp,     &
                                       isw_uiwp, isw_uwvp, isw_reliq, isw_reice, isw_o3mr,  &
                                       isw_co2, isw_ch4, isw_n2o, isw_tauaer, isw_z, isw_dz,&
                                       isw_dp
  use iso_c_binding,             only: c_float, c_null_char
  use inferof
  use netcdf

  implicit none

  real(kind_phys),dimension(1,127,1,28) :: npmatrix_offline

  type(infero_model) :: model_lw, model_sw

  logical :: do_debug_once

  ! Default effective radii, used when coupling radiation to single-moment cloud microphysics
  ! This is set by the host using, effr_in=F
  real (kind_phys), parameter :: &
       reliq_def  = 10.0, & ! Default liquid radius: 10 microns
       reice_def  = 50.0, & ! Default ice radius:  50 microns
       epsm1      = 0.999999, eps=0.000001
!       eps        = epsilon(1._kind_phys), &
!       epsm1      = 1._kind_phys - epsilon(1._kind_phys)
  integer :: cmode,ncid,status,varID0,varID1,varID2,varID3,varID4,varID5,varID6,varID7,     &
       varID8,varID9,varID10,varID11,varID12,varID13,varID14,varID15,varID16,varID17,       &
       varID18,varID19,varID20,varID21,dimID1,dimID2,record_counter(1),varID22,varID23,     &
       varID24,varID25,varID26,varID27,varID28

  
  public mlrad_driver_init, mlrad_driver_run

contains
! ###########################################################################################
!! \section arg_table_mlrad_driver_init
!! \htmlinclude mlrad_driver_init.html
!!
! ###########################################################################################
  subroutine mlrad_driver_init(do_mlrad, debug, infero_mpath_lw, infero_mtype_lw,           &
       infero_mpath_sw, infero_mtype_sw, mlrad_data, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         debug,           & ! Debug mode
         do_mlrad           ! Use ML emulator for radiation?
    character(len=128), intent(in) :: &
         infero_mpath_lw, & ! Infero model (LW) path
         infero_mtype_lw, & ! Infero model (LW) type
         infero_mpath_sw, & ! Infero model (SW) path
         infero_mtype_sw    ! Infero model (SW) type 
    type(ty_mlrad_data), intent(inout) :: &
         mlrad_data         ! DDT containing training data
    
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg          ! CCPP error message
    integer, intent(out) :: &
         errflg          ! CCPP error flag

    ! Locals
    character(1024) :: yaml_config_lw, yaml_config_sw
    integer :: iCol, iPred, iLay, iName, count

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_mlrad) return

    status = nf90_create('debug.mlrad_driver.fluxes.nc',cmode=nf90_clobber,ncid=ncid)
    status = nf90_def_dim(ncid,'lev',127,dimID1)
    status = nf90_def_dim(ncid,'time',nf90_unlimited,dimID2)
    status = nf90_def_var(ncid,'time_rad',     nf90_real, dimID2, varID0)
    status = nf90_def_var(ncid,'lwdnsfc_rrtmg',nf90_real, dimID2, varID1)
    status = nf90_def_var(ncid,'lwuptoa_rrtmg',nf90_real, dimID2, varID2)
    status = nf90_def_var(ncid,'lwdnsfc_mlrad',nf90_real, dimID2, varID3)
    status = nf90_def_var(ncid,'lwuptoa_mlrad',nf90_real, dimID2, varID4)
    status = nf90_def_var(ncid,'lwc',          nf90_real, (/dimID1,dimID2/), varID5)
    status = nf90_def_var(ncid,'iwc',          nf90_real, (/dimID1,dimID2/), varID6)
    status = nf90_def_var(ncid,'z',            nf90_real, (/dimID1,dimID2/), varID7)
    status = nf90_def_var(ncid,'dlwp',         nf90_real, (/dimID1,dimID2/), varID8)
    status = nf90_def_var(ncid,'diwp',         nf90_real, (/dimID1,dimID2/), varID9)
    status = nf90_def_var(ncid,'dwvp',         nf90_real, (/dimID1,dimID2/), varID10)
    status = nf90_def_var(ncid,'qv',           nf90_real, (/dimID1,dimID2/), varID11)
    status = nf90_def_var(ncid,'rho',          nf90_real, (/dimID1,dimID2/), varID12)
    status = nf90_def_var(ncid,'ulwp',         nf90_real, (/dimID1,dimID2/), varID13)
    status = nf90_def_var(ncid,'uiwp',         nf90_real, (/dimID1,dimID2/), varID14)
    status = nf90_def_var(ncid,'uwvp',         nf90_real, (/dimID1,dimID2/), varID15)
    status = nf90_def_var(ncid,'mr_liq',       nf90_real, (/dimID1,dimID2/), varID16)
    status = nf90_def_var(ncid,'mr_ice',       nf90_real, (/dimID1,dimID2/), varID17)
    status = nf90_def_var(ncid,'mr_rain',      nf90_real, (/dimID1,dimID2/), varID18)
    status = nf90_def_var(ncid,'mr_snow',      nf90_real, (/dimID1,dimID2/), varID19)
    status = nf90_def_var(ncid,'mr_grpl',      nf90_real, (/dimID1,dimID2/), varID20)
    status = nf90_def_var(ncid,'mr_tcld',      nf90_real, (/dimID1,dimID2/), varID21)
    status = nf90_def_var(ncid,'dp',           nf90_real, (/dimID1,dimID2/), varID22)
    status = nf90_def_var(ncid,'dz',  	       nf90_real, (/dimID1,dimID2/), varID23)
    status = nf90_def_var(ncid,'p',            nf90_real, (/dimID1,dimID2/), varID24)

    status = nf90_enddef(ncid)
    record_counter(1) = 1    
    if (debug) then
       open(93, file='debug.mlrad_driver.heating_rates.txt',   status='unknown')
       open(94, file='debug.mlrad_driver.fluxes.txt',          status='unknown')
       write(94,'(4a18)'  ) 'ML(LW down @ SFC)','G(LW down @ SFC)','ML(LW up @ TOA)', 'G(LW up @ TOA)'
       open(95, file='debug.mlrad_driver.upmatrix.example.txt',status='unknown')
       open(96, file='debug.mlrad_driver.npmatrix.example.txt',status='unknown')
       open(97, file='debug.mlrad_driver.upmatrix.inline.txt', status='unknown')
       open(98, file='debug.mlrad_driver.npmatrix.inline.txt', status='unknown')
       open(99, file='debug.mlrad_driver.breakpoints.txt',     status='unknown')
       do_debug_once = .true.
    endif

    ! #######################################################################################
    !
    ! Determine the mapping between the fields in training data and the order of expected 
    ! inputs to the emulator. Store in "ip2io"
    ! Also, store which data container the predictors are in (1D vs 2D)
    ! (The emulator expects the predictors in order. The inputs may not be in this order)
    !
    ! #######################################################################################
    ! Longwave
    count = 0
    ip2io_lw(:) = -999
    do iName = 1,size(pnames_lw)
       ! Scalars
       do iPred= 1,mlrad_data%lw%npreds
          if (trim(pnames_lw(iName)) == trim(mlrad_data%lw%scalar_pnames(iPred))) then
             ip2io_lw(iName) = iPred
             is2D_lw(iName)  = .false.
             count = count + 1
          endif
       enddo
       ! Vectors
       do iPred = 1,mlrad_data%lw%npredv
          if (trim(pnames_lw(iName)) == trim(mlrad_data%lw%vector_pnames(iPred))) then
             ip2io_lw(iName) = iPred
             is2D_lw(iName)  = .true.
             count = count + 1
          endif
       enddo
    enddo

    ! Shortwave
    count = 0
    ip2io_sw(:) = -999
    do iName = 1,size(pnames_sw)
       ! Scalars
       do iPred= 1,mlrad_data%sw%npreds
          if (trim(pnames_sw(iName)) == trim(mlrad_data%sw%scalar_pnames(iPred))) then
             ip2io_sw(iName) = iPred
             is2D_sw(iName)  = .false.
             count = count + 1
          endif
       enddo
       ! Vectors
       do iPred = 1,mlrad_data%sw%npredv
          if (trim(pnames_sw(iName)) == trim(mlrad_data%sw%vector_pnames(iPred))) then
             ip2io_sw(iName) = iPred
             is2D_sw(iName)  = .true.
             count = count + 1
          endif
       enddo
    enddo

    ! #######################################################################################
    !
    ! Initialize Infero
    !
    ! #######################################################################################
    call infero_check(infero_initialise())

    ! Longwave
    yaml_config_lw = "---"//NEW_LINE('A') &
         //"  path: "//TRIM(infero_mpath_lw)//NEW_LINE('A') &
         //"  type: "//TRIM(infero_mtype_lw)//c_null_char
    !
    call infero_check(model_lw%initialise_from_yaml_string(yaml_config_lw))

    ! Shortwave
    yaml_config_sw = "---"//NEW_LINE('A') &
         //"  path: "//TRIM(infero_mpath_sw)//NEW_LINE('A') &
         //"  type: "//TRIM(infero_mtype_sw)//c_null_char
    !
    call infero_check(model_sw%initialise_from_yaml_string(yaml_config_sw))

  end subroutine mlrad_driver_init
! ###########################################################################################
!! \section arg_table_mlrad_driver_run
!! \htmlinclude mlrad_driver_run.html
!!
! ###########################################################################################
  subroutine mlrad_driver_run(do_mlrad, effr_in, debug, do_norm, nCol, nLev, nDay, ntrac, i_cldliq,&
       i_cldice, i_ozone, ico2, isubc, iaermdl, iaerflg, icseed, idx, lsmask, semis, sfcalb,&
       coszen, lon, lat, prsl, tgrs, prslk, prsi, cld_reliq, cld_reice, qgrs, aerfld,       &
       mlrad_data, con_epsqs, con_eps, con_epsm1, con_rd, con_fvirt, con_g, con_pi, htrlw,  &
       htrsw, sfcflw, sfcfsw, topflw, topfsw, oro, fhour, errmsg, errflg, ref_data)
    use module_radsw_parameters, only: topfsw_type, sfcfsw_type
    use module_radlw_parameters, only: topflw_type, sfcflw_type

    ! Inputs
    type(ty_rad_ml_ref_data),intent(in) :: ref_data
    type(ty_mlrad_data), intent(in) :: &
         mlrad_data     ! DDT containing training data
    logical, intent(in) :: &
         do_mlrad,    & ! Use ML emulator for LW radiation?
         effr_in,     & ! Provide hydrometeor radii from macrophysics? 
         debug,       & ! Debug mode?
         do_norm        ! Normalize (LW) predictor matrix?
    integer, intent(in) ::  &
         nCol,        & ! Number of horizontal grid points
         nLev,        & ! Number of vertical layers
         nDay,        & ! Number of daylit columns
         ntrac,       & ! Number of tracers
         i_cldliq,    & ! Index into tracer array for cloud liquid.
         i_cldice,    & ! Index into tracer array for cloud ice.
         i_ozone,     & ! Index into tracer array for ozone concentration.
         ico2,        & ! Flag for co2 radiation scheme
         isubc,       & ! Flag for cloud-seeding (rng) for cloud-sampling
         iaermdl,     & ! Aerosol model scheme flag
         iaerflg        ! Aerosol effects to include
    integer,intent(in),dimension(:) :: &
         icseed,      & ! Seed for random number generation for longwave radiation
         idx            ! Index array for daytime points
    real(kind_phys), intent(in) :: &
         fhour          ! Forecast time in hours
    real(kind_phys), dimension(:), intent(in) :: &
         lsmask,      & ! Land/sea/sea-ice mask
         lon,         & ! Longitude
         lat,         & ! Latitude
         semis,       & ! Longwave surface emissivity
         coszen,      & ! Cosine(SZA)
         oro
    real(kind_phys), dimension(:,:), intent(in) :: &
         prsl,        & ! Pressure at model-layer centers (Pa)
         tgrs,        & ! Temperature at model-layer centers (K)
         prslk,       & ! Exner function at model layer centers (1)
         prsi,        & ! Pressure at model-interfaces (Pa)
         cld_reliq,   & ! Effective radius (m)
         cld_reice,   & ! Effective radius (m)
         sfcalb         ! Surface albedo
    real(kind_phys), dimension(:,:,:), intent(in) :: &
         qgrs           ! Tracer concentrations (kg/kg)
    real(kind_phys), dimension(:, :,:),intent(in) :: &
         aerfld         ! Aerosol input concentrations 
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
         htrlw,       & ! Longwave heating-rate       (K/s)
         htrsw          ! Shortwave heating-rate      (K/s)
    type(sfcflw_type), dimension(:), intent(inout) :: &
         sfcflw         ! Longwave fluxes at surface  (W/m2)
    type(sfcfsw_type), dimension(:), intent(inout) :: &
         sfcfsw         ! Shortwave fluxes at surface (W/m2)
    type(topflw_type), dimension(:), intent(inout) :: &
         topflw         ! Longwave fluxes at TOA      (W/m2)
    type(topfsw_type), dimension(:), intent(inout) :: &
         topfsw         ! Shortwave fluxes at TOA     (W/m2)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg         ! CCPP error message
    integer, intent(out) :: &
         errflg         ! CCPP error flag

    ! Locals
    logical :: top_at_1
    integer :: ipred, ilev, icase, iinf, iLay, iCol, iDay, ncol_pred, iBnd, nbpts, ibpt, itrac
    integer, dimension(nCol) :: ipseed
    real(kind_phys) :: es, qs, dp, tem1, tem2, pfac, ranku(nCol), rankn(nCol), rankk(1), &
         bin(1), edge(nLev+1), bot, top
    real(kind_phys), dimension(nLev) :: rho, tempVar
    real(kind_phys), dimension(nLev+1) :: hgtb, zo, zi
    real(kind_phys), dimension(nCol, nLev) :: o3_lay, tv, rh, aerodp, tau_aero, g_aero, ssa_aero
    real(kind_phys), dimension(nCol, nLev, NF_VGAS) :: gas_vmr
    real(kind_phys), dimension(nCol, nLev, 14, 3) :: aerosolssw
    real(kind_phys), dimension(nCol, nLev, 16, 3) :: aerosolslw
    real(kind_phys), dimension(nCol, 3) :: ext550
    real(c_float), allocatable :: it2f(:,:,:) ! data for inference in profile, height, value order
    real(c_float), allocatable :: ot2f(:,:)   ! data from inference in profile, height order
    real(c_float), dimension(nCol, nLev, size(pnames_lw)) :: predictor_matrix_lw, upredictor_matrix_lw,predictor_matrix_lw2
    real(c_float), dimension(nDay, nLev, size(pnames_sw)) :: predictor_matrix_sw, upredictor_matrix_sw
    real(c_float), dimension(nCol, nLev+2) :: target_matrix_lw
    real(c_float), dimension(nDay, nLev+2) :: target_matrix_sw
    real(kind_phys), dimension(nCol, nLev+1) :: plev, zlev
    real(kind_phys), dimension(nCol, nLev)   :: play, dprs, dz
    real(kind=kind_phys), dimension(nCol,nLev,2:ntrac)  :: tracer1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Are we supposed to be here?
    if (.not. do_mlrad) return

    ! What is vertical ordering of the host?
    top_at_1 = (prsi(1,1) .lt.  prsi(1, nLev))

    ! Change random number seed value for each radiation invocation (isubc = 1 or 2).
    if(isubc == 1) then      ! advance prescribed permutation seed
       do iCol=1,nCol
          ipseed(iCol) = 128 + iCol - 1
       enddo
    elseif (isubc == 2) then ! use input array of permutaion seeds
       do iCol=1,nCol
          ipseed(iCol) = icseed(iCol)
       enddo
    endif

    ! #######################################################################################
    ! Set up trace gas concentrations
    ! #######################################################################################

    ! Do we need to get climatological ozone? (only if not using prognostic ozone)
!    if (i_ozone .le. 0) then
!       call getozn (prslk, lat, nCol, nLev, top_at_1, o3_lay)
!    endif

    ! Get trace-gas concentrations.
    call getgases (prsi/100., lon, lat, nCol, nLev, ico2, top_at_1, con_pi, gas_vmr)

    do itrac = 2, ntrac
       do iLev = 1, nLev
          tracer1(:,iLev,itrac) = max(0.0, qgrs(:,iLev,itrac))
       enddo
    enddo
    
    ! #######################################################################################
    !
    ! Compute prediction matrices for longwave and shortwave
    !
    ! #######################################################################################

    !
    ! Compute Longwave predictor matrix...
    !
    tem1 = 1._kind_phys/con_g
    tem2 = con_rd/con_g
    do iCol=1,nCol
       ! Pressure business (DJS: Will revisit after retraining)
       predictor_matrix_lw(iCol,:,ilw_p)  = prsi(iCol,1:nLev) - 1._kind_phys
       do iLay=1,nLev
          predictor_matrix_lw(iCol,iLay,ilw_dp) = predictor_matrix_lw(iCol,iLay,ilw_p) - predictor_matrix_lw(iCol,iLay+1,ilw_p)
       enddo
       edge(1)      = predictor_matrix_lw(iCol,1,ilw_p)        - predictor_matrix_lw(iCol,1,ilw_dp)/2.0
       edge(2:nLev) = predictor_matrix_lw(iCol,1:nLev-1,ilw_p) + predictor_matrix_lw(iCol,1:nLev-1,ilw_dp)/2.0
       edge(nLev+1) = predictor_matrix_lw(iCol,nLev,ilw_p)     + predictor_matrix_lw(iCol,nLev,ilw_dp)/2.0
       predictor_matrix_lw(iCol,:,ilw_dp) = abs(edge(2:nlev+1)-edge(1:nlev))

       do iLay=1,nLev
          ! (DJS: Will revisit after retraining)
          ! Layer pressure (Pa)
!          predictor_matrix_lw(iCol,iLay,ilw_p)  = prsl(iCol,iLay)

          ! Pressure thickness
!          predictor_matrix_lw(iCol,iLay,ilw_dp) = abs((prsi(iCol,iLay+1))-(prsi(iCol,iLay)))

          ! Thermodynamic temporaries
          es            = min( predictor_matrix_lw(iCol,iLay,ilw_p),  fpvs( tgrs(iCol,iLay) ) )
          qs            = max( con_epsqs, con_eps * es / (predictor_matrix_lw(iCol,iLay,ilw_p) + con_epsm1*es) )
          tv(iCol,iLay) = tgrs(iCol,iLay) * (1._kind_phys + con_fvirt*qgrs(iCol,iLay,1))
          rho(iLay)     = predictor_matrix_lw(iCol,iLay,ilw_p)/(con_rd*tv(iCol,iLay))
          rh(iCol,iLay) = max( 0._kind_phys, min( 1._kind_phys, max(con_epsqs, qgrs(iCol,iLay,1))/qs ) )

          ! Layer temperature (K)
          predictor_matrix_lw(iCol,iLay,ilw_t)  = tgrs(iCol,iLay)

          ! Layer specific-humidity (kg/kg)
          predictor_matrix_lw(iCol,iLay,ilw_q)  = qgrs(iCol,iLay,1)

          ! Compute layer relative-humidity (kg/kg)->(1)
          predictor_matrix_lw(iCol,iLay,ilw_rh)  = rh(iCol,iLay)

          ! Compute layer liquid/Ice water content (kg/kg)->(g/m3).
          ! tracer1 in           (kg/kg)  *
          ! rho in               (kg/m3)  *
          ! conversion-factor    (1e3g/kg)*
          predictor_matrix_lw(iCol,iLay,ilw_lwc)  = max(0._kind_phys, tracer1(iCol,iLay,i_cldliq)*rho(iLay)*1.e3)
          predictor_matrix_lw(iCol,iLay,ilw_iwc)  = max(0._kind_phys, tracer1(iCol,iLay,i_cldice)*rho(iLay)*1.e3)

          ! Compute layer liquid/ice/vapor condensate path, from mixing ratios (kg/kg)->(g/m2).
          ! tracer1 in           (kg/kg)  *
          ! dp in                (kg/ms2) * - Pressure=Force/Area=mass*acceleration/Area=(kg*m*s-2)/(m2) => (kg/ms2)
          ! tem1 in              (s2/m)   * 
          ! conversion-factor in (1e3g/kg)*
          !
          !                      (kg/kg)  * (kg/ms2) * (s2/m) * (1e3g/kg) = > (g/m2)
          !
          ! *NOTE* Multiplying by 1e3 makes to answers worse!!!!!! 
          predictor_matrix_lw(iCol,iLay,ilw_dlwp)  = max(0._kind_phys, tracer1(iCol,iLay,i_cldliq) * &
               tem1 * predictor_matrix_lw(iCol,iLay,ilw_dp))
          predictor_matrix_lw(iCol,iLay,ilw_diwp)  = max(0._kind_phys, tracer1(iCol,iLay,i_cldice) * &
               tem1 * predictor_matrix_lw(iCol,iLay,ilw_dp))
          predictor_matrix_lw(iCol,iLay,ilw_dwvp)  = max(0._kind_phys, qgrs(iCol,iLay,1          ) * &
               tem1 * predictor_matrix_lw(iCol,iLay,ilw_dp))

          ! Cloud effective radii (m).
          if (effr_in) then
             predictor_matrix_lw(iCol,iLay,ilw_reliq) = cld_reliq(iCol,iLay)
             predictor_matrix_lw(iCol,iLay,ilw_reice) = cld_reice(iCol,iLay)
          else
             predictor_matrix_lw(iCol,iLay,ilw_reliq) = reliq_def*1e-6 ! (microns)->(meters)
             predictor_matrix_lw(iCol,iLay,ilw_reice) = reice_def*1e-6 ! (microns)->(meters)
          endif

          ! Ozone mixing-ratio (kg/kg).
          if (i_ozone .gt. 0) then
             predictor_matrix_lw(iCol,iLay,ilw_o3mr) = tracer1(iCol,iLay,i_ozone)
          else
             predictor_matrix_lw(iCol,iLay,ilw_o3mr) = o3_lay(iCol,iLay)
          endif

          ! Trace gases (ppmv)
          predictor_matrix_lw(iCol,iLay,ilw_co2) = gas_vmr(iCol,iLay,1)*1e6
          predictor_matrix_lw(iCol,iLay,ilw_ch4) = gas_vmr(iCol,iLay,3)*1e6
          predictor_matrix_lw(iCol,iLay,ilw_n2o) = gas_vmr(iCol,iLay,2)*1e6
 
       enddo ! END vertical loop
       
       ! Zenith-angle (1) Why do we need this for Longwave?
       predictor_matrix_lw(iCol,:,ilw_sza) = acos(coszen(iCol))

       ! Surface temperature (K)
       predictor_matrix_lw(iCol,:,ilw_sfct) = tgrs(iCol,1)

       ! Surface emissivity (1)
       predictor_matrix_lw(iCol,:,ilw_emiss) = semis(iCol)

       ! Compute vertically integrated (in both directions) liquid/ice/vapor condensate path (kg/m2).
       do iLay=1,nLev
          predictor_matrix_lw(iCol,iLay,ilw_ulwp) = sum(predictor_matrix_lw(iCol,1:iLay,ilw_dlwp))
          predictor_matrix_lw(iCol,iLay,ilw_uiwp) = sum(predictor_matrix_lw(iCol,1:iLay,ilw_diwp))
          predictor_matrix_lw(iCol,iLay,ilw_uwvp) = sum(predictor_matrix_lw(iCol,1:iLay,ilw_dwvp))
       enddo
       do iLay=1,nLev
          predictor_matrix_lw(iCol,iLay,ilw_dlwp)  = sum(predictor_matrix_lw(iCol,iLay:nLev,ilw_dlwp))
          predictor_matrix_lw(iCol,iLay,ilw_diwp)  = sum(predictor_matrix_lw(iCol,iLay:nLev,ilw_diwp))
          predictor_matrix_lw(iCol,iLay,ilw_dwvp)  = sum(predictor_matrix_lw(iCol,iLay:nLev,ilw_dwvp))
       enddo

       ! Layer thickness and height above ground (m) (assumes SFC->TOA vertical ordering)
       ! Layer thickness (m)
       do iLay=nLev,1,-1
          predictor_matrix_lw(iCol,iLay,ilw_dz) = tem2 * abs(log(prsi(iCol,iLay)) - log(prsi(iCol,iLay+1))) * tv(iCol,iLay)
       enddo

       ! Height at layer boundaries
       hgtb(1) = 0._kind_phys
       do iLay=1,nLev
          hgtb(iLay+1)= hgtb(iLay) + predictor_matrix_lw(iCol,iLay,ilw_dz)
          predictor_matrix_lw(iCol,iLay,ilw_z) = hgtb(iLay+1)
       enddo

    enddo    ! END column loop

    !
    ! Compute Shortwave predictor matrix...
    ! (Mostly the same as LW, plus aerosols)
    !
    if (nDay > 0) then

       ! SW Aerosol optics
       ! aerosolssw contains the MERRA aerosol optics properties for the RRTMG SW bands (14).
       call setaer(prsi*0.01, prsl*0.01, prslk, tv, rh, lsmask, tracer1, aerfld, lon, lat, &
            nCol, nLev, nLev+1, .true., .true., iaermdl, iaerflg, top_at_1, con_pi, con_rd,&
            con_g, aerosolssw, aerosolslw, aerodp, ext550, errflg, errmsg)

       ! setaer() provides spectrally defined aerosol optical properties. We need to sum
       ! them to get broadband values.
       tau_aero = 0._kind_phys
       ssa_aero = 0._kind_phys
       g_aero   = 0._kind_phys
       do iCol=1,nCol
          do iLay=1,nLev
             do iBnd=1,14
                tau_aero(iCol,iLay) = tau_aero(iCol,iLay) + aerosolssw(iCol,iLay,iBnd,1)
                ssa_aero(iCol,iLay) = ssa_aero(iCol,iLay) + aerosolssw(iCol,iLay,iBnd,1)* &
                                                            aerosolssw(iCol,iLay,iBnd,2)
                g_aero(iCol,iLay)   = g_aero(iCol,iLay)   + aerosolssw(iCol,iLay,iBnd,1)* &
                                                            aerosolssw(iCol,iLay,iBnd,2)* &
                                                            aerosolssw(iCol,iLay,iBnd,3)
             enddo
             ssa_aero(iCol,iLay) = ssa_aero(iCol,iLay) / max(con_eps, tau_aero(iCol,iLay))
             g_aero(iCol,iLay)   = g_aero(iCol,iLay)   / max(con_eps, ssa_aero(iCol,iLay)*tau_aero(iCol,iLay))
          enddo
          ! Compute vertical extinction.
          do iLay=1,nLev
             tau_aero(iCol,iLay) = sum(tau_aero(iCol,iLay:nLev))
          enddo
       enddo

       do iDay=1,nDay
          predictor_matrix_sw(iDay,iLay,isw_sza)    = predictor_matrix_lw(idx(iDay),iLay,ilw_sza)
          predictor_matrix_sw(iDay,iLay,isw_alb)    = sum(sfcalb(idx(iDay),:))/14.
          predictor_matrix_sw(iDay,iLay,isw_aerssa) = ssa_aero(idx(iDay),iLay)
          predictor_matrix_sw(iDay,iLay,isw_aerasy) = g_aero(idx(iDay),iLay)
          predictor_matrix_sw(iDay,iLay,isw_p)      = predictor_matrix_lw(idx(iDay),iLay,ilw_p)
          predictor_matrix_sw(iDay,iLay,isw_t)      = predictor_matrix_lw(idx(iDay),iLay,ilw_t)
          predictor_matrix_sw(iDay,iLay,isw_q)      = predictor_matrix_lw(idx(iDay),iLay,ilw_q)
          predictor_matrix_sw(iDay,iLay,isw_rh)     = predictor_matrix_lw(idx(iDay),iLay,ilw_rh)
          predictor_matrix_sw(iDay,iLay,isw_lwc)    = predictor_matrix_lw(idx(iDay),iLay,ilw_lwc)
          predictor_matrix_sw(iDay,iLay,isw_iwc)    = predictor_matrix_lw(idx(iDay),iLay,ilw_iwc)
          predictor_matrix_sw(iDay,iLay,isw_dlwp)   = predictor_matrix_lw(idx(iDay),iLay,ilw_dlwp)
          predictor_matrix_sw(iDay,iLay,isw_diwp)   = predictor_matrix_lw(idx(iDay),iLay,ilw_diwp)
          predictor_matrix_sw(iDay,iLay,isw_dwvp)   = predictor_matrix_lw(idx(iDay),iLay,ilw_dwvp)
          predictor_matrix_sw(iDay,iLay,isw_ulwp)   = predictor_matrix_lw(idx(iDay),iLay,ilw_ulwp)
          predictor_matrix_sw(iDay,iLay,isw_uiwp)   = predictor_matrix_lw(idx(iDay),iLay,ilw_uiwp)
          predictor_matrix_sw(iDay,iLay,isw_uwvp)   = predictor_matrix_lw(idx(iDay),iLay,ilw_uwvp)
          predictor_matrix_sw(iDay,iLay,isw_reliq)  = predictor_matrix_lw(idx(iDay),iLay,ilw_reliq)
          predictor_matrix_sw(iDay,iLay,isw_reice)  = predictor_matrix_lw(idx(iDay),iLay,ilw_reice)
          predictor_matrix_sw(iDay,iLay,isw_o3mr)   = predictor_matrix_lw(idx(iDay),iLay,ilw_o3mr)
          predictor_matrix_sw(iDay,iLay,isw_co2)    = predictor_matrix_lw(idx(iDay),iLay,ilw_co2)
          predictor_matrix_sw(iDay,iLay,isw_ch4)    = predictor_matrix_lw(idx(iDay),iLay,ilw_ch4)
          predictor_matrix_sw(iDay,iLay,isw_n2o)    = predictor_matrix_lw(idx(iDay),iLay,ilw_n2o)
          predictor_matrix_sw(iDay,iLay,isw_tauaer) = tau_aero(idx(iDay),iLay)
          predictor_matrix_sw(iDay,iLay,isw_z)      = predictor_matrix_lw(idx(iDay),iLay,ilw_z)
          predictor_matrix_sw(iDay,iLay,isw_dz)     = predictor_matrix_lw(idx(iDay),iLay,ilw_dz)
          predictor_matrix_sw(iDay,iLay,isw_dp)     = predictor_matrix_lw(idx(iDay),iLay,ilw_dp)
       enddo       ! END DO daylit columns
    endif          ! END IF daylit columns

    ! #######################################################################################
    !
    ! Normalize...
    !
    ! #######################################################################################
    
    !
    ! Longwave
    ! 
    if (debug .and. do_debug_once) upredictor_matrix_lw = predictor_matrix_lw
    do iPred = 1,size(pnames_lw)
       do iLay=1,nLev
          !
          ! Determine number of valid breakpoints for current predictor
          !
          if (is2D_lw(iPred)) then
             nbpts = size(mlrad_data%lw%vector_breakpoint(:, iLay, ip2io_lw(iPred)))
             if (any(isnan(mlrad_data%lw%vector_breakpoint(:, iLay, ip2io_lw(iPred))))) then
                do ibpt=1,nbpts
                   if (isnan(mlrad_data%lw%vector_breakpoint(ibpt, iLay, ip2io_lw(iPred)))) then
                      nbpts = ibpt-1
                      exit
                   endif
                enddo
             endif
          else
             nbpts = size(mlrad_data%lw%scalar_breakpoint(:, ip2io_lw(iPred)))
             if (any(isnan(mlrad_data%lw%scalar_breakpoint(:, ip2io_lw(iPred))))) then
                do ibpt=1,nbpts
                   if (isnan(mlrad_data%lw%scalar_breakpoint(ibpt, ip2io_lw(iPred)))) then
                      nbpts = ibpt-1
                      exit
                   endif
                enddo
             endif
          endif
          
          do iCol=1,nCol
             ! Compute the rank of the predictor variable with respect to the training data.
             ! *NOTE* The input data isn't in the order that the emulator requires. Here we 
             ! map the input data, mlrad_data%lw, using ip2io_lw determined during init.
             if (is2D_lw(iPred)) then
                if     (predictor_matrix_lw(iCol,iLay,iPred) <= mlrad_data%lw%vector_breakpoint(  1, iLay, ip2io_lw(iPred))) then
                   rankk(1) = 1
                   ranku(iCol) = eps
                elseif (predictor_matrix_lw(iCol,iLay,iPred) >= mlrad_data%lw%vector_breakpoint(nbpts, iLay, ip2io_lw(iPred))) then
                   rankk(1) = nbpts
                   ranku(iCol) = epsm1
                else
                   rankk = max(minloc(abs(predictor_matrix_lw(iCol,iLay,iPred) - mlrad_data%lw%vector_breakpoint(:,iLay, ip2io_lw(iPred)))) - 1, 1)
                   ranku(iCol) = mlrad_data%lw%vector_intercept(rankk(1), iLay, ip2io_lw(iPred)) + &
                        predictor_matrix_lw(iCol,iLay,iPred)*mlrad_data%lw%vector_slope(rankk(1), iLay, ip2io_lw(iPred))
                endif
                if (debug .and. do_debug_once) then
                   write(99,'(a40)') '------------------------------------------------------------------------------------------------V'
                   write(99,'(a20,i5)' ) trim(pnames_lw(iPred)),nbpts
                   write(99,'(a14,e14.6)')     "        val = ",predictor_matrix_lw(iCol,iLay,iPred)
                   write(99,'(a14,f14.6)')     "      rankk = ",rankk
                   write(99,'(a14,f14.6)')     "      ranku = ",ranku(iCol)
                   write(99,'(a14,i14  )')     "       iLay = ",iLay
                   write(99,'(a14,2e14.6)')    "breakpoints = ",mlrad_data%lw%vector_breakpoint(1, iLay,  ip2io_lw(iPred)),mlrad_data%lw%vector_breakpoint(nbpts, iLay,  ip2io_lw(iPred))
                   write(99,*)mlrad_data%lw%vector_breakpoint(:, iLay,  ip2io_lw(iPred))
                endif
             else
                if     (predictor_matrix_lw(iCol,iLay,iPred) <= mlrad_data%lw%scalar_breakpoint(  1, ip2io_lw(iPred))) then
                   rankk(1) = 1
                   ranku(iCol) = eps
                elseif (predictor_matrix_lw(iCol,iLay,iPred) >= mlrad_data%lw%scalar_breakpoint(nbpts, ip2io_lw(iPred))) then
                   rankk(1) = nbpts
                   ranku(iCol) = epsm1
                else
                   rankk = max(minloc(abs(predictor_matrix_lw(iCol,iLay,iPred) - mlrad_data%lw%scalar_breakpoint(:, ip2io_lw(iPred)))) - 1, 1)
                   ranku(iCol) = mlrad_data%lw%scalar_intercept(rankk(1), ip2io_lw(iPred)) + &
                        predictor_matrix_lw(iCol,iLay,iPred)*mlrad_data%lw%scalar_slope(rankk(1), ip2io_lw(iPred))
                endif
                if (debug .and. do_debug_once) then
                   write(99,'(a40)') '------------------------------------------------------------------------------------------------S'
                   write(99,'(a20)' ) trim(pnames_lw(iPred))
                   write(99,'(a14,e14.6)')     "        val = ",predictor_matrix_lw(iCol,iLay,iPred)
                   write(99,'(a14,f14.6)')     "      rankk = ",rankk
                   write(99,'(a14,f14.6)')     "      ranku = ",ranku(iCol)
                   write(99,'(a14,i14  )')     "       iLay = ",iLay
                   write(99,'(a14,2e14.6)')    "breakpoints = ",mlrad_data%lw%scalar_breakpoint(1, ip2io_lw(iPred)),mlrad_data%lw%scalar_breakpoint(nbpts,  ip2io_lw(iPred))
                endif
             endif
             
             ! Normalize.
             if (do_norm) then
                predictor_matrix_lw(iCol,iLay,iPred) = sqrt(2._kind_phys)*erfinv(2.*ranku(iCol)-1._kind_phys)
             endif

             if (debug .and. do_debug_once) then
                write(99,'(a14,f14.6)')     "       nval = ",predictor_matrix_lw(iCol,iLay,iPred)
             endif
          enddo
       enddo
    enddo

    !
    ! Shortwave
    !
    if (nDay > 0) then
       if (debug .and. do_debug_once) upredictor_matrix_sw = predictor_matrix_sw
       do iPred = 1,size(pnames_sw)
          do iLay=1,nLev
             !
             ! Determine number of valid breakpoints for current predictor
             !
             if (is2D_sw(iPred)) then
                nbpts = size(mlrad_data%sw%vector_breakpoint(:, iLay, ip2io_sw(iPred)))
                if (any(isnan(mlrad_data%sw%vector_breakpoint(:, iLay, ip2io_sw(iPred))))) then
                   do ibpt=1,nbpts
                      if (isnan(mlrad_data%sw%vector_breakpoint(ibpt, iLay, ip2io_sw(iPred)))) then
                         nbpts = ibpt-1
                         exit
                      endif
                   enddo
                endif
             else
                nbpts = size(mlrad_data%sw%scalar_breakpoint(:, ip2io_sw(iPred)))
                if (any(isnan(mlrad_data%sw%scalar_breakpoint(:, ip2io_sw(iPred))))) then
                   do ibpt=1,nbpts
                      if (isnan(mlrad_data%sw%scalar_breakpoint(ibpt, ip2io_sw(iPred)))) then
                         nbpts = ibpt-1
                         exit
                      endif
                   enddo
                endif
             endif

             do iCol=1,nDay
                if (is2D_sw(iPred)) then
                   if     (predictor_matrix_sw(iCol,iLay,iPred) <= mlrad_data%sw%vector_breakpoint(  1, iLay, ip2io_sw(iPred))) then
                      rankk(1) = 1
                      ranku(iCol) = eps
                   elseif (predictor_matrix_sw(iCol,iLay,iPred) >= mlrad_data%sw%vector_breakpoint(nbpts, iLay, ip2io_sw(iPred))) then
                      rankk(1) = nbpts
                      ranku(iCol) = epsm1
                   else
                      rankk = max(minloc(abs(predictor_matrix_sw(iCol,iLay,iPred) - mlrad_data%sw%vector_breakpoint(:,iLay, ip2io_sw(iPred)))) - 1, 1)
                      ranku(iCol) = mlrad_data%sw%vector_intercept(rankk(1), iLay, ip2io_sw(iPred)) + &
                           predictor_matrix_sw(iCol,iLay,iPred)*mlrad_data%sw%vector_slope(rankk(1), iLay, ip2io_sw(iPred))
                   endif
                   if (debug .and. do_debug_once) then
                      write(99,'(a40)') '------------------------------------------------------------------------------------------------V'
                      write(99,'(a20,i5)' ) trim(pnames_sw(iPred)),nbpts
                      write(99,'(a14,e14.6)')     "        val = ",predictor_matrix_sw(iCol,iLay,iPred)
                      write(99,'(a14,f14.6)')     "      rankk = ",rankk
                      write(99,'(a14,f14.6)')     "      ranku = ",ranku(iCol)
                      write(99,'(a14,i14  )')     "       iLay = ",iLay
                      write(99,'(a14,2e14.6)')    "breakpoints = ",mlrad_data%sw%vector_breakpoint(1, iLay,  ip2io_sw(iPred)),mlrad_data%sw%vector_breakpoint(nbpts, iLay,  ip2io_sw(iPred))
                      write(99,*)mlrad_data%sw%vector_breakpoint(:, iLay,  ip2io_sw(iPred))
                   endif
                else
                   if     (predictor_matrix_sw(iCol,iLay,iPred) <= mlrad_data%sw%scalar_breakpoint(  1, ip2io_sw(iPred))) then
                      rankk(1) = 1
                      ranku(iCol) = eps
                   elseif (predictor_matrix_sw(iCol,iLay,iPred) >= mlrad_data%sw%scalar_breakpoint(nbpts, ip2io_sw(iPred))) then
                      rankk(1) = nbpts
                      ranku(iCol) = epsm1
                   else
                      rankk = max(minloc(abs(predictor_matrix_sw(iCol,iLay,iPred) - mlrad_data%sw%scalar_breakpoint(:, ip2io_sw(iPred)))) - 1, 1)
                      ranku(iCol) = mlrad_data%sw%scalar_intercept(rankk(1), ip2io_sw(iPred)) + &
                           predictor_matrix_sw(iCol,iLay,iPred)*mlrad_data%sw%scalar_slope(rankk(1), ip2io_sw(iPred))
                   endif
                   if (debug .and. do_debug_once) then
                      write(99,'(a40)') '------------------------------------------------------------------------------------------------S'
                      write(99,'(a20)' ) trim(pnames_sw(iPred))
                      write(99,'(a14,e14.6)')     "        val = ",predictor_matrix_sw(iCol,iLay,iPred)
                      write(99,'(a14,f14.6)')     "      rankk = ",rankk
                      write(99,'(a14,f14.6)')     "      ranku = ",ranku(iCol)
                      write(99,'(a14,i14  )')     "       iLay = ",iLay
                      write(99,'(a14,2e14.6)')    "breakpoints = ",mlrad_data%sw%scalar_breakpoint(1, ip2io_sw(iPred)),mlrad_data%sw%scalar_breakpoint(nbpts,  ip2io_sw(iPred))
                   endif
                endif

                ! Normalize
                predictor_matrix_sw(iDay,iLay,iPred) = sqrt(2._kind_phys)*erfinv(2.*ranku(iDay)-1._kind_phys)

             enddo
          enddo
       enddo
    endif

    ! #######################################################################################
    ! Begin Inference...
    ! #######################################################################################
    ! Run inferences
    ! Longwave (all columns)
    call infero_check(model_lw%infer(predictor_matrix_lw, target_matrix_lw))

    ! Shortwave (daylit-columns only)
    if (nDay > 0) then
       call infero_check(model_sw%infer(predictor_matrix_sw, target_matrix_sw))
    endif

    if (debug) then
       do iCol=1,nCol
          status = nf90_put_var(ncid,varID0, fhour,                                start=record_counter)
          status = nf90_put_var(ncid,varID1, sfcflw(iCol)%dnfxc,                   start=record_counter)
          status = nf90_put_var(ncid,varID2, topflw(iCol)%upfxc,                   start=record_counter)
          status = nf90_put_var(ncid,varID3, target_matrix_lw(iCol,nLev+1),        start=record_counter)
          status = nf90_put_var(ncid,varID4, target_matrix_lw(iCol,nLev+2),        start=record_counter)
          status = nf90_put_var(ncid,varID5, predictor_matrix_lw(iCol,:,ilw_lwc),  start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID6, predictor_matrix_lw(iCol,:,ilw_iwc),  start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID7, predictor_matrix_lw(iCol,:,ilw_z),    start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID8, predictor_matrix_lw(iCol,:,ilw_dlwp), start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID9, predictor_matrix_lw(iCol,:,ilw_diwp), start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID10,predictor_matrix_lw(iCol,:,ilw_dwvp), start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID12, rho,                                 start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID13, predictor_matrix_lw(iCol,:,ilw_ulwp),start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID14, predictor_matrix_lw(iCol,:,ilw_uiwp),start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID15, predictor_matrix_lw(iCol,:,ilw_uwvp),start=(/1,record_counter/),count=(/127,1/))

          status = nf90_put_var(ncid,varID11, qgrs(iCol,:,       1),               start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID16, tracer1(iCol,:,i_cldliq),            start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID17, tracer1(iCol,:,i_cldice),            start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID18, tracer1(iCol,:,       3),            start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID19, tracer1(iCol,:,       5),            start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID20, tracer1(iCol,:,       6),            start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID21, tracer1(iCol,:,       9),            start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID22, predictor_matrix_lw(iCol,:,ilw_dp),  start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID23, predictor_matrix_lw(iCol,:,ilw_dz),  start=(/1,record_counter/),count=(/127,1/))
          status = nf90_put_var(ncid,varID24, predictor_matrix_lw(iCol,:,ilw_p),   start=(/1,record_counter/),count=(/127,1/))
       enddo
       record_counter(1) = record_counter(1) + 1
       do iCol=1,nCol
          !write(94,'(a20,4a18)'  ) '                  ','ML(LW down @ SFC)','G(LW down @ SFC)','ML(LW up @ TOA)', 'G(LW up @ TOA)'
          !write(94,'(a20,3f12.2)') 'LW surface(down): ',target_matrix_lw(iCol,nLev+1),sfcflw(iCol)%dnfxc,sfcflw(iCol)%dnfx0
          !write(94,'(a20,3f12.2)') 'LW toa(up):       ',target_matrix_lw(iCol,nLev+2),topflw(iCol)%upfxc,topflw(iCol)%upfx0
          write(94,'(4f18.2)') target_matrix_lw(iCol,nLev+1),sfcflw(iCol)%dnfxc,target_matrix_lw(iCol,nLev+2),topflw(iCol)%upfxc
          !write(94,'(a20,3f12.2)') 'SW surface(down): ',target_matrix_sw(iCol,nLev+1),sfcfsw(iCol)%dnfxc,sfcfsw(iCol)%dnfx0
          !write(94,'(a20,3f12.2)') 'SW toa(up):       ',target_matrix_sw(iCol,nLev+2),topfsw(iCol)%upfxc,topfsw(iCol)%upfx0
       enddo

       do iCol=1,nCol
          !write(93,'(a30)') 'heating profile (K/day): '
          !write(93,'(3a12)'  ) 'Layer','ML(all)','G(all)'
          !write(93,'(3a12)'  ) '     ','(LW)',   '(LW)'
          !do iLay=1,nLev
          !   write(93,'(a9,i3,2f12.2)') '',iLay,target_matrix_lw(iCol,iLay),htrlw(iCol,iLay)*3600.*24.
          !enddo
          write(93,'(127f12.2)') target_matrix_lw(iCol,:)
          write(93,'(127f12.2)') htrlw(iCol,:)*3600.*24.
       enddo
    end if
    
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! REMOVE WHEN WORKING
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    if (debug .and. do_debug_once) then
       ! ####################################################################################
       ! Compute fluxes from reference data (sanity check for proper infero implementation)
       ! ####################################################################################
       do iCol=1,nCol
          do iLay=1,nLev
             predictor_matrix_lw2(iCol,iLay,:) = ref_data%predictor_matrix(:,iLay,iCol)
          enddo
       enddo
       call infero_check(model_lw%infer(predictor_matrix_lw2, target_matrix_lw))

       ! ####################################################################################
       ! Raw predictor file (97)
       ! ####################################################################################
       write(97, '(a5, 24a10)')  'Layer', 'sza', 'sfcT', 'sfc_emiss', 'p', 'T', 'q', 'rh',  &
            'LWC', 'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                      &
            'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                             &
            'dp'
       write(97, '(a5, 24a10)')  '', '(1)', '(K)', '(1)','(Pa)', '(k)', '(kg/kg)', '(1)',   &
            '(kg/m3)', '(kg/m3)','(kg/m2)', '(kg/m2)', '(kg/m2)', '(kg/m2)', '(kg/m2)', '(kg/m2)',&
            '(liq)','(ice)', '(mg/kg)','(ppmv)', '(ppmv)', '(ppmv)', '(m)', '(m)',          &
            '(Pa)'
       do iCol=1,nCol
          do iLay=1,nLev
             write (97,'(i5,7f10.2,4e10.2,f10.2,2e10.2,10f10.4)')&
                  iLay,                                       &
                  upredictor_matrix_lw(iCol,iLay,ilw_sza),    &
                  upredictor_matrix_lw(iCol,iLay,ilw_sfct),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_emiss),  &
                  upredictor_matrix_lw(iCol,iLay,ilw_p),      &
                  upredictor_matrix_lw(iCol,iLay,ilw_t),      &
                  upredictor_matrix_lw(iCol,iLay,ilw_q),      &
                  upredictor_matrix_lw(iCol,iLay,ilw_rh),     &
                  upredictor_matrix_lw(iCol,iLay,ilw_lwc),    &
                  upredictor_matrix_lw(iCol,iLay,ilw_iwc),    &
                  upredictor_matrix_lw(iCol,iLay,ilw_dlwp),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_diwp),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_dwvp),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_ulwp),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_uiwp),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_uwvp),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_reliq),  &
                  upredictor_matrix_lw(iCol,iLay,ilw_reice),  &
                  1e9*upredictor_matrix_lw(iCol,iLay,ilw_o3mr),   &
                  upredictor_matrix_lw(iCol,iLay,ilw_co2),    &
                  upredictor_matrix_lw(iCol,iLay,ilw_ch4),    &
                  upredictor_matrix_lw(iCol,iLay,ilw_n2o),    &
                  upredictor_matrix_lw(iCol,iLay,ilw_z),      &
                  upredictor_matrix_lw(iCol,iLay,ilw_dz),     &
                  upredictor_matrix_lw(iCol,iLay,ilw_dp)
          enddo
       enddo

       ! ####################################################################################
       ! Normalizaed predictor file (98)
       ! ####################################################################################
       write(98, '(a5, 24a10)')  'Layer', 'sza', 'sfcT', 'sfc_emiss', 'p', 'T', 'q', 'rh', 'LWC',&
            'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                             &
            'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                             &
            'dp'

       do iCol=1,nCol
          do iLay=1,nLev
             write (98,'(i5,24f10.2)')                       &
                  iLay,                                      &
                  predictor_matrix_lw(iCol,iLay,ilw_sza),    &
                  predictor_matrix_lw(iCol,iLay,ilw_sfct),   &
                  predictor_matrix_lw(iCol,iLay,ilw_emiss),  &
                  predictor_matrix_lw(iCol,iLay,ilw_p),      &
                  predictor_matrix_lw(iCol,iLay,ilw_t),      &
                  predictor_matrix_lw(iCol,iLay,ilw_q),      &
                  predictor_matrix_lw(iCol,iLay,ilw_rh),     &
                  predictor_matrix_lw(iCol,iLay,ilw_lwc),    &
                  predictor_matrix_lw(iCol,iLay,ilw_iwc),    &
                  predictor_matrix_lw(iCol,iLay,ilw_dlwp),   &
                  predictor_matrix_lw(iCol,iLay,ilw_diwp),   &
                  predictor_matrix_lw(iCol,iLay,ilw_dwvp),   &
                  predictor_matrix_lw(iCol,iLay,ilw_ulwp),   &
                  predictor_matrix_lw(iCol,iLay,ilw_uiwp),   &
                  predictor_matrix_lw(iCol,iLay,ilw_uwvp),   &
                  predictor_matrix_lw(iCol,iLay,ilw_reliq),  &
                  predictor_matrix_lw(iCol,iLay,ilw_reice),  &
                  predictor_matrix_lw(iCol,iLay,ilw_o3mr),   &
                  predictor_matrix_lw(iCol,iLay,ilw_co2),    &
                  predictor_matrix_lw(iCol,iLay,ilw_ch4),    &
                  predictor_matrix_lw(iCol,iLay,ilw_n2o),    &
                  predictor_matrix_lw(iCol,iLay,ilw_z),      &
                  predictor_matrix_lw(iCol,iLay,ilw_dz),     &
                  predictor_matrix_lw(iCol,iLay,ilw_dp)
          enddo
       enddo

       ! ####################################################################################
       ! Raw predictor example (95)
       ! ####################################################################################
       write(95, '(a5, 25a10)')  'Layer', 'p', 'T', 'q', 'rh', 'LWC',                       &
            'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                             &
            'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                             &
            'dp', 'sza', 'Tsfc', 'sfc_emiss'
       write(95, '(a5, 25a10)')  '', '(Pa)', '(k)', '(g/kg)', '(1)', '(g/m3)',              &
            '(g/m3)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)', '(g/m2)',           &
            '(liq)','(ice)', '(mg/kg)','(ppmv)', '(ppmv)', '(ppmv)', '(m)', '(m)',          &
            '(Pa)', '(1)', '(K)', '(1)'
       do iCol=1,nCol
          do iLay=1,nLev
             write (95,'(i5,4f10.2,4e10.2,f10.2,2e10.2,17f10.4)')      &
                  iLay,                                                &
                  ref_data%unnorm_predictor_matrix(1,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(2,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(3,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(4,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(5,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(6,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(7,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(8,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(9,iLay,iCol),       &
                  ref_data%unnorm_predictor_matrix(10,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(11,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(12,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(13,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(14,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(15,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(16,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(17,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(18,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(19,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(20,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(21,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(22,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(23,iLay,iCol),      &
                  ref_data%unnorm_predictor_matrix(24,iLay,iCol)
          enddo
       enddo

       ! ####################################################################################
       ! Normalized predictor examples (96)
       ! ####################################################################################
       write(96, '(a5, 25a10)')  'Layer', 'p', 'T', 'q', 'rh', 'LWC',                       &
            'IWC', 'LWP', 'IWP', 'WVP', 'iLWP', 'iIWP', 'iWVP',                             &
            'Reff','Reff', 'o3','co2', 'ch4', 'n2o', 'z', 'dz',                             &
            'dp', 'sza', 'Tsfc', 'sfc_emiss'
       do iCol=1,nCol
          do iLay=1,nLev
             write (96,'(i5,25f10.2)')                     &
                  iLay,                                    &
                  ref_data%predictor_matrix(1,iLay,iCol),  &
                  ref_data%predictor_matrix(2,iLay,iCol),  &
                  ref_data%predictor_matrix(3,iLay,iCol),  &
                  ref_data%predictor_matrix(4,iLay,iCol),  &
                  ref_data%predictor_matrix(5,iLay,iCol),  &
                  ref_data%predictor_matrix(6,iLay,iCol),  &
                  ref_data%predictor_matrix(7,iLay,iCol),  &
                  ref_data%predictor_matrix(8,iLay,iCol),  &
                  ref_data%predictor_matrix(9,iLay,iCol),  &
                  ref_data%predictor_matrix(10,iLay,iCol), &
                  ref_data%predictor_matrix(11,iLay,iCol), &
                  ref_data%predictor_matrix(12,iLay,iCol), &
                  ref_data%predictor_matrix(13,iLay,iCol), &
                  ref_data%predictor_matrix(14,iLay,iCol), &
                  ref_data%predictor_matrix(15,iLay,iCol), &
                  ref_data%predictor_matrix(16,iLay,iCol), &
                  ref_data%predictor_matrix(17,iLay,iCol), &
                  ref_data%predictor_matrix(18,iLay,iCol), &
                  ref_data%predictor_matrix(19,iLay,iCol), &
                  ref_data%predictor_matrix(20,iLay,iCol), &
                  ref_data%predictor_matrix(21,iLay,iCol), &
                  ref_data%predictor_matrix(22,iLay,iCol), &
                  ref_data%predictor_matrix(23,iLay,iCol), &
                  ref_data%predictor_matrix(24,iLay,iCol)
          enddo
       enddo
    endif
    do_debug_once = .false.
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
    ! END DEBUG BLOCK
    ! #######################################################################################
    ! #######################################################################################
    ! #######################################################################################
 
    ! Copy from prediction-matrix to ccpp interstitials (for prognostic ML rad)
    do iCol=1,nCol
!       htrlw(iCol,1:nLev) = target_matrix_lw(iCol,1:nLev)/(3600.*24.) ! K/day -> K/sec
!       sfcflw(iCol)%dnfxc = target_matrix_lw(iCol,nLev+1)
!       sfcflw(iCol)%upfxc = 0.
!       topflw(iCol)%upfxc = target_matrix_lw(iCol,nLev+2)
    enddo
    do iDay=1,nDay
!       htrsw(idx(iDay),1:nLev) = target_matrix_sw(iDay,1:nLev)/(3600.*24.) ! K/day -> K/sec
!       sfcfsw(idx(iDay))%dnfxc = target_matrix_sw(iDay,nLev+1)
!       sfcfsw(idx(iDay))%upfxc = 0.
!       topfsw(idx(iDay))%upfxc = target_matrix_sw(iDay,nLev+2)
!       topfsw(idx(iDay))%dnfxc = 0.
    enddo

  end subroutine mlrad_driver_run

! ###########################################################################################
!! \section arg_table_mlrad_driver_finalize
!! \htmlinclude mlrad_driver_finalize.html
!!
! ###########################################################################################
  subroutine mlrad_driver_finalize(do_mlrad, debug, errmsg, errflg)
    ! Inputs
    logical,           intent(in) :: do_mlrad, debug
    ! Outputs
    character(len=*),  intent(out) :: errmsg
    integer,           intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_mlrad) return

    ! Free the model(s)
    call infero_check(model_lw%free())
    call infero_check(model_sw%free())

    ! Finalize
    call infero_check(infero_finalise())

    status = nf90_close(ncid)
    if (debug) then
       close(93)
       close(94)
       close(95)
       close(96)
       close(97)
       close(98)
       close(99)
    endif

  end subroutine mlrad_driver_finalize

  ! #########################################################################################
  ! ERFINV: Using Newton refinement method
  ! (http://www.mimirgames.com/articles/programming/approximations-of-the-inverse-error-function)
  ! #########################################################################################
  function erfinv(var_in) result(var_out)
    real(kind_phys),intent(in) :: var_in
    real(kind_phys) :: var_out
    real(kind_dbl_prec),parameter,dimension(4) :: &
         a = (/ 0.886226899, -1.645349621,  0.914624893, -0.140543331/), &
         b = (/-2.118377725,  1.442710462, -0.329097515,  0.012229801/),&
         c = (/-1.970840454, -1.624906493,  3.429567803,  1.641345311/)
    real(kind_dbl_prec),parameter,dimension(2) :: &
         d = (/ 3.543889200,  1.637067800/)
    real(kind_dbl_prec) :: x, z

    if (abs(var_in) == 1._kind_dbl_prec) then
       x = var_in
    elseif(var_in < -0.7_kind_dbl_prec) then
       z = sqrt(-log((1._kind_dbl_prec+var_in)/2._kind_dbl_prec))
       x = -(((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0)
    else
       if (var_in > 0.7) then
          z = var_in*var_in
          x = var_in*(((a(4)*z+a(3))*z+a(2))*z+a(1))/((((b(4)*z+b(4))*z+b(2))*z+b(1))*z+1.0)
       else
          z = sqrt(-Log((1._kind_dbl_prec-var_in)/2._kind_dbl_prec));
          x = (((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0)
       endif
    endif
    var_out = x - (erf(x) - var_in) / (2._kind_dbl_prec/sqrt(3.14_kind_dbl_prec) * exp(-x*x))

  end function erfinv

end module mlrad_driver
