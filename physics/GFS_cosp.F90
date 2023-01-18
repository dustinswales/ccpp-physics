!> \file GFS_cosp.F90
  ! #########################################################################################
module GFS_cosp
  use machine,              only: kind_phys
  use quickbeam,            only: radar_cfg
  use mod_quickbeam_optics, only: size_distribution, hydro_class_init, quickbeam_optics_init
  use mod_cosp,             only: cosp_outputs,cosp_optical_inputs,cosp_column_inputs,      &
       cosp_init
  use mod_cosp_config,      only: pres_binCenters, pres_binEdges, tau_binCenters,           &
       tau_binEdges, cloudsat_binCenters, cloudsat_binEdges, calipso_binCenters,            &
       calipso_binEdges, misr_histHgtCenters, misr_histHgtEdges,  PARASOL_SZA, R_UNDEF,     &
       PARASOL_NREFL, LIDAR_NCAT,SR_BINS, N_HYDRO, RTTOV_MAX_CHANNELS, numMISRHgtBins,      &
       CLOUDSAT_DBZE_BINS, LIDAR_NTEMP, calipso_histBsct, numMODISTauBins, numMODISPresBins,&
       numMODISReffIceBins, numMODISReffLiqBins, numISCCPTauBins, numISCCPPresBins,         &
       numMISRTauBins, reffICE_binEdges, reffICE_binCenters, reffLIQ_binEdges,              &
       reffLIQ_binCenters, LIDAR_NTYPE, nCloudsatPrecipClass, nsza_cosp => PARASOL_NREFL,   &
       nprs_cosp => npres, ntau_cosp => ntau, ntau_cosp_modis => ntau, nsr_cosp => SR_BINS, &
       nhtmisr_cosp => numMISRHgtBins, nhydro => N_HYDRO, cloudsat_preclvl, vgrid_zl,       &
       vgrid_zu, vgrid_z
    use mod_cosp_stats,     only: cosp_change_vertical_grid

  implicit none
  
  ! DO WE NEED THESE COPIES?
  real(kind_phys) :: prsmid_cosp(nprs_cosp)
  real(kind_phys) :: prslim_cosp(2,nprs_cosp)
  real(kind_phys) :: taumid_cosp(ntau_cosp)
  real(kind_phys) :: taulim_cosp(2,ntau_cosp)
  real(kind_phys) :: srmid_cosp(nsr_cosp)
  real(kind_phys) :: srlim_cosp(2,nsr_cosp)
  real(kind_phys) :: sza_cosp(nsza_cosp)
  real(kind_phys) :: dbzemid_cosp(CLOUDSAT_DBZE_BINS)
  real(kind_phys) :: dbzelim_cosp(2,CLOUDSAT_DBZE_BINS)
  real(kind_phys) :: htmisrmid_cosp(nhtmisr_cosp)
  real(kind_phys) :: htmisrlim_cosp(2,nhtmisr_cosp)
  real(kind_phys) :: taumid_cosp_modis(ntau_cosp_modis)
  real(kind_phys) :: taulim_cosp_modis(2,ntau_cosp_modis)
  real(kind_phys) :: reffICE_binEdges_cosp(2,numMODISReffIceBins)
  real(kind_phys) :: reffLIQ_binEdges_cosp(2,numMODISReffLiqBins)
  real(kind_phys) :: reffICE_binCenters_cosp(numMODISReffIceBins)
  real(kind_phys) :: reffLIQ_binCenters_cosp(numMODISReffLiqBins)

  ! Note: Unless otherwise specified, these are parameters that cannot be set by the CAM namelist.
  integer, parameter :: Npoints_it = 10000       ! Max # gridpoints to be processed in one iteration (10,000)
  integer :: ncolumns = 50                       ! Number of subcolumns in SCOPS (50), can be changed from default by CAM namelist
  integer :: nlr = 40                            ! Number of levels in statistical outputs 
                                                 ! (only used if USE_VGRID=.true.)  (40)
  logical :: use_vgrid = .true.                  ! Use fixed vertical grid for outputs? 
                                                 ! (if .true. then define # of levels with nlr)  (.true.)
  logical :: csat_vgrid = .true.                 ! CloudSat vertical grid? 
                                                 ! (if .true. then the CloudSat standard grid is used.
                                                 ! If set, overides use_vgrid.) (.true.)
  ! namelist variables for COSP input related to radar simulator
  real(kind_phys) :: radar_freq = 94.0_kind_phys ! CloudSat radar frequency (GHz) (94.0)
  integer :: surface_radar = 0                   ! surface=1, spaceborne=0 (0)
  integer :: use_mie_tables = 0                  ! use a precomputed lookup table? yes=1,no=0 (0)
  integer :: use_gas_abs = 1                     ! include gaseous absorption? yes=1,no=0 (1)
  integer :: do_ray = 0                          ! calculate/output Rayleigh refl=1, not=0 (0)
  integer :: melt_lay = 0                        ! melting layer model off=0, on=1 (0)
  real(kind_phys) :: k2 = -1                     ! |K|^2, -1=use frequency dependent default (-1)
  ! namelist variables for COSP input related to lidar simulator
  integer, parameter :: Nprmts_max_hydro = 12    ! Max # params for hydrometeor size distributions (12)
  integer, parameter :: Naero = 1                ! Number of aerosol species (Not used) (1)
  integer, parameter :: Nprmts_max_aero = 1      ! Max # params for aerosol size distributions (not used) (1)
  integer :: lidar_ice_type = 0                  ! Ice particle shape in lidar calculations
                                                 ! (0=ice-spheres ; 1=ice-non-spherical) (0)
  integer, parameter :: overlap = 3              ! overlap type: 1=max, 2=rand, 3=max/rand (3)

  !! namelist variables for COSP input related to ISCCP simulator
  integer :: isccp_topheight = 1                 ! 1 = adjust top height using both a computed infrared
                                                 ! brightness temperature and the visible
                                                 ! optical depth to adjust cloud top pressure.
                                                 ! Note that this calculation is most appropriate to compare
                                                 ! to ISCCP data during sunlit hours.
                                                 ! 2 = do not adjust top height, that is cloud top pressure
                                                 ! is the actual cloud top pressure in the model
                                                 ! 3 = adjust top height using only the computed infrared
                                                 ! brightness temperature. Note that this calculation is most
                                                 ! appropriate to compare to ISCCP IR only algortihm (i.e.
                                                 ! you can compare to nighttime ISCCP data with this option) (1)
  integer :: isccp_topheight_direction = 2       ! direction for finding atmosphere pressure level with
                                                 ! interpolated temperature equal to the radiance
                                                 ! determined cloud-top temperature
                                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                                 ! 1 = default setting in COSP v1.1, matches all versions of
                                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                                 ! 2 = default setting in COSP v1.3. default since V4.0 of ISCCP simulator

  type(radar_cfg)                      :: rcfg_cloudsat ! Radar configuration (Cloudsat)
  type(radar_cfg), allocatable         :: rcfg_cs(:)    ! chunked version of rcfg_cloudsat
  type(size_distribution)              :: sd            ! Size distribution used by radar simulator
  type(size_distribution), allocatable :: sd_cs(:)      ! chunked version of sd
  character(len=64)                    :: cloudsat_micro_scheme

contains

  ! #########################################################################################
!! \section arg_table_GFS_cosp_init
!! \htmlinclude GFS_cosp_init.html
!!
  subroutine GFS_cosp_init(do_cosp_isccp, do_cosp_modis, do_cosp_misr, do_cosp_cloudsat,    &
       do_cosp_calipso, do_cosp_grLidar532, do_cosp_atlid, do_cosp_parasol, imp_physics,    &
       imp_physics_thompson, imp_physics_gfdl, errmsg, errflg)
    USE mod_cosp_modis_interface,      ONLY: cosp_modis_init
    USE mod_cosp_misr_interface,       ONLY: cosp_misr_init
    USE mod_cosp_isccp_interface,      ONLY: cosp_isccp_init
    USE mod_cosp_calipso_interface,    ONLY: cosp_calipso_init
    USE mod_cosp_atlid_interface,      ONLY: cosp_atlid_init
    USE mod_cosp_grlidar532_interface, ONLY: cosp_grLidar532_init
    USE mod_cosp_parasol_interface,    ONLY: cosp_parasol_init
    USE mod_cosp_cloudsat_interface,   ONLY: cosp_cloudsat_init

    ! Inputs
    logical, intent(in)    :: &
         do_cosp_isccp,        & !
         do_cosp_modis,        & !
         do_cosp_misr,         & !
         do_cosp_cloudsat,     & !
         do_cosp_calipso,      & !
         do_cosp_grLidar532,   & !
         do_cosp_atlid,        & !
         do_cosp_parasol
    integer, intent(in)    ::  &
         imp_physics,          & ! Choice of microphysics scheme
         imp_physics_thompson, & ! Choice of Thompson
         imp_physics_gfdl        ! Choice of GFDL
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Initialize bin boundaries
    prsmid_cosp             = pres_binCenters
    prslim_cosp             = pres_binEdges
    taumid_cosp             = tau_binCenters
    taulim_cosp             = tau_binEdges
    srmid_cosp              = calipso_binCenters
    srlim_cosp              = calipso_binEdges
    sza_cosp                = parasol_sza
    dbzemid_cosp            = cloudsat_binCenters
    dbzelim_cosp            = cloudsat_binEdges
    htmisrmid_cosp          = misr_histHgtCenters
    htmisrlim_cosp          = misr_histHgtEdges
    taumid_cosp_modis       = tau_binCenters
    taulim_cosp_modis       = tau_binEdges
    reffICE_binCenters_cosp = reffICE_binCenters
    reffICE_binEdges_cosp   = reffICE_binEdges
    reffLIQ_binCenters_cosp = reffLIQ_binCenters
    reffLIQ_binEdges_cosp   = reffLIQ_binEdges

    ! Initialize the distributional parameters for hydrometeors in radar simulator
    if (imp_physics == imp_physics_thompson) then
       call hydro_class_init(.false., .true., sd)
       cloudsat_micro_scheme = 'MMF_v3.5_double_moment'
    endif
    if (imp_physics == imp_physics_gfdl) then
       call hydro_class_init(.true., .false., sd)
       cloudsat_micro_scheme = 'MMF_v3.5_single_moment'
    endif

    !
    if (do_cosp_cloudsat .or. do_cosp_grLidar532) then
       call quickbeam_optics_init()
    endif
    if (do_cosp_isccp) then
       call cosp_isccp_init(isccp_topheight,isccp_topheight_direction)
       endif
    if (do_cosp_modis) then
       call cosp_modis_init()
    endif
    if (do_cosp_misr) then
       call cosp_misr_init()
    endif
    if (do_cosp_cloudsat) then
       call cosp_cloudsat_init(radar_freq, k2, use_gas_abs, do_ray, R_UNDEF, N_HYDRO,      &
            surface_radar, rcfg_cloudsat, cloudsat_micro_scheme)
    endif
    if (do_cosp_calipso) then
       call cosp_calipso_init()
    endif
    if (do_cosp_grLidar532) then
       call cosp_grLidar532_init()
    endif
    if (do_cosp_atlid) then
       call cosp_atlid_init()
    endif
    if (do_cosp_parasol) then
       call cosp_parasol_init()
    endif

  end subroutine GFS_cosp_init

  ! #########################################################################################
!! \section arg_table_GFS_cosp_run
!! \htmlinclude GFS_cosp_run.html
!!
  subroutine GFS_cosp_run(errmsg, errflg)
    ! Inputs
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine GFS_cosp_run

end module GFS_cosp
