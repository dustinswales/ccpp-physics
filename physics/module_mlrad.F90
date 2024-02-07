module module_mlrad
!> \section arg_table_module_mlrad Argument table
!! \htmlinclude module_mlrad.html
!!
  use machine, only : kind_phys
  use netcdf
  implicit none

  ! #####################################################################################
  !
  ! Shortwave (SW)
  !
  ! ##################################################################################### 
  
  ! Number of SW predictors
  integer,parameter :: npred_sw = 26

  ! Predictor varaible names, from training data, in order expected by emulator.
  character(len=32),dimension(npred_sw) :: &
       pnames_sw =  &
       (/'zenith_angle_radians            ','albedo                          ',&
         'aerosol_single_scattering_albed ','aerosol_asymmetry_param         ',&
         'pressure_pascals                ','temperature_kelvins             ',&
         'specific_humidity_kg_kg01       ','relative_humidity_unitless      ',&
         'liquid_water_content_kg_m03     ','ice_water_content_kg_m03        ',&
         'liquid_water_path_kg_m02        ','ice_water_path_kg_m02           ',&
         'vapour_path_kg_m02              ','upward_liquid_water_path_kg_m02 ',&
         'upward_ice_water_path_kg_m02    ','upward_vapour_path_kg_m02       ',&
         'liquid_effective_radius_metres  ','ice_effective_radius_metres     ',&
         'o3_mixing_ratio_kg_kg01         ','co2_concentration_ppmv          ',&
         'ch4_concentration_ppmv          ','n2o_concentration_ppmv          ',&
         'aerosol_extinction_metres01     ','height_m_agl                    ',&
         'height_thickness_metres         ','pressure_thickness_pascals      '/)

  ! Indices into prediction matrix
  integer, parameter :: &
       isw_sza    = 1,  &
       isw_alb    = 2,  &
       isw_aerssa = 3,  &
       isw_aerasy = 4,  &
       isw_p      = 5,  &
       isw_t      = 6,  &
       isw_q      = 7,  &
       isw_rh     = 8,  &
       isw_lwc    = 9,  &
       isw_iwc    = 10, &
       isw_dlwp   = 11, &
       isw_diwp   = 12, &
       isw_dwvp   = 13, &
       isw_ulwp   = 14, &
       isw_uiwp   = 15, &
       isw_uwvp   = 16, &
       isw_reliq  = 17, &
       isw_reice  = 18, &
       isw_o3mr   = 19, &
       isw_co2    = 20, &
       isw_ch4    = 21, &
       isw_n2o    = 22, &
       isw_tauaer = 23, &
       isw_z      = 24, &
       isw_dz     = 25, &
       isw_dp     = 26

  ! ##################################################################################### 
  !
  ! Longwave (LW)
  !
  ! ##################################################################################### 

  ! Number of LW predictors.
  integer,parameter :: npred_lw = 24

  ! Predictor varaible names, from training data, in order expected by emulator.
!  character(len=32),dimension(npred_lw) :: &
!       pnames_lw =  &
!       (/'zenith_angle_radians            ','surface_temperature_kelvins     ',&
!         'surface_emissivity              ',                                   &
!         'pressure_pascals                ','temperature_kelvins             ',&
!         'specific_humidity_kg_kg01       ','relative_humidity_unitless      ',&
!         'liquid_water_content_kg_m03     ','ice_water_content_kg_m03        ',&
!         'liquid_water_path_kg_m02        ','ice_water_path_kg_m02           ',&
!         'vapour_path_kg_m02              ','upward_liquid_water_path_kg_m02 ',&
!         'upward_ice_water_path_kg_m02    ','upward_vapour_path_kg_m02       ',&
!         'liquid_effective_radius_metres  ','ice_effective_radius_metres     ',&
!         'o3_mixing_ratio_kg_kg01         ','co2_concentration_ppmv          ',&
!         'ch4_concentration_ppmv          ','n2o_concentration_ppmv          ',&
!         'height_m_agl                    ','height_thickness_metres         ',&
!         'pressure_thickness_pascals      '/)

  ! Indices into prediction matrix
!  integer, parameter :: &
!       ilw_sza    = 1,  &
!       ilw_sfct   = 2,  &
!       ilw_emiss  = 3,  &
!       ilw_p      = 4,  &
!       ilw_t      = 5,  &
!       ilw_q      = 6,  &
!       ilw_rh     = 7,  &
!       ilw_lwc    = 8,  &
!       ilw_iwc    = 9,  &
!       ilw_dlwp   = 10, &
!       ilw_diwp   = 11, &
!       ilw_dwvp   = 12, &
!       ilw_ulwp   = 13, &
!       ilw_uiwp   = 14, &
!       ilw_uwvp   = 15, &
!       ilw_reliq  = 16, &
!       ilw_reice  = 17, &
!       ilw_o3mr   = 18, &
!       ilw_co2    = 19, &
!       ilw_ch4    = 20, &
!       ilw_n2o    = 21, &
!       ilw_z      = 22, &
!       ilw_dz     = 23, &
!       ilw_dp     = 24

  integer, parameter :: &
       ilw_p      = 1,  &
       ilw_t      = 2,  &
       ilw_q      = 3,  &
       ilw_rh     = 4,  &
       ilw_lwc    = 5,  &
       ilw_iwc    = 6,  &
       ilw_dlwp   = 7, &
       ilw_diwp   = 8, &
       ilw_dwvp   = 9, &
       ilw_ulwp   = 10, &
       ilw_uiwp   = 11, &
       ilw_uwvp   = 12, &
       ilw_reliq  = 13, &
       ilw_reice  = 14, &
       ilw_o3mr   = 15, &
       ilw_co2    = 16, &
       ilw_ch4    = 17, &
       ilw_n2o    = 18, &
       ilw_z      = 19, &
       ilw_dz     = 20, &
       ilw_dp     = 21, &
       ilw_sza    = 22,  &
       ilw_sfct   = 23,  &
       ilw_emiss  = 24
  character(len=32),dimension(npred_lw) :: &
       pnames_lw =  &
       (/'pressure_pascals                ','temperature_kelvins             ',&
         'specific_humidity_kg_kg01       ','relative_humidity_unitless      ',&
         'liquid_water_content_kg_m03     ','ice_water_content_kg_m03        ',&
         'liquid_water_path_kg_m02        ','ice_water_path_kg_m02           ',&
         'vapour_path_kg_m02              ','upward_liquid_water_path_kg_m02 ',&
         'upward_ice_water_path_kg_m02    ','upward_vapour_path_kg_m02       ',&
         'liquid_effective_radius_metres  ','ice_effective_radius_metres     ',&
         'o3_mixing_ratio_kg_kg01         ','co2_concentration_ppmv          ',&
         'ch4_concentration_ppmv          ','n2o_concentration_ppmv          ',&
         'height_m_agl                    ','height_thickness_metres         ',&
         'pressure_thickness_pascals      ','zenith_angle_radians            ',&
         'surface_temperature_kelvins     ','surface_emissivity              '/)

  ! #####################################################################################
  !
  ! Bookeeping indices 
  !
  ! #####################################################################################
  integer, dimension(npred_lw) :: &
       ip2io_lw ! Index mapping from vector/scalar inputs into LW predictor matrix
  integer, dimension(npred_sw) :: &
       ip2io_sw ! Index mapping from vector/scalar inputs into SW predictor matrix
  logical, dimension(npred_lw) :: &
       is2D_lw  ! Are input LW predictors 2 dimensional?
  logical, dimension(npred_sw) :: &
       is2D_sw  ! Are input SW predictors 2 dimensional?
  
  ! #####################################################################################
  !
  ! Base type containing data needed by emulators.
  !
  ! #####################################################################################
  type ty_mlrad_base_data
     integer :: &
          nCol,              & ! Number of piecewise linear fits of training data(200)
          nLev,              & ! Number of vertical layers in training data (127)
          npreds,            & ! Number of predictor variables (1D)
          npredv               ! Number of predictor variables (2D)
     character(len=31),dimension(:), allocatable :: &
          scalar_pnames,     & ! Name for scalar predictors  [npreds]
          vector_pnames        ! Name for vector predictors  [npredv]
     real(kind_phys),dimension(:,:), allocatable :: &
          scalar_slope,      & ! [nCol,   npreds]
          scalar_intercept,  & ! [nCol,   npreds]
          scalar_breakpoint    ! [nCol+1, npreds]
     real(kind_phys),dimension(:,:,:), allocatable :: &
          vector_slope,      & ! [nCol,   nLev, npredv]
          vector_intercept,  & ! [nCol,   nLev, npredv]
          vector_breakpoint    ! [nCol+1, nLev, npredv]
  end type ty_mlrad_base_data

  public ty_mlrad_data, ty_rad_ml_ref_data

! ########################################################################################
! Type containing training data used for emulator.
!
!> \section arg_table_ty_mlrad_data Argument Table
!! \htmlinclude ty_mlrad_data.html
!!
! ########################################################################################
  type ty_mlrad_data
     type(ty_mlrad_base_data) :: &
          lw, & ! Data for LW emulator.
          sw    ! Data for SW emulator.
   contains
     generic,   public  :: load => load_training_data_piecewise
     procedure, private :: load_training_data_piecewise
  end type ty_mlrad_data

! ########################################################################################
! Type containing reference data for emulator.
!
!> \section arg_table_ty_rad_ml_ref_data Argument Table
!! \htmlinclude ty_rad_ml_ref_data.html
!!
! ########################################################################################
  type ty_rad_ml_ref_data
     integer :: &
          nCol,                   & ! Number of horizontal columns in reference data.
          nLev,                   & ! Number of vertical layers in reference data.
          npred,                  & ! Number of predictor variables.
          ntar_scalar,            & ! Number of targets (scalar)
          ntar_vector               ! Number of targets (vector)
     real(kind_phys),dimension(:,:,:), allocatable :: &
          predictor_matrix,       & ! Predictor matrix (normalized)
          unnorm_predictor_matrix   ! Predictor matrix (unnormalized)
     real(kind_phys),dimension(:,:), allocatable :: &
          scalar_prediction         ! Target predictions (scalar)
     real(kind_phys),dimension(:,:,:), allocatable :: &
          vector_prediction         ! Target predictions (vector)
   contains
     procedure, public :: load =>load_reference_data
  end type ty_rad_ml_ref_data

contains

  ! ######################################################################################
  !
  ! Type-bound procedure to load piecewise linear training data into ty_mlrad_data.
  !
  ! ######################################################################################
  function load_training_data_piecewise(this, file, spect_lw) result(err_message)
    class(ty_mlrad_data), intent(inout) :: this
    ! Inputs
    character(len=*), intent(in) :: file
    logical, intent(in) :: spect_lw
    ! Outputs
    character(len=128) :: err_message
    ! Locals
    integer :: ncid, dimid, varid, iLay, iPred, iName, count, stride

    ! Open file
    call check_netCDF(nf90_open(file, NF90_NOWRITE, ncid),err_message)

    if (spect_lw) then
       ! Get dimensions
       call check_netCDF(nf90_inq_dimid(ncid, 'linear_piece', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%lw%nCol),err_message)
       call check_netCDF(nf90_inq_dimid(ncid, 'height', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%lw%nLev),err_message)
       call check_netCDF(nf90_inq_dimid(ncid, 'scalar_predictor', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%lw%npreds),err_message)
       call check_netCDF(nf90_inq_dimid(ncid, 'vector_predictor', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%lw%npredv),err_message)

       ! Allocate space
       allocate(this%lw%scalar_pnames(this%lw%npreds))
       allocate(this%lw%scalar_slope(     this%lw%nCol,   this%lw%npreds))
       allocate(this%lw%scalar_intercept( this%lw%nCol,   this%lw%npreds))
       allocate(this%lw%scalar_breakpoint(this%lw%nCol+1, this%lw%npreds))
       allocate(this%lw%vector_pnames(this%lw%npredv))
       allocate(this%lw%vector_slope(     this%lw%nCol,   this%lw%nLev, this%lw%npredv))
       allocate(this%lw%vector_intercept( this%lw%nCol,   this%lw%nLev, this%lw%npredv))
       allocate(this%lw%vector_breakpoint(this%lw%nCol+1, this%lw%nLev, this%lw%npredv))

       ! Read in training data (piecewise uniform)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_predictor", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%scalar_pnames),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_predictor", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%vector_pnames),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_slope", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%scalar_slope),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_intercept", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%scalar_intercept),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_break_point_physical_units", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%scalar_breakpoint),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_slope", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%vector_slope),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_intercept", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%vector_intercept),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_break_point_physical_units", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%lw%vector_breakpoint),err_message)

       ! HACK UNTIL NETCDF INPUT FILES ARE FIXED!!
       this%lw%scalar_pnames(1)  = "zenith_angle_radians"
       this%lw%scalar_pnames(2)  = "latitude_deg_n"
       this%lw%scalar_pnames(3)  = "longitude_deg_e"
       this%lw%scalar_pnames(4)  = "column_liquid_water_path_kg_m02"
       this%lw%scalar_pnames(5)  = "column_ice_water_path_kg_m02"
       this%lw%scalar_pnames(6)  = "surface_temperature_kelvins"
       this%lw%scalar_pnames(7)  = "surface_emissivity"
       !
       this%lw%vector_pnames(1)  = "pressure_pascals"
       this%lw%vector_pnames(2)  = "temperature_kelvins"
       this%lw%vector_pnames(3)  = "specific_humidity_kg_kg01"
       this%lw%vector_pnames(4)  = "liquid_water_content_kg_m03"
       this%lw%vector_pnames(5)  = "ice_water_content_kg_m03"
       this%lw%vector_pnames(6)  = "o3_mixing_ratio_kg_kg01"
       this%lw%vector_pnames(7)  = "co2_concentration_ppmv"
       this%lw%vector_pnames(8)  = "ch4_concentration_ppmv"
       this%lw%vector_pnames(9)  = "n2o_concentration_ppmv"
       this%lw%vector_pnames(10) = "liquid_effective_radius_metres"
       this%lw%vector_pnames(11) = "ice_effective_radius_metres"
       this%lw%vector_pnames(12) = "height_m_agl"
       this%lw%vector_pnames(13) = "height_thickness_metres"
       this%lw%vector_pnames(14) = "pressure_thickness_pascals"
       this%lw%vector_pnames(15) = "liquid_water_path_kg_m02"
       this%lw%vector_pnames(16) = "ice_water_path_kg_m02"
       this%lw%vector_pnames(17) = "vapour_path_kg_m02"
       this%lw%vector_pnames(18) = "upward_liquid_water_path_kg_m02"
       this%lw%vector_pnames(19) = "upward_ice_water_path_kg_m02"
       this%lw%vector_pnames(20) = "upward_vapour_path_kg_m02"
       this%lw%vector_pnames(21) = "relative_humidity_unitless"
    else
       ! Get dimensions
       call check_netCDF(nf90_inq_dimid(ncid, 'linear_piece', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%sw%nCol),err_message)
       call check_netCDF(nf90_inq_dimid(ncid, 'height', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%sw%nLev),err_message)
       call check_netCDF(nf90_inq_dimid(ncid, 'scalar_predictor', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%sw%npreds),err_message)
       call check_netCDF(nf90_inq_dimid(ncid, 'vector_predictor', dimid),err_message)
       call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%sw%npredv),err_message)

       ! Allocate space
       allocate(this%sw%scalar_pnames(this%sw%npreds))
       allocate(this%sw%scalar_slope(     this%sw%nCol,   this%sw%npreds))
       allocate(this%sw%scalar_intercept( this%sw%nCol,   this%sw%npreds))
       allocate(this%sw%scalar_breakpoint(this%sw%nCol+1, this%sw%npreds))
       allocate(this%sw%vector_pnames(this%sw%npredv))
       allocate(this%sw%vector_slope(     this%sw%nCol,   this%sw%nLev, this%sw%npredv))
       allocate(this%sw%vector_intercept( this%sw%nCol,   this%sw%nLev, this%sw%npredv))
       allocate(this%sw%vector_breakpoint(this%sw%nCol+1, this%sw%nLev, this%sw%npredv))

       ! Read in training data (piecewise uniform)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_predictor", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%scalar_pnames),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_predictor", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%vector_pnames),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_slope", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%scalar_slope),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_intercept", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%scalar_intercept),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "scalar_break_point_physical_units", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%scalar_breakpoint),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_slope", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%vector_slope),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_intercept", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%vector_intercept),err_message)
       call check_netCDF(nf90_inq_varid(ncid, "vector_break_point_physical_units", varid),err_message)
       call check_netCDF(nf90_get_var(ncid, varid, this%sw%vector_breakpoint),err_message)

       ! HACK UNTIL NETCDF INPUT FILES ARE FIXED!! 
       this%sw%scalar_pnames(1)  = "zenith_angle_radians"
       this%sw%scalar_pnames(2)  = "albedo"
       this%sw%scalar_pnames(3)  = "latitude_deg_n"
       this%sw%scalar_pnames(4)  = "longitude_deg_e"
       this%sw%scalar_pnames(5)  = "column_liquid_water_path_kg_m02"
       this%sw%scalar_pnames(6)  = "column_ice_water_path_kg_m02"
       this%sw%scalar_pnames(7)  = "aerosol_single_scattering_albedo"
       this%sw%scalar_pnames(8)  = "aerosol_asymmetry_param"
       !
       this%sw%vector_pnames(1)  = "pressure_pascals"
       this%sw%vector_pnames(2)  = "temperature_kelvins"
       this%sw%vector_pnames(3)  = "specific_humidity_kg_kg01"
       this%sw%vector_pnames(4)  = "liquid_water_content_kg_m03"
       this%sw%vector_pnames(5)  = "ice_water_content_kg_m03"
       this%sw%vector_pnames(6)  = "o3_mixing_ratio_kg_kg01"
       this%sw%vector_pnames(7)  = "co2_concentration_ppmv"
       this%sw%vector_pnames(8)  = "ch4_concentration_ppmv"
       this%sw%vector_pnames(9)  = "n2o_concentration_ppmv"
       this%sw%vector_pnames(10) = "aerosol_extinction_metres01"
       this%sw%vector_pnames(11) = "liquid_effective_radius_metres"
       this%sw%vector_pnames(12) = "ice_effective_radius_metres"
       this%sw%vector_pnames(13) = "height_m_agl"
       this%sw%vector_pnames(14) = "height_thickness_metres"
       this%sw%vector_pnames(15) = "pressure_thickness_pascals"
       this%sw%vector_pnames(16) = "liquid_water_path_kg_m02"
       this%sw%vector_pnames(17) = "ice_water_path_kg_m02"
       this%sw%vector_pnames(18) = "vapour_path_kg_m02"
       this%sw%vector_pnames(19) = "upward_liquid_water_path_kg_m02"
       this%sw%vector_pnames(20) = "upward_ice_water_path_kg_m02"
       this%sw%vector_pnames(21) = "upward_vapour_path_kg_m02"
       this%sw%vector_pnames(22) = "relative_humidity_unitless"
    endif

    ! Close file
    call check_netCDF(nf90_close(ncid),err_message)

  end function load_training_data_piecewise

  ! ######################################################################################
  !
  ! Type-bound procedure to load reference data into ty_rad_ml_ref_data from input file.
  !
  ! ######################################################################################
  function load_reference_data(this,file) result(err_message)
    class(ty_rad_ml_ref_data), intent(inout) :: this
    ! Inputs
    character(len=*), intent(in) :: file
    ! Output
    character(len=128) :: err_message
    ! Locals
    integer :: ncid, dimid, varid

    ! Open file
    call check_netCDF(nf90_open(file, NF90_NOWRITE, ncid),err_message)

    ! Get dimensions
    call check_netCDF(nf90_inq_dimid(ncid, 'example', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%nCol),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'height', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%nLev),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'predictor_variable', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%npred),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'vector_target_variable', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%ntar_vector),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'scalar_target_variable', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%ntar_scalar),err_message)

    ! Allocate space
    allocate(this%predictor_matrix( this%npred,       this%nLev, this%nCol))
    allocate(this%vector_prediction(this%ntar_vector, this%nLev, this%nCol))
    allocate(this%scalar_prediction(this%ntar_scalar,            this%nCol))
    allocate(this%unnorm_predictor_matrix( this%npred,       this%nLev, this%nCol))

    ! Read in reference data
    call check_netCDF(nf90_inq_varid(ncid, "predictor_matrix", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%predictor_matrix),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "vector_prediction_matrix", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%vector_prediction),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "scalar_prediction_matrix", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%scalar_prediction),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "unnorm_predictor_matrix", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%unnorm_predictor_matrix),err_message)

    ! Close file
    call check_netCDF(nf90_close(ncid),err_message)

  end function load_reference_data

  ! ######################################################################################
  !
  ! Simple function to output netcdf error message and stop code.
  !
  ! ######################################################################################
  subroutine check_netCDF(status,err_message)
    use netcdf
    integer, intent ( in) :: status
    character(len=128),intent(out) :: err_message

    err_message = ''
    if(status /= nf90_noerr) then
       err_message = trim(nf90_strerror(status))
       print*,err_message
       stop
    end if

  end subroutine check_netCDF

end module module_mlrad
