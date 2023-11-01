module module_mlrad
!> \section arg_table_module_mlrad Argument table
!! \htmlinclude module_mlrad.html
!!
  use machine, only : kind_phys
  use netcdf
  implicit none

  ! #####################################################################################
  !
  ! Predictor varaible names, from training data, in order expected by emulator.
  !
  ! #####################################################################################
  integer,parameter :: nPredictors_lw = 24
  integer, dimension(nPredictors_lw) :: &
       ip2io_lw   ! Index mapping from vector/scalar inputs into predictor matrix
  logical, dimension(nPredictors_lw) :: &
       is2D_lw    ! Are input predictors 2 dimensional?
  character(len=32),dimension(nPredictors_lw) :: predictor_names_lw =  &
       (/'pressure_pascals                 ','temperature_kelvins             ',        &
         'specific_humidity_kg_kg01        ','relative_humidity_unitless      ',        &
         'liquid_water_content_kg_m03      ','ice_water_content_kg_m03        ',        &
         'liquid_water_path_kg_m02         ','ice_water_path_kg_m02           ',        &
         'vapour_path_kg_m02               ','upward_liquid_water_path_kg_m02 ',        &
         'upward_ice_water_path_kg_m02     ','upward_vapour_path_kg_m02       ',        &
         'liquid_effective_radius_metres   ','ice_effective_radius_metres     ',        &
         'o3_mixing_ratio_kg_kg01          ','co2_concentration_ppmv          ',        &
         'ch4_concentration_ppmv           ','n2o_concentration_ppmv          ',        &
         'height_m_agl                     ','height_thickness_metres         ',        &
         'pressure_thickness_pascals       ','zenith_angle_radians            ',        &
         'surface_temperature_kelvins      ','surface_emissivity              '/)
  public ty_rad_ml_data, ty_rad_ml_ref_data

! ########################################################################################
! Type containing training data used for emulator.
!
!> \section arg_table_ty_rad_ml_data Argument Table
!! \htmlinclude ty_rad_ml_data.html
!!
! ########################################################################################
  type ty_rad_ml_data
     integer :: &
          nCol,                   & ! Number of training columns OR piecewise linear fit
                                    ! of training data 
          nLev,                   & ! Number of vertical layers in training data
          npred_scalar,           & ! Number of predictor variables (1D)
          npred_vector,           & ! Number of predictor variables (2D)
          field_id_char             ! Character length
     character(len=31),dimension(:), allocatable :: &
          scalar_predictor_names, & ! Name for scalar predictors  [npred_scalar]
          vector_predictor_names    ! Name for vector predictors  [npred_vector]
     ! Used when reading in full training data
     real(kind_phys),dimension(:,:), allocatable :: &
          scalar_predictor          ! Scalar predictors           [npred_scalar, nCol]
     real(kind_phys),dimension(:,:,:), allocatable :: &
          vector_predictor          ! Vector predictors           [npred_vector, nLev, nCol]
     real(kind_phys),dimension(:), allocatable :: &
          scalar_predictor_mean,  & ! Mean of scalar predictors   [npred_scalar]
          scalar_predictor_stdev    ! Stdev of scalar predictors  [npred_scalar]
     real(kind_phys),dimension(:,:), allocatable :: &
          vector_predictor_mean,  & ! Mean of vector predictors   [npred_vector, nLev]
          vector_predictor_stdev    ! Stdev of vector predictors  [npred_vector, nLev]
     ! Used for piecewise discretization of traning data
     real(kind_phys),dimension(:,:), allocatable :: &
          scalar_slope,           & ! [nCol,   npred_scalar]
          scalar_intercept,       & ! [nCol,   npred_scalar]
          scalar_breakpoint         ! [nCol+1, npred_scalar]
     real(kind_phys),dimension(:,:,:), allocatable :: &
          vector_slope,           & ! [nCol,   nLev, npred_vector]
          vector_intercept,       & ! [nCol,   nLev, npred_vector]
          vector_breakpoint         ! [nCol+1, nLev, npred_vector]
   contains
     generic,   public  :: load => load_training_data_piecewise
     procedure, private :: load_training_data_piecewise
  end type ty_rad_ml_data

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
  ! Type-bound procedure to load piecewise linear training data into ty_rad_ml_data.
  !
  ! ######################################################################################
  function load_training_data_piecewise(this, file, spect_lw) result(err_message)
    class(ty_rad_ml_data), intent(inout) :: this
    ! Inputs
    character(len=*), intent(in) :: file
    logical, intent(in) :: spect_lw
    ! Outputs
    character(len=128) :: err_message
    ! Locals
    integer :: ncid, dimid, varid, iLay, iPred, iName, count, stride

    ! Open file
    call check_netCDF(nf90_open(file, NF90_NOWRITE, ncid),err_message)

    ! Get dimensions
    call check_netCDF(nf90_inq_dimid(ncid, 'linear_piece', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%nCol),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'height', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%nLev),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'scalar_predictor', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%npred_scalar),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'vector_predictor', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%npred_vector),err_message)

    ! Allocate space
    allocate(this%scalar_predictor_names(this%npred_scalar))
    allocate(this%scalar_slope(     this%nCol,   this%npred_scalar))
    allocate(this%scalar_intercept( this%nCol,   this%npred_scalar))
    allocate(this%scalar_breakpoint(this%nCol+1, this%npred_scalar))
    allocate(this%vector_predictor_names(this%npred_vector))
    allocate(this%vector_slope(     this%nCol,   this%nLev, this%npred_vector))
    allocate(this%vector_intercept( this%nCol,   this%nLev, this%npred_vector))
    allocate(this%vector_breakpoint(this%nCol+1, this%nLev, this%npred_vector))

    ! Read in training data (piecewise uniform)
    call check_netCDF(nf90_inq_varid(ncid, "scalar_predictor", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%scalar_predictor_names),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "vector_predictor", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%vector_predictor_names),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "scalar_slope", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%scalar_slope),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "scalar_intercept", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%scalar_intercept),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "scalar_break_point_physical_units", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%scalar_breakpoint),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "vector_slope", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%vector_slope),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "vector_intercept", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%vector_intercept),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "vector_break_point_physical_units", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%vector_breakpoint),err_message)

    ! HACK UNTIL NETCDF INPUT FILES ARE FIXED!!
    if (spect_lw) then
       this%scalar_predictor_names(1)  = "zenith_angle_radians"
       this%scalar_predictor_names(2)  = "latitude_deg_n"
       this%scalar_predictor_names(3)  = "longitude_deg_e"
       this%scalar_predictor_names(4)  = "column_liquid_water_path_kg_m02"
       this%scalar_predictor_names(5)  = "column_ice_water_path_kg_m02"
       this%scalar_predictor_names(6)  = "surface_temperature_kelvins"
       this%scalar_predictor_names(7)  = "surface_emissivity"
    else
       this%scalar_predictor_names(1)  = "zenith_angle_radians"
       this%scalar_predictor_names(2)  = "albedo"
       this%scalar_predictor_names(3)  = "latitude_deg_n"
       this%scalar_predictor_names(4)  = "longitude_deg_e"
       this%scalar_predictor_names(5)  = "column_liquid_water_path_kg_m02"
       this%scalar_predictor_names(6)  = "column_ice_water_path_kg_m02"
       this%scalar_predictor_names(7)  = "aerosol_single_scattering_albedo"
       this%scalar_predictor_names(8)  = "aerosol_asymmetry_param"
    endif
    if (spect_lw) then
       this%vector_predictor_names(1)  = "pressure_pascals"
       this%vector_predictor_names(2)  = "temperature_kelvins"
       this%vector_predictor_names(3)  = "specific_humidity_kg_kg01"
       this%vector_predictor_names(4)  = "liquid_water_content_kg_m03"
       this%vector_predictor_names(5)  = "ice_water_content_kg_m03"
       this%vector_predictor_names(6)  = "o3_mixing_ratio_kg_kg01"
       this%vector_predictor_names(7)  = "co2_concentration_ppmv"
       this%vector_predictor_names(8)  = "ch4_concentration_ppmv"
       this%vector_predictor_names(9)  = "n2o_concentration_ppmv"
       this%vector_predictor_names(10) = "liquid_effective_radius_metres"
       this%vector_predictor_names(11) = "ice_effective_radius_metres"
       this%vector_predictor_names(12) = "height_m_agl"
       this%vector_predictor_names(13) = "height_thickness_metres"
       this%vector_predictor_names(14) = "pressure_thickness_pascals"
       this%vector_predictor_names(15) = "liquid_water_path_kg_m02"
       this%vector_predictor_names(16) = "ice_water_path_kg_m02"
       this%vector_predictor_names(17) = "vapour_path_kg_m02"
       this%vector_predictor_names(18) = "upward_liquid_water_path_kg_m02"
       this%vector_predictor_names(19) = "upward_ice_water_path_kg_m02"
       this%vector_predictor_names(20) = "upward_vapour_path_kg_m02"
       this%vector_predictor_names(21) = "relative_humidity_unitless"
    else
       this%vector_predictor_names(1)  = "pressure_pascals"
       this%vector_predictor_names(2)  = "temperature_kelvins"
       this%vector_predictor_names(3)  = "specific_humidity_kg_kg01"
       this%vector_predictor_names(4)  = "liquid_water_content_kg_m03"
       this%vector_predictor_names(5)  = "ice_water_content_kg_m03"
       this%vector_predictor_names(6)  = "o3_mixing_ratio_kg_kg01"
       this%vector_predictor_names(7)  = "co2_concentration_ppmv"
       this%vector_predictor_names(8)  = "ch4_concentration_ppmv"
       this%vector_predictor_names(9)  = "n2o_concentration_ppmv"
       this%vector_predictor_names(10) = "aerosol_extinction_metres01"
       this%vector_predictor_names(11) = "liquid_effective_radius_metres"
       this%vector_predictor_names(12) = "ice_effective_radius_metres"
       this%vector_predictor_names(13) = "height_m_agl"
       this%vector_predictor_names(14) = "height_thickness_metres"
       this%vector_predictor_names(15) = "pressure_thickness_pascals"
       this%vector_predictor_names(16) = "liquid_water_path_kg_m02"
       this%vector_predictor_names(17) = "ice_water_path_kg_m02"
       this%vector_predictor_names(18) = "vapour_path_kg_m02"
       this%vector_predictor_names(19) = "upward_liquid_water_path_kg_m02"
       this%vector_predictor_names(20) = "upward_ice_water_path_kg_m02"
       this%vector_predictor_names(21) = "upward_vapour_path_kg_m02"
       this%vector_predictor_names(22) = "relative_humidity_unitless"
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
