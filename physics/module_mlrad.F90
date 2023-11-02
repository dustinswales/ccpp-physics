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
  integer,parameter :: nPredictors_lw = 24, nPredictors_sw = 26
  character(len=32),dimension(nPredictors_lw) :: &
       pnames_lw =  &
       (/'pressure_pascals                 ','temperature_kelvins             ',&
         'specific_humidity_kg_kg01        ','relative_humidity_unitless      ',&
         'liquid_water_content_kg_m03      ','ice_water_content_kg_m03        ',&
         'liquid_water_path_kg_m02         ','ice_water_path_kg_m02           ',&
         'vapour_path_kg_m02               ','upward_liquid_water_path_kg_m02 ',&
         'upward_ice_water_path_kg_m02     ','upward_vapour_path_kg_m02       ',&
         'liquid_effective_radius_metres   ','ice_effective_radius_metres     ',&
         'o3_mixing_ratio_kg_kg01          ','co2_concentration_ppmv          ',&
         'ch4_concentration_ppmv           ','n2o_concentration_ppmv          ',&
         'height_m_agl                     ','height_thickness_metres         ',&
         'pressure_thickness_pascals       ','zenith_angle_radians            ',&
         'surface_temperature_kelvins      ','surface_emissivity              '/)
  character(len=32),dimension(nPredictors_sw) :: &
       pnames_sw =  &
       (/'pressure_pascals                 ','temperature_kelvins             ',&
         'specific_humidity_kg_kg01        ','relative_humidity_unitless      ',&
         'liquid_water_content_kg_m03      ','ice_water_content_kg_m03        ',&
         'liquid_water_path_kg_m02         ','ice_water_path_kg_m02           ',&
         'vapour_path_kg_m02               ','upward_liquid_water_path_kg_m02 ',&
         'upward_ice_water_path_kg_m02     ','upward_vapour_path_kg_m02       ',&
         'liquid_effective_radius_metres   ','ice_effective_radius_metres     ',&
         'o3_mixing_ratio_kg_kg01          ','co2_concentration_ppmv          ',&
         'ch4_concentration_ppmv           ','n2o_concentration_ppmv          ',&
         'aerosol_extinction_metres01      ','height_m_agl                    ',&
         'height_thickness_metres          ','pressure_thickness_pascals      ',&
         'zenith_angle_radians             ','albedo                          ',&
         'aerosol_single_scattering_albedo ','aerosol_asymmetry_param         '/)

  ! Bookeeping indices 
  integer, dimension(nPredictors_lw) :: &
       ip2io_lw ! Index mapping from vector/scalar inputs into LW predictor matrix
  integer, dimension(nPredictors_sw) :: &
       ip2io_sw ! Index mapping from vector/scalar inputs into SW predictor matrix
  logical, dimension(nPredictors_lw) :: &
       is2D_lw  ! Are input LW predictors 2 dimensional?
  logical, dimension(nPredictors_sw) :: &
       is2D_sw  ! Are input SW predictors 2 dimensional?
  
  ! Base type containing data needed by emulators.
  type ty_mlrad_base_data
     integer :: &
          nCol,             & ! Number of piecewise linear fits of training data(200)
          nLev,             & ! Number of vertical layers in training data (127)
          npreds,           & ! Number of predictor variables (1D)
          npredv              ! Number of predictor variables (2D)
     character(len=31),dimension(:), allocatable :: &
          scalar_pnames, & ! Name for scalar predictors  [npreds]
          vector_pnames    ! Name for vector predictors  [npredv]
     real(kind_phys),dimension(:,:), allocatable :: &
          scalar_slope,           & ! [nCol,   npreds]
          scalar_intercept,       & ! [nCol,   npreds]
          scalar_breakpoint         ! [nCol+1, npreds]
     real(kind_phys),dimension(:,:,:), allocatable :: &
          vector_slope,           & ! [nCol,   nLev, npredv]
          vector_intercept,       & ! [nCol,   nLev, npredv]
          vector_breakpoint         ! [nCol+1, nLev, npredv]
  end type ty_mlrad_base_data

  public ty_mlrad_data

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
