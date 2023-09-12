module module_rad_ml
!> \section arg_table_module_rad_ml Argument table
!! \htmlinclude module_rad_ml.html
!!
  use machine, only : kind_phys
  use netcdf
  implicit none

  ! #####################################################################################
  !
  ! Predictor varaible names, from training data, in order expected by emulator.
  !
  ! #####################################################################################
  integer,parameter :: nPredictors = 24
  character(len=32),dimension(nPredictors) :: predictor_names =  &
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
  integer, dimension(nPredictors) :: &
       ip2io   ! Index mapping from vector/scalar inputs into predictor matrix
  logical, dimension(nPredictors) :: &
       is2D    ! Are input predictors 2 dimensional?

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
          nCol,                   & ! Number of training columns in training data 
          nLev,                   & ! Number of vertical layers in training data
          npred_scalar,           & ! Number of predictor variables (1D)
          npred_vector,           & ! Number of predictor variables (2D)
          field_id_char             ! Character length
     character(len=32),dimension(:), allocatable :: &
          scalar_predictor_names, & ! Name for scalar predictors  [npred_scalar]
          vector_predictor_names    ! Name for vector predictors  [npred_vector]
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
   contains
     procedure, public :: load => load_training_data
     procedure, public :: vert_int
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
  ! Type-bound procedure to load training data into ty_rad_ml_data from input file.
  !
  ! ######################################################################################
  function load_training_data(this, file, nPts) result(err_message)
    class(ty_rad_ml_data), intent(inout) :: this
    ! Inputs
    character(len=*), intent(in) :: file
    integer, intent(in), optional :: nPts
    ! Outputs
    character(len=128) :: err_message
    ! Locals
    integer :: ncid, dimid, varid, iLay, iPred, iName, count, stride
    
    ! Open file
    call check_netCDF(nf90_open(file, NF90_NOWRITE, ncid),err_message)

    ! Get dimensions
    call check_netCDF(nf90_inq_dimid(ncid, 'example', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%nCol),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'height', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%nLev),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'scalar_predictor', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%npred_scalar),err_message)
    call check_netCDF(nf90_inq_dimid(ncid, 'vector_predictor', dimid),err_message)
    call check_netCDF(nf90_inquire_dimension(ncid, dimid, len = this%npred_vector),err_message)

    ! If provided, use only a subset of the points
    if (present(nPts)) then
       stride    = max(1,this%nCol/nPts)
       this%nCol = nPts
    endif
    
    ! Allocate space
    allocate(this%scalar_predictor_names(this%npred_scalar))
    allocate(this%vector_predictor_names(this%npred_vector))
    allocate(this%scalar_predictor(this%npred_scalar, this%nCol))
    allocate(this%vector_predictor(this%npred_vector, this%nLev, this%nCol))
    allocate(this%scalar_predictor_mean(this%npred_scalar))
    allocate(this%vector_predictor_mean(this%npred_vector, this%nLev))
    allocate(this%scalar_predictor_stdev(this%npred_scalar))
    allocate(this%vector_predictor_stdev(this%npred_vector, this%nLev))

    ! Read in training data
    call check_netCDF(nf90_inq_varid(ncid, "scalar_predictor_names", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%scalar_predictor_names),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "vector_predictor_names", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%vector_predictor_names),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "scalar_predictor_matrix", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%scalar_predictor,count=(/this%npred_scalar, this%nCol/),stride=(/1,stride/)),err_message)
    call check_netCDF(nf90_inq_varid(ncid, "vector_predictor_matrix", varid),err_message)
    call check_netCDF(nf90_get_var(ncid, varid, this%vector_predictor,count=(/this%npred_vector, this%nLev, this%nCol/),stride=(/1,1,stride/)),err_message)
    
    ! HACK UNTIL NETCDF INPUT FILES ARE FIXED!!
    this%scalar_predictor_names(1)  = "zenith_angle_radians"
    this%scalar_predictor_names(2)  = "latitude_deg_n"
    this%scalar_predictor_names(3)  = "longitude_deg_e"
    this%scalar_predictor_names(4)  = "column_liquid_water_path_kg_m02"
    this%scalar_predictor_names(5)  = "column_ice_water_path_kg_m02"
    this%scalar_predictor_names(6)  = "surface_temperature_kelvins"
    this%scalar_predictor_names(7)  = "surface_emissivity"
    !
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

    ! Close file
    call check_netCDF(nf90_close(ncid),err_message)

    ! Compute scalar-predictor statistics.
    do iPred = 1,this%npred_scalar
       this%scalar_predictor_mean(iPred)  = sum(this%scalar_predictor(iPred,:))/this%nCol
       this%scalar_predictor_stdev(iPred) = sqrt(abs(sum(this%scalar_predictor(iPred,:)**2)-&
            sum(this%scalar_predictor(iPred,:))**2/this%nCol)/(this%nCol-1))
    enddo

    ! Compute vector-predictor statistics.
    do iPred = 1,this%npred_vector
       do iLay = 1,this%nLev
          this%vector_predictor_mean(iPred,iLay)  = sum(this%vector_predictor(iPred,iLay,:))/this%nCol
          this%vector_predictor_stdev(iPred,iLay) = sqrt(abs(sum(this%vector_predictor(iPred,iLay,:)**2)-&
               sum(this%vector_predictor(iPred,iLay,:))**2/this%nCol)/(this%nCol-1))
       enddo
    enddo

  end function load_training_data

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
  ! Procedure (type-bound) to interpolate model field onto emulators expected grid.
  ! ######################################################################################
  subroutine vert_int(this, var_in, zi, var_out, zo)
    class(ty_rad_ml_data), intent(in) :: this
    real(kind_phys),dimension(:),intent(in) :: &
         var_in, & ! Field to be regridded
         zi,     & ! Height of field, at interface (m)
         zo        ! Height of regridded field, at interface (m)
    real(kind_phys),dimension(:),intent(out) :: &
         var_out   ! Regriided field
    integer :: nlay_in, nlev_in, ilay_in, nlay_out, nlev_out, ilay_out, num_w
    real(kind_phys) :: wts, wt, dbb, dbt, dtb, dtt

    var_out(:) = 0._kind_phys
    nlay_in  = size(var_in)
    nlev_in  = size(zi)
    nlay_out = size(var_out)
    nlev_out = nlay_out+1

    do iLay_out = 1,nlay_out
       num_w = 0
       wts   = 0._kind_phys
       do iLay_in=1,nlay_in
          wt = 0._kind_phys
          if (ilay_in > nlev_in) exit
          ! Distances between edges of both grids
          dbb = zi(iLay_in)   - zo(iLay_out)
          dtb = zi(iLay_in+1) - zo(iLay_out)
          dbt = zi(iLay_in)   - zo(iLay_out+1)
          dtt = zi(iLay_in+1) - zo(iLay_out+1)
          if (dbt >= 0.0) exit ! Do next level in the new grid
          if (dtb > 0.0) then
             if (dbb <= 0.0) then
                if (dtt <= 0) then
                   wt = dtb
                else
                   wt = zo(iLay_out+1) - zo(iLay_out)
                endif
             else
                if (dtt <= 0) then
                   wt = zi(iLay_in+1) - zi(iLay_in)
                else
                   wt = -dbt
                endif
             endif
             ! If layers overlap (w/=0), then accumulate
             if (wt /= 0.0) then
                num_w = num_w + 1
                wts = wts + wt
                var_out(iLay_out) = var_out(iLay_out) + wt*var_in(iLay_in)
             endif
          endif
       enddo

       ! Calculate average in new grid
       if (num_w > 0) then
          var_out(iLay_out) = var_out(iLay_out)/wts
       endif
    enddo
    var_out(nlay_out) = var_in(nlay_in) - (zi(nLay_in) - zo(nLay_out))*(var_in(nLay_in) - var_in(nLay_in-1)) / (zi(nLay_in)-zi(nLay_in-1))


  end subroutine vert_int

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
       stop
    end if

  end subroutine check_netCDF

end module module_rad_ml
