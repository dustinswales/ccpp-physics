! ###########################################################################################
!> \file ml_rad_phys.F90
!!
! ###########################################################################################
module ml_rad_phys
  use machine, only: kind_phys
  use netcdf
  use inferof
  use iso_c_binding, only : c_double, c_int, c_float, c_char, c_null_char, c_ptr
  implicit none

  type(infero_model) :: model

  public ml_rad_phys_init, ml_rad_phys_run

contains
! ########################################################################################
!! \section arg_table_ml_rad_phys_init
!! \htmlinclude ml_rad_phys_init.html
!!
! #########################################################################################
  subroutine ml_rad_phys_init(do_ml_rad, errmsg, errflg)
    ! Inputs
    logical,           intent(in) :: do_ml_rad

    ! Outputs
    character(len=*),  intent(out) :: errmsg
    integer,           intent(out) :: errflg

    ! Localc
    character (len = *), parameter :: MODEL_PATH = "/scratch1/BMC/gmtb/Dustin.Swales/ML-radiation/libs/Jebb/infero/infero/gsl/model.onnx"
    character (len = *), parameter :: MODEL_TYPE = "onnx"
    character(1024) :: yaml_config

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ml_rad) return

    call infero_check(infero_initialise())

    ! YAML config string
    yaml_config = "---"//NEW_LINE('A') &
         //"  path: "//TRIM(MODEL_PATH)//NEW_LINE('A') &
         //"  type: "//TRIM(MODEL_TYPE)//c_null_char

    ! Get a infero model
    call infero_check(model%initialise_from_yaml_string(yaml_config))

  end subroutine ml_rad_phys_init

! #########################################################################################
!! \section arg_table_ml_rad_phys_run
!! \htmlinclude ml_rad_phys_run.html
!!
! #########################################################################################
  subroutine ml_rad_phys_run(do_ml_rad, predictor_matrix, vector_prediction_matrix,       &
       scalar_prediction_matrix, errmsg, errflg)
    ! Inputs
    logical,          intent(in) :: &
         do_ml_rad
    real(kind_phys),  intent(in), dimension(:,:,:) :: &
         predictor_matrix, &
         vector_prediction_matrix
    real(kind_phys),  intent(in), dimension(:,:) :: &
         scalar_prediction_matrix

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg
    integer,          intent(out) :: &
         errflg

    ! Locals
    integer :: npred_ml, nlev_ml, ncase_ml, ipred, ilev, icase, iinf
    real(c_float), allocatable :: it2f(:,:,:) ! data for inference in profile, height, value order
    real(c_float), allocatable :: ot2f(:,:)   ! data from inference in profile, height order

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. do_ml_rad) return

    print*,'predictor_matrix:         ', shape(predictor_matrix)
    print*,'vector_prediction_matrix: ', shape(vector_prediction_matrix)
    print*,'scalar_prediction_matrix: ', shape(scalar_prediction_matrix)

    ! Dimensions
    npred_ml = size(predictor_matrix(:,0,0))
    nlev_ml  = size(predictor_matrix(0,:,0))
    ncase_ml = size(predictor_matrix(0,0,:))

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

!    ! Check the data
!    print *, SHAPE(vector_prediction_matrix), "  ", SHAPE(scalar_prediction_matrix), "  ", SHAPE(ot2f)
!    do icase = 1, ncase_ml
!       do ilev = 1, nlev_ml + 2
!          do ipred = 1, 1
!             if (ilev <= nlev_ml) then
!                if (abs(ot2f(icase,ilev) - vector_prediction_matrix(ipred,ilev,icase)) .gt. 1e-3) then
!                   write(*,*) "ERROR: output element ",icase,ilev, " (", ot2f(icase,ilev) ,") ", &
!                        "is different from expected value ", vector_prediction_matrix(ipred,ilev,icase)
!                   stop 1
!                end if
!                print *, "( ",icase , ", ", ilev, ", ", ipred,") = ", ot2f(icase, ilev), vector_prediction_matrix(ipred, ilev, icase)
!             else
!                if (abs(ot2f(icase,ilev) - scalar_prediction_matrix(ilev-127,icase)) .gt. 1e-3) then
!                   write(*,*) "ERROR: output element ",icase,ilev, " (", ot2f(icase,ilev) ,") ", &
!                        "is different from expected value ", scalar_prediction_matrix(ilev-127,icase)
!                   stop 1
!                end if
!                print *,  "( ",icase , ", ", ilev, ", ", ipred,") = ", ot2f(icase, ilev), scalar_prediction_matrix(ilev-127, icase)
!             end if
!          end do
!       end do
!    end do

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
