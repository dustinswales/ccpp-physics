module rrtmgp_aux
  use machine, only: &
       kind_phys                   ! Working type
  implicit none

  real(kind_phys) :: &
       rrtmgp_minP, & ! Minimum pressure allowed in RRTMGP
       rrtmgp_minT    ! Minimum temperature allowed in RRTMGP
contains
  ! #########################################################################################
  ! SUBROUTINE check_error_msg
  ! #########################################################################################
  subroutine check_error_msg(routine_name, error_msg, error_flag, error_str)
    character(len=*), intent(in) :: &
         error_msg, routine_name
    integer,intent(out) :: error_flag
    character(len=*), intent(out) :: error_str
    
    error_str  = ""
    error_flag = 0
    if(error_msg /= "") then
       error_str  = "ERROR("//trim(routine_name)//"): "//trim(error_msg)
       error_flag = 1
       return
    end if
  end subroutine check_error_msg  
end module rrtmgp_aux
