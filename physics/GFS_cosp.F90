!> \file GFS_cosp.F90
module GFS_cosp
  use machine,  only: kind_phys
  use mod_cosp, only: cosp_init,cosp_optical_inputs,cosp_column_inputs,cosp_outputs,cosp_cleanUp,cosp_simulator

  implicit none
  public GFS_cosp_init, GFS_cosp_run
contains

!! \htmlinclude GFS_cosp_init.html
  subroutine GFS_cosp_init( errmsg, errflg)
    
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg             ! Error message
    integer, intent(out) :: &
         errflg             ! Error flag
  end subroutine GFS_cosp_init


!! \htmlinclude GFS_cosp_run.html
  subroutine GFS_cosp_run( errmsg, errflg)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg             ! Error message
    integer, intent(out) :: &
         errflg             ! Error flag

  end subroutine GFS_cosp_run


end module GFS_cosp
