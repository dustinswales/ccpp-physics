! ###########################################################################################
!> \file ml_rad_phys.F90
!!
! ###########################################################################################
module ml_rad_phys
  use machine, only: kind_phys
  use netcdf
  use inferof
  implicit none

  public ml_rad_phys_init, ml_rad_phys_run
contains
! ########################################################################################
!! \section arg_table_ml_rad_phys_init
!! \htmlinclude ml_rad_phys_init.html
!!
! #########################################################################################
  subroutine ml_rad_phys_init(errmsg,errflg)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
  end subroutine ml_rad_phys_init

! #########################################################################################
!! \section arg_table_ml_rad_phys_run
!! \htmlinclude ml_rad_phys_run.html
!!
! #########################################################################################
  subroutine ml_rad_phys_run(errmsg,errflg)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
  end subroutine ml_rad_phys_run
  
end module ml_rad_phys
