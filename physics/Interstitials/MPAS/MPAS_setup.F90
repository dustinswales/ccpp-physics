!> \file MPAS_setup.F90
!! This file initializes

module MPAS_setup
  use mpi_f08
  use ccpp_wp, only : kind_phys

  implicit none

  public MPAS_setup_init

  private

contains

!> \section arg_table_MPAS_setup_init Argument Table
!! \htmlinclude GFS_rrtmgp_setup_init.html
!!
  subroutine MPAS_setup_init(errmsg, errflg)
    ! Outputs
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine MPAS_setup_init

end module MPAS_setup
