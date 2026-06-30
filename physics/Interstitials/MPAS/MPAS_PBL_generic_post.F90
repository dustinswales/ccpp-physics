!> \file MPAS_PBL_generic_post.F90
!!

! #########################################################################################
!> Contains code related to PBL schemes to be called prior to PBL schemes within
!! MPAS-based physics suites.
! #########################################################################################
module MPAS_PBL_generic_post
  implicit none

contains

!! \htmlinclude MPAS_PBL_generic_post_timestep_init.html
!!
  subroutine MPAS_PBL_generic_post_timestep_init (errmsg, errflg)
    use ccpp_wp,  only : kind_phys

    implicit none
    ! CCPP error handling variables
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    print*,'CCPP: Calling MPAS_PBL_generic_post_timestep_init()'
  end subroutine MPAS_PBL_generic_post_timestep_init
  
! #########################################################################################
!> \brief This scheme sets up the vertically diffused tracer array for any PBL scheme based
!! on the microphysics scheme chosen
!! \section arg_table_MPAS_PBL_generic_post_run Argument Table
!!  
!! \htmlinclude MPAS_PBL_generic_post_run.html
!!
! #########################################################################################
  subroutine MPAS_PBL_generic_post_run (errmsg, errflg)
    use ccpp_wp,  only : kind_phys

    implicit none

    ! CCPP error handling variables
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    print*,'CCPP: Calling MPAS_PBL_generic_post_run()'
  end subroutine MPAS_PBL_generic_post_run

end module MPAS_PBL_generic_post
