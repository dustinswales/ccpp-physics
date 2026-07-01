!> \file MPAS_PBL_generic_pre.F90
!!

! #########################################################################################
!> Contains code related to PBL schemes to be called prior to PBL schemes within
!! MPAS-based physics suites.
! #########################################################################################
module MPAS_PBL_generic_pre
  implicit none

contains

! #########################################################################################
!> \brief This scheme sets up the vertically diffused tracer array for any PBL scheme based
!! on the microphysics scheme chosen
!! \section arg_table_MPAS_PBL_generic_pre_run Argument Table
!!  
!! \htmlinclude MPAS_PBL_generic_pre_run.html
!!
! #########################################################################################
  subroutine MPAS_PBL_generic_pre_run (ncols, nlevs, nvdiff, rtg_ozone_index,             &
       ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc, ntwa, ntia, ntgl, ntoz,         &
       imp_physics, imp_physics_thompson, ltaerosol, mraerosol,                           &
       qgrs, ugrs, vgrs, tgrs, vdftra, errmsg, errflg)
    use ccpp_wp,  only : kind_phys

    implicit none

    integer, intent(out) :: rtg_ozone_index
    integer, intent(in) :: ncols, nlevs, nvdiff
    integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntrnc
    integer, intent(in) :: ntwa, ntia, ntgl, ntoz
    integer, intent(in) :: imp_physics, imp_physics_thompson
    logical, intent(in) :: ltaerosol, mraerosol

    real(kind=kind_phys), dimension(:,:,:), intent(in) :: qgrs
    real(kind=kind_phys), dimension(:,:  ), intent(in) :: ugrs, vgrs, tgrs
    real(kind=kind_phys), dimension(:,:,:), intent(inout) :: vdftra

    ! CCPP error handling variables
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg


    ! Local variables
    integer :: i, k, kk, k1, n

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    rtg_ozone_index=-1

    ! Thompson MP
    if (imp_physics == imp_physics_thompson) then
       if(ltaerosol) then
          do k=1,nlevs
             do i=1,ncols
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntoz)
                vdftra(i,k,11) = qgrs(i,k,ntwa)
                vdftra(i,k,12) = qgrs(i,k,ntia)
             enddo
          enddo
          rtg_ozone_index = 10
       elseif(mraerosol) then
          do k=1,nlevs
             do i=1,ncols
                vdftra(i,k,1)  = qgrs(i,k,ntqv)
                vdftra(i,k,2)  = qgrs(i,k,ntcw)
                vdftra(i,k,3)  = qgrs(i,k,ntiw)
                vdftra(i,k,4)  = qgrs(i,k,ntrw)
                vdftra(i,k,5)  = qgrs(i,k,ntsw)
                vdftra(i,k,6)  = qgrs(i,k,ntgl)
                vdftra(i,k,7)  = qgrs(i,k,ntlnc)
                vdftra(i,k,8)  = qgrs(i,k,ntinc)
                vdftra(i,k,9)  = qgrs(i,k,ntrnc)
                vdftra(i,k,10) = qgrs(i,k,ntoz)
             enddo
          enddo
          rtg_ozone_index = 10
       else
          do k=1,nlevs
             do i=1,ncols
                vdftra(i,k,1) = qgrs(i,k,ntqv)
                vdftra(i,k,2) = qgrs(i,k,ntcw)
                vdftra(i,k,3) = qgrs(i,k,ntiw)
                vdftra(i,k,4) = qgrs(i,k,ntrw)
                vdftra(i,k,5) = qgrs(i,k,ntsw)
                vdftra(i,k,6) = qgrs(i,k,ntgl)
                vdftra(i,k,7) = qgrs(i,k,ntinc)
                vdftra(i,k,8) = qgrs(i,k,ntrnc)
                vdftra(i,k,9) = qgrs(i,k,ntoz)
             enddo
          enddo
          rtg_ozone_index = 9
       endif
    endif ! END MP

    end subroutine MPAS_PBL_generic_pre_run

  end module MPAS_PBL_generic_pre
