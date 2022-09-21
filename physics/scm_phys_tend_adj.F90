! ########################################################################################
! ########################################################################################
module scm_phys_tend_adj
  use machine, only: kind_phys
  use netcdf
#ifdef MPI
  use mpi
#endif
  
  public scm_phys_tend_adj_run
contains

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_scm_phys_tend_adj_init
!! \htmlinclude scm_phys_tend_adj_init.html
!!
  subroutine scm_phys_tend_adj_init(errmsg, errflg)
    ! Outputs
    character(len=*),     intent(out) :: errmsg
    integer,              intent(out) :: errflg

  end subroutine scm_phys_tend_adj_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_scm_phys_tend_adj_run
!! \htmlinclude scm_phys_tend_adj_run.html
!!
  subroutine scm_phys_tend_adj_run(dtp, dtf, tgrs, ugrs, vgrs, qgrs, dTdt_lwrad, dTdt_swrad, &
       dTdt_pbl, dqdt_pbl, dudt_pbl, dvdt_pbl, dTdt_gwd, dudt_gwd, dvdt_gwd, gt0, gu0, gv0, gq0, errmsg, errflg)

    ! Inputs
    real(kind_phys), intent(in   )                   :: dtp, dtf
    real(kind_phys), intent(in   ), dimension(:,:)   :: tgrs, ugrs, vgrs
    real(kind_phys), intent(in   ), dimension(:,:,:) :: qgrs
    real(kind_phys), intent(in   ), dimension(:,:)   :: dTdt_lwrad, dTdt_swrad, dTdt_pbl, &
                                                        dqdt_pbl, dudt_pbl, dvdt_pbl,     &
                                                        dTdt_gwd, dudt_gwd, dvdt_gwd
    ! Outputs
    real(kind_phys), intent(inout), dimension(:,:)   :: gt0, gu0, gv0
    real(kind_phys), intent(inout), dimension(:,:,:) :: gq0
    character(len=*),intent(out  )                   :: errmsg
    integer,         intent(out  )                   :: errflg

    ! Locals
    integer :: iCol, iLay

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    do iCol=1,size(tgrs(:,1))
       do iLay=1,size(tgrs(1,:))
          print*,tgrs(iCol,iLay),dTdt_lwrad(iCol,iLay)*dtp + dTdt_swrad(iCol,iLay)*dtp
       enddo
    enddo
    !gt0(:,:)   = tgrs(:,:)   + dtdt(:,:)   * dtp
    !gu0(:,:)   = ugrs(:,:)   + dudt(:,:)   * dtp
    !gv0(:,:)   = vgrs(:,:)   + dvdt(:,:)   * dtp
    !gq0(:,:,:) = qgrs(:,:,:) + dqdt(:,:,:) * dtp

  end subroutine scm_phys_tend_adj_run

end module scm_phys_tend_adj
