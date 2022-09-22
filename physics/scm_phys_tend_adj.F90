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
  subroutine scm_phys_tend_adj_run(dtp, dtf, tgrs, ugrs, vgrs, qgrs, dTdt_lwrad,         &
       dTdt_swrad, dTdt_pbl, dqdt_pbl, dudt_pbl, dvdt_pbl, dTdt_gwd, dudt_gwd, dvdt_gwd, &
       dTdt_SCNV, dqdt_SCNV, dudt_SCNV, dvdt_SCNV, dTdt_DCNV, dqdt_DCNV, dudt_DCNV,      &
       dvdt_DCNV, gt0, gu0, gv0, gq0, errmsg, errflg)

    ! Inputs
    real(kind_phys), intent(in   )                   :: dtp, dtf
    real(kind_phys), intent(in   ), dimension(:,:)   :: tgrs, ugrs, vgrs
    real(kind_phys), intent(in   ), dimension(:,:,:) :: qgrs, dqdt_DCNV, dqdt_SCNV
    real(kind_phys), intent(in   ), dimension(:,:)   :: dTdt_lwrad, dTdt_swrad, dTdt_pbl, &
                                                        dqdt_pbl,   dudt_pbl,   dvdt_pbl, &
                                                        dTdt_gwd,   dudt_gwd,   dvdt_gwd, &
                                                        dTdt_SCNV,  dudt_SCNV,  dvdt_SCNV,&
                                                        dTdt_DCNV,  dudt_DCNV,  dvdt_DCNV
    ! Outputs
    real(kind_phys), intent(inout), dimension(:,:)   :: gt0, gu0, gv0
    real(kind_phys), intent(inout), dimension(:,:,:) :: gq0
    character(len=*),intent(out  )                   :: errmsg
    integer,         intent(out  )                   :: errflg

    ! Locals
    integer :: iCol, iLay, iTracer
    real(kind_phys), dimension(:,:), allocatable :: gt1, gu1, gv1
    real(kind_phys), dimension(:,:,:), allocatable :: gq1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    allocate(gt1(size(gq0(:,1,1)),size(gq0(1,:,1))), &
             gu1(size(gq0(:,1,1)),size(gq0(1,:,1))), &
             gv1(size(gq0(:,1,1)),size(gq0(1,:,1))), &
             gq1(size(gq0(:,1,1)),size(gq0(1,:,1)),size(gq0(1,1,:))))
    !
    ! Reconstruct state using scheme tendencies...
    !
    gt1(:,:)   = tgrs(:,:)   + (dTdt_pbl(:,:) + dTdt_gwd(:,:) + dTdt_SCNV(:,:)   + &
                                dTdt_DCNV(:,:) + dTdt_lwrad(:,:)+ dTdt_swrad(:,:)) * dtp
    gu1(:,:)   = ugrs(:,:)   + (dudt_pbl(:,:) + dudt_gwd(:,:) + dudt_SCNV(:,:)   + &
                                dudt_DCNV(:,:)) * dtp 
    gv1(:,:)   = vgrs(:,:)   + (dvdt_pbl(:,:) + dvdt_gwd(:,:) + dvdt_SCNV(:,:)   + &
                                dvdt_DCNV(:,:)) * dtp 
    gq1(:,:,1) = qgrs(:,:,1) + (dqdt_pbl(:,:) + dqdt_SCNV(:,:,1) + dqdt_DCNV(:,:,1)) * dtp
    do iTracer=2,size(gq0(1,1,:))
       gq1(:,:,iTracer) = qgrs(:,:,iTracer) + (dqdt_SCNV(:,:,iTracer) + dqdt_DCNV(:,:,iTracer)) * dtp
    enddo

    do iCol=1,size(gq0(:,1,1))
       do iLay=1,size(gq0(1,:,1))
          write(*,'(i5,3f8.3)') iLay,gt0(iCol,iLay),gt1(iCol,iLay),gt0(iCol,iLay)-gt1(iCol,iLay)
       enddo
    enddo

  end subroutine scm_phys_tend_adj_run

end module scm_phys_tend_adj
