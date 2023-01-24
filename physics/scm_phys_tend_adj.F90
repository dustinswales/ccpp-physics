! ########################################################################################
! ########################################################################################
module scm_phys_tend_adj
  use machine, only: kind_phys
  use netcdf
#ifdef MPI
  use mpi
#endif

  !
  ! Data driven phsyics tendencies
  !
  real(kind_phys), allocatable, dimension(:)   :: time_data
  real(kind_phys), allocatable, dimension(:,:) :: dTdt_lwrad_data, dTdt_swrad_data,      &
       dTdt_pbl_data, dqdt_pbl_data, dudt_pbl_data, dvdt_pbl_data, dTdt_gwd_data,        &
       dudt_gwd_data, dvdt_gwd_data, dTdt_SCNV_data, dudt_SCNV_data, dvdt_SCNV_data,     &
       dqdt_SCNV_data, dTdt_DCNV_data, dudt_DCNV_data, dvdt_DCNV_data, dqdt_DCNV_data,   &
       dTdt_cldMP_data, dqdt_cldMP_data
  !
  ! Logical switches for each scheme and state.
  !
  logical :: have_dTdt_lwrad_data = .false.,    use_dTdt_lwrad_data = .false.,           &
             have_dTdt_swrad_data = .false.,    use_dTdt_swrad_data = .false.,           &
             have_dTdt_pbl_data   = .false.,    use_dTdt_pbl_data   = .false.,           &
             have_dqdt_pbl_data   = .false.,    use_dqdt_pbl_data   = .false.,           &
             have_dudt_pbl_data   = .false.,    use_dudt_pbl_data   = .false.,           &
             have_dvdt_pbl_data   = .false.,    use_dvdt_pbl_data   = .false.,           &
             have_dTdt_gwd_data   = .false.,    use_dTdt_gwd_data   = .false.,           &
             have_dudt_gwd_data   = .false.,    use_dudt_gwd_data   = .false.,           &
             have_dvdt_gwd_data   = .false.,    use_dvdt_gwd_data   = .false.,           &
             have_dTdt_SCNV_data  = .false.,    use_dTdt_SCNV_data  = .false.,           &
             have_dudt_SCNV_data  = .false.,    use_dudt_SCNV_data  = .false.,           &
             have_dvdt_SCNV_data  = .false.,    use_dvdt_SCNV_data  = .false.,           &
             have_dqdt_SCNV_data  = .false.,    use_dqdt_SCNV_data  = .false.,           &
             have_dTdt_DCNV_data  = .false.,    use_dTdt_DCNV_data  = .false.,           &
             have_dudt_DCNV_data  = .false.,    use_dudt_DCNV_data  = .false.,           &
             have_dvdt_DCNV_data  = .false.,    use_dvdt_DCNV_data  = .false.,           &
             have_dqdt_DCNV_data  = .false.,    use_dqdt_DCNV_data  = .false.,           &
             have_dTdt_cldMP_data = .false.,    use_dTdt_cldMP_data = .false.,           &
             have_dqdt_cldMP_data = .false.,    use_dqdt_cldMP_data = .false.

  public scm_phys_tend_adj_init, scm_phys_tend_adj_run
contains

  ! ######################################################################################
  !
  ! SUBROUTINE scm_phys_tend_adj_init
  !
  ! ######################################################################################
!! \section arg_table_scm_phys_tend_adj_init
!! \htmlinclude scm_phys_tend_adj_init.html
!!
  subroutine scm_phys_tend_adj_init(me, master, nlunit, nml_file, errmsg, errflg)

    ! Inputs
    integer,          intent (in) :: me, master, nlunit
    character(len=*), intent (in) :: nml_file

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    integer :: ncid, dimID, varID, status, nlon, nlat, nlev, ntime, ios, init_year,     &
         init_month, init_day, init_hour
    character(len=128) :: fileIN
    logical :: exists

    ! Namelist
    namelist / scm_data_nml / &
         fileIN, use_dTdt_lwrad_data, use_dTdt_swrad_data, use_dTdt_pbl_data,           &
         use_dqdt_pbl_data, use_dudt_pbl_data, use_dvdt_pbl_data, use_dTdt_gwd_data,    &
         use_dudt_gwd_data, use_dvdt_gwd_data, use_dTdt_SCNV_data, use_dudt_SCNV_data,  &
         use_dvdt_SCNV_data, use_dqdt_SCNV_data, use_dTdt_DCNV_data, use_dudt_DCNV_data,&
         use_dvdt_DCNV_data, use_dqdt_DCNV_data, use_dTdt_cldMP_data,                   &
         use_dqdt_cldMP_data

    !
    ! Read in namelist
    !
    inquire (file = trim (nml_file), exist = exists)
    if (.not. exists) then
        errmsg = 'SCM data tendency :: namelist file: '//trim(nml_file)//' does not exist'
        errflg = 1
        return
    else
        open (unit = nlunit, file = nml_file, action = 'read', status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read (nlunit, nml = scm_data_nml)
    close (nlunit)

    !
    if (me == 0) then
       print*, "--- Using SCM data tendencies ---"
       print*, "year:            ", init_year
       print*, "month:           ", init_month
       print*, "day:             ", init_day
       print*, "hour:            ", init_hour
       print*, "---------------------------------"
       print*, "                 AV RQ"
       print*, "dTdt_lwrad_data: ", have_dTdt_lwrad_data,  use_dTdt_lwrad_data
       print*, "dTdt_swrad_data: ", have_dTdt_swrad_data,  use_dTdt_swrad_data
       print*, "dTdt_pbl_data:   ", have_dTdt_pbl_data,    use_dTdt_pbl_data
       print*, "dqdt_pbl_data:   ", have_dqdt_pbl_data,    use_dqdt_pbl_data
       print*, "dudt_pbl_data:   ", have_dudt_pbl_data,    use_dudt_pbl_data
       print*, "dvdt_pbl_data:   ", have_dvdt_pbl_data,    use_dvdt_pbl_data
       print*, "dTdt_gwd_data:   ", have_dTdt_gwd_data,    use_dTdt_gwd_data
       print*, "dudt_gwd_data:   ", have_dudt_gwd_data,    use_dudt_gwd_data
       print*, "dvdt_gwd_data:   ", have_dvdt_gwd_data,    use_dvdt_gwd_data
       print*, "dTdt_SCNV_data:  ", have_dTdt_SCNV_data,   use_dTdt_SCNV_data
       print*, "dudt_SCNV_data:  ", have_dudt_SCNV_data,   use_dudt_SCNV_data
       print*, "dvdt_SCNV_data:  ", have_dvdt_SCNV_data,   use_dvdt_SCNV_data
       print*, "dqdt_SCNV_data:  ", have_dqdt_SCNV_data,   use_dqdt_SCNV_data
       print*, "dTdt_DCNV_data:  ", have_dTdt_DCNV_data,   use_dTdt_DCNV_data
       print*, "dudt_DCNV_data:  ", have_dudt_DCNV_data,   use_dudt_DCNV_data
       print*, "dvdt_DCNV_data:  ", have_dvdt_DCNV_data,   use_dvdt_DCNV_data
       print*, "dqdt_DCNV_data:  ", have_dqdt_DCNV_data,   use_dqdt_DCNV_data
       print*, "dTdt_cldMP_data: ", have_dTdt_cldMP_data,  use_dTdt_cldMP_data
       print*, "dqdt_cldMP_data: ", have_dqdt_cldMP_data,  use_dqdt_cldMP_data
       print*, "---------------------------------"
    endif

  end subroutine scm_phys_tend_adj_init

  ! ######################################################################################
  !
  ! SUBROUTINE scm_phys_tend_adj_run
  !
  ! ######################################################################################
!! \section arg_table_scm_phys_tend_adj_run
!! \htmlinclude scm_phys_tend_adj_run.html
!!
  subroutine scm_phys_tend_adj_run(solhr, kdt, dtp, dtf, tgrs, ugrs, vgrs, qgrs,         &
       dTdt_lwrad, dTdt_swrad, dTdt_pbl, dqdt_pbl, dudt_pbl, dvdt_pbl, dTdt_gwd,         &
       dudt_gwd, dvdt_gwd, dTdt_SCNV, dqdt_SCNV, dudt_SCNV, dvdt_SCNV, dTdt_DCNV,        &
       dqdt_DCNV, dudt_DCNV, dvdt_DCNV, dTdt_cldMP, dqdt_cldMP,                          &
       gt0, gu0, gv0, gq0, errmsg, errflg)

    ! Inputs
    integer,         intent(in   ) :: kdt
    real(kind_phys), intent(in   ) :: dtp, dtf, solhr
    real(kind_phys), intent(in   ), dimension(:,:) :: tgrs, ugrs, vgrs
    real(kind_phys), intent(in   ), dimension(:,:) :: dTdt_lwrad, dTdt_swrad, dTdt_pbl,  &
         dudt_pbl, dvdt_pbl, dTdt_gwd, dudt_gwd, dvdt_gwd, dTdt_SCNV, dudt_SCNV,         &
         dvdt_SCNV, dTdt_DCNV, dudt_DCNV, dvdt_DCNV, dTdt_cldMP
    real(kind_phys), intent(in   ), dimension(:,:,:) :: qgrs, dqdt_DCNV, dqdt_SCNV,      &
         dqdt_cldMP, dqdt_pbl

    ! Outputs
    real(kind_phys), intent(inout), dimension(:,:) :: gt0, gu0, gv0
    real(kind_phys), intent(inout), dimension(:,:,:) :: gq0
    character(len=*),intent(out  ) :: errmsg
    integer,         intent(out  ) :: errflg

    ! Locals
    integer :: iCol, iLay, iTracer, nCol, nLay, nTracer, ti(1), tf(1)
    real(kind_phys) :: w1, w2
    real(kind_phys), dimension(:,:),   allocatable :: gt1, gu1, gv1
    real(kind_phys), dimension(:,:,:), allocatable :: gq1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Dimensions
    nCol = size(gq0(:,1,1))
    nLay = size(gq0(1,:,1))
    nTrc = size(gq0(1,1,:))
    
    ! Allocate temporaries
    allocate(gt1(nCol,nLay), gu1(nCol,nLay), gv1(nCol,nLay),gq1(nCol,nLay,nTrc))

    !
    ! Determine temporal interpolation weights for data-tendecies.
    !
    ti = findloc(abs(time_data-kdt*dtp),minval(abs(time_data-kdt*dtp)))
    if (kdt*dtp - time_data(ti(1)) .le. 0) ti = ti-1
    tf = ti + 1
    w1 = (time_data(tf(1))-kdt*dtp) / (time_data(tf(1)) - time_data(ti(1)))
    w2 = 1 - w1

    do iCol = 1,nCol
       !
       ! Temperature
       !
       gt1(iCol,:) = tgrs(iCol,:)
       ! Longwave radiation
       if (have_dTdt_lwrad_data .and. use_dTdt_lwrad_data) then
          gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_lwrad_data(ti(1),:) + (1-w1)*dTdt_lwrad_data(tf(1),:)) * dtp
       else
          gt1(iCol,:) = gt1(iCol,:) + dTdt_lwrad(iCol,:) * dtp
       endif
       ! Shortwave radiation
       if (have_dTdt_swrad_data .and. use_dTdt_swrad_data) then
          gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_swrad_data(ti(1),:) + (1-w1)*dTdt_swrad_data(tf(1),:)) * dtp
       else
          gt1(iCol,:) = gt1(iCol,:) + dTdt_swrad(iCol,:) * dtp
       endif
       ! PBL physics
       if (have_dTdt_pbl_data .and. use_dTdt_pbl_data) then
          gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_pbl(ti(1),:) + (1-w1)*(dTdt_pbl(tf(1),:))) * dtp
       else
          gt1(iCol,:) = gt1(iCol,:) + dTdt_pbl(iCol,:) * dtp
       endif
       ! Gravity-wave drag (DJS Need to break into oro/non-oro components)
       if (have_dTdt_gwd_data .and. use_dTdt_gwd_data) then
          gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_gwd(ti(1),:) + (1-w1)*(dTdt_gwd(tf(1),:))) * dtp
       else
          gt1(iCol,:) = gt1(iCol,:) + dTdt_gwd(iCol,:) * dtp
       endif
       ! Shallow convection
       if (have_dTdt_SCNV_data .and. use_dTdt_SCNV_data) then
          gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_SCNV(ti(1),:) + (1-w1)*(dTdt_SCNV(tf(1),:))) * dtp
       else
          gt1(iCol,:) = gt1(iCol,:) + dTdt_SCNV(iCol,:) * dtp
       endif
       ! Deep convection
       if (have_dTdt_DCNV_data .and. use_dTdt_DCNV_data) then
          gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_DCNV(ti(1),:) + (1-w1)*(dTdt_DCNV(tf(1),:) )) * dtp
       else
          gt1(iCol,:) = gt1(iCol,:) + dTdt_DCNV(iCol,:) * dtp
       endif
       ! Cloud macrophysics
       if (have_dTdt_cldMP_data .and. use_dTdt_cldMP_data) then
          gt1(iCol,:) = gt1(iCol,:) + (w1*dTdt_cldMP(ti(1),:) + (1-w1)*(dTdt_cldMP(tf(1),:))) * dtp
       else
          gt1(iCol,:) = gt1(iCol,:) + dTdt_cldMP(iCol,:) * dtp
       endif
       
       !
       ! u-momentum
       !
       gu1(iCol,:) = ugrs(iCol,:)
       ! PBL physics
       if (have_dudt_pbl_data .and. use_dudt_pbl_data) then
          gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_pbl(ti(1),:) + (1-w1)*(dudt_pbl(tf(1),:))) * dtp
       else
          gu1(iCol,:) = gu1(iCol,:) + dudt_pbl(iCol,:) * dtp
       endif
       ! Gravity-wave drag (DJS Need to break into oro/non-oro components)
       if (have_dudt_gwd_data .and. use_dudt_gwd_data) then
          gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_gwd(ti(1),:) + (1-w1)*(dudt_gwd(tf(1),:))) * dtp
       else
          gu1(iCol,:) = gu1(iCol,:) + dudt_gwd(iCol,:) * dtp
       endif
       ! Shallow convection
       if (have_dudt_SCNV_data .and. use_dudt_SCNV_data) then
          gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_SCNV(ti(1),:) + (1-w1)*(dudt_SCNV(tf(1),:))) * dtp
       else
          gu1(iCol,:) = gu1(iCol,:) + dudt_SCNV(iCol,:) * dtp
       endif
       ! Deep convection 
       if (have_dudt_DCNV_data .and. use_dudt_DCNV_data) then
          gu1(iCol,:) = gu1(iCol,:) + (w1*dudt_DCNV(ti(1),:) + (1-w1)*(dudt_DCNV(tf(1),:))) * dtp
       else
          gu1(iCol,:) = gu1(iCol,:) + dudt_DCNV(iCol,:) * dtp
       endif
       
       !
       ! v-momentum
       !
       gv1(iCol,:) = vgrs(iCol,:)
       ! PBL physics 
       if (have_dvdt_pbl_data .and. use_dvdt_pbl_data) then
          gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_pbl(ti(1),:) + (1-w1)*(dvdt_pbl(tf(1),:))) * dtp
       else
          gv1(iCol,:) = gv1(iCol,:) + dvdt_pbl(iCol,:) * dtp
       endif
       ! Gravity-wave drag (DJS Need to break into oro/non-oro components)
       if (have_dvdt_gwd_data .and. use_dvdt_gwd_data) then
          gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_gwd(ti(1),:) + (1-w1)*(dvdt_gwd(tf(1),:))) * dtp
       else
          gv1(iCol,:) = gv1(iCol,:) + dvdt_gwd(iCol,:) * dtp
       endif
       ! Shallow convection
       if (have_dvdt_SCNV_data .and. use_dvdt_SCNV_data) then
          gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_SCNV(ti(1),:) + (1-w1)*(dvdt_SCNV(tf(1),:))) * dtp
       else
          gv1(iCol,:) = gv1(iCol,:) + dvdt_SCNV(iCol,:) * dtp
       endif
       ! Deep convection 
       if (have_dvdt_DCNV_data .and. use_dvdt_DCNV_data) then
          gv1(iCol,:) = gv1(iCol,:) + (w1*dvdt_DCNV(ti(1),:) + (1-w1)*(dvdt_DCNV(tf(1),:) )) * dtp
       else
          gv1(iCol,:) = gv1(iCol,:) + dvdt_DCNV(iCol,:) * dtp
       endif
       
       !
       ! Moisture
       !
!       gq1(iCol,:,1) = qgrs(iCol,:,1)
!       ! PBL physics
!       if (have_dqdt_pbl_data .and. use_dqdt_pbl_data) then
!          gq1(iCol,:,1) = gq1(iCol,:,1) + (w1*dqdt_pbl(ti(1),:) + (1-w1)*(dqdt_pbl(tf(1),:))) * dtp
!       else
!          gq1(iCol,:,1) = gq1(iCol,:,1) + dqdt_pbl(iCol,:) * dtp
!       endif
!       ! Shallow convection
!       if (have_dqdt_SCNV_data .and. use_dqdt_SCNV_data) then
!          gq1(iCol,:,1) = gq1(iCol,:,1) + (w1*dqdt_SCNV(ti(1),:,1) + (1-w1)*(dqdt_SCNV(tf(1),:,1))) * dtp
!       else
!          gq1(iCol,:,1) = gq1(iCol,:,1) + dqdt_SCNV(iCol,:,1) * dtp
!       endif
!       ! Deep convection 
!       if (have_dqdt_DCNV_data .and. use_dqdt_DCNV_data) then
!          gq1(iCol,:,1) = gq1(iCol,:,1) + (w1*dqdt_DCNV(ti(1),:,1) + (1-w1)*(dqdt_DCNV(tf(1),:,1))) * dtp
!       else
!          gq1(iCol,:,1) = gq1(iCol,:,1) + dqdt_DCNV(iCol,:,1) * dtp
!       endif
       
       !
       ! Moisture(1) + Tracers(2:ntracer)
       !
       do iTracer=1,nTrc
          gq1(iCol,:,iTracer) = qgrs(iCol,:,iTracer)
          ! PBL
          if (have_dqdt_PBL_data .and. use_dqdt_PBL_data) then
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + (w1*dqdt_PBL(ti(1),:,iTracer) + (1-w1)*(dqdt_PBL(tf(1),:,iTracer))) * dtp
          else
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + dqdt_PBL(iCol,:,iTracer) * dtp
          endif
          ! Shallow convection 
          if (have_dqdt_SCNV_data .and. use_dqdt_SCNV_data) then
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + (w1*dqdt_SCNV(ti(1),:,iTracer) + (1-w1)*(dqdt_SCNV(tf(1),:,iTracer))) * dtp
          else
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + dqdt_SCNV(iCol,:,iTracer) * dtp
          endif
          ! Deep convection 
          if (have_dqdt_DCNV_data .and. use_dqdt_DCNV_data) then
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + (w1*dqdt_DCNV(ti(1),:,iTracer) + (1-w1)*(dqdt_DCNV(tf(1),:,iTracer))) * dtp
          else
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + dqdt_DCNV(iCol,:,iTracer) * dtp
          endif
          ! Cloud macrophysics
          if (have_dqdt_cldMP_data .and. use_dqdt_cldMP_data) then
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + (w1*dqdt_cldMP(ti(1),:,iTracer) + (1-w1)*(dqdt_cldMP(tf(1),:,iTracer))) * dtp
          else
             gq1(iCol,:,iTracer) = gq1(iCol,:,iTracer) + dqdt_cldMP(iCol,:,iTracer) * dtp
          endif
       enddo  ! tracers
    enddo     ! columns

    !
    do iCol=1,size(gq0(:,1,1))
       do iLay=1,size(gq0(1,:,1))
          write(*,'(i5,4f8.3)') iLay,tgrs(iCol,iLay),gt0(iCol,iLay),gt1(iCol,iLay),gt0(iCol,iLay)-gt1(iCol,iLay)
       enddo
    enddo

  end subroutine scm_phys_tend_adj_run

  subroutine nf90_error_reporting(status, errmsg, errflg)
    integer,          intent(in ) :: status
    integer,          intent(out) :: errflg
    character(len=*), intent(out) :: errmsg
    !
    if (status /= nf90_noerr) then
       errflg = 1
       errmsg = trim(nf90_strerror(status))
    endif
    !
  end subroutine nf90_error_reporting

end module scm_phys_tend_adj
