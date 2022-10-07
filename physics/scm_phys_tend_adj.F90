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
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dTdt_lwrad_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dTdt_swrad_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dTdt_pbl_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dqdt_pbl_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dudt_pbl_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dvdt_pbl_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dTdt_gwd_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dudt_gwd_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dvdt_gwd_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dTdt_SCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dudt_SCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dvdt_SCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dqdt_SCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dTdt_DCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dudt_DCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dvdt_DCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dqdt_DCNV_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dTdt_cldMP_data
  real(kind_phys), allocatable, dimension(:,:,:,:) :: dqdt_cldMP_data
  !
  logical :: have_dTdt_lwrad_data = .false.
  logical :: have_dTdt_swrad_data = .false.
  logical :: have_dTdt_pbl_data   = .false.
  logical :: have_dqdt_pbl_data   = .false.
  logical :: have_dudt_pbl_data   = .false.
  logical :: have_dvdt_pbl_data   = .false.
  logical :: have_dTdt_gwd_data   = .false.
  logical :: have_dudt_gwd_data   = .false.
  logical :: have_dvdt_gwd_data   = .false.
  logical :: have_dTdt_SCNV_data  = .false.
  logical :: have_dudt_SCNV_data  = .false.
  logical :: have_dvdt_SCNV_data  = .false.
  logical :: have_dqdt_SCNV_data  = .false.
  logical :: have_dTdt_DCNV_data  = .false.
  logical :: have_dudt_DCNV_data  = .false.
  logical :: have_dvdt_DCNV_data  = .false.
  logical :: have_dqdt_DCNV_data  = .false.
  logical :: have_dTdt_cldMP_data = .false.
  logical :: have_dqdt_cldMP_data = .false.

  public scm_phys_tend_adj_run
contains

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_scm_phys_tend_adj_init
!! \htmlinclude scm_phys_tend_adj_init.html
!!
  !subroutine scm_phys_tend_adj_init(fileIN, errmsg, errflg)
  subroutine scm_phys_tend_adj_init(errmsg, errflg)

    ! Inputs
!    character(len=*),     intent(in ) :: fileIN

    ! Outputs
    character(len=*),     intent(out) :: errmsg
    integer,              intent(out) :: errflg

    ! Local variables
    integer :: ncid, dimID, varID, status, nlon, nlat, nlev, ntime

    !
    ! Open file
    !
    status = nf90_open('/glade/u/home/dswales/testfile.nc', NF90_NOWRITE, ncid)

    !
    ! Dimensions
    !
    status = nf90_inq_dimid(ncid, 'lon', dimid)
    status = nf90_inquire_dimension(ncid, dimid, len = nlon)
    status = nf90_inq_dimid(ncid, 'lat', dimid)
    status = nf90_inquire_dimension(ncid, dimid, len = nlat)
    status = nf90_inq_dimid(ncid, 'lev', dimid)
    status = nf90_inquire_dimension(ncid, dimid, len = nlev)
    status = nf90_inq_dimid(ncid, 'time', dimid)
    status = nf90_inquire_dimension(ncid, dimid, len = ntime)

    !
    ! Read in physics data tendencies
    !
    status = nf90_inq_varid(ncid, 'dtend_temp_lw', varID)
    if (status == nf90_noerror) then
       allocate(dTdt_lwrad_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_lwrad_data)
       have_dTdt_lwrad_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_temp_sw', varID)
    if (status == nf90_noerror) then
       allocate(dTdt_swrad_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_swrad_data)
       have_dTdt_swrad_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_temp_pbl', varID)
    if (status == nf90_noerror) then
       allocate(dTdt_pbl_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_pbl_data)
       have_dTdt_pbl_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_qv_pbl', varID)
    if (status == nf90_noerror) then
       allocate(dqdt_pbl_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dqdt_pbl_data)
       have_dqdt_pbl_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_u_pbl', varID)
    if (status == nf90_noerror) then
       allocate(dudt_pbl_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_pbl_data)
       have_dudt_pbl_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_v_pbl', varID)
    if (status == nf90_noerror) then
       allocate(dvdt_pbl_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_pbl_data)
       have_dvdt_pbl_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_temp_cnvgwd', varID)
    if (status == nf90_noerror) then
       allocate(dTdt_gwd_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_gwd_data)
       have_dTdt_gwd_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_u_cnvgwd', varID)
    if (status == nf90_noerror) then
       allocate(dudt_gwd_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_gwd_data)
       have_dudt_gwd_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_v_cnvgwd', varID)
    if (status == nf90_noerror) then
       allocate(dvdt_gwd_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_gwd_data)
       have_dvdt_gwd_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_temp_shalcnv', varID)
    if (status == nf90_noerror) then
       allocate(dTdt_SCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_SCNV_data)
       have_dTdt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_u_shalcnv', varID)
    if (status == nf90_noerror) then
       allocate(dudt_SCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_SCNV_data)
       have_dudt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_v_shalcnv', varID)
    if (status == nf90_noerror) then
       allocate(dvdt_SCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_SCNV_data)
       have_dvdt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_qv_shalcnv', varID)
    if (status == nf90_noerror) then
       allocate(dqdt_SCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dqdt_SCNV_data)
       have_dqdt_SCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_temp_deepcnv', varID)
    if (status == nf90_noerror) then
       allocate(dTdt_DCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_DCNV_data)
       have_dTdt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_u_deepcnv', varID)
    if (status == nf90_noerror) then
       allocate(dudt_DCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dudt_DCNV_data)
       have_dudt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_v_deepcnv', varID)
    if (status == nf90_noerror) then
       allocate(dvdt_DCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dvdt_DCNV_data)
       have_dvdt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_qv_deepcnv', varID)
    if (status == nf90_noerror) then
       allocate(dqdt_DCNV_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dqdt_DCNV_data)
       have_dqdt_DCNV_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_temp_mp', varID)
    if (status == nf90_noerror) then
       allocate(dTdt_cldMP_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dTdt_cldMP_data)
       have_dTdt_cldMP_data = .true.
    endif
    !
    status = nf90_inq_varid(ncid, 'dtend_qv_mp', varID)
    if (status == nf90_noerror) then
       allocate(dqdt_cldMP_data(nlon, nlat, nlev, ntime))
       status = nf90_get_var(  ncid, varID, dqdt_cldMP_data)
       have_dqdt_cldMP_data = .true.
    endif

    !
    ! Close file
    !
    status = nf90_close(ncid)

  end subroutine scm_phys_tend_adj_init

  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_scm_phys_tend_adj_run
!! \htmlinclude scm_phys_tend_adj_run.html
!!
  subroutine scm_phys_tend_adj_run(dtp, dtf, tgrs, ugrs, vgrs, qgrs, dTdt_lwrad,         &
       dTdt_swrad, dTdt_pbl, dqdt_pbl, dudt_pbl, dvdt_pbl, dTdt_gwd, dudt_gwd, dvdt_gwd, &
       dTdt_SCNV, dqdt_SCNV, dudt_SCNV, dvdt_SCNV, dTdt_DCNV, dqdt_DCNV, dudt_DCNV,      &
       dvdt_DCNV, dTdt_cldMP, dqdt_cldMP, gt0, gu0, gv0, gq0, errmsg, errflg)

    ! Inputs
    real(kind_phys), intent(in   )                   :: dtp, dtf
    real(kind_phys), intent(in   ), dimension(:,:)   :: tgrs, ugrs, vgrs
    real(kind_phys), intent(in   ), dimension(:,:)   :: dTdt_lwrad, dTdt_swrad, dTdt_pbl,&
         dqdt_pbl, dudt_pbl, dvdt_pbl, dTdt_gwd, dudt_gwd, dvdt_gwd, dTdt_SCNV,          &
         dudt_SCNV, dvdt_SCNV, dTdt_DCNV, dudt_DCNV, dvdt_DCNV, dTdt_cldMP
    real(kind_phys), intent(in   ), dimension(:,:,:) :: qgrs, dqdt_DCNV, dqdt_SCNV,      &
         dqdt_cldMP

    ! Outputs
    real(kind_phys), intent(inout), dimension(:,:)   :: gt0, gu0, gv0
    real(kind_phys), intent(inout), dimension(:,:,:) :: gq0
    character(len=*),intent(out  )                   :: errmsg
    integer,         intent(out  )                   :: errflg

    ! Locals
    integer :: iCol, iLay, iTracer
    real(kind_phys), dimension(:,:),   allocatable   :: gt1, gu1, gv1
    real(kind_phys), dimension(:,:,:), allocatable   :: gq1

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
    gt1(:,:)   = tgrs(:,:)   + (dTdt_pbl(:,:) + dTdt_gwd(:,:) + dTdt_SCNV(:,:)    +      &
                                dTdt_DCNV(:,:) + dTdt_lwrad(:,:)+ dTdt_swrad(:,:) +      &
                                dTdt_cldMP(:,:)) * dtp
    gu1(:,:)   = ugrs(:,:)   + (dudt_pbl(:,:) + dudt_gwd(:,:) + dudt_SCNV(:,:)    +      &
                                dudt_DCNV(:,:)) * dtp 
    gv1(:,:)   = vgrs(:,:)   + (dvdt_pbl(:,:) + dvdt_gwd(:,:) + dvdt_SCNV(:,:)    +      &
                                dvdt_DCNV(:,:)) * dtp 
    gq1(:,:,1) = qgrs(:,:,1) + (dqdt_pbl(:,:) + dqdt_SCNV(:,:,1) + dqdt_DCNV(:,:,1)) * dtp
    do iTracer=2,size(gq0(1,1,:))
       gq1(:,:,iTracer) = qgrs(:,:,iTracer) + (dqdt_SCNV(:,:,iTracer)             +      &
            dqdt_DCNV(:,:,iTracer) + dqdt_cldMP(:,:,iTracer)) * dtp
    enddo

    do iCol=1,size(gq0(:,1,1))
       do iLay=1,size(gq0(1,:,1))
          write(*,'(i5,4f8.3)') iLay,tgrs(iCol,iLay),gt0(iCol,iLay),gt1(iCol,iLay),gt0(iCol,iLay)-gt1(iCol,iLay)
       enddo
    enddo

  end subroutine scm_phys_tend_adj_run

end module scm_phys_tend_adj
