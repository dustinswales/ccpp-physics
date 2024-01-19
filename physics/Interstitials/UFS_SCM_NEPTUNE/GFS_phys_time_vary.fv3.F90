! #########################################################################################
!> \file GFS_phys_time_vary.fv3.F90
!!  Contains code related to GFS physics suite setup (physics part of time_vary_step)

!>\defgroup mod_GFS_phys_time_vary GFS Physics Time Update
!! This module contains GFS physics time vary subroutines including stratospheric water vapor,
!! aerosol, IN&CCN and surface properties updates.
! #########################################################################################
module GFS_phys_time_vary
#ifdef _OPENMP
  use omp_lib
#endif
  use machine,                    only : kind_phys, kind_dbl_prec, kind_sngl_prec
  use mersenne_twister,           only : random_setseed, random_number, random_index,       &
                                         random_stat
  use module_radiation_astronomy, only : sol_init
  use module_radiation_aerosols,  only : aer_init
  use module_radiation_gases,     only : gas_init
  use module_radiation_clouds,    only : cld_init
  use rrtmg_lw,                   only : rlwinit
  use rrtmg_sw,                   only : rswinit
  use module_ozphys,              only : ty_ozphys
  use h2o_def,                    only : levh2o, h2o_coeff, h2o_lat, h2o_pres, h2o_time,    &
                                         h2oplin
  use h2ointerp,                  only : read_h2odata, setindxh2o, h2ointerpol
  use aerclm_def,                 only : aerin, aer_pres, ntrcaer, ntrcaerm, iamin, iamax,  &
                                         jamin, jamax
  use aerinterp,                  only : read_aerdata, setindxaer, aerinterpol, read_aerdataf
  use iccn_def,                   only : ciplin, ccnin, ci_pres
  use iccninterp,                 only : read_cidata, setindxci, ciinterpol
  use gcycle_mod,                 only : gcycle
  use funcphys,                   only : gfuncphys
  use cires_tauamf_data,          only : cires_indx_ugwp,  read_tau_amf, tau_amf_interp
  use cires_tauamf_data,          only : tau_limb,  days_limb, ugwp_taulat
  use namelist_soilveg,           only : salp_data, snupx
  use set_soilveg_mod,            only : set_soilveg
  use noahmp_tables,              only : read_mp_table_parameters, laim_table, saim_table,  &
                                         sla_table, bexp_table, smcmax_table,smcwlt_table,  &
                                         dwsat_table, dksat_table, psisat_table,            &
                                         isurban_table, isbarren_table, isice_table,        &
                                         iswater_table

  implicit none

  private copy_error

  public GFS_phys_time_vary_init, GFS_phys_time_vary_timestep_init,                          &
       GFS_phys_time_vary_timestep_finalize, GFS_phys_time_vary_finalize

  ! Module paramaters
  integer ::  &
       month0 = 0, &
       iyear0 = 0, &
       monthd = 0
  logical :: is_initialized = .false.
  ! Control flag for the first time of reading climatological ozone data
  ! (set/reset in subroutines GFS_phys_time_vary_init/GFS_phys_time_vary_timestep_init,
  ! it is used only if the control parameter ntoz=0)
  logical :: loz1st = .true.

  !
  real(kind=kind_phys), parameter :: con_hr        =  3600.0_kind_phys
  real(kind=kind_phys), parameter :: con_24        =   24.0_kind_phys
  real(kind=kind_phys), parameter :: con_99        =    99.0_kind_phys
  real(kind=kind_phys), parameter :: con_100       =   100.0_kind_phys
  real(kind=kind_phys), parameter :: missing_value = 9.99e20_kind_phys
  real(kind=kind_phys), parameter :: drythresh     =   1.e-4_kind_phys
  real(kind=kind_phys), parameter :: zero          =     0.0_kind_phys
  real(kind=kind_phys), parameter :: one           =     1.0_kind_phys

contains
  ! #########################################################################################
  ! SUBROUTINE GFS_phys_time_vary_init
  ! #########################################################################################
!> \section arg_table_GFS_phys_time_vary_init Argument Table
!! \htmlinclude GFS_phys_time_vary_init.html
!!
!>\section gen_GFS_phys_time_vary_init GFS_phys_time_vary_init General Algorithm
!> @{
  subroutine GFS_phys_time_vary_init (me, master, ntoz, h2o_phys, iccn, iflip, im, levs,    &
       nx, ny, idate, xlat_d, xlon_d,                                                       &
       jindx1_o3, jindx2_o3, ddy_o3, jindx1_h, jindx2_h, ddy_h, h2opl,fhour,                &
       jindx1_aer, jindx2_aer, ddy_aer, iindx1_aer, iindx2_aer, ddx_aer, aer_nm,            &
       jindx1_ci, jindx2_ci, ddy_ci, iindx1_ci, iindx2_ci, ddx_ci, imap, jmap,              &
       do_ugwp_v1, jindx1_tau, jindx2_tau, ddy_j1tau, ddy_j2tau,                            &
       isot, ivegsrc, nlunit, sncovr, sncovr_ice, lsm, lsm_noahmp, lsm_ruc, min_seaice,     &
       fice, landfrac, vtype, weasd, lsoil, zs, dzs, lsnow_lsm_lbound, lsnow_lsm_ubound,    &
       tvxy, tgxy, tahxy, canicexy, canliqxy, eahxy, cmxy, chxy, fwetxy, sneqvoxy, alboldxy,&
       qsnowxy, wslakexy, albdvis_lnd, albdnir_lnd, albivis_lnd, albinir_lnd, albdvis_ice,  &
       albdnir_ice, albivis_ice, albinir_ice, emiss_lnd, emiss_ice, taussxy, waxy, wtxy,    &
       zwtxy, xlaixy, xsaixy, lfmassxy, stmassxy, rtmassxy, woodxy, stblcpxy, fastcpxy,     &
       smcwtdxy, deeprechxy, rechxy, snowxy, snicexy, snliqxy, tsnoxy , smoiseq, zsnsoxy,   &
       slc, smc, stc, tsfcl, snowd, canopy, tg3, stype, lsm_cold_start, nthrds, lkm,        &
       use_lake_model, lakefrac, lakedepth, iopt_lake, iopt_lake_clm, iopt_lake_flake,      &
       lakefrac_threshold, lakedepth_threshold, ozphys, iaermdl, iaerflg, si, levr, ictm,   &
       isol, solar_file, ico2, iaer, ntcw, num_p3d, npdf3d, iovr, iovr_rand, iovr_maxrand,  &
       iovr_max, iovr_dcorr, iovr_exp, iovr_exprand, icliq_sw, lcrick, lcnorm, imp_physics, &
       lnoprec, do_RRTMGP, lalw1bd, aeros_file, con_pi, con_t0c, con_c,                 &
       con_boltz, con_plnk, con_solr_2008, con_solr_2002, con_g, con_rd, co2usr_file,       &
       co2cyc_file, rad_hr_units, inc_minor_gas, icliq_lw, isubcsw, isubclw, iswmode,       &
       ipsd0, ltp, lextop, iaerclm, errmsg, errflg)

    implicit none

    ! Interface variables
    integer,              intent(in)    :: me, master, ntoz, iccn, iflip, im, nx, ny, levs
    logical,              intent(in)    :: h2o_phys, lsm_cold_start
    integer,              intent(in)    :: idate(:), iopt_lake, iopt_lake_clm, iopt_lake_flake
    real(kind_phys),      intent(in)    :: fhour, lakefrac_threshold, lakedepth_threshold
    real(kind_phys),      intent(in)    :: xlat_d(:), xlon_d(:)
    integer,              intent(in)    :: lkm
    integer,              intent(inout) :: use_lake_model(:)
    real(kind=kind_phys), intent(in   ) :: lakefrac(:), lakedepth(:)
    integer,              intent(inout) :: jindx1_o3(:), jindx2_o3(:), jindx1_h(:), jindx2_h(:)
    real(kind_phys),      intent(inout) :: ddy_o3(:),  ddy_h(:)
    real(kind_phys),      intent(in)    :: h2opl(:,:,:)
    integer,              intent(inout) :: jindx1_aer(:), jindx2_aer(:), iindx1_aer(:), iindx2_aer(:)
    real(kind_phys),      intent(inout) :: ddy_aer(:), ddx_aer(:)
    real(kind_phys),      intent(out)   :: aer_nm(:,:,:)
    integer,              intent(inout) :: jindx1_ci(:), jindx2_ci(:), iindx1_ci(:), iindx2_ci(:)
    real(kind_phys),      intent(inout) :: ddy_ci(:), ddx_ci(:)
    integer,              intent(inout) :: imap(:), jmap(:)
    logical,              intent(in)    :: do_ugwp_v1
    real(kind_phys),      intent(inout) :: ddy_j1tau(:), ddy_j2tau(:)
    integer,              intent(inout) :: jindx1_tau(:), jindx2_tau(:)
    integer,              intent(in)    :: isot, ivegsrc, nlunit
    real(kind_phys),      intent(inout) :: sncovr(:), sncovr_ice(:)
    integer,              intent(in)    :: lsm, lsm_noahmp, lsm_ruc, vtype(:)
    real(kind_phys),      intent(in)    :: min_seaice, fice(:)
    real(kind_phys),      intent(in)    :: landfrac(:)
    real(kind_phys),      intent(inout) :: weasd(:)
    type(ty_ozphys),      intent(in)    :: ozphys
    ! NoahMP - only allocated when NoahMP is used
    integer,              intent(in)    :: lsoil, lsnow_lsm_lbound, lsnow_lsm_ubound
    real(kind_phys),      intent(in)    :: zs(:)
    real(kind_phys),      intent(in)    :: dzs(:)
    real(kind_phys),      intent(inout) :: tvxy(:)
    real(kind_phys),      intent(inout) :: tgxy(:)
    real(kind_phys),      intent(inout) :: tahxy(:)
    real(kind_phys),      intent(inout) :: canicexy(:)
    real(kind_phys),      intent(inout) :: canliqxy(:)
    real(kind_phys),      intent(inout) :: eahxy(:)
    real(kind_phys),      intent(inout) :: cmxy(:)
    real(kind_phys),      intent(inout) :: chxy(:)
    real(kind_phys),      intent(inout) :: fwetxy(:)
    real(kind_phys),      intent(inout) :: sneqvoxy(:)
    real(kind_phys),      intent(inout) :: alboldxy(:)
    real(kind_phys),      intent(inout) :: qsnowxy(:)
    real(kind_phys),      intent(inout) :: wslakexy(:)
    real(kind_phys),      intent(inout) :: albdvis_lnd(:)
    real(kind_phys),      intent(inout) :: albdnir_lnd(:)
    real(kind_phys),      intent(inout) :: albivis_lnd(:)
    real(kind_phys),      intent(inout) :: albinir_lnd(:)
    real(kind_phys),      intent(inout) :: albdvis_ice(:)
    real(kind_phys),      intent(inout) :: albdnir_ice(:)
    real(kind_phys),      intent(inout) :: albivis_ice(:)
    real(kind_phys),      intent(inout) :: albinir_ice(:)
    real(kind_phys),      intent(inout) :: emiss_lnd(:)
    real(kind_phys),      intent(inout) :: emiss_ice(:)
    real(kind_phys),      intent(inout) :: taussxy(:)
    real(kind_phys),      intent(inout) :: waxy(:)
    real(kind_phys),      intent(inout) :: wtxy(:)
    real(kind_phys),      intent(inout) :: zwtxy(:)
    real(kind_phys),      intent(inout) :: xlaixy(:)
    real(kind_phys),      intent(inout) :: xsaixy(:)
    real(kind_phys),      intent(inout) :: lfmassxy(:)
    real(kind_phys),      intent(inout) :: stmassxy(:)
    real(kind_phys),      intent(inout) :: rtmassxy(:)
    real(kind_phys),      intent(inout) :: woodxy(:)
    real(kind_phys),      intent(inout) :: stblcpxy(:)
    real(kind_phys),      intent(inout) :: fastcpxy(:)
    real(kind_phys),      intent(inout) :: smcwtdxy(:)
    real(kind_phys),      intent(inout) :: deeprechxy(:)
    real(kind_phys),      intent(inout) :: rechxy(:)
    real(kind_phys),      intent(inout) :: snowxy(:)
    real(kind_phys),      intent(inout) :: snicexy(:,lsnow_lsm_lbound:)
    real(kind_phys),      intent(inout) :: snliqxy(:,lsnow_lsm_lbound:)
    real(kind_phys),      intent(inout) :: tsnoxy (:,lsnow_lsm_lbound:)
    real(kind_phys),      intent(inout) :: smoiseq(:,:)
    real(kind_phys),      intent(inout) :: zsnsoxy(:,lsnow_lsm_lbound:)
    real(kind_phys),      intent(inout) :: slc(:,:)
    real(kind_phys),      intent(inout) :: smc(:,:)
    real(kind_phys),      intent(inout) :: stc(:,:)
    real(kind_phys),      intent(in)    :: tsfcl(:)
    real(kind_phys),      intent(in)    :: snowd(:)
    real(kind_phys),      intent(in)    :: canopy(:)
    real(kind_phys),      intent(in)    :: tg3(:)
    integer,              intent(in)    :: stype(:)
    integer,              intent(in)    :: nthrds
    real (kind_phys),     intent(in)    :: si(:)
    integer,              intent(in)    :: levr, ictm, isol, ico2, iaer, ntcw, &
         num_p3d, ltp, npdf3d, iovr, iovr_rand, iovr_maxrand, iovr_max,        &
         iovr_dcorr, iovr_exp, iovr_exprand, icliq_sw, imp_physics,            &
         rad_hr_units, icliq_lw, isubcsw, isubclw, iswmode
    logical,              intent(in)    :: lcrick, lcnorm, lnoprec, do_RRTMGP, &
         lalw1bd, inc_minor_gas, lextop, iaerclm
    character(len=26),    intent(in)    :: aeros_file, solar_file, co2usr_file,&
         co2cyc_file
    real(kind_phys),      intent(in)    :: con_pi, con_t0c, con_c, con_boltz,  &
         con_plnk, con_solr_2008, con_solr_2002, con_g, con_rd
    integer,              intent(inout) :: ipsd0
    integer,              intent(out)   :: iaermdl, iaerflg
    character(len=*),     intent(out)   :: errmsg
    integer,              intent(out)   :: errflg

    ! Local variables
    integer :: i, j, ix, vegtyp
    real(kind_phys) :: rsnow
    integer         :: soiltyp, isnow, is, imn
    real(kind_phys) :: masslai, masssai, snd
    real(kind_phys) :: bexp, ddz, smcmax, smcwlt, dwsat, dksat, psisat
    real(kind_phys), dimension(:), allocatable :: dzsno
    real(kind_phys), dimension(:), allocatable :: dzsnso
    integer :: myerrflg
    character(len=255) :: myerrmsg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (is_initialized) return

    !> Call gfuncphys (funcphys.f) to compute all physics function tables.
    call gfuncphys ()

    !>
    !> Radiation setup
    !>
    if ( ictm==0 .or. ictm==-2 ) then
       iaerflg = mod(iaer, 100)        ! no volcanic aerosols for clim hindcast
    else
       iaerflg = mod(iaer, 1000)
    endif
    iaermdl = iaer/1000               ! control flag for aerosol scheme selection
    if ( iaermdl < 0 .or.  (iaermdl>2 .and. iaermdl/=5) ) then
       print *, ' Error -- IAER flag is incorrect, Abort'
       errflg = 1
       errmsg = 'ERROR(GFS_phys_time_vary): IAER flag is incorrect'
       return
    endif

    !>
    !> Assign initial permutation seed for mcica cloud-radiation.
    !>
    if ( isubcsw>0 .or. isubclw>0 ) then
       ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4)
    endif

    !>
    !> Call initialization routines for radiaiton modules.
    !>
    call sol_init ( me, isol, solar_file, con_solr_2008,con_solr_2002, &
         con_pi )
    call aer_init ( levr, me, iaermdl, iaerflg, lalw1bd, aeros_file,   &
         con_pi, con_t0c, con_c, con_boltz, con_plnk, errflg, errmsg)
    call gas_init ( me, co2usr_file, co2cyc_file, ico2, ictm, con_pi,  &
         errflg, errmsg )
    call cld_init ( si, levr, imp_physics, me, con_g, con_rd, errflg,  &
         errmsg)
    call rlwinit ( me, rad_hr_units, inc_minor_gas, icliq_lw, isubcsw, &
         iovr, iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr,          &
         iovr_exp, iovr_exprand, errflg, errmsg )
    call rswinit ( me, rad_hr_units, inc_minor_gas, icliq_sw, isubclw, &
         iovr, iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr,          &
         iovr_exp, iovr_exprand,iswmode, errflg, errmsg )
  
    if ( me == 0 ) then
       print *,' In radiation initialize (GFS_phys_time_vary)'
       print *,' si =',si
       print *,' levr=',levr,' ictm=',ictm,' isol=',isol,' ico2=',ico2,&
               ' iaermdl=',iaermdl,' iaerflg=',iaerflg
       print *,' np3d=',num_p3d,' ntoz=',ntoz,                         &
               ' iovr=',iovr,' isubcsw=',isubcsw,                      &
               ' isubclw=',isubclw,' icliq_sw=',icliq_sw,              &
               ' iflip=',iflip,'  me=',me
       print *,' lcrick=',lcrick,                                      &
               ' lcnorm=',lcnorm,' lnoprec=',lnoprec
       print *, 'lextop=',lextop, ' ltp=',ltp
       print *,' Radiation sub-cloud initial seed =',ipsd0,            &
               ' IC-idate =',idate
    endif

    !>
    !> - Call read_h2odata() to read stratospheric water vapor data
    !>
    need_h2odata: if(h2o_phys) then
       call read_h2odata (h2o_phys, me, master)

       ! Consistency check that the hardcoded values for levh2o and
       ! h2o_coeff in GFS_typedefs.F90 match what is set by read_h2odata
       ! in GFS_typedefs.F90: allocate (Tbd%h2opl (IM,levh2o,h2o_coeff))
       if (size(h2opl, dim=2).ne.levh2o) then
          write(myerrmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",     &
               "levh2o from read_h2odata does not match value in GFS_typedefs.F90: ", &
               levh2o, " /= ", size(h2opl, dim=2)
          myerrflg = 1
          call copy_error(myerrmsg, myerrflg, errmsg, errflg)
       end if
       if (size(h2opl, dim=3).ne.h2o_coeff) then
          write(myerrmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",       &
               "h2o_coeff from read_h2odata does not match value in GFS_typedefs.F90: ", &
               h2o_coeff, " /= ", size(h2opl, dim=3)
          myerrflg = 1
          call copy_error(myerrmsg, myerrflg, errmsg, errflg)
       end if
    endif need_h2odata

    !>
    !> - Call read_aerdata() to read aerosol climatology, Anning added coupled
    !>  added coupled gocart and radiation option to initializing aer_nm
    !>
    if (iaerclm) then
       ntrcaer = ntrcaerm
       myerrflg = 0
       myerrmsg = 'read_aerdata failed without a message'
       call read_aerdata (me,master,iflip,idate,myerrmsg,myerrflg)
       call copy_error(myerrmsg, myerrflg, errmsg, errflg)
    else if(iaermdl ==2 ) then
       do ix=1,ntrcaerm
          do j=1,levs
             do i=1,im
                aer_nm(i,j,ix) = 1.e-20_kind_phys
             end do
          end do
       end do
       ntrcaer = ntrcaerm
    else
       ntrcaer = 1
    endif

    !> - Call read_cidata() to read IN and CCN data
    if (iccn == 1) then
       call read_cidata (me,master)
       ! No consistency check needed for in/ccn data, all values are
       ! hardcoded in module iccn_def.F and GFS_typedefs.F90
    endif

    !> - Call tau_amf dats for  ugwp_v1
    if (do_ugwp_v1) then
       myerrflg = 0
       myerrmsg = 'read_tau_amf failed without a message'
       call read_tau_amf(me, master, myerrmsg, myerrflg)
       call copy_error(myerrmsg, myerrflg, errmsg, errflg)
    endif
    
    !> - Initialize soil vegetation (needed for sncovr calculation further down)
    myerrflg = 0
    myerrmsg = 'set_soilveg failed without a message'
    call set_soilveg(me, isot, ivegsrc, nlunit, myerrmsg, myerrflg)
    call copy_error(myerrmsg, myerrflg, errmsg, errflg)
    
    !> - read in NoahMP table (needed for NoahMP init)
    if(lsm == lsm_noahmp) then
       myerrflg = 0
       myerrmsg = 'read_mp_table_parameters failed without a message'
       call read_mp_table_parameters(myerrmsg, myerrflg)
       call copy_error(myerrmsg, myerrflg, errmsg, errflg)
    endif

    ! Need an OpenMP barrier here (implicit in "end sections")

    !> - Setup spatial interpolation indices for ozone physics.
    if (ntoz > 0) then
       call ozphys%setup_o3prog(xlat_d, jindx1_o3, jindx2_o3, ddy_o3)
    endif

    !> - Call setindxh2o() to initialize stratospheric water vapor data
    if (h2o_phys) then
       call setindxh2o (im, xlat_d, jindx1_h, jindx2_h, ddy_h)
    endif

    !> - Call setindxaer() to initialize aerosols data
    if (iaerclm) then
       call setindxaer (im, xlat_d, jindx1_aer, jindx2_aer, ddy_aer, xlon_d,     &
            iindx1_aer, iindx2_aer, ddx_aer,  me, master)
       iamin = 999
       iamax = -999
       jamin = 999
       jamax = -999
       iamin = min(minval(iindx1_aer), iamin)
       iamax = max(maxval(iindx2_aer), iamax)
       jamin = min(minval(jindx1_aer), jamin)
       jamax = max(maxval(jindx2_aer), jamax)
    endif

    !> - Call setindxci() to initialize IN and CCN data
    if (iccn == 1) then
       call setindxci (im, xlat_d, jindx1_ci,jindx2_ci, ddy_ci, xlon_d,  &
            iindx1_ci, iindx2_ci, ddx_ci)
    endif

    !> - Call  cires_indx_ugwp to read monthly-mean GW-tau diagnosed from FV3GFS-runs that can resolve GWs
    if (do_ugwp_v1) then
       call cires_indx_ugwp (im, me, master, xlat_d, jindx1_tau, jindx2_tau,  &
            ddy_j1tau, ddy_j2tau)
    endif

    !--- initial calculation of maps local ix -> global i and j
    ix = 0
    do j = 1,ny
       do i = 1,nx
          ix = ix + 1
          jmap(ix) = j
          imap(ix) = i
       enddo
    enddo
    
    !--- if sncovr does not exist in the restart, need to create it
    if (all(sncovr < zero)) then
       if (me == master ) write(*,'(a)') 'GFS_phys_time_vary_init: compute sncovr from weasd and soil vegetation parameters'
       !--- compute sncovr from existing variables
       !--- code taken directly from read_fix.f
       sncovr(:) = zero
       do ix=1,im
          if (landfrac(ix) >= drythresh .or. fice(ix) >= min_seaice) then
             vegtyp = vtype(ix)
             if (vegtyp == 0) vegtyp = 7
             rsnow  = 0.001_kind_phys*weasd(ix)/snupx(vegtyp)
             if (0.001_kind_phys*weasd(ix) < snupx(vegtyp)) then
                sncovr(ix) = one - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
             else
                sncovr(ix) = one
             endif
          endif
       enddo
    endif
    
    !--- For RUC LSM: create sncovr_ice from sncovr
    if (lsm == lsm_ruc) then
       if (all(sncovr_ice < zero)) then
          if (me == master ) write(*,'(a)') 'GFS_phys_time_vary_init: fill sncovr_ice with sncovr for RUC LSM'
          sncovr_ice(:) = sncovr(:)
       endif
    endif
    
    if (errflg/=0) return
    
    if (iaerclm) then
       ! This call is outside the OpenMP section, so it should access errmsg & errflg directly.
       call read_aerdataf (me, master, iflip, idate, fhour, errmsg, errflg)
       ! If it is moved to an OpenMP section, it must use myerrmsg, myerrflg, and copy_error.
       if (errflg/=0) return
    end if
    
    !--- For Noah MP or RUC LSMs: initialize four components of albedo for
    !--- land and ice - not for restart runs
    lsm_init: if (lsm_cold_start) then
       if (lsm == lsm_noahmp .or. lsm == lsm_ruc) then
          if (me == master ) write(*,'(a)') 'GFS_phys_time_vary_init: initialize albedo for land and ice'
          do ix=1,im
             albdvis_lnd(ix)  = 0.2_kind_phys
             albdnir_lnd(ix)  = 0.2_kind_phys
             albivis_lnd(ix)  = 0.2_kind_phys
             albinir_lnd(ix)  = 0.2_kind_phys
             emiss_lnd(ix)    = 0.95_kind_phys
          enddo
       endif
       if (lsm == lsm_ruc) then
          do ix=1,im
             albdvis_ice(ix)  = 0.6_kind_phys
             albdnir_ice(ix)  = 0.6_kind_phys
             albivis_ice(ix)  = 0.6_kind_phys
             albinir_ice(ix)  = 0.6_kind_phys
             emiss_ice(ix)    = 0.97_kind_phys
          enddo
       endif
       
       noahmp_init: if (lsm == lsm_noahmp) then
          allocate(dzsno (lsnow_lsm_lbound:lsnow_lsm_ubound))
          allocate(dzsnso(lsnow_lsm_lbound:lsoil)           )
          dzsno(:)    = missing_value
          dzsnso(:)   = missing_value
          
          tvxy(:)     = missing_value
          tgxy(:)     = missing_value
          tahxy(:)    = missing_value
          canicexy(:) = missing_value
          canliqxy(:) = missing_value
          eahxy(:)    = missing_value
          cmxy(:)     = missing_value
          chxy(:)     = missing_value
          fwetxy(:)   = missing_value
          sneqvoxy(:) = missing_value
          alboldxy(:) = missing_value
          qsnowxy(:)  = missing_value
          wslakexy(:) = missing_value
          taussxy(:)  = missing_value
          waxy(:)     = missing_value
          wtxy(:)     = missing_value
          zwtxy(:)    = missing_value
          xlaixy(:)   = missing_value
          xsaixy(:)   = missing_value
          
          lfmassxy(:)   = missing_value
          stmassxy(:)   = missing_value
          rtmassxy(:)   = missing_value
          woodxy(:)     = missing_value
          stblcpxy(:)   = missing_value
          fastcpxy(:)   = missing_value
          smcwtdxy(:)   = missing_value
          deeprechxy(:) = missing_value
          rechxy(:)     = missing_value

          snowxy (:)    = missing_value
          snicexy(:,:)  = missing_value
          snliqxy(:,:)  = missing_value
          tsnoxy (:,:)  = missing_value
          smoiseq(:,:)  = missing_value
          zsnsoxy(:,:)  = missing_value

          imn           = idate(2)

!$OMP parallel do num_threads(nthrds) default(none)                     &
!$OMP          shared(im,lsoil,con_t0c,landfrac,tsfcl,tvxy,tgxy,tahxy)  &
!$OMP          shared(snowd,canicexy,canliqxy,canopy,eahxy,cmxy,chxy)   &
!$OMP          shared(fwetxy,sneqvoxy,weasd,alboldxy,qsnowxy,wslakexy)  &
!$OMP          shared(taussxy)                                          &
!$OMP          shared(waxy,wtxy,zwtxy,imn,vtype,xlaixy,xsaixy,lfmassxy) &
!$OMP          shared(stmassxy,rtmassxy,woodxy,stblcpxy,fastcpxy)       &
!$OMP          shared(isbarren_table,isice_table,isurban_table)         &
!$omp          shared(iswater_table,laim_table,sla_table,bexp_table)    &
!$omp          shared(stc,smc,slc,tg3,snowxy,tsnoxy,snicexy,snliqxy)    &
!$omp          shared(zsnsoxy,stype,smcmax_table,smcwlt_table,zs,dzs)   & 
!$omp          shared(dwsat_table,dksat_table,psisat_table,smoiseq)     &
!$OMP          shared(smcwtdxy,deeprechxy,rechxy,errmsg,errflg)         &
!$OMP          private(vegtyp,masslai,masssai,snd,dzsno,dzsnso,isnow)   &
!$OMP          private(soiltyp,bexp,smcmax,smcwlt,dwsat,dksat,psisat)   &
!$OMP          private(myerrmsg,myerrflg,ddz)
          do ix=1,im
             if (landfrac(ix) >= drythresh) then
                tvxy(ix)     = tsfcl(ix)
                tgxy(ix)     = tsfcl(ix)
                tahxy(ix)    = tsfcl(ix)
                
                if (snowd(ix) > 0.01_kind_phys .and. tsfcl(ix) > con_t0c ) then
                   tvxy(ix)  = con_t0c
                   tgxy(ix)  = con_t0c
                   tahxy(ix) = con_t0c
                end if
                
                canicexy(ix) = 0.0_kind_phys
                canliqxy(ix) = canopy(ix)

                eahxy(ix)    = 2000.0_kind_phys

                cmxy(ix)     = zero
                chxy(ix)     = zero
                fwetxy(ix)   = zero
                sneqvoxy(ix) = weasd(ix)     ! mm
                alboldxy(ix) = 0.65_kind_phys
                qsnowxy(ix)  = zero

!                 if (srflag(ix) > 0.001) qsnowxy(ix) = tprcp(ix)/dtp
                ! already set to 0.0
                wslakexy(ix) = zero
                taussxy(ix)  = zero

                waxy(ix)     = 4900.0_kind_phys
                wtxy(ix)     = waxy(ix)
                zwtxy(ix)    = (25.0_kind_phys + 2.0_kind_phys) - waxy(ix) / 1000.0_kind_phys / 0.2_kind_phys

                vegtyp       = vtype(ix)
                if (vegtyp == 0) vegtyp = 7

                if ((vegtyp == isbarren_table) .or. (vegtyp == isice_table) .or. (vegtyp == isurban_table) .or. (vegtyp == iswater_table)) then
                   xlaixy(ix)   = zero
                   xsaixy(ix)   = zero
                   lfmassxy(ix) = zero
                   stmassxy(ix) = zero
                   rtmassxy(ix) = zero
                   woodxy   (ix) = zero
                   stblcpxy (ix) = zero
                   fastcpxy (ix) = zero
                else
                   xlaixy(ix)   = max(laim_table(vegtyp, imn),0.05_kind_phys)
                   xsaixy(ix)   = max(xlaixy(ix)*0.1_kind_phys,0.05_kind_phys)
                   masslai      = 1000.0_kind_phys / max(sla_table(vegtyp),one)
                   lfmassxy(ix) = xlaixy(ix)*masslai
                   masssai      = 1000.0_kind_phys / 3.0_kind_phys
                   stmassxy(ix) = xsaixy(ix)* masssai
                   rtmassxy(ix) = 500.0_kind_phys
                   woodxy(ix)   = 500.0_kind_phys
                   stblcpxy(ix) = 1000.0_kind_phys
                   fastcpxy(ix) = 1000.0_kind_phys
                endif  ! non urban ...

                if (vegtyp == isice_table) then
                   do is = 1,lsoil
                      stc(ix,is) = min(stc(ix,is),min(tg3(ix),263.15_kind_phys))
                      smc(ix,is) = one
                      slc(ix,is) = zero
                   enddo
                endif
                
                snd = snowd(ix)/1000.0_kind_phys  ! go to m from snwdph
                
                if (weasd(ix) /= zero .and. snd == zero ) then
                   snd = weasd(ix)/1000.0
                endif
                
                if (vegtyp == 15) then                      ! land ice in MODIS/IGBP
                   if (weasd(ix) < 0.1_kind_phys) then
                      weasd(ix) = 0.1_kind_phys
                      snd       = 0.01_kind_phys
                   endif
                endif
                
                if (snd < 0.025_kind_phys ) then
                   snowxy(ix)   = zero
                   dzsno(-2:0)  = zero
                elseif (snd >= 0.025_kind_phys .and. snd <= 0.05_kind_phys ) then
                   snowxy(ix)   = -1.0_kind_phys
                   dzsno(0)     = snd
                elseif (snd > 0.05_kind_phys .and. snd <= 0.10_kind_phys ) then
                   snowxy(ix)   = -2.0_kind_phys
                   dzsno(-1)    = 0.5_kind_phys*snd
                   dzsno(0)     = 0.5_kind_phys*snd
                elseif (snd > 0.10_kind_phys .and. snd <= 0.25_kind_phys ) then
                   snowxy(ix)   = -2.0_kind_phys
                   dzsno(-1)    = 0.05_kind_phys
                   dzsno(0)     = snd - 0.05_kind_phys
                elseif (snd > 0.25_kind_phys .and. snd <= 0.45_kind_phys ) then
                   snowxy(ix)   = -3.0_kind_phys
                   dzsno(-2)    = 0.05_kind_phys
                   dzsno(-1)    = 0.5_kind_phys*(snd-0.05_kind_phys)
                   dzsno(0)     = 0.5_kind_phys*(snd-0.05_kind_phys)
                elseif (snd > 0.45_kind_phys) then
                   snowxy(ix)   = -3.0_kind_phys
                   dzsno(-2)    = 0.05_kind_phys
                   dzsno(-1)    = 0.20_kind_phys
                   dzsno(0)     = snd - 0.05_kind_phys - 0.20_kind_phys
                else
                   myerrmsg = 'Error in GFS_phys_time_vary.fv3.F90: Problem with the logic assigning snow layers in Noah MP initialization'
                   myerrflg = 1
                   call copy_error(myerrmsg, myerrflg, errmsg, errflg)
                endif
                
                ! Now we have the snowxy field
                ! snice + snliq + tsno allocation and compute them from what we have
                
                tsnoxy(ix,:)  = zero
                snicexy(ix,:) = zero
                snliqxy(ix,:) = zero
                zsnsoxy(ix,:) = zero
                
                isnow = nint(snowxy(ix))+1 ! snowxy <=0.0, dzsno >= 0.0
                
                do is = isnow,0
                   tsnoxy(ix,is)  = tgxy(ix)
                   snliqxy(ix,is) = zero
                   snicexy(ix,is) = one * dzsno(is) * weasd(ix)/snd
                enddo
                !
                !zsnsoxy, all negative ?
                !
                do is = isnow,0
                   dzsnso(is) = -dzsno(is)
                enddo
                
                do is = 1,4
                   dzsnso(is) = -dzs(is)
                enddo
                !
                ! Assign to zsnsoxy
                !
                zsnsoxy(ix,isnow) = dzsnso(isnow)
                do is = isnow+1,4
                   zsnsoxy(ix,is) = zsnsoxy(ix,is-1) + dzsnso(is)
                enddo
                !
                ! smoiseq
                ! Init water table related quantities here
                !
                soiltyp  = stype(ix)
                if (soiltyp /= 0) then
                   bexp   = bexp_table(soiltyp)
                   smcmax = smcmax_table(soiltyp)
                   smcwlt = smcwlt_table(soiltyp)
                   dwsat  = dwsat_table(soiltyp)
                   dksat  = dksat_table(soiltyp)
                   psisat = -psisat_table(soiltyp)
                endif
                
                if (vegtyp == isurban_table) then
                   smcmax = 0.45_kind_phys
                   smcwlt = 0.40_kind_phys
                endif
                
                if ((bexp > zero) .and. (smcmax > zero) .and. (-psisat > zero)) then
                   do is = 1, lsoil
                      if ( is == 1 )then
                         ddz = -zs(is+1) * 0.5_kind_phys
                      elseif ( is < lsoil ) then
                         ddz = ( zs(is-1) - zs(is+1) ) * 0.5_kind_phys
                      else
                         ddz = zs(is-1) - zs(is)
                      endif
                      smoiseq(ix,is) = min(max(find_eq_smc(bexp, dwsat, dksat, ddz, smcmax),1.e-4_kind_phys),smcmax*0.99_kind_phys)
                   enddo
                else                                    ! bexp <= 0.0
                   smoiseq(ix,1:4) = smcmax
                endif                                   ! end the bexp condition
                
                smcwtdxy(ix)   = smcmax
                deeprechxy(ix) = zero
                rechxy(ix)     = zero
                
             endif
             
          enddo ! ix
!$OMP end parallel do

          if (errflg/=0) return
          
          deallocate(dzsno)
          deallocate(dzsnso)

       endif noahmp_init
    endif lsm_init

    !
    ! Lake model
    !
    if(lkm>0 .and. iopt_lake>0) then
       ! A lake model is enabled.
       do i = 1, im
          !if (lakefrac(i) > 0.0 .and. lakedepth(i) > 1.0 ) then

          ! The lake data must say there's a lake here (lakefrac) with a depth (lakedepth)
          if (lakefrac(i) > lakefrac_threshold .and. lakedepth(i) > lakedepth_threshold ) then
             ! This is a lake point. Inform the other schemes to use a lake model, and possibly nsst (lkm)
             use_lake_model(i) = lkm
             cycle
          else
             ! Not a valid lake point.
             use_lake_model(i) = 0
          endif
       enddo
    else
       ! Lake model is disabled or settings are invalid.
       use_lake_model = 0
    endif

    is_initialized = .true.
    
  contains
    !
    ! Use newton-raphson method to find eq soil moisture
    !
    function find_eq_smc(bexp, dwsat, dksat, ddz, smcmax) result(smc)
      implicit none
      real(kind=kind_phys), intent(in) :: bexp, dwsat, dksat, ddz, smcmax
      real(kind=kind_phys) :: smc
      real(kind=kind_phys) :: expon, aa, bb, func, dfunc, dx
      integer :: iter
      !
      expon = bexp + 1.
      aa    = dwsat / ddz
      bb    = dksat / smcmax ** expon
      smc = 0.5 * smcmax
      !
      do iter = 1,100
         func  = (smc - smcmax) * aa +  bb * smc ** expon
         dfunc = aa + bb * expon * smc ** bexp
         dx    = func / dfunc
         smc   = smc - dx
         if ( abs (dx) < 1.e-6_kind_phys) return
      enddo
    end function find_eq_smc

  end subroutine GFS_phys_time_vary_init
!> @}

  ! #########################################################################################
  ! SUBROUTINE GFS_phys_time_vary_timestep_init
  ! #########################################################################################
!> \section arg_table_GFS_phys_time_vary_timestep_init Argument Table
!! \htmlinclude GFS_phys_time_vary_timestep_init.html
!!
!>\section gen_GFS_phys_time_vary_timestep_init GFS_phys_time_vary_timestep_init General Algorithm
!> @{
  subroutine GFS_phys_time_vary_timestep_init (me, master, cnx, cny, isc, jsc, nrcm, im,    &
       levs, nsswr, nslwr, fhswr, imfdeepcnv, cal_pre, random_clds, nscyc, ntoz,            &
       h2o_phys, iaerclm, iccn, clstp, jindx1_o3, jindx2_o3, ddy_o3, ozpl, jindx1_h,        &
       jindx2_h, ddy_h, h2opl, iflip, jindx1_aer, jindx2_aer, ddy_aer, iindx1_aer,          &
       iindx2_aer, ddx_aer, aer_nm, jindx1_ci, jindx2_ci, ddy_ci, iindx1_ci, iindx2_ci,     &
       ddx_ci, in_nm, ccn_nm, fn_nml, imap, jmap, prsl, seed0, rann, nthrds, nx, ny, nsst,  &
       tile_num, nlunit, lsoil, lsoil_lsm, kice, ialb, isot, ivegsrc, input_nml_file,       &
       use_ufo, nst_anl, frac_grid, fhcyc, lakefrac, min_seaice, min_lakeice, smc, slc,     &
       stc, smois, sh2o, tslb, tiice, tg3, tref, tsfc, tsfco, tisfc, hice, fice, facsf,     &
       facwf, alvsf, alvwf, alnsf, alnwf, zorli, zorll, zorlo, weasd, slope, snoalb, canopy,&
       vfrac, vtype, stype,scolor, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf, xlat_d,&
       xlon_d, slmsk, landfrac, ozphys, do_ugwp_v1, jindx1_tau, jindx2_tau, ddy_j1tau,      &
       ddy_j2tau, tau_amf, lrseeds, rseeds, isubc_lw, isubc_sw, icsdsw, icsdlw, imp_physics,&
       imp_physics_zhao_carr, ipsd0, ipsdlim, ps_2delt, ps_1delt, t_2delt, t_1delt,         &
       qv_2delt, qv_1delt, t, qv, ps, jdat, idat, nhfrad, debug, dtp, kdt, yearlen, ipt,    &
       lprnt, lssav, lslwr, lsswr, sec, phour, zhour, fhour, julian, solhr, deltim,         &
       iaermdl, aeros_file, isol, slag, sdec, cdec, solcon, con_pi, co2dat_file,            &
       co2gbl_file, ictm, ico2, errmsg, errflg)
    
    implicit none

    ! Interface variables
    integer,              intent(in)    :: me, master, cnx, cny, isc, jsc, nrcm, im, levs
    integer,              intent(in)    :: nsswr, nslwr, imfdeepcnv, iccn, nscyc, ntoz, iflip
    real(kind_phys),      intent(in)    :: fhswr
    logical,              intent(in)    :: cal_pre, random_clds, h2o_phys, iaerclm
    real(kind_phys),      intent(out)   :: clstp
    integer,              intent(in)    :: jindx1_o3(:), jindx2_o3(:), jindx1_h(:), jindx2_h(:)
    real(kind_phys),      intent(in)    :: ddy_o3(:),  ddy_h(:)
    real(kind_phys),      intent(inout) :: ozpl(:,:,:), h2opl(:,:,:)
    integer,              intent(in)    :: jindx1_aer(:), jindx2_aer(:), iindx1_aer(:), iindx2_aer(:)
    real(kind_phys),      intent(in)    :: ddy_aer(:), ddx_aer(:)
    real(kind_phys),      intent(inout) :: aer_nm(:,:,:)
    integer,              intent(in)    :: jindx1_ci(:), jindx2_ci(:), iindx1_ci(:), iindx2_ci(:)
    real(kind_phys),      intent(in)    :: ddy_ci(:), ddx_ci(:)
    real(kind_phys),      intent(inout) :: in_nm(:,:), ccn_nm(:,:)
    integer,              intent(in)    :: imap(:), jmap(:)
    real(kind_phys),      intent(in)    :: prsl(:,:)
    integer,              intent(in)    :: seed0
    real(kind_phys),      intent(inout) :: rann(:,:)
    logical,              intent(in)    :: do_ugwp_v1
    integer,              intent(in)    :: jindx1_tau(:), jindx2_tau(:)
    real(kind_phys),      intent(in)    :: ddy_j1tau(:), ddy_j2tau(:)
    real(kind_phys),      intent(inout) :: tau_amf(:)
    type(ty_ozphys),      intent(inout) :: ozphys
    ! For gcycle only
    integer,              intent(in)    :: nthrds, nx, ny, nsst, tile_num, nlunit, lsoil
    integer,              intent(in)    :: lsoil_lsm, kice, ialb, isot, ivegsrc
    character(len=*),     intent(in)    :: input_nml_file(:)
    character(len=*),     intent(in)    :: fn_nml
    logical,              intent(in)    :: use_ufo, nst_anl, frac_grid
    real(kind_phys),      intent(in)    :: fhcyc, lakefrac(:), min_seaice, min_lakeice,  &
                                           xlat_d(:), xlon_d(:), landfrac(:)
    real(kind_phys),      intent(inout) :: smc(:,:), slc(:,:), stc(:,:), smois(:,:), sh2o(:,:), &
                                           tslb(:,:), tiice(:,:), tg3(:), tref(:),                        &
                                           tsfc(:), tsfco(:), tisfc(:), hice(:), fice(:),                 &
                                           facsf(:), facwf(:), alvsf(:), alvwf(:), alnsf(:), alnwf(:),    &
                                           zorli(:), zorll(:), zorlo(:), weasd(:), snoalb(:),             &
                                           canopy(:), vfrac(:), shdmin(:), shdmax(:),                     &
                                           snowd(:), cv(:), cvb(:), cvt(:), oro(:), oro_uf(:), slmsk(:)
    integer,              intent(inout) :: vtype(:), stype(:),scolor(:), slope(:) 
    ! For radiation
    logical,              intent(in)    :: lrseeds
    integer,              intent(in)    :: rseeds(:,:)
    integer,              intent(in)    :: isubc_lw, isubc_sw
    integer,              intent(in)    :: imp_physics, imp_physics_zhao_carr, ipsd0, ipsdlim
    logical,              intent(out)   :: lslwr, lsswr
    integer,              intent(inout) :: icsdsw(:), icsdlw(:)
    real(kind_phys),      intent(out)   :: sec
    real(kind_phys),      intent(in)    :: deltim
    real(kind_phys),      intent(in)    :: con_pi
    integer,              intent(in)    :: iaermdl, isol, ictm, ico2
    character(len=26),    intent(in)    :: aeros_file, co2dat_file, co2gbl_file
    real(kind_phys),      intent(out)   :: slag
    real(kind_phys),      intent(out)   :: sdec
    real(kind_phys),      intent(out)   :: cdec
    real(kind_phys),      intent(out)   :: solcon
    ! For Zhao-Carr Microphysics
    real(kind_phys),      intent(inout) :: ps_2delt(:)
    real(kind_phys),      intent(inout) :: ps_1delt(:)
    real(kind_phys),      intent(inout) :: t_2delt(:,:)
    real(kind_phys),      intent(inout) :: t_1delt(:,:)
    real(kind_phys),      intent(inout) :: qv_2delt(:,:)
    real(kind_phys),      intent(inout) :: qv_1delt(:,:)
    real(kind_phys),      intent(in)    :: t(:,:), qv(:,:), ps(:)
    ! For GFS calendar
    integer,              intent(in)    :: jdat(:), idat(:)
    integer,              intent(in)    :: nhfrad
    logical,              intent(in)    :: debug
    real(kind_phys),      intent(in)    :: dtp
    real(kind_phys),      intent(out)   :: phour, zhour, fhour, julian, solhr
    integer,              intent(out)   :: kdt, yearlen, ipt
    logical,              intent(out)   :: lprnt, lssav
    !
    character(len=*),     intent(out)   :: errmsg
    integer,              intent(out)   :: errflg

    ! Local variables
    integer :: i, j, k, iseed, iskip, ix, idat_local(8), jdat_local(8), j1, j2, &
         nc, n1, n2, jdow, jdoy, jday, w3kindreal, w3kindint
    real(kind_phys) :: wrk(1), tem, tx1, tx2, rjday
    real(kind_phys) :: rannie(cny)
    real(kind_phys) :: rndval(cnx*cny*nrcm)
    real(kind_dbl_prec)  :: rinc8(5)
    real(kind_sngl_prec) :: rinc4(5)
    type (random_stat) :: stat
    integer :: ipseed
    integer :: numrdm(cnx*cny*2)
    integer :: iw3jdn
    integer :: jd0, jd1
    real    :: fjd
    integer :: iyear, imon, iday, ihour
    integer :: kyear, kmon, kday, khour
    logical :: lmon_chg       ! month change flag
    logical :: lco2_chg       ! cntrl flag for updating co2 data
    logical :: lsol_chg       ! cntrl flag for updating solar constant

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Check initialization status
    if (.not.is_initialized) then
       write(errmsg,'(*(a))') "Logic error: GFS_phys_time_vary_timestep_init called before GFS_phys_time_vary_init"
       errflg = 1
       return
    end if

    !
    ! Set up time stamp at fcst time and that for green house gases.
    !
    iyear = jdat(1)
    imon  = jdat(2)
    iday  = jdat(3)
    ihour = jdat(5)

    ! Set up time stamp used for green house gases (** currently co2 only)
    ! get external data at initial condition time
    if ( ictm==0 .or. ictm==-2 ) then
       kyear = idat(1)
       kmon  = idat(2)
       kday  = idat(3)
       khour = idat(5)
    ! get external data at fcst or specified time
    else
       kyear = iyear
       kmon  = imon
       kday  = iday
       khour = ihour
    endif

    if ( month0 /= imon ) then
       lmon_chg = .true.
       month0   = imon
    else
       lmon_chg = .false.
    endif

    !
    ! Update solar forcing...
    !
    if (lsswr) then
       if ( isol == 0 .or. isol == 10 ) then
          lsol_chg = .false.
       elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
       else
          lsol_chg = ( isol==4 .and. lmon_chg )
       endif
       iyear0 = iyear
       call sol_update(jdat, kyear, fhswr, deltim, lsol_chg, me, slag, sdec, cdec, solcon, con_pi, errmsg, errflg)
    endif

    !
    ! Update aerosols...
    !
    if ( lmon_chg ) then
       call aer_update ( iyear, imon, me, iaermdl, aeros_file, errflg, errmsg)
    endif

    !
    ! Update trace gases (co2 only)...
    !
    if ( monthd /= kmon ) then
       monthd = kmon
       lco2_chg = .true.
    else
       lco2_chg = .false.
    endif
    call gas_update (kyear, kmon, kday, khour, lco2_chg, me, co2dat_file, co2gbl_file, ictm,&
         ico2, errflg, errmsg )
    
    !
    ! Update ozone concentration...
    !
    if (ntoz == 0) then
       call ozphys%update_o3clim(kmon, kday, khour, loz1st)
    endif
    if ( loz1st ) loz1st = .false.

!$OMP parallel num_threads(nthrds) default(none)                                         &
!$OMP          shared(kdt,nsswr,nslwr,lsswr,lslwr,clstp,imfdeepcnv,cal_pre,random_clds)  &
!$OMP          shared(fhswr,fhour,seed0,cnx,cny,nrcm,wrk,rannie,rndval)                  &
!$OMP          shared(rann,im,isc,jsc,imap,jmap,ntoz,me,jindx1_o3,jindx2_o3)             &
!$OMP          shared(ozpl,ddy_o3,h2o_phys,jindx1_h,jindx2_h,h2opl,ddy_h,iaerclm,master) &
!$OMP          shared(levs,prsl,iccn,jindx1_ci,jindx2_ci,ddy_ci,iindx1_ci,iindx2_ci)     &
!$OMP          shared(ddx_ci,in_nm,ccn_nm,do_ugwp_v1,jindx1_tau,jindx2_tau,ddy_j1tau)    &
!$OMP          shared(ddy_j2tau,tau_amf,iflip,ozphys,rjday,n1,n2,idat_local,jdat_local)  &
!$OMP          shared(rinc8,rinc4,w3kindreal,w3kindint,jdow,jdoy,jday,jdat,idat,sec)     &
!$OMP          shared(zhour,phour,yearlen,julian,solhr,ipt,lprnt,lssav,dtp,nhfrad)&
!$OMP          private(iseed,iskip,i,j,k,jd0,jd1,fjd)

!$OMP sections

!$OMP section
    !
    ! Update GFS calendars.
    !
    call w3kind(w3kindreal,w3kindint)
    if (w3kindreal == 8) then
       rinc8(1:5) = 0
       call w3difdat(jdat,idat,4,rinc8)
       sec = rinc8(4)
    else if (w3kindreal == 4) then
       rinc4(1:5) = 0
       call w3difdat(jdat,idat,4,rinc4)
       sec = rinc4(4)
    else
       write(0,*)' FATAL ERROR: Invalid w3kindreal'
       call abort
    endif
    phour = sec/con_hr
    
    ! Set current bucket hour
    zhour = phour
    fhour = (sec + dtp)/con_hr
    kdt   = nint((sec + dtp)/dtp)
    
    !GJF* These calculations were originally in GFS_physics_driver.F90 for
    !     NoahMP. They were moved to this routine since they only depend
    !     on time (not space). Note that this code is included as-is from
    !     GFS_physics_driver.F90, but it may be simplified by using more
    !     NCEP W3 library calls (e.g., see W3DOXDAT, W3FS13 for Julian day
    !     of year and W3DIFDAT to determine the integer number of days in
    !     a given year). *GJF
    ! Julian day calculation (fcst day of the year)
    ! we need yearln and julian to
    ! pass to noah mp sflx, idat is init, jdat is fcst;idat = jdat when kdt=1
    ! jdat is changing
    !
    jd1    = iw3jdn(jdat(1),jdat(2),jdat(3))
    jd0    = iw3jdn(jdat(1),1,1)
    fjd    = float(jdat(5))/24.0 + float(jdat(6))/1440.0
    julian = float(jd1-jd0) + fjd
    
    !
    ! Year length
    !
    ! what if the integration goes from one year to another?
    ! iyr or jyr ? from 365 to 366 or from 366 to 365
    !
    ! is this against model's noleap yr assumption?
    if (mod(jdat(1),4) == 0) then
       yearlen = 366
       if (mod(jdat(1),100) == 0) then
          yearlen = 365
          if (mod(jdat(1),400) == 0) then
             yearlen = 366
          endif
       endif
    endif

    ipt    = 1
    lprnt  = .false.
    lssav  = .true.

    !
    ! Setup radiation triggers
    !
    lsswr  = (mod(kdt, nsswr) == 1)
    lslwr  = (mod(kdt, nslwr) == 1)
    ! allow for radiation to be called on every physics time step, if needed
    if (nsswr == 1)  lsswr = .true.
    if (nslwr == 1)  lslwr = .true.
    ! allow for radiation to be called on every physics time step
    ! for the first nhfrad timesteps (for spinup, coldstarts only)
    if (kdt <= nhfrad) then
       lsswr = .true.
       lslwr = .true.
    end if

    !
    ! Set the solar hour based on a combination of phour and time initial hour
    !
    solhr  = mod(phour+idat(1),con_24)

!$OMP section

    !--- switch for saving convective clouds - cnvc90.f
    !--- aka Ken Campana/Yu-Tai Hou legacy
    if ((mod(kdt,nsswr) == 0) .and. (lsswr)) then
       !--- initialize,accumulate,convert
       clstp = 1100 + min(fhswr/con_hr,fhour,con_99)
    elseif (mod(kdt,nsswr) == 0) then
       !--- accumulate,convert
       clstp = 0100 + min(fhswr/con_hr,fhour,con_99)
    elseif (lsswr) then
       !--- initialize,accumulate
       clstp = 1100
    else
       !--- accumulate
       clstp = 0100
    endif

!$OMP section

    !--- random number needed for RAS and old SAS and when cal_pre=.true.
    !    imfdeepcnv < 0 when ras = .true.
    if ( (imfdeepcnv <= 0 .or. cal_pre) .and. random_clds ) then
       iseed = mod(con_100*sqrt(fhour*con_hr),1.0d9) + seed0
       call random_setseed(iseed)
       call random_number(wrk)
       do i = 1,cnx*nrcm
          iseed = iseed + nint(wrk(1)*1000.0) * i
          call random_setseed(iseed)
          call random_number(rannie)
          rndval(1+(i-1)*cny:i*cny) = rannie(1:cny)
       enddo

       do k = 1,nrcm
          iskip = (k-1)*cnx*cny
          do ix=1,im
             j = jmap(ix)
             i = imap(ix)
             rann(ix,k) = rndval(i+isc-1 + (j+jsc-2)*cnx + iskip)
          enddo
       enddo
    endif  ! imfdeepcnv, cal_re, random_clds

!$OMP section
    !> - Compute temporal interpolation indices for updating gas concentrations.
    idat_local=0
    idat_local(1)=idat(4)
    idat_local(2)=idat(2)
    idat_local(3)=idat(3)
    idat_local(5)=idat(1)
    rinc8=0.
    rinc8(2)=fhour
    call w3kind(w3kindreal,w3kindint)
    if(w3kindreal==4) then
       rinc4=rinc8
       CALL w3movdat(rinc4,idat_local,jdat_local)
    else
       CALL w3movdat(rinc8,idat_local,jdat_local)
    endif
    jdow = 0
    jdoy = 0
    jday = 0
    call w3doxdat(jdat_local,jdow,jdoy,jday)
    rjday = jdoy + jdat_local(5) / 24.
    if (rjday < ozphys%time(1)) rjday = rjday + 365.
    
    n2 = ozphys%ntime + 1
    do j=2,ozphys%ntime
       if (rjday < ozphys%time(j)) then
          n2 = j
          exit
       endif
    enddo
    n1 = n2 - 1
    if (n2 > ozphys%ntime) n2 = n2 - ozphys%ntime

    !
    !> - Update ozone concentration.
    !
    if (ntoz > 0) then
       call ozphys%update_o3prog(jindx1_o3, jindx2_o3, ddy_o3, rjday, n1, n2, ozpl)
    endif

!$OMP section
    !
    !> - Update stratospheric water vapor data.
    !
    if (h2o_phys) then
       call h2ointerpol (me, im, idat, fhour, jindx1_h, jindx2_h, h2opl, ddy_h)
    endif

!$OMP section
    !
    !> - Call ciinterpol() to make IN and CCN data interpolation
    !
    if (iccn == 1) then
       call ciinterpol (me, im, idat, fhour, jindx1_ci, jindx2_ci, ddy_ci, iindx1_ci,&
            iindx2_ci, ddx_ci, levs, prsl, in_nm, ccn_nm)
    endif

!$OMP section
    !
    !> - Call  cires_indx_ugwp to read monthly-mean GW-tau diagnosed from FV3GFS-runs that resolve GW-activ
    !
    if (do_ugwp_v1) then
       call tau_amf_interp(me, master, im, idat, fhour, jindx1_tau, jindx2_tau,       &
            ddy_j1tau, ddy_j2tau, tau_amf)
    endif

!$OMP end sections
!$OMP end parallel

    !
    !> - Call aerinterpol() to make aerosol interpolation
    !
    if (iaerclm) then
       ! aerinterpol is using threading inside, don't
       ! move into OpenMP parallel section above
       call aerinterpol (me, master, nthrds, im, idat, fhour, iflip, jindx1_aer, jindx2_aer, &
            ddy_aer, iindx1_aer, iindx2_aer, ddx_aer, levs, prsl, aer_nm, errmsg, errflg)
       if(errflg /= 0) then
          return
       endif
    endif

    !
    !> - Call gcycle() to repopulate specific time-varying surface properties for AMIP/forecast runs
    !
    if (nscyc >  0) then
       if (mod(kdt,nscyc) == 1) THEN
          call gcycle (me, nthrds, nx, ny, isc, jsc, nsst, tile_num, nlunit, fn_nml,       &
               input_nml_file, lsoil, lsoil_lsm, kice, idat, ialb, isot, ivegsrc,          &
               use_ufo, nst_anl, fhcyc, phour, landfrac, lakefrac, min_seaice, min_lakeice,&
               frac_grid, smc, slc, stc, smois, sh2o, tslb, tiice, tg3, tref, tsfc,        &
               tsfco, tisfc, hice, fice, facsf, facwf, alvsf, alvwf, alnsf, alnwf,         &
               zorli, zorll, zorlo, weasd, slope, snoalb, canopy, vfrac, vtype,            &
               stype, scolor, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf,            &
               xlat_d, xlon_d, slmsk, imap, jmap, errmsg, errflg)
       endif
    endif

    !
    !> - Set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
    !
    if (lsswr .or. lslwr) then
       if ((isubc_lw==2) .or. (isubc_sw==2)) then
          !NRL If random seeds supplied by NEPTUNE
          if(lrseeds) then
             do ix=1,size(jmap)
                icsdsw(ix) = rseeds(ix,1)
                icsdlw(ix) = rseeds(ix,2)
             enddo
          else
             ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
             call random_setseed (ipseed, stat)
             call random_index (ipsdlim, numrdm, stat)
             do ix=1,size(jmap)
                j = jmap(ix)
                i = imap(ix)
                !--- for testing purposes, replace numrdm with '100'
                icsdsw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx)
                icsdlw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx + cnx*cny)
             enddo
          end if !lrseeds
       endif  ! isubc_lw and isubc_sw
    endif

    !
    !> - Store state needed by Zhao-Carr Microphysics
    !
    if (imp_physics == imp_physics_zhao_carr) then
       if (kdt == 1) then
          t_2delt  = t
          t_1delt  = t
          qv_2delt = qv
          qv_1delt = qv
          ps_2delt = ps
          ps_1delt = ps
       endif
    endif
  end subroutine GFS_phys_time_vary_timestep_init
!> @}

  ! #########################################################################################
  ! SUBROUTINE GFS_phys_time_vary_timestep_finalize
  ! #########################################################################################
!> \section arg_table_GFS_phys_time_vary_timestep_finalize Argument Table
!! \htmlinclude GFS_phys_time_vary_timestep_finalize.html
!!
!>\section gen_GFS_phys_time_vary_timestep_finalize GFS_phys_time_vary_timestep_finalize General Algorithm
!> @{
  subroutine GFS_phys_time_vary_timestep_finalize (errmsg, errflg)

    implicit none

    ! Interface variables
    character(len=*),                 intent(out)   :: errmsg
    integer,                          intent(out)   :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine GFS_phys_time_vary_timestep_finalize
!> @}

  ! #########################################################################################
  ! SUBROUTINE GFS_phys_time_vary_finalize
  ! #########################################################################################
!> \section arg_table_GFS_phys_time_vary_finalize Argument Table
!! \htmlinclude GFS_phys_time_vary_finalize.html
!!
  subroutine GFS_phys_time_vary_finalize(errmsg, errflg)

    implicit none

    ! Interface variables
    character(len=*),                 intent(out)   :: errmsg
    integer,                          intent(out)   :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    ! Deallocate h2o arrays
    if (allocated(h2o_lat) ) deallocate(h2o_lat)
    if (allocated(h2o_pres)) deallocate(h2o_pres)
    if (allocated(h2o_time)) deallocate(h2o_time)
    if (allocated(h2oplin) ) deallocate(h2oplin)

    ! Deallocate aerosol arrays
    if (allocated(aerin)   ) deallocate(aerin)
    if (allocated(aer_pres)) deallocate(aer_pres)

    ! Deallocate IN and CCN arrays
    if (allocated(ciplin)  ) deallocate(ciplin)
    if (allocated(ccnin)   ) deallocate(ccnin)
    if (allocated(ci_pres) ) deallocate(ci_pres)

    ! Deallocate UGWP-input arrays
    if (allocated(ugwp_taulat)) deallocate(ugwp_taulat)
    if (allocated(tau_limb   )) deallocate(tau_limb)
    if (allocated(days_limb  )) deallocate(days_limb)

    ! DH* this is the place to deallocate whatever is allocated by gfuncphys() in GFS_time_vary_pre_init

    is_initialized = .false.

  end subroutine GFS_phys_time_vary_finalize
  subroutine copy_error(myerrmsg, myerrflg, errmsg, errflg)
    implicit none
    character(*), intent(in) :: myerrmsg
    integer, intent(in) :: myerrflg
    character(*), intent(out) :: errmsg
    integer, intent(inout) :: errflg
    if(myerrflg /= 0 .and. errflg == 0) then
       !$OMP CRITICAL
       if(errflg == 0) then
          errmsg = myerrmsg
          errflg = myerrflg
       endif
       !$OMP END CRITICAL
    endif
  end subroutine copy_error

end module GFS_phys_time_vary
