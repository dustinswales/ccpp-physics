!> \file cosp_simulator.F90
  ! #########################################################################################
module cosp_simulator
  use machine,              only: kind_phys
  use quickbeam,            only: radar_cfg
  use mod_quickbeam_optics, only: size_distribution, hydro_class_init, quickbeam_optics_init
  use mod_cosp_config,      only: pres_binCenters, pres_binEdges, tau_binCenters,           &
       tau_binEdges, cloudsat_binCenters, cloudsat_binEdges, calipso_binCenters,            &
       calipso_binEdges, misr_histHgtCenters, misr_histHgtEdges,  PARASOL_SZA, R_UNDEF,     &
       PARASOL_NREFL, LIDAR_NCAT,SR_BINS, N_HYDRO, RTTOV_MAX_CHANNELS, numMISRHgtBins,      &
       CLOUDSAT_DBZE_BINS, LIDAR_NTEMP, calipso_histBsct, numMODISTauBins, numMODISPresBins,&
       numMODISReffIceBins, numMODISReffLiqBins, numISCCPTauBins, numISCCPPresBins,         &
       numMISRTauBins, reffICE_binEdges, reffICE_binCenters, reffLIQ_binEdges,              &
       reffLIQ_binCenters, LIDAR_NTYPE, nCloudsatPrecipClass, nsza_cosp => PARASOL_NREFL,   &
       nprs_cosp => npres, ntau_cosp => ntau, ntau_cosp_modis => ntau, nsr_cosp => SR_BINS, &
       nhtmisr_cosp => numMISRHgtBins, nhydro => N_HYDRO, cloudsat_preclvl, vgrid_zl,       &
       vgrid_zu, vgrid_z
    use mod_cosp_stats,     only: cosp_change_vertical_grid
    use module_cosp,        only: subsample_and_optics_CAM6, subsample_and_optics_UFS
  implicit none

  ! #########################################################################################
  ! module types
  ! #########################################################################################
  type(radar_cfg) :: &
       csat_rcfg          !< Radar configuration (Cloudsat)
  type(size_distribution) :: &
       csat_sd            !< Size distribution used by radar simulator
  character(len=64) :: &
       csat_micro_scheme  !<

contains

!! \section arg_table_cosp_simulator_init
!! \htmlinclude cosp_simulator_init.html
!!
  ! #########################################################################################
  ! 
  ! #########################################################################################
  subroutine cosp_simulator_init(mpirank, mpiroot, do_cosp, do_isccp, do_misr, do_modis,    &
       do_cloudsat, do_calipso, do_grLidar532, do_atlid, do_parasol, cosp_nsubcol,          &
       imp_physics, imp_physics_thompson, imp_physics_gfdl, isccp_topht, isccp_topht_dir,   &
       errmsg, errflg)
    USE mod_cosp_modis_interface,      ONLY: cosp_modis_init
    USE mod_cosp_misr_interface,       ONLY: cosp_misr_init
    USE mod_cosp_isccp_interface,      ONLY: cosp_isccp_init
    USE mod_cosp_calipso_interface,    ONLY: cosp_calipso_init
    USE mod_cosp_atlid_interface,      ONLY: cosp_atlid_init
    USE mod_cosp_grlidar532_interface, ONLY: cosp_grLidar532_init
    USE mod_cosp_parasol_interface,    ONLY: cosp_parasol_init
    USE mod_cosp_cloudsat_interface,   ONLY: cosp_radar_init => cosp_cloudsat_init
    implicit none

    ! Inputs
    logical, intent(in)    :: &
         do_cosp, do_isccp, do_misr, do_modis, do_cloudsat, do_calipso, do_grLidar532,      &
         do_atlid, do_parasol
    integer, intent(in)    ::  &
         mpirank,              & ! Current MPI rank 
         mpiroot,              & ! Master MPI
         imp_physics,          & ! Choice of microphysics scheme
         imp_physics_thompson, & ! Choice of Thompson
         imp_physics_gfdl,     & ! Choice of GFDL
         isccp_topht,          & !
         isccp_topht_dir,      & !
         cosp_nsubcol
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg              !
    integer, intent(out) :: &
         errflg              !

    ! Local
    logical :: exists
    integer :: ios
    integer, parameter :: spaceborne_radar  = 0
    integer, parameter :: groundbased_radar = 1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Has COSP been requested?
    if (.not. do_cosp) return

    ! ######################################################################################
    ! Initialize COSP
    ! ######################################################################################

    ! Initialize the distributional parameters for hydrometeors in radar simulator
    ! DJS: This needs to be revisited.
    if (imp_physics == imp_physics_thompson) then
       call hydro_class_init(.false., .true., csat_sd)
       csat_micro_scheme = 'MMF_v3.5_double_moment'
    endif
    if (imp_physics == imp_physics_gfdl) then
       call hydro_class_init(.true., .false., csat_sd)
       csat_micro_scheme = 'MMF_v3.5_single_moment'
    endif

    ! Initialize requested simulators
    ! DJS: In CAM cosp_init(), see src/cosp.F90, is called instead of the indivdual simulators.
    if (do_cloudsat .or. do_grLidar532) then
       call quickbeam_optics_init()
    endif
    if (do_isccp) then
       call cosp_isccp_init(isccp_topht, isccp_topht_dir)
       endif
    if (do_modis) then
       call cosp_modis_init()
    endif
    if (do_misr) then
       call cosp_misr_init()
    endif
    if (do_cloudsat) then
!       call cosp_radar_init(csat_freq, csat_k2, csat_gas_abs, csat_do_ray, R_UNDEF, N_HYDRO,      &
!            spaceborne_radar, csat_rcfg, csat_micro_scheme)
    endif
    if (do_calipso) then
       call cosp_calipso_init()
    endif
    if (do_grLidar532) then
       call cosp_grLidar532_init()
    endif
    if (do_atlid) then
       call cosp_atlid_init()
    endif
    if (do_parasol) then
       call cosp_parasol_init()
    endif

    !
    if (mpirank .eq. mpiroot) then
       if (do_cosp) then 
          print*,'COSP enabled:'
          print*,'  Number of COSP subcolumns                   = ', cosp_nsubcol
          print*,'  Enable Cloudsat RADAR simulator             = ', do_cloudsat
          print*,'  Enable Calipso LIDAR simulator              = ', do_calipso
          print*,'  Enable EarthCare LIDAR simulator            = ', do_atlid
          print*,'  Enable Ground-based (532nm) LIDAR simulator = ', do_grLidar532
          print*,'  Enable ISCCP simulator                      = ', do_isccp
          print*,'  Enable MISR simulator                       = ', do_misr
          print*,'  Enable MODIS simulator                      = ', do_modis
          print*,'  RADAR_SIM microphysics scheme               = ', trim(csat_micro_scheme)
       else
          print*, 'COSP not enabled'
       endif
    endif

  end subroutine cosp_simulator_init

  ! #########################################################################################
!! \section arg_table_cosp_simulator_run
!! \htmlinclude cosp_simulator_run.html
!!
  subroutine cosp_simulator_run(nCol, nLev, cosp_nlvgrid, cosp_nsubcol, tsfc, coszen, slmsk,&
       prsl, prsi, phil, phii, tgrs, qgrs, cldtau_lw, cldtau_sw, cld_frac, ccld_frac,       &
       top_at_1, con_g, iSFC, iTOA, do_cosp, do_isccp, do_misr, do_modis, do_cloudsat,      &
       do_calipso, do_grLidar532, do_atlid, do_parasol, overlap,                            &
       cosp_mp, cosp_mp_cam6, cosp_mp_ufs,                                                  &
       f1isccp_cosp, cldtot_isccp, meancldalb_isccp, meanptop_isccp, meantau_isccp,         &
       meantb_isccp, meantbclr_isccp, tau_isccp, cldptop_isccp, errmsg, errflg)
    use mod_cosp,  only: cosp_outputs, cosp_optical_inputs, cosp_column_inputs, cosp_simulator
    implicit none

    ! Inputs
    logical, intent(in) :: &
         do_cosp, do_isccp, do_misr, do_modis, do_cloudsat, do_calipso, do_grLidar532,      &
         do_atlid, do_parasol, &
         top_at_1              ! Vertical ordering flag
    integer, intent(in) :: &
         nCol,               & ! Number of horizontal grid points
         nLev,               & ! Number of vertical layers
         cosp_nlvgrid,       & !
         cosp_nsubcol,       & ! Number of COSP subcolumns
         overlap,            & ! Cloud overlap assumption
         iSFC,               & ! Vertical index for surface
         iTOA,               & ! Vertical index for TOA
         cosp_mp,            & ! Choice of subsampling and optics.
         cosp_mp_cam6,       & ! Choice of subsampling and optics CAM6 method.
         cosp_mp_ufs           ! Choice of subsampling and optics UFS method.
    real(kind_phys), intent(in) :: &
         con_g                 ! Physical constant: gravitational constant
    real(kind_phys), dimension(:), intent(in) :: & 
         tsfc,               & ! Surface skin temperature (K)
         coszen,             & ! Cosine of SZA
         slmsk                 ! Area type
    real(kind_phys), dimension(:,:), intent(in) :: & 
         prsl,               & ! Pressure at model-layer centers (Pa)
         tgrs,               & ! Temperature at model-layer centers (K)
         prsi,               & ! Pressure at model-interfaces (Pa)
         phii,               & ! Geopotential at model-interface (m2/s2)
         phil,               & ! Geopotential at model-layer centers
         cld_frac,           & ! Total cloud fraction
         ccld_frac,          & ! Convective cloud fraction
         cldtau_lw,          & ! In-cloud 10 micron optical depth
         cldtau_sw             ! In-cloud 0.67 micron optical depth
    real(kind_phys), dimension(:,:,:), intent(in) :: & 
         qgrs                  ! Tracer concentrations (kg/kg)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg
    real(kind_phys), dimension(:,:,:), intent(inout) :: &
         f1isccp_cosp        ! ISCCP CFAD
    real(kind_phys), dimension(:,:), intent(inout) :: &
         tau_isccp,        & ! ISCCP subcolumn optical-depth
         cldptop_isccp       ! ISCCP subcolumn cloud-top pressure
    real(kind_phys), dimension(:), intent(inout) :: &
         cldtot_isccp,     & ! ISCCP mean cloud-fraction
         meancldalb_isccp, & ! ISCCP mean cloud albedo
         meanptop_isccp,   & ! ISCCP mean cloud-top pressure
         meantau_isccp,    & ! ISCCP mean optical-depth
         meantb_isccp,     & ! ISCCP mean brightness temperature
         meantbclr_isccp     ! ISCCP mean brightness temperature (clear-sky)

    ! Local
    type(cosp_outputs)        :: cospOUT
    type(cosp_optical_inputs) :: cospIN
    type(cosp_column_inputs)  :: cospstateIN
    integer, dimension(nCol)  :: sunlit
    integer :: iCol, nerror, iErr, vs, iprs, itau, iSubCol
    character(len=256),dimension(100) :: cosp_status
    integer :: iSFCa = 1, iTOAa = 127
    logical :: top_at_1a = .false.

    if (.not. do_cosp) return

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Vertical stride direction
    if (top_at_1a)       vs = 1
    if (.not. top_at_1a) vs = -1
    
    ! Compute sunlit flag.
    sunlit(:) = 0
    do iCol = 1, nCol
       if (coszen(iCol) >= 0.0001) then
          sunlit(iCol) = 1
       endif
    enddo

    ! Type containing COSP outputs.
    call construct_cosp_outputs(do_isccp, do_modis, do_misr, do_cloudsat, do_calipso,    &
         do_grLidar532, do_atlid, do_parasol, nCol, cosp_nsubcol, nLev, cosp_nlvgrid,    &
         cospOUT)

    ! Host-model state for COSP (toa-2-sfc vertical ordering).
    call construct_cospstateIN(nCol, nLev, cospstateIN)
    cospstateIN%sunlit          = sunlit(:)
    cospstateIN%skt             = tsfc(:)
    cospstateIN%land            = slmsk(:)
    cospstateIN%at              = tgrs(:,iTOAa:iSFCa:vs)
    cospstateIN%pfull           = prsl(:,iTOAa:iSFCa:vs)
    cospstateIN%phalf           = prsi(:,iTOAa+1:iSFCa:vs)
    cospstateIN%qv              = qgrs(:,iTOAa:iSFCa:vs,1)
    cospstateIN%hgt_matrix      = phil/con_g
    cospstateIN%hgt_matrix_half = phii/con_g

    ! Derived (optical) inputs for COSP.
    call construct_cospIN(do_isccp, do_modis, do_misr, do_cloudsat, do_calipso,          &
         do_grLidar532, do_atlid, do_parasol, nCol, cosp_nsubcol, nLev, cospIN)

    !
    ! Call subsample_and_optics
    !
    ! Couple COSP to CAM6 microphysics.
    if (cosp_mp == cosp_mp_cam6) then
       call subsample_and_optics_CAM6(nCol, cosp_nsubcol, nLev, do_isccp, do_misr,       &
            do_modis, prsi(:,iSFCa), cld_frac, ccld_frac, overlap, cldtau_lw, cldtau_sw, &
            cospIN)
    endif
    ! Couple COSP to UFS microphysics.
    if (cosp_mp == cosp_mp_ufs) then
       call subsample_and_optics_UFS(cospIN)
    endif

    !
    ! Call COSP
    !
    cosp_status = cosp_simulator(cospIN, cospstateIN, cospOUT, start_idx=1,              &
         stop_idx=nCol, debug=.false.)

    ! Error checking
    nerror = 0
    do iErr = 1, ubound(cosp_status, 1)
       if (len_trim(cosp_status(iErr)) > 0) then
          errmsg = "cosp_simulator: ERROR: "//trim(cosp_status(iErr))
          nerror = nerror + 1
       end if
    end do
    if (nerror > 0) errflg = -1

    ! Set dark-scenes to fill value. Only done for passive simulators
    if (do_isccp) then
       ! 1D
       where(sunlit(1:nCol) .eq. 0)
          cospOUT%isccp_totalcldarea(1:nCol)  = R_UNDEF
          cospOUT%isccp_meanptop(1:nCol)      = R_UNDEF
          cospOUT%isccp_meantaucld(1:nCol)    = R_UNDEF
          cospOUT%isccp_meanalbedocld(1:nCol) = R_UNDEF
          cospOUT%isccp_meantb(1:nCol)        = R_UNDEF
          cospOUT%isccp_meantbclr(1:nCol)     = R_UNDEF
       end where
       ! 2D
       do iSubCol=1,cosp_nsubcol
          where (sunlit(1:nCol) .eq. 0)
             cospOUT%isccp_boxtau(1:nCol,iSubCol)  = R_UNDEF
             cospOUT%isccp_boxptop(1:nCol,iSubCol) = R_UNDEF
          end where
       enddo
       ! 3D
       do iprs=1,nprs_cosp
          do itau=1,ntau_cosp
             where(sunlit(1:nCol) .eq. 0)
                cospOUT%isccp_fq(1:nCol,iprs,itau) = R_UNDEF
             end where
          end do
       end do
    endif

    ! Copy COSP outputs to host interstitials
    f1isccp_cosp     = cospOUT%isccp_fq
    tau_isccp        = cospOUT%isccp_boxtau
    cldptop_isccp    = cospOUT%isccp_boxptop
    cldtot_isccp     = cospOUT%isccp_totalcldarea
    meanptop_isccp   = cospOUT%isccp_meanptop
    meantau_isccp    = cospOUT%isccp_meantaucld
    meancldalb_isccp = cospOUT%isccp_meanalbedocld
    meantb_isccp     = cospOUT%isccp_meantb
    meantbclr_isccp  = cospOUT%isccp_meantbclr

  end subroutine cosp_simulator_run
  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  ! ######################################################################################
  subroutine construct_cosp_outputs(do_isccp, do_modis, do_misr, do_cloudsat, do_calipso,&
       do_grLidar532, do_atlid, do_parasol, nCol, nSubCol, nLev, Nlvgrid, x)
    use mod_cosp,  only: cosp_outputs
    implicit none

    ! Inputs
    logical, intent(in) :: &
         do_isccp, do_modis, do_misr, do_cloudsat, do_calipso, do_grLidar532, do_atlid,  &
         do_parasol
    integer, intent(in) :: &
         nCol, nSubCol, nLev, Nlvgrid
    
    ! Outputs
    type(cosp_outputs),intent(out) :: &
         x           ! COSP output structure  
  
     ! ISCCP simulator outputs
    if (do_isccp) then
       allocate(x%isccp_boxtau(nCol,nSubCol)) 
       allocate(x%isccp_boxptop(nCol,nSubCol))
       allocate(x%isccp_fq(nCol,numISCCPTauBins,numISCCPPresBins))
       allocate(x%isccp_totalcldarea(nCol))
       allocate(x%isccp_meanptop(nCol))
       allocate(x%isccp_meantaucld(nCol))
       allocate(x%isccp_meantb(nCol))
       allocate(x%isccp_meantbclr(nCol))
       allocate(x%isccp_meanalbedocld(nCol))
    endif

    ! MISR simulator
    if (do_misr) then 
       allocate(x%misr_fq(nCol,numMISRTauBins,numMISRHgtBins))
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.
       allocate(x%misr_dist_model_layertops(nCol,numMISRHgtBins))
       allocate(x%misr_meanztop(nCol))
       allocate(x%misr_cldarea(nCol))    
    endif
    
    ! MODIS simulator
    if (do_modis) then
       allocate(x%modis_Cloud_Fraction_Total_Mean(nCol))
       allocate(x%modis_Cloud_Fraction_Water_Mean(nCol))
       allocate(x%modis_Cloud_Fraction_Ice_Mean(nCol))
       allocate(x%modis_Cloud_Fraction_High_Mean(nCol))
       allocate(x%modis_Cloud_Fraction_Mid_Mean(nCol))
       allocate(x%modis_Cloud_Fraction_Low_Mean(nCol))
       allocate(x%modis_Optical_Thickness_Total_Mean(nCol))
       allocate(x%modis_Optical_Thickness_Water_Mean(nCol))
       allocate(x%modis_Optical_Thickness_Ice_Mean(nCol))
       allocate(x%modis_Optical_Thickness_Total_LogMean(nCol))
       allocate(x%modis_Optical_Thickness_Water_LogMean(nCol))
       allocate(x%modis_Optical_Thickness_Ice_LogMean(nCol))
       allocate(x%modis_Cloud_Particle_Size_Water_Mean(nCol))
       allocate(x%modis_Cloud_Particle_Size_Ice_Mean(nCol))
       allocate(x%modis_Cloud_Top_Pressure_Total_Mean(nCol))
       allocate(x%modis_Liquid_Water_Path_Mean(nCol))
       allocate(x%modis_Ice_Water_Path_Mean(nCol))
       allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nCol,numModisTauBins,numMODISPresBins))
       allocate(x%modis_Optical_thickness_vs_ReffLIQ(nCol,numMODISTauBins,numMODISReffLiqBins))   
       allocate(x%modis_Optical_Thickness_vs_ReffICE(nCol,numMODISTauBins,numMODISReffIceBins))
    endif
    
    ! CALIPSO simulator
    if (do_calipso) then
       allocate(x%calipso_beta_mol(nCol,nLev))
       allocate(x%calipso_beta_tot(nCol,nSubCol,nLev))
       allocate(x%calipso_srbval(SR_BINS+1))
       allocate(x%calipso_cfad_sr(nCol,SR_BINS,Nlvgrid))
       allocate(x%calipso_betaperp_tot(nCol,nSubCol,nLev))  
       allocate(x%calipso_lidarcld(nCol,Nlvgrid))
       allocate(x%calipso_cldlayer(nCol,LIDAR_NCAT))        
       allocate(x%calipso_lidarcldphase(nCol,Nlvgrid,6))
       allocate(x%calipso_lidarcldtmp(nCol,LIDAR_NTEMP,5))
       allocate(x%calipso_cldlayerphase(nCol,LIDAR_NCAT,6))     
       ! These 2 outputs are part of the calipso output type, but are not controlled by an 
       ! logical switch in the output namelist, so if all other fields are on, then allocate
       allocate(x%calipso_tau_tot(nCol,nSubCol,nLev))       
       allocate(x%calipso_temp_tot(nCol,nLev))               
       ! Calipso opaque cloud diagnostics
       allocate(x%calipso_cldtype(nCol,LIDAR_NTYPE))
       allocate(x%calipso_cldtypetemp(nCol,LIDAR_NTYPE))  
       allocate(x%calipso_cldtypemeanz(nCol,2)) 
       allocate(x%calipso_cldtypemeanzse(nCol,3)) 
       allocate(x%calipso_cldthinemis(nCol))
       allocate(x%calipso_lidarcldtype(nCol,Nlvgrid,LIDAR_NTYPE+1))
    endif 
      
    ! PARASOL
    if (do_parasol) then
       allocate(x%parasolPix_refl(nCol,nSubCol,PARASOL_NREFL))
       allocate(x%parasolGrid_refl(nCol,PARASOL_NREFL))
    endif

    ! Cloudsat simulator
    if (do_cloudsat) then
       allocate(x%cloudsat_Ze_tot(nCol,nSubCol,nLev))
       allocate(x%cloudsat_cfad_ze(nCol,CLOUDSAT_DBZE_BINS,Nlvgrid))
       allocate(x%lidar_only_freq_cloud(nCol,Nlvgrid))
       allocate(x%radar_lidar_tcc(nCol))
       allocate(x%cloudsat_precip_cover(nCol,nCloudsatPrecipClass))
       allocate(x%cloudsat_pia(nCol))
    endif

  end subroutine construct_cosp_outputs

  ! ######################################################################################
  ! SUBROUTINE construct_cospstate
  ! ######################################################################################
  subroutine construct_cospstateIN(nCol, nLev, y)
    use mod_cosp,  only: cosp_column_inputs
    implicit none

    ! Inputs
    integer,intent(in) :: &
         nCol, nLev
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y
    
    allocate(y%sunlit(nCol),y%skt(nCol),y%land(nCol),y%at(nCol,nLev), y%pfull(nCol,nLev),&
         y%phalf(nCol,nLev+1),y%qv(nCol,nLev), y%hgt_matrix(nCol,nLev),                  &
         y%hgt_matrix_half(nCol,nLev+1))

  end subroutine construct_cospstateIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(do_isccp, do_modis, do_misr, do_cloudsat, do_calipso,      &
       do_grLidar532, do_atlid, do_parasol, nCol, nSubCol, nLev, y)
    use mod_cosp,  only: cosp_optical_inputs
    implicit none

    ! Inputs
    integer,intent(in) :: &
         nCol, nSubCol, nLev
    logical, intent(in) :: &
         do_isccp, do_modis, do_misr, do_cloudsat, do_calipso, do_grLidar532, do_atlid,  &
         do_parasol
    ! Outputs 
    type(cosp_optical_inputs),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = nCol
    y%Ncolumns = nSubCol
    y%Nlevels  = nLev

    if (do_isccp) then
       if (.not. allocated(y%frac_out)) allocate(y%frac_out(nCol, nSubCol, nLev))
       if (.not. allocated(y%emiss_11)) allocate(y%emiss_11(nCol, nSubCol, nLev))
       if (.not. allocated(y%tau_067))  allocate(y%tau_067( nCol, nSubCol, nLev))
    endif
    if (do_misr) then
       if (.not. allocated(y%tau_067))  allocate(y%tau_067( nCol, nSubCol, nLev))
    endif
    if (do_modis) then
       if (.not. allocated(y%fracLiq))  allocate(y%fracLiq(nCol, nSubCol, nLev))
       if (.not. allocated(y%tau_067))  allocate(y%tau_067(nCol, nSubCol, nLev))
       if (.not. allocated(y%asym))     allocate(y%asym(   nCol, nSubCol, nLev))
       if (.not. allocated(y%ss_alb))   allocate(y%ss_alb( nCol, nSubCol, nLev))
    endif
    if (do_cloudsat) then
       if (.not. allocated(y%z_vol_cloudsat))  allocate(y%z_vol_cloudsat(  nCol, nSubCol, nLev))
       if (.not. allocated(y%kr_vol_cloudsat)) allocate(y%kr_vol_cloudsat( nCol, nSubCol, nLev))
       if (.not. allocated(y%g_vol_cloudsat))  allocate(y%g_vol_cloudsat(  nCol, nSubCol, nLev))
       if (.not. allocated(y%fracPrecipIce))   allocate(y%fracPrecipIce(   nCol, nSubCol))
    endif
    if (do_calipso) then
       if (.not. allocated(y%betatot_calipso))     allocate(y%betatot_calipso(    nCol, nSubCol, nLev))
       if (.not. allocated(y%betatot_ice_calipso)) allocate(y%betatot_ice_calipso(nCol, nSubCol, nLev))
       if (.not. allocated(y%betatot_liq_calipso)) allocate(y%betatot_liq_calipso(nCol, nSubCol, nLev))
       if (.not. allocated(y%tautot_calipso))      allocate(y%tautot_calipso(     nCol, nSubCol, nLev))
       if (.not. allocated(y%tautot_ice_calipso))  allocate(y%tautot_ice_calipso( nCol, nSubCol, nLev))
       if (.not. allocated(y%tautot_liq_calipso))  allocate(y%tautot_liq_calipso( nCol, nSubCol, nLev))
    endif
    if (do_parasol) then
       y%Nrefl = PARASOL_NREFL
       if (.not. allocated(y%tautot_S_ice)) allocate(y%tautot_S_ice(nCol, nSubCol))
       if (.not. allocated(y%tautot_S_liq)) allocate(y%tautot_S_liq(nCol, nSubCol))
    endif

  end subroutine construct_cospIN

end module cosp_simulator
