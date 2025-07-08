! ###########################################################################################
!> \file GFS_cosp.F90
!!
!> \defgroup GFS_cosp GFS_cosp.F90
!!
!! \brief This module contains the CCPP interface to the Cloud-Feedback Model Intercomparison
!! Project (CFMIP) Observational Simulator Package Version 2.0 (COSP)
!!
! ###########################################################################################
module GFS_cosp
  use machine,  only: kind_phys
  use mod_cosp, only: cosp_outputs, cosp_optical_inputs, cosp_column_inputs

  implicit none
  real(kind_phys), parameter :: emsfc_lw = 0.99_kind_phys ! longwave emissivity of surface at 10.5 microns

contains

! ###########################################################################################
!! \section arg_table_GFS_cosp_init
!! \htmlinclude GFS_cosp_init.html
!!
!!
!> \ingroup GFS_cosp
!!
!! \brief
!!
!! \section GFS_cosp_init
!> @{
! ###########################################################################################
  subroutine GFS_cosp_init(mpirank, mpiroot, do_cosp, do_isccp, do_misr, do_modis,          &
       cosp_nsubcol, imp_physics, imp_physics_thompson, imp_physics_gfdl, isccp_topht,      &
       isccp_topht_dir, errmsg, errflg)
    use mod_cosp_config,          only: modis_histTau, modis_histTauEdges, modis_histTauCenters
    use mod_cosp_config,          only: numMODISTauBins, ntau
    use mod_cosp_config,          only: tau_binBounds, tau_binEdges, tau_binCenters
    use mod_cosp_modis_interface, only: cosp_modis_init
    use mod_cosp_misr_interface,  only: cosp_misr_init
    use mod_cosp_isccp_interface, only: cosp_isccp_init

    ! Inputs
    logical, intent(in)    :: &
         do_cosp,              & ! Flag for COSP diagnostics
	 do_isccp,             & ! Flag for COSP ISCCP diagnostics
	 do_misr,              & ! Flag for COSP MISR diagnostics
	 do_modis                ! Flag for COSP MODIS diagnostics
    integer, intent(in)    ::  &
         mpirank,              & ! Current MPI rank 
         mpiroot,              & ! Master MPI
         imp_physics,          & ! Choice of microphysics scheme
         imp_physics_thompson, & ! Choice of Thompson
         imp_physics_gfdl,     & ! Choice of GFDL
         isccp_topht,          & ! Cloud top height adjustment in cosp isccp simulator
         isccp_topht_dir,      & ! Cloud top height direction in cosp isccp simulator
         cosp_nsubcol            ! Number of COSP subcolumns.

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                ! CCPP error message
    integer, intent(out) :: &
         errflg                ! CCPP error flag

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Has COSP been requested?
    if (.not. do_cosp) return

    ! Initialize requested simulators
    if (do_isccp) then
       call cosp_isccp_init(isccp_topht, isccp_topht_dir)
    endif
    if (do_modis) then
       ! Initialize MODIS optical-depth bin boundaries for joint-histogram. (defined in cosp_config.F90)
       if (.not. allocated(modis_histTau)) then
          allocate(modis_histTau(ntau+1),modis_histTauEdges(2,ntau),modis_histTauCenters(ntau))
          numMODISTauBins      = ntau
          modis_histTau        = tau_binBounds
          modis_histTauEdges   = tau_binEdges
          modis_histTauCenters = tau_binCenters
    endif
       call cosp_modis_init()
    endif
    if (do_misr) then
       call cosp_misr_init()
    endif

    !
    if (mpirank .eq. mpiroot) then
       print*,'CFMIP Observational Simulator Package (COSP) enabled:'
       print*,'  Number of COSP subcolumns = ', cosp_nsubcol
       print*,'  Enable ISCCP simulator    = ', do_isccp
       if (do_isccp) then
          print*,'     ISCCP top height           = ',isccp_topht
	  print*,'     ISCCP top height direction = ',isccp_topht_dir
       endif
       print*,'  Enable MISR simulator     = ', do_misr
       print*,'  Enable MODIS simulator    = ', do_modis
    endif

  end subroutine GFS_cosp_init
!> @}

! ###########################################################################################
!! \section arg_table_GFS_cosp_run
!! \htmlinclude GFS_cosp_run.html
!!
!> \ingroup GFS_cosp
!!
!! \brief
!!
!! \section GFS_cosp_run
!> @{
! ###########################################################################################
  subroutine GFS_cosp_run(nCol, nLay, cosp_nlvgrid, cosp_nsubcol, tsfc, coszen, slmsk,      &
       prsl, prsi, phil, phii, tgrs, qgrs, cldtau_lw, cldtau_sw, cld_frac, ccld_frac,       &
       top_at_1, con_g, cld_liq, cld_ice, cld_rain, cld_snow, cld_graupel, ccld_liq,        &
       cld_reliq, cld_reice, cld_rerain, cld_resnow, iSFC, iTOA,                            &
       n_isccp_pres_bins,  isccp_pres_bins,  n_isccp_tau_bins,   isccp_tau_bins,            &
       n_modis_pres_bins,  modis_pres_bins,  n_modis_tau_bins,   modis_tau_bins,            &
       n_modis_reffi_bins, modis_reffi_bins, n_modis_reffl_bins, modis_reffl_bins,          &
       n_misr_hgt_bins,    misr_hgt_bins,    n_misr_tau_bins,    misr_tau_bins,             &
       doSWrad, doLWrad, do_cosp, do_isccp, do_misr, do_modis, overlap,                     &
       f1isccp_cosp, cldtot_isccp, meancldalb_isccp, meanptop_isccp, meantau_isccp,         &
       meantb_isccp, meantbclr_isccp, tau_isccp, cldptop_isccp,                             &
       clt_modis, clw_modis, cli_modis, clh_modis, clm_modis, cll_modis, taut_modis,        &
       tauw_modis, taui_modis, tautlog_modis, tauwlog_modis, tauilog_modis, reffclw_modis,  &
       reffcli_modis, pct_modis, lwp_modis, iwp_modis, cl_modis, clri_modis, clrl_modis,    &
       errmsg, errflg)
    use mod_cosp,        only: cosp_simulator
    use mod_cosp_config, only: R_UNDEF
    real(kind_phys), parameter :: missing_value = 0._kind_phys!9.99e20_kind_phys

    ! Inputs
    logical, intent(in) :: &
         doSWrad,            & ! Logical flags for sw radiation calls
         doLWrad,            & ! Logical flags for lw radiation calls
         do_cosp,            & ! Flag for COSP diagnostics
	 do_isccp,           & ! Flag for COSP ISCCP diagnostics
	 do_misr,            & ! Flag for COSP MISR diagnostics
	 do_modis,           & ! Flag for COSP MODIS diagnostics
         top_at_1              ! Vertical ordering flag
    integer, intent(in) :: &
         nCol,               & ! Number of horizontal grid points
         nLay,               & ! Number of vertical layers
         cosp_nlvgrid,       & ! Number of vertical layers in COSP statistical grid.
         cosp_nsubcol,       & ! Number of COSP subcolumns
         n_isccp_pres_bins,  & ! Number of pressure      bins in ISCCP CFAD.
	 n_isccp_tau_bins,   & ! Number of optical-depth bins in ISCCP CFAD.
         n_modis_pres_bins,  & ! Number of pressure      bins in MODIS CFAD.
         n_modis_tau_bins,   & ! Number of optical-depth bins in MODIS CFAD.
         n_modis_reffi_bins, & ! Number of ice-radii     bins in MODIS CFAD.
         n_modis_reffl_bins, & ! Number of liquid-radii  bins in MODIS CFAD.
         n_misr_hgt_bins,    & ! Number of height        bins in MISR CFAD.
         n_misr_tau_bins,    & ! Number of optical-depth bins in MISR CFAD.
         overlap,            & ! Cloud overlap assumption
         iSFC,               & ! Vertical index for surface
         iTOA                  ! Vertical index for TOA
    real(kind_phys), intent(in) :: &
         con_g                 ! Physical constant: gravitational constant
    real(kind_phys), dimension(:), intent(in) :: & 
         tsfc,               & ! Surface skin temperature (K)
         coszen,             & ! Cosine of SZA
         slmsk,              & ! Area type
	 isccp_pres_bins,    & ! Pressure bin boundaries for ISCCP CFAD.
	 isccp_tau_bins,     & ! Optical-depth bin boundaries for ISCCP CFAD.
         modis_pres_bins,    & ! Pressure bin boundaries for MODIS CFAD.
         modis_tau_bins,     & ! Optical-depth bin boundaries for  MODIS CFAD.
         modis_reffi_bins,   & ! Ice-radii bin boundaries for  MODIS CFAD.
         modis_reffl_bins,   & ! Liquid-radii bin boundaries for  MODIS CFAD.
         misr_hgt_bins,      & ! Pressure bin boundaries for MISR CFAD.
         misr_tau_bins         ! Optical-depth bin boundaries for MISR CFAD.
    real(kind_phys), dimension(:,:), intent(in) :: & 
         prsl,               & ! Pressure at model-layer centers (Pa)
         tgrs,               & ! Temperature at model-layer centers (K)
         prsi,               & ! Pressure at model-interfaces (Pa)
         phii,               & ! Geopotential at model-interface (m2/s2)
         phil,               & ! Geopotential at model-layer centers
         cld_frac,           & ! Total cloud fraction
         cld_liq,            & ! Liquid cloud water mixing ratio (kg/kg)
         cld_ice,            & ! Ice cloud water mixing ratio (kg/kg)
         cld_rain,           & ! Rain cloud water mixing ratio (kg/kg)
         cld_snow,           & ! Snow cloud water mixing ratio (kg/kg)
         cld_graupel,        & ! Graupel cloud water mixing ratio (kg/kg)
         cldtau_lw,          & ! In-cloud 10 micron optical depth
         cldtau_sw             ! In-cloud 0.67 micron optical depth
    real(kind_phys), dimension(:,:), intent(in), optional :: &
         ccld_liq,           & ! Convective cloud water mixing ratio (kg/kg)
         ccld_frac,          & ! Convective cloud fraction
         cld_reliq,          &
         cld_reice,          &
         cld_rerain,         &
         cld_resnow
    real(kind_phys), dimension(:,:,:), intent(in) :: & 
         qgrs                  ! Tracer concentrations (kg/kg)

    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                ! CCPP error message
    integer, intent(out) :: &
         errflg                ! CCPP error flag
    real(kind_phys), dimension(:,:,:), intent(out) :: &
         f1isccp_cosp,       & ! ISCCP CFAD
         cl_modis,           & ! MODIS CFAD
         clri_modis,         & ! MODIS CFAD
         clrl_modis            ! MODIS CFAD
    real(kind_phys), dimension(:,:), intent(out) :: &
         tau_isccp,          & ! ISCCP subcolumn optical-depth
         cldptop_isccp         ! ISCCP subcolumn cloud-top pressure
    real(kind_phys), dimension(:), intent(out) :: &
         cldtot_isccp,       & ! ISCCP mean cloud-fraction
         meancldalb_isccp,   & ! ISCCP mean cloud albedo
         meanptop_isccp,     & ! ISCCP mean cloud-top pressure
         meantau_isccp,      & ! ISCCP mean optical-depth
         meantb_isccp,       & ! ISCCP mean brightness temperature
         meantbclr_isccp,    & ! ISCCP mean brightness temperature (clear-sky)
         clt_modis,          & ! MODIS
         clw_modis,          & ! MODIS
         cli_modis,          & ! MODIS
         clh_modis,          & ! MODIS
         clm_modis,          & ! MODIS
         cll_modis,          & ! MODIS
         taut_modis,  	     & ! MODIS
         tauw_modis,         & ! MODIS
         taui_modis,         & ! MODIS
         tautlog_modis,      & ! MODIS
         tauwlog_modis,      & ! MODIS
         tauilog_modis,      & ! MODIS
         reffclw_modis,      & ! MODIS
         reffcli_modis,      & ! MODIS
         pct_modis,          & ! MODIS
         lwp_modis,          & ! MODIS
         iwp_modis             ! MODIS

    ! Local
    type(cosp_outputs)        :: cospOUT
    type(cosp_optical_inputs) :: cospIN
    type(cosp_column_inputs)  :: cospstateIN
    integer, dimension(nCol)  :: sunlit
    integer :: iCol, nerror, iErr, vs, iprs, itau, iSubCol
    character(len=256),dimension(100) :: cosp_status

    if (.not. do_cosp) return

    ! Only call COSP on radiation time-step.
    if (.not. (doLWrad .or. doSWrad)) return
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Initialize.
!    f1isccp_cosp     = 0._kind_phys
!    tau_isccp        = 0._kind_phys
!    cldptop_isccp    = 0._kind_phys
!    cldtot_isccp     = 0._kind_phys
!    meanptop_isccp   = 0._kind_phys
!    meantau_isccp    = 0._kind_phys
!    meancldalb_isccp = 0._kind_phys
!    meantb_isccp     = 0._kind_phys
!    meantbclr_isccp  = 0._kind_phys

    ! Vertical stride direction
    if (top_at_1)       vs = 1
    if (.not. top_at_1) vs = -1
    
    ! Compute sunlit flag.
    sunlit(:) = 0
    do iCol = 1, nCol
       if (coszen(iCol) > 0._kind_phys) then
          sunlit(iCol) = 1
       endif
    enddo

    ! Type containing COSP outputs.
    call construct_cosp_outputs(do_isccp, do_modis, do_misr, nCol, cosp_nsubcol, nLay,      &
    	 cosp_nlvgrid, n_isccp_pres_bins, n_isccp_tau_bins, n_modis_pres_bins,              &
	 n_modis_tau_bins, n_modis_reffl_bins, n_modis_reffi_bins, n_misr_hgt_bins,         &
         n_misr_tau_bins, cospOUT)

    ! Host-model state for COSP (toa-2-sfc vertical ordering).
    call construct_cospstateIN(nCol, nLay, cospstateIN)
    cospstateIN%sunlit          = sunlit(:)
    cospstateIN%skt             = tsfc(:)
    cospstateIN%land            = slmsk(:)
    cospstateIN%at              = tgrs(:,iTOA:iSFC:vs)
    cospstateIN%pfull           = prsl(:,iTOA:iSFC:vs)
    cospstateIN%phalf           = prsi(:,iTOA-vs:iSFC:vs)
    cospstateIN%qv              = qgrs(:,iTOA:iSFC:vs,1)
    cospstateIN%hgt_matrix      = phil(:,iTOA:iSFC:vs)/con_g
    cospstateIN%hgt_matrix_half = phii(:,iTOA-vs:iSFC:vs)/con_g
    where(cospstateIN%qv .lt. 1.e-6) cospstateIN%qv = 1.e-6

    ! Derived (optical) inputs for COSP.
    call construct_cospIN(do_isccp, do_modis, do_misr, nCol, cosp_nsubcol, nLay, cospIN)
    cospIN%emsfc_lw = emsfc_lw

    !
    ! Call subsample_and_optics
    !
    call subsample_and_optics(nCol, cosp_nsubcol, nLay, do_isccp, do_misr, do_modis,     &
    	 prsi(:,iSFC), cld_frac, ccld_frac, overlap, cldtau_lw, cldtau_sw, cld_liq,&
         ccld_liq, cld_ice, cld_rain, cld_snow, cld_graupel, cld_reliq, cld_reice, cld_rerain, cld_resnow, cospIN)

    !
    ! Call COSP
    !
    cosp_status = cosp_simulator(cospIN, cospstateIN, cospOUT, start_idx=1, stop_idx=nCol,debug=.false.)

    ! Error checking
    nerror = 0
    do iErr = 1, ubound(cosp_status, 1)
       if (len_trim(cosp_status(iErr)) > 0) then
          errmsg = "cosp: ERROR: "//trim(cosp_status(iErr))
          nerror = nerror + 1
       end if
    end do
    if (nerror > 0) errflg = -1

    !
    ! Replace COSP masking.
    !
    if (do_isccp) then
       where(cospOUT%isccp_totalcldarea  .eq. R_UNDEF) cospOUT%isccp_totalcldarea  = missing_value
       where(cospOUT%isccp_meanptop      .eq. R_UNDEF) cospOUT%isccp_meanptop      = missing_value
       where(cospOUT%isccp_meantaucld    .eq. R_UNDEF) cospOUT%isccp_meantaucld    = missing_value
       where(cospOUT%isccp_meanalbedocld .eq. R_UNDEF) cospOUT%isccp_meanalbedocld = missing_value
       where(cospOUT%isccp_meantb        .eq. R_UNDEF) cospOUT%isccp_meantb        = missing_value
       where(cospOUT%isccp_meantbclr     .eq. R_UNDEF) cospOUT%isccp_meantbclr     = missing_value
    endif
    if (do_modis) then
       where(cospOUT%modis_Cloud_Fraction_Total_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Fraction_Total_Mean       = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Fraction_Water_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Fraction_Water_Mean       = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Fraction_Ice_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Fraction_Ice_Mean         = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Fraction_High_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Fraction_High_Mean        = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Fraction_Mid_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Fraction_Mid_Mean         = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Fraction_Low_Mean  .eq. R_UNDEF)
          cospOUT%modis_Cloud_Fraction_Low_Mean         = missing_value
       endwhere
       where(cospOUT%modis_Optical_Thickness_Total_Mean .eq. R_UNDEF)
          cospOUT%modis_Optical_Thickness_Total_Mean    = missing_value
       endwhere
       where(cospOUT%modis_Optical_Thickness_Water_Mean .eq. R_UNDEF)
          cospOUT%modis_Optical_Thickness_Water_Mean    = missing_value
       endwhere
       where(cospOUT%modis_Optical_Thickness_Ice_Mean .eq. R_UNDEF)
          cospOUT%modis_Optical_Thickness_Ice_Mean      = missing_value
       endwhere
       where(cospOUT%modis_Optical_Thickness_Total_LogMean .eq. R_UNDEF)
          cospOUT%modis_Optical_Thickness_Total_LogMean = missing_value
       endwhere
       where(cospOUT%modis_Optical_Thickness_Water_LogMean .eq. R_UNDEF)
          cospOUT%modis_Optical_Thickness_Water_LogMean = missing_value
       endwhere
       where(cospOUT%modis_Optical_Thickness_Ice_LogMean .eq. R_UNDEF)
          cospOUT%modis_Optical_Thickness_Ice_LogMean   = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Particle_Size_Water_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Particle_Size_Water_Mean  = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Particle_Size_Ice_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Particle_Size_Ice_Mean    = missing_value
       endwhere
       where(cospOUT%modis_Cloud_Top_Pressure_Total_Mean .eq. R_UNDEF)
          cospOUT%modis_Cloud_Top_Pressure_Total_Mean   = missing_value
       endwhere
       where(cospOUT%modis_Liquid_Water_Path_Mean .eq. R_UNDEF)
          cospOUT%modis_Liquid_Water_Path_Mean          = missing_value
       endwhere
       where(cospOUT%modis_Ice_Water_Path_Mean .eq. R_UNDEF)
          cospOUT%modis_Ice_Water_Path_Mean             = missing_value
       endwhere
    endif
    
    ! Set dark-scenes to fill value. Only done for passive simulators
    if (do_isccp) then
       ! 1D
       where(sunlit(1:nCol) .eq. 0)
          cospOUT%isccp_totalcldarea(1:nCol)  = missing_value
          cospOUT%isccp_meanptop(1:nCol)      = missing_value
          cospOUT%isccp_meantaucld(1:nCol)    = missing_value
          cospOUT%isccp_meanalbedocld(1:nCol) = missing_value
          cospOUT%isccp_meantb(1:nCol)        = missing_value
          cospOUT%isccp_meantbclr(1:nCol)     = missing_value
       end where
       ! 2D
       do iSubCol=1,cosp_nsubcol
          where (sunlit(1:nCol) .eq. 0)
             cospOUT%isccp_boxtau(1:nCol,iSubCol)  = missing_value
             cospOUT%isccp_boxptop(1:nCol,iSubCol) = missing_value
          end where
       enddo
       ! 3D
       do iprs=1,n_isccp_pres_bins
          do itau=1,n_isccp_tau_bins
             where(sunlit(1:nCol) .eq. 0)
                cospOUT%isccp_fq(1:nCol,iprs,itau) = missing_value
             end where
          end do
       end do
    endif
    if (do_misr) then
       do iprs=1,n_misr_hgt_bins
          do itau=1,n_misr_tau_bins
             where(sunlit(1:ncol) .eq. 0)
                cospOUT%misr_fq(1:ncol,itau,iprs) = missing_value
             end where
          end do
       end do
    end if
    if (do_modis) then
       ! 1D
       where(sunlit(1:nCol) .eq. 0)
          cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol)       = missing_value
          cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol)       = missing_value
          cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol)         = missing_value
          cospOUT%modis_Cloud_Fraction_High_Mean(1:ncol)        = missing_value
          cospOUT%modis_Cloud_Fraction_Mid_Mean(1:ncol)         = missing_value
          cospOUT%modis_Cloud_Fraction_Low_Mean(1:ncol)         = missing_value
          cospOUT%modis_Optical_Thickness_Total_Mean(1:ncol)    = missing_value
          cospOUT%modis_Optical_Thickness_Water_Mean(1:ncol)    = missing_value
          cospOUT%modis_Optical_Thickness_Ice_Mean(1:ncol)      = missing_value
          cospOUT%modis_Optical_Thickness_Total_LogMean(1:ncol) = missing_value
          cospOUT%modis_Optical_Thickness_Water_LogMean(1:ncol) = missing_value
          cospOUT%modis_Optical_Thickness_Ice_LogMean(1:ncol)   = missing_value
          cospOUT%modis_Cloud_Particle_Size_Water_Mean(1:ncol)  = missing_value
          cospOUT%modis_Cloud_Particle_Size_Ice_Mean(1:ncol)    = missing_value
          cospOUT%modis_Cloud_Top_Pressure_Total_Mean(1:ncol)   = missing_value
          cospOUT%modis_Liquid_Water_Path_Mean(1:ncol)          = missing_value
          cospOUT%modis_Ice_Water_Path_Mean(1:ncol)             = missing_value
       end where
       ! 3D
       do iprs=1,n_modis_pres_bins
          do itau=1,n_modis_tau_bins
             where(sunlit(1:ncol) .eq. 0)
                cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(1:ncol,itau,iprs) = missing_value
             end where
          enddo
       enddo
       do iprs=1,n_modis_reffi_bins
          do itau=1,n_modis_tau_bins
             where(sunlit(1:ncol) .eq. 0)
                cospOUT%modis_Optical_Thickness_vs_ReffICE(1:ncol,itau,iprs) = missing_value
             end where
          end do
       enddo
       do iprs=1,n_modis_reffl_bins
          do itau=1,n_modis_tau_bins
             where(sunlit(1:ncol) .eq. 0)
                cospOUT%modis_Optical_Thickness_vs_ReffLIQ(1:ncol,itau,iprs) = missing_value
             end where
          enddo
       enddo
    end if
    
    ! Copy COSP outputs to host interstitials.
    if (do_isccp) then
       f1isccp_cosp     = cospOUT%isccp_fq
       tau_isccp        = cospOUT%isccp_boxtau
       cldptop_isccp    = cospOUT%isccp_boxptop
       cldtot_isccp     = cospOUT%isccp_totalcldarea
       meanptop_isccp   = cospOUT%isccp_meanptop
       meantau_isccp    = cospOUT%isccp_meantaucld
       meancldalb_isccp = cospOUT%isccp_meanalbedocld
       meantb_isccp     = cospOUT%isccp_meantb
       meantbclr_isccp  = cospOUT%isccp_meantbclr
    endif
    if (do_misr) then

    endif
    if (do_modis) then
       clt_modis        = cospOUT%modis_Cloud_Fraction_Total_Mean
       clw_modis        = cospOUT%modis_Cloud_Fraction_Water_Mean
       cli_modis        = cospOUT%modis_Cloud_Fraction_Ice_Mean
       clh_modis        = cospOUT%modis_Cloud_Fraction_High_Mean
       clm_modis        = cospOUT%modis_Cloud_Fraction_Mid_Mean
       cll_modis        = cospOUT%modis_Cloud_Fraction_Low_Mean
       taut_modis       = cospOUT%modis_Optical_Thickness_Total_Mean
       tauw_modis       = cospOUT%modis_Optical_Thickness_Water_Mean
       taui_modis       = cospOUT%modis_Optical_Thickness_Ice_Mean
       tautlog_modis    = cospOUT%modis_Optical_Thickness_Total_LogMean
       tauwlog_modis    = cospOUT%modis_Optical_Thickness_Water_LogMean
       tauilog_modis    = cospOUT%modis_Optical_Thickness_Ice_LogMean
       reffclw_modis    = cospOUT%modis_Cloud_Particle_Size_Water_Mean
       reffcli_modis    = cospOUT%modis_Cloud_Particle_Size_Ice_Mean
       pct_modis        = cospOUT%modis_Cloud_Top_Pressure_Total_Mean
       lwp_modis        = cospOUT%modis_Liquid_Water_Path_Mean
       iwp_modis        = cospOUT%modis_Ice_Water_Path_Mean
       cl_modis         = cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure
       clri_modis       = cospOUT%modis_Optical_Thickness_vs_ReffICE
       clrl_modis       = cospOUT%modis_Optical_Thickness_vs_ReffLIQ
    endif

    ! Clean up
    call destroy_cospIN(cospIN)
    call destroy_cospstateIN(cospstateIN)
    call destroy_cosp_outputs(cospOUT) 

  end subroutine GFS_cosp_run
!> @}

  ! #########################################################################################
  ! SUBROUTINE subsample_and_optics
  !
  ! This routine contains the radiaiton to cloud coupling needed by COSP. Changes to host
  ! cloud and radiative configurations need to also occur here.
  !
  ! #########################################################################################
  subroutine subsample_and_optics(nCol, nSubCol, nLay, do_isccp, do_misr, do_modis,         &
       sfcP, cld_frac, ccld_frac, overlap, cldtau_lw, cldtau_sw, cld_liq, ccld_liq, cld_ice,&
       cld_rain, cld_snow, cld_graupel, cld_reliq, cld_reice, cld_rerain, cld_resnow, cospIN)
    use cosp_optics,     only: modis_optics, modis_optics_partition
    use mod_scops,       only: scops
    use mod_prec_scops,  only: prec_scops
    use mod_rng,         only: rng_state, init_rng
    use mod_cosp_config, only: nHydro => N_HYDRO
    ! Inputs
    logical, intent(in) :: &
    	 do_isccp,  & ! Flag for COSP ISCCP diagnostics
	 do_misr,   & ! Flag for COSP MISR diagnostics
	 do_modis     ! Flag for COSP MODIS diagnostics
    integer, intent(in) :: &
    	 nCol,      & ! Number of horizontal gridpoints
	 nSubCol,   & ! Number of COSP subcolumns
	 nLay,      & ! Number of vertical layers
	 overlap      ! Cloud overlap assumption
    real(kind_phys), dimension(nCol), intent(in) :: &
         sfcP         ! Pressure @ surface (Pa)
    real(kind_phys), dimension(nCol,nLay), intent(in) :: &
         cld_frac,  & ! Cloud-fraction from cloud-mp
	 cldtau_lw, & ! In-cloud 10 micron optical depth
         cldtau_sw, & ! In-cloud 0.67 micron optical depth
         cld_liq,   & ! Liquid cloud water mixing ratio (kg/kg)
         cld_ice,   & ! Ice cloud water mixing ratio (kg/kg)
         cld_rain,  & ! Rain cloud water mixing ratio (kg/kg)
         cld_snow,  & ! Snow cloud water mixing ratio (kg/kg)
         cld_graupel  ! Graupel cloud water mixing ratio (kg/kg)
    real(kind_phys), dimension(nCol,nLay), intent(in), optional :: &
         ccld_frac,  & ! Convective cloud fraction
         ccld_liq,   & ! Convective cloud water mixing ratio (kg/kg)
         cld_reliq,  & !
         cld_reice,  & !
         cld_rerain, & !
         cld_resnow    !
    type(cosp_optical_inputs), intent(inout) :: &
         cospIN       ! DDT containing optical inputs needed by COSP.

    ! Locals
    type(rng_state), dimension(nCol) :: rngs
    integer,         dimension(nCol) :: seed
    integer :: i, j, k, iSub, istat
    real(kind_phys), dimension(nCol,nLay) :: cldemis_lw_strat, cldemis_lw_conv
    real(kind_phys), dimension(nCol,nLay) :: cldtau_sw_conv, cldtau_sw_strat
    real(kind_phys), dimension(nCol,nLay) :: column_frac_out, column_prec_out
    real(kind_phys), dimension(nCol,nLay) :: ls_p_rate, cv_p_rate
    real(kind_phys), dimension(nCol,nLay) :: cld_reliq_local, cld_reice_local, cld_rerain_local, cld_resnow_local
    real(kind_phys), dimension(nCol,nLay) :: ccld_frac_local, ccld_liq_local
    real(kind_phys),dimension(:,:),  allocatable  :: frac_ls, prec_ls, frac_cv, prec_cv
    real(kind_phys),dimension(:,:,:),  allocatable  :: frac_prec, &
         MODIS_cloudWater,MODIS_cloudIce, MODIS_watersize,MODIS_iceSize,          &
         MODIS_snowSize,MODIS_cloudSnow,  MODIS_opticalThicknessLiq,              &
         MODIS_opticalThicknessSnow,      MODIS_opticalThicknessIce
    real(kind_phys),dimension(:,:,:,:),  allocatable  :: &
         mr_hydro, Reff

    !
    integer,parameter :: &
         I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
         I_LSCICE = 2, & ! Large-scale (stratiform) ice
         I_LSRAIN = 3, & ! Large-scale (stratiform) rain
         I_LSSNOW = 4, & ! Large-scale (stratiform) snow
         I_CVCLIQ = 5, & ! Convective liquid
         I_CVCICE = 6, & ! Convective ice
         I_CVRAIN = 7, & ! Convective rain
         I_CVSNOW = 8, & ! Convective snow
         I_LSGRPL = 9    ! Large-scale (stratiform) groupel

    ! Stratiform and convective clouds in frac_out (scops output).
    integer, parameter :: &
         I_LSC = 1,    & ! Large-scale clouds
         I_CVC = 2       ! Convective clouds

    ! #####################################################################################
    !
    ! Set default values for optional MP arguments
    !
    ! #####################################################################################
    
    ! We may/maynot have convective cloud-condensate (depends on MP choice).
    ! SCOPS applies sub subsampling to the gridbox, assuming a convective fraction of 5%.
    ! This is satisfactory(?) for the resolution and microphysics used by GCMs for CMIP
    ! experiments, but this is Not scale aware and may cause problems at hi-res.
    ! Need to explore more.
    
    ccld_frac_local(:,:) = 0._kind_phys
    if (present(ccld_frac)) ccld_frac_local = ccld_frac
    ccld_liq_local(:,:) = 0._kind_phys
    if (present(ccld_liq)) ccld_liq_local = ccld_liq

    ! Ditto for hydrometeors sizes.
    cld_reliq_local  = 0._kind_phys
    cld_reice_local  = 0._kind_phys
    cld_rerain_local = 0._kind_phys
    cld_resnow_local = 0._kind_phys
    if (present(cld_reliq))  cld_reliq_local = cld_reliq
    if (present(cld_reice))  cld_reice_local = cld_reice
    if (present(cld_rerain)) cld_rerain_local = cld_rerain
    if (present(cld_resnow)) cld_resnow_local = cld_resnow
    
    ! #####################################################################################
    !
    ! Begin sub-sampling step using SCOPS
    !
    ! #####################################################################################
    if (nSubCol .gt. 1) then
    
       ! RNG used for subcolumn generation
       seed = int(sfcP)
       if (nCol .gt. 1) seed=(sfcP-int(sfcP))*1000000
       call init_rng(rngs, seed)

       ! Call scops
       call scops(nCol, nLay, nSubCol, rngs, cld_frac, ccld_frac_local, overlap, cospIN%frac_out, 0)

       ! Sum up precipitation rates
       ls_p_rate(:,1:nLay) = cld_rain + cld_snow + cld_graupel
       cv_p_rate(:,1:nLay) = 0 

       ! Call Prec_scops
       allocate(frac_prec(nCol, nSubCol, nLay))
       call prec_scops(nCol, nLay, nSubCol, ls_p_rate, cv_p_rate, cospIN%frac_out, frac_prec)
       
       ! ##################################################################################
       ! Compute precipitation fraction in each gridbox
       ! ################################################################################## 
       allocate(frac_ls(nCol, nLay),prec_ls(nCol, nLay), frac_cv(nCol, nLay), &
            prec_cv(nCol, nLay), stat=istat)

       ! Initialize
       frac_ls(1:nCol,1:nLay) =	0._kind_phys
       prec_ls(1:nCol,1:nLay) =	0._kind_phys
       frac_cv(1:nCol,1:nLay) =	0._kind_phys
       prec_cv(1:nCol,1:nLay) = 0._kind_phys
       do j=1,nCol
          do k=1,nLay
             do i=1,nSubCol
                if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._kind_phys
                if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._kind_phys
                if (frac_prec(j,i,k) .eq. 1)         prec_ls(j,k) = prec_ls(j,k)+1._kind_phys
                if (frac_prec(j,i,k) .eq. 2)         prec_cv(j,k) = prec_cv(j,k)+1._kind_phys
                if (frac_prec(j,i,k) .eq. 3)         prec_cv(j,k) = prec_cv(j,k)+1._kind_phys
                if (frac_prec(j,i,k) .eq. 3)         prec_ls(j,k) = prec_ls(j,k)+1._kind_phys
             enddo
             frac_ls(j,k)=frac_ls(j,k)/nSubCol
             frac_cv(j,k)=frac_cv(j,k)/nSubCol
             prec_ls(j,k)=prec_ls(j,k)/nSubCol
             prec_cv(j,k)=prec_cv(j,k)/nSubCol
          enddo
       enddo

       ! ################################################################################## 
       ! Compute mixing ratios, effective radii and precipitation fluxes for clouds
       ! and precipitation
       ! ################################################################################## 
       allocate(mr_hydro(nCol, nSubCol, nLay, nHydro), Reff(nCol, nSubCol, nLay, nHydro))
       mr_hydro(:,:,:,:) = 0._kind_phys
       Reff(:,:,:,:)     = 0._kind_phys
       do iSub=1,nSubCol
          ! Subcolumn clouds
          column_frac_out = cospIN%frac_out(:,iSub,:)

          ! LS clouds
          where (column_frac_out == I_LSC)
             mr_hydro(:,iSub,:,I_LSCLIQ) = cld_liq
             mr_hydro(:,iSub,:,I_LSCICE) = cld_ice
             Reff(:,iSub,:,I_LSCLIQ)     = cld_reliq_local
             Reff(:,iSub,:,I_LSCICE)     = cld_reice_local
          ! CONV clouds
          elsewhere (column_frac_out == I_CVC)
             mr_hydro(:,iSub,:,I_CVCLIQ) = ccld_liq_local
             mr_hydro(:,iSub,:,I_CVCICE) = cld_ice
             Reff(:,iSub,:,I_CVCLIQ)     = cld_reliq_local
             Reff(:,iSub,:,I_CVCICE)     = cld_reice_local
          end where
          ! Subcolumn precipitation
          column_prec_out = frac_prec(:,iSub,:)

          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3) )
             Reff(:,iSub,:,I_LSRAIN) = cld_rerain_local
             Reff(:,iSub,:,I_LSSNOW) = cld_resnow_local
             Reff(:,iSub,:,I_LSGRPL) = cld_resnow_local
          ! CONV precipitation
          elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
             Reff(:,iSub,:,I_CVRAIN) = cld_rerain_local
             Reff(:,iSub,:,I_CVSNOW) = cld_resnow_local
          end where
       enddo

       ! ##################################################################################
       ! Convert the mixing ratio and precipitation fluxes from gridbox mean to
       ! the fraction-based values
       ! ##################################################################################
       do k=1,nLay
          do j=1,nCol
             ! In-cloud mixing ratios.
             if (frac_ls(j,k) .ne. 0._kind_phys) then
                mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
             endif
             if (frac_cv(j,k) .ne. 0._kind_phys) then
                mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
             endif
             
             ! Precipitation
             if (prec_ls(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_LSRAIN) = mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                mr_hydro(j,:,k,I_LSSNOW) = mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                mr_hydro(j,:,k,I_LSGRPL) = mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
             endif
             if (prec_cv(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_CVRAIN) = mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                mr_hydro(j,:,k,I_CVSNOW) = mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
             endif
          enddo
       enddo
       
    else
       cospIN%frac_out = 1
       allocate(mr_hydro(nCol,1,nLay,nHydro),Reff(nCol,1,nLay,nHydro))
       mr_hydro(:,1,:,I_LSCLIQ) = cld_liq
       mr_hydro(:,1,:,I_LSCICE) = cld_ice
       mr_hydro(:,1,:,I_CVCLIQ) = ccld_liq_local
       mr_hydro(:,1,:,I_CVCICE) = cld_ice
       Reff(:,1,:,I_LSRAIN)     = cld_reliq_local
       Reff(:,1,:,I_LSSNOW)   	= cld_resnow_local
       Reff(:,1,:,I_LSGRPL)   	= cld_resnow_local
       Reff(:,1,:,I_CVRAIN)     = cld_rerain_local
       Reff(:,1,:,I_CVSNOW)     = cld_resnow_local
    endif
    
    ! ##################################################################################
    ! 11-micron emissivity (in-cloud), needed by ISCCP simulator.
    ! ##################################################################################
    if (do_isccp) then
       ! Assume same radiative properties for stratiform and convective clouds
       ! True in current RRTMG and RRTMGP implementations
       cldemis_lw_strat = 1._kind_phys - exp(-cldtau_lw)
       cldemis_lw_conv  = cldemis_lw_strat
       !
       cospIN%emiss_11(:,:,:) = 0._kind_phys
       do iSub=1,nSubCol
          where(cospIN%frac_out(:,iSub,:) .eq. 1)
             cospIN%emiss_11(:,iSub,:) = cldemis_lw_conv
          endwhere
          where(cospIN%frac_out(:,iSub,:) .eq. 2)
             cospIN%emiss_11(:,iSub,:) = cldemis_lw_strat
          endwhere
       enddo
    endif

    ! ##################################################################################
    ! 0.67 micron optical-depth (in-cloud), needed by ISCCP, MISR and MODIS simulators.
    ! ##################################################################################
    if (do_isccp .or. do_modis .or. do_misr) then
       ! Assume same radiative properties for stratiform and convective clouds
       ! True in current RRTMG and RRTMGP implementations
       cldtau_sw_strat = cldtau_sw
       cldtau_sw_conv  = cldtau_sw
       !
       cospIN%tau_067(:,:,:) = 0._kind_phys
       do iSub=1,nSubCol
          where(cospIN%frac_out(:,iSub,:) .eq. 1)
             cospIN%tau_067(:,iSub,:) = cldtau_sw_conv
          endwhere
          where(cospIN%frac_out(:,iSub,:) .eq. 2)
             cospIN%tau_067(:,iSub,:) = cldtau_sw_strat
          endwhere
       enddo
    endif

    ! ##################################################################################
    ! MODIS optics
    ! ##################################################################################
    if (do_modis) then
       allocate(MODIS_cloudWater(nCol,nSubCol,nLay),                                   &
                MODIS_cloudIce(nCol,nSubCol,nLay),                                     &
                MODIS_waterSize(nCol,nSubCol,nLay),                                    &
                MODIS_iceSize(nCol,nSubCol,nLay),                                      &
                MODIS_opticalThicknessLiq(nCol,nSubCol,nLay),                          &
                MODIS_opticalThicknessIce(nCol,nSubCol,nLay), stat=istat)

       ! Sample (stratiform/convective) cloud properties.
       ! Liquid cloud-water
       MODIS_cloudWater(:,:,:) = 0._kind_phys
       do iSub=1,nSubCol
          where(cospIN%frac_out(:,iSub,:) .eq. 1)
             MODIS_cloudWater(:,iSub,:) = mr_hydro(:,iSub,:,I_CVCLIQ)
          endwhere
          where(cospIN%frac_out(:,iSub,:) .eq. 2)
             MODIS_cloudWater(:,iSub,:) = mr_hydro(:,iSub,:,I_LSCLIQ)
          endwhere
       enddo
       
       ! Ice cloud-water
       MODIS_cloudIce(:,:,:) = 0._kind_phys
       do iSub=1,nSubCol
          where(cospIN%frac_out(:,iSub,:) .eq. 1)
             MODIS_cloudIce(:,iSub,:) = mr_hydro(:,iSub,:,I_CVCICE)
          endwhere
          where(cospIN%frac_out(:,iSub,:) .eq. 2)
             MODIS_cloudIce(:,iSub,:) = mr_hydro(:,iSub,:,I_LSCICE)
          endwhere
       enddo
       
       ! Liquid cloud particle size
       MODIS_waterSize(:,:,:) = 0._kind_phys
       do iSub=1,nSubCol
          where(cospIN%frac_out(:,iSub,:) .eq. 1)
             MODIS_waterSize(:,iSub,:) = Reff(:,iSub,:,I_CVCLIQ)
          endwhere
          where(cospIN%frac_out(:,iSub,:) .eq. 2)
             MODIS_waterSize(:,iSub,:) = Reff(:,iSub,:,I_LSCLIQ)
          endwhere
       enddo
       
       ! Ice cloud particle size
       MODIS_iceSize(:,:,:) =	0._kind_phys
       do iSub=1,nSubCol
          where(cospIN%frac_out(:,iSub,:) .eq. 1)
             MODIS_iceSize(:,iSub,:) = Reff(:,iSub,:,I_CVCICE)
          endwhere
          where(cospIN%frac_out(:,iSub,:) .eq. 2)
             MODIS_iceSize(:,iSub,:) = Reff(:,iSub,:,I_LSCICE)
          endwhere
       enddo

       ! Partition optical thickness into liquid and ice parts
       call modis_optics_partition(nCol, nLay, nSubCol, MODIS_cloudWater,                &
            MODIS_cloudIce, MODIS_waterSize, MODIS_iceSize, cospIN%tau_067,              &
            MODIS_opticalThicknessLiq, MODIS_opticalThicknessIce)
       
       ! Compute assymetry parameter and single scattering albedo 
       call modis_optics(nCol, nLay, nSubCol, MODIS_opticalThicknessLiq,                 &
            MODIS_waterSize*1.0e6_kind_phys, MODIS_opticalThicknessIce,                  &
            MODIS_iceSize*1.0e6_kind_phys, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)
       
       deallocate (MODIS_cloudWater, MODIS_cloudIce, MODIS_waterSize)
       deallocate (MODIS_iceSize,  MODIS_opticalThicknessLiq)
       deallocate (MODIS_opticalThicknessIce)
       deallocate (mr_hydro, Reff)
    endif

  end subroutine subsample_and_optics

  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  ! ######################################################################################
  subroutine construct_cosp_outputs(do_isccp, do_modis, do_misr, nCol, nSubCol, nLay,    &
  	     Nlvgrid, n_isccp_pres_bins, n_isccp_tau_bins, n_modis_pres_bins,            &
	     n_modis_tau_bins, n_modis_reffl_bins, n_modis_reffi_bins, n_misr_hgt_bins,  &
             n_misr_tau_bins, x)

    ! Inputs
    logical, intent(in) ::   &
         do_isccp,           & ! Flag for COSP ISCCP diagnostics
         do_misr,            & ! Flag for COSP MISR diagnostics
         do_modis              ! Flag for COSP MODIS diagnostics
    integer, intent(in) ::   &
         nCol,               & ! Number of horizontal gridpoints.
	 nSubCol,            & ! Number of COSP subcolumns.
	 nLay,               & ! Number of vertical layers.
	 Nlvgrid,            & ! Number of vertical layers in COSP statistical grid.
         n_isccp_pres_bins,  & ! Number of pressure      bins in ISCCP CFAD.
         n_isccp_tau_bins,   & ! Number of optical-depth bins in ISCCP CFAD.
         n_modis_pres_bins,  & ! Number of pressure      bins in MODIS CFAD.
         n_modis_tau_bins,   & ! Number of optical-depth bins in MODIS CFAD.
         n_modis_reffi_bins, & ! Number of ice-radii     bins in MODIS CFAD.
         n_modis_reffl_bins, & ! Number of liquid-radii  bins in MODIS CFAD
         n_misr_hgt_bins,    & ! Number of height        bins in MISR CFAD.
         n_misr_tau_bins       ! Number of optical-depth bins in MISR CFAD.
    
    ! Outputs
    type(cosp_outputs),intent(out) :: &
         x                    ! COSP output structure  
  
     ! ISCCP simulator outputs
    if (do_isccp) then
       allocate(x%isccp_boxtau(nCol, nSubCol)) 
       allocate(x%isccp_boxptop(nCol, nSubCol))
       allocate(x%isccp_fq(nCol, n_isccp_tau_bins, n_isccp_pres_bins))
       allocate(x%isccp_totalcldarea(nCol))
       allocate(x%isccp_meanptop(nCol))
       allocate(x%isccp_meantaucld(nCol))
       allocate(x%isccp_meantb(nCol))
       allocate(x%isccp_meantbclr(nCol))
       allocate(x%isccp_meanalbedocld(nCol))
    endif

    ! MISR simulator
    if (do_misr) then 
       allocate(x%misr_fq(nCol, n_misr_tau_bins, n_misr_hgt_bins))
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
       allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nCol, n_modis_tau_bins, n_modis_pres_bins))
       allocate(x%modis_Optical_thickness_vs_ReffLIQ(nCol, n_modis_tau_bins, n_modis_reffl_bins))
       allocate(x%modis_Optical_Thickness_vs_ReffICE(nCol, n_modis_tau_bins, n_modis_reffi_bins))
    endif

  end subroutine construct_cosp_outputs

  ! ######################################################################################
  ! SUBROUTINE construct_cospstate
  ! ######################################################################################
  subroutine construct_cospstateIN(nCol, nLay, y)

    ! Inputs
    integer,intent(in) :: &
         nCol, & ! Number of horizontal gridpoints
	 nLay    ! Number of vertical layers
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y
    
    allocate(y%sunlit(nCol),y%skt(nCol),y%land(nCol),y%at(nCol,nLay), y%pfull(nCol,nLay),&
         y%phalf(nCol,nLay+1),y%qv(nCol,nLay), y%hgt_matrix(nCol,nLay),                  &
         y%hgt_matrix_half(nCol,nLay+1))

  end subroutine construct_cospstateIN

  ! ###################################################################################### 
  ! SUBROUTINE construct_cospIN
  ! ######################################################################################
  subroutine construct_cospIN(do_isccp, do_modis, do_misr, nCol, nSubCol, nLay, y)

    ! Inputs
    logical, intent(in) :: &
         do_isccp, & ! Flag for COSP ISCCP diagnostics
         do_misr,  & ! Flag for COSP MISR diagnostics
         do_modis    ! Flag for COSP MODIS diagnostics
    integer,intent(in) :: &
         nCol,     & ! Number of horizontal gridpoints
	 nSubCol,  & ! Number of COSP subcolumns
	 nLay        ! Number of vertical layers
    ! Outputs 
    type(cosp_optical_inputs),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = nCol
    y%Ncolumns = nSubCol
    y%Nlevels  = nLay

    if (do_isccp) then
       if (.not. allocated(y%frac_out)) allocate(y%frac_out(nCol, nSubCol, nLay))
       if (.not. allocated(y%emiss_11)) allocate(y%emiss_11(nCol, nSubCol, nLay))
       if (.not. allocated(y%tau_067))  allocate(y%tau_067( nCol, nSubCol, nLay))
    endif
    if (do_misr) then
       if (.not. allocated(y%tau_067))  allocate(y%tau_067( nCol, nSubCol, nLay))
    endif
    if (do_modis) then
       if (.not. allocated(y%fracLiq))  allocate(y%fracLiq(nCol, nSubCol, nLay))
       if (.not. allocated(y%tau_067))  allocate(y%tau_067(nCol, nSubCol, nLay))
       if (.not. allocated(y%asym))     allocate(y%asym(   nCol, nSubCol, nLay))
       if (.not. allocated(y%ss_alb))   allocate(y%ss_alb( nCol, nSubCol, nLay))
    endif

  end subroutine construct_cospIN
  
  ! ######################################################################################
  ! SUBROUTINE destroy_cosp_outputs
  ! ######################################################################################
  subroutine destroy_cosp_outputs(y)
    type(cosp_outputs),intent(inout) :: y

    ! Deallocate and nullify
    if (associated(y%isccp_totalcldarea))        then
       deallocate(y%isccp_totalcldarea) 
       nullify(y%isccp_totalcldarea)  
    endif
    if (associated(y%isccp_meantb))              then
       deallocate(y%isccp_meantb) 
       nullify(y%isccp_meantb)     
    endif
    if (associated(y%isccp_meantbclr))           then
       deallocate(y%isccp_meantbclr)
       nullify(y%isccp_meantbclr)  
    endif
    if (associated(y%isccp_meanptop))            then
       deallocate(y%isccp_meanptop)
       nullify(y%isccp_meanptop)     
    endif
    if (associated(y%isccp_meantaucld))          then
       deallocate(y%isccp_meantaucld) 
       nullify(y%isccp_meantaucld)       
    endif
    if (associated(y%isccp_meanalbedocld))       then
       deallocate(y%isccp_meanalbedocld)
       nullify(y%isccp_meanalbedocld)     
    endif
    if (associated(y%isccp_boxtau))              then
       deallocate(y%isccp_boxtau)
       nullify(y%isccp_boxtau)       
    endif
    if (associated(y%isccp_boxptop))             then
       deallocate(y%isccp_boxptop)
       nullify(y%isccp_boxptop)     
    endif
    if (associated(y%isccp_fq))                  then
       deallocate(y%isccp_fq)
       nullify(y%isccp_fq)       
    endif
    ! MISR
    if (associated(y%misr_fq))                   then
       deallocate(y%misr_fq) 
       nullify(y%misr_fq)     
    endif
    if (associated(y%misr_dist_model_layertops)) then
       deallocate(y%misr_dist_model_layertops)
       nullify(y%misr_dist_model_layertops)       
    endif
    if (associated(y%misr_meanztop))             then
       deallocate(y%misr_meanztop)
       nullify(y%misr_meanztop)     
    endif
    if (associated(y%misr_cldarea))              then
       deallocate(y%misr_cldarea)
       nullify(y%misr_cldarea)      
    endif

    ! MODIS
    if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
       deallocate(y%modis_Cloud_Fraction_Total_Mean)       
       nullify(y%modis_Cloud_Fraction_Total_Mean)       
    endif
    if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
       deallocate(y%modis_Cloud_Fraction_Ice_Mean)     
       nullify(y%modis_Cloud_Fraction_Ice_Mean)     
    endif
    if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
       deallocate(y%modis_Cloud_Fraction_Water_Mean)           
       nullify(y%modis_Cloud_Fraction_Water_Mean)           
    endif
    if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
       deallocate(y%modis_Cloud_Fraction_High_Mean)     
       nullify(y%modis_Cloud_Fraction_High_Mean)     
    endif
    if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
       deallocate(y%modis_Cloud_Fraction_Mid_Mean)       
       nullify(y%modis_Cloud_Fraction_Mid_Mean)       
    endif
    if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
       deallocate(y%modis_Cloud_Fraction_Low_Mean)     
       nullify(y%modis_Cloud_Fraction_Low_Mean)     
    endif
    if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
       deallocate(y%modis_Optical_Thickness_Total_Mean)  
       nullify(y%modis_Optical_Thickness_Total_Mean)  
    endif
    if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
       deallocate(y%modis_Optical_Thickness_Water_Mean)     
       nullify(y%modis_Optical_Thickness_Water_Mean)     
    endif
    if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
       deallocate(y%modis_Optical_Thickness_Ice_Mean)       
       nullify(y%modis_Optical_Thickness_Ice_Mean)       
    endif
    if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
       deallocate(y%modis_Optical_Thickness_Total_LogMean)    
       nullify(y%modis_Optical_Thickness_Total_LogMean)    
    endif
    if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
       deallocate(y%modis_Optical_Thickness_Water_LogMean)     
       nullify(y%modis_Optical_Thickness_Water_LogMean)     
    endif
    if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
       deallocate(y%modis_Optical_Thickness_Ice_LogMean)     
       nullify(y%modis_Optical_Thickness_Ice_LogMean)     
    endif
    if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
       deallocate(y%modis_Cloud_Particle_Size_Water_Mean)       
       nullify(y%modis_Cloud_Particle_Size_Water_Mean)       
    endif
    if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
       deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)     
       nullify(y%modis_Cloud_Particle_Size_Ice_Mean)     
    endif
    if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
       deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)           
       nullify(y%modis_Cloud_Top_Pressure_Total_Mean)           
    endif
    if (associated(y%modis_Liquid_Water_Path_Mean))                         then
       deallocate(y%modis_Liquid_Water_Path_Mean)     
       nullify(y%modis_Liquid_Water_Path_Mean)     
    endif
    if (associated(y%modis_Ice_Water_Path_Mean))                            then
       deallocate(y%modis_Ice_Water_Path_Mean)       
       nullify(y%modis_Ice_Water_Path_Mean)       
    endif
    if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
       deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
       nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
    endif
    if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
       deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
       nullify(y%modis_Optical_thickness_vs_ReffLIQ)
    endif
    if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
       deallocate(y%modis_Optical_thickness_vs_ReffICE)
       nullify(y%modis_Optical_thickness_vs_ReffICE)
    endif
  end subroutine destroy_cosp_outputs
  
  ! ######################################################################################
  ! SUBROUTINE destroy_cospIN
  ! ######################################################################################
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y
    
    if (allocated(y%tau_067))             deallocate(y%tau_067)
    if (allocated(y%emiss_11))            deallocate(y%emiss_11)
    if (allocated(y%frac_out))            deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))        deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))        deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))                deallocate(y%asym)
    if (allocated(y%ss_alb))              deallocate(y%ss_alb)
    if (allocated(y%fracLiq))             deallocate(y%fracLiq)
    if (allocated(y%fracPrecipIce))       deallocate(y%fracPrecipIce)
  end subroutine destroy_cospIN

  ! ######################################################################################
  ! SUBROUTINE destroy_cospstateIN     
  ! ###################################################################################### 
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%surfelev))        deallocate(y%surfelev)
    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)    
  end subroutine destroy_cospstateIN
end module GFS_cosp
