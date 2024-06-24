! ###########################################################################################
!> \file cosp_simulator.F90
!!
!> \defgroup cosp_simulator cosp_simulator.F90
!!
!! \brief This module contains the CCPP interface to the Cloud-Feedback Model Intercomparison
!! Project (CFMIP) Observational Simulator Package Version 2.0 (COSP)
!!
! ###########################################################################################
module cosp_simulator
  use machine,                  only: kind_phys
  use mod_cosp,                 only: cosp_outputs, cosp_optical_inputs, cosp_column_inputs,&
      				      cosp_simulator
  use mod_cosp_modis_interface, only: cosp_modis_init
  use mod_cosp_misr_interface,  only: cosp_misr_init
  use mod_cosp_isccp_interface, only: cosp_isccp_init
  implicit none

contains

! ###########################################################################################
!! \section arg_table_cosp_simulator_init
!! \htmlinclude cosp_simulator_init.html
!!
!!
!> \ingroup cosp_simulator
!!
!! \brief
!!
!! \section cosp_simulator_init
!> @{
! ###########################################################################################
  subroutine cosp_simulator_init(mpirank, mpiroot, do_cosp, do_isccp, do_misr, do_modis,    &
       cosp_nsubcol, imp_physics, imp_physics_thompson, imp_physics_gfdl, isccp_topht,      &
       isccp_topht_dir, top_at_1, iSFC, iTOA, errmsg, errflg)

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
    real(kind_phys), dimension(:,:), intent(in) :: &
         prsi                    ! Pressure at model-interfaces (Pa)

    ! Outputs
    logical, intent(out), :: &
         top_at_1              ! Vertical ordering flag
    integer, intent(out) :: &
         iSFC,               & ! Vertical index for surface
         iTOA                  ! Vertical index for TOA
    character(len=*), intent(out) :: &
         errmsg                ! CCPP error message
    integer, intent(out) :: &
         errflg                ! CCPP error flag

    ! Local
    integer :: nLev, nLay

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Has COSP been requested?
    if (.not. do_cosp) return

    ! What is vertical ordering of the host?
    nLay = size(prsl,2)
    nLev = size(prsi,2)
    top_at_1 = (prsi(1,1) .lt.  prsi(1, nLev))
    if (top_at_1) then
       iSFC = nLay
       iTOA = 1
    else
       iSFC = 1
       iTOA = nLay
    endif

    ! Initialize requested simulators
    if (do_isccp) then
       call cosp_isccp_init(isccp_topht, isccp_topht_dir)
    endif
    if (do_modis) then
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

  end subroutine cosp_simulator_init
!> @}

! ###########################################################################################
!! \section arg_table_cosp_simulator_run
!! \htmlinclude cosp_simulator_run.html
!!
!> \ingroup cosp_simulator
!!
!! \brief
!!
!! \section cosp_simulator_run
!> @{
! ###########################################################################################
  subroutine cosp_simulator_run(nCol, nLay, cosp_nlvgrid, cosp_nsubcol, tsfc, coszen, slmsk,&
       prsl, prsi, phil, phii, tgrs, qgrs, cldtau_lw, cldtau_sw, cld_frac, ccld_frac,       &
       top_at_1, con_g, iSFC, iTOA, n_isccp_pres_bins, isccp_pres_bins, n_isccp_tau_bins,   &
       isccp_tau_bins, n_modis_pres_bins, modis_pres_bins, n_modis_tau_bins, modis_tau_bins,&
       n_misr_pres_bins, misr_pres_bins, n_misr_tau_bins, misr_tau_bins, do_cosp, do_isccp, &
       do_misr, do_modis, overlap, cosp_mp, cosp_mp_ufs,                                    &
       f1isccp_cosp, cldtot_isccp, meancldalb_isccp, meanptop_isccp, meantau_isccp,         &
       meantb_isccp, meantbclr_isccp, tau_isccp, cldptop_isccp, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
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
         n_misr_pres_bins,   & ! Number of pressure      bins in MISR CFAD.
         n_misr_tau_bins,    & ! Number of optical-depth bins in MISR CFAD.
         overlap,            & ! Cloud overlap assumption
         iSFC,               & ! Vertical index for surface
         iTOA,               & ! Vertical index for TOA
         cosp_mp,            & ! Choice of subsampling and optics.
         cosp_mp_ufs           ! Choice of subsampling and optics UFS method.
    real(kind_phys), intent(in) :: &
         con_g                 ! Physical constant: gravitational constant
    real(kind_phys), dimension(:), intent(in) :: & 
         tsfc,               & ! Surface skin temperature (K)
         coszen,             & ! Cosine of SZA
         slmsk,              & ! Area type
	 isccp_pres_bins,    & ! Pressure bin boundaries for ISCCP CFAD.
	 isccp_tau_bins        ! Optical-depth bin boundaries for ISCCP CFAD.
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
         errmsg                ! CCPP error message
    integer, intent(out) :: &
         errflg                ! CCPP error flag
    real(kind_phys), dimension(:,:,:), intent(inout) :: &
         f1isccp_cosp          ! ISCCP CFAD
    real(kind_phys), dimension(:,:), intent(inout) :: &
         tau_isccp,          & ! ISCCP subcolumn optical-depth
         cldptop_isccp         ! ISCCP subcolumn cloud-top pressure
    real(kind_phys), dimension(:), intent(inout) :: &
         cldtot_isccp,       & ! ISCCP mean cloud-fraction
         meancldalb_isccp,   & ! ISCCP mean cloud albedo
         meanptop_isccp,     & ! ISCCP mean cloud-top pressure
         meantau_isccp,      & ! ISCCP mean optical-depth
         meantb_isccp,       & ! ISCCP mean brightness temperature
         meantbclr_isccp       ! ISCCP mean brightness temperature (clear-sky)

    ! Local
    type(cosp_outputs)        :: cospOUT
    type(cosp_optical_inputs) :: cospIN
    type(cosp_column_inputs)  :: cospstateIN
    integer, dimension(nCol)  :: sunlit
    integer :: iCol, nerror, iErr, vs, iprs, itau, iSubCol
    character(len=256),dimension(100) :: cosp_status

    if (.not. do_cosp) return

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Vertical stride direction
    if (top_at_1)       vs = 1
    if (.not. top_at_1) vs = -1
    
    ! Compute sunlit flag.
    sunlit(:) = 0
    do iCol = 1, nCol
       if (coszen(iCol) >= 0.0001) then
          sunlit(iCol) = 1
       endif
    enddo

    ! Type containing COSP outputs.
    call construct_cosp_outputs(do_isccp, do_modis, do_misr, nCol, cosp_nsubcol, nLay,      &
    	 cosp_nlvgrid, n_isccp_pres_bins, n_isccp_tau_bins, n_modis_pres_bins,              &
	 n_modis_tau_bins, n_misr_pres_bins, n_misr_tau_bins, cospOUT)

    ! Host-model state for COSP (toa-2-sfc vertical ordering).
    call construct_cospstateIN(nCol, nLay, cospstateIN)
    cospstateIN%sunlit          = sunlit(:)
    cospstateIN%skt             = tsfc(:)
    cospstateIN%land            = slmsk(:)
    cospstateIN%at              = tgrs(:,iTOA:iSFC:vs)
    cospstateIN%pfull           = prsl(:,iTOA:iSFC:vs)
    cospstateIN%phalf           = prsi(:,iTOA-vs:iSFC:vs)
    cospstateIN%qv              = qgrs(:,iTOA:iSFC:vs,1)
    cospstateIN%hgt_matrix      = phil/con_g
    cospstateIN%hgt_matrix_half = phii/con_g

    ! Derived (optical) inputs for COSP.
    call construct_cospIN(do_isccp, do_modis, do_misr, nCol, cosp_nsubcol, nLay, cospIN)

    !
    ! Call subsample_and_optics
    !
    call subsample_and_optics(nCol, nSubCol, nLay, do_isccp, do_misr, do_modis,          &
    	 prsi(:,iSFC), cld_frac, ccld_frac, overlap, cldtau_lw, cldtau_sw, cospIN)

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
       do iprs=1,n_isccp_pres_bins
          do itau=1,n_isccp_tau_bins
             where(sunlit(1:nCol) .eq. 0)
                cospOUT%isccp_fq(1:nCol,iprs,itau) = R_UNDEF
             end where
          end do
       end do
    endif

    ! Copy COSP outputs to host interstitials.
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
!> @}

  ! ######################################################################################
  ! SUBROUTINE subsample_and_optics
  !
  ! This routine contains the radiaiton to cloud coupling needed by COSP. Changes to host
  ! cloud and radiative configurations need to also occur here.
  !
  ! ######################################################################################
  subroutine subsample_and_optics(nCol, nSubCol, nLay, do_isccp, do_misr, do_modis,      &
       sfcP, cld_frac, ccld_frac, overlap, cldtau_lw, cldtau_sw, cospIN)

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
	 ccld_frac, & ! Convective cloud fraction
	 cldtau_lw, & ! In-cloud 10 micron optical depth
         cldtau_sw    ! In-cloud 0.67 micron optical depth
    type(cosp_optical_inputs), intent(inout) :: &
         cospIN       ! DDT containing optical inputs needed by COSP.

    ! Locals
    type(rng_state), dimension(nCol) :: rngs
    integer,         dimension(nCol) :: seed
    integer :: iSub

    ! RNG used for subcolumn generation
    seed = int(sfcP)
    if (nCol .gt. 1) seed=(sfcP-int(sfcP))*1000000
    call init_rng(rngs, seed)

    ! Call scops
    call scops(nCol, nLay, nSubCol, rngs, cld_frac, ccld_frac, overlap, cospIN%frac_out, 0)

    ! 11-micron emissivity (in-cloud), needed by ISCCP simulator.
    if (do_isccp) then
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

    ! 0.67 micron optical-depth (in-cloud), needed by ISCCP, MISR and MODIS simulators.
    if (do_isccp .or. do_modis .or. do_misr) then
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

  end subroutine subsample_and_optics

  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  ! ######################################################################################
  subroutine construct_cosp_outputs(do_isccp, do_modis, do_misr, nCol, nSubCol, nLay,    &
  	     Nlvgrid, n_isccp_pres_bins, n_isccp_tau_bins, n_modis_pres_bins,            &
	     n_modis_tau_bins, n_misr_pres_bins, n_misr_tau_bins, x)

    ! Inputs
    logical, intent(in) :: &
         do_isccp,          & ! Flag for COSP ISCCP diagnostics
         do_misr,           & ! Flag for COSP MISR diagnostics
         do_modis             ! Flag for COSP MODIS diagnostics
    integer, intent(in) :: &
         nCol,              & ! Number of horizontal gridpoints.
	 nSubCol,           & ! Number of COSP subcolumns.
	 nLay,              & ! Number of vertical layers.
	 Nlvgrid,           & ! Number of vertical layers in COSP statistical grid.
         n_isccp_pres_bins, & ! Number of pressure      bins in ISCCP CFAD.
         n_isccp_tau_bins,  & ! Number of optical-depth bins in ISCCP CFAD.
         n_modis_pres_bins, & ! Number of pressure      bins in MODIS CFAD.
         n_modis_tau_bins,  & ! Number of optical-depth bins in MODIS CFAD.
         n_misr_pres_bins,  & ! Number of pressure      bins in MISR CFAD.
         n_misr_tau_bins      ! Number of optical-depth bins in MISR CFAD.
    
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
       allocate(x%misr_fq(nCol, n_misr_tau_bins, n_modis_pres_bins))
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
       allocate(x%modis_Optical_thickness_vs_ReffLIQ(nCol, n_modis_tau_bins, n_modis_pres_bins))
       allocate(x%modis_Optical_Thickness_vs_ReffICE(nCol, n_modis_tau_bins, n_modis_pres_bins))
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

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

end module cosp_simulator
