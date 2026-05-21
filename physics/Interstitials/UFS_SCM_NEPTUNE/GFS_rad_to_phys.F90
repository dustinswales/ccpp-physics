! #############################################################################################
!>\file GFS_rad_to_phys.F90
!!
!! This module contains the CCPP-compliant codes that fits radiative fluxes and heating rates
!! from a coarse radiation calculation time interval into model's more frequent time steps.
!!
!! Solar heating rates and fluxes are scaled by the ratio of cosine  of zenith angle at the
!! current time to the mean value used in radiation calculation. Surface downward LW flux is
!! scaled by the  ratio of current surface air temperature to the corresponding temperature
!! saved during LW radiation calculation. Upward LW flux at the surface is computed by current
!! ground surface temperature.
!! Surface emissivity effect will be taken in other part of the model.
!!
!! program history:
!!-          198?  nmc mrf    - created, similar as treatment in gfdl radiation treatment
!!-          1994  y. hou     - modified solar zenith angle calculation
!!-     nov  2004  x. wu      - add sfc sw downward flux to the variable list for sea-ice model
!!-     mar  2008  y. hou     - add cosine of zenith angle as output for sunshine duration time calc.
!!-     sep  2008  y. hou     - separate net sw and downward lw in slrad,
!!                              changed the sign of sfc net sw to consistent with
!!                              other parts of the mdl (positive value defines from
!!                              atmos to the ground). rename output fluxes as adjusted
!!                              fluxes. other minor changes such as renaming some of
!!                              passing argument names to be consistent with calling
!!                              program.
!!-     apr  2009  y. hou     - integrated with the new parallel model
!!                              along with other modifications
!!-     mar  2011  y. hou     - minor modification including rearrange
!!                              loop orders and loop structures to improve efficiency
!!-     mar  2014  x. wu      - add sfc nir/vis bm/df to the variable
!!                              list for the coupled model input
!!-     jul  2014  s moorthi  - merge gfs and nems versions
!!-     jun  2014  y. hou     - revised to include both up and down sw
!!                              spectral component fluxes
!!-     Oct  2014  y. hous s. moorthi - add emissivity contribution to
!!                                      upward longwave flux
!!-     Mar  2019  s. moorthi - modify xmu calculation in a time centered
!!                              way and add more accuracy when physics
!!                              time step is close to radiation time step
! #############################################################################################
module GFS_rad_to_phys
  implicit none
  private
  public :: GFS_rad_to_phys_run
contains

!> \section arg_table_GFS_rad_to_phys_run Argument Table
!! \htmlinclude GFS_rad_to_phys_run.html
!!
!!\section GFS_rad_to_phys_general General Algorithm
!> @{
  subroutine GFS_rad_to_phys_run( solhr, slag, sdec, cdec, sinlat, coslat,                     &
       con_g, con_cp, con_pi, con_sbc, xlon, coszen, tsfc_lnd, tsfc_ice, tsfc_wat, tf, tsflw,  &
       tsfc, sfcemis_lnd, sfcemis_ice, sfcemis_wat, sfcdsw, sfcdswc, sfcnsw, sfcdlw, swh, swhc,&
       hlw,hlwc, sfcnirbmu, sfcnirdfu, sfcvisbmu, sfcvisdfu, sfcnirbmd, sfcnirdfd, sfcvisbmd,  &
       sfcvisdfd, im, levs, deltim, fhswr, dry, icy, wet, damp_LW_fluxadj, lfnc_k, lfnc_p0,    &
       use_LW_jacobian, sfculw, use_med_flux, sfculw_med, fluxlwUP_jac, t_lay, p_lay, p_lev,   &
       flux2D_lwUP, flux2D_lwDOWN,pert_radtend,do_sppt,ca_global,tsfc_radtime, dtdt, dtdtnp,   &
       htrlw, adjsfcdsw, adjsfcdswc, adjsfcnsw ,adjsfcdlw, adjsfculw_lnd, adjsfculw_ice,       &
       adjsfculw_wat, xmu, xcosz, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                  &
       adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, errmsg, errflg)
    use machine,         only : kind_phys
    implicit none

    ! Inputs
    integer, intent(in) :: im, levs
    logical, intent(in) :: dry(:), icy(:), wet(:), use_LW_jacobian, damp_LW_fluxadj,           &
         pert_radtend, use_med_flux, do_sppt, ca_global
    real(kind=kind_phys),   intent(in) :: solhr, slag, cdec, sdec, deltim, fhswr, lfnc_k,      &
         lfnc_p0,  con_g, con_cp, con_pi, con_sbc
    real(kind=kind_phys), dimension(:), intent(in) :: sinlat, coslat, xlon, coszen, tf, tsflw, &
         sfcdlw, sfcdsw, sfcdswc, sfcnsw, sfculw, tsfc, tsfc_lnd, tsfc_ice, tsfc_wat,          &
         sfcemis_lnd, sfcemis_ice, sfcemis_wat, sfcnirbmu, sfcnirdfu, sfcvisbmu, sfcvisdfu,    &
         sfcnirbmd, sfcnirdfd, sfcvisbmd, sfcvisdfd
    real(kind=kind_phys), dimension(:,:), intent(in) :: swh, hlw, swhc, hlwc, p_lay, t_lay, p_lev
    real(kind=kind_phys), dimension(:),   intent(in), optional :: sfculw_med, tsfc_radtime
    real(kind=kind_phys), dimension(:,:), intent(in), optional :: flux2D_lwUP, flux2D_lwDOWN,  &
         fluxlwUP_jac

    ! Input/output:
    real(kind=kind_phys), dimension(:,:), intent(inout) :: dtdt
    real(kind=kind_phys), dimension(:,:), intent(inout), optional :: dtdtnp, htrlw
    
    ! Outputs:
    real(kind=kind_phys), dimension(:), intent(out) :: adjsfcdsw, adjsfcnsw, adjsfcdlw, xmu,   &
         xcosz, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu, adjnirbmd, adjnirdfd, adjvisbmd,   &
         adjvisdfd, adjsfcdswc, adjsfculw_lnd, adjsfculw_ice, adjsfculw_wat
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Locals:
    integer :: i, k, nstp, nstl, it, istsun(im),iSFC,iTOA
    real(kind=kind_phys) :: cns,  coszn, tem1, tem2, anginc, rstl, solang, dT, pid12
    real(kind=kind_phys), dimension(im,levs+1) :: flxlwup_adj, flxlwdn_adj
    real(kind=kind_phys) :: fluxlwnet_adj,fluxlwnet,dT_sfc, fluxlwDOWN_jac,lfnc,c1

    ! Parameters
    real(kind=kind_phys), parameter :: L      = 1.  ! Length scale for flux-adjustment scalin
    real(kind=kind_phys), parameter :: gamma  = 0.2 ! Scaling factor for downwelling LW Jacobian profile.
    real(kind=kind_phys), parameter :: f_eps  = 0.0001_kind_phys
    real(kind=kind_phys), parameter :: zero   = 0.0d0, one = 1.0d0
    real(kind=kind_phys), parameter :: hour12 = 12.0_kind_phys
    real(kind=kind_phys), parameter :: f3600  = one/3600.0_kind_phys
    real(kind=kind_phys), parameter :: f7200  = one/7200.0_kind_phys
    real(kind=kind_phys), parameter :: czlimt = 0.0001_kind_phys  

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Vertical ordering?
    if (p_lev(1,1) .lt.  p_lev(1, levs)) then 
       iSFC = levs + 1
       iTOA = 1
    else
       iSFC = 1
       iTOA = levs + 1
    endif

    tem1 = fhswr / deltim
    nstp = max(6, nint(tem1))
    nstl = max(1, nint(nstp/tem1))
    pid12  = con_pi / hour12

    ! -------------------------------------------------------------------------- 
    !  --- ...  sw time-step adjustment for current cosine of zenith angle
    ! --------------------------------------------------------------------------
    if (nstl == 1) then
       cns = pid12 * (solhr + deltim*f7200 - hour12) + slag
       do i = 1, IM
          xcosz(i) = sdec*sinlat(i) + cdec*coslat(i)*cos(cns+xlon(i))
       enddo
    elseif (nstl == nstp) then
       do i = 1, IM
          xcosz(i) = coszen(i)
       enddo
    else
       rstl = one / float(nstl)
       solang = pid12 * (solhr - hour12)         
       anginc = pid12 * deltim * f3600 * rstl
       do i = 1, im
          xcosz(i)  = zero
          istsun(i) = zero
       enddo
       do it=1,nstl
          cns = solang + (float(it)-0.5_kind_phys)*anginc + slag
          do i = 1, IM
             coszn    = sdec*sinlat(i) + cdec*coslat(i)*cos(cns+xlon(i))
             xcosz(i) = xcosz(i) + max(zero, coszn)
             if (coszn > czlimt) istsun(i) = istsun(i) + 1
          enddo
       enddo
       do i = 1, IM
          if (istsun(i) > 0) xcosz(i) = xcosz(i) / istsun(i)  ! mean cosine of solar zenith angle at current time
       enddo
    endif

    do i = 1, im
       ! --------------------------------------------------------------------------
       !> - LW time-step adjustment:
       ! -------------------------------------------------------------------------- 
       tem1 = tf(i) / tsflw(i)
       tem2 = tem1 * tem1
       adjsfcdlw(i) = sfcdlw(i) * tem2 * tem2
       !!  - adjust \a sfc downward LW flux to account for t changes in the lowest model layer.
       !! compute 4th power of the ratio of \c tf in the lowest model layer over the mean value \c tsflw.
       if (dry(i)) then
          tem2 = tsfc_lnd(i) * tsfc_lnd(i)
          adjsfculw_lnd(i) =  sfcemis_lnd(i) * con_sbc * tem2 * tem2 &
                           + (one - sfcemis_lnd(i)) * adjsfcdlw(i)
       endif
       if (icy(i)) then
          tem2 = tsfc_ice(i) * tsfc_ice(i)
          adjsfculw_ice(i) =  sfcemis_ice(i) * con_sbc * tem2 * tem2 &
                           + (one - sfcemis_ice(i)) * adjsfcdlw(i)
       endif
       if (wet(i)) then
          tem2 = tsfc_wat(i) * tsfc_wat(i)
          adjsfculw_wat(i) =  sfcemis_wat(i) * con_sbc * tem2 * tem2 &
                           + (one - sfcemis_wat(i)) * adjsfcdlw(i)
          !>  - replace upward longwave flux provided by the mediator (zero over lakes)
          if (use_med_flux) then
             if (sfculw_med(i) > f_eps) then
                adjsfculw_wat(i) = sfculw_med(i)
             end if
          end if
       endif

       !>  - normalize by average value over radiation period for daytime.
       if ( xcosz(i) > f_eps .and. coszen(i) > f_eps ) then
          xmu(i) = xcosz(i) / coszen(i)
       else
          xmu(i) = zero
       endif

       ! --------------------------------------------------------------------------
       !>  - adjust \a sfc net and downward SW fluxes for zenith angle changes.
       !      note: sfc emiss effect will not be appied here
       ! -------------------------------------------------------------------------- 
       adjsfcnsw(i) = sfcnsw(i)    * xmu(i)
       adjsfcdsw(i) = sfcdsw(i)    * xmu(i)
       adjsfcdswc(i)= sfcdswc(i)   * xmu(i)

       adjnirbmu(i) = sfcnirbmu(i) * xmu(i)
       adjnirdfu(i) = sfcnirdfu(i) * xmu(i)
       adjvisbmu(i) = sfcvisbmu(i) * xmu(i)
       adjvisdfu(i) = sfcvisdfu(i) * xmu(i)

       adjnirbmd(i) = sfcnirbmd(i) * xmu(i)
       adjnirdfd(i) = sfcnirdfd(i) * xmu(i)
       adjvisbmd(i) = sfcvisbmd(i) * xmu(i)
       adjvisdfd(i) = sfcvisdfd(i) * xmu(i)
    enddo

    ! Adjust the LW and SW heating-rates.
    ! For LW, optionally scale using the Jacobian of the upward LW flux. *RRTMGP ONLY*
    ! For SW, adjust heating rates with zenith angle change.
    if (use_LW_jacobian) then
       ! Compute adjusted net LW flux foillowing Hogan and Bozzo 2015 (10.1002/2015MS000455)
       ! Here we assume that the profile of the downwelling LW Jaconiam has the same shape
       ! as the upwelling, but scaled and offset.
       ! The scaling factor is 0.2
       ! The profile of the downwelling Jacobian (J) is offset so that
       !     J_dn_sfc / J_up_sfc = scaling_factor
       !     J_dn_toa / J_up_sfc = 0
       !
       ! Optionally, the flux adjustment can be damped with height using a logistic function
       ! fx ~ L / (1 + exp(-k*dp)), where dp = p - p0
       ! L  = 1, fix scale between 0-1.      - Fixed
       ! k  = 1 / pressure decay length (Pa) - Controlled by namelist
       ! p0 = Transition pressure (Pa)       - Controlled by namelsit
       do i = 1, im
          c1 = fluxlwUP_jac(i,iTOA) / fluxlwUP_jac(i,iSFC)
          dT_sfc = tsfc(i) - tsfc_radtime(i)
          do k = 1, levs
             ! LW net flux
             fluxlwnet = (flux2D_lwUP(i,  k+1) - flux2D_lwUP(i,  k) - &
                          flux2D_lwDOWN(i,k+1) + flux2D_lwDOWN(i,k))
             ! Downward LW Jacobian (Eq. 9)
             fluxlwDOWN_jac = gamma *                                 &
                  (fluxlwUP_jac(i,k)/fluxlwUP_jac(i,iSFC) - c1) /     &
                  (1 - c1)
             ! Adjusted LW net flux(Eq. 10)
             fluxlwnet_adj = fluxlwnet + dT_sfc*                      &
                  (fluxlwUP_jac(i,k)/fluxlwUP_jac(i,iSFC) -           &
                  fluxlwDOWN_jac)
             ! Adjusted LW heating rate
             htrlw(i,k) = fluxlwnet_adj * con_g /                     &
                  (con_cp * (p_lev(i,k+1) - p_lev(i,k)))

             ! Add radiative heating rates to physics heating rate. Optionally, scaled w/ height
             ! using a logistic function
             if (damp_LW_fluxadj) then
                lfnc = L / (1+exp(-(p_lev(i,k) - lfnc_p0)/lfnc_k))
             else
                lfnc = 1.
             endif
             dtdt(i,k) = dtdt(i,k) + swh(i,k)*xmu(i) +                &
                  htrlw(i,k)*lfnc + (1.-lfnc)*hlw(i,k)
          enddo
       enddo
    else
       do k = 1, levs
          do i = 1, im
             dtdt(i,k)  = dtdt(i,k)  + swh(i,k)*xmu(i)  + hlw(i,k)
          enddo
       enddo
    endif

    if (do_sppt .or. ca_global) then
       if (pert_radtend) then
          ! clear sky
          do k = 1, levs
             do i = 1, im
                dtdtnp(i,k) = dtdtnp(i,k) + swhc(i,k)*xmu(i) + hlwc(i,k)  
             enddo
          enddo
       else
          ! all sky
          do k = 1, levs
             do i = 1, im
                dtdtnp(i,k) = dtdtnp(i,k) + swh(i,k)*xmu(i) + hlw(i,k)
             enddo
          enddo
       endif
    endif

    return

  end subroutine GFS_rad_to_phys_run
!> @}
!-----------------------------------
end module GFS_rad_to_phys
