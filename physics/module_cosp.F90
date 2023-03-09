module module_cosp
  use machine, only: kind_phys
  implicit none
  public subsample_and_optics_CAM6, subsample_and_optics_UFS
contains
  ! ######################################################################################
  ! SUBROUTINE subsample_and_optics_CAM6
  ! ######################################################################################
  subroutine subsample_and_optics_CAM6(nCol, nSubCol, nLev, do_isccp, do_misr, do_modis, &
       sfcP, cld_frac, ccld_frac, overlap, cldtau_lw, cldtau_sw, cospIN)
    use mod_scops,      only: scops
    use mod_prec_scops, only: prec_scops
    use mod_rng,        only: rng_state, init_rng
    use mod_cosp,       only: cosp_optical_inputs
    implicit none

    ! Inputs
    logical, intent(in) :: do_isccp, do_misr, do_modis
    integer, intent(in) :: nCol, nSubCol, nLev, overlap
    real(kind_phys), dimension(nCol), intent(in) :: sfcP
    real(kind_phys), dimension(nCol,nLev), intent(in) :: cld_frac, ccld_frac, cldtau_lw, &
         cldtau_sw
    type(cosp_optical_inputs) :: cospIN

    ! Outputs

    ! Locals
    integer :: iSub
    type(rng_state), dimension(nCol)      :: rngs
    integer,         dimension(nCol)      :: seed
    real(kind_phys), dimension(nCol,nLev) :: cldemis_lw_strat, cldemis_lw_conv,          &
         cldtau_sw_strat, cldtau_sw_conv

    ! RNG used for subcolumn generation
    seed = int(sfcP)
    if (nCol .gt. 1) seed=(sfcP-int(sfcP))*1000000 
    call init_rng(rngs, seed)

    ! Call scops
    call scops(nCol, nLev, nSubCol, rngs, cld_frac, ccld_frac, overlap, cospIN%frac_out, 0)
    
    ! Call prec_scops

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

  end subroutine subsample_and_optics_CAM6

  ! ######################################################################################
  ! SUBROUTINE subsample_and_optics_UFS
  ! ######################################################################################
  subroutine subsample_and_optics_UFS(cospIN)
    use mod_cosp,       only: cosp_optical_inputs
    implicit none
    type(cosp_optical_inputs), intent(inout) :: cospIN

  end subroutine subsample_and_optics_UFS

end module module_cosp
