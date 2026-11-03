!>\file rrtmgp_sampling.F90
!! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2019,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> This module provides a simple implementation of sampling for the
!!   Monte Carlo Independent Pixel Approximation (McICA, doi:10.1029/2002jd003322)
! Cloud optical properties, defined by band and assumed homogenous within each cell (column/layer),
!   are randomly sampled to preserve the mean cloud fraction and one of several possible overlap assumptions
! Users supply random numbers with order ngpt,nlay,ncol
!   These are only accessed if cloud_fraction(icol,ilay) > 0 so many values don't need to be filled in
!
! Adapted by Dustin Swales on 8/11/2020 for use in UFS (NOAA-PSL/CU-CIRES)
!
! -------------------------------------------------------------------------------------------------
module rrtmgp_sampling
  use machine,          only: kind_dbl_prec
  use mersenne_twister, only: random_setseed, random_number, random_stat
  use mo_rte_kind,      only: wp, wl
  use mo_optical_props, only: ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  implicit none
  private
  public :: draw_samples, sampled_mask, sample_cloud_type
contains
!> Apply a T/F sampled cloud mask to cloud optical properties defined by band to produce
!!  McICA-sampled cloud optical properties
  function draw_samples(cloud_mask,do_twostream,clouds,clouds_sampled) result(error_msg)
	! Inputs
    logical, dimension(:,:,:),      intent(in   ) :: cloud_mask     ! Dimensions ncol,nlay,ngpt
    logical,                        intent(in   ) :: do_twostream   ! Do two-stream?
    class(ty_optical_props_arry),   intent(in   ) :: clouds         ! Defined by band
 
 	! Outputs
    class(ty_optical_props_arry),   intent(inout) :: clouds_sampled ! Defined by g-point
    character(len=128)                            :: error_msg

	! Local variables
    integer :: ncol,nlay,nbnd,ngpt
    integer :: imom
    
    error_msg = ""

    ! Array extents
    ncol = clouds%get_ncol()
    nlay = clouds%get_nlay()
    nbnd = clouds%get_nband()
    ngpt = clouds_sampled%get_ngpt()

    ! Optical depth assignment works for 1scl, 2str (also nstr)
    call apply_cloud_mask(ncol,nlay,nbnd,ngpt,clouds_sampled%get_band_lims_gpoint(),cloud_mask,clouds%tau,clouds_sampled%tau)
    !
    ! For 2-stream
    !
    select type(clouds)
    type is (ty_optical_props_2str)
      select type(clouds_sampled)
      type is (ty_optical_props_2str)
         if (do_twostream) then
            call apply_cloud_mask(ncol,nlay,nbnd,ngpt,clouds_sampled%get_band_lims_gpoint(),cloud_mask,clouds%ssa,clouds_sampled%ssa)
            call apply_cloud_mask(ncol,nlay,nbnd,ngpt,clouds_sampled%get_band_lims_gpoint(),cloud_mask,clouds%g,  clouds_sampled%g  )
         endif
      class default
          error_msg = "draw_samples: by-band and sampled cloud properties need to be the same variable type"
      end select
    end select
  end function draw_samples

!> Generate a McICA-sampled cloud mask
  subroutine sampled_mask(randoms, cloud_frac, cloud_mask, overlap_param, randoms2)
    ! Inputs
    real(wp), dimension(:,:,:),  intent(in )           :: randoms       ! ngpt,nlay,ncol
    real(wp), dimension(:,:),    intent(in )           :: cloud_frac    ! ncol,nlay

    ! Outputs
    logical,  dimension(:,:,:),  intent(out)           :: cloud_mask    ! ncol,nlay,ngpt

    ! Inputs (optional)
    real(wp), dimension(:,:),    intent(in ), optional :: overlap_param ! ncol,nlay-1
    real(wp), dimension(:,:,:),  intent(in ), optional :: randoms2      ! ngpt,nlay,ncol
    
    ! Local variables
    integer                              :: ncol, nlay, ngpt, icol, ilay, igpt
    integer                              :: cloud_lay_fst, cloud_lay_lst
    real(wp)                             :: rho
    real(wp), dimension(size(randoms,1)) :: local_rands
    logical,  dimension(size(randoms,2)) :: cloud_mask_layer
    ! NOTE: no initializers on these flags; an initializer implies SAVE, so a flag set by one
    !       call (e.g. dcorr) used to persist into later calls (e.g. max-rand).
    logical                              :: l_use_overlap_param, l_use_second_rng
    character(len=128)                   :: error_msg

    ! Array dimensions
    ncol = size(randoms, 3)
    nlay = size(randoms, 2)
    ngpt = size(randoms, 1)
    
    ! Using cloud-overlap parameter (alpha)?
    l_use_overlap_param = present(overlap_param)
    
    ! Using a second RNG?
    l_use_second_rng    = present(randoms2)
    
    ! Construct the cloud mask for each column
    do icol = 1, ncol
       cloud_mask_layer(1:nlay) = cloud_frac(icol,1:nlay) > 0._wp
       if(.not. any(cloud_mask_layer)) then
          cloud_mask(icol,1:nlay,1:ngpt) = .false.
          cycle
       end if
       cloud_lay_fst = findloc(cloud_mask_layer, .true., dim=1)
       cloud_lay_lst = findloc(cloud_mask_layer, .true., dim=1, back = .true.)
       cloud_mask(icol,1:cloud_lay_fst,1:ngpt) = .false.
       
       ilay = cloud_lay_fst
       local_rands(1:ngpt) = randoms(1:ngpt,cloud_lay_fst,icol)
       cloud_mask(icol,ilay,1:ngpt) = local_rands(1:ngpt) > (1._wp - cloud_frac(icol,ilay))
       do ilay = cloud_lay_fst+1, cloud_lay_lst
          ! ################################################################################
          ! Exponential-decorrelation overlap
          !   new  random deviates if the adjacent layer isn't cloudy
          !   correlated  deviates if the adjacent layer is    cloudy and decorrelation-length
          ! ################################################################################ 
          if (l_use_second_rng) then
             where (randoms2(1:ngpt,ilay,icol) > overlap_param(icol,ilay-1)) &
                  local_rands(1:ngpt) = randoms(1:ngpt,ilay,icol)
             cloud_mask(icol,ilay,1:ngpt) = local_rands(1:ngpt) > (1._wp - cloud_frac(icol,ilay))
          ! ################################################################################
          ! Exponential-random overlap
          !   new  random deviates if the adjacent layer isn't cloudy
          !   correlated  deviates if the adjacent layer is    cloudy
          ! ################################################################################
          elseif (l_use_overlap_param) then
             if(cloud_mask_layer(ilay)) then
                if(cloud_mask_layer(ilay-1)) then
                   ! Create random deviates correlated between this layer and the previous layer
                   !    (have to remove mean value before enforcing correlation).
                   rho = overlap_param(icol,ilay-1)
                   local_rands(1:ngpt) =  rho*(local_rands(1:ngpt)      -0.5_wp) + &
                        sqrt(1._wp-rho*rho)*(randoms(1:ngpt,ilay,icol)-0.5_wp) + 0.5_wp
                else
                   local_rands(1:ngpt) = randoms(1:ngpt,ilay,icol)
                end if
                cloud_mask(icol,ilay,1:ngpt) = local_rands(1:ngpt) > (1._wp - cloud_frac(icol,ilay))
             endif
          ! ################################################################################
          ! Max-random overlap
          !   new  random deviates if the adjacent layer isn't cloudy
          !   same random deviates if the adjacent layer is    cloudy
          ! ################################################################################
          else
             if(cloud_mask_layer(ilay)) then
                if(.not. cloud_mask_layer(ilay-1)) local_rands(1:ngpt) = randoms(1:ngpt,ilay,icol)
                cloud_mask(icol,ilay,1:ngpt) = local_rands(1:ngpt) > (1._wp - cloud_frac(icol,ilay))
             else
                cloud_mask(icol,ilay,1:ngpt) = .false.
             end if
          end if   ! END COND: Maximum-random overlap
       end do    ! END LOOP: Layers
       
       ! Set cloud-mask in layer below clouds to false      
       cloud_mask(icol,cloud_lay_lst+1:nlay, 1:ngpt) = .false.
    end do		! END LOOP: Columns
    
  end subroutine sampled_mask

  ! -------------------------------------------------------------------------------------------------
  !> Generate a McICA cloud mask for a single cloud type (e.g. subgrid convective cloud) with
  !! its own cloud fraction, overlap assumption, and random-number stream.
  !!
  !! Random numbers are drawn from the Mersenne twister seeded with seeds(icol). Callers should
  !! offset the seeds from those used for the stratiform cloud so the two masks are statistically
  !! independent, i.e. the two cloud types overlap randomly within a layer and the expected layer
  !! cover is cf_a + cf_b*(1 - cf_a).
  !!
  !! For decorrelation-length overlap the second set of random numbers continues the same stream
  !! (no re-seeding), so it is independent of the first set.
  ! -------------------------------------------------------------------------------------------------
  function sample_cloud_type(iovr_type, iovr_max, iovr_maxrand, iovr_rand, iovr_dcorr,       &
       iovr_exp, iovr_exprand, seeds, cloud_frac, cloud_mask, overlap_param) result(error_msg)
    integer,                    intent(in )           :: iovr_type     ! Overlap assumption for this cloud type
    integer,                    intent(in )           :: iovr_max, iovr_maxrand, iovr_rand, &
                                                         iovr_dcorr, iovr_exp, iovr_exprand
    integer,  dimension(:),     intent(in )           :: seeds         ! ncol
    real(wp), dimension(:,:),   intent(in )           :: cloud_frac    ! ncol,nlay
    logical,  dimension(:,:,:), intent(out)           :: cloud_mask    ! ncol,nlay,ngpt
    real(wp), dimension(:,:),   intent(in ), optional :: overlap_param ! ncol,nlay-1
    character(len=128)                                :: error_msg

    integer :: ncol, nlay, ngpt, icol, ilay
    type(random_stat) :: rng_stat
    real(kind_dbl_prec), dimension(size(cloud_mask,3))                    :: rng1D
    real(kind_dbl_prec), dimension(size(cloud_mask,3)*size(cloud_mask,2)) :: rng2D
    real(wp), dimension(size(cloud_mask,3),size(cloud_mask,2),size(cloud_mask,1)) :: rng3D, rng3D2

    error_msg = ""
    ncol = size(cloud_mask,1)
    nlay = size(cloud_mask,2)
    ngpt = size(cloud_mask,3)
    cloud_mask = .false.

    if (iovr_type == iovr_dcorr .or. iovr_type == iovr_exp .or. iovr_type == iovr_exprand) then
       if (.not. present(overlap_param)) then
          error_msg = "sample_cloud_type: overlap_param required for dcorr/exp/exprand overlap"
          return
       endif
    endif

    ! Skip columns with no cloud of this type
    if (.not. any(cloud_frac > 0._wp)) return

    do icol = 1, ncol
       call random_setseed(seeds(icol), rng_stat)
       if (iovr_type == iovr_max) then
          call random_number(rng1D, rng_stat)
          do ilay = 1, nlay
             rng3D(:,ilay,icol) = real(rng1D, wp)
          enddo
       else
          do ilay = 1, nlay
             call random_number(rng1D, rng_stat)
             rng3D(:,ilay,icol) = real(rng1D, wp)
          enddo
       endif
       if (iovr_type == iovr_dcorr) then
          ! Continue the stream: independent of rng3D
          call random_number(rng2D, rng_stat)
          rng3D2(:,:,icol) = reshape(real(rng2D, wp), shape=[ngpt,nlay])
       endif
    enddo

    if (iovr_type == iovr_maxrand .or. iovr_type == iovr_rand .or. iovr_type == iovr_max) then
       call sampled_mask(rng3D, cloud_frac, cloud_mask)
    elseif (iovr_type == iovr_dcorr) then
       call sampled_mask(rng3D, cloud_frac, cloud_mask, overlap_param=overlap_param, randoms2=rng3D2)
    elseif (iovr_type == iovr_exp .or. iovr_type == iovr_exprand) then
       call sampled_mask(rng3D, cloud_frac, cloud_mask, overlap_param=overlap_param)
    else
       error_msg = "sample_cloud_type: unsupported overlap assumption"
    endif
  end function sample_cloud_type

!> Apply a true/false cloud mask to a homogeneous field
!!   This could be a kernel
  subroutine apply_cloud_mask(ncol,nlay,nbnd,ngpt,band_lims_gpt,cloud_mask,input_field,sampled_field)
    integer,                                intent(in ) :: ncol,nlay,nbnd,ngpt
    integer,     dimension(2,nbnd),         intent(in ) :: band_lims_gpt
    logical,     dimension(ncol,nlay,ngpt), intent(in ) :: cloud_mask
    real(wp),    dimension(ncol,nlay,nbnd), intent(in ) :: input_field
    real(wp),    dimension(ncol,nlay,ngpt), intent(out) :: sampled_field

    integer :: icol,ilay,ibnd,igpt

    do ibnd = 1, nbnd
      do igpt = band_lims_gpt(1,ibnd), band_lims_gpt(2,ibnd)
        do ilay = 1, nlay
          sampled_field(1:ncol,ilay,igpt) = merge(input_field(1:ncol,ilay,ibnd), 0._wp, cloud_mask(1:ncol,ilay,igpt))
        end do
      end do
    end do
  end subroutine apply_cloud_mask

end module rrtmgp_sampling
