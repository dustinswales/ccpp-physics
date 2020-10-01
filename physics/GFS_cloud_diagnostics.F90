! ########################################################################################
! This module contains code to produce the UFS High/Mid/Low cloud-diagnostics. 
! This was bundled together with the prognostic cloud modules within the RRTMG implementation.
! For the RRTMGP implementation we propose to keep these diagnostics independent.
! ########################################################################################
module GFS_cloud_diagnostics
  use machine,    only: kind_phys
  use rrtmgp_aux, only: check_error_msg

  ! Module parameters
  integer, parameter :: &
       NK_CLDS     = 3,         & ! Number of cloud vertical domains [low,middle,high]
       ik_lowclds  = 1,         & ! Index for low clouds
       ik_midclds  = 2,         & ! Index for middle clouds
       ik_highclds = 3,         & ! Index for high clouds
       ik_totclds  = NK_CLDS+1, & ! Index for total clouds
       ik_blclds   = NK_CLDS+2    ! Index for boundary-layer clouds
  ! Low/Middle/High pressure boundaries
  real(kind_phys), parameter, dimension(NK_CLDS+1,2) :: &
       ptopc = reshape(source=(/ 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /),  &
                       shape=(/NK_CLDS+1,2/))
  real(kind_phys), parameter :: &
       climit = 0.001,      & ! Lowest allowable cloud-fraction
       ovcst  = 1.0 - 1.0e-8  ! Overcast cloud-fraction 0.999999999 
                      
  ! Version tag and last revision date
  character(40), parameter :: VTAGCLD='UFS-cloud-diagnostics    vX.x May 2020 '
  
  ! Module variables
  integer :: &
       llyr = 2            ! Upper limit of boundary layer clouds
     
  public GFS_cloud_diagnostics_run, GFS_cloud_diagnostics_init,&
       GFS_cloud_diagnostics_finalize, hml_cloud_diagnostics_init
contains
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_cloud_diagnostics_init()  
  end subroutine GFS_cloud_diagnostics_init
  
  ! ######################################################################################
  ! ######################################################################################
!! \section arg_table_GFS_cloud_diagnostics_run
!! \htmlinclude GFS_cloud_diagnostics_run.html
!!  
  subroutine GFS_cloud_diagnostics_run(nCol, nLev, lsswr, lslwr, lat, de_lgth, p_lay,    &
       cld_frac, p_lev, deltaZ, cloud_overlap_param, precip_overlap_param, con_pi,       &
       iovr_lw, iovr_sw, iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr, iovr_exp,        &
       iovr_exprand, ivflip, mbota, mtopa, cldsa, errmsg, errflg)
    implicit none
     
    ! Inputs 
    integer, intent(in) :: &
         nCol,                & ! Number of horizontal grid-points
         nLev,                & ! Number of vertical-layers
         iovr_lw,             & ! Choice of LW cloud-overlap method
         iovr_sw,             & ! Choice of SW cloud-overlap method
         iovr_rand,           & ! Flag for random cloud overlap method
         iovr_maxrand,        & ! Flag for maximum-random cloud overlap method
         iovr_max,            & ! Flag for maximum cloud overlap method
         iovr_dcorr,          & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,            & ! Flag for exponential cloud overlap method
         iovr_exprand,        & ! Flag for exponential-random cloud overlap method  
         ivflip                 ! Flag for vertical grid order (0=toa-2-sfc;1=sfc-2-toa)       
    logical, intent(in) :: &
    	 lsswr,               & ! Call SW radiation?
    	 lslwr                  ! Call LW radiation 
    real(kind_phys), intent(in) :: &
         con_pi                 ! Physical constant: pi  
    real(kind_phys), dimension(nCol), intent(in) :: &
         lat,                 & ! Latitude       
         de_lgth                ! Decorrelation length     
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
         p_lay,               & ! Pressure at model-layer
         cld_frac               ! Total cloud fraction
    real(kind_phys), dimension(nCol,nLev+1), intent(in) :: &
         p_lev                  ! Pressure at model interfaces         
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
    	 deltaZ,              & ! Layer thickness (km)
         cloud_overlap_param, & ! Cloud-overlap parameter
         precip_overlap_param   ! Precipitation overlap parameter
    
    ! Outputs
    character(len=*), intent(out) :: &
         errmsg                 ! Error message
    integer, intent(out) :: &  
         errflg                 ! Error flag
    integer,dimension(nCol,NK_CLDS),intent(out) :: &
         mbota,               & ! Vertical indices for cloud tops (low, middle, high)
         mtopa                  ! Vertical indices for cloud bases (low, middle, high)
    real(kind_phys), dimension(nCol,NK_CLDS+2), intent(out) :: &
         cldsa                  ! Fraction of clouds for low, middle, high (3), total and BL (2)
    
    ! Local variables
    integer i,id,iCol,iLay,icld,iovr
    real(kind_phys) :: tem1
    real(kind_phys),dimension(nCol,NK_CLDS+1) :: ptop1
    real(kind_phys),dimension(nCol) :: rlat
    real(kind_phys),dimension(nCol,nLev) :: cldcnv
	
    if (.not. (lsswr .or. lslwr)) return
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    iovr = max(iovr_lw,iovr_sw)	
    
    ! This is set to zero in all of the progcld() routines and passed to gethml().
    cldcnv(:,:) = 0._kind_phys
    
    do icld = 1, NK_CLDS+1
       tem1 = ptopc(icld,2) - ptopc(icld,1)
       do i=1,nCol
          rlat(i) = abs(lat(i) / con_pi )
          ptop1(i,icld) = ptopc(icld,1) + tem1*max( 0.0, 4.0*rlat(i)-1.0 )
       enddo
    enddo
	
    ! Compute low, mid, high, total, and boundary layer cloud fractions and clouds top/bottom 
    ! layer indices for low, mid, and high clouds. The three cloud domain boundaries are 
    ! defined by ptopc. The cloud overlapping method is defined by control flag 'iovr', which may
    ! be different for lw and sw radiation programs, but the same overlap method is assumed here.
    call gethml(p_lay/100., ptop1, cld_frac, cldcnv, deltaZ, de_lgth, cloud_overlap_param,&
         nCol, nLev, iovr, iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr, iovr_exp,       &
         iovr_exprand, cldsa, mtopa, mbota, errmsg, errflg)	
    
  end subroutine GFS_cloud_diagnostics_run
  
  ! ######################################################################################
  ! ######################################################################################
  subroutine GFS_cloud_diagnostics_finalize()
  end subroutine GFS_cloud_diagnostics_finalize
  
  ! ######################################################################################
  ! Initialization routine for High/Mid/Low cloud diagnostics. 
  ! ######################################################################################
  subroutine hml_cloud_diagnostics_initialize(imp_physics, imp_physics_fer_hires,        &
          imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,                      &
          imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, imp_physics_mg, nLev,        &
          mpi_rank, sigmainit, errflg)
    implicit none
    ! Inputs
    integer, intent(in) :: &
          imp_physics,               & ! Flag for MP scheme
          imp_physics_fer_hires,     & ! Flag for fer-hires scheme
          imp_physics_gfdl,          & ! Flag for gfdl scheme
          imp_physics_thompson,      & ! Flag for thompsonscheme
          imp_physics_wsm6,          & ! Flag for wsm6 scheme
          imp_physics_zhao_carr,     & ! Flag for zhao-carr scheme
          imp_physics_zhao_carr_pdf, & ! Flag for zhao-carr+PDF scheme
          imp_physics_mg               ! Flag for MG scheme
    integer, intent(in) :: &
         nLev,        & ! Number of vertical-layers
         mpi_rank 
    real(kind_phys), dimension(nLev+1), intent(in) :: &
         sigmainit
    ! Outputs
    integer, intent(out) :: &
         errflg
    
    ! Local variables
    integer :: iLay, kl
    
    ! Initialize error flag
    errflg = 0
    
    if (mpi_rank == 0) then
       print *, VTAGCLD      !print out version tag
       print *,' - Using Prognostic Cloud Method'
       if (imp_physics == imp_physics_zhao_carr) then
          print *,'   --- Zhao/Carr/Sundqvist microphysics'
       elseif (imp_physics == imp_physics_zhao_carr_pdf) then
          print *,'   --- zhao/carr/sundqvist + pdf cloud'
       elseif (imp_physics == imp_physics_gfdl) then
          print *,'   --- GFDL Lin cloud microphysics'
       elseif (imp_physics == imp_physics_thompson) then
          print *,'   --- Thompson cloud microphysics'
       elseif (imp_physics == imp_physics_wsm6) then
          print *,'   --- WSM6 cloud microphysics'
       elseif (imp_physics == imp_physics_mg) then
          print *,'   --- MG cloud microphysics'
       elseif (imp_physics == imp_physics_fer_hires) then
          print *,'   --- Ferrier-Aligo cloud microphysics'
       else
          print *,'  !!! ERROR in cloud microphysc specification!!!', &
               '  imp_physics (NP3D) =',imp_physics
          errflg = 1
       endif
    endif
    
    ! Compute the top of BL cld (llyr), which is the topmost non cld(low) layer for 
    ! stratiform (at or above lowest 0.1 of the atmosphere).
    lab_do_k0 : do iLay = nLev, 2, -1
       kl = iLay
       if (sigmainit(iLay) < 0.9e0) exit lab_do_k0
    enddo  lab_do_k0
    llyr = kl      
    
    return
  end subroutine hml_cloud_diagnostics_initialize
  
  ! #########################################################################################
  ! Compute high, mid, low, total, and boundary cloud fractions and cloud top/bottom layer 
  ! indices for model diagnostic output. The three cloud domain boundaries are defined by
  ! ptopc.
  ! #########################################################################################
  subroutine gethml(plyr, ptop1, cldtot, cldcnv, dz, de_lgth, alpha, nCol, nLev, iovr,      &
       iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr, iovr_exp, iovr_exprand,               &
       clds, mtop, mbot, errmsg, errflg)
    implicit none
    
    ! Inputs
    integer, intent(in) :: &
         nCol,         & ! Horizontal dimension
         nLev,         & ! Vertical dimension
         iovr,         & ! Choice of cloud-overlap method
         iovr_rand,    & ! Flag for random cloud overlap method
         iovr_maxrand, & ! Flag for maximum-random cloud overlap method
         iovr_max,     & ! Flag for maximum cloud overlap method
         iovr_dcorr,   & ! Flag for decorrelation-length cloud overlap method
         iovr_exp,     & ! Flag for exponential cloud overlap method
         iovr_exprand    ! Flag for exponential-random cloud overlap method
    real(kind_phys), dimension(nCol,NK_CLDS+1),intent(in) :: &
    	 ptop1           ! ???
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
         plyr,         & ! Pressure at model-layer centers (Pa)
         cldtot,       & ! Layer cloud-fraction
         cldcnv,       & ! Layer convective cloud-fraction
         dz,           & ! Layer thickness (km)
         alpha           ! Cloud overlap parameter
    real(kind_phys), dimension(nCol),   intent(in) :: &
         de_lgth         ! Decorrelation length
    ! Outputs
    integer,dimension(nCol,NK_CLDS),intent(out) :: &
         mbot,         & ! Vertical indices for cloud tops
         mtop            ! Vertical indices for cloud bases
    real(kind_phys), dimension(nCol,NK_CLDS+2), intent(out) :: &
         clds            ! Fraction of clouds for low, middle, high (3), total and BL (2) 
    character(len=*), intent(out) :: &
         errmsg          ! Error message
    integer, intent(out) :: &  
         errflg          ! Error flag
             
    ! Local
    real(kind_phys),dimension(nCol) :: cl1, cl2, dz1
    real(kind_phys) :: pcur, pnxt, ccur, cnxt, alfa
    integer, dimension(nCol):: idom, kbt1, kth1, kbt2, kth2
    integer :: i, k, id, id1, kSFC, kTOA, kInc
    logical :: top_at_1
    
    ! Initialize
    errmsg = ''
    errflg = 0    

    ! Vertical ordering?    
    top_at_1 = (plyr(1,1) .lt.  plyr(1, nLev))
    if (top_at_1) then 
       kSFC = nLev
       kTOA = 1
       kInc = -1
    else
       kSFC = 1
       kTOA = nLev
       kInc = 1
    endif
     
    !
    clds(:,:) = 0.0
    cl1(:)    = 1.0
    cl2(:)    = 1.0
    !
    ! random overlap
    !
    if ( iovr == iovr_rand ) then                     
       do k = kSFC,kTOA,kInc
          do i = 1, nCol
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) cl1(i) = cl1(i) * (1.0 - ccur)
          enddo
          ! save bl cloud
          if (k == llyr) then
             clds(:,ik_blclds) = 1.0 - cl1(:)          
          endif
       enddo
       ! save total cloud
       clds(:,ik_totclds) = 1.0 - cl1(:)              
    !
    ! max/ran overlap
    !
    elseif ( iovr == iovr_maxrand ) then                 
       do k = kSFC,kTOA,kInc
          do i = 1, nCol
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) then             ! cloudy layer
                cl2(i) = min( cl2(i), (1.0 - ccur) )
             else                                ! clear layer
                cl1(i) = cl1(i) * cl2(i)
                cl2(i) = 1.0
             endif
          enddo
          ! save bl cloud
          if (k == llyr) then
             clds(:,ik_blclds) = 1.0 - cl1(:) * cl2(:) 
          endif
       enddo
       ! save total cloud
       clds(:,ik_totclds) = 1.0 - cl1(:) * cl2(:)        
    !
    ! maximum overlap all levels
    ! 
    elseif ( iovr == iovr_max ) then                
       cl1(:) = 0.0
       do k = kSFC,kTOA,kInc
          do i = 1, nCol
             ccur = min( ovcst,  max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) cl1(i) = max( cl1(i), ccur )
          enddo
          ! save bl cloud
          if (k == llyr) then
             clds(:,ik_blclds) = cl1(:)
          endif
       enddo
       ! save total cloud          
       clds(:,ik_totclds) = cl1(:)
    !
    ! random if clear-layer divided, otherwise de-corrlength method
    !
    elseif ( iovr == iovr_dcorr ) then
       dz1(:) = - dz(:,kSFC)       
       do k = kSFC,kTOA,kInc
          do i = 1, nCol
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) then                               ! cloudy layer
                alfa = exp( -0.5*((dz1(i)+dz(i,k)))/de_lgth(i) )
                dz1(i) = dz(i,k)
                cl2(i) =    alfa      * min(cl2(i), (1.0 - ccur)) & ! maximum part
                     + (1.0 - alfa) * (cl2(i) * (1.0 - ccur))       ! random part
             else                                                   ! clear layer
                cl1(i) = cl1(i) * cl2(i)
                cl2(i) = 1.0
                if (k /= 1) dz1(i) = -dz(i,k-1)
             endif
          enddo
          ! save bl cloud
          if (k == llyr) then
             clds(:,ik_blclds) = 1.0 - cl1(:) * cl2(:) 
          endif
       enddo
       ! save total cloud
       clds(:,ik_totclds) = 1.0 - cl1(:) * cl2(:)     
    !
    ! Exponential or exponential-random (alpha is computed differently beforehand)
    !
    elseif ( iovr == iovr_exp .or. iovr == iovr_exprand ) then 
       do k = kSFC,kTOA,kInc
          do i = 1, nCol
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             if (ccur >= climit) then                                   ! cloudy layer
                cl2(i) =   alpha(i,k) * min(cl2(i), (1.0 - ccur))     & ! maximum part
                     + (1.0 - alpha(i,k)) * (cl2(i) * (1.0 - ccur))     ! random part
             else                                                       ! clear layer
                cl1(i) = cl1(i) * cl2(i)
                cl2(i) = 1.0
             endif
          enddo
          ! save bl cloud
          if (k == llyr) then
             clds(:,ik_blclds) = 1.0 - cl1(:) * cl2(:) 
          endif
       enddo
       ! save total cloud
       clds(:,ik_totclds) = 1.0 - cl1(:) * cl2(:)     
    endif
    
    ! Calculte high, mid, low cloud fractions and vertical indices of cloud tops/bases.       
    if (top_at_1) then
       cl1 (1:nCol) = 0.0
       cl2 (1:nCol) = 0.0
       kbt1(1:nCol) = nLev
       kbt2(1:nCol) = nLev
       kth1(1:nCol) = 0
       kth2(1:nCol) = 0
       idom(1:nCol) = 1
       mbot(1:nCol,ik_lowclds)  = nLev
       mtop(1:nCol,ik_lowclds)  = nLev
       mbot(1:nCol,ik_midclds)  = nLev - 1
       mtop(1:nCol,ik_midclds)  = nLev - 1
       mbot(1:nCol,ik_highclds) = nLev - 1
       mtop(1:nCol,ik_highclds) = nLev - 1
       do k = kSFC,kTOA,kInc
          do i = 1, nCol
             id = idom(i)
             id1= id + 1
             
             pcur = plyr(i,k)
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             
             if (k > 1) then
                pnxt = plyr(i,k-1)
                cnxt = min( ovcst, max( cldtot(i,k-1), cldcnv(i,k-1) ))
             else
                pnxt = -1.0
                cnxt = 0.0
             endif
             
             if (pcur < ptop1(i,id1)) then
                id = id + 1
                id1= id1 + 1
                idom(i) = id
             endif
             
             if (ccur >= climit) then
                if (kth2(i) == 0) kbt2(i) = k
                kth2(i) = kth2(i) + 1
                
                if ( iovr == iovr_rand ) then
                   cl2(i) = cl2(i) + ccur - cl2(i)*ccur
                else
                   cl2(i) = max( cl2(i), ccur )
                endif
                
                if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                   kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i) )      &
                        / (cl1(i) + cl2(i)) )
                   kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i) )      &
                        / (cl1(i) + cl2(i)) )
                   cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)
                   
                   kbt2(i) = k - 1
                   kth2(i) = 0
                   cl2 (i) = 0.0
                endif   ! end_if_cnxt_or_pnxt
             endif      ! end_if_ccur
             
             if (pnxt < ptop1(i,id1)) then
                clds(i,id) = cl1(i)
                mtop(i,id) = min( kbt1(i), kbt1(i)-kth1(i)+1 )
                mbot(i,id) = kbt1(i)
                
                cl1 (i) = 0.0
                kbt1(i) = k - 1
                kth1(i) = 0
                
                if (id1 <= NK_CLDS) then
                   mbot(i,id1) = kbt1(i)
                   mtop(i,id1) = kbt1(i)
                endif
             endif     ! end_if_pnxt             
          enddo        ! end_do_i_loop
       enddo           ! end_do_k_loop
    else
       cl1 (1:nCol) = 0.0
       cl2 (1:nCol) = 0.0
       kbt1(1:nCol) = 1
       kbt2(1:nCol) = 1
       kth1(1:nCol) = 0
       kth2(1:nCol) = 0
       idom(1:nCol) = 1
       mbot(1:nCol,ik_lowclds)  = 1
       mtop(1:nCol,ik_lowclds)  = 1
       mbot(1:nCol,ik_midclds)  = 2
       mtop(1:nCol,ik_midclds)  = 2
       mbot(1:nCol,ik_highclds) = 2
       mtop(1:nCol,ik_highclds) = 2
       do k = 1, nLev
          do i = 1, nCol
             id = idom(i)
             id1= id + 1
             
             pcur = plyr(i,k)
             ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
             
             if (k < nLev) then
                pnxt = plyr(i,k+1)
                cnxt = min( ovcst, max( cldtot(i,k+1), cldcnv(i,k+1) ))
             else
                pnxt = -1.0
                cnxt = 0.0
             endif
             
             if (pcur < ptop1(i,id1)) then
                id = id + 1
                id1= id1 + 1
                idom(i) = id
             endif
             
             if (ccur >= climit) then
                if (kth2(i) == 0) kbt2(i) = k
                kth2(i) = kth2(i) + 1
                
                if ( iovr == 0 ) then
                   cl2(i) = cl2(i) + ccur - cl2(i)*ccur
                else
                   cl2(i) = max( cl2(i), ccur )
                endif
                
                if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                   kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
                        / (cl1(i) + cl2(i)) )
                   kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
                        / (cl1(i) + cl2(i)) )
                   cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)
                   
                   kbt2(i) = k + 1
                   kth2(i) = 0
                   cl2 (i) = 0.0
                endif     ! end_if_cnxt_or_pnxt
             endif       ! end_if_ccur
             
             if (pnxt < ptop1(i,id1)) then
                clds(i,id) = cl1(i)
                mtop(i,id) = max( kbt1(i), kbt1(i)+kth1(i)-1 )
                mbot(i,id) = kbt1(i)
                
                cl1 (i) = 0.0
                kbt1(i) = min(k+1, nLev)
                kth1(i) = 0
                
                if (id1 <= NK_CLDS) then
                   mbot(i,id1) = kbt1(i)
                   mtop(i,id1) = kbt1(i)
                endif
             endif     ! end_if_pnxt
          enddo        ! end_do_i_loop
       enddo           ! end_do_k_loop
    endif              ! End vertical ordering

    return
  end subroutine gethml
end module GFS_cloud_diagnostics
