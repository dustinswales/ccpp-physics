!>\file ozinterp.f90
!! This file contains ozone climatology interpolation subroutines.

!>\ingroup mod_GFS_phys_time_vary
!! This module contains subroutines of reading and interpolating ozone coefficients.
module ozinterp
  
  implicit none

  private

  public :: read_o3data, setindxoz, ozinterpol

contains
  !
  !**********************************************************************
  !
  SUBROUTINE read_o3data (ntoz, kozpl, latsozp, levozp, timeoz,oz_coeff,&
       oz_lat, oz_pres, oz_time, ozplin)
    use machine, only: kind_phys
    
    ! Inputs
    integer, intent(in) :: ntoz, kozpl, latsozp, levozp, timeoz, oz_coeff

    ! In/outs
    real(kind_phys), dimension(latsozp), intent(inout) :: oz_lat(latsozp),&
         oz_pres(levozp), oz_time(timeoz), ozplin(latsozp,levozp,oz_coeff,timeoz)
  
    ! Locals
    integer :: i, n, k, a1, a2, a3, a4
    real(kind=4), dimension(latsozp)  :: oz_lat4, tempin
    real(kind=4), dimension(levozp)   :: oz_pres4
    real(kind=4), dimension(timeoz+1) :: oz_time4

    if (ntoz <= 0) return

    ! Open file
    open(unit=kozpl,file='global_o3prdlos.f77', form='unformatted', convert='big_endian')

    ! Read in indices and axis
    read (kozpl) a1, a2, a3, a4, oz_lat4, oz_pres4, oz_time4

    ! Store
    oz_pres(:) = oz_pres4(:)
    oz_pres(:) = log(100.0*oz_pres(:)) ! from mb to ln(Pa) 
    oz_lat(:)  = oz_lat4(:)
    oz_time(:) = oz_time4(:)

    ! Read in ozplin which is in order of (latitudes, ozone levels, coeff number, time)
    ! assume latitudes is on a uniform gaussian grid
    DO i=1,timeoz
       DO n=1,oz_coeff
          DO k=1,levozp
             READ(kozpl) tempin
             ozplin(:,k,n,i) = tempin(:)
          ENDDO
       ENDDO
    ENDDO

    ! Close file
    close(kozpl)

  END SUBROUTINE read_o3data
  !
  !**********************************************************************
  !
  SUBROUTINE setindxoz(npts, jo3, oz_lat, dlat, jindx1, jindx2, ddy)
    USE MACHINE,  ONLY: kind_phys
    implicit none

    ! Inputs
    integer,         intent(in)  :: npts, jo3
    real(kind_phys), intent(in)  :: oz_lat(jo3), dlat(npts)
    ! Outputs
    integer,         intent(out) :: JINDX1(npts),JINDX2(npts)
    real(kind_phys), intent(out) :: DDY(npts)
    ! Local
    integer i,j,lat
    
    DO J=1,npts
       jindx2(j) = jo3 + 1
       do i=1,jo3
          if (dlat(j) < oz_lat(i)) then
             jindx2(j) = i
             exit
          endif
       enddo
       jindx1(j) = max(jindx2(j)-1,1)
       jindx2(j) = min(jindx2(j),jo3)
       if (jindx2(j) .ne. jindx1(j)) then
          DDY(j) = (dlat(j)           - oz_lat(jindx1(j))) &
                 / (oz_lat(jindx2(j)) - oz_lat(jindx1(j)))
       else
          ddy(j) = 1.0
       endif
!       print *,' j=',j,' dlat=',dlat(j),' jindx12=',jindx1(j), &
!         jjindx2(j),' oz_lat=',oz_lat(jindx1(j)),              &
!         oz_lat(jindx2(j)),' ddy=',ddy(j)
    ENDDO
 
    RETURN
  END SUBROUTINE setindxoz
  !
  !**********************************************************************
  !
  SUBROUTINE ozinterpol(me, npts, IDATE, fhour, jindx1, jindx2, latsozp,&
       levozp, oz_coeff, timeoz, ozplin, oz_time, oz_pres, oz_lat,      &
       ddy, ozplout)
    USE MACHINE,  ONLY : kind_phys
    implicit none

    ! Inputs
    integer,         intent(in)  :: me, idate(4), latsozp, levozp,      &
         oz_coeff, timeoz, JINDX1(npts), JINDX2(npts)
    real(kind_phys), intent(in)  :: fhour,  ddy(npts), oz_lat(latsozp), &
         oz_pres(levozp), oz_time(timeoz), ozplin(latsozp,levozp,oz_coeff,timeoz)

    ! Outputs
    real(kind_phys), intent(out) :: ozplout(npts,levozp,oz_coeff)

    ! Local
    integer :: IDAT(8),JDAT(8),iday,j,j1,j2,l,npts,nc,n1,n2,jdow,jdoy,&
         jday,w3kindreal,w3kindint
    real(kind_phys) :: tem, tx1, tx2, rjday
    real(8) :: rinc(5)
    real(4) :: rinc4(5)

    IDAT=0
    IDAT(1)=IDATE(4)
    IDAT(2)=IDATE(2)
    IDAT(3)=IDATE(3)
    IDAT(5)=IDATE(1)
    RINC=0.
    RINC(2)=FHOUR
    call w3kind(w3kindreal,w3kindint)
    if(w3kindreal==4) then
       rinc4=rinc
       CALL W3MOVDAT(RINC4,IDAT,JDAT)
    else
       CALL W3MOVDAT(RINC,IDAT,JDAT)
    endif
    !
    jdow = 0
    jdoy = 0
    jday = 0
    call w3doxdat(jdat,jdow,jdoy,jday)
    rjday = jdoy + jdat(5) / 24.
    IF (RJDAY < oz_time(1)) RJDAY = RJDAY + 365.
    !
    n2 = timeoz + 1
    do j=2,timeoz
       if (rjday < oz_time(j)) then
          n2 = j
          exit
       endif
    enddo
    n1 = n2 - 1
    !
    !     if (me == 0) print *,' n1=',n1,' n2=',n2,' rjday=',rjday
    !    &,'oz_time=',oz_time(n1),oz_time(n2)
    !

    tx1 = (oz_time(n2) - rjday) / (oz_time(n2) - oz_time(n1))
    tx2 = 1.0 - tx1
    
    if (n2 > timeoz) n2 = n2 - timeoz
    !
    do nc=1,oz_coeff
       DO L=1,levozp
          DO J=1,npts
             J1  = JINDX1(J)
             J2  = JINDX2(J)
             TEM = 1.0 - ddy(J)
             ozplout(j,L,nc) = & 
                  tx1*(TEM*ozplin(J1,L,nc,n1)+ddy(J)*ozplin(J2,L,nc,n1)) & 
                  + tx2*(TEM*ozplin(J1,L,nc,n2)+ddy(J)*ozplin(J2,L,nc,n2))
          ENDDO
       ENDDO
    enddo
    !
    RETURN
  END SUBROUTINE ozinterpol
  
end module ozinterp
