module LUdcmpModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ONE, ZERO, CONE, CZERO
!
public :: LUdcmp
    interface LUdcmp
       module procedure LUdcmp_r, LUdcmp_c
    end interface
!
private
!
contains
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    SUBROUTINE LUdcmp_r(a,n,np,indx,d)
!   ==================================================================
    implicit none
!
    INTEGER (kind=IntKind), intent(in) :: n,np
    INTEGER (kind=IntKind), intent(out) :: indx(n)
    integer (kind=IntKind), parameter :: NMAX = 500
    INTEGER (kind=IntKind) :: i,imax,j,k
!
    REAL (kind=RealKind), intent(inout) :: a(np,np)
    REAL (kind=RealKind), intent(out) :: d
    REAL (kind=RealKind), parameter :: TINY=1.0d-20
    REAL (kind=RealKind) :: aamax,dum,sum,vv(NMAX)
!
    d=ONE
    do i=1,n
       aamax=ZERO
       do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       enddo
       if (aamax.eq.ZERO) stop 'Singular matrix in LUdcmp'
       vv(i)=ONE/aamax
    enddo
    do j=1,n
       do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
       enddo
       aamax=ZERO
       do i=j,n
          sum=a(i,j)
          do k=1,j-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
             imax=i
             aamax=dum
          endif
       enddo
       if (j.ne.imax)then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if (a(j,j).eq.ZERO)a(j,j)=TINY
       if (j.ne.n)then
          dum=ONE/a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          enddo
       endif
    enddo
    END subroutine LUdcmp_r
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    SUBROUTINE LUdcmp_c(a,n,np,indx,d)
!   ==================================================================
    implicit none
!
    INTEGER (kind=IntKind), intent(in) :: n,np
    INTEGER (kind=IntKind), intent(out) :: indx(n)
    integer (kind=IntKind), parameter :: NMAX = 500
    INTEGER (kind=IntKind) :: i,imax,j,k
!
    complex (kind=CmplxKind), intent(out) :: d
    REAL (kind=RealKind), parameter :: TINY=1.0d-40
    real (kind=RealKind) :: aamax,vv(NMAX)
!
    complex (kind=CmplxKind), intent(inout) :: a(np,np)
    complex (kind=CmplxKind) :: dum,sum
!
    d=ONE
    do i=1,n
       aamax=ZERO
       do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       enddo
       if (aamax.eq.ZERO) stop 'Singular matrix in LUdcmp'
       vv(i)=ONE/aamax
    enddo
    do j=1,n
       do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
       enddo
       aamax=ZERO
       do i=j,n
          sum=a(i,j)
          do k=1,j-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (abs(dum).ge.aamax) then
             imax=i
             aamax=abs(dum)
          endif
       enddo
       if (j.ne.imax)then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if (abs(a(j,j)).eq.ZERO)a(j,j)=TINY
       if (j.ne.n)then
          dum=CONE/a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          enddo
       endif
    enddo
    END subroutine LUdcmp_c
!
end module LUdcmpModule
