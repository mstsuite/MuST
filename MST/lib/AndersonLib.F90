!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dgaleq(a,y,n,ipits,info)
!     ================================================================
      use KindParamModule, only : IntKind, RealKind
!
      use MathParamModule, only : ONE, TEN2m12
!
      use ErrorHandlermodule, only : ErrorHandler, WarningHandler
!
      implicit   none
!
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind), intent(in) :: ipits
      integer (kind=IntKind), intent(out) :: info
      integer (kind=IntKind) :: i
      integer (kind=IntKind) :: j
      integer (kind=IntKind) :: k
      integer (kind=IntKind) :: ij
!
      real (kind=RealKind), intent(inout) :: a(ipits+1,ipits+1)
      real (kind=RealKind), intent(inout) :: y(ipits+1)
      real (kind=RealKind) :: f1
      real (kind=RealKind) :: f2
!
!     ****************************************************************
!     Important note: the following codes may will result a problem:
!             a(i-1,i-1) becomes zero!!!
!     It needs to be fixed.
!     ****************************************************************
      if (n > ipits+1) then
         call ErrorHandler('dgaleq','n > ipits+1',n,ipits+1)
      endif
!
      info = 0
      do i=1,n-1
         if (abs(a(i,i)) < TEN2m12) then
!           call WarningHandler('dgaleq','Ill conditioned a < 10^{-12}',a(i,i))
            info = 1
            return
         endif
      enddo
!
      do i=2,n
         f1=-ONE/a(i-1,i-1)
         do j=i,n
            f2=f1*a(j,i-1)
            do k=1,n
               a(j,k)=a(j,k)+f2*a(i-1,k)
            enddo
            y(j)=y(j)+f2*y(i-1)
         enddo
      enddo
      y(n)=y(n)/a(n,n)
      do ij=1,n-1
         i=n-ij
         do j=1,ij
            y(i)=y(i)-y(i+j)*a(i,i+j)
         enddo
         y(i)=y(i)/a(i,i)
      enddo
!
      end subroutine dgaleq
!     ================================================================
!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function dptc(tcini,tcouti,tcinj,tcoutj,ndim) result(d)
!     ================================================================
      use KindParamModule, only : IntKind, RealKind
! 
      use MathParamModule, only : ZERO
!
      implicit   none
!
      integer (kind=IntKind), intent(in) :: ndim
      integer (kind=IntKind) :: k
!
      real (kind=RealKind) :: d
      real (kind=RealKind), intent(in) :: tcini(ndim)
      real (kind=RealKind), intent(in) :: tcouti(ndim)
      real (kind=RealKind), intent(in) :: tcinj(ndim)
      real (kind=RealKind), intent(in) :: tcoutj(ndim)
!
      d=ZERO
      do k=1,ndim
         d=d+(tcini(k)-tcouti(k))*(tcinj(k)-tcoutj(k))
      enddo
!
      end function dptc
!     ================================================================
