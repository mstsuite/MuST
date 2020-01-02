!
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zeroout(x,nx)
!     =================================================================
      use KindParamModule, only : IntKind, RealKind
      use MathParamModule, only : ZERO
!
      implicit   none
!
      integer (kind=IntKind), intent(in) :: nx
      integer (kind=IntKind) :: n
!
      real (kind=RealKind), intent(out) :: x(nx)
!
      do n=1,nx
         x(n)=ZERO
      enddo
!
      end subroutine zeroout
