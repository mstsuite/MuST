!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genGaussianPoints(ng,x1,x2,xg,wg)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : HALF, ONE
!
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ng
   integer (kind=IntKind) :: ig
!
   real (kind=RealKind), intent(in)  :: x1, x2
   real (kind=RealKind), intent(out) :: xg(ng)
   real (kind=RealKind), intent(out) :: wg(ng)
!
   if (ng < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('genGaussianPoints',                          &
                        'invalid number of Gaussian Points',ng)
!     ----------------------------------------------------------------
   else if (x1 > x2) then
!     ----------------------------------------------------------------
      call ErrorHandler('genGaussianPoints',                          &
                        'invalid interval: x1 > x2',x1,x2)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   call gauleg(-ONE,ONE,xg,wg,ng)
!  -------------------------------------------------------------------
!
   do ig = 1, ng
      xg(ig) = HALF*(x2-x1)*xg(ig) + HALF*(x2+x1)
      wg(ig) = HALF*(x2-x1)*wg(ig)
   enddo
!
   end subroutine genGaussianPoints
