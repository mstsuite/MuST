!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine beschb(x,gam1,gam2,gampl,gammi)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : one
   use MathParamModule, only : eight
!
   implicit none
!
   integer (kind=IntKind), parameter :: NUSE1=7
   integer (kind=IntKind), parameter :: NUSE2=8
!
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind), intent(out) :: gam1,gam2,gammi,gampl
!
   real (kind=RealKind) :: xx,c1(7),c2(8),chebev
   SAVE c1,c2
   data c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,    &
           -3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
   data c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,  &
           -4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
   xx=eight*x*x-one
   gam1=chebev(-one,one,c1,NUSE1,xx)
   gam2=chebev(-one,one,c2,NUSE2,xx)
   gampl=gam2-x*gam1
   gammi=gam2+x*gam1
!
   end subroutine beschb
