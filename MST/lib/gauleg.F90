!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine gauleg(x1,x2,x,w,n)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : zero, fourth, half, one, two, three, pi
   implicit   none
!
   integer (kind=IntKind), intent(in) :: n
!
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
!
   real (kind=RealKind), intent(in) :: x1
   real (kind=RealKind), intent(in) :: x2
   real (kind=RealKind), intent(out) :: x(n)
   real (kind=RealKind), intent(out) :: w(n)
!
   real (kind=RealKind) :: xm
   real (kind=RealKind) :: xl
   real (kind=RealKind) :: z
   real (kind=RealKind) :: z1
   real (kind=RealKind) :: p1
   real (kind=RealKind) :: p2
   real (kind=RealKind) :: p3
   real (kind=RealKind) :: pp
   real (kind=RealKind), parameter :: eps = 3.0d-14
!
!  *******************************************************************
!  generates gaussian points and weights to be used in ...............
!  Gauss-Legendre quadrature..........................................
!  See:: Numerical Recipes, page 125..................................
!  *******************************************************************
!
   m=(n+1)/2
   xm=half*(x2+x1)
   xl=half*(x2-x1)
   do i=1,m
      z=cos(pi*(i-fourth)/(n+half))
      LOOP_do: do
         p1=one
         p2=zero
         do j=1,n
            p3=p2
            p2=p1
            p1=((two*j-one)*z*p2-(j-one)*p3)/j
         enddo
         pp=n*(z*p1-p2)/(z*z-one)
         z1=z
         z=z1-p1/pp
         if(abs(z-z1) <= eps) then
            exit LOOP_do
         endif
      enddo LOOP_do
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=two*xl/((one-z*z)*pp*pp)
      w(n+1-i)=w(i)
   enddo
!
   end subroutine gauleg
