module DerivativeModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
public :: derv5
   interface derv5
      module procedure derv5r, derv5c
   end interface
!
private
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine derv5r(y,dy,x,nn)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: nn

   integer (kind=IntKind) :: i,i1

   real (kind=RealKind), intent(in) :: x(nn)
   real (kind=RealKind), intent(in) :: y(nn)
   real (kind=RealKind), intent(out) :: dy(nn)

   real (kind=RealKind) :: a,b,c,d,e

!  ===================================================================
!  lagrangian fit of y=f(x)
!  derivative is then calculated
!  ===================================================================
   do i=1,nn
      i1=max(i-2,1)
      i1=min(i1,nn-4)
      a=y(i1)
      b=(y(i1+1)-a)/(x(i1+1)-x(i1))
      c=(y(i1+2)-a)/(x(i1+2)-x(i1))
      c=(c-b)/(x(i1+2)-x(i1+1))
      d=(y(i1+3)-a)/(x(i1+3)-x(i1))
      d=(d-b)/(x(i1+3)-x(i1+1))
      d=(d-c)/(x(i1+3)-x(i1+2))
      e=(y(i1+4)-a)/(x(i1+4)-x(i1))
      e=(e-b)/(x(i1+4)-x(i1+1))
      e=(e-c)/(x(i1+4)-x(i1+2))
      e=(e-d)/(x(i1+4)-x(i1+3))
      dy(i)=b+c*(  (x(i)-x(i1)) + (x(i)-x(i1+1))  )                 &
             +d*(  (x(i)-x(i1)) * (x(i)-x(i1+1))                    &
                 + (x(i)-x(i1+1)) * (x(i)-x(i1+2))                  &
                 + (x(i)-x(i1+2)) * (x(i)-x(i1))  )                 &
             +e*(  (x(i)-x(i1)) * (x(i)-x(i1+1)) * (x(i)-x(i1+2))   &
                 + (x(i)-x(i1+3)) * (x(i)-x(i1+1)) * (x(i)-x(i1+2)) &
                 + (x(i)-x(i1)) * (x(i)-x(i1+1)) * (x(i)-x(i1+3))   &
                 + (x(i)-x(i1)) * (x(i)-x(i1+2)) * (x(i)-x(i1+3)) )
   enddo
   dy(1)=dy(2)
   end subroutine derv5r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine derv5c(y,dy,x,nn)
!  ===================================================================
   implicit   none

   integer (kind=IntKind), intent(in) :: nn
   integer (kind=IntKind) :: i,i1

   real (kind=RealKind), intent(in) :: x(nn)

   complex (kind=CmplxKind), intent(in) :: y(nn)
   complex (kind=CmplxKind), intent(out) :: dy(nn)
   complex (kind=CmplxKind) :: a,b,c,d,e

!  ===================================================================
!  lagrangian fit of y=f(x)
!  derivative is then calculated
!  ===================================================================
   do i=1,nn
      i1=max(i-2,1)
      i1=min(i1,nn-4)
      a=y(i1)
      b=(y(i1+1)-a)/(x(i1+1)-x(i1))
      c=(y(i1+2)-a)/(x(i1+2)-x(i1))
      c=(c-b)/(x(i1+2)-x(i1+1))
      d=(y(i1+3)-a)/(x(i1+3)-x(i1))
      d=(d-b)/(x(i1+3)-x(i1+1))
      d=(d-c)/(x(i1+3)-x(i1+2))
      e=(y(i1+4)-a)/(x(i1+4)-x(i1))
      e=(e-b)/(x(i1+4)-x(i1+1))
      e=(e-c)/(x(i1+4)-x(i1+2))
      e=(e-d)/(x(i1+4)-x(i1+3))
      dy(i)=b+c*(  (x(i)-x(i1)) + (x(i)-x(i1+1))  )                 &
             +d*(  (x(i)-x(i1)) * (x(i)-x(i1+1))                    &
                 + (x(i)-x(i1+1)) * (x(i)-x(i1+2))                  &
                 + (x(i)-x(i1+2)) * (x(i)-x(i1))  )                 &
             +e*(  (x(i)-x(i1)) * (x(i)-x(i1+1)) * (x(i)-x(i1+2))   &
                 + (x(i)-x(i1+3)) * (x(i)-x(i1+1)) * (x(i)-x(i1+2)) &
                 + (x(i)-x(i1)) * (x(i)-x(i1+1)) * (x(i)-x(i1+3))   &
                 + (x(i)-x(i1)) * (x(i)-x(i1+2)) * (x(i)-x(i1+3)) )
   enddo
   dy(1)=dy(2)
   end subroutine derv5c
!  ===================================================================
end module DerivativeModule
