      SUBROUTINE ratintr(n,xa,ya,x,y,dy)
!
      use KindParamModule, only : IntKind, RealKind
      use MathParamModule, only : ZERO
      use ErrorHandlerModule, only : ErrorHandler
!
      INTEGER (kind=IntKind), intent(in) :: n
      INTEGER (kind=IntKind) :: i,m,ns
      INTEGER (kind=IntKind), parameter :: NMAX = 10
!
      REAL (kind=RealKind), intent(in) :: xa(n),ya(n),x
      REAL (kind=RealKind), intent(out) :: dy,y
      REAL (kind=RealKind), parameter :: TINY = 1.0d-25
      REAL (kind=RealKind) :: dd,h,hh,t,w,c(NMAX),d(NMAX)
!
      ns=1
      hh=abs(x-xa(1))
!
      do i=1,n
         h=abs(x-xa(i))
         if (h < TINY)then
            y=ya(i)
            dy=ZERO
            return
         else if (h < hh) then
            ns=i
            hh=h
         endif
         c(i)=ya(i)
         d(i)=ya(i)+TINY
      enddo
!
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            w=c(i+1)-d(i)
            h=xa(i+m)-x
            t=(xa(i)-x)*d(i)/h
            dd=t-c(i+1)
            if(abs(dd) < TINY) then
               call ErrorHandler('ratintr','failure in ratint')
            endif
            dd=w/dd
            d(i)=c(i+1)*dd
            c(i)=t*dd
         enddo
         if (2*ns < n-m)then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
      enddo
!
      END SUBROUTINE ratintr
