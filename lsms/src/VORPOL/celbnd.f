      subroutine celbnd(x0,y0,z0,
     >                  xp,nbnd,
     >                  x,y,z,n,incr)
c     ================================================================
c
      implicit   none
c
      integer    n
      integer    nbnd
      integer    incr
      integer    m
      integer    sigma
c
      real*8     x0
      real*8     y0
      real*8     z0
      real*8     x(n+1)
      real*8     y(n+1)
      real*8     z(n+1)
      real*8     xp(3,nbnd)
      real*8     tol
c
      parameter  (tol=1.0d-10)
c
      incr=0
      if (sigma(x0,y0,z0,xp,nbnd,1).eq.1) then
         do m=1,n
            if (abs(x0-x(m)).lt.tol .and.
     >          abs(y0-y(m)).lt.tol .and.
     >          abs(z0-z(m)).lt.tol) then
                return
            end if
         end do
         incr=1
         n=n+1
         x(n)=x0
         y(n)=y0
         z(n)=z0
      end if
c
      return
      end
