      function sigma(x0,y0,z0,xp,nbnd,mode)
c     ================================================================
c
      implicit   none
c
      integer    nbnd
      integer    sigma
      integer    mode
      integer    i
c
      real*8     x0
      real*8     y0
      real*8     z0
      real*8     xp(3,nbnd)
      real*8     tol
c
      parameter  (tol=1.0d-10) 
c
      do i=1,nbnd
         if((x0-xp(1,i))*xp(1,i)+(y0-xp(2,i))*xp(2,i)+
     &      (z0-xp(3,i))*xp(3,i)
     &      .gt.tol*mode) then
            sigma=0
            return
         end if
      end do
c
      sigma=1
      return
      end
