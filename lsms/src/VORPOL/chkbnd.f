      subroutine chkbnd(x0,y0,z0,r02,i0,
     >                  xp,dp2,nbnd,isigma)
c     ================================================================
c
      implicit   none
c
      integer    i0
      integer    nbnd
      integer    isigma
      integer    j
      integer    jp
      integer    k
      integer    kp
      integer    njoint
      integer    ierr
      integer    nbm1
      integer    nbm2
      integer    incr
c
      real*8     x0
      real*8     y0
      real*8     z0
      real*8     r02
      real*8     xp(3,nbnd)
      real*8     dp2(nbnd)
      real*8     xm(9)
      real*8     xminv(9)
      real*8     xc(4)
      real*8     yc(4)
      real*8     zc(4)
c
c     ****************************************************************
c     check if the plane defined by (x0,y0,z0) is an actual boundary
c     plane...........................................................
c     ****************************************************************
c
c     ================================================================
c     The idea is to find the number of joint points with any other
c     two planes not being outside the polyhedron. If it is more than 
c     three, the plane is a boundary plane............................
c     ================================================================
      xm(1)=x0/r02
      xm(2)=y0/r02
      xm(3)=z0/r02
      njoint=0
      nbm2=nbnd-2
      nbm1=nbnd-1
      kp=0
      do while(njoint.le.2 .and. kp.le.nbm2)
         kp=kp+1
         k=mod(kp+i0-1,nbnd)+1
         xm(4)=xp(1,k)/dp2(k)
         xm(5)=xp(2,k)/dp2(k)
         xm(6)=xp(3,k)/dp2(k)
         do jp=kp+1,nbm1
            j=mod(jp+i0-1,nbnd)+1
c           write(6,'(''i0,k,j ='',3i5)')i0,k,j
            xm(7)=xp(1,j)/dp2(j)
            xm(8)=xp(2,j)/dp2(j)
            xm(9)=xp(3,j)/dp2(j)
c           ----------------------------------------------------------
            call invm3(xm,xminv,ierr)
c           ----------------------------------------------------------
            if(ierr.eq.0) then
c              -------------------------------------------------------
               call celbnd(xminv(1)+xminv(2)+xminv(3),
     >                     xminv(4)+xminv(5)+xminv(6),
     >                     xminv(7)+xminv(8)+xminv(9),
     >                     xp,nbnd,
     >                     xc,yc,zc,njoint,incr)
c              -------------------------------------------------------
            endif
         enddo
      enddo
c
      if(njoint.ge.3) then
         isigma=1
      else
         isigma=0
      endif
c
      return
      end
