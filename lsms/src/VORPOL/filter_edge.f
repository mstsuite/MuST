      subroutine filter_edge(x0,y0,z0,r02,i0,
     >                       xp,dp2,nbnd,isigma,i)
c     ================================================================
c
      implicit   none
c
      integer    i0
      integer    nbnd
      integer    isigma
      integer    i
      integer    ip
c
      real*8     x0
      real*8     y0
      real*8     z0
      real*8     r02
      real*8     xp(3,nbnd)
      real*8     dp2(nbnd)
      real*8     xm(9)
      real*8     ptest
      real*8     tol
c
      parameter  (tol=1.0d-10)
c
c     ****************************************************************
c     check if the plane defined by (x0,y0,z0) is an possible boundary
c     plane...........................................................
c
c     The idea is to determine the edge formed by this plane and any
c     other one is partially inside the polyhedron...................
c     ****************************************************************
c
c     ================================================================
c     Look for edges. edges are represented by two vectors and one 
c     parameters as:
c          ->       ->
c          p  + t * e
c           n        n
c     where:
c          ->
c          p  = the vector starting from the origin and endind at and
c           n   perpendicular to the edge     
c          ->
c          e  = a vector parallel to the edge.
c           n
c          t  = a real parameter.
c     ================================================================
c
      xm(1)=x0/r02
      xm(2)=y0/r02
      xm(3)=z0/r02
      do ip=1,nbnd-1
	 i=mod(ip+i0-1,nbnd)+1
         ptest = sqrt(dp2(i)*r02)-abs(xp(1,i)*x0+xp(2,i)*y0+xp(3,i)*z0)
         if (ptest.gt.tol) then
c           ==========================================================
c           check if there is any portion of the edge inside the 
c           voronoi polyhedron.....................................
c           ----------------------------------------------------------
            call chkedge(x0,y0,z0,i,xm,
     >                   xp,dp2,nbnd,isigma)
c           ----------------------------------------------------------
            if(isigma.eq.1) return
         endif
      enddo
c
      return
      end
