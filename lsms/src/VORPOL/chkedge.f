      subroutine chkedge(x0,y0,z0,i,xm,
     >                   xp,dp2,nbnd,isigma)
c     ================================================================
c returns pij in xm(4..6), eij in xm(8..9)
c
      implicit   none
c
      integer    nbnd
      integer    isigma
      integer    i,j,ierr
c
      real*8 x0,y0,z0
      real*8 xm(9),xminv(9)
      real*8 pxy,pyz,pzx
      real*8     eijx
      real*8     eijy
      real*8     eijz
      real*8     xp(3,nbnd)
      real*8     dp2(nbnd)
      real*8     p2mpr
      real*8     pe
      real*8     t
      real*8     tlt
      real*8     tgt
      real*8     zero,one
      real*8     tol
c
      parameter  (zero=0.d0)
      parameter  (one=1.d0)
      parameter  (tol=1.0d-10)
c
         isigma=0
c           ==========================================================
c           look for vector eij first by solving equations:
c
c              ->   ->
c              e  * R = 0
c               i    i
c
c              ->   ->
c              e  * R = 0
c               i    0
c           ==========================================================
            pxy=xp(1,i)*y0-x0*xp(2,i)
            pyz=xp(2,i)*z0-y0*xp(3,i)
            pzx=xp(3,i)*x0-z0*xp(1,i)
            if(abs(pxy) .gt. tol) then
               eijz=one
               eijx=(xp(2,i)*z0-y0*xp(3,i))/pxy
               eijy=(x0*xp(3,i)-xp(1,i)*z0)/pxy
            else if(abs(pyz) .gt. tol) then
               eijx=one
               eijy=(xp(3,i)*x0-z0*xp(1,i))/pyz
               eijz=(y0*xp(1,i)-xp(2,i)*x0)/pyz
            else if(abs(pzx) .gt. tol) then
               eijy=one
               eijz=(xp(1,i)*y0-x0*xp(2,i))/pzx
               eijx=(xp(3,i)*y0-z0*xp(2,i))/pzx
            else
c              =======================================================
c              plane i is parallel to plane j.
c              =======================================================
               return
            endif
c           ==========================================================
c           look for vector pij by solving equations:
c
c              ->   ->   2
c              p  * R = R
c               ij   j   j
c
c              ->   ->   2
c              p  * R = R
c               ij   i   i
c
c              ->   ->
c              p  * e  = 0
c               ij   ij
c           ==========================================================
            xm(4)=xp(1,i)/dp2(i)
            xm(5)=xp(2,i)/dp2(i)
            xm(6)=xp(3,i)/dp2(i)
            xm(7)=eijx
            xm(8)=eijy
            xm(9)=eijz
c           ----------------------------------------------------------
            call invm3(xm,xminv,ierr)
c           ----------------------------------------------------------
            if(ierr.ne.0) then
               write(6,'(''chkedge:  Cannot find pij.'')')
               stop'chkedge'
            endif
            xm(4)=xminv(1)+xminv(2)
            xm(5)=xminv(4)+xminv(5)
            xm(6)=xminv(7)+xminv(8)
c           ==========================================================
c           check if there is any portion of the edge inside the 
c           voronoi polyhedron.....................................
c     ****************************************************************
c     check if the edge is outside the polyhedron, returns isigma:
c
c     isigma = 0, if outside
c            = 1, otherwise
c     ****************************************************************
c
      tlt= 1.0d+30
      tgt=-1.0d+30
      do j=1,nbnd
         p2mpr=dp2(j)-xp(1,j)*xm(4)-xp(2,j)*xm(5)-xp(3,j)*xm(6)
         pe=xp(1,j)*eijx+xp(2,j)*eijy+xp(3,j)*eijz
	 if(abs(pe).le.tol) then
	   if(p2mpr.lt.-tol) return
         else
            t=p2mpr/pe
            if(pe.gt.zero) then
	      if(t.lt.tlt) tlt=t
            else
	      if(t.gt.tgt) tgt=t
            endif
	 endif
      enddo
c
      if(tlt-tgt .ge. tol) isigma=1
c
      return
      end
