c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine polyhedron(xp,dp2,nbnd,
     >                      corn,dc2,ipcorn,ncorn,indxc,
     >                      edge,edgp,ipedge,nedge)
c     ================================================================
c
c     ****************************************************************
c     Look for those corners and edges of the polyhedron, consists of
c     boundary planes given by xp(1), xp(2), xp(3)
c     ****************************************************************
c
      implicit   none
c
      integer    nbnd
      integer    ipcorn,ncorn
      integer    ipedge,nedge
      integer    n
      integer    k
      integer    i
      integer    j
      integer    ierr
      integer    incr
      integer    isigma
c
      real*8     corn(ipcorn,3)
      real*8     dc2(ipcorn)
      integer    indxc(ipcorn)
      real*8     edge(3,ipedge)
      real*8     edgp(3,ipedge)
      real*8     xm(3,3)
      real*8     xminv(3,3)
      real*8     xp(3,nbnd),dp2(nbnd)
      real*8     x0
      real*8     y0
      real*8     z0
      real*8     ptest
      real*8     tol
c
      parameter  (tol=1.0d-10)
      real*8 x,y,z,length2
      length2(x,y,z)=x*x+y*y+z*z
c
c     ================================================================
c     The equation of a plane is:
c         xp * x + yp * y + zp * z = d2,
c     where:
c         d2 = xp**2+yp**2+zp**2
c
c     Look for corners by solving the equations:
c         xp1 * x + yp1 * y + zp1 * z = d12
c         xp2 * x + yp2 * y + zp2 * z = d22
c         xp3 * x + yp3 * y + zp3 * z = d32
c     ================================================================
      do k=1,nbnd
	dp2(k)=length2(xp(1,k),xp(2,k),xp(3,k))
      enddo
      n=0
      do k=1,nbnd-2
         xm(1,1)=xp(1,k)/dp2(k)
         xm(2,1)=xp(2,k)/dp2(k)
         xm(3,1)=xp(3,k)/dp2(k)
         do j=k+1,nbnd-1
            xm(1,2)=xp(1,j)/dp2(j)
            xm(2,2)=xp(2,j)/dp2(j)
            xm(3,2)=xp(3,j)/dp2(j)
            do i=j+1,nbnd
               xm(1,3)=xp(1,i)/dp2(i)
               xm(2,3)=xp(2,i)/dp2(i)
               xm(3,3)=xp(3,i)/dp2(i)
c              -------------------------------------------------------
               call invm3(xm,xminv,ierr)
c              -------------------------------------------------------
               if(ierr.eq.0) then
                  x0=xminv(1,1)+xminv(2,1)+xminv(3,1)
                  y0=xminv(1,2)+xminv(2,2)+xminv(3,2)
                  z0=xminv(1,3)+xminv(2,3)+xminv(3,3)
c                 ----------------------------------------------------
                  call celbnd(x0,y0,z0,
     >                        xp,nbnd,
     >                        corn(1,1),corn(1,2),corn(1,3),n,incr)
c                 ----------------------------------------------------
                  if(incr.eq.1) then
		    if(n.gt.ipcorn) then
		      write(6,'(''polyhedron: n > ipcorn='',i3)')
     &                      ipcorn
		      stop'polyhedron'
		    endif
                     dc2(n)=length2(corn(n,1),corn(n,2),corn(n,3))
                  endif
               endif
            enddo
         enddo
      enddo
      ncorn=n
      if(ncorn.gt.1) then
c        -------------------------------------------------------------
         call sortidx(ncorn,dc2,indxc)
c        -------------------------------------------------------------
      endif
c
c     ================================================================
c     Look for edges. edges are represented by two vectors and one 
c     parameters as:
c          ->       ->
c          p  + t * e
c           n        n
c     where:
c          ->
c          p  = the vector starting from the origin and ending 
c           n   perpendicular to the edge     
c          ->
c          e  = a vector parallel to the edge.
c           n
c          t  = a real parameter.
c     ================================================================
      n=0
      do j=1,nbnd-1
         xm(1,1)=xp(1,j)/dp2(j)
         xm(2,1)=xp(2,j)/dp2(j)
         xm(3,1)=xp(3,j)/dp2(j)
         do 10 i=j+1,nbnd
            ptest = sqrt(dp2(i)*dp2(j))
     >             -abs(xp(1,i)*xp(1,j)+xp(2,i)*xp(2,j)+xp(3,i)*xp(3,j))
            if(ptest.gt.tol) then
c              =======================================================
c              look for vectors eij and pij
c              check if there is any portion of the edge inside the 
c              voronoi polyhedron.....................................
c              -------------------------------------------------------
               call chkedge(xp(1,j),xp(2,j),xp(3,j),i,xm,
     >                      xp,dp2,nbnd,isigma)
c              -------------------------------------------------------
               if(isigma.eq.1) then
                  n=n+1
		  if(n.gt.ipedge) then
		    write(6,'(''polyhedron: n > ipedge='',i3)') ipedge
		    stop'polyhedron'
		  endif
                  edgp(1,n)=xm(1,2)
                  edgp(2,n)=xm(2,2)
                  edgp(3,n)=xm(3,2)
                  edge(1,n)=xm(1,3)
                  edge(2,n)=xm(2,3)
                  edge(3,n)=xm(3,3)
               endif
            endif
10       continue 
      enddo
      nedge=n
c
      if(nbnd.gt.1) then
c        -------------------------------------------------------------
         call sort(nbnd,dp2)
c        -------------------------------------------------------------
      endif

      return 
      end 
