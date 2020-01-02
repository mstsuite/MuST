      subroutine calsig(lmax,r,
     >                  xp,nbnd,
     >                  tnode,node,wgq,xgq,npts,
     >                  wylm,plm,sumfi)
c     ================================================================
c
c*********************************************************************
c*                                                                   *
c*    This is the new program to expand sigma(r) on complex spherical*
c*    harmonics. It only requires the position of boundary planes.   *
c*                                                                   *
c*    the integration :                                              *
c*                       /                                           *
c*                       |  _          _      *   _                  *
c*                       | do * sigma( r ) * y  ( r ) = wylm         *
c*                       |                    l,m                    *
c*                       /                                           *
c*                                                                   *
c*    for certain r and m => 0, is carried through this subroutine.  *
c*                                                                   *
c*                             m        *                            *
c*    For m < 0, w(l,-m) = (-1) * w(l,m)                             *
c*                                                                   *
c*                               *                                   *
c*                       = w(l,m)         if m = even                *
c*                                                                   *
c*    The array is stored in the following way:                      *
c*                                                                   *
c*        jl       l       m              jl      l      m           *
c*                                                                   *
c*      w( 1 )     0       0            w( 7 )    3      0           *
c*      w( 2 )     1       0            w( 8 )    3      1           *
c*      w( 3 )     1       1            w( 9 )    3      2           *
c*      w( 4 )     2       0            w(10 )    3      3           *
c*      w( 5 )     2       1               .      .      .           *
c*      w( 6 )     2       2               .      .      .           *
c*                                         .      .      .           *
c*                                                                   *
c*    input:                                                         *
c*                                                                   *
c*          lmax     = the maximum value of l.                       *
c*          xgq      = the Gaussian quadrature points.               *
c*          wgq      = the Gaussian quadrature weight associating    *
c*                     each value of xgq.                            *
c*          npts     = the number of Gaussian qurdrature points.     *
c*          r        = the ratial value.                             *
c*          xp       = the vector normal to a boundary plane and     *
c*                     ending on the plane (from origin).            *
c*          nbnd     = the number of boundary planes.                *
c*          tnode    = the nodes along theta integration which is the*
c*                     non-regular points of function: f(theta).     *
c*                     The integration over theta is broken into     *
c*                     several pieces, each of which has two nodes as*
c*                     the terminal points.                          *
c*          node     = the number of tnode.                          *
c*                                                                   *
c*                                                                   *
c*    output:                                                        *
c*                                                                   *
c*          wylm     = step function expansion value for current     *
c*                     value of r.                                   *
c*                                                                   * 
c*********************************************************************
c
      implicit   none
c
      integer    lmax
      integer    nbnd
      integer    node
      integer    npts
      integer    i
      integer    l,m,jl
      integer    n
c
      real*8     wgq(npts)
      real*8     xgq(npts)
      real*8     tnode(node)
      real*8     wt
      real*8     plm((lmax+1)*(lmax+2)/2)
      real*8     xp(3,nbnd)
      real*8     cosrth,offset
      real*8     half 
      real*8     one
      real*8     two
      real*8     wtol
      real*8     d,halfd
      real*8     h0,sh0
      real*8     r,u
      real*8 t(2)
c
      complex*16 sumfi(0:lmax,3)
      complex*16 wylm((lmax+1)*(lmax+2)/2)
      complex*16 af1
      complex*16 czero
c
      parameter  (half=0.5d0)
      parameter  (one=1.0d0)
      parameter  (two=2.0d0)
      parameter  (wtol=1.0d-14)
      parameter  (czero=(0.d0,0.d0))
c
c     ================================================================
c     start to integrate over xn=cos(theta), using gaussian method. 
c
c             2pi              ->   -i*m*phi
c     sumfi = int d(phi) sigma(r )*e        ,   for certain cos(theta)
c              0
c
c     ================================================================
      do i=1,node-1
         d=tnode(i+1)-tnode(i)
         h0=d*1.d-6
	 if(h0.lt.1.d-8) h0=min(1.d-8,d*1.d-3)
         sh0=sqrt(h0)

         call intphi(lmax,tnode(i),r,xp,nbnd,sumfi(0,3))
         call intphi(lmax,tnode(i)+h0,r,xp,nbnd,sumfi(0,1))
         do m=0,lmax
	   sumfi(m,1)=(sumfi(m,1)-sumfi(m,3))/sh0
	 enddo
c
         call intphi(lmax,tnode(i+1),r,xp,nbnd,sumfi(0,3))
         call intphi(lmax,tnode(i+1)-h0,r,xp,nbnd,sumfi(0,2))
         do m=0,lmax
	   sumfi(m,2)=(sumfi(m,2)-sumfi(m,3))/sh0
	 enddo
c
	 offset=half*(tnode(i)+tnode(i+1))
	 halfd=half*d
         do n=1,npts
            cosrth=offset+halfd*xgq(n)
c
c           ==========================================================
c           generate the Legendre functions up to
c           l = lmax for each cos(theta) value.
c           ==========================================================
c           ----------------------------------------------------------
! meis: changed to normalized associated Legendre functions
            call plm_normalized(lmax,cosrth,plm)
c           ----------------------------------------------------------
            call intphi(lmax,cosrth,r,xp,nbnd,sumfi(0,3))
c           ----------------------------------------------------------
c
	    t(1)=sqrt(halfd+halfd*xgq(n))
	    t(2)=sqrt(halfd-halfd*xgq(n))
            wt=wgq(n)*halfd
	    do m=0,lmax
              af1= wt*(sumfi(m,3) -
     >                    sumfi(m,1)*t(1)-sumfi(m,2)*t(2))
	      jl=m*(m+1)/2+1
	    do l=m,lmax
	      jl=jl+l
               wylm(jl) = wylm(jl)+af1*plm(jl)
            enddo
            enddo
c
	    if(xgq(n).ge.0.d0) then
            u=sqrt(d)*xgq(n)
	    wt=two*sqrt(d)*wgq(n)*u*u
c           ----------------------------------------------------------
! meis: changed to normalized associated Legendre functions
            call plm_normalized(lmax,u*u+tnode(i),plm)
c           ----------------------------------------------------------
	    do m=0,lmax
	      af1=wt*sumfi(m,1)
	      jl=m*(m+1)/2+1
	    do l=m,lmax
	      jl=jl+l
               wylm(jl)=wylm(jl)+af1*plm(jl)
            enddo
            enddo
c           ----------------------------------------------------------
! meis: changed to normalized associated Legendre functions
            call plm_normalized(lmax,tnode(i+1)-u*u,plm)
c           ----------------------------------------------------------
	    do m=0,lmax
	      af1=wt*sumfi(m,2)
	      jl=m*(m+1)/2+1
	    do l=m,lmax
	      jl=jl+l
               wylm(jl)=wylm(jl)+af1*plm(jl)
            enddo
            enddo
	    endif
         enddo
      enddo
c
      do jl=1,(lmax+1)*(lmax+2)/2
         if(abs(wylm(jl)) .lt. wtol) then
            wylm(jl)=czero
         endif
      enddo
c
      return
      end
