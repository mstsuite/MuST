c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stepyll(r,rmt,wylm,
     >                   lmax,clm,
     >                   ngauss,xgq,wgq,
     >                   xp,xpt,nbnd,
     >                   rcs2,
     >                   edge,edgp,nedge,
     >                   np_index,rp_index,tnode,plm,sumfi)
c     ================================================================
c
      implicit   none
c
      integer    lmax
      integer    ngauss
      integer    nbnd
      integer    nedge
      integer    np_index(2)
      integer    ip
      integer    np
      integer    node
      integer    it
      integer    jl
      integer    l
      integer    sigma
c
      real*8     r
      real*8     xgq(ngauss)
      real*8     wgq(ngauss)
      real*8     xp(3,nbnd)
      real*8     xpt(3,nbnd)
      real*8     rmt
      real*8     rcs2
      real*8     edge(3,nedge)
      real*8     edgp(3,nedge)
      real*8     rp_index(2)
      real*8     tnode(nbnd+2*nedge+2)
      real*8     clm((lmax+1)*(lmax+2)/2)
      real*8     plm((lmax+1)*(lmax+2)/2)
      real*8     pi
      parameter (pi=3.141592653589793d0)
      real*8     begth
      real*8     endth
      real*8     one
      real*8     two,twopi
      parameter  (two=2.0d0)
      parameter  (twopi=two*pi)
      real*8     zero
      real*8     tol
c
      complex*16 wylm((lmax+1)*(lmax+2)/2),sumfi(0:lmax,3)
c
      parameter  (one=1.0d0)
      parameter  (zero=0.0d0)
      parameter  (tol=1.0d-10)
c
c     ****************************************************************
c     this subroutine calculates the expansion coefficients ,        *
c                                       _                            *
c     sigma (r),  of stepfunction sigma(r) on complex spherical      *
c          l                                                         *
c                m _                                                 *
c     harmonics Y (r)                                                *
c                l                                                   *
c                                                                    *
c                 ___                                                *
c            _    \                   m _                            *
c     sigma( r ) = >   sigma ( r ) * Y (r)                           *
c                 /__       l         l                              *
c                 l,l'                                               *
c                                                                    *
c     it calls calsig for calculating the expansion coefficients     *
c     of sigma ( r ) through complex spherical harmonics.            *
c             l                                                      *
c     yang wang:: dept. of physics, FAU, Boca raton, FL 33431        *
c                    ver 2.1 written in Aug. 1992                    *
c     ****************************************************************
c
c     ----------------------------------------------------------------
      call zeroout(wylm,(lmax+1)*(lmax+2))
c     ----------------------------------------------------------------
c
c     ================================================================
c     special case of r being inside the muffin-tin radius............
c     ================================================================
      if ( r .le. rmt ) then
         wylm(1)=two*sqrt(pi)
         return
      endif
c
c     ================================================================
c     special case of r being outside the poyhedron...................
c     ================================================================
      if (r*r .ge. rcs2) then
         return
      end if
c
c     ================================================================
c     case of r being greater than the muffin-tin radius..............
c     ----------------------------------------------------------------
      if(sigma(zero,zero,r,xpt,nbnd,1).eq.1) then
         endth=+one
      else
         endth=two
      endif
c     ----------------------------------------------------------------
      if(sigma(zero,zero,-r,xpt,nbnd,1).eq.1) then
         begth=-one
      else
         begth=-two
      endif
      ip=0
      do np=1,nbnd
         if( (np_index(1).eq.np .and. r.lt.rp_index(1)) .or. 
     >       (np_index(2).eq.np .and. r.lt.rp_index(2)) ) then
            if (xpt(3,np).gt.zero .and. xpt(3,np).le.r) then
               endth=xpt(3,np)/r
            else if(xpt(3,np).lt.zero .and. -xpt(3,np).le.r) then
               begth=xpt(3,np)/r
            endif
         else if((np_index(1).ne.np) .and. (np_index(2).ne.np)) then
c           ==========================================================
c           only store those boundaries planes not perpendicular to 
c           the z-axis.
c           ==========================================================
            ip=ip+1
            xp(1,ip)=xpt(1,np)
            xp(2,ip)=xpt(2,np)
            xp(3,ip)=xpt(3,np)
         endif
      enddo
c     ----------------------------------------------------------------
      call caltnode(r,xpt,nbnd,
     >              edge,edgp,nedge,
     >              begth,endth,
     >              np_index,tnode,node)
c     ----------------------------------------------------------------
c
      if(abs(tnode(1)+one) .le. tol) then
c        -------------------------------------------------------------
         call intpl0(-one,tnode(2),lmax,plm)
c        -------------------------------------------------------------
         do l=0,lmax
            jl=(l+1)*(l+2)/2-l
            wylm(jl)=wylm(jl)+twopi*plm(l+1)
         enddo
         do it=1,node-1
            tnode(it)=tnode(it+1)
         enddo
         node=node-1
      endif
c
      if(abs(tnode(node)-one) .le. tol.and.node.gt.1) then
c        -------------------------------------------------------------
         call intpl0(tnode(node-1),one,lmax,plm)
c        -------------------------------------------------------------
         do l=0,lmax
            jl=(l+1)*(l+2)/2-l
            wylm(jl)=wylm(jl)+twopi*plm(l+1)
         enddo
         node=node-1
      endif
c
c     ================================================================
c     it calls calsig to calculate  sigma   ( r ) ,
c                                        l,m
c     which are the expansion coefficients of step function on the 
c     complex spherical harmonics. the data are stored in wylm(jl).
c     ----------------------------------------------------------------
      call calsig(lmax,r,xp,ip,
     >            tnode,node,wgq,xgq,ngauss,
     >            wylm,plm,sumfi)
c     ----------------------------------------------------------------
      do jl=1,(lmax+1)*(lmax+2)/2
	wylm(jl)=wylm(jl)*clm(jl)
      enddo
c
c     ================================================================
c
      return
      end
