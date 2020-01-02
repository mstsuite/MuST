      subroutine intphi(lmax,cosrth,r,
     >                  xp,nbnd,
     >                  sumfi)
c     ================================================================
c
      implicit   real*8 (a-h,o-z)
c
      integer    sigma
c
      real*8     cosrth
      real*8     ford(50)
      real*8     xp(3,nbnd)
      real*8     dp2,rhop
      real*8     zero
      real*8     one
      real*8     two
      real*8     half
      real*8     macherr
c
      complex*16 sumfi(0:lmax)
      complex*16 sqrtm1
      complex*16 czero
      complex*16 dummy1,dummy2,exp1,exp2
c
      parameter  (zero=0.0d0)
      parameter  (one=1.0d0)
      parameter  (two=2.0d0)
      parameter  (half=0.5d0)
      parameter  (sqrtm1=(0.0d0,1.0d0))
      parameter  (czero=(0.0d0,0.0d0))
      parameter (pi=3.141592653589793d0)
      parameter (halfpi=half*pi)
      parameter (twopi=two*pi)
      parameter (macherr=1.d-12)
c
      sinrth=sqrt(one-cosrth*cosrth)
      z=r*cosrth
      rho=r*sinrth
c
c     ================================================================
c     find those circle portions within the w-s boundary, and only
c     integrate over phi through them. the value of phi for both
c     ending points of any portion on the boundary are stored in
c     ford.
c     ================================================================
c
      mp=1
      ford(mp)=zero
      do ip=1,nbnd
c
c        =============================================================
c        rp = the distance of the plane from the origin.
c        if r < rp the curve does not touch the plane.
c        if r = rp the curve touches the plane.
c        if r > rp the curve intersects with the plane.
c        mp = the number of intersection points.
c        =============================================================
c
	 rhop=xp(1,ip)*xp(1,ip)+xp(2,ip)*xp(2,ip)
	 dp2=rhop+xp(3,ip)*xp(3,ip)
         if(r*r.gt.dp2.and.rhop.gt.macherr*macherr) then
	    rhop=sqrt(rhop)
	    a=(dp2-xp(3,ip)*z)/(rho*rhop)
            if(abs(a) .le. one+macherr) then
	    if(abs(a).gt.one) a=a/abs(a)
              if(abs(xp(2,ip)) .le. macherr) then
               fi1=halfpi*(one-sign(one,xp(1,ip)))
              else if(xp(2,ip) .ge. zero) then
               fi1=acos(xp(1,ip)/rhop)
              else
               fi1=twopi-acos(xp(1,ip)/rhop)
              endif
	       fi2=acos(a)
               mp=mp+1
               f=fi1-fi2
               if (f.ge.zero) then
                  ford(mp)=f
               else
                  ford(mp)=twopi+f
               end if
               mp=mp+1
               f=fi1+fi2
               if (f.lt.twopi) then
                  ford(mp)=f
               else
                  ford(mp)=f-twopi
               endif
            endif
         endif
      enddo
      mp=mp+1
      ford(mp)=twopi
c
c     ================================================================
c     sort ford so that :   ford(1) < ford(2) < ... < ford(mp).
c     ================================================================
c        
c     ----------------------------------------------------------------
      call sort(mp,ford)
c     ----------------------------------------------------------------
c        write(6,'(/,''  the number of phi = '',i4)')mp
c        write(6,'(  ''  ford:: '',i4,2x,1p1e16.9)')
c    >                             (ip,ford(ip),ip=1,mp)
c
c     ================================================================
c     start the phi integration.
c     ================================================================
c
      do m=0,lmax
         sumfi(m)=czero
      enddo
      exp1=exp(-sqrtm1*ford(1))
      exp2=one
      do ip=1,mp-1
c           phi1=ford(ip)
c           phi2=ford(ip+1)
         phi=ford(ip+1)-ford(ip)
         exp1=exp1*conjg(exp2)
	 exp2=exp(sqrtm1*half*phi)
c exp1=exp(-sqrtm1*half*(phi1+phi2))
	 exp1=exp1*conjg(exp2)
c
c     ================================================================
c     sigma:    = 1       if curve between ford(i) and
c                         ford(i+1) is within the cell.
c               = 0       otherwise.
c
c     i = 1,2,...,mp-1
c     ================================================================
c
         fh=half*(ford(ip)+ford(ip+1))
         x=rho*cos(fh)
         y=rho*sin(fh)
c
c        =============================================================
c        check if the point (x,y,z) is inside or outside the 
c        polyhedron.
c        =============================================================
c
	if(sigma(x,y,z,xp,nbnd,0).eq.1) then
c
         sumfi(0)=sumfi(0)+phi
	 dummy1=one
	 dummy2=one
         do m=1,lmax
	    dummy1=dummy1*exp1
	    dummy2=dummy2*exp2
            sumfi(m) = sumfi(m) + dummy1*dimag(dummy2)
         enddo
	endif  ! sigma=1
      enddo
      do m=1,lmax
	 sumfi(m)=two*sumfi(m)/m
      enddo
c
      return
      end
