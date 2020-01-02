c     ====================================================================
      subroutine interfsmr(plmofz,
     >                     eta,rm,lmax,ndlm,mofj,lofj,clm,
     >                     xgs,wgs,ngauss_mad,ngaussr_mad,rdif,dqintl)
c     ====================================================================
c
c     ********************************************************************
c     integrate gaussian from 0 to  xi0
c     to get( error function)  divid by r to get
c     potential from gaussian adistance rdif from origin on z-axis
c     integrate over cos(theta) at a distance rm from origin to get dqintl
c     ********************************************************************
c
      implicit none
c
      integer lmax
c
      integer ndlm
      integer lofj(ndlm)
      integer mofj(ndlm)
      integer l
      integer m
      integer n
      integer j
      integer iz
      integer ngauss_mad
      integer ngaussr_mad
c
      real*8 clm(ndlm)
      real*8 plmofz(ndlm)
      real*8 fnpi
      real*8 rm
c
      real*8 xgs(ngauss_mad)
      real*8 wgs(ngauss_mad)
c
      real*8 dqintl(0:lmax)
      real*8 rssq
      real*8 eta
      real*8 xi
      real*8 xisq
      real*8 xi0
      real*8 arg
      real*8 sum
      real*8 zero
      real*8 half
      real*8 one
      real*8 rho2
      real*8 z
      real*8 rdif
      real*8 pi
      real*8 twopi
      real*8 two
      real*8 xg(ngaussr_mad),wg(ngaussr_mad)
      parameter (zero=0.0d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.0d0)
c    -------------------------------------------------------------------
      pi=fnpi()
      twopi=two*pi
      call zeroout(dqintl(0),lmax+1)
      xi0 = sqrt(eta)*half
!      call gauleg(zero,xi0,xg,wg,ngaussr_mad)
      call gauss_legendre_points(zero,xi0,xg,wg,ngaussr_mad)
c         loop over cos(theta)
      do iz=1,ngauss_mad
        z=xgs(iz)
            rho2=one-z*z
c           actual x*x+y*y
            rho2=rm*rm*rho2
c           distance from sphere to gaussian center at z=rdif
c           actual z determined by cos(theta) and rm
            rssq=(rdif-rm*z)**2+rho2
c           calculate potential at this distance
            sum=zero
            do n=1,ngaussr_mad
              xi=xg(n)
              xisq = xi * xi
              arg = - rssq*xisq
              sum = sum+exp(arg)*wg(n)
            enddo
c           potential for negative gaussian
c Don't divide by r because the integral was rescaled from (0,r*xi0)
c to (0,xi0) thus missing a factor of r.
            sum=-sum*two/sqrt(pi)
	    sum=wgs(iz)*sum*twopi
c     get legendre polynomials on z=cos(theta) grid
! meis: changed to normalized associated Legendre functions
        call plm_normalized(lmax,z,plmofz)
c    -------------------------------------------------------------------
c     loop over L
        do j=1,ndlm
          m=mofj(j)
c       only m=0 contributes
          if(m.eq.0)then
            l=lofj(j)
c           potential for negative gaussian plus delta funtion
            dqintl(l) = dqintl(l)+clm(j)*plmofz(j)*sum
          endif
        enddo
      enddo
      return
c
      end
