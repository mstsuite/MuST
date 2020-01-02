      subroutine interf(plmofz,
     >eta,rm,lmax,ndlm,mofj,lofj,clm,xgs,wgs,ngauss_mad,rdif,dqintl)
c     integrate gaussian from xi0 to infinity
c     to get(one/r- error function/r) subtract one/r to get 
c     potential from -gaussian a distance rdif from origin on z-axis
c     integrate over cos(theta) at a distance rm from origin to get dqintl
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
      integer j
      integer iz
      integer ngauss_mad
c
      real*8 erfc
      real*8 clm(ndlm)
      real*8 plmofz(ndlm)
      real*8 fnpi
      real*8 rm
c
      real*8 xgs(ngauss_mad)
      real*8 wgs(ngauss_mad)
c
      real*8 dqintl(0:lmax)
      real*8 r
      real*8 eta
      real*8 xi0
      real*8 sum
      real*8 half
      real*8 one
      real*8 rho2
      real*8 z
      real*8 rdif
      real*8 pi
      real*8 twopi
      real*8 two
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.0d0)
c    -------------------------------------------------------------------
      pi=fnpi()
      twopi=two*pi
      call zeroout(dqintl(0),lmax+1)
      xi0 = sqrt(eta)*half
c         loop over cos(theta)
      do iz=1,ngauss_mad
        z=xgs(iz)
c           sin(theta)**2
            rho2=one-z*z
c           actual x*x+y*y
            rho2=rm*rm*rho2
c           distance from sphere to gaussian center at z=rdif
c           actual z determined by cos(theta) and rm
	    r=sqrt((rdif-rm*z)**2+rho2)
c           calculate potential at this distance
c    -------------------------------------------------------------------
c erfc is defined as
c  2/sqrt(pi) int_x^infty exp(-x*x)dx
c   --- xgz
              sum=(erfc(r*xi0)-one)/r
c             sum contribution of potential to particular L
          sum=wgs(iz)*sum*twopi
c             it doesn't matter if we use ylm or cylm because m=0
c             we could save time with a legengre polynomial routine for m=0
c     get legendre polynomials on z=cos(theta) grid
! meis: changed to normalized associated Legendre functions
        call plm_normalized(lmax,z,plmofz)
c     loop over L
        do j=1,ndlm
          m=mofj(j)
c       only m=0 contributes
          if(m.eq.0)then
            l=lofj(j)
            dqintl(l) = dqintl(l)+clm(j)*plmofz(j)*sum
          endif
        enddo
      enddo
c
      return
c
      end
