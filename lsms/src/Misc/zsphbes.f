      subroutine zsphbes(lmax,x,bj,hpm,dj,dh,job)
c_______________________________________________________________________
c calculate spherical bessel and hankel functions and their derivatives
c input: lmax integer scalar, max l
c        x    complex*16 scalar, argument of the bessel functions
c        job  integer scalar. job=+1: retarded hankel function is calculated
c             job=-1: advanced hankel function is calculated
c returns: bj complex*16 array of (0:lmax), j_l
c          hpm complex*16 array of (0:lmax), h(+-)_l (hankel function)
c          dj complex*16 array of (0:lmax), derivative of bj
c          dh complex*16 array of (0:lmax), derivative of hpm
c calls zsphbesjh
c  computes spherical bessel/hankel functions j and h+-
c    h+ = jl + i nl = sinx/x - i cosx/x
c    h- = jl - i nl = sinx/x + i cosx/x
c
c   xgz ornl 1994
      implicit real*8 (a-h,o-z)
      complex*16 ci
      parameter(ci=(0.0d0,1.0d0))
      complex*16 bj(0:lmax),hpm(0:lmax)
      complex*16 dj(0:lmax),dh(0:lmax)
      complex*16 x,eiz,fll
c
      if(x.eq.0.0d0) goto 200

c zsphbesjh returns bj, hpm and eiz
c eiz = job*exp(job*ci*x)
      call zsphbesjh(lmax,x,bj,hpm,eiz,job)
c put back factor of e^ix/x to hpm
      eiz=eiz/x
      do l=0,lmax
      hpm(l)=eiz*hpm(l)
      enddo

      fll=2.d0/x
      dj(0) = -bj(1)
      dj(1) = bj(0)-fll*bj(1)
      dh(0) = -hpm(1)
      dh(1) = hpm(0)-fll*hpm(1)
c
      do l=2,lmax
	 fll=l/x
         dj(l) = bj(l-1)-fll*bj(l)
         dh(l) = hpm(l-1)-fll*hpm(l)
      enddo
      return

200   write(6,*)'error bessel: zero argument.'
      stop
      end
