      subroutine zsphbesjh(lmax,x,bj,hm,eiz,job)
c_______________________________________________________________________
c calculate spherical bessel and hankel functions
c input: lmax integer scalar, max l
c        x    complex*16 scalar, argument of the bessel functions
c        job  integer scalar. job=+1: retarded hankel function is calculated
c             job=-1: advanced hankel function is calculated
c returns: bj complex*16 array of (0:lmax), j_l
c          hm complex*16 array of (0:lmax), h(+)_l (retarded hankel function)
c             this is scaled out of a factor of exp(i*x)/x
c          eiz complex*16 scalar, factor exp(i*x)
c calls zsphbesj
c  computes spherical bessel/hankel functions j and h+
c    h+ = jl + i nl = sinx/x - i cosx/x
c    h- = jl - i nl = sinx/x + i cosx/x
c
c   xgz ornl 1994
      implicit real*8 (a-h,o-z)
      complex*16 ci
      parameter(ci=(0.0d0,1.0d0))
      complex*16 bj(0:lmax),hm(0:lmax)
      complex*16 x,x1,eiz,fll,ac,as
c
      if(x.eq.0.0d0) goto 200

c use hm for scratch in zsphbesj
c zsphbesj returns x1=1/x, as=sin(x),ac=cos(x)
      call zsphbesj(lmax,x,as,ac,x1,bj,hm)
c pull out factor of e^ix/x from h
c eiz = job*exp(job*ci*x)
      eiz=job*exp(job*ci*x)

c now use hm for hankel hm = jl + i nl
      hm(0)=-ci
c hm(1) = -1.d0 - ci*x1
      hm(1)=-dcmplx(1.0d0-dimag(x1),dreal(x1))

      do l=2,lmax
        hm(l)=(2*l-1)*hm(l-1)*x1-hm(l-2)
      enddo
      return

200   write(6,*)'error bessel: zero argument.'
      stop
      end
