c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine csbf(n,p,r,fb,fn,fh)
c     ================================================================
c
c     ****************************************************************
c     spherical bessel and neuman functions
c     ****************************************************************
c
      implicit none
c
      integer     n
      integer     l
c 
      real*8      r
c
      complex*16  fb(0:n)
      complex*16  fn(0:n)
      complex*16  fh(0:n)
      complex*16  p
      complex*16  x1
      complex*16  x2
      complex*16  sqrtm1
c
      parameter (sqrtm1=(0.d0,1.d0))
c
      x1=r*p
      x2=x1*x1
c
      fb(0)= sin(x1)/x1
      fb(1)= sin(x1)/x2 - cos(x1)/x1
      fn(0)=-cos(x1)/x1
      fn(1)=-cos(x1)/x2 - sin(x1)/x1
      fh(0)=fb(0)+sqrtm1*fn(0)
      fh(1)=fb(1)+sqrtm1*fn(1)
c
      do l=2,n
        fb(l)=(2*l-1)*fb(l-1)/x1-fb(l-2)
      end do
      do l=2,n
        fn(l)=(2*l-1)*fn(l-1)/x1-fn(l-2)
      end do
      do l=2,n
        fh(l)=fb(l)+sqrtm1*fn(l)
      end do
c
      return
      end
