c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine intpl0(x1,x2,lmax,yl)
c     ================================================================
c
c     ****************************************************************
c                                  x2       0
c     calculate the integration:  inte {dx P (x)}
c                                  x1       l
c
c     and store in intpl0.
c     ****************************************************************
c
      implicit   none
c
      integer    lmax
      integer    l
c
      real*8     x1
      real*8     x2
      real*8     p0x1
      real*8     p0x2
      real*8     p1x1
      real*8     p1x2
      real*8     p2x1
      real*8     p2x2
      real*8     yl(0:lmax)
      real*8     one, pi
c
      parameter  (one=1.0d0)
      parameter (pi=3.141592653589793d0)
c
      p0x1=one
      p0x2=one
      p1x1=x1
      p1x2=x2
      yl(0)=x2-x1
c
      do l=2,lmax+1
         p2x1=((2*l-1)*x1*p1x1-(l-1)*p0x1)/dble(l)
         p2x2=((2*l-1)*x2*p1x2-(l-1)*p0x2)/dble(l)
         yl(l-1)=(p2x2-p2x1-x2*p1x2+x1*p1x1)/dble(l-1)
         p0x1=p1x1
         p0x2=p1x2
         p1x1=p2x1
         p1x2=p2x2
      enddo
c
! meis: transform to normalized Legendre functions
      do l=0,lmax
        yl(l)=yl(l)*sqrt(dble(2*l+1)/(4.0d0*pi))
      end do

      return
      end
