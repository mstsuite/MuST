c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zeroout(x,nx)
c     =================================================================
c
      implicit   none
c
      integer    nx,n
c
      real*8     x(nx)
      real*8     zero
c
      parameter  (zero=0.0)
c
      do n=1,nx
         x(n)=zero
      enddo
c
      return
      end
