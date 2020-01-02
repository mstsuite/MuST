c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mbeqa(a,b,n)
c     ================================================================
c
      integer    i,n
c
      real*8     a(n)
      real*8     b(n)
c
      do i=1,n
         b(i)=a(i)
      enddo
c
      return
      end
