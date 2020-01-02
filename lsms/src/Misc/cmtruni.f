c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmtruni(a,n)
c     ================================================================
c
      implicit   none
c
      integer    n
      integer    i
c
      complex*16 a(n,n)
      complex*16 cone
c
      parameter  (cone=(1.0d0,0.0d0))
c
c     ****************************************************************
c     set up a unit matrix and store it in a.
c     ****************************************************************
c
c     ----------------------------------------------------------------
      call zeroout(a,2*n*n)
c     ----------------------------------------------------------------
c
      do i=1,n
         a(i,i)=cone
      enddo
c
      return
      end
