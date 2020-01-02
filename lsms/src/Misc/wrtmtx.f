c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wrtmtx(x,n,istop)
c     ================================================================
c
      implicit   none
c
      character  sname*32
      character  istop*32
c
      integer    n
      integer    i
      integer    j
c
      real*8     tol
c
      complex*16 x(n,n)
c
      parameter  (sname='wrtmtx')
      parameter  (tol=1.0d-8)
c
c     ****************************************************************
c     writes out the non-zero elements of a NxN complex matrix........
c     ****************************************************************
c
      do j=1,n
         do i=1,n
            if(abs(x(i,j)).gt.tol) then
               write(6,'(2i4,1p3d23.13)') i,j,x(i,j),abs(x(i,j))
            endif
         enddo
      enddo
c
      if (istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
