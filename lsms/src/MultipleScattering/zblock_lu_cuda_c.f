! Notes
! 

      subroutine zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)

      implicit none
      integer lda,mp,nblk,k
      integer ipvt(mp),blk_sz(*)
      integer idcol(blk_sz(1))
      complex*16 a(lda,*)

      external zblock_lu_cuda_c

!      write (*,*), "Inside zblock_lu_cuda_c.f !!! "

      call zblock_lu_cuda_c
     &     ( a, lda, blk_sz, nblk, ipvt, mp, idcol, k)

      return
      
      end

