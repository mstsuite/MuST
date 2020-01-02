#if defined(ACCELERATOR_LIBSCI)
!#warning "Using LIBSCI"
#include "zblock_lu_cray.f"
#elif defined(ACCELERATOR_CUDA_C)
#include "zblock_lu_cuda_c.f"
#else
!#warning "Using CPU"
#include "zblock_lu_CPU.f"
c      subroutine zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
c      implicit none
c      write(*,*) "unknown accelerator in zblock_lu: STOP!"
c      stop
c      end subroutine
#endif
c$$$      subroutine zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
c$$$c does a partitioning of a matrix and a partial inversion to
c$$$c get the upper subblock of the inverse
c$$$c   a -- complex*16 of (lda,size) the matrix to be inverted
c$$$c        where size is the sum of the all elements of blk_sz
c$$$c        on return, the upper subblock of a is untouched, and
c$$$c        postprocessing is needed to obtain the inverse
c$$$c   blk_sz -- integer of (nblk), the size of each subblock of matrix a
c$$$c   nblk -- number of blocks (number of elements in blk_sz)
c$$$c   ipvt -- integer of (mp), work array
c$$$c   idcol -- integer of (blk_sz(1)), if idcol(i)=idcol(j), then
c$$$c            the two columns are equivalent by symmetry, and only
c$$$c            the min(i,j) th column of the inverse is calculated.
c$$$c   k -- returns actual number of columns in the calculated inverse

