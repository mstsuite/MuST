#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_CUDA_C)
#include "Accelerator/DeviceStorage.hpp"
#endif

extern "C" {
  void zblock_lu_(Complex *a, int *lda, int *blk_sz, int *nblk, int *ipvt, int *mp, int *idcol, int *k);
  void  zblock_lu_cuda_c_(Complex *a, int *lda, int *blk_sz, int *nblk, int *ipvt, int *mp, int *idcol, int *k);
}

int zblock_lu_cpp(Matrix<Complex> &a, int *blk_sz, int nblk, int *ipvt, int *idcol);

#if defined(ACCELERATOR_CUBLAS)
int zblock_lu_cublas(cublasHandle_t handle, Matrix<Complex> &a, int *blk_sz, int nblk, int *ipvt, int *idcol);
#endif

void block_inverse(Matrix<Complex> &a, int *blk_sz, int nblk, Matrix<Complex> &delta, int *ipvt, int *idcol)
{
  int k;
  
  for(int i=0; i<blk_sz[0]; i++)
    for(int j=0; j<blk_sz[0]; j++)
      a(j,i)=0.0;

#if defined(ACCELERATOR_CUBLAS)
  zblock_lu_cublas(DeviceStorage::getCublasHandle(), a, blk_sz, nblk, ipvt, idcol);
#elif  defined(ACCELERATOR_CUDA_C)
  int mp=a.l_dim();
  int lda=a.l_dim();
  zblock_lu_cuda_c_( &a(0,0), &lda, blk_sz, &nblk, ipvt, &mp, idcol, &k);
#else
  k=zblock_lu_cpp(a, blk_sz, nblk, ipvt, idcol);
  // int mp=a.l_dim();
  // int lda=a.l_dim();
  // zblock_lu_(&a(0,0), &lda, blk_sz, &nblk, ipvt, &mp, idcol, &k);
#endif
  
  for(int i=0; i<blk_sz[0]; i++)
    for(int j=0; j<blk_sz[0]; j++)
      delta(j,i) = -a(j,i);
  
}

/*
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine block_inv(a,lda,na,mp,ipvt,blk_sz,nblk,delta,
     &                     idcol,iprint)
c     ================================================================
c
c     ****************************************************************
c     PURPOSE:   inverse the first block of a complex matrix: a
c                (a^{-1})_00=(a_00-delta)^{-1}
c
c     INPUT:     a,      the complex matrix to be inverted
c     OUTPUT:    a,   contains the first block of the inverted
c                            matrix
c                delta
c     ****************************************************************

      implicit none
      integer lda,na,mp,nblk
      integer ipvt(mp),blk_sz(nblk)
      integer i,j,k
      integer idcol(blk_sz(1))
      integer  iprint

      real*8 time,clock_time

      complex*16 a(lda,na)
      complex*16 delta(blk_sz(1),blk_sz(1))
      complex*16 czero
      parameter (czero=(0.d0,0.d0))
c

      time = clock_time()
      do i=1,blk_sz(1)
      do j=1,blk_sz(1)
        a(j,i)=czero
      enddo
      enddo


c        =============================================================
c        Use the  LU algorithm........................................
c        =============================================================
         call zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
c        -------------------------------------------------------------


        do i=1,blk_sz(1)
        do j=1,blk_sz(1)
          delta(j,i)=-a(j,i)
        enddo
        enddo

      time = clock_time() - time

      if( iprint .ge. 1 ) then
         write(6,'(''BLOCK_INV:: Using LU: Block size='',
     &     i4,'', Time='',f10.3,
     >   '' Sec'')')  blk_sz(2),time
         call flush(6)
      endif

      return
      end
*/
