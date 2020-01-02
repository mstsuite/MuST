#include "DeviceStorage.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include "cuda_runtime.h"
#include "cublas_v2.h"

#ifdef _OPENMP
#include <omp.h>
#else
#ifndef LSMS_DUMMY_OPENMP
#define LSMS_DUMMY_OPENMP
inline int omp_get_max_threads() {return 1;}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
#endif
#endif

#include "cudaCheckError.hpp"
#include "cudaDoubleComplex.hpp"
#include <assert.h>
#include <lapack.h>
//#include <mpi.h>

extern "C" int zmatinv_prep1_ (void **a, void **b, int *n, int *lda, cudaStream_t thisstream);
extern "C" int zmatinv_batch_ (cuDoubleComplex **A, cuDoubleComplex **Ainv, int *n, int *batch, cudaStream_t thisstream);
extern "C" int ilaenv_(int*,const char*,const char*,int*,int*,int*,int*);

void handle_cuda_error ( cudaError_t cerr, const char *errmsg ) 
{
  if ( cerr ) {
    printf ("CUDA ERROR (%d) %s \n", cerr, errmsg);
    abort();
    //MPI_Abort(MPI_COMM_WORLD, 100);
  }
  else {
    //printf ("SUCCESS !!! %s \n", errmsg);
  }
}

void handle_cublas_error ( cublasStatus_t cs, const char *errmsg ) 
{
  if ( cs ) {
    printf ("cuBLAS ERROR (%d) %s \n", cs, errmsg);
    abort();
    //MPI_Abort(MPI_COMM_WORLD, 101);
  }
  else {
    //printf ("SUCCESS !!! %s \n", errmsg);
  }
}


//TODO call directly from calculateTauMatrix (don't route through fortran)
extern "C"
void zblock_lu_cuda_c_ ( std::complex<double> *a, int *lda, int *blk_sz, int *nblk, int *ipvt, int *mp, int *idcol, int *k)
  //===================================================================================================================
/*
  Performs a partial inversion of the a matrix to return the inverse of the upper diagonal
  subblock.  

  a      : input matrix - double complex
  blk_sz : integer array giving the size of each subblock
  nblk   : the number of subblocks
  ipvt   : integer work array (not tested in c version)
  idcol  : integer array specifying symmetry (not tested in c version)
  k      : returns the actual number of columns in the calculated inverse

*/
{
  //TODO:  
      //  adjust allocation sizes
      //  dynamically choose hybrid or not
      //  validate flop count

  unsigned long long flops=0;
  /********************paramters for zgemm rank maximization*******************/
  int zgemm_rank = 600;
  int gpu_only_blks = 0;

  int remaining=0;
  for(int i=0;i<gpu_only_blks;i++) {
    (*nblk)--;
    remaining+=blk_sz[*nblk];
  }

  while(remaining>0) {
    blk_sz[*nblk]=min(55,remaining);
    remaining-=55;
    (*nblk)++;
  }
  
  int currentRank=0;

  int m, n;
  int ioff, joff;
  int info;
  cudaError_t ce;
  cublasStatus_t cublasStat;

  // set constants
  const cuDoubleComplex cone  = make_cuDoubleComplex( 1.0, 0.0);
  const cuDoubleComplex cmone = make_cuDoubleComplex(-1.0, 0.0);
  const cuDoubleComplex czero = make_cuDoubleComplex( 0.0, 0.0);
  
  // get the thread number
  int threadId = omp_get_thread_num();

  /***************************One time initialization, should be moved outside******************************************/
  int max_blk_sz=blk_sz[0];
  if(*nblk>1)
  {
    for(int i=1; i<*nblk; i++)
      max_blk_sz=max(max_blk_sz,blk_sz[i]);
  }

  const int MAX_THREADS=16;
  //TODO dynamically size
  static bool initialized[MAX_THREADS] = {false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false};
  static cuDoubleComplex *vdevWork[MAX_THREADS];
  static cuDoubleComplex *vdevInv[MAX_THREADS];
  static cuDoubleComplex *vdevA2[MAX_THREADS];
  static std::complex<double> *vdevHostDiag[MAX_THREADS];
  static std::complex<double> *vwork[MAX_THREADS];
  static int *vhostIPVT[MAX_THREADS];
  static int lwork;
  if ( ! initialized[threadId] ) {

    //calculate optimial work size for zgetri
    int one=1;  int mone=-1;
    int NB = ilaenv_(&one,"ZGETRI","",&max_blk_sz,&mone,&mone,&mone);
    lwork= max_blk_sz * NB; 

    // allocate space on device
    ce = cudaMalloc ( &vdevWork[threadId], max_blk_sz*max_blk_sz*sizeof(cuDoubleComplex));
    handle_cuda_error (ce, "cudaMalloc devWork");
    ce = cudaMalloc ( &vdevInv[threadId], max_blk_sz*max_blk_sz*sizeof(cuDoubleComplex));
    handle_cuda_error (ce, "cudaMalloc devInv");
  

    int LDA= *lda;
    ce = cudaMallocHost ( &vwork[threadId], lwork *sizeof(std::complex<double>));
    handle_cuda_error (ce, "cudaMallocHost vwork");
    ce = cudaMalloc ( &vdevA2[threadId], max_blk_sz * LDA *sizeof(cuDoubleComplex));
    handle_cuda_error (ce, "cudaMalloc devA2");
    ce = cudaMallocHost ( &vdevHostDiag[threadId], max_blk_sz * max_blk_sz *sizeof(std::complex<double>));
    handle_cuda_error (ce, "cudaMallocHost vdevHostDiag");
    ce = cudaMallocHost((void**)&vhostIPVT[threadId], max_blk_sz*sizeof(int));
    handle_cuda_error (ce, "cudaMallocHost vhostIPVT");
    
    //this speeds up the small block inverse
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    initialized[threadId] = true;
  }
  /**********************************************************************************************************************/

  /********************assign thread private variables********************************/
  cudaStream_t stream1=get_stream_(0);
  cudaStream_t stream2=get_stream_(1);
  cudaEvent_t done_event=get_cuda_event_();
  cuDoubleComplex *devWork = vdevWork[threadId];
  cuDoubleComplex *devInv = vdevInv[threadId];
  cuDoubleComplex *devA=(cuDoubleComplex*)get_dev_m_();
  cuDoubleComplex *devA2 = vdevA2[threadId];
  std::complex<double> *work = vwork[threadId];
  cublasHandle_t cublasHandle = get_cublas_handle_();
  int *hostIPVT=vhostIPVT[threadId];
  Complex *hostAdiag = (Complex*)vdevHostDiag[threadId];
  /***********************************************************************************/

  // add up the sizes of the subblocks to get the size of the entire matrix
  int na;
  na = 0;
  for ( int i=0; i<abs(*nblk); i++ ) {
    na += blk_sz[i];
  }

  // eliminate columns that are equivalent due to symmetry
  if ( idcol[0] == 0 ) {
    *k = 1;
  }
  else {
    *k = blk_sz[0]+1;
    for ( int i=blk_sz[0]-1; i>=0; i-- ) {
      if ( idcol[0] == 0 || idcol[i] == i ) {
        *k -= 1;
        if ( *k != i ) {
          printf ("Eliminate columns that are equivalent due to symmetry section in zblock_lu_cuda_c not tested\n");
          abort();
          // zcopy ( na-blk_sz[0], a[i*lda+blk_sz[0]], 1, a[*k*lda+blk_sz[0]], 1 );
        }
      }
    }
  }

#ifndef BUILDKKRMATRIX_GPU
  // copy matrix to device 
  cublasStat = cublasSetMatrix ( na, na, sizeof(cuDoubleComplex), a, *lda, devA, *lda); 
  handle_cublas_error ( cublasStat, "cublasSetMatrix devA ");
#endif
  if ( *nblk > 0 ) {

    n = blk_sz[*nblk-1];
    joff = na - n;

    // loop over sub-blocks
    for ( int iblk=*nblk-1; iblk>0; iblk-- ) {  
      m = n;
      ioff = joff;
      n = blk_sz[iblk-1];
      joff = joff-n;

      //TODO update condition to chose branch, should do this branch when remaining size is small...

      // HPL factorization and left propagation
      if ( m<56 ) { //CUDA only version
        //A^-1 // invert the clique
        // re-package the diagonal block into a dense matrix suitable for sending to zmatinv
        cuDoubleComplex *devAdiag;
        devAdiag = &devA[ioff* *lda + ioff];
        info = zmatinv_prep1_ ( (void**)&devAdiag, (void**)&devWork, &m, lda, stream1 );
        if ( info ) { printf (" zmatinv_prep1 returned error code %d \n", info); abort();}

        int one = 1;
        info = zmatinv_batch_ ( &devWork, &devInv, &m, &one, stream1 );
        if ( info ) { printf (" zmatinv_batch returned error code %d \n", info); printf (" m = %d, one = %d \n", m, one ); abort(); }
        flops += m * m * m;

      }
      else { //HYBRID version, do small inverse on the host. This works well.

        cuDoubleComplex *devAdiag = (cuDoubleComplex*)&devA[ioff* *lda + ioff];

        cublasSetStream ( cublasHandle, stream1 );
        
        cublasStat = cublasGetMatrixAsync ( m, m, sizeof(cuDoubleComplex), devAdiag, *lda, hostAdiag, m, stream1 );

        cudaEventRecord(done_event,stream1);
        //wait for transfers to the host to finish
        cudaEventSynchronize(done_event);
        
        int info;
        //zgetrf on host
        zgetrf_(&m, &m, hostAdiag, &m, hostIPVT, &info);

        //zgetri on host
        zgetri_(&m, hostAdiag, &m, hostIPVT, (Complex*)work, &lwork, &info);
        
        flops += m * m * m;

        //copy_async down to device
        cublasStat = cublasSetMatrixAsync ( m, m, sizeof(cuDoubleComplex), hostAdiag, m, devInv, m, stream1 );

        cudaEventRecord(done_event,stream1);
        //wait for transfers to the host to finish
        cudaEventSynchronize(done_event);
      }

      //CA^-1
#ifdef PRINT_ZGEMM
      fprintf(stderr, "m: %d n: %d k: %d lda: %d ldb: %d ldc: %d\n", m, ioff, m, m, *lda, max_blk_sz);
#endif
      cublasStat = cublasZgemm ( cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, m, ioff, m, &cone, devInv, m, &devA[ioff], *lda, &czero, devA2, max_blk_sz );
      handle_cublas_error ( cublasStat, "Error in cublasZgemm #1\n" );
      flops+= m * ioff * m;

      //Mark end of small zgemm in stream1
      cudaEventRecord(done_event,stream1);  

      //stream 2 must wait for the small zgemm to finish
      cudaStreamWaitEvent(stream2,done_event,0);

      // Trailing matrix update 
      currentRank+=m;

      if ( currentRank<zgemm_rank && iblk>1)  {

        // only update the next block row
        // little chance for hybrid acceleration here - so ignore for now.

        cublasSetStream ( cublasHandle, stream1 );

        // need to place A2 back into A
        ce = cudaMemcpy2DAsync ( &devA[ioff], *lda*sizeof(cuDoubleComplex), devA2, max_blk_sz*sizeof(cuDoubleComplex), m*sizeof(cuDoubleComplex), ioff, cudaMemcpyDeviceToDevice, stream1 );

#ifdef PRINT_ZGEMM
        fprintf(stderr, "m: %d n: %d k: %d lda: %d ldb: %d ldc: %d\n", n, ioff, currentRank, *lda, *lda, *lda);
#endif
        cublasStat = cublasZgemm ( cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, ioff, currentRank, &cmone,
            &devA[ioff* *lda +ioff-n], *lda,
            //devA2, max_blk_sz, &cone,
            &devA[ioff], *lda, &cone,
            &devA[ioff-n], *lda ); 
        
        flops += n * ioff * currentRank;

#ifdef PRINT_ZGEMM
        fprintf(stderr, "m: %d n: %d k: %d lda: %d ldb: %d ldc: %d\n", ioff-n, n, currentRank, *lda, *lda, *lda);
#endif
        cublasStat = cublasZgemm ( cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, ioff-n, n, currentRank, &cmone,
            &devA[ioff * *lda], *lda,
            &devA[(ioff-n) * *lda + ioff], *lda, &cone,
            &devA[(ioff-n) * *lda], *lda ); 
        
        flops += (ioff-n)*n*currentRank;
      }
      else {
        // update the full trailing matrix 

        cublasSetStream ( cublasHandle, stream1 );
        // perform a portion of the zgemm on the gpu

        // first need to place A2 back into A
        ce = cudaMemcpy2DAsync ( &devA[ioff], *lda*sizeof(cuDoubleComplex), devA2, max_blk_sz*sizeof(cuDoubleComplex), m*sizeof(cuDoubleComplex), ioff, cudaMemcpyDeviceToDevice, stream1 );

        //D=CA^-1B
#ifdef PRINT_ZGEMM
        fprintf(stderr, "m: %d n: %d k: %d lda: %d ldb: %d ldc: %d\n", ioff, ioff, currentRank, *lda, *lda, *lda);
#endif
        cublasStat = cublasZgemm ( cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, ioff, ioff, currentRank, &cmone, 
            &devA[ioff* *lda], *lda, 
            &devA[ioff], *lda , &cone, 
            devA, *lda);

        flops += ioff * ioff * currentRank;

        // just did a full trailing submatrix update, so reset block row delay counter
        currentRank=0;
      }

    } // end for

    cublasStat = cublasGetMatrixAsync ( blk_sz[0], blk_sz[0], sizeof(cuDoubleComplex), devA, *lda, a, *lda, stream1 );
  } // end if ( *nblk > 0 )

  *k = blk_sz[0];

  cudaEventRecord(done_event,stream1);
  //wait for last transfer to finish
  cudaEventSynchronize(done_event);

  // clean up
  //cudaFree (devWork);
  //cudaFree (devInv);
  //cudaFree (devA);
  //cudaFree (devA2);
  //cudaFreeHost (ap);
#ifdef PRINT_FLOPS
  printf("BLOCK_INV ZGEMM FLOPS: %llu\n", flops*4*2); 
#endif
}
