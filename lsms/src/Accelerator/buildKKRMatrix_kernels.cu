// -*- mode: c++; -*-

#include <stdlib.h>
#include <Complex.hpp>
#include "cudaDoubleComplex.hpp"

#include <iostream>
#include <cublas_v2.h>
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() {return 1;}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
#endif
#include <DeviceMatrix.hpp>
#include <DeviceVector.hpp>
#include <DeviceArray3d.hpp>
using namespace std;

#define MAX_BLOCKS 512  //This gives 32 blocks an SM which should be plenty to keep batched the SM's busy
#define NUM_BLOCKS(N,NUM_THREADS) min( (N+NUM_THREADS-1)/NUM_THREADS, MAX_BLOCKS )

#include "cudaCheckError.hpp"

template<class T> 
__global__ void cudaMemset_kernel(T* mem, T val, int N) {
  for(int idx=blockIdx.x*blockDim.x+threadIdx.x; idx<N; idx+=blockDim.x*gridDim.x) {
    mem[idx]=val;
  }
}

template<class T> 
__inline__ void cudaMemset_custom(T* mem, const T val, int N, cudaStream_t s) {
  int num_threads=256;
  int num_blocks=NUM_BLOCKS(N,num_threads);
  cudaMemset_kernel<<<num_blocks,num_threads,0,s>>>(mem,val,N);
}

//initialize the diagonal of a matrix
__global__
void setDiagonal_kernel(cudaDoubleComplex *m, cudaDoubleComplex v, int rows, int lda) {
  for(int idx=blockIdx.x*blockDim.x+threadIdx.x;idx<rows;idx+=blockDim.x*gridDim.x) {
    m[idx*lda+idx]=v;
  }
}

//initialize the diagonal of a matrix
void setDiagonal_cuda(Complex *m, int rows, int lda,  cudaStream_t s ) {
 
  int num_threads=256;
  int num_blocks=NUM_BLOCKS(rows,num_threads);

  //write 1's to the diagonal
  cudaDoubleComplex v(1.,0.);
  
  setDiagonal_kernel<<<num_blocks,num_threads,0,s>>>((cudaDoubleComplex*)m,v,rows,lda);
  cudaCheckError();
}

__global__ void clearM00_kernel(cudaDoubleComplex *m, cudaDoubleComplex v, int blk_sz, int lda)
{
  int off=blockIdx.x*lda;
  for(int i=0; i<blk_sz; i++)
    m[off+i]=v;
}

// clear m00 block
void clearM00(Complex *m, int blk_sz, int lda, cudaStream_t s)
{
  cudaDoubleComplex v(0.0,0.0);
  clearM00_kernel<<<blk_sz,1,0,s>>>((cudaDoubleComplex *)m,v,blk_sz,lda);
  cudaCheckError();
}

//copy the matrix from tmat_store into tmat_n
__global__
void makeTmat_batched_kernel(int* atom_liz_store_idx, int* LIZlmax, cudaDoubleComplex* tmat_store, int lDimTmatStore, int n_spin_cant, int kkrsz, cudaDoubleComplex *tmat_n) {
  //x,y,z,w
  int ir1=blockIdx.y;
  int kkr1=(LIZlmax[ir1]+1)*(LIZlmax[ir1]+1);
  int xdim=kkr1;
  int ydim=xdim*n_spin_cant;
  int zdim=ydim*kkr1;
  int N=zdim*n_spin_cant;

  tmat_n+=4*kkrsz*kkrsz*ir1;

  for(int idx=blockIdx.x*blockDim.x+threadIdx.x;idx<N;idx+=gridDim.x*blockDim.x) {
    int js=idx/zdim;
    int j=(idx%zdim)/ydim;
    int is=(idx%ydim)/xdim;
    int i=idx%xdim;

    int jsm=kkrsz*kkrsz*n_spin_cant*js;
    int jm=jsm+kkrsz*n_spin_cant*j+kkrsz*is;

    tmat_n[idx]=tmat_store[ atom_liz_store_idx[ir1]*lDimTmatStore + jm + i];
  }
}

__global__
void makeTmat_batched_kernel(int* atom_liz_store_idx, int* LIZlmax, cudaDoubleComplex* tmat_store, int lDimTmatStore, int n_spin_cant, int kkr1, int kkrsz, cudaDoubleComplex *tmat_n) {
  int ld1=kkr1;
  int ld2=ld1*n_spin_cant;
  int ld3=ld2*kkr1;

  int ir1=blockIdx.y;
  int js=blockIdx.x;
  int j=threadIdx.z;
  int is=threadIdx.y;
  int k=threadIdx.x;

  int im=js*ld3+j*ld2+is*ld1;
   
  //jump to the current tmat
  tmat_n+=4*kkrsz*kkrsz*ir1;

  int jsm = kkrsz*kkrsz*n_spin_cant*js;  //linear through matrix of size kkrsz*kkrsz_ns
  int jm=jsm+kkrsz*n_spin_cant*j+kkrsz*is;
  
  tmat_n[im+k]=tmat_store[ atom_liz_store_idx[ir1]*lDimTmatStore + jm + k];
}

void makeTmat_cuda_batched(DeviceVector<int> &LIZStoreIdx, DeviceVector<int> &LIZlmax, int iie, int blkSizeTmatStore, DeviceMatrix<Complex> &tmat_store, int n_spin_cant, int numLIZ, int kkr1, int kkrsz, Complex *tmat_n, cudaStream_t s) {
  dim3 threads,blocks;
  threads.x=kkrsz;
  threads.y=n_spin_cant;
  threads.z=kkrsz;
  blocks.x=n_spin_cant;
  blocks.y=numLIZ;
  int blk_offset=iie*blkSizeTmatStore;
  
  makeTmat_batched_kernel<<<blocks,threads,0,s>>>(LIZStoreIdx.raw(), LIZlmax.raw(), (cudaDoubleComplex*)tmat_store.raw()+blk_offset, tmat_store.l_dim(), n_spin_cant, kkr1, kkrsz, (cudaDoubleComplex*)tmat_n);
  cudaCheckError();
}

//helper class specifying the shared memory size we need to for makeGBijs
//Note: the order of these appear to affect the compiler.  This order appears to work but other orders may crash! 
class SM_sizes {
  public:
  int hfn_off;
  int dlm_off;
  int sinmp_off;
  int cosmp_off;
  int plm_off;
  int total;
  SM_sizes (int lmax, int ndlm, int ndlj) {
    int hfn_size=(2*lmax+1)*sizeof(cudaDoubleComplex)*2;
    int dlm_size=ndlj*sizeof(cudaDoubleComplex);
    int sinmp_size=(2*lmax+1)*sizeof(double);
    int cosmp_size=(2*lmax+1)*sizeof(double);
    int plm_size=ndlm*sizeof(double);

    hfn_off=0;
    dlm_off=hfn_off+hfn_size;
    sinmp_off=dlm_off+dlm_size;
    cosmp_off=sinmp_off+sinmp_size;
    plm_off=cosmp_off+cosmp_size;
    total=plm_off+plm_size;
  }
};

#include "makebgij_device.hpp"

//constructs a batch of B matrices.  
__global__
void 
__launch_bounds__(128, 8) 
makeBGijs_kernel_batched(SM_sizes sm, int kkrsz, int maxlmax, int ndlj, int ndlm, int nspin, cudaDoubleComplex prel,  cudaDoubleComplex energy,
                     DeviceMatrix<double> LIZPos, DeviceVector<int> LIZlmax, DeviceVector<double> clm, DeviceArray3d<double> cgnt, 
                     int gntlmax, DeviceVector<int> lofk, DeviceVector<int> mofk, cudaDoubleComplex* ilp1, DeviceMatrix<cudaDoubleComplex> illp, cudaDoubleComplex *bgij) {

  extern char __shared__ sm_mem[];
  double* sinmp=(double*)(sm_mem+sm.sinmp_off);
  double* cosmp=(double*)(sm_mem+sm.cosmp_off);
  double* plm=(double*)(sm_mem+sm.plm_off);
  cudaDoubleComplex* hfn=(cudaDoubleComplex*)(sm_mem+sm.hfn_off);
  cudaDoubleComplex* dlm=(cudaDoubleComplex*)(sm_mem+sm.dlm_off);
  
  int ir1=blockIdx.y;
  int ir2=blockIdx.x;

  if(ir1!=ir2) {
    //advance bgij to this block
    bgij+=gridDim.x*4*kkrsz*kkrsz*ir1 + 4*kkrsz*kkrsz*ir2;

    double rij[3];

    rij[0]=LIZPos(0,ir1)-LIZPos(0,ir2);
    rij[1]=LIZPos(1,ir1)-LIZPos(1,ir2);
    rij[2]=LIZPos(2,ir1)-LIZPos(2,ir2);
    
    int kkr1=(LIZlmax[ir1]+1)*(LIZlmax[ir1]+1);
    int kkr2=(LIZlmax[ir2]+1)*(LIZlmax[ir2]+1);
    Real pi4=4.0*2.0*std::asin(1.0);

    makebgij_device(LIZlmax[ir1],kkr1,LIZlmax[ir2],kkr2,maxlmax,kkrsz,ndlj,ndlm, nspin,
        prel,rij,sinmp,cosmp,clm.raw(),plm,cgnt,gntlmax,lofk.raw(),mofk.raw(), ilp1,illp,hfn,dlm,bgij, pi4);
  }
}

void makeBGijs_cuda_batched(int numLIZ, int ndlj, int ndlm, int lmax, int gntlmax, int kkrsz, int nspin, Complex prel, Complex energy,
    DeviceVector<int> &LIZlmax, DeviceMatrix<Real> &LIZPos, DeviceVector<Real> &clm, DeviceArray3d<Real> &cgnt, 
    DeviceVector<int> &lofk, DeviceVector<int> &mofk, DeviceVector<Complex> &ilp1, DeviceMatrix<Complex> &illp, 
    Complex *bgij, cudaStream_t s) {
 
  cudaDoubleComplex pr=cudaDoubleComplex(prel.real(),prel.imag());
  cudaDoubleComplex e=cudaDoubleComplex(energy.real(),energy.imag());
  int num_threads=128; 
  dim3 num_blocks;
  //1 block per atom
  num_blocks.x=numLIZ; 
  num_blocks.y=numLIZ; 

  SM_sizes sm(lmax,ndlm,ndlj);
  DeviceMatrix<cudaDoubleComplex> *illp_hack=(DeviceMatrix<cudaDoubleComplex>*)&illp;
    
  const cuDoubleComplex czero=make_cuDoubleComplex(0.0,0.0);
  cudaMemset_custom((cuDoubleComplex*)bgij, czero , 4*kkrsz*kkrsz*numLIZ*numLIZ, s);
  cudaCheckError();

  //cudaFuncSetCacheConfig(makeBGijs_kernel,cudaFuncCachePreferL1);
  makeBGijs_kernel_batched<<<num_blocks,num_threads,sm.total,s>>>(sm,kkrsz,lmax,ndlj,ndlm,nspin,pr,e,LIZPos,LIZlmax,clm,cgnt,gntlmax,lofk,mofk,(cudaDoubleComplex*)ilp1.raw(),*illp_hack, (cudaDoubleComplex*)bgij);
  cudaCheckError();
}

//computes C = A * B, with N=kkr1_ns M=kkr2_ns, K=nrmat_ns
void zgemm_cuda_batched(int numLIZ, int kkrsz, int nspin, int nrmat_ns, cublasHandle_t &cublas_h, Complex* a, Complex *b, Complex *c, cudaStream_t s) {

  int kkrsz_ns=kkrsz*nspin;
  //save old stream
  cudaStream_t old;
  cublasGetStream(cublas_h,&old);
  //set cublas stream
  cublasSetStream(cublas_h,s);
  const cuDoubleComplex cmone=make_cuDoubleComplex(-1.0,0.0);
  const cuDoubleComplex czero=make_cuDoubleComplex(0.0,0.0);

  //create batches
  static vector<cuDoubleComplex*> h_av[16], h_bv[16], h_cv[16];

  vector<cuDoubleComplex*> &h_a=h_av[omp_get_thread_num()];
  vector<cuDoubleComplex*> &h_b=h_bv[omp_get_thread_num()];
  vector<cuDoubleComplex*> &h_c=h_cv[omp_get_thread_num()];

  static bool initialized[16]={0};

  if(!initialized[omp_get_thread_num()]) {
    h_a.resize(numLIZ);
    h_b.resize(numLIZ);
    h_c.resize(numLIZ);

    cudaHostRegister(&h_a[0],h_a.size()*sizeof(cuDoubleComplex*),0);
    cudaHostRegister(&h_b[0],h_b.size()*sizeof(cuDoubleComplex*),0);
    cudaHostRegister(&h_c[0],h_c.size()*sizeof(cuDoubleComplex*),0);
    initialized[omp_get_thread_num()]=true;
  }
  
  for(int i=0;i<numLIZ;i++) {
    h_a[i]=reinterpret_cast<cuDoubleComplex*>(a)+4*kkrsz*kkrsz*i;
    h_b[i]=reinterpret_cast<cuDoubleComplex*>(b)+4*kkrsz*kkrsz*numLIZ*i;
    h_c[i]=reinterpret_cast<cuDoubleComplex*>(c)+kkrsz*nspin*i;
  }

  static DeviceVector<cuDoubleComplex*> d_av[16], d_bv[16], d_cv[16];

  DeviceVector<cuDoubleComplex*> &d_a=d_av[omp_get_thread_num()];
  DeviceVector<cuDoubleComplex*> &d_b=d_bv[omp_get_thread_num()];
  DeviceVector<cuDoubleComplex*> &d_c=d_cv[omp_get_thread_num()];

  //copy batch vectors to the device
  d_a.copy_async(h_a,s);
  d_b.copy_async(h_b,s);
  d_c.copy_async(h_c,s);

  cublasCheckError(
      cublasZgemmBatched(cublas_h, 
        CUBLAS_OP_N, 
        CUBLAS_OP_N, 
        kkrsz_ns, 
        kkrsz_ns*numLIZ, 
        kkrsz_ns, 
        &cmone, 
        (const cuDoubleComplex**)d_a.raw(), 
        kkrsz_ns, 
        (const cuDoubleComplex**)d_b.raw(), 
        kkrsz_ns,
        &czero,
        d_c.raw(),
        nrmat_ns,
        numLIZ)
      );

  //  cudaStreamSynchronize(s);

  //restore old stream
  cublasSetStream(cublas_h,old);
}
