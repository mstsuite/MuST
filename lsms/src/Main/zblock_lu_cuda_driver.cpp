#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SystemParameters.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "Accelerator/Accelerator.hpp"

//#ifdef BUILDKKRMATRIX_GPU
#include "Accelerator/DeviceStorage.hpp"
void *allocateDConst(void);
void freeDConst(void *);

std::vector<void *> deviceConstants;
// std::vector<void *> deviceStorage;
void * deviceStorage;
//#endif
extern "C"
void zblock_lu_cuda_c_ ( std::complex<double>*, int *, int *, int *, int *, int *, int *, int *);

int main(int argc, char *argv[])
{
  LSMSSystemParameters lsms;
  LocalTypeInfo local;

  lsms.global.setIstop("main");
  lsms.global.iprint=0;
  lsms.global.default_iprint=-1;
#ifdef _OPENMP
  lsms.global.GPUThreads=std::min(12,omp_get_max_threads());
#else
  lsms.global.GPUThreads=1;
#endif
  int nLIZ   = 113;
  int kkrsz_max = 16;
  int nspin = 2;
  int kkrsz_ns = kkrsz_max * nspin;
  int nrmat_ns = nLIZ * kkrsz_ns;

// initialize the potential accelerators (GPU)
// we need to know the max. size of the kkr matrix to invert: lsms.n_spin_cant*local.maxNrmat()
// which is only available after building the LIZ

  acceleratorInitialize(nrmat_ns,lsms.global.GPUThreads);
  local.tmatStore.pinMemory();

#ifdef BUILDKKRMATRIX_GPU
  deviceConstants.resize(local.num_local);
  for(int i=0; i<local.num_local; i++) deviceConstants[i]=allocateDConst();
  // deviceStorage.resize(lsms.global.GPUThreads);
  // for(int i=0; i<lsms.global.GPUThreads; i++) deviceStorage[i]=allocateDStore();
  deviceStorage=allocateDStore();
#endif
  initDStore(deviceStorage,kkrsz_max,nspin,nLIZ,lsms.global.GPUThreads);

  //FIXME cublas not getting initialized to DeviceStorage.allocate not getting
  //called
  Matrix<Complex> m(nrmat_ns,nrmat_ns);
  // nrmat_ns is sum of all block sizes (atoms in liz * 32)
  // Synthetic matrix from calculateTauMatrix.cpp
  for(int i=0; i<nrmat_ns; i++)
    for(int j=0; j<nrmat_ns; j++)
    {
      m(i,j)=Complex((11.0*(i+1)-3.0*(j+1))/Real(i+j+2),(5.0*(i+1)-2.0*(j+1))/Real(i+j+2));
    }

  // setup blocks (from calculateTauMatrices)
    int max_blk_sz=175;
    int nblk;
    //assign blocks in a load balanced way
    if((nrmat_ns-kkrsz_ns)%max_blk_sz==0)
      nblk=(nrmat_ns-kkrsz_ns)/max_blk_sz+1;
    else
    {
      nblk=(nrmat_ns-kkrsz_ns)/(max_blk_sz-1)+1;
      if((nrmat_ns-kkrsz_ns)%(max_blk_sz-1) > max_blk_sz/2)
        nblk++;
    }

    int blk_sz[1000];

    blk_sz[0]=kkrsz_ns;
    if(nblk==1)
      blk_sz[0]=nrmat_ns;
    else if(nblk==2)
      blk_sz[1]=nrmat_ns-blk_sz[0];
    else if(nblk>2)
    {
//  {
//    int min_sz=(nrmat_ns-blk_sz[0])/(nblk-1);
//    for(int i=1; i<nblk; i++) blk_sz[i]=min_sz;
//    blk_sz[nblk-1]=nrmat_ns-blk_sz[0]-(nblk-2)*min_sz;
//  }
      int min_sz=(nrmat_ns-blk_sz[0])/(nblk-1);
      int rem=(nrmat_ns-blk_sz[0])%(nblk-1);
      int i=1;
      for(i;i<=rem;i++)
        blk_sz[i]=min_sz+1;
      for(i;i<nblk;i++)
        blk_sz[i]=min_sz;
    }

    Complex *vecs;
    Complex tmp_vecs; vecs=&tmp_vecs;
    int *ipvt = new int[nrmat_ns];
    Complex *delta = new Complex[kkrsz_ns*kkrsz_ns];
    int *iwork = new int[nrmat_ns];
    //int *iwork = new int[kkrsz_ns*atom.numLIZ];
    Real *rwork = new Real[nrmat_ns];
    //Real *rwork = new Real[kkrsz_ns*atom.numLIZ];
    Complex *work1 = new Complex[nrmat_ns];
    //Complex *work1 = new Complex[kkrsz_ns*atom.numLIZ];
    int alg=2;
    if(alg>2) vecs=(Complex *)malloc((nrmat_ns*(kkrsz_ns*6+6))*sizeof(Complex));
    int *idcol = new int[blk_sz[0]]; idcol[0]=0;
#ifdef BUILDKKRMATRIX_GPU
    Complex *dev_m=get_dev_m_();
    clearM00(dev_m, blk_sz[0], nrmat_ns,get_stream_(0));
#endif
#ifdef _OPENMP
    double st = omp_get_wtime();
#endif
    {
      int k;
//      block_inv_(&m(0,0),vecs,&nrmat_ns,&nrmat_ns,&nrmat_ns,ipvt,
//          blk_sz,&nblk,delta,
//          iwork,rwork,work1,&alg,idcol,&lsms.global.iprint);
       zblock_lu_cuda_c_(&m(0,0),&nrmat_ns, blk_sz, &nblk, ipvt, &nrmat_ns,idcol,&k);
    }
#ifdef _OPENMP
    double et = omp_get_wtime();
    printf("Time: %lf (s)\n", et - st);
#endif

  local.tmatStore.unpinMemory();
#ifdef BUILDKKRMATRIX_GPU
  for(int i=0; i<local.num_local; i++) freeDConst(deviceConstants[i]);
  freeDStore(deviceStorage);
#endif

  acceleratorFinalize();
  return 0;
}
