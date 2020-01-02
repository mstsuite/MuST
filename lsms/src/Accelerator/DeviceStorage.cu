// -*- mode: c++; -*-

#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "DeviceMatrix.hpp"
#include "DeviceArray3d.hpp"
#include "DeviceVector.hpp"
#include "Main/SystemParameters.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <iostream>

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

// #include "cuda_error.h"
#include "cudaCheckError.hpp"

using namespace std;

//TODO move inside DeviceStorage?
//allocate a thread specific matrix on the host and pin its memory
extern "C"
Complex* get_host_m_(const int &max_nrmat_ns) {
  static Complex * m_v=0;
  static int cur_size=0;
  static cudaError_t pinned;

  if(cur_size<max_nrmat_ns) {

    //release previously allocated memory
    if(m_v!=0) {
      if(pinned) cudaFreeHost(m_v);
      else free(m_v);
    }

    //allocate new memory
    pinned = cudaMallocHost((void**)&m_v,max_nrmat_ns*max_nrmat_ns*sizeof(Complex)*omp_get_max_threads());

    if ( pinned != cudaSuccess )
    {
      fprintf(stderr, "Matrix not pinned\n");
      m_v = (Complex*)malloc(max_nrmat_ns*max_nrmat_ns*sizeof(Complex)*omp_get_max_threads());
    }
    cur_size=max_nrmat_ns;
  }
  return m_v; 
}

static const int MAX_THREADS=16;
class DeviceStorage {
private:
  static int nThreads;
  static Complex *dev_m[MAX_THREADS], *dev_bgij[MAX_THREADS], *dev_tmat_n[MAX_THREADS];
  static int *dev_ipvt[MAX_THREADS];
  static cublasHandle_t cublas_h[MAX_THREADS];
  static cudaEvent_t event[MAX_THREADS];
  static cudaStream_t stream[MAX_THREADS][2];
  static DeviceMatrix<Complex> dev_tmat_store;
  static bool initialized;
public:
  static int allocate(int kkrsz_max,int nspin, int numLIZ, int _nThreads)
  {
    if(!initialized)
    {
      //printf("*************************************MEMORY IS BEING ALLOCATED\n");
      if(_nThreads>MAX_THREADS)
      {
        printf("nThreads (%d) in DeviceStorage::allocate exceeds MAX_THREADS (%d)\n",_nThreads,MAX_THREADS);
        printf("  change MAX_THREADS in src/Accelerator/DeviceStorage.cu and recompile!\n");
        exit(1);
      }
      nThreads=_nThreads;
      int N=kkrsz_max*nspin*numLIZ;
      for(int i=0;i<nThreads;i++)
      {
        cudaMalloc((void**)&dev_m[i],N*N*sizeof(Complex));
        cudaMalloc((void**)&dev_ipvt[i],N*sizeof(int));
#ifdef BUILDKKRMATRIX_GPU
        cudaMalloc((void**)&dev_bgij[i],4*kkrsz_max*kkrsz_max*numLIZ*numLIZ*sizeof(Complex));
        cudaMalloc((void**)&dev_tmat_n[i],4*kkrsz_max*kkrsz_max*numLIZ*sizeof(Complex)); 
#endif
        cudaStreamCreate(&stream[i][0]);
        cudaStreamCreate(&stream[i][1]);
        cudaEventCreateWithFlags(&event[i],cudaEventDisableTiming);
        cublasCreate(&cublas_h[i]);
      }
      cudaCheckError();
      initialized=true;
    }
    return 0;
  }
  
  static void free()
  {
    if(initialized) {
   //     printf("*************************************MEMORY IS BEING FREED\n");
      // for(int i=0;i<omp_get_max_threads();i++)
      for(int i=0; i<nThreads; i++)
      {
        cudaFree(dev_m[i]);
        cudaFree(dev_ipvt[i]);
#ifdef BUILDKKRMATRIX_GPU
        cudaFree(dev_bgij[i]);
        cudaFree(dev_tmat_n[i]);
#endif
        cudaStreamDestroy(stream[i][0]);
        cudaStreamDestroy(stream[i][1]);
        cudaEventDestroy(event[i]);
        cublasDestroy(cublas_h[i]);
      }
      dev_tmat_store.clear();
      cudaCheckError();
      initialized=false;
    }
  }

  static Complex* getDevM() { return dev_m[omp_get_thread_num()]; } 
  static Complex* getDevBGij() { if(!initialized) {printf("DeviceStorage not initialized\n"); exit(1);}
                                 return dev_bgij[omp_get_thread_num()]; } 
  static Complex* getDevTmatN() { return dev_tmat_n[omp_get_thread_num()]; } 
  static int* getDevIpvt() { return dev_ipvt[omp_get_thread_num()]; } 
  static cudaStream_t getStream(int i) { return stream[omp_get_thread_num()][i]; }
  static cudaEvent_t getEvent() { return event[omp_get_thread_num()]; }
  static cublasHandle_t getCublasHandle() { return cublas_h[omp_get_thread_num()]; }
  static DeviceMatrix<Complex>* getDevTmatStore() { return &dev_tmat_store; }
};

bool DeviceStorage::initialized = false;
Complex *DeviceStorage::dev_m[MAX_THREADS], *DeviceStorage::dev_bgij[MAX_THREADS], *DeviceStorage::dev_tmat_n[MAX_THREADS];
int *DeviceStorage::dev_ipvt[MAX_THREADS];
cublasHandle_t DeviceStorage::cublas_h[MAX_THREADS];
cudaEvent_t DeviceStorage::event[MAX_THREADS];
cudaStream_t DeviceStorage::stream[MAX_THREADS][2];
DeviceMatrix<Complex> DeviceStorage::dev_tmat_store;
int DeviceStorage::nThreads=1;
bool initialized;


/****************Fortran Interfaces*********************/
extern "C"
Complex* get_dev_m_() {
  return DeviceStorage::getDevM();
}

extern "C"
Complex* get_dev_bgij_() {
  return DeviceStorage::getDevBGij();
}

extern "C"
Complex* get_dev_tmat_n_() {
  return DeviceStorage::getDevTmatN();
}

extern "C"
int* get_dev_ipvt_() {
  return DeviceStorage::getDevIpvt();
}

extern "C"
cudaStream_t get_stream_(const int &id) {
  return DeviceStorage::getStream(id);
}

extern "C"
cublasHandle_t get_cublas_handle_() {
  return DeviceStorage::getCublasHandle();
}

//allocate a thread specific event
extern "C"
cudaEvent_t get_cuda_event_() {
  return DeviceStorage::getEvent();
}
/********************************************************/

DeviceMatrix<Complex>* get_dev_tmat_store() {
  return DeviceStorage::getDevTmatStore();
}

void *allocateDStore(void)
{
  return static_cast<void *>(new DeviceStorage);
}

void freeDStore(void * d_store)
{
  static_cast<DeviceStorage*>(d_store)->free();
  delete static_cast<DeviceStorage*>(d_store);
}

int initDStore(void * d_store,int kkrsz_max, int nspin, int numLIZ, int nthreads)
{
  return (*static_cast<DeviceStorage*>(d_store)).allocate(kkrsz_max,nspin,numLIZ,nthreads);
}

void copyTmatStoreToDevice(LocalTypeInfo &local) {
  DeviceMatrix<Complex> &d_tmat_store=*get_dev_tmat_store();
  d_tmat_store.copy_async(local.tmatStore,0);
}
