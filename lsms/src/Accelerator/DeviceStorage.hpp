#ifndef LSMS_DEVICE_STORAGE_HPP
#define LSMS_DEVICE_STORAGE_HPP

#include "Real.hpp"
#include "Complex.hpp"
#include <cublas_v2.h>

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

// #include "DeviceMatrix.hpp"

template <class T> class DeviceMatrix;

extern "C" Complex* get_dev_m_();
extern "C" Complex* get_dev_bgij_();
extern "C" Complex* get_dev_tmat_n_();
extern "C" int* get_dev_ipvt_();
extern "C" cudaStream_t get_stream_(const int &id);
extern "C" cublasHandle_t get_cublas_handle_();
extern "C" cudaEvent_t get_cuda_event_();
extern "C" Complex* get_host_m_(const int &max_nrmat_ns);

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
  static int allocate(int kkrsz_max,int nspin, int numLIZ, int _nThreads);
  static void free();

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

DeviceMatrix<Complex>* get_dev_tmat_store();

void *allocateDStore(void);
void freeDStore(void * d_store);
int initDStore(void * d_store,int kkrsz_max, int nspin, int numLIZ, int nthreads);

#endif
