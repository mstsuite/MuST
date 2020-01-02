#include "Real.hpp"
#include "Complex.hpp"
#include <cublas_v2.h>

template <class T> class DeviceMatrix;

extern "C" Complex* get_dev_m_();
extern "C" Complex* get_dev_bgij_();
extern "C" Complex* get_dev_tmat_n_();
extern "C" int* get_dev_ipvt_();
extern "C" cudaStream_t get_stream_(const int &id);
extern "C" cublasHandle_t get_cublas_handle_();
extern "C" cudaEvent_t get_cuda_event_();
extern "C" Complex* get_host_m_(const int &max_nrmat_ns);



DeviceMatrix<Complex>* get_dev_tmat_store();

void *allocateDStore(void);
void freeDStore(void * d_store);
int initDStore(void * d_store,int kkrsz_max, int nspin, int numLIZ, int nthreads);
