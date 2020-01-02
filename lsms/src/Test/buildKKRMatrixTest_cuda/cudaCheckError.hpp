#ifndef ME_CUDACHECKERROR_H
#define ME_CUDACHECKERROR_H

// inline int cudaCheckError(void) {return 0;}
#include <cstdio>

#define cudaCheckError() {                                          \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}

#define cublasCheckError(e) {                                          \
 if(e!=CUBLAS_STATUS_SUCCESS) {                                              \
   printf("CUBLAS failure %s:%d:\n",__FILE__,__LINE__);           \
   exit(0); \
 }                                                                 \
}

#endif
