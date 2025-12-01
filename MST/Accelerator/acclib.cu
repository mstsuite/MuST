#include <unistd.h>
#include <iostream>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "acclib.hpp"

void check_cuda_error(cudaError_t result, char const *const func, const char *const file, int const line) {
   if (result != cudaSuccess) {
      fprintf(stderr, "\nCUDA Runtime Error at %s:%d code=%d(%s) \"%s\"\n", file, line,
                      static_cast<unsigned int>(result), cudaGetErrorName(result), func);
      exit(EXIT_FAILURE);
   }
}

void check_cublas_error(cublasStatus_t status, char const *const func, const char *const file, int const line) {
    if (status != CUBLAS_STATUS_SUCCESS) {
       fprintf(stderr, "\ncuBLAS Runtime Error at %s:%d code=%d(%s) \"%s\"\n", file, line,
                       static_cast<unsigned int>(status), cublasGetStatusString(status), func);
       exit(EXIT_FAILURE);
    }
}

extern "C"
void get_node_resources_(int *num_cpu_cores, int *num_gpu_cards){
    *num_cpu_cores = sysconf(_SC_NPROCESSORS_ONLN);    // available cores
    // *num_cpu_cores = sysconf(_SC_NPROCESSORS_CONF); // configured cores

    int deviceCount = 0;
    checkCudaErrors(cudaGetDeviceCount(&deviceCount));
    *num_gpu_cards = deviceCount;
}
