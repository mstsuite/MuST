#include "Accelerator.hpp"
#if defined(ACCELERATOR_CUDA_C)
#include <iostream>
#include <cuda_runtime.h>
#endif

void acceleratorInitialize(int sz, int nthreads)
{
  accelerator_initialize_(&sz);
}

void acceleratorFinalize(void)
{
  accelerator_finalize_();
#if defined(ACCELERATOR_CUDA_C)
  cudaDeviceReset();
#endif
}

void acceleratorPrint(void)
{
#if defined(ACCELERATOR_CUDA_C)
  cudaError_t cuda_error;
  int cuda_count;
  cudaDeviceProp cuda_prop;
  cuda_error = cudaGetDeviceCount(&cuda_count);
  std::cout << "Found " << cuda_count << " CUDA GPUs." << std::endl;
  for (int i = 0; i < cuda_count; i++)
  {
    cuda_error = cudaGetDeviceProperties(&cuda_prop, i);
    std::cout << "Device " << i << ": " << cuda_prop.name << std::endl;
  }
#else
  ;
#endif
}
