#ifndef ACCLIB_HPP
#define ACCLIB_HPP

// CUDA error-checking macro for reliable memory allocation and transfers
#define checkCudaErrors(val) check_cuda_error((val), #val, __FILE__, __LINE__)
#define checkCublasErrors(val) check_cublas_error((val), #val, __FILE__, __LINE__)

void check_cuda_error(cudaError_t result, char const *const func, const char *const file, int const line);
void check_cublas_error(cublasStatus_t status, char const *const func, const char *const file, int const line);

#endif
