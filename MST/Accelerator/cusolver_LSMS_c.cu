#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <complex.h>
//#include "cudaDoubleComplex.hpp"
//#include "DeviceStorage.hpp"

// CUDA error-checking macro for reliable memory allocation and transfers
#define cu_checkCudaErrors(val) cu_check_cuda((val), #val, __FILE__, __LINE__)
void cu_check_cuda(cudaError_t result, char const *const func, const char *const file, int const line) {
    if (result) {
        fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\"\n", file, line, 
                        static_cast<unsigned int>(result), cudaGetErrorName(result), func);
        exit(EXIT_FAILURE);
    }
}

// CUDA kernel to create a unit matrix
__global__ void cu_createUnitMatrixKernel(cuDoubleComplex *d_matrix, int n) {
    // Calculate the unique thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Boundary check to ensure the thread does not go out of bounds
    if (idx < n * n) {
        // Calculate the row and column from the linear index
        int row = idx / n;
        int col = idx % n;

        // If the thread is on the main diagonal, set the value to 1.0 + 0.0i
        if (row == col) {
            d_matrix[idx] = make_cuDoubleComplex(1.0, 0.0);
        }
        // Otherwise, set the value to 0.0 + 0.0i
        else {
            d_matrix[idx] = make_cuDoubleComplex(0.0, 0.0);
        }
    }
}

bool cu_initialized = false;
cuDoubleComplex  *aDev;
cuDoubleComplex  *aInvDev;
cuDoubleComplex  *tau00Dev;
int *cu_pivotArray;
int *cu_infoArray;
int bssq = 0;

// CUDA kernel to copy the upper-left nB x nB block of matrix A of size nA x nA into matrix B
__global__ void cu_copyBlockKernel(cuDoubleComplex *A, int nA, cuDoubleComplex *B, int nB) {
   int row = blockIdx.y * blockDim.y + threadIdx.y;
   int col = blockIdx.x * blockDim.x + threadIdx.x;

   // Check if the current thread is within the n x n block
   if (row < nB && col < nB) {
      B[row * nB + col] = A[row * nA + col];
   }
}

extern "C"
void cu_init_lsms_gpu_(int *m, int *block_size) {
   if (!cu_initialized) {
      cudaError_t error;

      printf("CUDA memory assigned \n");
      error=cudaMalloc((void**)&aDev,  sizeof(cuDoubleComplex)* *m * *m);
      if (error != cudaSuccess) fprintf(stderr,"\nError1: %s\n",cudaGetErrorString(error));

      error=cudaMalloc((void**)&cu_pivotArray,  sizeof(int) * *m);
      if (error != cudaSuccess) fprintf(stderr,"\nError2: %s\n",cudaGetErrorString(error));

      error=cudaMalloc((void**)&cu_infoArray,  sizeof(int));
      if (error != cudaSuccess) fprintf(stderr,"\nError3: %s\n",cudaGetErrorString(error));

      error=cudaMalloc(&aInvDev, sizeof(cuDoubleComplex)* *m * *m);
      if (error != cudaSuccess) fprintf(stderr,"\nError4: %s\n",cudaGetErrorString(error));

      bssq = *block_size * *block_size;
      error=cudaMalloc(&tau00Dev, sizeof(cuDoubleComplex)*bssq);
      if (error != cudaSuccess) fprintf(stderr,"\nError5: %s\n",cudaGetErrorString(error));

      cu_initialized = true;
   }
}

extern "C"
void cusolver_lsms_c_(int *m, double _Complex *a, int *block_size, double _Complex *b) {
   // Added on 9/16/2025 ===
   // ======================
   static cudaError_t error;
   static int Lwork;
   // float time_copyin=0;
   // float time_copyout=0;
   // float time_compute=0;

   if (!cu_initialized) {
      fprintf(stderr,"\nError: %s\n","lsms_gpu_init() needs to be called first");
      exit(1);
   }
   else if (bssq != *block_size * *block_size) {
      fprintf(stderr,"\nError: bssq <> block_size**2, %d,%d\n",bssq,*block_size * *block_size);
      exit(1);
   }

   //cudaEvent_t start, stop;
   //cudaEventCreate(&start);
   //cudaEventCreate(&stop);
   cusolverDnHandle_t cusolverHandle;
   cusolverStatus_t cusolverStatus;
   cusolverDnCreate(&cusolverHandle);
   cusolverDnZgetrf_bufferSize(cusolverHandle, *m, *m, aDev, *m, &Lwork);
   //printf("Lwork is %d\n",Lwork);

   cuDoubleComplex  *workArray;
   error=cudaMalloc((void**)&workArray, Lwork*sizeof(cuDoubleComplex));
   if (error != cudaSuccess) fprintf(stderr,"\nError5: %s\n",cudaGetErrorString(error));

   //cudaEventRecord(start); 
   error = cudaMemcpy(aDev, a, sizeof(cuDoubleComplex)* *m * *m, cudaMemcpyHostToDevice);
   if (error != cudaSuccess) fprintf(stderr,"\nError6: %s\n",cudaGetErrorString(error));
   //cudaEventRecord(stop);
   //cudaEventSynchronize(stop);
   //cudaEventElapsedTime(&time_copyin, start, stop);

   // We are create a unit matrix on device, instead of copying it from the host
   // ========================================
   //* error = cudaMemcpy(aInvDev, tau00, sizeof(cuDoubleComplex)* *m * *m, cudaMemcpyHostToDevice);
   //* if (error != cudaSuccess) fprintf(stderr,"\nError7: %s\n",cudaGetErrorString(error));
   // ----------------------------------------
   // Define kernel launch parameters
   int threads_per_block = 512;
   int num_blocks = (*m * *m + threads_per_block - 1) / threads_per_block;
   cu_createUnitMatrixKernel<<<num_blocks, threads_per_block>>>(aInvDev, *m);
   cu_checkCudaErrors(cudaPeekAtLastError());
   // ========================================

   //cudaEventRecord(start);
   cusolverStatus = cusolverDnZgetrf(cusolverHandle, *m, *m, aDev, *m, workArray, cu_pivotArray, cu_infoArray);
   //if (cusolverStatus == CUSOLVER_STATUS_SUCCESS)
   //  printf("cuSOLVER ZGETRF SUCCESSFUL! \n");
   //else
   //  printf("cuSOLVER ZGETRF UNSUCCESSFUL! \n");

   cusolverStatus = cusolverDnZgetrs(cusolverHandle,CUBLAS_OP_N,*m,*m,aDev,*m, cu_pivotArray,aInvDev,*m,cu_infoArray); 
   //if (cusolverStatus == CUSOLVER_STATUS_SUCCESS)
   //  printf("cuSOLVER ZGETRS SUCCESSFUL! \n");
   //else
   //  printf("cuSOLVER ZGETRS UNSUCCESSFUL! \n");
   //cudaEventRecord(stop);
   //cudaEventSynchronize(stop);
   //cudaEventElapsedTime(&time_compute, start, stop);     

   // --------------------------------------------
   // Aggregate the data into tau00Dev array and then make one cudaMemcpy call
   // --------------------------------------------
   // Define grid and block dimensions for the kernel
   dim3 threadsPerBlock(32, 32); // Using a 32x32 block
   dim3 blocksPerGrid((*block_size + threadsPerBlock.x - 1) / threadsPerBlock.x,
                      (*block_size + threadsPerBlock.y - 1) / threadsPerBlock.y);

   // Launch the kernel
   cu_copyBlockKernel<<<blocksPerGrid, threadsPerBlock>>>(aInvDev, *m, tau00Dev, *block_size);
   cudaMemcpy(b, tau00Dev, sizeof(cuDoubleComplex)*bssq, cudaMemcpyDeviceToHost);
   // ============================================
   //cudaEventRecord(stop);
   //cudaEventSynchronize(stop);
   //cudaEventElapsedTime(&time_copyout, start, stop);

   //Print the time (in ms) for GPU data transfer and GPU compute
   //printf("Time for copyin: %f\tfor copyout: %f\tfor compute inverse: %f\n",
   //                       time_copyin*0.001,time_copyout*0.001,time_compute*0.001);
   //clean up
   cusolverDnDestroy(cusolverHandle);
   cudaFree(workArray);
}
