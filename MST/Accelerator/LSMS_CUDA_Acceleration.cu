#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <complex.h>
//#include "cudaDoubleComplex.hpp"
//#include "DeviceStorage.hpp"

// CUDA error-checking macro for reliable memory allocation and transfers
#define checkCudaErrors(val) check_cuda_error((val), #val, __FILE__, __LINE__)
void check_cuda_error(cudaError_t result, char const *const func, const char *const file, int const line) {
    if (result) {
        fprintf(stderr, "\nCUDA error at %s:%d code=%d(%s) \"%s\"\n", file, line, 
                        static_cast<unsigned int>(result), cudaGetErrorName(result), func);
        exit(EXIT_FAILURE);
    }
}

void checkCuda(cudaError_t result) {
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
        exit(EXIT_FAILURE);
    }
}

void checkCublas(cublasStatus_t status) {
    if (status != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "cuBLAS Error\n");
        exit(EXIT_FAILURE);
    }
}

// CUDA kernel to create a unit matrix
__global__ void createUnitMatrixKernel(cuDoubleComplex *d_matrix, int n) {
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

// CUDA kernel to copy the upper-left nB x nB block of matrix A of size nA x nA into matrix B
__global__ void copyBlockKernel(cuDoubleComplex *A, int nA, cuDoubleComplex *B, int nB) {
   int row = blockIdx.y * blockDim.y + threadIdx.y;
   int col = blockIdx.x * blockDim.x + threadIdx.x;

   // Check if the current thread is within the n x n block
   if (row < nB && col < nB) {
      B[row * nB + col] = A[row * nA + col];
   }
}

bool initialized = false;
int tau_size = 0;
int mmat_size = 0;

double _Complex  *sineHost;  // mat_id = 1
double _Complex  *jinvHost;  // mat_id = 2
double _Complex  *gijHost;   // mat_id = 3

cuDoubleComplex  *sineDev;
cuDoubleComplex  *jinvDev;
cuDoubleComplex  *gijDev;
cuDoubleComplex  *BigMatDev;
cuDoubleComplex  *BigMatInvDev;
cuDoubleComplex  *block00Dev;
int *pivotArray;
int *infoArray;

void copySubBlockToMatrix(double _Complex *a, int *size_a, int *row, int *col, double _Complex *b, int *size_b) {
   if (*row < 0 || *col < 0) {
      fprintf(stderr, "\nError in copySubBlockToMatrix: invalid row or col value, %d,%d\n",*row,*col); 
      exit(EXIT_FAILURE);
   }
   else if (*row * *size_a > *size_b || *col * *size_a > *size_b) {
      fprintf(stderr, "\nError in copySubBlockToMatrix: row or col value could cause out-of-bound %d,%d\n",
                      *row,*col); 
      exit(EXIT_FAILURE);
   }
   else {
      for (int i=0; i<*size_a; i++) {
         for (int j=0; j<*size_a; j++) {
            b[(*col * *size_a+i) * *size_b + *row * *size_a + j] = a[i * *size_a + j];
         }
      }
   }
}

extern "C"
void init_lsms_gpu_(int *dsize, int *block_size, int *my_pe) {
   if (!initialized) {
      mmat_size = *dsize;
      tau_size = *block_size;

      // This code assumes that each atom in the LIZ has the same scattering matrix size
      // ==============================================================

      if (*my_pe == 0) printf("CUDA memory assigned \n");

      cudaError_t error;
      size_t size_block = tau_size * tau_size * sizeof(cuDoubleComplex);
      size_t size_bigmat = mmat_size * mmat_size * sizeof(cuDoubleComplex);

      checkCuda(cudaMalloc((void**)&BigMatDev, size_bigmat));

      checkCuda(cudaMalloc((void**)&pivotArray, sizeof(int)*mmat_size));

      checkCuda(cudaMalloc((void**)&infoArray, sizeof(int)));

      checkCuda(cudaMalloc((void**)&BigMatInvDev, size_bigmat));

      checkCuda(cudaMalloc((void**)&block00Dev, size_block));

      checkCuda(cudaMalloc((void**)&sineDev, size_bigmat));

      checkCuda(cudaMalloc((void**)&jinvDev, size_bigmat));

      checkCuda(cudaMalloc((void**)&gijDev, size_bigmat));

      sineHost = (double _Complex *) malloc(size_bigmat);
      memset(sineHost, 0, size_bigmat);

      jinvHost = (double _Complex *) malloc(size_bigmat);
      memset(jinvHost, 0, size_bigmat);

      gijHost = (double _Complex *) malloc(size_bigmat);
      memset(gijHost, 0, size_bigmat);

      initialized = true;
   }
}

extern "C"
void finalize_lsms_gpu_() {
   if (initialized) {
      free(sineHost);
      free(jinvHost);
      free(gijHost);
      cudaFree(BigMatDev);
      cudaFree(pivotArray);
      cudaFree(infoArray);
      cudaFree(BigMatInvDev);
      cudaFree(block00Dev);
      cudaFree(sineDev);
      cudaFree(jinvDev);
      cudaFree(gijDev);
   }
   initialized = false;
}

extern "C"
void init_bigmatrix_gpu_(int *b_size) {
   if (*b_size != mmat_size) {
      fprintf(stderr,"\nError: b_size <> mmat_size, %d,%d\n",*b_size,mmat_size);
      exit(EXIT_FAILURE);
   }

   // Define kernel launch parameters and set BigMatDev to be a unit matrix
   // ========================================
   int threads_per_block = 512;
   int num_blocks = (mmat_size*mmat_size + threads_per_block - 1) / threads_per_block;
   createUnitMatrixKernel<<<num_blocks, threads_per_block>>>(BigMatDev, mmat_size);
   checkCudaErrors(cudaPeekAtLastError());
}

extern "C"
void push_submatrix_gpu_(int *mat_id, int *row, int *col, double _Complex *mat, int *msize ) {
   if (*msize != tau_size) {
      fprintf(stderr,"\nError: msize <> tau_size, %d,%d\n",*msize,tau_size);
      exit(EXIT_FAILURE);
   }
   else if (*mat_id == 1) {
      copySubBlockToMatrix(mat,msize,row,col,sineHost,&mmat_size);
   }
   else if (*mat_id == 2) {
      copySubBlockToMatrix(mat,msize,row,col,jinvHost,&mmat_size);
   }
   else if (*mat_id == 3) {
      copySubBlockToMatrix(mat,msize,row,col,gijHost,&mmat_size);
   }
   else {
      fprintf(stderr,"\nError: invalid matrix ID, %d\n",*mat_id);
      exit(EXIT_FAILURE);
   }
}

extern "C"
void push_gij_matrix_gpu_(int *cant, int *row, int *col, double _Complex *gij, int *kkri) {
   if (*cant != 1 && *cant != 2) {
      fprintf(stderr,"\nError: cant <> 1 and 2, %d\n",*cant);
      exit(EXIT_FAILURE);
   } 
   else if (*kkri * *cant != tau_size) {
      fprintf(stderr,"\nError: kkri*cant <> tau_size, %d,%d\n",*kkri * *cant,tau_size);
      exit(EXIT_FAILURE);
   }

   if (*cant == 1) {
      copySubBlockToMatrix(gij,kkri,row,col,gijHost,&mmat_size);
   }
   else {
      double _Complex *dijHost;
      dijHost = (double _Complex *) malloc(sizeof(double _Complex)*tau_size*tau_size);
      memset(dijHost, 0, sizeof(double _Complex)*tau_size*tau_size);
      int diag = 0;
      copySubBlockToMatrix(gij,kkri,&diag,&diag,dijHost,&tau_size);
      diag = 1;
      copySubBlockToMatrix(gij,kkri,&diag,&diag,dijHost,&tau_size);
      copySubBlockToMatrix(dijHost,&tau_size,row,col,gijHost,&mmat_size);
      free(dijHost);
   }
}

extern "C"
void commit_to_gpu_(int *mat_id) {
   size_t size = sizeof(cuDoubleComplex)*mmat_size*mmat_size;
   if (*mat_id == 1) {
      checkCuda(cudaMemcpy(sineDev, sineHost, size, cudaMemcpyHostToDevice));
   }
   else if (*mat_id == 2) {
      checkCuda(cudaMemcpy(jinvDev, jinvHost, size, cudaMemcpyHostToDevice));
   }
   else if (*mat_id == 3) {
      checkCuda(cudaMemcpy(gijDev, gijHost, size, cudaMemcpyHostToDevice));
   }
   else {
      fprintf(stderr,"\nError: invalid matrix ID, %d\n",*mat_id);
      exit(EXIT_FAILURE);
   }
}

extern "C"
void construct_bigmatrix_gpu_(double _Complex *neg_kappa_inv) {
    cuDoubleComplex *jig_d;
    cuDoubleComplex minus_one = make_cuDoubleComplex(-1.0, 0.0);
    cuDoubleComplex one = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex zero = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex alpha = make_cuDoubleComplex(creal(*neg_kappa_inv),cimag(*neg_kappa_inv));

    size_t size = mmat_size * mmat_size * sizeof(cuDoubleComplex);

    // Allocate device memory
    checkCuda(cudaMalloc((void**)&jig_d, size));

    // Create cuBLAS handle
    cublasHandle_t handle;
    checkCublas(cublasCreate(&handle));

    // Compute jig = jinv * gij
    checkCublas(cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                            mmat_size, mmat_size, mmat_size, &one,
                            jinvDev, mmat_size,
                            gijDev, mmat_size,
                            &zero, jig_d, mmat_size));

    // Compute BigMat = 1 - jig * sine / kappa
    checkCublas(cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                            mmat_size, mmat_size, mmat_size, &alpha,
                            jig_d, mmat_size,
                            sineDev, mmat_size,
                            &one, BigMatDev, mmat_size));

    // Cleanup
    cudaFree(jig_d);
    cublasDestroy(handle);
}

extern "C"
void invert_bigmatrix_gpu_(double _Complex *block00, int *block_size) {
   static cudaError_t error;
   static int Lwork;

   if (tau_size != *block_size) {
      fprintf(stderr,"\nError: tau_size <> block_size, %d,%d\n",tau_size,*block_size);
      exit(EXIT_FAILURE);
   }

   cusolverDnHandle_t cusolverHandle;
   cusolverStatus_t cusolverStatus;
   cusolverDnCreate(&cusolverHandle);
   cusolverDnZgetrf_bufferSize(cusolverHandle, mmat_size, mmat_size, BigMatDev, mmat_size, &Lwork);

   // Create a working array on the device
   // ========================================
   cuDoubleComplex  *workArray;
   checkCuda(cudaMalloc((void**)&workArray, Lwork*sizeof(cuDoubleComplex)));

   // Define kernel launch parameters and create a unit matrix on device
   // ========================================
   int threads_per_block = 512;
   int num_blocks = (mmat_size*mmat_size + threads_per_block - 1) / threads_per_block;
   createUnitMatrixKernel<<<num_blocks, threads_per_block>>>(BigMatInvDev, mmat_size);
   checkCudaErrors(cudaPeekAtLastError());

   // Perform matrix inverse on the device
   // ============================================
   cusolverStatus = cusolverDnZgetrf(cusolverHandle, mmat_size, mmat_size, BigMatDev, mmat_size, 
                                     workArray, pivotArray, infoArray);
   if (cusolverStatus != CUSOLVER_STATUS_SUCCESS) {
      fprintf(stderr,"\nError: cuSOLVER ZGETRF UNSUCCESSFUL! \n");
   }
   cusolverStatus = cusolverDnZgetrs(cusolverHandle, CUBLAS_OP_N, mmat_size, mmat_size, BigMatDev, mmat_size, 
                                     pivotArray, BigMatInvDev, mmat_size, infoArray); 
   if (cusolverStatus != CUSOLVER_STATUS_SUCCESS) {
      fprintf(stderr,"\nError: cuSOLVER ZGETRS UNSUCCESSFUL! \n");
   }

   // --------------------------------------------
   // Aggregate the data into block00Dev array and then make one cudaMemcpy call
   // --------------------------------------------
   // Define grid and block dimensions for the kernel
   dim3 threadsPerBlock(32, 32); // Using a 32x32 block
   dim3 blocksPerGrid((tau_size + threadsPerBlock.x - 1) / threadsPerBlock.x,
                      (tau_size + threadsPerBlock.y - 1) / threadsPerBlock.y);

   // Launch the kernel
   // ============================================
   copyBlockKernel<<<blocksPerGrid, threadsPerBlock>>>(BigMatInvDev, mmat_size, block00Dev, tau_size);
   checkCudaErrors(cudaPeekAtLastError());
   checkCuda(cudaMemcpy(block00, block00Dev, sizeof(cuDoubleComplex)*tau_size*tau_size, cudaMemcpyDeviceToHost));

   // clean up
   // ============================================
   cusolverDnDestroy(cusolverHandle);
   cudaFree(workArray);
}
