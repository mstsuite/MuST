#include <math.h>
#include <cmath>
#include <complex.h>
#include "cuComplex.h"

// CUDA kernel to create a unit matrix
__global__ void createUnitMatrixKernel(cuDoubleComplex *matrix_d, int n) {
    // Calculate the unique thread index
    // int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Boundary check to ensure the thread does not go out of bounds
    // if (idx < n * n) {
    for (size_t idx = tid; idx < n*n; idx += blockDim.x*gridDim.x) {
        // Calculate the row and column from the linear index
        int row = idx / n;
        int col = idx % n;

        // If the thread is on the main diagonal, set the value to 1.0 + 0.0i
        if (row == col) {
            matrix_d[idx] = make_cuDoubleComplex(1.0, 0.0);
        }
        // Otherwise, set the value to 0.0 + 0.0i
        else {
            matrix_d[idx] = make_cuDoubleComplex(0.0, 0.0);
        }
    }  
    __syncthreads();
}

// CUDA kernel to copy the upper-left nB x nB block of matrix A of size nA x nA into matrix B
__global__ void copyBlockKernel(cuDoubleComplex *A, int nA, cuDoubleComplex *B, int nB) {
   // int row = blockIdx.y * blockDim.y + threadIdx.y;
   // int col = blockIdx.x * blockDim.x + threadIdx.x;
   int tid_row = blockIdx.y * blockDim.y + threadIdx.y;
   int tid_col = blockIdx.x * blockDim.x + threadIdx.x;

   // Check if the current thread is within the n x n block
   // if (row < nB && col < nB) {
   for (size_t row = tid_row; row < nB; row += blockDim.y*gridDim.y) {
      for (size_t col = tid_col; col < nB; col += blockDim.x*gridDim.x) {
         B[row * nB + col] = A[row * nA + col];
      }
   }
}

// CUDA kernel to copy the upper-left nB x nB block of matrix A of size nA x nA into matrix B
__global__ void copySubBlockKernel(cuDoubleComplex *A, int nA, int i, int j, cuDoubleComplex *B, int nB) {
   int col = blockIdx.y * blockDim.y + threadIdx.y;
   int row = blockIdx.x * blockDim.x + threadIdx.x;

   // Check if the current thread is within the n x n block
   if (row < nB && col < nB) {
      B[col * nB + row] = A[j*nA*nB + col*nA + i*nB + row];
   }
}

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
