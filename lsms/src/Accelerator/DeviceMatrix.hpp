#ifndef LSMS_DEVICE_MATRIX_HPP
#define LSMS_DEVICE_MATRIX_HPP

#include "Matrix.hpp"
#include "cudaCheckError.hpp"
#include <cstdio>
  template <class T>
  struct DeviceMatrix {
    public:
      typedef size_t       size_type;
      __inline__ DeviceMatrix() : nRow(0), nCol(0), lDim(0), data(0), owner(0)  {}
      
      __inline__ DeviceMatrix(Matrix<T>& in) : data(0) {
        *this=in;
      }
      
      __inline__ DeviceMatrix(size_type nRows,size_type nCols,size_type ldim=0) : nRow(nRows), nCol(nCols), lDim(ldim)  {
        if(lDim<nRow)
          ldim=nRow;
        size_type num_bytes=sizeof(T)*lDim*nCol;
        allocate(num_bytes);
      }

      __inline__ ~DeviceMatrix() {
        free();
      }

      __inline__ void clear() { free(); }

      // <DeviceMatrix Access@>;
      // We provide two methods to access the lements of a matrix. The most natural one uses |operator()| with row and collumn arguments and looks simimar to a matrix access in {\bf Fortran}.
      __inline__ __device__ T& operator() (size_type i, size_type j) {
        return data[j*lDim+i];
      }

      // The second way to access the elements of a matrix uses the familair C |[]| operator and vies the |data| array as a continous block of memory.
      __inline__ __device__ T& operator[](size_type i) {
        return data[i];
      }

      __inline__ DeviceMatrix<T> &operator=(Matrix<T>& mat) {
        //TODO BACK TO ASYNC
        copy_async(mat,0);
        cudaCheckError();
        return *this;
      }

      __inline__ void copy(Matrix<T> &mat) {
        size_type curSize=lDim*nCol;
        nRow=mat.n_row(); nCol=mat.n_col(); lDim=mat.l_dim(); 
        size_type num_bytes=sizeof(T)*lDim*nCol;
        
        printf("lDim: %d, nCol: %d, nRow: %d, num_bytes: %d, curSize: %d\n",lDim,nCol,nRow,num_bytes,curSize);
        if(lDim*nCol!=curSize) {
          if((nRow*nCol)>0) {
            allocate(num_bytes);
          }
          else { free(); nRow=0; nCol=0; lDim=0; } 
        }
        cudaMemcpy(data,&mat[0],num_bytes,cudaMemcpyHostToDevice);
      }

      __inline__ void copy_async(Matrix<T> &mat, cudaStream_t s) {
        size_type curSize=lDim*nCol;
        nRow=mat.n_row(); nCol=mat.n_col(); lDim=mat.l_dim(); 
        size_type num_bytes=sizeof(T)*lDim*nCol;

        if(lDim*nCol!=curSize) {
          if(data!=0) free();
          if((nRow*nCol)>0) {
            allocate(num_bytes);
          }
          else { nRow=0; nCol=0; lDim=0; data=0;} 
        }
        cudaMemcpyAsync(data,&mat[0],num_bytes,cudaMemcpyHostToDevice,s);
      }


      // We provide a few functions to return information about the matrix.
      __inline__ __host__ __device__ size_type size() const { return lDim*nCol; }
      __inline__ __host__ __device__ T* raw() const { return data; }
      __inline__ __host__ __device__ size_type n_row() const { return nRow; }
      __inline__ __host__ __device__ size_type n_col() const { return nCol; }
      __inline__ __host__ __device__ size_type l_dim() const { return lDim; }

    private:
      size_type lDim;
      size_type nRow,nCol;

      T* data;
      DeviceMatrix<T> *owner;
      
      __inline__ void allocate(size_type num_bytes) {
        if(num_bytes>0) {
          free();
          owner=this;
          cudaMalloc(&data,num_bytes);
          cudaCheckError();
        }
        else {      // DANGEROUS!! Might cause memory leak!
          owner = 0;
          data = 0;
          nRow = 0;
          nCol = 0;
          lDim = 0;
        }
      }
      __inline__ void free() {
        if(owner==this && data!=0) {
          cudaFree(data);
        }
        owner=0;
        data=0;
      }
  };

#endif

