#ifndef LSMS_DEVICE_ARRAY3D_HPP
#define LSMS_DEVICE_ARRAY3D_HPP

#include "Array3d.hpp"

  template <class T>
  class DeviceArray3d {
    public:
      typedef size_t       size_type;

      // The basic constructor builds an empty matrix, that needs to be resized before use.
      __inline__ DeviceArray3d() : nRow(0), nCol(0), nSlice(0), lDim1(0), lDim2(0), lDim12(0), data(0), owner(0) {}

      __inline__ DeviceArray3d(Array3d<T>& in): data(0) {
        *this=in;
      }

      __inline__ DeviceArray3d(size_type nRows,size_type nCols, size_type nSlices, size_type ldim1=0, size_type ldim2=0) {
        nRow=nRows; nCol=nCols; nSlice=nSlices; lDim1=ldim1; lDim2=ldim2; lDim12=lDim1*lDim2;
        if(ldim1<nRows)
          ldim1=nRows;
        if(ldim2<nCols)
          ldim2=nCols;
        size_type num_bytes=sizeof(T)*lDim12*nSlice;
        allocate(num_bytes);
      }

      // <DeviceArray3d Destructor@>;
      // Finaly we need a destructor to free the memory allocated by the array.
      __inline__ ~DeviceArray3d() {
        free();
      }

      // <Matrix Access@>;
      // We provide two methods to access the elements of a matrix. The most natural one uses |operator()| with row and collumn arguments and looks simimar to a matrix access in {\bf Fortran}.
      __inline__ __device__ T& operator() (size_type i, size_type j, size_type k) {
        return data[k*lDim12+j*lDim1+i];
      }

      // The second way to access the elements of a matrix uses the familair C |[]| operator and vies the |data| array as a continous block of memory.
      __inline__ __device__ T& operator[](size_type i) {
        return data[i];
      }

      // We provide a few functions to return information about the matrix.
      __inline__ __host__ __device__ size_type size() const { return lDim12*nSlice; }
      __inline__ __host__ __device__ T* raw() const { return data; }
      __inline__ __host__ __device__ size_type n_row() const { return nRow; }
      __inline__ __host__ __device__ size_type n_col() const { return nCol; }
      __inline__ __host__ __device__ size_type n_slice() const {return nSlice;}
      __inline__ __host__ __device__ size_type l_dim1() const { return lDim1; }
      __inline__ __host__ __device__ size_type l_dim2() const { return lDim2;}

      // Assignments and copy:
      __inline__ DeviceArray3d<T> &operator=(Array3d<T>& a)
      {
        copy_async(a,0);
        return *this;
      }
      __inline__ void copy(Array3d<T> &a) {
        size_type curSize=lDim12*nSlice;
        nRow=a.n_row(); nCol=a.n_col(); nSlice=a.n_slice(); lDim1=a.l_dim1(); lDim2=a.l_dim2(); lDim12=lDim1*lDim2; 
        size_type num_bytes=sizeof(T)*lDim12*nSlice;

        if(curSize!=lDim12*nSlice) {
          if(nRow*nCol*nSlice>0) {
            allocate(num_bytes);
          }
          else { nRow=0; nCol=0; nSlice=0; lDim1=0; lDim2=0; lDim12=0; data=0; }
        }
        cudaMemcpy(data,&a[0],num_bytes,cudaMemcpyHostToDevice);
      }
      
      __inline__ void copy_async(Array3d<T> &a, cudaStream_t s) {
        size_type curSize=lDim12*nSlice;
        nRow=a.n_row(); nCol=a.n_col(); nSlice=a.n_slice(); lDim1=a.l_dim1(); lDim2=a.l_dim2(); lDim12=lDim1*lDim2; 
        size_type num_bytes=sizeof(T)*lDim12*nSlice;

        if(curSize!=lDim12*nSlice) {
          if(nRow*nCol*nSlice>0) {
            allocate(num_bytes);
          }
          else { nRow=0; nCol=0; nSlice=0; lDim1=0; lDim2=0; lDim12=0; data=0; }
        }
        cudaMemcpyAsync(data,&a[0],num_bytes,cudaMemcpyHostToDevice,s);
      }

    private:
      size_type nRow,nCol,nSlice,lDim1,lDim2,lDim12;
      T* data;
      DeviceArray3d<T> *owner;

      __inline__ void allocate(size_type num_bytes) {
        if(num_bytes>0) {
          free();
          owner=this;
          cudaMalloc(&data,num_bytes);
          cudaCheckError();
        }
        else {      // DANGEROUS!! Might cause memory leak!
          nRow = 0;
          nCol = 0;
          nSlice = 0;
          lDim1 = 0;
          lDim2 = 0;
          lDim12 = 0;
          data = 0;
          owner = 0;
        }
      }
      __inline__ void free() {
        if(owner==this && data!=0) {
          cudaFree(data);
          owner=0;
          data=0;
        }
      }
  };
	 
#endif
 
