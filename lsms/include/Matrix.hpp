// \section{Matrix Class}
// We use an addapted version of the Matrix class from {\tt psimag}.
// The memory layout of this matrix class is designed to be the same in {\bf Fortran}.

#ifndef LSMS_MATRIX_H
#define LSMS_MATRIX_H

#include <string.h>
#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <typeinfo>
// #include "PSIMAGAssert.h"

#ifdef BUILDKKRMATRIX_GPU
#include <cuda_runtime.h>
#endif

  template <typename T>
  class Matrix {
  public:
    typedef Matrix<T>    ThisType;
    typedef size_t       size_type;
    typedef T            value_type;
    typedef T*           iterator;
    typedef const T*     const_iterator;
    typedef T&           reference;
    typedef const T&     const_reference;

    // The basic constructor builds an empty matrix, that needs to be resized before use.
    Matrix() : nRow(0), nCol(0), lDim(0), owner(true), data(0), physicalSize(0) {;}

    // We also can construct a matrix with |nRows| rows and |nCols| the leading dimension |ldim| defaults to |nRows|, but can be set to any larger value to get the desired memory alignement.
    Matrix(size_type nRows,size_type nCols,size_type ldim=0,const T& val=T()) {
      if(ldim<nRows) 
	ldim=nRows;
      if(bool(ldim*nCols)) 
	{
	  nRow=nRows; nCol=nCols; lDim=ldim;
	  owner=true;
          physicalSize=lDim*nCol;
	  data = new T[physicalSize];
	}
      else { 
	nRow=0; nCol=0; lDim=0; owner=true; data=0; physicalSize=0;
      } 
      if(val==T(0))
	memset(data,0,sizeof(T)*lDim*nCol);
      else
	for(size_type i=0;i<nCol*lDim;i++) 
	  data[i]=val;
    }

    // Also a special constructor is provided to map a memory region to this matrix class. The user is responsible for the proper allocation and dealocation of memory and the proper dimensions.
    /* matrix that does not own data (user is responsible for proper dim. */
    Matrix(size_type nRows, size_type nCols, T* dat, size_type ldim=0) 
      : data(dat) {
      if(ldim<nRows) ldim=nRows;
      if(bool(nRows*nCols))
	{
	  nRow=nRows; nCol=nCols; lDim=ldim; physicalSize=lDim*nCol;
	  owner=false;
	}
      else { nRow=0; nCol=0; lDim=0; owner=false; physicalSize=0;}
    }

// The final constructor makes a deep copy of another matrix.
    /* makes a deep copy of mat */
    Matrix(const Matrix<T> & mat) 
      : nRow(mat.nRow), nCol(mat.nCol), lDim(mat.nRow),owner(true) {
      if(bool(nRow*nCol)) {
	owner=true;
        physicalSize=lDim*nCol;
	data = new T[physicalSize];
	memcpy(data,mat.data,sizeof(T)*physicalSize);
      }
      else { nRow=0; nCol=0; lDim=0; owner=true; data=0; physicalSize=0;} 
    }

    // <Resizing a Matrix@>;
    void resize(size_type m,size_type n,size_type ldim=0) {
      if(ldim<m) ldim=m;
      if(!owner && (ldim*n>physicalSize)){std::logic_error("matrix not locally owned in T& Matrix<T>::resize");}
      if(ldim*n <= physicalSize) {
        if(bool(ldim*n) && data) {
          nCol = n; nRow = m; lDim = ldim;
        }
        else {
          nRow=0; nCol=0; lDim=0; if(owner) {delete [] data; data=0; physicalSize=0;}
        }
      }
      else {
        if(owner && data) delete [] data;
        nCol = n; nRow = m; lDim = ldim; owner=true;
        physicalSize=lDim*nCol;
        data = new T[physicalSize];
      }
    }

  void retarget(size_type nRows, size_type nCols, T* dat, size_type ldim=0) {
      if(owner && data) delete [] data;
      data=dat;
      if(ldim<nRows) ldim=nRows;
      if(bool(nRows*nCols))
        {
          nRow=nRows; nCol=nCols; lDim=ldim; physicalSize=lDim*nCol;
          owner=false;
        }
      else { nRow=0; nCol=0; lDim=0; owner=false; physicalSize=0;}
    }


    // <Matrix Destructor@>;
    // Finaly we need a destructor to free the memory allocated by the matrix.
    ~Matrix() {
      if(owner && data) delete [] data;
    }

    // <Matrix Access@>;
    // We provide two methods to access the lements of a matrix. The most natural one uses |operator()| with row and collumn arguments and looks simimar to a matrix access in {\bf Fortran}.
    inline T& operator() (size_type i, size_type j) {
      // ASSERT(i<nRow,
      //        std::range_error("i>=n_row in T& Matrix<T>::operator()(size_type i,size_type j"));
      // ASSERT(j<nCol,
      //        std::range_error("j>=n_col in T& Matrix<T>::operator()(size_type i,size_type j"));
      return data[j*lDim+i];
    }

  // The second way to access the elements of a matrix uses the familair C |[]| operator and vies the |data| array as a continous block of memory.
    inline T& operator[](size_type i) {
      return data[i];
    }

  // We provide a few functions to return information about the matrix.
    size_type size() const { return lDim*nCol; }
    size_type n_row() const { return nRow; }
    size_type n_col() const { return nCol; }
    size_type l_dim() const { return lDim; }
    
  // \subsection{Operations on Matrices}

  // Assignments and copy:
    Matrix<T> &operator=(const Matrix<T>& mat)
    {
      if (this == &mat) return *this;   // Gracefully handle self assignment[12.1]
      if(owner && data) delete [] data;
      nRow=mat.n_row(), nCol=mat.n_col(), lDim=mat.n_row(),owner=true;
      if(bool(nRow*nCol)) {
        physicalSize=lDim*nCol;
	data=new T[physicalSize];
	if(mat.lDim==lDim)
	  memcpy(data,mat.data,sizeof(T)*nCol*lDim);
	else
	  for(size_type j=0;j<nCol;j++) for(size_type i=0;i<nRow;i++)
//					  data[j+i*lDim] = mat(i,j);
                                          data[i+j*lDim] = mat.data[i+j*mat.l_dim()];
      }
      else { nRow=0; nCol=0; lDim=0; owner=true; data=0; physicalSize=0;}
      return *this;
    }

    void copy(const Matrix<T>& A) {
      if(A.l_dim()==lDim && A.n_col()==nCol) {
        for(size_t i=0; i<nCol*lDim; i++) data[i]=A[i];
      }
      else if (A.n_row()==nRow && A.n_col()==nCol){
        for(size_t j=0; j<nCol; j++)
          for(size_t i=0; i<nRow; i++)
//            data[j+lDim*i]=A(i,j);
            data[i+j*lDim] = A.data[i+j*A.l_dim()];
      }
      else std::range_error("matrix sizes don't match in Matrix<T>::copy");
    }

/*
    bool operator==(const Matrix<T>&a)
    {
      bool r=true;
      if(a.n_row()!=nRow || a.n_col()!=nCol) return false;
      for(size_t i=0; i<a.n_row(); i++)
        for(size_t j=0; j<a.n_col(); j++)
          r=r&&(a(i,j)==data[i+j*lDim]);
      return r;
    }

    bool operator!=(const Matrix<T>&a)
    {
      return !(*this==a);
    }
*/

    inline
    ThisType& operator = (const T& val)
    {
      for(size_type i=0;i<nCol*lDim;i++) 
	data[i]=val;
      return *this;
    }

    // Scaling and adding:
    void scale(T a) {
      for(size_t i=0; i<nCol*lDim; i++) data[i]*=a;
    }

    void add(Matrix<T> &A) {
      if(A.l_dim()==lDim && A.n_col()==nCol) {
        for(size_t i=0; i<nCol*lDim; i++) data[i]+=A[i];
      }
      else if (A.n_row()==nRow && A.n_col()==nCol){
        for(size_t j=0; j<nCol; j++)
            for(size_t i=0; i<nRow; i++)
            data[j*lDim+i]+=A(i,j);
      }
      else std::range_error("matrix sizes don't match in Matrix<T>::add");
    }

    void addScaled(Matrix<T> &A, T a) {
      if(A.l_dim()==lDim && A.n_col()==nCol) {
        for(size_t i=0; i<nCol*lDim; i++) data[i]+=a*A[i];
      }
      else if (A.n_row()==nRow && A.n_col()==nCol){
        for(size_t j=0; j<nCol; j++)
          for(size_t i=0; i<nRow; i++)
            data[j*lDim+i]+=a*A(i,j);
      }
      else std::range_error("matrix sizes don't match in Matrix<T>::addScaled");
    }

    void pinMemory() {
#ifdef BUILDKKRMATRIX_GPU
      cudaHostRegister(data,physicalSize*sizeof(T),0);
#endif
    }

    void unpinMemory() {
#ifdef BUILDKKRMATRIX_GPU
      cudaHostUnregister(data);
#endif
    }

  private:
    size_type nRow,nCol,lDim, physicalSize;
    bool owner;
    T* data;
    // T** col;
  };

    // \subsection{Matrix functions}
  // Here we use the {\tt BLAS} and {\tt LAPACK} interfaces from PsiMag.

#include "BLAS.hpp"
#include "LAPACK.hpp"

// Matrix multiplication. Here we use the {\tt BLAS} |xGEMM| interface. The operation performed is
// $$ C \leftarrow \alpha AB + \beta C$$

template<class T>
void gemm(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, 
          const T& alpha=1, const T& beta=0, char opA = 'N', char opB = 'N') 
{
  // ASSERT(A.n_row() == C.n_row(),
  //        std::range_error("A.n_row() != C.n_row() in GEMM for Matrix<T>"));
  // ASSERT(A.n_col() == C.n_col(),
  //        std::range_error("A.n_col() != C.n_col() in GEMM for Matrix<T>"));
  // ASSERT(A.n_col() == B.n_row(),
  //        std::range_error("A.n_col() != B.n_row() in GEMM for Matrix<T>"));
  BLAS::GEMM(char(opA),
             char(opB),
             static_cast<int>(A.n_row()),
             static_cast<int>(B.n_col()),
             static_cast<int>(A.n_col()),
             alpha,
             &A(0,0),
             static_cast<int>(A.l_dim()),
             &B(0,0),
             static_cast<int>(B.l_dim()),
             beta,
             &C(0,0),
             static_cast<int>(C.l_dim()));
}
	  
#endif
