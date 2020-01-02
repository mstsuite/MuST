// \section{3 dim. Array Class}
// We use an adapted the Matrix class to an additional dimension.
// The memory layout of this class is designed to be the same in {\bf Fortran}.

#ifndef LSMS_ARRAY3D_H
#define LSMS_ARRAY3D_H

#include <string.h>
#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <typeinfo>
// #include "PSIMAGAssert.h"

  template <typename T>
  class Array3d {
  public:
    typedef Array3d<T>    ThisType;
    typedef size_t       size_type;
    typedef T            value_type;
    typedef T*           iterator;
    typedef const T*     const_iterator;
    typedef T&           reference;
    typedef const T&     const_reference;

    // The basic constructor builds an empty matrix, that needs to be resized before use.
    Array3d() : nRow(0), nCol(0), nSlice(0), lDim1(0), lDim2(0), lDim12(0), owner(true), data(0) {;}

    // We also can construct an array with |nRows| rows and |nCols| collumns and |nSlice| as third dimension
    // the leading dimension |ldim1| defaults to |nRows| and |ldim2| to |nCol|,
    // but can be set to any larger value to get the desired memory alignement.
    Array3d(size_type nRows,size_type nCols, size_type nSlices, size_type ldim1=0, size_type ldim2=0, const T& val=T()) {
      if(ldim1<nRows) 
	ldim1=nRows;
      if(ldim2<nCols)
        ldim2=nCols;
      if(bool(ldim1*ldim2*nSlices)) 
	{
	  nRow=nRows; nCol=nCols; nSlice=nSlices; lDim1=ldim1; lDim2=ldim2; lDim12=lDim1*lDim2;
	  owner=true;
	  data = new T[lDim12*nSlice];
	}
      else { 
	nRow=0; nCol=0; nSlice=0; lDim1=0; lDim2=0; lDim12=0; owner=true; data=0; 
      } 
      if(val==T(0))
	memset(data,0,sizeof(T)*lDim12*nSlice);
      else
	for(size_type i=0;i<lDim12*nSlice;i++) 
	  data[i]=val;
    }

    // Also a special constructor is provided to map a memory region to this matrix class. The user is responsible for the proper allocation and dealocation of memory and the proper dimensions.
    /* matrix that does not own data (user is responsible for proper dim. */
    Array3d(size_type nRows, size_type nCols, size_type nSlices, T* dat, size_type ldim1=0, size_type ldim2=0) 
      : data(dat) {
      if(ldim1<nRows) ldim1=nRows;
      if(ldim2<nCols) ldim2=nCols;
      if(bool(nRows*nCols*nSlices))
	{
	  nRow=nRows; nCol=nCols; nSlice=nSlices; lDim1=ldim1; lDim2=ldim2; lDim12=lDim1*lDim2;
	  owner=false;
	}
      else { nRow=0; nCol=0; nSlice=0; lDim1=0; lDim2=0; lDim12=0; owner=false; }
    }

// The final constructor makes a deep copy of another array.
    /* makes a deep copy of a */
    Array3d(const Array3d<T> & a) 
      : nRow(a.nRow), nCol(a.nCol), nSlice(a.nSlice), lDim1(a.nRow), lDim2(a.nCol),owner(true) {
      lDim12=lDim1*lDim2;
      if(bool(nRow*nCol*nSlice)) {
	owner=true;
	data = new T[lDim12*nSlice];
	memcpy(data,a.data,sizeof(T)*lDim12*nSlice);
      }
      else { nRow=0; nCol=0; nSlice=0; lDim1=0; lDim2=0; lDim12=0; owner=true; data=0; } 
    }

    // <Resizing an Array@>;
    void resize(size_type m,size_type n, size_type k, size_type ldim1=0, size_type ldim2=0) {
      if(!owner){std::logic_error("matrix not locally owned in T& Array3d<T>::resize");}
      if(ldim1<m) ldim1=m;
      if(ldim2<n) ldim2=n;
      if(ldim1*ldim2*k <= lDim1*lDim2*nSlice) {
        if(bool(ldim1*ldim2*k) && data) {
          nCol = n; nRow = m; nSlice=k; lDim1 = ldim1; lDim2 = ldim2; lDim12=lDim1*lDim2;
        }
        else {
          nRow=0; nCol=0; nSlice=0; lDim1=0; lDim2=0; lDim12=0;
        }
      }
      else {
        if(owner && data) delete [] data;
        nCol = n; nRow = m; nSlice = k; lDim1 = ldim1; lDim2 = ldim2; lDim12=lDim1*lDim2; owner=true;
        data = new T[lDim1*lDim2*nSlice];
      }
    }

    // <Array3d Destructor@>;
    // Finaly we need a destructor to free the memory allocated by the array.
    ~Array3d() {
      if(owner && data) delete [] data;
    }

    // <Matrix Access@>;
    // We provide two methods to access the lements of a matrix. The most natural one uses |operator()| with row and collumn arguments and looks simimar to a matrix access in {\bf Fortran}.
    inline T& operator() (size_type i, size_type j, size_type k) {
      // ASSERT(i<nRow,
      //        std::range_error("i>=n_row in T& Matrix<T>::operator()(size_type i,size_type j"));
      // ASSERT(j<nCol,
      //        std::range_error("j>=n_col in T& Matrix<T>::operator()(size_type i,size_type j"));
      return data[k*lDim12+j*lDim1+i];
    }

  // The second way to access the elements of a matrix uses the familair C |[]| operator and vies the |data| array as a continous block of memory.
    inline T& operator[](size_type i) {
      return data[i];
    }

  // We provide a few functions to return information about the matrix.
    size_type size() const { return lDim12*nSlice; }
    size_type n_row() const { return nRow; }
    size_type n_col() const { return nCol; }
    size_type n_slice() const {return nSlice;}
    size_type l_dim1() const { return lDim1; }
    size_type l_dim2() const { return lDim2;}
    
  // \subsection{Operations on 3 dim. arrays}

  // Assignments and copy:
    Array3d<T> &operator=(const Array3d<T>& a)
    {
      if (this == &a) return *this;   // Gracefully handle self assignment[12.1]
      if(owner && data) delete [] data;
      nRow=a.n_row(); nCol=a.n_col(); nSlice=a.n_slice(); lDim1=a.n_row(); lDim2=a.n_col(); lDim12=lDim1*lDim2; owner=true;
      if(bool(nRow*nCol*nSlice)) {
	data=new T[lDim12*nSlice];
	if(a.lDim12==lDim12)
	  memcpy(data,a.data,sizeof(T)*lDim12*nSlice);
	else
	  for(size_type k=0; k<nSlice; k++)
            for(size_type j=0;j<nCol;j++)
              for(size_type i=0;i<nRow;i++)
//					  data[j+i*lDim] = mat(i,j);
                data[i+j*lDim1+k*lDim12] = a.data[i+j*a.lDim1+k*lDim12];
      }
      else { nRow=0; nCol=0; nSlice=0; lDim1=0; lDim2=0; lDim12=0; owner=true; data=0; }
      return *this;
    }

    void copy(const Array3d<T>& A) {
      if(A.l_dim1()==lDim1 && A.l_dim2()==lDim2 && A.n_slice()==nSlice) {
        for(size_t i=0; i<lDim12*nSlice; i++) data[i]=A.data[i];
      }
      else if (A.n_row()==nRow && A.n_col()==nCol && A.n_slice()==nSlice){
        for(size_t k=0; k<nSlice; k++)        
          for(size_t j=0; j<nCol; j++)
            for(size_t i=0; i<nRow; i++)
//            data[j+lDim*i]=A(i,j);
              data[i+j*lDim1+k*lDim12] = A.data[i+j*A.l_dim1()+k*A.lDim12];
      }
      else std::range_error("matrix sizes don't match in Array3d<T>::copy");
    }

    inline
    ThisType& operator = (const T& val)
    {
      for(size_type i=0;i<lDim12*nSlice;i++) 
	data[i]=val;
      return *this;
    }

    // Scaling and adding:
    void scale(T a) {
      for(size_t i=0; i<lDim12*nSlice; i++) data[i]*=a;
    }

    void add(Array3d<T> &A) {
      if(A.l_dim1()==lDim1 && A.l_dim2()==lDim2 && A.n_slice()==nSlice) {
        for(size_t i=0; i<lDim12*nSlice; i++) data[i]+=A[i];
      }
      else if (A.n_row()==nRow && A.n_col()==nCol && A.n_slice()==nSlice){
        for(size_t k=0; k<nSlice; k++)
          for(size_t j=0; j<nCol; j++)
            for(size_t i=0; i<nRow; i++)
              data[i+j*lDim1+k*lDim12]+=A(i,j,k);
      }
      else std::range_error("matrix sizes don't match in Array3d<T>::add");
    }

    void addScaled(Array3d<T> &A, T a) {
      if(A.l_dim1()==lDim1 && A.l_dim2()==lDim2 && A.n_slice()==nSlice) {
        for(size_t i=0; i<lDim12*nSlice; i++) data[i]+=a*A[i];
      }
      else if (A.n_row()==nRow && A.n_col()==nCol && A.n_slice()==nSlice){
        for(size_t k=0; k<nSlice; k++)
          for(size_t j=0; j<nCol; j++)
            for(size_t i=0; i<nRow; i++)
              data[i+j*lDim1+k*lDim12]+=a*A(i,j,k);
      }
      else std::range_error("matrix sizes don't match in Array3d<T>::addScaled");
    }

   
  private:
    size_type nRow,nCol,nSlice,lDim1,lDim2,lDim12;
    bool owner;
    T* data;
    // T** col;
  };

template <class T>
inline  bool operator==(Array3d<T>& A, Array3d<T>& B)
{
  if(A.n_row()!=B.n_row() || A.n_col()!= B.n_col() || A.n_slice()!=B.n_slice())
    return false;
  for(size_t k=0; k<A.n_slice(); k++)
    for(size_t j=0; j<A.n_col(); j++)
      for(size_t i=0; i<A.n_row(); i++)
        if(A(i,j,k)!=B(i,j,k)) return false;
  return true;
}

template <class T>
inline  bool operator!=(Array3d<T>& A, Array3d<T>& B)
{
  if(A.n_row()!=B.n_row() || A.n_col()!= B.n_col() || A.n_slice()!=B.n_slice())
    return true;
  for(size_t k=0; k<A.n_slice(); k++)
    for(size_t j=0; j<A.n_col(); j++)
      for(size_t i=0; i<A.n_row(); i++)
        if(A(i,j,k)==B(i,j,k)) return false;
  return true;
}
	  
#endif
