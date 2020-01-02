#ifndef LSMS_REAL_H
#define LSMS_REAL_H

#include <cmath>

typedef double Real;

/*

#include <hdf5.h>
#include "TypeTraits.hpp"

#ifndef __CUDACC__
#include <mpi.h>


template<>
class TypeTraits<double>
{
//  static hid_t hdf5_type;
public:
  inline static MPI_Datatype mpiType(void) {return MPI_DOUBLE;}
//  inline static hid_t hdf5Type(void) { return hdf5_type;}
  inline static hid_t hdf5Type(void) { return H5T_NATIVE_DOUBLE;}
};
// hid_t TypeTraits<double>::hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);

template<>
class TypeTraits<int>
{
//   static hid_t hdf5_type;
public:
  inline static MPI_Datatype mpiType(void) {return MPI_INT;}
  inline static hid_t hdf5Type(void) { return H5T_NATIVE_INT;}
};
// hid_t TypeTraits<int>::hdf5_type = H5Tcopy(H5T_NATIVE_INT);

template<>
class TypeTraits<char>
{
//   static hid_t hdf5_type;
public:
  inline static MPI_Datatype mpiType(void) {return MPI_CHAR;}
  inline static hid_t hdf5Type(void) { return H5T_NATIVE_CHAR;}
};
// hid_t TypeTraits<char>::hdf5_type = H5Tcopy(H5T_NATIVE_CHAR);

#endif

*/
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510 
#endif
