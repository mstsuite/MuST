// -*- mode: c++; -*-

#ifndef BUILDKKRMATRIX_GPU_H
#define BUILDKKRMATRIX_GPU_H

#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "DeviceMatrix.hpp"
#include "DeviceArray3d.hpp"
#include "DeviceVector.hpp"

// #include "TestStructures.hpp"
#include "Misc/Indices.hpp"
#include "Main/SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"
#include "Misc/Coeficients.hpp"

//Instantiate this once, use the = operator from their host counterparts after the host has valid data. Then leave them on the device, pass into kernels as necessary
class DeviceConstants {
  public:
//  DeviceConstants() { }
//  ~DeviceConstants() { }

  DeviceVector<int> LIZStoreIdx;
  DeviceVector<int> LIZlmax;
  DeviceVector<Real> clm;
  DeviceVector<int> lofk;
  DeviceVector<int> mofk;
  DeviceVector<Complex> ilp1;
  DeviceMatrix<Real> LIZPos;
  DeviceMatrix<Complex> illp;
  DeviceArray3d<Real> cgnt;
};


void *allocateDConst(void);
void freeDConst(void * d_store);

void setupForBuildKKRMatrix_gpu(LSMSSystemParameters &lsms, AtomData &atom,DeviceConstants &d_const);
void setupForBuildKKRMatrix_gpu_opaque(LSMSSystemParameters &lsms, AtomData &atom, void *d_const);

void buildKKRMatrix_gpu(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, Complex energy, Complex prel, int iie, Matrix<Complex> &m, DeviceConstants &d_const);
void buildKKRMatrix_gpu_opaque(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, int ispin, Complex energy, Complex prel, int iie, Matrix<Complex> &m, void *d_const);

#endif
