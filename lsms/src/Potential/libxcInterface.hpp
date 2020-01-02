#ifndef LSMS_LIBXCINTERFACE_HPP
#define LSMS_LIBXCINTERFACE_HPP

#include "Main/SystemParameters.hpp"
#include "Real.hpp"

#ifdef USE_LIBXC
#include <xc.h>

// LibxcInterface is a singleton class to provide an interface to libxc
class LibxcInterface {
public:
  xc_func_type functional[numFunctionalIndices-1];
  int numFunctionals;
  bool needGradients; // the functional needs gradients of the density (for GGAs)
  bool needLaplacian; // need laplacians of the density (for MetaGGAs)
  bool needKineticEnergyDensity; // for MetaGGAs
  bool needExactExchange; // for Hybrid Functionals

  int init(int nSpin, int *xcFunctional);
  void evaluate(std::vector<Real> &rMesh, Matrix<Real> &rhoIn, int jmt, int nSpin, Matrix<Real> &xcEnergyOut, Matrix<Real> &xcPotOut);
  void evaluateSingle(Real *rhoIn, int nSpin, Real *xcEnergyOut, Real *xcPotOut);
};

#else // we don't link with libxc

class LibxcInterface {
public:
  int init(int nSpin, int *xcFunctional);
  void evaluate(std::vector<Real> &rMesh, Matrix<Real> &rhoIn, int jmt, int nSpin, std::vector<Real> &xcEnergyOut, Matrix<Real> &xcPotOut);
  void evaluateSingle(Real *rhoIn, int nSpin, Real *xcEnergyOut, Real *xcPotOut);
};
#endif
#endif
