#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "TestStructures.hpp"
#include "Misc/Coeficients.hpp"

// [eric] complex conversion
#define CXX2C(a) (*((double _Complex*)(&(a))))


extern "C"
void c_kkrmat(int, int, int, int, int, int, int*, int*, double*, int*, int*,
              double _Complex*, double _Complex*, double*, double*,
              double _Complex, double _Complex*, double _Complex*);


void buildKKRMatrix_nrel_ns2_acc(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, Complex energy, Complex prel, Matrix<Complex> &m)
{
  // make sure that lsms.n_spin_cant == 2!
  if(lsms.n_spin_cant!=2)
  {
    printf("lsms.n_spin_cant!=2 in buildKKRMatrix_nrel_ns2\n");
    exit(1);
  }

  c_kkrmat(lsms.maxlmax,lsms.angularMomentumIndices.ndlm,lsms.angularMomentumIndices.ndlj,
    2*atom.nrmat, atom.kkrsz, atom.numLIZ,
    &lsms.angularMomentumIndices.lofk[0], &lsms.angularMomentumIndices.mofk[0],
    &atom.LIZPos(0,0), &atom.LIZlmax[0], &atom.LIZStoreIdx[0],
    &CXX2C(iFactors.ilp1[0]), &CXX2C(iFactors.illp(0,0)),
    &gauntCoeficients.cgnt(0,0,0),
    &sphericalHarmonicsCoeficients.clm[0],
    CXX2C(prel), &CXX2C(m(0,0)), &CXX2C(local.tmatStore(0,0)));
}
