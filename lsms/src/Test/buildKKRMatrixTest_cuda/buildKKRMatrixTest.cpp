#include <stdlib.h>
#include <omp.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "TestStructures.hpp"
#include "Misc/Coeficients.hpp"

#include "LIZ_pos.h"

#include "buildKKRMatrix_gpu.hpp"
extern std::vector<DeviceConstants> deviceConstants;
//extern void setupForBuildKKRMatrix_gpu(LSMSSystemParameters &lsms, AtomData &atom,DeviceConstants &d_const);
//extern void buildKKRMatrix_gpu(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, Complex energy, Complex prel, int iie, Matrix<Complex> &m, DeviceConstants &d_const);
//#include <cuda_runtime.h>
//void clearM00(Complex *m, int blk_sz, int lda, cudaStream_t s);
#include "Accelerator/DeviceStorage.hpp"
extern DeviceStorage *deviceStorage;

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

void buildKKRMatrix_orig(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom,
                         Complex energy, Complex prel, Matrix<Complex> &m);
void buildKKRMatrix_nrel_ns2(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom,
                             Complex energy, Complex prel, Matrix<Complex> &m);

const double cphot=274.072;
const double c2inv=1.0/(cphot*cphot);

void set_tmat(Complex energy,Matrix<Complex> &tmat)
{
  energy=std::complex<Real>(-0.2,0.1);

  for(int i=0; i<32; i++)
    for(int j=0; j<32; j++)
      tmat(i,j)=0.0;

// spin up
// l=0
  tmat(0,0)=std::complex<Real>( 2.72989453e-01, 4.75691304e-02);
// l=1
  tmat(1,1)=std::complex<Real>(-1.53433276e-02, 9.50739941e-03);
  tmat(3,3)=tmat(2,2)=tmat(1,1);
// l=2
  tmat(4,4)=std::complex<Real>(-1.23189965e-02, 1.45766627e-02);
  tmat(8,8)=tmat(7,7)=tmat(6,6)=tmat(5,5)=tmat(4,4);
// l=3
  tmat(9,9)=std::complex<Real>( 1.53919816e-05,-1.07135429e-04);
  tmat(15,15)=tmat(14,14)=tmat(13,13)=tmat(12,12)=tmat(11,11)=tmat(10,10)=tmat(9,9);
// spin down
// l=0
  tmat(16,16)=std::complex<Real>( 4.51280708e-01, 2.16447659e-02);
// l=1
  tmat(17,17)=std::complex<Real>(-3.04540765e-02, 1.86224585e-02);
  tmat(19,19)=tmat(18,18)=tmat(17,17);
// l=2
  tmat(20,20)=std::complex<Real>(-1.00226432e-02, 1.24349804e-02);
  tmat(24,24)=tmat(23,23)=tmat(22,22)=tmat(21,21)=tmat(20,20);
// l=3
  tmat(25,25)=std::complex<Real>( 1.49731043e-05,-1.06213516e-04);
  tmat(31,31)=tmat(30,30)=tmat(29,29)=tmat(28,28)=tmat(27,27)=tmat(26,26)=tmat(25,25);
}

void set_atom(AtomData &atom)
{
  atom.kkrsz=(atom.lmax+1)*(atom.lmax+1);
  set_LIZ(atom);
  atom.nrmat=0;
  atom.LIZStoreIdx.resize(atom.numLIZ);
  atom.LIZlmax.resize(atom.numLIZ);
  for(int i=0; i<atom.numLIZ; i++)
  {
    atom.LIZStoreIdx[i]=0;
    atom.LIZlmax[i]=3;
    atom.nrmat+=(atom.LIZlmax[i]+1)*(atom.LIZlmax[i]+1);
  }
}

int main(int argc, char *argv[])
{
  Complex energy;
  Matrix<Complex> tmat(32,32);
  LSMSSystemParameters lsms;
  LocalTypeInfo local;
  AtomData atom;
  Matrix<Complex> m;

  atom.lmax=lsms.maxlmax=3;
  lsms.global.iprint=0;
  lsms.global.setIstop("main");
  lsms.nrel_rel=0;
  lsms.n_spin_cant=2;

  lsms.angularMomentumIndices.init(2*lsms.maxlmax);
  sphericalHarmonicsCoeficients.init(2*lsms.maxlmax);
  gauntCoeficients.init(lsms,lsms.angularMomentumIndices,sphericalHarmonicsCoeficients);
  iFactors.init(lsms,lsms.maxlmax);

  local.lDimTmatStore=32;
  local.tmatStore.resize(local.lDimTmatStore*local.lDimTmatStore,1);
  set_tmat(energy,tmat);
  for(int i=0; i<32*32; i++)
    local.tmatStore(i,0)=tmat[i];
  Complex prel=std::sqrt(energy*(1.0+energy*c2inv));

  set_atom(atom);

  m.resize(atom.nrmat*lsms.n_spin_cant,atom.nrmat*lsms.n_spin_cant);

  if(argc<2)
  {
    printf("Usage: buildKKRMatrixTest c|n|o [loop_count]\n");
    printf("  c - compare old and new construction results\n");
    printf("  n - time new matrix construction\n");
    printf("  o - time old matrix construction\n");
    exit(0); 
  }

  if(argc>1 && ((*argv[1]=='c')||(*argv[1]='n')))
  {
    #pragma omp parallel for default(none) shared(lsms,local,deviceConstants)
    for(int i=0; i<local.num_local; i++)
    {
      setupForBuildKKRMatrix_gpu(lsms,local.atom[i],deviceConstants[i]);
      // setupForBuildKKRMatrix_gpu_opaque(lsms,local.atom[i],deviceConstants[i]);
    }
    deviceStorage->allocate(maxkkrsz,lsms.n_spin_cant,maxNumLIZ,lsms.global.GPUThreads);
    copyTmatStoreToDevice(local);
  }
  double t0, t1;
  if(argc>1 && *argv[1]=='c')
  {
    Matrix<Complex> m1;
    printf("Comparing matrices\n");
    m1.resize(atom.nrmat*lsms.n_spin_cant,atom.nrmat*lsms.n_spin_cant);
    printf("Running new version ");
    t0 = omp_get_wtime();
    buildKKRMatrix_gpu(lsms, local, atom, energy, prel,m);
    t1 = omp_get_wtime();
    printf(" time: %lf seconds.\nRunning old version ", t1-t0);
    t0 = omp_get_wtime();
    buildKKRMatrix_orig(lsms, local, atom, energy, prel,m1);
    t1 = omp_get_wtime();
    printf(" time: %lf seconds.\n", t1-t0);
    for(int i=0; i<atom.nrmat*lsms.n_spin_cant; i++)
      for(int j=0; j<atom.nrmat*lsms.n_spin_cant; j++)
      {
        Real d=std::abs(m(i,j)-m1(i,j));
        Real rd=0.0;
        if(m(i,j)!=m1(i,j)) Real rd=d/std::abs(0.5*(m(i,j)+m1(i,j)));
        if(rd>1.0e-10) printf("m(i,j)!=m1(i,j): i=%d, j=%d, m(i,j)=(%lg,%lg), m1(i,j)=(%lg,%lg) d=%lg rd=%lg\n",
                                   i,j,std::real(m(i,j)),std::imag(m(i,j)),std::real(m1(i,j)),std::imag(m1(i,j)),
                                   d,rd);
      }
  } else if(argc>1 && *argv[1]=='n') {
    int loop_count=10;
    if(argc>2) loop_count=atoi(argv[2]);
    printf("Timing new matrix construction\n  loop_count=%d ",loop_count);
    t0 = omp_get_wtime();
    for(int i=0; i<loop_count; i++)
    {
      buildKKRMatrix_gpu(lsms, local, atom, energy, prel,m);
    }
    t1 = omp_get_wtime();
    printf(" time: %lf seconds.\n", t1-t0);
  } else {
// timing loop
    int loop_count=10;
    if(argc>2) loop_count=atoi(argv[2]);
    printf("Timing old matrix construction\n  loop_count=%d ",loop_count);
    t0 = omp_get_wtime();
    for(int i=0; i<loop_count; i++)
    {
      buildKKRMatrix_orig(lsms, local, atom, energy, prel,m);
    }
    t1 = omp_get_wtime();
    printf(" time: %lf seconds.\n", t1-t0);
  }

  return 0;
}
