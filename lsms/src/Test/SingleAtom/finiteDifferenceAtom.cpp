/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// Non Relativistic Renormalized Atom

#include <tuple>
#include <algorithm>

// typedef double Real;
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "../../Misc/rationalFit.hpp"
#include "../../Misc/integrateOneDim.cpp"

#include "calculateXC.cpp"

// solve the time independent Schroedinger equation
// using finite difference method on a log-mesh
//
// H^l f_i^l = E_i^l f_i^l
//
// where H^l = -d^2/dr^2 + V(r) + l(l+1)/r^2
//
// on a regular r mesh with spacing dr the form of H^l
// is a tridiagonal symmetric matrix of the form
// h_{i i} = - 2/dr^2 - V(r_i) - l(l+1)/r_i^2
// h_{i-1 i} = h_{i i-1} = 1/dr^2
//

void buildMatrix(std::vector<Real> &vr, int l, Real dr,
  std::vector<Real> &hDiagonal, std::vector<Real> &hSubdiagonal)
{
  for(int i=0; i<hDiagonal.size(); i++)
  {
    hDiagonal[i] = 2.0/(dr*dr) + vr[i]/(Real(i+1)*dr) + Real(l*(l+1))/((Real(i+1)*dr)*(Real(i+1)*dr));
  }
  for(int i=0; i<hSubdiagonal.size(); i++)
  {
    hSubdiagonal[i] = -1.0/(dr*dr);
  }
}

Real generateRadialMesh(std::vector<Real> &r_mesh, int N, Real r0, Real rN)
{
    if (N != r_mesh.size()) r_mesh.resize(N);

    Real dr = (rN-r0) / (N-1);
    for(int j=0; j<N; j++)
    {
      r_mesh[j] = (Real(j+1)*dr);
    }
    return dr;
}

int countNodes(Matrix<Real> &y, int c=0)
{
  int n=0;
  bool s = std::signbit(y(0,c));
  for(int i=1; i<y.n_row()-1; i++)
  {
    if(s != std::signbit(y(i,c)))
    {
      n++;
      s = std::signbit(y(i,c));
    }
  }
  return n;
}

void calculateRadialDensity(std::vector<Real> &r_mesh, Matrix<Real> &f, int solution, std::vector<Real> &rho)
{
  std::vector<Real> rhoIntegrated(r_mesh.size());

  for(int i=0; i<r_mesh.size(); i++)
  {
    rho[i] = std::abs(f(i,solution))*std::abs(f(i,solution));
  }

  // integrateOneDimSpherical(r_mesh, rho, rhoIntegrated);
  // the p already contains a factor of r,
  // thus the rho above is actually r^2 rho and a standard integration
  // is to be performed to yield the spherical integral of |psi|^2
  // we only need a facto of 4 pi
  integrateOneDim(r_mesh, rho, rhoIntegrated);

  Real normalization = 1.0/(4.0 * M_PI * rhoIntegrated[rhoIntegrated.size()-1]);

  for(int i=0; i<r_mesh.size(); i++)
  {
    rho[i] *= normalization;
  }
}

extern "C" {
void dstevx_(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, double *vu,
             int *il, int *iu, double *abstol,
             int *m, double *w, double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);
}

void findSchroedingerEigenvalues(std::vector<Real> &r_mesh, std::vector<Real> &vr, int l, int numSolutions,
                                 std::vector<Real> &energy, Matrix<Real> &f)
{
  std::vector<Real> hDiagonal(r_mesh.size()), hSubdiagonal(r_mesh.size());
  buildMatrix(vr, l, r_mesh[1]-r_mesh[0], hDiagonal, hSubdiagonal);

// finde eigenvalues and eigenvectors using LAPACK dstevx:
// SUBROUTINE DSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL,
//                          M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
//
//       .. Scalar Arguments ..
//       CHARACTER          JOBZ, RANGE
//       INTEGER            IL, INFO, IU, LDZ, M, N
//       DOUBLE PRECISION   ABSTOL, VL, VU
//       ..
//       .. Array Arguments ..
//       INTEGER            IFAIL( * ), IWORK( * )
//       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )

  // jobz = 'V' : find eigenvalues and eigenvectors
  // range = 'I' : find IL through IU eigenvalues
  // N = dim of Hamiltonian
  // D[N] = diagonal elements, might multiplied by a factor on exit
  // E[N-1] the subdiagonal elements, might multiplied by a factor on exit
  // VU, VL used only if range = 'V'
  // IL, IU idex of the smalest and largest eigenvalues, 1 <= IL <= IU <= N, used only if range = 'I'
  // ABSTOL = 2*DLAMCH('S')
  // out:
  // M = The total number of eigenvalues found.  0 <= M <= N. if RANGE = 'I', M = IU-IL+1
  // W = the eigenvalues
  // Z(LDZ, M) = eigenvectors
  // WORK[5*N], IWORK[5*N]
  // INFO = 0 success, <0 illegal value in Argument -INFO, >0 INFO eigenvectors didn't converge
  // IFAIL[N] indices of eigenvalues that didn't converge

  int n=hDiagonal.size();
  Real zero=0.0;
  Real abstol = 0.0;
  int one=1, numSolutionsFound, info;
  std::vector<Real> work(5*n);
  std::vector<int> iwork(5*n);
  std::vector<int> ifail(n);
  dstevx_("V", "I", &n, &hDiagonal[0], &hSubdiagonal[0], &zero, &zero, &one, &numSolutions, &abstol,
          &numSolutionsFound, &energy[0], &f(0,0), &n, &work[0], &iwork[0], &ifail[0], &info);
}

class AtomOrbital {
public:
  int n,l;
  Real energy;
  std::vector<Real> rho; // radial charge density contribution from an electron in this level
  // char type; // 'C' for core lectron, 'S' for semi-core and 'V' for valence electrons
};

void accumulateDensities(std::vector<AtomOrbital> &orbitalEnergiesAndDensities, int atomicNumber,
                         std::vector<Real> &rhotot)
{
  for(int i=0; i<rhotot.size(); i++) rhotot[i]=0;
  int electronsMissing = atomicNumber; // still need Z electrons
  int orbitalIdx=0; // the index of the current orbital
  while(electronsMissing>0)
  {
    if(electronsMissing >= 2*(2*std::abs(orbitalEnergiesAndDensities[orbitalIdx].l)+1)) // a filled orbital
    {
      for(int i=0; i<rhotot.size(); i++)
	rhotot[i] += orbitalEnergiesAndDensities[orbitalIdx].rho[i]
	  * 2 * (2 * std::abs(orbitalEnergiesAndDensities[orbitalIdx].l) + 1);
      electronsMissing -= 2 * (2 * std::abs(orbitalEnergiesAndDensities[orbitalIdx].l) + 1);
      printf("%d %2d %lg Ry: filled (%2d electrons)\n",
	     orbitalEnergiesAndDensities[orbitalIdx].n,
	     orbitalEnergiesAndDensities[orbitalIdx].l,
	     orbitalEnergiesAndDensities[orbitalIdx].energy,
	     2 * (2 * std::abs(orbitalEnergiesAndDensities[orbitalIdx].l) + 1));
    } else {
      for(int i=0; i<rhotot.size(); i++)
	rhotot[i] += orbitalEnergiesAndDensities[orbitalIdx].rho[i]
	  * electronsMissing;
      
      printf("%d %2d %lg Ry: partially filled (%2d electrons)\n",
	     orbitalEnergiesAndDensities[orbitalIdx].n,
	     orbitalEnergiesAndDensities[orbitalIdx].l,
	     orbitalEnergiesAndDensities[orbitalIdx].energy,
	     electronsMissing);

      electronsMissing = 0;
    }
    orbitalIdx++;
  }
}

// generate lmax for fixed atom
void lmaxAndNMaxFromZ(int Z, int &lmax, int &nMax)
{
  // He: 1s^2
  lmax = 0; nMax = 1;
  if(Z<3) return;
  // Ne: [He] 2s^2 2p^6
  lmax = 1; nMax = 2;
  if(Z<11) return;
  // Ar: [Ne] 3s^2 3p^6
  nMax = 3;
  if(Z<19) return;
  // Kr: [Ar] 3d^10 4s^2 4p^6
  lmax = 2; nMax =4;
  if(Z<37) return;
  // Xe: [Kr] 4d^10 5s^2 5p^6
  nMax = 5;
  if(Z<55) return;
  // Rn: [Xe] 4f^14 5d^10 6s^2 6p^6
  lmax = 3; nMax = 6;
  if(Z<87) return;
  // Og: [Rn] 5f^14 6d^10 7s^2 7p^6
  nMax = 7;
}

void calculateOrbitals(std::vector<Real> &r_mesh, std::vector<Real> &vr, int atomicNumber,
                       std::vector<AtomOrbital> &orbitalEnergiesAndDensities)
{
  int lmax, nMax;
  int orbitalIdx =0;
  lmaxAndNMaxFromZ(atomicNumber, lmax, nMax);
  int totalNumberOrbitals = 0;

  for(int l=0; l<=lmax; l++)
    totalNumberOrbitals +=  nMax-l;
  if(orbitalEnergiesAndDensities.size() != totalNumberOrbitals)
    orbitalEnergiesAndDensities.resize(totalNumberOrbitals);
  
  printf(" n  l  energy\n");
  for(int l=0; l<=lmax; l++)
  {
    std::vector<Real> energies(20);
    Matrix<Real> f(r_mesh.size(),20);

    int numStates = nMax-l;
    int principalQuantumNumber;
  
    findSchroedingerEigenvalues(r_mesh, vr, l, numStates, energies, f);

    for(int i=0; i<numStates; i++)
    {
      principalQuantumNumber = i + l + 1;
      orbitalEnergiesAndDensities[orbitalIdx].n = principalQuantumNumber;
      orbitalEnergiesAndDensities[orbitalIdx].l = l;
      orbitalEnergiesAndDensities[orbitalIdx].energy = energies[i];

      printf(" %d %2d %lg Ry",principalQuantumNumber, l, orbitalEnergiesAndDensities[orbitalIdx].energy);

      orbitalEnergiesAndDensities[orbitalIdx].rho.resize(r_mesh.size());
      calculateRadialDensity(r_mesh, f, i, orbitalEnergiesAndDensities[orbitalIdx].rho);

      if(orbitalEnergiesAndDensities[orbitalIdx].rho[r_mesh.size()-1] > 0.0001)
      {
        printf(" !!\n");
      } else {
        printf("\n");
      }
      orbitalIdx++;
    }
  }
  
  std::sort(orbitalEnergiesAndDensities.begin(), orbitalEnergiesAndDensities.end(),
	    [](AtomOrbital const & a, AtomOrbital const &b){return a.energy < b.energy;});

}

// assume that rho is actually r^2 rho, calulate r v
void sphericalPoisson(std::vector<Real> &rho, std::vector<Real> &r_mesh, std::vector<Real> &vr)
{
  // V_l(r) = 4pi/(2l+1) [1/r^(l+1) int_0^r rho_l(r') r'^(l+2) dr'
  //                      + r^l     int_r^R rho_l(r') / r'^(l-1) dr']
  // l = 0
  // vr <- int_0^r rho(r') r'^2 dr'
  integrateOneDim(r_mesh, rho, vr);
  // vr <- 1/r  int_0^r rho(r') r'^2 dr' THIS IS NOT NEEDED AS WE STORE V*r
  // for(int i=0; i<r_mesh.size(); i++)
  //  vr[i] = vr[i] / r_mesh[i];

  //  int_r^R rho(r') / r^(-1) dr'
  //  = int_0^R rho(r') r' dr' - int_0^r rho(r') r' dr'
  std::vector<Real> integral;
  integral.resize(r_mesh.size());
  //
  
  integrateOneDimRPower(r_mesh, rho, integral, -1); // p=-1 because rho is actually rho * r^2
  for(int i=0; i<r_mesh.size(); i++)
    vr[i] = 2.0*4.0*M_PI*(vr[i] + (integral[integral.size()-1] - integral[i]) * r_mesh[i]); 
}

void calculateSingeOrbital(int atomicNumber, int l, int principleQuantumNumber, Real atomRadius)
{
  std::vector<Real> r_mesh, vr, rho;
  generateRadialMesh(r_mesh, 6000, 1.5e-10, atomRadius);
  vr.resize(r_mesh.size());
  rho.resize(r_mesh.size());
  std::vector<Real> energies(20);
  Matrix<Real> f(r_mesh.size(),20);

  for(int i=0; i<r_mesh.size(); i++) vr[i] = -2.0*Real(atomicNumber);
  
  findSchroedingerEigenvalues(r_mesh, vr, l, principleQuantumNumber-l, energies, f);
  calculateRadialDensity(r_mesh, f, principleQuantumNumber-l-1, rho);

  printf("number of nodes: %d\n",countNodes(f,principleQuantumNumber-l-1));
  
  FILE *outf=fopen("orbital.out","w");
  fprintf(outf,"# Atomic Number: %d\n",atomicNumber);
  fprintf(outf,"# energy: %lg\n", energies[principleQuantumNumber-l-1]);
  fprintf(outf,"# i r[i] rho\n");
  Real exc;
  for(int i=0; i<r_mesh.size(); i++)
  {
    fprintf(outf,"%5d %lg %lg\n",i, r_mesh[i], rho[i]);
  }
  fclose(outf);
}

void printUsage(const char *name, Real R)
{
  printf("Usage: %s Z [R]\n",name);
  printf("       Z: atomic number\n");
  printf("       R: atomic sphere radius (optional, default=%lf)\n",R);
  printf("       %s o Z l principleQuantumNumber [R]\n",name);
}

int main(int argc, char *argv[])
{
  int atomicNumber = 29; // test copper
  int principalQuantumNumber = 1;
  int l=0;
  Real atomRadius = 3.0;
  Real rmsRho;
  char mode = 's'; // scf atom
  Real energy;

  if(argc != 2 && argc != 3 && argc != 5 && argc != 6)
  {
    printUsage(argv[0], atomRadius);
    exit(1);
  }
  if(argv[1][0]=='o')
  {
    mode = 'o'; // one orbital
    atomicNumber = atoi(argv[2]);
    l = atoi(argv[3]);
    principalQuantumNumber = atoi(argv[4]);
    if (argc == 6)
      atomRadius = atof(argv[2]);
    calculateSingeOrbital(atomicNumber, l, principalQuantumNumber, atomRadius);
    return 0;
  } else {
    atomicNumber = atoi(argv[1]);
    if(argc == 3)
      atomRadius = atof(argv[2]);
  }

  int maxPrincipalQuantumNumber;
  
  std::vector<Real> r_mesh, vr, vXC;
  // std::vector<Matrix<Real> > pqs;
  std::vector<AtomOrbital> orbitalEnergiesAndDensities;
  // pqs.resize(orbitals.size());
  orbitalEnergiesAndDensities.resize(50);
  
// initialize r_mesh (unit of length is the Bohr radius)
  generateRadialMesh(r_mesh, 1500, 1.5e-10, atomRadius);
  std::vector<Real> rhotot, rhonew;
  rhotot.resize(r_mesh.size());
  rhonew.resize(r_mesh.size());
  
// initialize Z/r potential e^2=2
  vr.resize(r_mesh.size());
  vXC.resize(r_mesh.size());
  for(int i=0; i<r_mesh.size(); i++) vr[i] = -2.0*Real(atomicNumber);

  calculateOrbitals(r_mesh, vr, atomicNumber, orbitalEnergiesAndDensities);
  // accumulate the charge from all the electrons in the atom
  accumulateDensities(orbitalEnergiesAndDensities, atomicNumber, rhotot);
  sphericalPoisson(rhotot, r_mesh, vr);
  // add nuclear charge
  for(int i=0; i<r_mesh.size(); i++) vr[i] += -2.0*Real(atomicNumber);
  
  // iterate on the charge density with mixing
  Real mixing = 0.05;
  Real rms=1.0;
  for(int iter=0; iter < 50 && rms > 0.0001; iter++)
  {
    calculateOrbitals(r_mesh, vr, atomicNumber, orbitalEnergiesAndDensities);
    // accumulate the charge from all the electrons in the atom
    accumulateDensities(orbitalEnergiesAndDensities, atomicNumber, rhonew);

    rms=0.0;
    for(int i=0; i<rhotot.size(); i++)
      rms += (rhotot[i]-rhonew[i])*(rhotot[i]-rhonew[i]);
    rms = std::sqrt(rms/Real(rhotot.size()));
    printf("iter %3d: rms = %lg\n",iter, rms);

    // mixing
    for(int i=0; i<rhotot.size(); i++)
      rhotot[i] = (1.0-mixing)*rhotot[i] + mixing*rhonew[i];

    sphericalPoisson(rhotot, r_mesh, vr);
    // calculate exchange-correltaion potential
    exchangeCorrelationPotentialLDA(rhotot, r_mesh, vXC);
    // add nuclear charge and exchange-correlation potential
    for(int i=0; i<r_mesh.size(); i++) vr[i] += -2.0*Real(atomicNumber) + vXC[i]*r_mesh[i];
  }
    
  ///*
  FILE *outf=fopen("rho_vr.out","w");
  fprintf(outf,"# Atomic Number: %d\n",atomicNumber);
  fprintf(outf,"# i r[i] rho vr vXC\n");
  Real exc;
  for(int i=0; i<r_mesh.size(); i++)
  {
    fprintf(outf,"%5d %lg %lg %lg %lg %lg\n",i, r_mesh[i], rhotot[i], vr[i], vXC[i],
            chachiyo2016(rhotot[i]/(4.0*M_PI*r_mesh[i]*r_mesh[i]), exc));
  }
  fclose(outf);
  //*/

  return 0;
}
