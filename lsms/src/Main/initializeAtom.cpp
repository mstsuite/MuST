/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <vector>
#include <algorithm>
#include <cstdio>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "SingleSite/AtomData.hpp"
#include "initializeAtom.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/integrateOneDim.cpp"
#include "Potential/calculateChargesPotential.hpp"

/* initializeAtom(AtomData &a)
   initialize an atom using the following information in AtomData

    Mesh information: xstart, rmt, jmt, jws
    Atom information: ztotss, zcorss, zsemss, zvalss
    (Note that ztotss=zcorss+zsemss+zvalss)
    lmax to limit core search
*/
class InitialAtomLevel {
public:
  Real energy;
  int n,l;
  std::vector<Real> rho; // radial charge density contribution from an electron in this level
  char type; // 'C' for core lectron, 'S' for semi-core and 'V' for valence electrons
} ;

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

Real generateLinearRadialMesh(std::vector<Real> &r_mesh, int N, Real r0, Real rN)
{
    if (N != r_mesh.size()) r_mesh.resize(N);

    Real dr = (rN-r0) / (N-1);
    for(int j=0; j<N; j++)
    {
      r_mesh[j] = (Real(j+1)*dr);
    }
    return dr;
}

void calculateRadialDensity(std::vector<Real> &r_mesh, Matrix<Real> &f, int solution, Real energy, std::vector<Real> &rho)
{
  const bool constantDensityForPositiveEnergies=false;
  std::vector<Real> rhoIntegrated(r_mesh.size());

  for(int i=0; i<r_mesh.size(); i++)
  {
    rho[i] = std::abs(f(i,solution))*std::abs(f(i,solution));
  }
  if(constantDensityForPositiveEnergies && energy > 0.0)
  {
    for(int i=0; i<r_mesh.size(); i++)
    {
      rho[i] = r_mesh[i]*r_mesh[i];
    } 
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
                       std::vector<InitialAtomLevel> &orbitalEnergiesAndDensities)
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
      calculateRadialDensity(r_mesh, f, i, energies[i], orbitalEnergiesAndDensities[orbitalIdx].rho);

      /*
      if(orbitalEnergiesAndDensities[orbitalIdx].rho[r_mesh.size()-1] > 0.0001)
      {
        printf(" !!\n");
      } else {
        printf("\n");
      }
      */
      printf("\n");

      orbitalIdx++;
    }
  }
  
  std::sort(orbitalEnergiesAndDensities.begin(), orbitalEnergiesAndDensities.end(),
	    [](InitialAtomLevel const & a, InitialAtomLevel const &b){return a.energy < b.energy;});

}

void accumulateDensities(std::vector<InitialAtomLevel> &orbitalEnergiesAndDensities, int atomicNumber,
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

// correlation energy from T. Chachiyo, J. Chem. Phys. 145, 021101 (2016).
Real chachiyo2016(Real rho, Real &exc)
{
  // total XC energy = int rho*exc dr
  // exchange energy density:
  //
  //Real ex = 0;
  Real ex = -6.0*std::pow(rho * 3.0/(64.0*M_PI), 1.0/3.0);
  // correlation energy density:
  // e_c = a ln(1 + b/r_s + b/r_s^2)
  // constatnts from eq. 3 in T. Chachiyo, J. Chem. Phys. 145, 021101 (2016)
  Real const a = (std::log(2) - 1.0)/(M_PI * M_PI); // converted from Hartree in Chachiyo to Rydberg
  Real const b = 20.4562557;
  // r_s = (4 pi rho / 3) ^ -1/3 ->
  // 1/rs = (4 pi rho / 3) ^ 1/3
  Real rsInv = std::pow(4.0 * M_PI * rho / 3.0, 1.0/3.0);
  Real ec = a*std::log(1.0 + b*rsInv + b*rsInv*rsInv);
  exc = ex + ec;
  
  // exchange correlation potential:
  // v_xc(r) = e_xc(rho(r)) + rho(r) * d e_xc(rho)/d rho
  // rho_dex = rho * d e_x / d rho = 1/3 e_x
  Real rho_dex = ex/3.0;
  //
  Real rho_dec = a*b*(rsInv + 2.0*rsInv*rsInv)/(3.0*(1.0 + b*rsInv + b*rsInv*rsInv));

  return exc + rho_dex + rho_dec;
}
    
Real exchangeCorrelationPotentialLDA(Real rho, Real r)
{
  Real eXC;
  /*
  if(rho < 1.0e-9)
    return 0.0;
  // return alpha2(std::pow(3.0/rho , 1.0/3.0), 0.0, 1.0, eXC);
  return alpha2(std::pow(3.0*r*r/rho , 1.0/3.0), 0.0, 1.0, eXC);
  //*/
  return chachiyo2016(rho/(4.0*M_PI*r*r), eXC);
}

void exchangeCorrelationPotentialLDA(std::vector<Real> &rho, std::vector<Real> &r_mesh, std::vector<Real> &vXC)
{
  for(int i=0; i<rho.size(); i++)
  {
    vXC[i] = exchangeCorrelationPotentialLDA(rho[i], r_mesh[i]);
  }
}


void iterateOrbitals(int atomicNumber, Real atomRadius, std::vector<Real> &r_mesh,
                     std::vector<InitialAtomLevel> &orbitalEnergiesAndDensities,
                     std::vector<Real> &vr)
{
  int linearMeshSize = 6000;
  Real meshOrigin = 1.5e-10;
  int maxIterations = 100;
  Real mixingDensity = 1.0;
  Real mixingPot = 0.1;
  Real rmsTarget = 0.0001;
  
  std::vector<Real> vXC, vrNew;
  // std::vector<InitialAtomLevel> orbitalEnergiesAndDensities;
  
// initialize r_mesh (unit of length is the Bohr radius)
  generateLinearRadialMesh(r_mesh, linearMeshSize, meshOrigin, atomRadius);
  std::vector<Real> rhotot, rhonew;
  rhotot.resize(r_mesh.size());
  rhonew.resize(r_mesh.size());
  
// initialize Z/r potential e^2=2
  vr.resize(r_mesh.size());
  vrNew.resize(r_mesh.size());
  vXC.resize(r_mesh.size());
  for(int i=0; i<r_mesh.size(); i++) vr[i] = -2.0*Real(atomicNumber);

  calculateOrbitals(r_mesh, vr, atomicNumber, orbitalEnergiesAndDensities);
  // accumulate the charge from all the electrons in the atom
  accumulateDensities(orbitalEnergiesAndDensities, atomicNumber, rhotot);
  sphericalPoisson(rhotot, r_mesh, vr);
  // add nuclear charge
  for(int i=0; i<r_mesh.size(); i++) vr[i] += -2.0*Real(atomicNumber);
  
  // iterate on the charge density with mixing
  
  Real rms=1.0;
  for(int iter=0; iter < maxIterations && rms > rmsTarget; iter++)
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
      rhotot[i] = (1.0-mixingDensity)*rhotot[i] + mixingDensity*rhonew[i];

    sphericalPoisson(rhotot, r_mesh, vrNew);
    // calculate exchange-correltaion potential
    exchangeCorrelationPotentialLDA(rhotot, r_mesh, vXC);
    // add nuclear charge and exchange-correlation potential
    for(int i=0; i<r_mesh.size(); i++) vrNew[i] += -2.0*Real(atomicNumber) + vXC[i]*r_mesh[i];

    // mixing
    for(int i=0; i<vr.size(); i++)
      vr[i] = (1.0-mixingPot)*vr[i] + mixingPot*vrNew[i];
  }
}

int calculateAtomLevels(AtomData &a, std::vector<InitialAtomLevel> &atomLevels,
			Matrix<Real> &rhoValence, bool writeOrbitals=false)
{
  // find the lowest zcorss+zsemss levels
  // maximum number of states to test: (lmax+1)*(lmax+2)/2
  // assuming kappa degeneracy, iterating over principal quantum number and l
  // we treat core and semi-core states the same for the initialization,
  // using finite difference atom

  // iterate the finite difference atom
  std::vector<Real> linear_mesh;
  std::vector<Real> vrLinearMesh;
  iterateOrbitals(a.ztotss, a.r_mesh[a.r_mesh.size()-1], linear_mesh, atomLevels, vrLinearMesh);

  char fname[256]; 
  FILE *orbitalFile;
  if(writeOrbitals)
  {
    for(int nl = 0; nl<atomLevels.size(); nl++)
    {
      sprintf(fname,"orbital_density_n%d_l%d",atomLevels[nl].n,atomLevels[nl].l);
      orbitalFile=fopen(fname,"w");

      for(int ir=0; ir<linear_mesh.size(); ir++)
      {
        fprintf(orbitalFile,"%4d %g %g\n",ir,linear_mesh[ir],atomLevels[nl].rho[ir]);
      }
      fclose(orbitalFile);
    }
  }
      
// fill core states:
  int coreTarget = a.zcorss + a.zsemss;
  int coreElectrons = 0;
  a.numc = 0;
  for (int i=0; i<atomLevels.size(); i++)
  {
    if (coreElectrons < coreTarget)
    {
      if(coreElectrons<a.zcorss)
      {
        printf("C ");
        atomLevels[i].type='C';
      } else {
        printf("S ");
        atomLevels[i].type='S';
      }
      coreElectrons += 2*(2*atomLevels[i].l+1);
      if(atomLevels[i].l == 0)
        a.numc += 1;
      else
        a.numc += 2;
    }
    else
    {
      printf("V ");
      atomLevels[i].type='V';
    }
    printf("n=%2d  l=%2d  :  Energy=%12.6lf\n",atomLevels[i].n,atomLevels[i].l,atomLevels[i].energy);
  }

  a.resizeCore(a.numc);
  int j = 0;
  for (int i=0; i<a.numc; i++)
  {
    a.ec(i,0) = a.ec(i,1) = atomLevels[j].energy;
    a.nc(i,0) = a.nc(i,1) = atomLevels[j].n;
    a.lc(i,0) = a.lc(i,1) = atomLevels[j].l;
    a.kc(i,0) = a.kc(i,1) = -atomLevels[j].l-1;
    if(atomLevels[j].type=='C') // accumulate the density for core levels
    {
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
        a.corden(ir,0) += Real((3-a.nspin)*(-a.kc(i,0))) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
        a.corden(ir,1) += Real((3-a.nspin)*(-a.kc(i,1))) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
      }
    } else if(atomLevels[j].type=='S') { // accumulate the density for semi-core levels
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
        a.semcor(ir,0) += Real((3-a.nspin)*(-a.kc(i,0))) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
        a.semcor(ir,1) += Real((3-a.nspin)*(-a.kc(i,1))) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
      }
    }
    if(atomLevels[j].l != 0)
    {
      i++;
      a.ec(i,0) = a.ec(i,1) = atomLevels[j].energy;
      a.nc(i,0) = a.nc(i,1) = atomLevels[j].n;
      a.lc(i,0) = a.lc(i,1) = atomLevels[j].l;
      a.kc(i,0) = a.kc(i,1) = atomLevels[j].l;
      if(atomLevels[j].type=='C') // accumulate the density for core levels
      {
        for(int ir=0; ir<a.r_mesh.size(); ir++)
        {
          a.corden(ir,0)+=Real((3-a.nspin)*a.kc(i,0)) *
            interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
          a.corden(ir,1)+=Real((3-a.nspin)*a.kc(i,1)) *
            interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
        }
      } else if(atomLevels[j].type=='S') { // accumulate the density for semi-core levels
        for(int ir=0; ir<a.r_mesh.size(); ir++)
        {
          a.semcor(ir,0)+=Real((3-a.nspin)*a.kc(i,0)) *
            interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
          a.semcor(ir,1)+=Real((3-a.nspin)*a.kc(i,1)) *
            interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
        }
      }
    }
    j++;
    // printf("%d %d\n",i,j);
  }
  // accumulate valence density
  for(int ir=0; ir<a.r_mesh.size(); ir++)
  {
    rhoValence(ir,0)=rhoValence(ir,1)=0.0;
  }
  int numValenceRemaining = a.zvalss;
  j=0; while(atomLevels[j].type != 'V') j++;
  while(numValenceRemaining>0)
  {
    int multiplicity = 2*atomLevels[j].l+1;
    if(2*multiplicity <= numValenceRemaining)
    {
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
	rhoValence(ir,0) += Real(multiplicity*(3-a.nspin)) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
	rhoValence(ir,1) += Real(multiplicity*(3-a.nspin)) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
      }
      numValenceRemaining -= 2*multiplicity;
      j++;
    } else {
      for(int ir=0; ir<a.r_mesh.size(); ir++)
      {
	rhoValence(ir,0) += 0.5*Real(numValenceRemaining*(3-a.nspin)) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
	rhoValence(ir,1) += 0.5*Real(numValenceRemaining*(3-a.nspin)) *
          interpolate(linear_mesh, atomLevels[j].rho, a.r_mesh[ir]);
      }
      numValenceRemaining = 0;
      j++;
    }
  }

  for(int ir=0; ir<a.r_mesh.size(); ir++)
  {
    a.vr(ir,0) = a.vrNew(ir,0) = a.vr(ir,1) = a.vrNew(ir,1) =
      interpolate(linear_mesh, vrLinearMesh, a.r_mesh[ir]);
  }
  
  return coreElectrons;
}

void initializeAtom(AtomData &a)
{
  a.generateRadialMesh();

  // for analysis purposes gnerate gnuplot files for the atom
  bool generatePlot=true;
  FILE * plotFile;
  if(generatePlot) plotFile=fopen("initializeAtom.plot","w");
  
  // inititalize potential to be V(r) = -2Z/r
  // add a potential corresponding to a gaussian charge distribution
  // for the core charge with \sigma=0.2*rmt 
  // (vr stores V(r) * r)
  Real sigmaSqrt2Inv = 1.0 / (0.2 * std::sqrt(2.0) * a.rmt);
  Real q = a.zcorss + a.zsemss;
  for (int ir=0; ir<a.r_mesh.size(); ir++)
  {
    a.vr(ir,0) = a.vr(ir,1) = -2.0*a.ztotss; // +q*approxErfc(a.r_mesh[ir]*sigmaSqrt2Inv);
    // clear core and semi-core densities
    a.corden(ir,0) = a.corden(ir,1) = a.semcor(ir,0) = a.semcor(ir,1) = 0.0;
    a.rhotot(ir,0) = a.rhotot(ir,1) = a.rhoNew(ir,0) = a.rhoNew(ir,1) = 0.0;
  }
  // find the lowest zcorss+zsemss levels
  // maximum number of states to test: (lmax+1)*(lmax+2)/2
  
  // assuming kappa degeneracy, iterating over principal quantum number and l
  // we treat core and semi-core states the same for the initialization


  std::vector<InitialAtomLevel> atomLevels(1);
  Matrix<Real> rhoValence(a.r_mesh.size(),2);
  
  int coreTarget = a.zcorss + a.zsemss;
  int coreElectrons = calculateAtomLevels(a, atomLevels, rhoValence, true);
  
  if (coreElectrons != coreTarget)
  {
    printf("Warning: initializeAtom can't satisfy the core electron requirement:\n  Target: %d (%lf + %lf)\n  Actual: %d\n",
           coreTarget, a.zcorss, a.zsemss, coreElectrons);
  }

  Real atomVolume=(4.0/3.0)*M_PI*a.rws*a.rws*a.rws;  // a.rmt*a.rmt*a.rmt;
 
  for(int ir=0; ir<a.r_mesh.size(); ir++)
  {
    a.rhotot(ir,0)=a.rhoNew(ir,0)=a.corden(ir,0)+a.semcor(ir,0)+rhoValence(ir,0);
    a.rhotot(ir,1)=a.rhoNew(ir,1)=a.corden(ir,1)+a.semcor(ir,1)+rhoValence(ir,1);
  }
  
  if(generatePlot)
  {
    fprintf(plotFile,"set term pdf\nset outp 'initialPotential.pdf'\n");
    fprintf(plotFile,"set xrange [0:%lf]\n",a.r_mesh[a.r_mesh.size()-1]);

    fprintf(plotFile,"plot '-' with lines title 'Vr spin up'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.vr(ir,0));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'Vr spin down'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.vr(ir,1));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'rho spin up'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.rhotot(ir,0));
    }
    fprintf(plotFile,"e\n");
    
    fprintf(plotFile,"plot '-' with lines title 'rho spin down'\n");
    for(int ir=0; ir<a.r_mesh.size(); ir++)
    {
      fprintf(plotFile,"%18.12lf   %18.12lf\n",a.r_mesh[ir],a.rhotot(ir,1));
    }
    fprintf(plotFile,"e\n");
  }
  
  
  if(generatePlot) fclose(plotFile);
  // exit(1);
}


int initializeNewPotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local)
{
  for (int i=0; i<local.num_local; i++)
  {
    printf("Initializing potential %d.\n", i);
    snprintf(local.atom[i].header, 40, "Potential Initialized by LSMS_3");
    for (int j=31; j<80; j++)
      local.atom[i].header[j] = ' ';
    local.atom[i].resizePotential(lsms.global.iprpts);
    local.atom[i].ztotss = (Real)crystal.types[local.global_id[i]].Z;
    local.atom[i].zcorss = (Real)crystal.types[local.global_id[i]].Zc;
    local.atom[i].zsemss = (Real)crystal.types[local.global_id[i]].Zs;
    local.atom[i].zvalss = (Real)crystal.types[local.global_id[i]].Zv;
    local.atom[i].forceZeroMoment = crystal.types[local.global_id[i]].forceZeroMoment;
    local.atom[i].vdif = 0.0;
    local.atom[i].vdifNew = 0.0;
    local.atom[i].spinFlipped = false;
    local.atom[i].xvalws[0] = local.atom[i].xvalws[1] = 0.5 * local.atom[i].zvalss;
    local.atom[i].lmax = crystal.types[local.global_id[i]].lmax;
    local.atom[i].kkrsz = (local.atom[i].lmax + 1) * (local.atom[i].lmax + 1);

    local.atom[i].xstart = -11.1309674;
    local.atom[i].jmt = 1001;
    local.atom[i].jws = 1012;
    printf("rmt=%lf\n",local.atom[i].rmt);
    local.atom[i].generateRadialMesh();
    local.atom[i].nspin = 2;
    local.atom[i].evec[0] = local.atom[i].evec[1] = 0.0;
    local.atom[i].evec[2] = 1.0;

    initializeAtom(local.atom[i]);
    local.atom[i].alat = local.atom[i].rmt;
    local.atom[i].efermi = 0.5;
  }

  lsms.chempot = 0.5;
  calculateCoreStates(comm, lsms, local);
  return 0;
}

int initializeNewAlloyBank(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyMixingDesc &alloyDesc, AlloyAtomBank &alloyBank)
{
  
  printf("Initializing new alloy bank\n");

  for(int id = 0, i = 0; i < alloyBank.size(); i++)
  for(int j = 0; j < alloyBank[i].size(); j++, id++) {

    AtomData &atom = alloyBank[i][j];

    printf("Initializing potential %d.\n",id);
    snprintf(atom.header, 40, "Potential Initialized by LSMS_3");

    for(int j=31;j<80;j++) atom.header[j]=' ';
    atom.resizePotential(lsms.global.iprpts);
    atom.ztotss=(Real)alloyDesc[i][j].Z;
    atom.zcorss=(Real)alloyDesc[i][j].Zc;
    atom.zsemss=(Real)alloyDesc[i][j].Zs;
    atom.zvalss=(Real)alloyDesc[i][j].Zv;
    atom.alloy_class=(Real)alloyDesc[i][j].alloy_class;
    atom.vdif=0.0;
    atom.xvalws[0]=atom.xvalws[1]=0.5*atom.zvalss;
    atom.lmax=alloyDesc[i][j].lmax;
    atom.kkrsz=(atom.lmax+1)*(atom.lmax+1);

    atom.xstart=-11.1309674;
    atom.jmt=1001;
    atom.jws=1001;
    atom.nspin=2;
    atom.evec[0]=atom.evec[1]=0.0; atom.evec[2]=1.0;

    initializeAtom(atom);
    atom.alat=atom.rmt;
    atom.efermi=0.5;
  }
  calculateCoreStates(comm, lsms, alloyBank);
  return 0;
}
