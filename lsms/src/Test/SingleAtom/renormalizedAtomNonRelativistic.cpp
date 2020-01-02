/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// Non Relativistic Renormalized Atom

#include <tuple>
#include <algorithm>

// typedef double Real;
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "../../Misc/bulirschStoerIntegrator.hpp"
#include "../../Misc/rationalFit.hpp"
#include "../../Misc/integrateOneDim.cpp"


#include "calculateXC.cpp"

/*
// integrate a system of n coupled ordinary differential equations
// from x0 to x1 with initial values in y0
// return values at x1 in y1
// rhs is a function that can evaluate the right hand side of the system of ODEs
// anywhere in the interval [x0,x1]
// eps is the desired max error
// returns an error code: 0 -> success
// templated on the type of real numbers
template <typename Rx, typename Ry>
int bulirschStoerIntegrator(Rx x0, Rx x1, Ry *y0,
                            Ry *y1, int n, std::function<void(Rx x, Ry* y, Ry* dy)> rhs,
                            Rx eps=1.0e-12)

/// Given a table of function values r[i] -> f(r[i])
/// find the interpolated value of f(x)
template<typename T>
T interpolate(std::vector<T> &r, std::vector<T> &f, T x)
*/

// Schroedinger RHS for radial electrostatic potential without magnetic field:
// in: r, p[0] , p[1] = dp[0]/dr
// out: dp[0] = p[1], dp[1] = (v + l*(l+1)/(r*r) - e) * p[0]
void schroedingerRHS(Real r, Real *p, Real *dp, std::vector<Real> &r_mesh, std::vector<Real> &vr,
              Real e, Real l)
{
  Real v = interpolate<Real>(r_mesh, vr, r)/r;
  dp[0] = p[1];
  dp[1] = (v + l*(l+1)/(r*r) - e) * p[0];
}


// p ~ r^(l+1); dp/dr ~ (l+1) r^l
template<typename PQType>
void schroedingerBoundaryConditionNearOrigin(Real r, PQType *pdp, Real Z, Real l)
{
  const Real cLight = 274.072;
  Real zeta = 2.0*Z/cLight; // Z * alpha

  pdp[0] = std::pow(r, l+1);
  pdp[1] = (l+1) * std::pow(r, l);
}

template<typename PQType>
void schroedingerBoundaryConditionAtBoundary(Real r, PQType *pdp, int numberOfNodes)
{
  pdp[0] = 0.0; // psi(R) = 0
  pdp[1] = 1.0; // 0.0001; //
  if(numberOfNodes % 2) // for an odd number of nodes change initital slope
    pdp[1] = -pdp[1];
}

Real energyGuess(int n, Real Z, int l)
{
  const Real cLight = 274.072;
  Real zeta = 2.0*Z/cLight; // Z * alpha

  return (2.0*Z*Z)*(-1.0/(2.0*Real(n*n)));
}

void generateRadialMesh(std::vector<Real> &r_mesh, int N, Real r0, Real rN)
{
    if (N != r_mesh.size()) r_mesh.resize(N);
    Real x0 = std::log(r0);
    Real xN = std::log(rN);
    Real h = (xN-x0) / (N-1);
    for(int j=0; j<N; j++)
    {
      Real x_mesh = x0 + (Real)j*h;
      r_mesh[j] = std::exp(x_mesh);
    }
}

int countNodes(Matrix<Real> &y, int c=0)
{
  int n=0;
  bool s = std::signbit(y(c,0));
  for(int i=1; i<y.n_col()-1; i++)
  {
    if(s != std::signbit(y(c,i)))
    {
      n++;
      s = std::signbit(y(c,i));
    }
  }
  return n;
}


void integrateSchroedinger(std::vector<Real> &r_mesh, Real atomicNumber, std::vector<Real> &vr,
                    Real energy, Real l, Matrix<Real> &pdp)
{
  // printf("integrateDirac: energy=%lg\n",energy);
// outward integration
  schroedingerBoundaryConditionNearOrigin(r_mesh[0], &pdp(0,0), Real(atomicNumber), l);

  // for(int i=1; i<iFitting; i++)
  for(int i=1; i<r_mesh.size(); i++)
  {
    /*
    modifiedMidpoint<Real, Real>(r_mesh[i-1], r_mesh[i], &pdp(0,i-1), &pdp(0,i), 2,
                                 [&](Real r, Real *y, Real *dy)
                                 {schroedingerRHS(r, y, dy, r_mesh, vr, energy, l);},
                                 10);
    //*/
    //*
    if(bulirschStoerIntegrator<Real, Real>(r_mesh[i-1], r_mesh[i], &pdp(0,i-1), &pdp(0,i), 2,
                           [&](Real r, Real *y, Real *dy)
                           {schroedingerRHS(r, y, dy, r_mesh, vr, energy, l);},
                           1.0e-12))
    {
      printf("outward integration did not succeed: %d:%lg -> %d:%lg!\n",i-1,r_mesh[i-1],i,r_mesh[i]);
      exit(0);
    }
    //*/
  }
}

void integrateSchroedingerInward(std::vector<Real> &r_mesh, Real atomicNumber, std::vector<Real> &vr,
                                 Real energy, Real l, int principleQuantumNumber, Matrix<Real> &pdp)
{
  // printf("integrateDirac: energy=%lg\n",energy);
// inward integration
  int numberOfNodes = principleQuantumNumber - l - 1;
  schroedingerBoundaryConditionAtBoundary(r_mesh[r_mesh.size()-1], &pdp(0,r_mesh.size()-1), numberOfNodes);

  // for(int i=1; i<iFitting; i++)
  for(int i=r_mesh.size()-1; i>1; i--)
  {
    /*
    modifiedMidpoint<Real, Real>(r_mesh[i], r_mesh[i-1], &pdp(0,i), &pdp(0,i-1), 2,
                                 [&](Real r, Real *y, Real *dy)
                                 {schroedingerRHS(r, y, dy, r_mesh, vr, energy, l);},
                                 10);
    //*/   
    //*
    if(bulirschStoerIntegrator<Real, Real>(r_mesh[i], r_mesh[i-1], &pdp(0,i), &pdp(0,i-1), 2,
                           [&](Real r, Real *y, Real *dy)
                           {schroedingerRHS(r, y, dy, r_mesh, vr, energy, l);},
                           1.0e-12))
    {
      printf("inward integration did not succeed: %d:%lg -> %d:%lg!\n",i-1,r_mesh[i-1],i,r_mesh[i]);
      exit(1);
    }
    //*/
  }
}

void calculateRadialDensity(std::vector<Real> &r_mesh, Matrix<Real> &pdp, std::vector<Real> &rho)
{
  std::vector<Real> rhoIntegrated(r_mesh.size());

  for(int i=0; i<r_mesh.size(); i++)
  {
    rho[i] = std::abs(pdp(0,i))*std::abs(pdp(0,i));
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

Real minimizationQuantityForMatching(Matrix<Real> &pdpOut, Matrix<Real> &pdpIn, int matchingPoint, Real energy)
{
  Real logDerivativeOut = pdpOut(1,matchingPoint)/(std::sqrt(std::abs(energy))*pdpOut(0,matchingPoint));
  Real logDerivativeIn = pdpIn(1,matchingPoint)/(std::sqrt(std::abs(energy))*pdpIn(0,matchingPoint));

  if(std::sqrt(std::abs(energy))<1)
  {
    logDerivativeOut = pdpOut(1,matchingPoint)/(pdpOut(0,matchingPoint));
    logDerivativeIn = pdpIn(1,matchingPoint)/(pdpIn(0,matchingPoint));
  }
  
  // return logDerivativeOut - logDerivativeIn; // log derivative difference
  return (1.0+logDerivativeIn)/(1.0-logDerivativeIn) -  (1.0+logDerivativeOut)/(1.0-logDerivativeOut);
}

Real findSchroedingerEigenvalue(std::vector<Real> &r_mesh, std::vector<Real> &vr, Real atomicNumber,
			 int principalQuantumNumber, int l, Real energy,
			 Matrix<Real> &pdp)
{
  Real energyUpper, energyLower;

  Real matchingQuantity;
  int matchingPoint = 95*r_mesh.size()/100;
  Matrix<Real> pdpInward;
  pdpInward.resize(2,r_mesh.size());
  
  int targetNumberOfNodes = principalQuantumNumber - l - 1;

// find energy bounds

  // printf("find energy bounds\n");
  integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
  integrateSchroedingerInward(r_mesh, atomicNumber, vr, energy, l, principalQuantumNumber, pdpInward);
  matchingQuantity = minimizationQuantityForMatching(pdp,pdpInward,matchingPoint,energy);

  int numNodesP = countNodes(pdp,0);
  int numNodesDP = countNodes(pdp,1);
  int numNodesPUpper, numNodesPLower;

  /*
  energyUpper = energy;
  numNodesPUpper = numNodesP;
  // find lower energy
  while(numNodesP>targetNumberOfNodes)
  {
    energy -= 0.1;;
    integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
    numNodesP = countNodes(pdp,0);
    printf(">nodes: %d target: %d energy: %.15lg\n",numNodesP,targetNumberOfNodes,energy);
  }
  energyLower = energy;
  numNodesPLower = numNodesP;
  // find upper energy
  // energy = energyUpper;
  // numNodesP =  numNodesPUpper;
  while(numNodesP<targetNumberOfNodes)
  {
    energy += 0.1;;
    integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
    numNodesP = countNodes(pdp,0);
    printf("<nodes: %d target: %d energy: %.15lg\n",numNodesP,targetNumberOfNodes,energy);
  }
  energyUpper = energy;
  numNodesPUpper = numNodesP;
  */
  //*
  if(numNodesP>targetNumberOfNodes)
  {
    energyUpper = energy;
    numNodesPUpper = numNodesP;
    while(numNodesP>targetNumberOfNodes)
    {
      energy -= 0.5;
      integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
      numNodesP = countNodes(pdp,0);
    }
    energyLower = energy;
    numNodesPLower = numNodesP;
  } else {
    energyLower = energy;
    numNodesPLower = numNodesP;
    while(numNodesP<=targetNumberOfNodes)
    {
      energy += 0.5;
      integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
      numNodesP = countNodes(pdp,0);
    }
    energyUpper = energy;
    numNodesPUpper = numNodesP;
  }
  //*/
  
  Real energyEpsilon = 1.0e-3; // 1.0e-15
  /*
  while(std::abs((energyUpper - energyLower)/energy) > energyEpsilon)
  {
    energy = energyLower + 0.5*(energyUpper-energyLower);
    integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
    numNodesP = countNodes(pdp,0);
    if(numNodesP>targetNumberOfNodes)
    {
      energyUpper = energy;
      numNodesPUpper = numNodesP;
    } else {
      energyLower = energy;
      numNodesPLower = numNodesP;
    }
    // printf("%lf %lf  %lg\n",energyLower, energyUpper, energyUpper - energyLower);
  }
  //*/

  integrateSchroedinger(r_mesh, atomicNumber, vr, energyUpper, l, pdp);
  integrateSchroedingerInward(r_mesh, atomicNumber, vr, energyUpper, l, principalQuantumNumber, pdpInward);
  Real matchingQuantityUpper = minimizationQuantityForMatching(pdp,pdpInward,matchingPoint,energyUpper);

  integrateSchroedinger(r_mesh, atomicNumber, vr, energyLower, l, pdp);
  integrateSchroedingerInward(r_mesh, atomicNumber, vr, energyLower, l, principalQuantumNumber, pdpInward);
  Real matchingQuantityLower = minimizationQuantityForMatching(pdp,pdpInward,matchingPoint,energyLower);

  // printf("***\n");
  // while(std::signbit(logDerivativeDifferenceLower) == std::signbit(logDerivativeDifferenceUpper))
  // {

  //   Real d = (logDerivativeDifferenceUpper - logDerivativeDifferenceLower); // /(energyUpper - energyLower);
  //   printf("eL lDDL  eU lDDU -> d: %.15lg %.15lg   %.15lg %.15lg  -> %.15lg\n",energyLower, logDerivativeDifferenceLower, energyUpper, logDerivativeDifferenceUpper, d);

  //   if(logDerivativeDifferenceUpper > 0.0)
  //   {
  //     if(d > 0.0)
  //     {
  //       energyLower -= logDerivativeDifferenceLower/d;
  //       integrateSchroedinger(r_mesh, atomicNumber, vr, energyLower, l, pdp);
  //       integrateSchroedingerInward(r_mesh, atomicNumber, vr, energyLower, l, principalQuantumNumber, pdpInward);
  //       logDerivativeDifferenceLower = pdp(1,matchingPoint)/pdp(0,matchingPoint)
  //         - pdpInward(1,matchingPoint)/pdpInward(0,matchingPoint);
  //     } else {
  //       energyUpper += logDerivativeDifferenceUpper/d;
  //       integrateSchroedinger(r_mesh, atomicNumber, vr, energyUpper, l, pdp);
  //       integrateSchroedingerInward(r_mesh, atomicNumber, vr, energyUpper, l, principalQuantumNumber, pdpInward);
  //       logDerivativeDifferenceUpper = pdp(1,matchingPoint)/pdp(0,matchingPoint)
  //         - pdpInward(1,matchingPoint)/pdpInward(0,matchingPoint);
  //     }
  //   } else {
  //     if(d < 0.0)
  //     {
  //       energyLower -= logDerivativeDifferenceLower/d;
  //       integrateSchroedinger(r_mesh, atomicNumber, vr, energyLower, l, pdp);
  //       integrateSchroedingerInward(r_mesh, atomicNumber, vr, energyLower, l, principalQuantumNumber, pdpInward);
  //       logDerivativeDifferenceLower = pdp(1,matchingPoint)/pdp(0,matchingPoint)
  //         - pdpInward(1,matchingPoint)/pdpInward(0,matchingPoint);
  //     } else {
  //       energyUpper += logDerivativeDifferenceUpper/d;
  //       integrateSchroedinger(r_mesh, atomicNumber, vr, energyUpper, l, pdp);
  //       integrateSchroedingerInward(r_mesh, atomicNumber, vr, energyUpper, l, principalQuantumNumber, pdpInward);
  //       logDerivativeDifferenceUpper = pdp(1,matchingPoint)/pdp(0,matchingPoint)
  //         - pdpInward(1,matchingPoint)/pdpInward(0,matchingPoint);
  //     }
  //   }
  // }

  energy = energyLower + 0.5*(energyUpper-energyLower);
  while(std::abs((energyUpper - energyLower)/energy) > energyEpsilon)
  {
    energy = energyLower + 0.5*(energyUpper-energyLower);
    integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
    integrateSchroedingerInward(r_mesh, atomicNumber, vr, energy, l, principalQuantumNumber, pdpInward);
    matchingQuantity = minimizationQuantityForMatching(pdp,pdpInward,matchingPoint,energy);
    
    Real f = pdp(0,matchingPoint)/pdpInward(0,matchingPoint);
    for(int i=matchingPoint; i<r_mesh.size(); i++)
    {
      pdp(0,i) = f*pdpInward(0,i);
      pdp(1,i) = f*pdpInward(1,i);
    }
    numNodesP = countNodes(pdp,0);
    if((std::signbit(matchingQuantityLower) == std::signbit(matchingQuantity))
       || (numNodesP<targetNumberOfNodes))
    {
      matchingQuantityLower = matchingQuantity;
      energyLower = energy;
    } else {
      matchingQuantityUpper = matchingQuantity;
      energyUpper = energy;
    }
  }

  energy = energyLower + 0.00001; //energyUpper;
  Real energyOld = energyLower;
  matchingQuantity = matchingQuantityUpper;
  Real matchingQuantityOld = matchingQuantityLower;


  while(std::abs(matchingQuantity) > 1.0e-8)
  {
    Real d = (energy - energyOld)/(matchingQuantity - matchingQuantityOld);
    Real tEnergy = energy -  matchingQuantity*d;
  
    energyOld = energy;
    matchingQuantityOld = matchingQuantity;
     energy = tEnergy;

    integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
    integrateSchroedingerInward(r_mesh, atomicNumber, vr, energy, l, principalQuantumNumber, pdpInward);
    matchingQuantity = minimizationQuantityForMatching(pdp,pdpInward,matchingPoint,energy);
  //   Real f = pdp(0,matchingPoint)/pdpInward(0,matchingPoint);
  //   for(int i=matchingPoint; i<r_mesh.size(); i++)
  //   {
  //     pdp(0,i) = f*pdpInward(0,i);
  //     pdp(1,i) = f*pdpInward(1,i);
  //   }
  //   numNodesP = countNodes(pdp,0);
  //   if((std::signbit(logDerivativeDifferenceLower) == std::signbit(logDerivativeDifference))
  //      || (numNodesP<targetNumberOfNodes))
  //   {
  //     logDerivativeDifferenceLower = logDerivativeDifference;
  //     energyLower = energy;
  //   } else {
  //     logDerivativeDifferenceUpper = logDerivativeDifference;
  //     energyUpper = energy;
  //   }
    // printf("energy = %.15lg  logDerivativeDifference = %.15lg\n",energy,matchingQuantity);
  }
  

  Real f = pdp(0,matchingPoint)/pdpInward(0,matchingPoint);
  for(int i=matchingPoint; i<r_mesh.size(); i++)
  {
    pdp(0,i) = f*pdpInward(0,i);
    pdp(1,i) = f*pdpInward(1,i);
  }

  // exit(0);
  return energy;
}

class AtomOrbital {
public:
  int n,l;
  Real energy;
  std::vector<Real> rho; // radial charge density contribution from an electron in this level
  // char type; // 'C' for core lectron, 'S' for semi-core and 'V' for valence electrons
} ;


// generate a list of (n, l) orbitals
void initOrbitals(int Z, std::vector<std::tuple<int, int>> &orbitals)
{
  orbitals.clear();
  // He: 1s^2
  orbitals.push_back(std::make_tuple<int,int>(1, 0));
  if(Z<3) return;
  // Ne: [He] 2s^2 2p^6
  orbitals.push_back(std::make_tuple<int,int>(2, 0));
  orbitals.push_back(std::make_tuple<int,int>(2, 1));
  if(Z<11) return;
  // Ar: [Ne] 3s^2 3p^6
  orbitals.push_back(std::make_tuple<int,int>(3, 0));
  orbitals.push_back(std::make_tuple<int,int>(3, 1));
  if(Z<19) return;
  // Kr: [Ar] 3d^10 4s^2 4p^6
  orbitals.push_back(std::make_tuple<int,int>(3, 2));
  orbitals.push_back(std::make_tuple<int,int>(4, 0));
  orbitals.push_back(std::make_tuple<int,int>(4, 1));
  if(Z<37) return;
  // Xe: [Kr] 4d^10 5s^2 5p^6
  orbitals.push_back(std::make_tuple<int,int>(4, 2));
  orbitals.push_back(std::make_tuple<int,int>(5, 0));
  orbitals.push_back(std::make_tuple<int,int>(5, 1));
  if(Z<55) return;
  // Rn: [Xe] 4f^14 5d^10 6s^2 6p^6
  orbitals.push_back(std::make_tuple<int,int>(4, 3));
  orbitals.push_back(std::make_tuple<int,int>(5, 2));
  orbitals.push_back(std::make_tuple<int,int>(6, 0));
  orbitals.push_back(std::make_tuple<int,int>(6, 1));
  if(Z<87) return;
  // Og: [Rn] 5f^14 6d^10 7s^2 7p^6
  orbitals.push_back(std::make_tuple<int,int>(5, 3));
  orbitals.push_back(std::make_tuple<int,int>(6, 2));
  orbitals.push_back(std::make_tuple<int,int>(7, 0));
  orbitals.push_back(std::make_tuple<int,int>(7, 1));
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

void calculateOrbitals(std::vector<Real> &r_mesh, std::vector<Real> &vr, int atomicNumber,
                       std::vector<std::tuple<int, int>> &orbitals,
                       std::vector<AtomOrbital> &orbitalEnergiesAndDensities)
{
  printf(" n  l  energy\n");
  for(int i=0; i<orbitals.size(); i++)
  {
    Matrix<Real> pdp;
    // int principalQuantumNumber = orbitalEnergiesAndDensities[i].n; //std::get<0>(orbitals[i]);
    // int l = orbitalEnergiesAndDensities[i].l; //std::get<1>(orbitals[i]);
    int principalQuantumNumber = std::get<0>(orbitals[i]);
    int l = std::get<1>(orbitals[i]);
    pdp.resize(2,r_mesh.size());
  
    Real energy = energyGuess(principalQuantumNumber, atomicNumber, l);
    // if(orbitalEnergiesAndDensities[i].energy != 0.0)
    //  energy = orbitalEnergiesAndDensities[i].energy;
  // printf("# energy: %lg Ry\n",energy);

    orbitalEnergiesAndDensities[i].n = principalQuantumNumber;
    orbitalEnergiesAndDensities[i].l = l;
    orbitalEnergiesAndDensities[i].energy = findSchroedingerEigenvalue(r_mesh, vr, atomicNumber,
					     principalQuantumNumber, l, energy,
					     pdp);

    printf(" %d %2d %lg Ry",principalQuantumNumber, l,orbitalEnergiesAndDensities[i].energy);

    orbitalEnergiesAndDensities[i].rho.resize(r_mesh.size());
    calculateRadialDensity(r_mesh, pdp, orbitalEnergiesAndDensities[i].rho);

    if(orbitalEnergiesAndDensities[i].rho[r_mesh.size()-1] > 0.0001)
    {
      printf(" !!\n");
    } else {
      printf("\n");
    }
    
    int targetNumberOfNodes = principalQuantumNumber - l - 1;
    int numNodesP = countNodes(pdp,0);
    int numNodesQ = countNodes(pdp,1);
  }
  
  std::sort(orbitalEnergiesAndDensities.begin(), orbitalEnergiesAndDensities.end(),
	    [](AtomOrbital const & a, AtomOrbital const &b){return a.energy < b.energy;});

}

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

void calculateSingeOrbitalLogDerivative(int atomicNumber, int l, Real energy,
                                        Real energyTop, Real dEnergy, Real atomRadius)
{
  std::vector<Real> r_mesh, vr, rho, rhoInward;
  Matrix<Real> pdp, pdpInward;
  generateRadialMesh(r_mesh, 1500, 1.5e-10, atomRadius);
  vr.resize(r_mesh.size());
  rho.resize(r_mesh.size());
  rhoInward.resize(r_mesh.size());
  pdp.resize(2,r_mesh.size());
  pdpInward.resize(2,r_mesh.size());
  for(int i=0; i<r_mesh.size(); i++) vr[i] = -2.0*Real(atomicNumber);

  int matchingPoint = 95*r_mesh.size()/100;
  
  FILE *outf=fopen("logDerivative.out","w");
  fprintf(outf,"# Atomic Number: %d\n",atomicNumber);
  fprintf(outf,"# energy derivOut derivIn difference nodesOut nodesIn\n");
  for(; energy<energyTop; energy += dEnergy)
  {
    integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
    calculateRadialDensity(r_mesh, pdp, rho);
    integrateSchroedingerInward(r_mesh, atomicNumber, vr, energy, l, 1, pdpInward);
    calculateRadialDensity(r_mesh, pdpInward, rhoInward);
 
    fprintf(outf,"%lg %lg %lg %lg %d %d\n",energy, pdp[1,matchingPoint]/pdp[0,matchingPoint],
            pdpInward[1,matchingPoint]/pdpInward[0,matchingPoint],
            minimizationQuantityForMatching(pdp,pdpInward,matchingPoint,energy),
            countNodes(pdp,0), countNodes(pdpInward,0));
  }
  fclose(outf);
}

void calculateSingeOrbital(int atomicNumber, int l, Real energy, Real atomRadius)
{
  std::vector<Real> r_mesh, vr, rho, rhoInward;
  Matrix<Real> pdp, pdpInward;
  generateRadialMesh(r_mesh, 1500, 1.5e-10, atomRadius);
  vr.resize(r_mesh.size());
  rho.resize(r_mesh.size());
  rhoInward.resize(r_mesh.size());
  pdp.resize(2,r_mesh.size());
  pdpInward.resize(2,r_mesh.size());
  for(int i=0; i<r_mesh.size(); i++) vr[i] = -2.0*Real(atomicNumber);
  
  integrateSchroedinger(r_mesh, atomicNumber, vr, energy, l, pdp);
  calculateRadialDensity(r_mesh, pdp, rho);
  integrateSchroedingerInward(r_mesh, atomicNumber, vr, energy, l, 1, pdpInward);
  calculateRadialDensity(r_mesh, pdpInward, rhoInward);

  printf("number of nodes: %d\n",countNodes(pdp,0));
  
  FILE *outf=fopen("orbital.out","w");
  fprintf(outf,"# Atomic Number: %d\n",atomicNumber);
  fprintf(outf,"# i r[i] rho rhoInward\n");
  Real exc;
  for(int i=0; i<r_mesh.size(); i++)
  {
    fprintf(outf,"%5d %lg %lg %lg\n",i, r_mesh[i], rho[i], rhoInward[i]);
  }
  fclose(outf);
}

void printUsage(const char *name, Real R)
{
  printf("Usage: %s Z [R]\n",name);
  printf("       Z: atomic number\n");
  printf("       R: atomic sphere radius (optional, default=%lf)\n",R);
  printf("       %s o Z l energy [R]\n",name);
  printf("       %s d Z l energy energyTop dEnergy [R]\n",name);
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

  if(argc != 2 && argc != 3 && argc != 5 && argc != 6 && argc != 7 && argc != 8)
  {
    printUsage(argv[0], atomRadius);
    exit(1);
  }
  if(argv[1][0]=='o')
  {
    mode = 'o'; // one orbital
    atomicNumber = atoi(argv[2]);
    l = atoi(argv[3]);
    energy = atof(argv[4]);
    if (argc == 6)
      atomRadius = atof(argv[2]);
    calculateSingeOrbital(atomicNumber, l, energy, atomRadius);
    return 0;
  } else if(argv[1][0] =='d') {
    mode = 'd';
    atomicNumber = atoi(argv[2]);
    l = atoi(argv[3]);
    energy = atof(argv[4]);
    Real energyTop = atof(argv[5]);
    Real dEnergy =  atof(argv[6]);
    if (argc == 8)
      atomRadius = atof(argv[7]);
    calculateSingeOrbitalLogDerivative(atomicNumber, l, energy,
                                       energyTop, dEnergy, atomRadius);
    return 0;
  } else {
    atomicNumber = atoi(argv[1]);
    if(argc == 3)
      atomRadius = atof(argv[2]);
  }

  int maxPrincipalQuantumNumber;

  std::vector<std::tuple<int, int>> orbitals;
  initOrbitals(atomicNumber, orbitals);

  printf("number of orbitals to be computed: %d\n",orbitals.size());
  
  std::vector<Real> r_mesh, vr, vXC;
  // std::vector<Matrix<Real> > pqs;
  std::vector<AtomOrbital> orbitalEnergiesAndDensities;
  // pqs.resize(orbitals.size());
  orbitalEnergiesAndDensities.resize(orbitals.size());
  for(int i=0; i<orbitals.size(); i++)
  {
    orbitalEnergiesAndDensities[i].n = std::get<0>(orbitals[i]);
    orbitalEnergiesAndDensities[i].l = std::get<1>(orbitals[i]);
    orbitalEnergiesAndDensities[i].energy = 0.0;
  }

// initialize r_mesh (unit of length is the Bohr radius)
  generateRadialMesh(r_mesh, 1500, 1.5e-10, atomRadius);
  std::vector<Real> rhotot, rhonew;
  rhotot.resize(r_mesh.size());
  rhonew.resize(r_mesh.size());
  
// initialize Z/r potential e^2=2
  vr.resize(r_mesh.size());
  vXC.resize(r_mesh.size());
  for(int i=0; i<r_mesh.size(); i++) vr[i] = -2.0*Real(atomicNumber);

  calculateOrbitals(r_mesh, vr, atomicNumber, orbitals, orbitalEnergiesAndDensities);
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
    calculateOrbitals(r_mesh, vr, atomicNumber, orbitals, orbitalEnergiesAndDensities);
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
