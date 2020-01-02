#include "Misc/integrateOneDim.cpp"
#include <cmath>

extern "C"
{
  void zeropt_(Real *ezpt, Real *tpzpt, Real *atvol, Real *ztotss);
}

void localTotalEnergy(LSMSSystemParameters &lsms, AtomData &atom, 
//output:
Real &energy, Real &pressure)
{
  Real eigenvalueSum=0.0;
  Real kineticEnergy=0.0;
  Real coulombEnergy=0.0;
  Real xcEnergy=0.0;
  Real ezpt=0.0;
  Real tpzpt=0.0;
  std::vector<Real> integral(atom.r_mesh.size()+1);
  std::vector<Real> integrand(atom.r_mesh.size()+1);
  std::vector<Real> grid0(atom.r_mesh.size()+1);
  std::vector<Real> gridSqrt(atom.r_mesh.size()+1);
  std::vector<Real> gridCbrt(atom.r_mesh.size()+1);
  energy =0.0;
  pressure=0.0;
  
  Real rSphere;
  switch (lsms.mtasa)
  {
  case 1:
    rSphere = atom.rws;
    break;
  case 2:
    rSphere = atom.rws;
    break;
  default:
    rSphere = atom.rInscribed;
  }

//        ============================================================
//        calculate the zeropoint energy...........................
//        ============================================================
//        ------------------------------------------------------------
  zeropt_(&ezpt,&tpzpt,&atom.omegaWS,&atom.ztotss);

// calculate kinetic energy T:
// T = \sum_{Core} e_i                            -- (1)
//   + \int^{E_F} e n(e) de                       -- (2)
//   - \int \rho(r) (v_{Coulomb} + v_{xc}) d^3r   -- (3)
//   + \int m(r) (B_{xc} + B_{external}) d^3      -- (4)

  if(lsms.n_spin_pola==1)
  {
    eigenvalueSum=atom.evalsum[0]+atom.esemv[0];

    kineticEnergy = atom.ecorv[0]+atom.esemv[0]; // (1)
    kineticEnergy += atom.evalsum[0]; // (2)

    // set up integrand for (3)
    for(int i=0; i<atom.r_mesh.size(); i++)
    {
      grid0[i+1]=atom.r_mesh[i];
      gridSqrt[i+1] = std::sqrt(atom.r_mesh[i]);
      gridCbrt[i+1] = std::cbrt(atom.r_mesh[i]);
      integrand[i+1]=(atom.rhoNew(i,0)*atom.vr(i,0))/(atom.r_mesh[i]);
      // integrand[i+1]=2.0*gridSqrt[i+1]*(atom.rhoNew(i,0)*atom.vr(i,0))/(atom.r_mesh[i]);
      // integrand[i+1]=3.0*gridCbrt[i+1]*gridCbrt[i+1]*(atom.rhoNew(i,0)*atom.vr(i,0))/(atom.r_mesh[i]);
    }
  } else {  // spin polarized
    eigenvalueSum=atom.evalsum[0]+atom.evalsum[1]+atom.esemv[0]+atom.esemv[1];

    kineticEnergy = atom.ecorv[0]+atom.ecorv[1]+atom.esemv[0]+atom.esemv[1]; // (1)
    kineticEnergy += atom.evalsum[0]+atom.evalsum[1]; // (2)

    // set up integrand for (3)
    for(int i=0; i<atom.r_mesh.size(); i++)
    {
      integrand[i+1]=(atom.rhoNew(i,0)*atom.vr(i,0)+atom.rhoNew(i,1)*atom.vr(i,1))/(atom.r_mesh[i]);
      grid0[i+1]=atom.r_mesh[i];
      gridSqrt[i+1] = std::sqrt(atom.r_mesh[i]);
    }
  }

  grid0[0]=gridSqrt[0]=0.0;
  RationalFit<Real> fit;
  fit.set(grid0,integrand,2);
  // fit.set(gridSqrt,integrand,2);
  // fit.set(gridCbrt, integrand, 2);
  integrand[0]=fit(0.0);
  kineticEnergy -= integrateOneDim(grid0, integrand, integral, rSphere); // (3)
  // kineticEnergy -= integrateOneDim<0>(gridSqrt, integrand, integral, std::sqrt(rSphere)); // (3)
  // kineticEnergy -= integrateOneDim<11>(gridCbrt, integrand, integral, std::cbrt(rSphere)); // (3)
  if (lsms.global.iprint >= 0)
  {
    printf("evssum                      = %35.25lf Ry\n",eigenvalueSum);
    printf("kinetic Energy              = %35.25lf Ry\n",kineticEnergy);
  }

// calculate Coulomb Energy E_C:
// E_C = \int \rho(r) v_{Coulomb} d^3r          -- (5)
//     + E_{Madelung}                           -- (6) // external to this routine
// break the Coulomb integral in two parts:
//  \int \rho(r) v_{Coulomb} d^3r
//      = \int \rho(r) \int^r rho(r')/r' dr' dr -- (5a)
//      + \int rho(r) Z/r dr                    -- (5b)

  if(lsms.n_spin_pola==1)
  {
    for(int i=0; i<atom.r_mesh.size(); i++)
    {
      integrand[i+1]=atom.rhoNew(i,0);
    }
  } else { //spin polarized
    for(int i=0; i<atom.r_mesh.size(); i++)
    {
      integrand[i+1]=(atom.rhoNew(i,0)+atom.rhoNew(i,1));
    }
  }
  fit.set(grid0,integrand,2);
  integrand[0]=fit(0.0);
  integrateOneDim(grid0, integrand, integral);
  if(lsms.n_spin_pola==1)
  {
    for(int i=0; i<atom.r_mesh.size(); i++)
      integral[i+1]=integral[i+1]*2.0*(atom.rhoNew(i,0))/(atom.r_mesh[i]);
  } else { // spin polarized
    for(int i=0; i<atom.r_mesh.size(); i++)
      integral[i+1]=integral[i+1]*2.0*(atom.rhoNew(i,0)+atom.rhoNew(i,1))/(atom.r_mesh[i]);
  }
  fit.set(grid0,integral,2);
  integral[0]=fit(0.0);
  Real erho = integrateOneDim(grid0, integral, integrand, rSphere); // (5a)
  if (lsms.global.iprint >= 0)
    printf("erho                        = %35.25lf Ry\n",erho);

  if(lsms.n_spin_pola==1)
  {
    for(int i=0; i<atom.r_mesh.size(); i++)
    {
      integrand[i+1]=2.0*(atom.rhoNew(i,0))*atom.ztotss/(atom.r_mesh[i]);
    }
  } else { // spin polarized
    for(int i=0; i<atom.r_mesh.size(); i++)
    {
      integrand[i+1]=2.0*(atom.rhoNew(i,0)+atom.rhoNew(i,1))*atom.ztotss/(atom.r_mesh[i]);
    }
  }
  fit.set(grid0,integrand,2);
  integrand[0]=fit(0.0);
  Real ezrho=-integrateOneDim(grid0, integrand, integral, rSphere); // (5b)
  if (lsms.global.iprint >= 0)
    printf("ezrho                       = %35.25lf Ry\n",ezrho);

  coulombEnergy = erho+ezrho; // (5)
  if (lsms.global.iprint >= 0)
    printf("Coulomb Energy              = %35.25lf Ry\n",coulombEnergy);

// Exchange-Correlation energy                  -- (7)

  if(lsms.xcFunctional[0]==0)  // for the built in xc functionals:
  {
    if(lsms.n_spin_pola==1)
    {
      for(int i=0; i<atom.r_mesh.size(); i++)
      { 
        integrand[i+1]=(atom.rhoNew(i,0)*atom.exchangeCorrelationEnergy(i,0));
      }
    } else { // spin polarized
      for(int i=0; i<atom.r_mesh.size(); i++)
      { 
        integrand[i+1]=(atom.rhoNew(i,0)*atom.exchangeCorrelationEnergy(i,0)
                        +atom.rhoNew(i,1)*atom.exchangeCorrelationEnergy(i,1));
      }
    }
    fit.set(grid0,integrand,2);
    integrand[0]=fit(0.0);
    xcEnergy=integrateOneDim(grid0, integrand, integral, rSphere); // (7)
  } else if(lsms.xcFunctional[0]==1) { // for libxc functionals
    if(lsms.n_spin_pola==1)
    {
      for(int i=0; i<atom.r_mesh.size(); i++)
      { 
        integrand[i+1]=atom.exchangeCorrelationEnergy(i,0)*(atom.rhoNew(i,0));
      }
    } else { // spin polarized
      for(int i=0; i<atom.r_mesh.size(); i++)
      { 
        integrand[i+1]=atom.exchangeCorrelationEnergy(i,0)*(atom.rhoNew(i,0)+atom.rhoNew(i,1));
      }
    }
    fit.set(grid0,integrand,2);
    integrand[0]=fit(0.0);
    xcEnergy=integrateOneDim(grid0, integrand, integral, rSphere); // (7)
  } else {  // unknown functional (we never should arrive here!
    printf("Unknown xc function in localTotalEnergy!\n");
    exit(1);
  }
  if (lsms.global.iprint >= 0)
  {
    printf("Exchange-Correlation Energy = %35.25lf Ry\n", xcEnergy);

    printf("ezpt                        = %35.25lf Ry\n\n",ezpt);
  }
// add all energy contributions:

  energy += kineticEnergy + coulombEnergy + xcEnergy + ezpt;
}
