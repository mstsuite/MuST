// routines to calculate the associated Legendre functions P_{lm} needed to calculate 
// spherical harmonics
// This also provides the coeficients c_{lm} to calculate spherical harmonics:
// Y_{lm}(\theta, \phi) = c_{lm} P_{lm}(\cos\theta) e^{im\phi}
//
// The new routine calculates normalized versions of the Legendre polynomials
// such that
// Y_{lm}(\theta, \phi) = \sqrt{1/2} \bar{P}_{lm}(\cos\theta) e^{im\phi}
//
// i.e. c_{lm} = 1
//

#include <cmath>
#include "associatedLegendreFunction.hpp"

extern "C" {
void plm_normalized_(int *lmax, double *x, double *plm)
{
  associatedLegendreFunctionNormalized<double>(*x, *lmax, plm);
}

void ylm_coefficients_(int *lmax, double *clm)
{
  for(int i=0; i<((*lmax) +1)*((*lmax) +2)/2; i++)
    clm[i]=1.0;
    // clm[i]=std::sqrt(0.5);
}
}

