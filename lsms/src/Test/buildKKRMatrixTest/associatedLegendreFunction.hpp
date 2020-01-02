// calculate the associated Legendre Function
// P_{lm}(x) = \sqrt{2 (2l+1)\frac{(l-m)!}{(l+m)!}} x^m \frac{d^m}{dx^m} P_l(x)
// where P_l(x) is the Legendre Polynomial of degree l
//
#ifndef ASSOCIATED_LEGENDRE_FUNCTIONS_HPP
#define ASSOCIATED_LEGENDRE_FUNCTIONS_HPP
#include <cmath>

// for the normalized associated Legendre functions \bar{P}_{lm}
// such that the spherical harmonics are:
// Y_lm(\theta, \phi) = \bar{P}_{lm}(\cos \theta) e^{i m \phi}
// use the recursion relation:
// P_{00}(x) = \sqrt{1/4\pi}
//
// i.e \bar{P}_{lm}=\sqrt{\frac{(2l+1)(l-m)!}{4\pi(l+m)!}} P_{lm}
//
// Plm is a 1-d array that will contain all the values of P_{lm}(x) from P_{00} to P_{l_{max} l_{max}}
// the index into this array is Plm[l*(l+1)/2 + m]
//

inline int plmIdx(int l, int m)
{ return l*(l+1)/2+m; }

template<typename R> // R has to be real data type (default double)
inline void associatedLegendreFunctionNormalized(R x, int lmax, R *Plm)
{
  const R pi = std::acos(-R(1));
  // y = \sqrt{1-x^2}
  R y = std::sqrt(R(1)-x*x);
  // initialize the first entry
  // Plm[0]=std::sqrt(R(1)/(R(2)*pi));
  Plm[0]=std::sqrt(R(1)/(R(4)*pi));

  if(lmax<1) return;

  for(int m=1; m<=lmax; m++)
  {
    // \bar{P}_{mm} = - \sqrt{\frac{2m+1}{2m}} y \bar{P}_{m-1, m-1}
    Plm[plmIdx(m,m)] = - std::sqrt(R(2*m+1)/R(2*m)) * y * Plm[plmIdx(m-1,m-1)];
    // \bar{P}_{mm-1} = \sqrt{2 m + 1} x \bar{P}_{m-1, m-1}
    Plm[plmIdx(m,m-1)] = std::sqrt(R(2*m+1)) * x * Plm[plmIdx(m-1,m-1)]; 
  }

  for(int m=0; m<lmax; m++)
  {
    for(int l=m+2; l<=lmax; l++)
    {
      // \bar{P}_{lm} = a_{lm} (x \bar{P}_{l-1. m} - b_{lm} \bar{P}_{l-2, m})
      // a_{lm} = \sqrt{\frac{(4 l^2 - 1)(l^2 - m^2)}}
      // b_{lm} = \sqrt{\frac{(l -1)^2 - m^2}{4 (l-1)^2 -1}}
      R a_lm = std::sqrt(R(4*l*l-1)/R(l*l - m*m));
      R b_lm = std::sqrt(R((l-1)*(l-1) - m*m)/R(4*(l-1)*(l-1)-1));
      Plm[plmIdx(l,m)] = a_lm * (x * Plm[plmIdx(l-1,m)] - b_lm * Plm[plmIdx(l-2,m)]);
    }
  }
}

// recursion for the unnormalized associated Legendre polynomials
template<typename R> // R has to be real data type (default double)
void associatedLegendreFunction(R x, int lmax, R *Plm)
{
  // y = \sqrt{1-x^2}
  R y = std::sqrt(R(1)-x*x);
  // initialize the first entry
  Plm[0]=R(1);
  if(lmax<1) return;
  Plm[1]=x;
  Plm[2]=-y;
  for(int m=2; m<lmax+1; m++)
  {
    // P_{mm}(x)= -(2m-1) y P_{m-1, m-1}(x)
    Plm[plmIdx(m,m)] = -R(2*m-1) * y * Plm[plmIdx(m-1,m-1)];
    // P_{mm-1} = (2m-1) x P_{m-1,m-1}
    Plm[plmIdx(m,m-1)] = R(2*m-1) * x * Plm[plmIdx(m-1,m-1)];
//    // P_{mm}(x)= -(2m+1) y P_{m-1, m-1}(x)
//    Plm[plmIdx(m,m)] = -R(2*m+1) * y * Plm[plmIdx(m-1,m-1)];
//    // P_{mm-1} = (2m+1) x P_{m-1,m-1}
//    Plm[plmIdx(m,m-1)] = R(2*m+1) * x * Plm[plmIdx(m-1,m-1)];
  }
  for(int m=0; m<lmax; m++)
  {
    for(int l=m+2; l<=lmax; l++)
    {
      // (l - m) P_{lm} = (2l - 1) x P_{l-1. m} - (l + m - 1) P_{l-2, m}
      Plm[plmIdx(l,m)] = (R(2*l-1) * x * Plm[plmIdx(l-1,m)] - R(l+m-1)*Plm[plmIdx(l-2,m)])/R(l-m);
    }
  }
}

#endif
