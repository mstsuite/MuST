// calculateGaussLegendreNodes calculates the nodes and weights for 
// the Gauss-Legendre quadrature formula of degree n>0
// the algorithm follows
// 
// David H. Bailey, Karthik Jeyabalan and Xiaoye S. Li
// "A Comparison of Three High-Precision Quadrature Schemes"
// Experimental Mathematics 14, 317-329 (2005).
//

#include <cmath>

template<typename R>  // R is a real type e.g. double
void calculateGaussLegendrePoints(R *x, R *w, int n, R xLeft = -1.0, R xRight = +1.0)
{
  const R pi = std::acos(-R(1));

  R Pk[n], PkMinus1[n], PkPlus1[n]; // the Legendre polynomials at points x[i]
  R Ppn[n]; // derivative of P_n(x_i)
  R dX[n]; // correction values for the Newton iteration

  R xHalfWidth = 0.5*(xRight - xLeft);
  R xMidpoint = 0.5*(xRight+xLeft);

  // the x[i] are the roots of the nth degree Legendre Polynomial P_n(x)
  // we search these using a Newton iteration with starting guess
  // x[i] ~ cos(\pi (j - 1/4) / (n + 1/2))

  for(int i=0; i<n; i++)
  {
    x[i]=std::cos(pi*(R(i+1)-0.25)/(R(n)+0.5));
    dX[i]=1000.0;

    while(x[i]+dX[i]/4.0 != x[i])
    {
      Pk[i] = 0.0; // P_{-1}(x_i) = 0
      PkPlus1[i] = 1.0; // P_0(x_i) = 1
      // iterate using the recurence relation for Legendre polynomials:
      // (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}
      // until k=n. I.e we have P_{n+1}, P_n and P_{n-1}
      for(int k=0; k<=n; k++)
      {
        PkMinus1[i]=Pk[i];
        Pk[i]=PkPlus1[i];
        PkPlus1[i]=(R(2*k+1) * x[i] * Pk[i] - R(k) * PkMinus1[i])/R(k+1);
      }

      // derivative of Pk:
      // P'_n(x) = n (x P_n(x) - P_{n-1}(x))/(x^2 -1)
      Ppn[i]=R(n) * (x[i] * Pk[i] - PkMinus1[i])/(x[i]*x[i]-R(1));

      // Newton step x_old <- x, x <- x_old - P_n(x)/P'_n(x)
      dX[i] = Pk[i]/Ppn[i];
      x[i] = x[i] - dX[i];
    }

    // calculate P'_n(x) and P_{n+1}(x) at the final points
      Pk[i] = 0.0; // P_{-1}(x_i) = 0
      PkPlus1[i] = 1.0; // P_0(x_i) = 1
      // iterate using the recurence relation for Legendre polynomials:
      // (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}
      // until k=n. I.e we have P_{n+1}, P_n and P_{n-1}
      for(int k=0; k<=n; k++)
      {
        PkMinus1[i]=Pk[i];
        Pk[i]=PkPlus1[i];
        PkPlus1[i]=(R(2*k+1) * x[i] * Pk[i] - R(k) * PkMinus1[i])/R(k+1);
      }
      // derivative of Pk:
      // P'_n(x) = n (x P_n(x) - P_{n-1}(x))/(x^2 -1)
      Ppn[i]=R(n) * (x[i] * Pk[i] - PkMinus1[i])/(x[i]*x[i]-R(1));

    // calculate the weights:
    w[i] = -R(2)/(R(n+1)*Ppn[i]*PkPlus1[i]);
    // adjust the points and weights for the integration interval
    x[i] = xMidpoint - xHalfWidth*x[i];
    w[i] = xHalfWidth * w[i];
  }
}


