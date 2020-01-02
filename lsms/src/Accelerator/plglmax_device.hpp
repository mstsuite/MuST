#include <stdio.h>

/*     subroutine plglmax(lmax,x,plm)
c     ****************************************************************
c     Associated Legendre function....................................
c     Calclates all the p(l,m)'s up to lmax...........................
c     based on the formulae given in "Numerical Recipes" pages 180-183
c     (Equations 6.6.7, 6.6.8 and 6.6.9...............................
c     W. H. Press, B. P. Flannery, S A Teukolsky and W. T. Vetterling.
c     Cambridge Univ Press 1986.......................................
c     ****************************************************************
c
c Inputs:   lmax    integer scalar, max l for Plm
c           x       real*8 scalar, argument of Plm(x)
c Returns:  plm     real*8 array of ((lmax+1)*(lmax+2)/2), Plm(x)
c
*/

__device__
__inline__ void plglmax_device(int lmax, double x, double *plm)
{

  double zero=0;
  double one=1;
  double two=2;
  double tol=1e-12;

  if(one-fabs(x) <= tol) {
    int lend=(lmax+1)*(lmax+2)/2;
    for(int l=threadIdx.x;l<lend;l+=blockDim.x) //Parallel
    {
      plm[l]=0;
    }
    if(x < zero) {
      for(int l=threadIdx.x;l<=lmax;l+=blockDim.x) //Parallel
      {
        int i=(l+1)*l/2; 
        plm[i]=one-two*(l%2);
      }
    } else {
      for(int l=threadIdx.x;l<=lmax;l+=blockDim.x) //Parallel
      {
        int i=(l+1)*l/2; 
        plm[i]=one;
      }
    }
    return;
  }
  if(threadIdx.x!=0) return; //do rest in serial
  //begin calculation of p(l,m)'s...................................
  if(lmax == zero) {
    //special case lmax=0..........................................
    plm[0]=one;
  } else {

    //        =============================================================
    //       minus sign added to be consistant with Numerical Recipes
    //       which has (-1)^m factor in plm :  July 97  by xgz.............
    //        =============================================================
    double somx2=-sqrt((one-x)*(one+x));
    if(lmax == one ) {
      //         ==========================================================
      //         special case lmax=1.......................................
      plm[0]=one;
      plm[1]=x;
      plm[2]=somx2;
    } else {
      for(int m=0;m<=lmax;m++) {

        //              =======================================================
        //                                       m       m
        //              calculate the first two P   and P
        //                                       m       m+1
        //              =======================================================

        if(m==zero) {
          plm[0]=one;
          plm[1]=x;
        } else {
          double pmm=somx2;
          double fact=one;
          for(int i=2;i<=m;i++) {
            fact=fact+two;
            pmm=pmm*fact*somx2;
          }
          int mm=(m+1)*(m+2)/2;
          plm[mm-1]=pmm;

          if( mm+m+1 <= (lmax+1)*(lmax+2)/2 ) {
            plm[mm+m]=x*(2*m+1)*pmm;
          }

        } //end else
        //              =======================================================
        //                                 m        m
        //              calculate the rest P     to P
        //                                  m+2      lmax
        //              =======================================================
        int ll=(m+2)*(m+1)/2;
        double fact=(two*m+one)*x;
        for(int l=m+2; l<=lmax; l++) {
          double pmm=(l+m-1)*plm[ll-1];
          fact=fact+two*x;
          ll=ll+l-1;
          plm[ll+l-1]=( fact*plm[ll-1] - pmm )/double(l-m);
        }

      } //end for m
    } //end else
  }
}
//
// for the normalized associated Legendre functions \bar{P}_{lm}
// such that the spherical harmonics are:
// Y_lm(\theta, \phi) = \bar{P}_{lm}(\cos \theta) e^{i m \phi}
// use the recursion relation:
// P_{00}(x) = \sqrt{1/2\pi}
//
// i.e \bar{P}_{lm}=\sqrt{\frac{(2l+1)(l-m)!}{4\pi(l+m)!}} P_{lm}
//
// Plm is a 1-d array that will contain all the values of P_{lm}(x) from P_{00} to P_{l_{max} l_{max}}
// the index into this array is Plm[l*(l+1)/2 + m]
//

__device__ __inline__ int plmIdxDev(int l, int m)
{ return l*(l+1)/2+m; }

__device__ __inline__
void associatedLegendreFunctionNormalizedDevice(double x, int lmax, double *Plm)
{
  const double pi = std::acos(-1.0);
  // y = \sqrt{1-x^2}
  double y = std::sqrt(1.0-x*x);
  // initialize the first entry
  // Plm[0]=std::sqrt(R(1)/(R(2)*pi));
  if ( threadIdx.x == 0 ) Plm[0]=std::sqrt(1.0/(4.0*pi));
  __syncthreads();

  if(lmax<1) return;

  //FIXME - Parallel versions of these currently broken. Running in parallel
  //        (incorrectly) gave an O(35%) speed-up on the kernel. -JLarkin 9/2017
  //for(int m=threadIdx.x+1; m<=lmax; m+=blockDim.x) //Parallel
  if ( threadIdx.x == 0 )
  {
  for(int m=1; m<=lmax; m++) //Serial
  {
    // \bar{P}_{mm} = - \sqrt{\frac{2m+1}{2m}} y \bar{P}_{m-1, m-1}
    Plm[plmIdxDev(m,m)] = - std::sqrt(((double)(2*m+1)/(double)(2*m))) * y * Plm[plmIdxDev(m-1,m-1)];

    // \bar{P}_{mm-1} = \sqrt{2 m + 1} x \bar{P}_{m-1, m-1}
    Plm[plmIdxDev(m,m-1)] = std::sqrt((double)(2*m+1)) * x * Plm[plmIdxDev(m-1,m-1)]; 
  }

  //for(int m=threadIdx.x; m<lmax; m+=blockDim.x) //Parallel
  for(int m=0; m<lmax; m++) //Serial
  {
    for(int l=m+2; l<=lmax; l++) //Serial
    {
      // \bar{P}_{lm} = a_{lm} (x \bar{P}_{l-1. m} - b_{lm} \bar{P}_{l-2, m})
      // a_{lm} = \sqrt{\frac{(4 l^2 - 1)(l^2 - m^2)}}
      // b_{lm} = \sqrt{\frac{(l -1)^2 - m^2}{4 (l-1)^2 -1}}
      double a_lm = std::sqrt((double)(4*l*l-1)/(double)(l*l - m*m));
      double b_lm = std::sqrt((double)(((l-1)*(l-1) - m*m)/(double)(4*(l-1)*(l-1)-1)));
      Plm[plmIdxDev(l,m)] = a_lm * (x * Plm[plmIdxDev(l-1,m)] - b_lm * Plm[plmIdxDev(l-2,m)]);
    }
  }
  }
  __syncthreads(); // Wait for serial section to complete
}
