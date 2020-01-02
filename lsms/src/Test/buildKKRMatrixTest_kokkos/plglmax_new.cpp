#include "Real.hpp"
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
#define LMAX_GE_2

inline void plglmax_new(int lmax, Real x, Real*plm)
{
  const Real tol=1.0e-12;
// removed checks from original plglmax.f

  if((1.0-std::abs(x))<tol)
  {
    for(int i=0; i<(lmax+1)*(lmax+2)/2;i++) plm[i]=0.0;
    if(x<0.0)
    {
      for(int l=0; l<=lmax; l++)
      {
        int i=l*(l+1)/2;
        plm[i]=1.0-2.0*(l%2);
      }
    } else {
      for(int l=0; l<=lmax; l++)
      {
        int i=l*(l+1)/2;
        plm[i]=1.0;
      }
    }
  } else {
    Real fact;
    int ll;
/*
c     ================================================================
c     begin calculation of p(l,m)'s...................................
c     ================================================================
c        special case lmax=0..........................................
*/
      plm[0]=1.0;
#ifndef LMAX_GE_2
      if(lmax==0) return;
#endif

//      else
//        =============================================================
//        minus sign added to be consistant with Numerical Recipes
//        which has (-1)^m factor in plm :  July 97  by xgz.............
//        =============================================================
      Real somx2=-std::sqrt((1.0-x)*(1.0+x));
//         if(lmax.eq.1) then
//           ==========================================================
//           special case lmax=1.......................................
      plm[1]=x;
      plm[2]=somx2;
#ifndef LMAX_GE_2
      if(lmax==1) return;
#endif
// m==0 special case
      fact=x;
      ll=0;
      for(int l=2; l<=lmax; l++)
      {
        Real pmm=(l-1.0)*plm[ll];
        fact+=2.0*x;
        ll+=l-1;
        plm[ll+l]=(fact*plm[ll]-pmm)/Real(l);
      }

      for(int m=1; m<lmax; m++)
      {
/*
c              =======================================================
c                                       m       m
c              calculate the first two P   and P
c                                       m       m+1
c              =======================================================
*/
        Real pmm=somx2;
        fact=1.0;
        for(int i=2; i<=m; i++)
        {
          fact+=2.0;
          pmm*=fact*somx2;
        }
        int mm=m*(m+1)/2+m;
        plm[mm]=pmm;
        if( mm+m+1<(lmax+1)*(lmax+2)/2 )
          plm[mm+m+1]=x*(2.0*m+1.0)*pmm;
/*
c              =======================================================
c                                  m        m
c              calculate the rest P     to P
c                                  m+2      lmax
c              =======================================================
*/
        ll=(m+1)*m/2+m;
        fact=(2.0*m+1.0)*x;
        for(int l=m+2; l<=lmax; l++)
        {
          pmm=(l+m-1.0)*plm[ll];
          fact+=2.0*x;
          ll=ll+l-1;
          plm[ll+l]=( fact*plm[ll] - pmm )/Real(l-m);
        }
      }
// m==lmax
      {
        int m=lmax;
        Real pmm=somx2;
        fact=1.0;
        for(int i=2; i<=m; i++)
        {
          fact+=2.0;
          pmm*=fact*somx2;
        }
        int mm=m*(m+1)/2+m;
        plm[mm]=pmm;
      }
  }
}
