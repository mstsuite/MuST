#include "Complex.hpp"
#include <vector>
#include <cmath>
#include <cstdlib>
#include "Array3d.hpp"
#include "Indices.hpp"

#ifdef USE_GMP
#include "Gaunt_gmp.hpp"
#else
#include "Gaunt.hpp"
#endif

/*
c
c     ****************************************************************
c     generate:
c
c               4pi  _  m1 _     m2* _     m3 _
c        clll = int do Y  (o) * Y   (o) * Y  (o)
c                0      l1       l2        l3
c
c                L3
c             = C
c                L1,L2
c
c     L1, L2: l = 0,1,..,lmax;  L3: l=0,1,...,2*lmax
c     ****************************************************************
*/ 

void calculateGauntCoeficients(int lmax, Array3d<Real> &cgnt, AngularMomentumIndices &a)
{
#ifdef USE_GMP
  mpf_set_default_prec(128);
#endif
  cgnt.resize(lmax+1,(lmax+1)*(lmax+1),(lmax+1)*(lmax+1));
  for(int j2=0; j2<(lmax+1)*(lmax+1); j2++)
  {
    int l2=a.lofk[j2];
    int m2=a.mofk[j2];
    for(int j1=0; j1<(lmax+1)*(lmax+1); j1++)
    {
      int l1=a.lofk[j1];
      int m1=a.mofk[j1];
      int m3=m2-m1;
      int lb=std::max(abs(m3),abs(l1-l2));
      for(int l3=l1+l2; l3>=lb; l3-=2)
      {
        cgnt(l3/2,j1,j2)=gaunt<Real,long int>(l2,l1,l3,m2,m1,m3);
      }
    }
  }
}

