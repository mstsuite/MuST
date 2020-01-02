#include <cstdlib>
#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "Array3d.hpp"


// #define INLINE_PLGLMAX

extern "C"
{
  void plglmax_(int *lend, Real*costheta, Real*plm);
}
// void plglmax_new(int lmax, Real x, Real*plm);

// #include "plglmax_new.cpp"
#include "associatedLegendreFunction.hpp"

inline void makegij_new(int lmaxi, int kkri, int lmaxj, int kkrj,
                 int lmax, int kkrsz, int ndlj, int ndlm,
                 Complex prel, Real *rij,Real *sinmp, Real *cosmp,
                 Real *clm, Real *plm, Array3d<Real> &cgnt,
                 int *lofk, int *mofk,
                 Complex *ilp1, Matrix<Complex> &illp,
                 Complex *hfn, Complex *dlm,
                 Complex *gij, Real pi4)
{
  const Complex sqrtm1=std::complex<Real>(0.0,1.0);
  const Real ptol=1.0e-6;
/*
c     *****************************************************************
c    
c
c      ij         l+1                                m ->  *
c     D  [E]  = -i    * sqrt(E) * h [sqrt(E)*R  ] * Y [R  ]
c      L                           l          ij     l  ij
c    
c
c                 -l+1                               m  ->  *
c             = -i    * sqrt(E) * h [sqrt(E)*R  ] * Y [-R  ]
c                                  l          ij     l   ij
c    
c The extra factor (-1)^m removed by xgz July 97 from plm (plglmax_c.f)
c to make Ylm consistent with Numerical Recipes.
c
c      m       [ (2*l+1)*(l-|m|)!]  m
c     Y  = sqrt[-----------------]*P [cos(theta)]*exp(i*m*phi)
c      l       [   4*pi*(l+|m|)! ]  l
c
c     for m>=0 
c
c      m                  m  -m           *
c     Y [theta,phi] = (-1)  Y  [theta,phi]      for m<0
c      l                     l
c
c     ->    ->   ->
c     R   = R  - R  ==> [theta,phi]
c      ij    j    i
c
c     cos(theta)=Rij(3)/sqrt(Rij(1)**2+Rij(2)**2+Rij(3)**2)
c
c     cos(phi)  =Rij(1)/sqrt(Rij(1)**2+Rij(2)**2)
c
c                m     m
c     Indexing: P  :  P(l*(l+1)/2+m+1)  only +m are calculated
c                l     l
c
c                m     m
c     Indexing: C  :  C(l*(l+1)/2+m+1)  only +m are calculated
c                l     l
c
c                m     m
c     Indexing: D  :  D(l*(l+1)+m+1)    all   m are calculated
c                l     l
c                    
c     Now for the real space structure contant we have :
c                    
c      ij             l-l'         L       ij
c     G   (E) = 4*pi*i    * SUM { C     * D  (E) }
c      L,L'                  L"    L',L"   L"
c
c     *****************************************************************
*/

        int lend=lmaxi+lmaxj;
/*
c     calculate the hankel function.[dangerous code if z is close to 0]
c     hankel function is hfn*exp(i*z)/z
*/
      Real rmag=std::sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
//      if(rmag.lt.ptol) stop 'Error in makegij: rmag = 0.'
      if(std::abs(prel)==0.0)
      {
// if prel is zero then we calculate Gij for multipole Coulomb interaction.
	hfn[0]=1.0/rmag;
	for(int l=1; l<=lend; l++)
        {
	  hfn[l]=sqrtm1*(2.0*l-1)*hfn[l-1]/rmag;
        }
      } else {
        Complex z=prel*rmag;
        hfn[0]=-sqrtm1;
        hfn[1]=-1.0-sqrtm1/z;
        for(int l=1; l<lend; l++)
        {
          hfn[l+1]=(2.0*l+1)*hfn[l]/z - hfn[l-1];
        }
/*
c             l+1
c     hfn = -i   *h (k*R  )*sqrt(E)
c                  l    ij
*/
        z=std::exp(sqrtm1*z)/rmag;
        for(int l=0; l<=lend;l++)
        {
          hfn[l]=-hfn[l]*z*ilp1[l];     
        }
      }

//     calculate p(l,m)'s...............................................
      Real costheta=rij[2]/rmag;
//      plglmax_(&lend,&costheta,plm);
#ifdef INLINE_PLGLMAX
{
  const Real pi = std::acos(-Real(1));
  // y = \sqrt{1-x^2}
  Real y = std::sqrt(Real(1)-costheta*costheta);
  // initialize the first entry
  // normalized such that clm[i]=1.0
  plm[0]=std::sqrt(Real(1)/(Real(4)*pi));

  if(lend>0)
  {
    for(int m=1; m<=lend; m++)
    {
      // \bar{P}_{mm} = - \sqrt{\frac{2m+1}{2m}} y \bar{P}_{m-1, m-1}
      plm[plmIdx(m,m)] = - std::sqrt(Real(2*m+1)/Real(2*m)) * y * plm[plmIdx(m-1,m-1)];
      // \bar{P}_{mm-1} = \sqrt{2 m + 1} x \bar{P}_{m-1, m-1}
      plm[plmIdx(m,m-1)] = std::sqrt(Real(2*m+1)) * costheta * plm[plmIdx(m-1,m-1)]; 
    }

    for(int m=0; m<lend; m++)
    {
      for(int l=m+2; l<=lend; l++)
      {
        // \bar{P}_{lm} = a_{lm} (x \bar{P}_{l-1. m} - b_{lm} \bar{P}_{l-2, m})
        // a_{lm} = \sqrt{\frac{(4 l^2 - 1)(l^2 - m^2)}}
        // b_{lm} = \sqrt{\frac{(l -1)^2 - m^2}{4 (l-1)^2 -1}}
        Real a_lm = std::sqrt(Real(4*l*l-1)/Real(l*l - m*m));
        Real b_lm = std::sqrt(Real((l-1)*(l-1) - m*m)/Real(4*(l-1)*(l-1)-1));
        plm[plmIdx(l,m)] = a_lm * (costheta * plm[plmIdx(l-1,m)] - b_lm * plm[plmIdx(l-2,m)]);
      }
    }
  }
}

#else
//      plglmax_new(lend,costheta,plm);
      associatedLegendreFunctionNormalized<Real>(costheta, lend, plm);
#endif
//     multiply be the normalization constant...........................
      int ndlm_local=(lend+1)*(lend+2)/2;
      if(ndlm_local>ndlm)
      {
	printf("MAKEGIJ:: ndlm incorrect!\n");
        printf("ndlm=%d\nndlm_local=%d\n",ndlm,ndlm_local);
        exit(1);
      }
//      for(int j=0;j<ndlm_local;j++)
//        plm[j]=clm[j]*plm[j];
//     =================================================================
//     calculate cos(phi) and sin(phi) .................................
      Real pmag=std::sqrt(rij[0]*rij[0]+rij[1]*rij[1]);
      cosmp[0]=1.0;
      sinmp[0]=0.0;
      if(pmag>ptol)
      {
        cosmp[1]=rij[0]/pmag;
        sinmp[1]=rij[1]/pmag;
      } else {
        cosmp[1]=0.0;
        sinmp[1]=0.0;
      }
      for(int m=2; m<=lend; m++)
      {
        cosmp[m]=cosmp[m-1]*cosmp[1]-sinmp[m-1]*sinmp[1];
        sinmp[m]=sinmp[m-1]*cosmp[1]+cosmp[m-1]*sinmp[1];
      }
  
      int j=0;
      for(int l=0; l<=lend; l++)
      {
        int ll=l*(l+1);
        j=ll;
        ll=ll/2;
        Real m1m=1.0;
        dlm[j]= hfn[l]*plm[ll];
        for(int m=1; m<=l; m++)
        {
          m1m=-m1m;
          Complex fac=plm[ll+m]*std::complex<Real>(cosmp[m],sinmp[m]);
          dlm[j-m]= hfn[l]*m1m*fac;
          dlm[j+m]= hfn[l]*std::conj(fac);
        }
      }

//     ================================================================
//     calculate g(R_ij)...............................................
      for(int i=0; i<kkri*kkrj; i++) gij[i]=0.0;

//     loop over l1,m1............................................
      for(int lm1=0; lm1<kkrj; lm1++)
      {
        int l1=lofk[lm1];
        int m1=mofk[lm1];

//        loop over l2,m2..............................................
        for(int lm2=0; lm2<kkri; lm2++)
        {
          int l2=lofk[lm2];
          int m2=mofk[lm2];
/*
c           ==========================================================
c                            l2-l1
c           illp(lm2,lm1) = i
c
c           perform sum over l3 with gaunt # ......................
c           ==========================================================
*/
          int m3=m2-m1;
          int llow=std::max(std::abs(m3),std::abs(l1-l2));
          if(std::abs(prel)==0.0) llow=l1+l2;
          for(int l3=l1+l2; l3>=llow;l3-=2)
          {
            int j=l3*(l3+1)+m3;
            gij[lm2+lm1*kkri] = gij[lm2+lm1*kkri]+cgnt(l3/2,lm1,lm2)*dlm[j];
          }
          gij[lm2+lm1*kkri]=pi4*illp(lm2,lm1)*gij[lm2+lm1*kkri];
        }
      }
}
