#ifndef LSMS_COEFICIENTS_HPP
#define LSMS_COEFICIENTS_HPP

#include "Complex.hpp"
#include <vector>
#include <cmath>
#include "Array3d.hpp"
#include "Indices.hpp"

class SphericalHarmonicsCoeficients {
public:
  int lmax;
  std::vector<Real> clm;
  void init(int _lmax)
  {
    Real pi=2.0*std::asin(1.0);
    Real fpi=4.0*pi;
    lmax=_lmax;
    clm.resize((lmax+1)*(lmax+2)/2);
// coefficients for normalized associated Legendre functions (P_00 = \sqrt{1/4pi})
    for(int lm=0; lm<(lmax+1)*(lmax+2)/2; lm++)
      clm[lm]=1.0;
      // clm[lm]=std::sqrt(0.5);
// coefficients for standard form of the assocoiated Legendere functions
/*
    clm[0]=std::sqrt(1.0/fpi);
    for(int l=1; l<=lmax; l++)
    {
      Real xfact = std::sqrt(Real(2*l+1)/fpi);
      for(int m=0; m<=l; m++)
      {
        clm[(l*(l+1))/2 + m]=1.0;
        for(int i=1; i<=2*m; i++)
          clm[(l*(l+1))/2 + m]*=Real(l-m+i);
        clm[(l*(l+1))/2 + m]=xfact*std::sqrt(1.0/clm[(l*(l+1))/2 + m]);
      }
    }
*/
  }
};

extern SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;

extern "C"
{
void cgaunt_(int * lmax, Real *clm,
             Real *plmg, Real *tg, Real *wg,
             Real *cgnt,
             int *lofk, int *mofk,
             int *iprint, char *istop);
}

void calculateGauntCoeficients(int lmax, Array3d<Real> &cgnt, AngularMomentumIndices &a);

class GauntCoeficients {
public:
  int lmax;
  Array3d<Real> cgnt;
  
  void init(LSMSSystemParameters &lsms, AngularMomentumIndices &a, SphericalHarmonicsCoeficients &s)
  {
    const bool useNewGaunt = false;
    const bool testNewGaunt = true;

    lmax=a.lmax/2;
    // lmax=lsms.maxlmax;

    if(useNewGaunt)
    {
      calculateGauntCoeficients(lmax, cgnt, a);
      if(testNewGaunt)
      {
        const Real tol=1.0e-12;
        std::vector<Real> tg,wg;
        Matrix<Real> plmg;
        Array3d<Real> cgntTest;
        tg.resize(2*(2*lmax+1)); wg.resize(2*(2*lmax+1));
        cgntTest.resize(lmax+1,(lmax+1)*(lmax+1),(lmax+1)*(lmax+1));
        plmg.resize(((2*lmax+1)*(2*lmax+2))/2,2*lmax+1);

        cgaunt_(&lmax,&s.clm[0],&plmg(0,0),&tg[0],&wg[0],
                &cgntTest(0,0,0),
                &a.lofk[0],&a.mofk[0],
                &lsms.global.iprint,lsms.global.istop);

        for(int l1=0; l1<lmax+1; l1++)
          for(int j2=0; j2<(lmax+1)*(lmax+1); j2++)
            for(int j3=0; j3<(lmax+1)*(lmax+1); j3++)
            {
              if(std::abs(cgnt(l1,j2,j3)-cgntTest(l1,j2,j3))>tol)
              {
                printf("Difference in Gaunt Coeficients: %d %d %d : new=%.12lf old=%.12lf\n",
                       l1,j2,j3,cgnt(l1,j2,j3),cgntTest(l1,j2,j3));
                exit(1);
              }
            }
      }
    } else {
      std::vector<Real> tg,wg;
      Matrix<Real> plmg;
      tg.resize(2*(2*lmax+1)); wg.resize(2*(2*lmax+1));
      cgnt.resize(lmax+1,(lmax+1)*(lmax+1),(lmax+1)*(lmax+1));
      plmg.resize(((2*lmax+1)*(2*lmax+2))/2,2*lmax+1);

      cgaunt_(&lmax,&s.clm[0],&plmg(0,0),&tg[0],&wg[0],
              &cgnt(0,0,0),
              &a.lofk[0],&a.mofk[0],
              &lsms.global.iprint,lsms.global.istop);
    }
  }
};

extern GauntCoeficients gauntCoeficients;

extern "C"
{
void ifacts_(int *lmax,Complex *illp,Complex *ilp1,int *iprint, char *istop, int istop_len);
}

class IFactors {
public:
  Matrix<Complex> illp;
  std::vector<Complex> ilp1;
  void init(LSMSSystemParameters &lsms, int _lmax)
  {
    ilp1.resize(2*_lmax+1);
    illp.resize((_lmax+1)*(_lmax+1),(_lmax+1)*(_lmax+1));
    ifacts_(&_lmax,&illp(0,0),&ilp1[0],&lsms.global.iprint,lsms.global.istop,32);
  }
};

extern IFactors iFactors;

#endif
