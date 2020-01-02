#ifndef BKKRM_TEST_STRUCTURES_HPP
#define BKKRM_TEST_STRUCTURES_HPP

#include "Misc/Indices.hpp"
/*
// from Misc/Indices.hpp
class AngularMomentumIndices {
public:
  int lmax,ndlj,ndlm;
  std::vector<int> lofk,mofk,lofj,mofj;

  void init(int _lmax)
  {
    lmax=_lmax;
    ndlj=(lmax+1)*(lmax+1);
    ndlm=((lmax+1)*(lmax+2))/2;
    mofk.resize(ndlj); lofk.resize(ndlj);
    mofj.resize(ndlm); lofj.resize(ndlm);
    int j=0;
    int k=0;
    for(int l=0; l<=lmax; l++)
    {
      for(int m=0; m<=l; m++)
      {
        lofj[j]=l;
        mofj[j++]=m;
      }
      for(int m=-l; m<=l; m++)
      {
        
        lofk[k]=l;
        mofk[k]=m;
        // printf("k,l,m: %d %d %d, lofk,mofk:%d %d\n",k,l,m,lofk[k],mofk[k]);
        k++;
      }
    }
    // printf("k=%d\n",k);
  }
};
*/

// only the needed parts from Main/SystemParameters.hpp:
class LSMSGlobals {
public:
  void setIstop(const char *c){strncpy(istop,c,32); for(int i=strlen(c); i<32; i++) istop[i]=' ';}
  bool checkIstop(const char *c){return (strncmp(istop,c,32)==0);}
// ...
  int iprint;
// ...
  char istop[32];
};

class LSMSSystemParameters {
public:
// ...
  int nrel_rel;
  int n_spin_cant;
// ...
  int maxlmax;
  LSMSGlobals global;
  AngularMomentumIndices angularMomentumIndices;
// ...
};

class LocalTypeInfo {
public:
// ...
  int lDimTmatStore;
  Matrix<Complex> tmatStore;
// ...
};


// only the needed parts from SingleSite/AtomData.hpp:
class AtomData {
public:
// ...
// Local Interaction Zone
  int numLIZ;
  std::vector<int> LIZGlobalIdx, LIZStoreIdx, LIZlmax;
  int nrmat; // sum (LIZlmax+1)^2
  std::vector<Real> LIZDist;
  Matrix<Real> LIZPos;
// ...
// General Data
  int lmax,kkrsz;
// ...
};

#endif
