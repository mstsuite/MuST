#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "TestStructures.hpp"
#include "Misc/Coeficients.hpp"

#include "associatedLegendreFunction.hpp"

// #define USE_ALLOCA
// #define INLINE_PLGLMAX
 #define INLINE_MAKEGIJ

void plglmax_new(int lmax, Real x, Real*plm);

#include "makegij_new.cpp" 

void buildKKRMatrix_nrel_ns2(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom,
                             Complex energy, Complex prel, Matrix<Complex> &m)
{
// make sure that lsms.n_spin_cant == 2!
  if(lsms.n_spin_cant!=2)
  {
    printf("lsms.n_spin_cant!=2 in buildKKRMatrix_nrel_ns2\n");
    exit(1);
  }

  Real rij[3];

  int lmax=lsms.maxlmax;
  int kkrsz=(lmax+1)*(lmax+1);
  int kkrsz_ns=2*kkrsz; // lsms.n_spin_cant
  int nrmat_ns=2*atom.nrmat;  // lsms.n_spin_cant
  int nrst,ncst;

#ifdef USE_ALLOCA
  Complex *gij=(Complex *)alloca(kkrsz*kkrsz*sizeof(Complex));
  Real *sinmp = (Real *)alloca(sizeof(Real)*(2*lmax+1));
  Real *cosmp = (Real *)alloca(sizeof(Real)*(2*lmax+1));
  Real *plm = (Real *)alloca(sizeof(Real)*lsms.angularMomentumIndices.ndlm);
  Complex *hfn = (Complex *)alloca(sizeof(Complex)*(2*lmax+1));
  Complex *dlm = (Complex *)alloca(sizeof(Complex)*lsms.angularMomentumIndices.ndlj);
  Complex *bgij = (Complex *)alloca(sizeof(Complex)*4*kkrsz*kkrsz);
  Complex *tmat_n = (Complex *)alloca(sizeof(Complex)*atom.kkrsz*atom.kkrsz*4);
#else
  Complex *gij = new Complex[kkrsz*kkrsz];
  Real *sinmp = new Real[2*lmax+1];
  Real *cosmp = new Real[2*lmax+1];
  Real *plm = new Real[lsms.angularMomentumIndices.ndlm];
  Complex *hfn = new Complex[2*lmax+1];
  Complex *dlm = new Complex[lsms.angularMomentumIndices.ndlj];
  Complex *bgij = new Complex[4*kkrsz*kkrsz];
  Complex *tmat_n = new Complex[atom.kkrsz*atom.kkrsz*4];
#endif

  const Complex cmone=-1.0;
  const Complex czero=0.0;

  Real pi4=4.0*2.0*std::asin(1.0);

  for(int i=0; i<nrmat_ns*nrmat_ns; i++) m[i]=0.0;
  for(int i=0; i<nrmat_ns; i++) m(i,i)=1.0;

  nrst=0;
  for(int ir1=0; ir1<atom.numLIZ; ir1++)
  {
    int kkr1=(atom.LIZlmax[ir1]+1)*(atom.LIZlmax[ir1]+1);
    int kkr1_ns=2*kkr1; // lsms.n_spin_cant;

    ncst=0;
// build local t_mat:
// the dimension is (2*kkr1,2*kkr1)
// if kkrsz == kkr1 then tmat_n=tmatStore
    if(kkrsz==kkr1)
    {
      // clbas_zcopy(kkrsz*kkrsz*4,&local.tmatStore(0,atom.LIZStoreIdx[ir1]),1,&tmat_n[0],1);
      for(int i=0; i<kkrsz*kkrsz*4; i++) tmat_n[i]=local.tmatStore(i,atom.LIZStoreIdx[ir1]);
    } else {
      int im=0;
      for(int js=0; js<lsms.n_spin_cant; js++)
      {
        int jsm = kkrsz*kkrsz_ns*js;
        for(int j=0; j<kkr1; j++)
        {
          for(int is=0; is<lsms.n_spin_cant; is++)
          {
            int jm=jsm+kkrsz_ns*j+kkrsz*is;
            // BLAS::zcopy_(&kkr1,&local.tmatStore(jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
            for(int k=0; k<kkr1; k++) tmat_n[im+k]=local.tmatStore(jm+k,atom.LIZStoreIdx[ir1]);
            im+=kkr1;
          }
        }
      }
    }

    for(int ir2=0; ir2<atom.numLIZ; ir2++)
    {
      int kkr2=(atom.LIZlmax[ir2]+1)*(atom.LIZlmax[ir2]+1);
      int kkr2_ns=2*kkr2; // *lsms.n_spin_cant;
      if(ir1!=ir2)
      {
        rij[0]=atom.LIZPos(0,ir1)-atom.LIZPos(0,ir2);
        rij[1]=atom.LIZPos(1,ir1)-atom.LIZPos(1,ir2);
        rij[2]=atom.LIZPos(2,ir1)-atom.LIZPos(2,ir2);

#ifndef INLINE_MAKEGIJ
        makegij_new(atom.LIZlmax[ir1], kkr1, atom.LIZlmax[ir2], kkr2,
                    lsms.maxlmax, kkrsz, lsms.angularMomentumIndices.ndlj, lsms.angularMomentumIndices.ndlm,
                    prel, rij, sinmp, cosmp,
                    &sphericalHarmonicsCoeficients.clm[0], plm, gauntCoeficients.cgnt,
                    &lsms.angularMomentumIndices.lofk[0], &lsms.angularMomentumIndices.mofk[0],
                    &iFactors.ilp1[0], iFactors.illp,
                    hfn, dlm,
                    gij, pi4);

#else
{
  int lmaxi=atom.LIZlmax[ir1];
  int kkri=kkr1;
  int lmaxj=atom.LIZlmax[ir2];
  int kkrj=kkr2;
  int lmax=lsms.maxlmax;
  int ndlj=lsms.angularMomentumIndices.ndlj;
  int ndlm=lsms.angularMomentumIndices.ndlm;
  Real *clm=&sphericalHarmonicsCoeficients.clm[0];
  int *lofk=&lsms.angularMomentumIndices.lofk[0];
  int *mofk=&lsms.angularMomentumIndices.mofk[0];

  const Complex sqrtm1=std::complex<Real>(0.0,1.0);
  const Real ptol=1.0e-6;

        int lend=lmaxi+lmaxj;

      Real rmag=std::sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);

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
          hfn[l]=-hfn[l]*z*iFactors.ilp1[l];     
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
      // plglmax_new(lend,costheta,plm);
      associatedLegendreFunctionNormalized<Real>(costheta, lend, plm);
#endif

      int ndlm_local=(lend+1)*(lend+2)/2;
      if(ndlm_local>ndlm)
      {
	printf("MAKEGIJ:: ndlm incorrect!\n");
        printf("ndlm=%d\nndlm_local=%d\n",ndlm,ndlm_local);
        exit(1);
      }
//      for(int j=0;j<ndlm_local;j++)
//        plm[j]=clm[j]*plm[j];

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
      for(int i=0; i<kkri*kkrj; i++) gij[i]=0.0;

      for(int lm1=0; lm1<kkrj; lm1++)
      {
        int l1=lofk[lm1];
        int m1=mofk[lm1];

        for(int lm2=0; lm2<kkri; lm2++)
        {
          int l2=lofk[lm2];
          int m2=mofk[lm2];
          int m3=m2-m1;
          int llow=std::max(std::abs(m3),std::abs(l1-l2));
          if(std::abs(prel)==0.0) llow=l1+l2;
          for(int l3=l1+l2; l3>=llow;l3-=2)
          {
            int j=l3*(l3+1)+m3;
            gij[lm2+lm1*kkri] = gij[lm2+lm1*kkri]+ gauntCoeficients.cgnt(l3/2,lm1,lm2)*dlm[j];
          }
          gij[lm2+lm1*kkri]=pi4*iFactors.illp(lm2,lm1)*gij[lm2+lm1*kkri];
        }
      }
}
#endif

        Complex psq=prel*prel;

//        for(int i=0; i<kkr1_ns*kkr2_ns; i++) bgij[i]=0.0;
        for(int j=0; j<kkr2; j++)
        {
          for(int i=0; i<kkr1; i++)
          {
            bgij[i+j*kkr1_ns]=gij[i+j*kkr1];
            bgij[i+kkr1+(j+kkr2)*kkr1_ns]=gij[i+j*kkr1];
            bgij[i+kkr1+j*kkr1_ns]=0.0;
            bgij[i+(j+kkr2)*kkr1_ns]=0.0;
          }
        }
//*
        BLAS::zgemm_("n","n",&kkr1_ns,&kkr2_ns,&kkr1_ns,&cmone,
                    tmat_n,&kkr1_ns,bgij,&kkr1_ns,&czero,
                    &m(nrst,ncst),&nrmat_ns);
//*/
/*
        for(int j=0; j<kkr2_ns; j++)
          for(int i=0; i<kkr1_ns; i++)
          {
            m(nrst+i,ncst+j)=0.0;
            for(int k=0; k<kkr1_ns; k++)
              m(nrst+i,ncst+j)-=tmat_n[i+k*kkr1_ns]*bgij[k+j*kkr1_ns];
          }
*/
      }
      ncst+=kkr2_ns;
    }
    nrst+=kkr1_ns;
  }

#ifndef USE_ALLOCA
  delete [] gij;
  delete [] sinmp;
  delete [] cosmp;
  delete [] plm;
  delete [] hfn;
  delete [] dlm;
  delete [] bgij;
  delete [] tmat_n;
#endif
}
