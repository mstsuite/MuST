#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "TestStructures.hpp"
#include "Misc/Coeficients.hpp"

#ifdef CRAYPAT
#include <pat_api.h>
#endif

extern "C"
{
  void write_kkrmat_(Complex *a,int *n,int *lda,Complex *e);

  void makegij_(int *lmaxi,int *kkri,int *lmaxj,int *kkrj,
                int *lmax,int *kkrsz,int *ndlj,int *ndlm,
                Complex *prel,double *rij,double *sinmp,double *cosmp,
                double *clm,double *plm,double *cgnt,int *lmax_cg,int *lofk,int *mofk,
                Complex *ilp1,Complex *illp,
                Complex *hfn,Complex *dlm,Complex *gij,
                double *pi4,int *iprint,char *istop,int len_sitop);

  void setgij_(Complex *gij,Complex *bgij,int *kkr1,int *kkr1_ns,int *kkr2,int *kkr2_ns,
               int *n_spin_cant,int *nrel_rel,Complex *psq,Complex *energy);
};

// #define SYNTHETIC_MATRIX
// #define WRITE_GIJ

#define USE_ALLOCA

void buildKKRMatrix_orig(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, Complex energy, Complex prel, Matrix<Complex> &m)
{
#ifdef CRAYPAT
PAT_region_begin(1,"buildKKRMatrix init");
#endif
  Real rij[3];

  int lmax=lsms.maxlmax;
  int kkrsz=(lmax+1)*(lmax+1);
  int kkrsz_ns=kkrsz*lsms.n_spin_cant;
  int nrmat_ns=lsms.n_spin_cant*atom.nrmat;
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
#ifdef CRAYPAT
PAT_region_end(1);
#endif
/*
  // print first t matrix
  if(lsms.global.iprint>=1)
  {
    printf("first t Matrix:\n");
    for(int i=0; i<kkrsz_ns; i++)
      for(int j=0; j<kkrsz_ns; j++)
      {
        printf("%d %d %.8le %.8le\n",i,j,std::real(local.tmatStore(i+j*kkrsz_ns,0)),std::imag(local.tmatStore(i+j*kkrsz_ns,0)));
      }
// write LIZ positions for first atom
    FILE *of=fopen("LIZ_pos.h","w");
    fprintf(of,"atom.numLIZ=%d;\n",atom.numLIZ);
    fprintf(of,"atom.LIZPos.resize(3,atom.numLIZ);\n");
    for(int i=0; i<atom.numLIZ; i++)
    {
      fprintf(of,"atom.LIZPos(0,%d)=%lg;\n",i,atom.LIZPos(0,i));
      fprintf(of,"atom.LIZPos(1,%d)=%lg;\n",i,atom.LIZPos(1,i));
      fprintf(of,"atom.LIZPos(2,%d)=%lg;\n",i,atom.LIZPos(2,i));
    }
    fclose(of);
    exit(0);
  }
*/

#ifdef WRITE_GIJ
  Matrix<Complex> Gij_full(nrmat_ns,nrmat_ns);
  for(int i=0; i<nrmat_ns*nrmat_ns; i++) Gij_full[i]=0.0;
#endif

  nrst=0;
  for(int ir1=0; ir1<atom.numLIZ; ir1++)
  {
    int kkr1=(atom.LIZlmax[ir1]+1)*(atom.LIZlmax[ir1]+1);
    int kkr1_ns=kkr1*lsms.n_spin_cant;

    ncst=0;
// build local t_mat:
#ifdef CRAYPAT
PAT_region_begin(2,"build local t_mat");
#endif
    int im=0;
    for(int js=0; js<lsms.n_spin_cant; js++)
    {
      int jsm = kkrsz*kkrsz_ns*js;
      for(int j=0; j<kkr1; j++)
      {
        for(int is=0; is<lsms.n_spin_cant; is++)
        {
          int jm=jsm+kkrsz_ns*j+kkrsz*is;
          int one=1;
          BLAS::zcopy_(&kkr1,&local.tmatStore(jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
          im+=kkr1;
        }
      }
    }
#ifdef CRAYPAT
PAT_region_end(2);
#endif
/*
    if(ir1==1)
    {
      printf("atom.LIZlmax[ir1]=%d\n",atom.LIZlmax[ir1]);
      for(int i=0; i<atom.kkrsz*atom.kkrsz*4; i++)
      {
        Complex t=local.tmatStore(i,atom.LIZStoreIdx[ir1]);
        if(std::norm(t)!=0.0)
        {
          printf("tmat_n[%d]=(%lf, %lf)\n",i,std::real(tmat_n[i]),std::imag(tmat_n[i]));
          printf("tmatStore[%d,ir]=(%lf, %lf)\n",i,std::real(t),std::imag(t));
        }
      }
    }
*/
//
#ifdef CRAYPAT
PAT_region_begin(3,"buildKKRMatrix ir2 loop");
#endif
    for(int ir2=0; ir2<atom.numLIZ; ir2++)
    {
      int kkr2=(atom.LIZlmax[ir2]+1)*(atom.LIZlmax[ir2]+1);
      int kkr2_ns=kkr2*lsms.n_spin_cant;
      if(ir1!=ir2)
      {
        rij[0]=atom.LIZPos(0,ir1)-atom.LIZPos(0,ir2);
        rij[1]=atom.LIZPos(1,ir1)-atom.LIZPos(1,ir2);
        rij[2]=atom.LIZPos(2,ir1)-atom.LIZPos(2,ir2);

        makegij_(&atom.LIZlmax[ir1],&kkr1,&atom.LIZlmax[ir2],&kkr2,
                 &lsms.maxlmax,&kkrsz,&lsms.angularMomentumIndices.ndlj,&lsms.angularMomentumIndices.ndlm,
                 &prel,&rij[0],sinmp,cosmp,
                 &sphericalHarmonicsCoeficients.clm[0],plm,
                 &gauntCoeficients.cgnt(0,0,0),&gauntCoeficients.lmax,
                 &lsms.angularMomentumIndices.lofk[0],&lsms.angularMomentumIndices.mofk[0],
                 &iFactors.ilp1[0],&iFactors.illp(0,0),
                 hfn,dlm,gij,
                 &pi4,&lsms.global.iprint,lsms.global.istop,32);
        Complex psq=prel*prel;
        setgij_(gij,bgij,&kkr1,&kkr1_ns,&kkr2,&kkr2_ns,
                &lsms.n_spin_cant,&lsms.nrel_rel,&psq,&energy);
#ifdef WRITE_GIJ
        for(int ii=nrst; ii<nrst+kkr1_ns; ii++)
          for(int jj=ncst; jj<ncst+kkr2_ns; jj++)
            Gij_full(ii,jj)=bgij[ii-nrst+(jj-ncst)*kkr2_ns];
#endif

        BLAS::zgemm_("n","n",&kkr1_ns,&kkr2_ns,&kkr1_ns,&cmone,
                    tmat_n,&kkr1_ns,bgij,&kkr1_ns,&czero,
                    &m(nrst,ncst),&nrmat_ns);
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
#ifdef CRAYPAT
PAT_region_end(3);
#endif
    nrst+=kkr1_ns;
  }

#ifdef SYNTHETIC_MATRIX
  std::cout<<"WARNING: USING SYNTHETIC MATRIX!!!!\n";
  for(int i=0; i<nrmat_ns; i++)
    for(int j=0; j<nrmat_ns; j++)
    {
      m(i,j)=Complex((11.0*(i+1)-3.0*(j+1))/Real(i+j+2),(5.0*(i+1)-2.0*(j+1))/Real(i+j+2));
    }
#endif

/*
#ifdef WRITE_GIJ
  write_kkrmat_(&Gij_full(0,0),&nrmat_ns,&nrmat_ns, &energy);
#else
  write_kkrmat_(&m(0,0),&nrmat_ns,&nrmat_ns, &energy);
#endif
  exit(0);
*/

/*
  if(lsms.global.checkIstop("buildKKRMatrix"))
  {
    if(lsms.global.iprint>=0)
    {
      FILE *f1=fopen("kkrmat.out","w");
      FILE *f2=fopen("kkrmat.pattern","w");
      for(int i=0; i<nrmat_ns; i++)
      {
        fprintf(f2,"%5d ",i);
        for(int j=0; j<nrmat_ns; j++)
        {
          fprintf(f1,"%5d %5d %lf %lf\n",i,j,std::real(m(i,j)),std::imag(m(i,j)));
          if(std::abs(m(i,j))<0.000001) fprintf(f2,".");
          else fprintf(f2,"x");
        }
        fprintf(f2,"\n");
      }
      fclose(f1); fclose(f2);
    }
// calculate matrix norm and condition number
    double work[2*nrmat_ns];
    int ipiv[nrmat_ns];
    int info;
    double norm=zlange_("1",&nrmat_ns,&nrmat_ns,&m(0,0),&nrmat_ns,work);
    printf("Matrix 1-norm of the kkr matrix=%lf\n",norm);
    zgetrf_(&nrmat_ns, &nrmat_ns, &m(0,0),&nrmat_ns,ipiv,&info);
    printf("ZGETRF on kkr-matrix returned info=%d\n",info);
    if(info==0)
    {
      zgetri_(&nrmat_ns,&m(0,0),&nrmat_ns,ipiv,(Complex *)work,&nrmat_ns,&info);
      printf("ZGETRI returned info=%d\n",info);
      double norm_i=zlange_("1",&nrmat_ns,&nrmat_ns,&m(0,0),&nrmat_ns,work);
      printf("Matrix 1-norm of the kkr matrix inverse=%lf\n",norm_i);
      printf("1-norm condition number of the kkr-matrix=%lf\n",norm*norm_i);
    }
    // exit(0);
  }
*/
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
