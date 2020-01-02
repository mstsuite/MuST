#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

//[eric] OpenACC has no support of cabs() and cexp()

void c_kkrmat(int lmax, int ndlm, int ndlj, int nrmat_ns, int akkrsz, int anumLIZ, 
    int *lofk,            // [ndlj]   Misc/Indices.hpp
    int *mofk,            // [ndlj]   Misc/Indices.hpp
    double *LIZPos,       // [3,atom.numLIZ] LIZ_pos.h
    int *LIZlmax,         // [atom.numLIZ] main()
    int *LIZStoreIdx,     // [atom.numLIZ] main()
    double complex *ilp1, // [2*lmax+1]   Misc/Coeficients.hpp
    double complex *illp, // [kkrsz,kkrsz]    Misc/Coeficients.hpp, kkrsz=(lmax+1)^2
    double *cgnt,         // [lmax+1,kkrsz,kkrsz]   Misc/Coeficients.hpp
    double *clm,          // [(2*lmax+1)*(2*lmax+2)/2]
    double complex prel, double complex *m, double complex *tmat0)
{

  int kkrsz=(lmax+1)*(lmax+1);
  int szBlk_tmat = kkrsz*kkrsz*4;
  int kkr1=(LIZlmax[0]+1)*(LIZlmax[0]+1);
  int kkr1_ns=2*kkr1; // lsms.n_spin_cant;
  int kkr2=(LIZlmax[0]+1)*(LIZlmax[0]+1);
  int kkr2_ns=2*kkr2; // lsms.n_spin_cant;

  double sinmp[2*lmax+1];              // 7
  double cosmp[2*lmax+1];              // 7
  double plm[ndlm];                    // 28
  double complex hfn[2*lmax+1];        // 7
  double complex dlm[ndlj*anumLIZ*anumLIZ];  // 49
  double complex gij[kkrsz*kkrsz];     // 256

#define MPdlm(i,ir1,ir2) (i+ir1*ndlj+ir2*ndlj*anumLIZ)

  double complex tmat_n[akkrsz*akkrsz*4*anumLIZ];

  const double complex J=_Complex_I;
  const double ptol=1.0e-6;
  const double pi4=4.0*2.0*asin(1.0);
  const double tol=1.0e-12;

  int prel0 = (creal(prel)==0.0 && cimag(prel)==0.0);

  // Preload tmat_n
  for(int ir1=0; ir1<anumLIZ; ir1++)
    for(int i=0; i<szBlk_tmat; i++)
      tmat_n[i+szBlk_tmat*ir1]=tmat0[i+LIZStoreIdx[ir1]*szBlk_tmat];

  #pragma acc data \
    copyin( tmat_n[0:akkrsz*akkrsz*4*anumLIZ],                                \
            lofk[0:ndlj],mofk[0:ndlj], ilp1[0:2*lmax+1], illp[0:kkrsz*kkrsz], \
            LIZPos[0:3*anumLIZ], LIZlmax[0:anumLIZ], LIZStoreIdx[0:anumLIZ],  \
            cgnt[0:(lmax+1)*kkrsz*kkrsz], clm[:(2*lmax+1)*(2*lmax+2)/2] )     \
    copyout( m[0:kkr1_ns*anumLIZ*kkr2_ns*anumLIZ] ) \
    create( dlm[0:ndlj*anumLIZ*anumLIZ] )
  {  // Data region start

  #pragma acc kernels loop independent collapse(2) vector
  for(int j=0; j<nrmat_ns; j++)
  for(int i=0; i<nrmat_ns; i++)
    m[i+j*nrmat_ns]=(i==j?1.0:0.0);

  #pragma acc kernels loop independent collapse(2) vector                     \
    private( sinmp,cosmp,plm,hfn )
  for(int ir1=0; ir1<anumLIZ; ir1++)
  for(int ir2_0=0; ir2_0<anumLIZ-1; ir2_0++)
  {
    // Skip the diagonals
    int ir2 = ir2_0 + (ir2_0>=ir1?1:0);

    int lmaxi=LIZlmax[ir1];
    int lmaxj=LIZlmax[ir2];
    int lend=lmaxi+lmaxj;

    double rij_1,rij_2,rij_0;

    rij_0=LIZPos[0+3*ir1]-LIZPos[0+3*ir2];
    rij_1=LIZPos[1+3*ir1]-LIZPos[1+3*ir2];
    rij_2=LIZPos[2+3*ir1]-LIZPos[2+3*ir2];

    double rmag=sqrt(rij_0*rij_0+rij_1*rij_1+rij_2*rij_2);

    //if(cabs(prel)==0.0)
    if (prel0)
    {
      // if prel is zero then we calculate Gij for multipole Coulomb interaction.
      hfn[0]=1.0/rmag;
      for(int l=1; l<=lend; l++)
      {
        // [eric] the commented one is correct on CPU but incorrect on GPU
        //hfn[l]=J*(2.0*l-1.0)*hfn[l-1]/rmag;
        hfn[l]=J*(2.0*l-1.0)*hfn[l-1]*hfn[0];
      }
    }
    else
    {
      double complex z=prel*rmag;
      hfn[0]=-J;
      hfn[1]=-1.0-J/z;
      for(int l=1; l<lend; l++) hfn[l+1]=(2.0*l+1)*hfn[l]/z - hfn[l-1];
      /*
                 l+1
         hfn = -i   *h (k*R)*sqrt(E)
                 l    ij
      */

      //z=cexp(J*z)/rmag;
      z = (cos(creal(z))+J*sin(creal(z)))/exp(cimag(z))/rmag;

      for(int l=0; l<=lend;l++) hfn[l]=-hfn[l]*z*ilp1[l];
    }

    //     calculate p(l,m)'s...............................................
    double costheta=rij_2/rmag;

    if((1.0-abs(costheta))<tol)
    {
      for(int i=0; i<(lend+1)*(lend+2)/2;i++) plm[i]=0.0;
      if(costheta<0.0)
      {
        for(int l=0; l<=lend; l++)
        {
          int i=l*(l+1)/2;
          plm[i]=1.0-2.0*(l%2);
        }
      }
      else
      {
        for(int l=0; l<=lend; l++)
        {
          int i=l*(l+1)/2;
          plm[i]=1.0;
        }
      }
    }
    else
    {
      double fact;
      int ll;
      plm[0]=1.0;
      /* #ifndef LMAX_GE_2 if(lend==0) return; #endif */

      double somx2=-sqrt((1.0-costheta)*(1.0+costheta));

      plm[1]=costheta;
      plm[2]=somx2;
      /* #ifndef LMAX_GE_2 if(lend==1) return; #endif */
      // m==0 special case
      fact=costheta;
      ll=0;
      for(int l=2; l<=lend; l++)
      {
        double pmm=(l-1.0)*plm[ll];
        fact+=2.0*costheta;
        ll+=l-1;
        plm[ll+l]=(fact*plm[ll]-pmm)/l;
      }

      for(int m=1; m<lend; m++)
      {
        double pmm=somx2;
        fact=1.0;
        for(int i=2; i<=m; i++)
        {
          fact+=2.0;
          pmm*=fact*somx2;
        }
        int mm=m*(m+1)/2+m;
        plm[mm]=pmm;
        if( mm+m+1<(lend+1)*(lend+2)/2 )
          plm[mm+m+1]=costheta*(2.0*m+1.0)*pmm;

        ll=(m+1)*m/2+m;
        fact=(2.0*m+1.0)*costheta;
        for(int l=m+2; l<=lend; l++)
        {
          pmm=(l+m-1.0)*plm[ll];
          fact+=2.0*costheta;
          ll=ll+l-1;
          plm[ll+l]=( fact*plm[ll] - pmm )/(l-m);
        }
      }
      // m==lend
      {
        int m=lend;
        double pmm=somx2;
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

    int ndlm_local=(lend+1)*(lend+2)/2;

    //// [eric] disabled for GPU offloading
    //if(ndlm_local>ndlm)
    //{
    //  printf("MAKEGIJ:: ndlm incorrect!\n");
    //  printf("ndlm=%d\nndlm_local=%d\n",ndlm,ndlm_local);
    //  exit(1);
    //}

    for(int j=0;j<ndlm_local;j++) plm[j]=clm[j]*plm[j];

    double pmag=sqrt(rij_0*rij_0+rij_1*rij_1);
    cosmp[0]=1.0;
    sinmp[0]=0.0;
    if(pmag>ptol)
    {
      cosmp[1]=rij_0/pmag;
      sinmp[1]=rij_1/pmag;
    }
    else
    {
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
      double m1m=1.0;
      dlm[MPdlm(j,ir1,ir2)]= hfn[l]*plm[ll];
      for(int m=1; m<=l; m++)
      {
        m1m=-m1m;
        double complex fac=plm[ll+m]*(cosmp[m]+J*sinmp[m]);
        dlm[MPdlm(j-m,ir1,ir2)]= hfn[l]*m1m*fac;
        dlm[MPdlm(j+m,ir1,ir2)]= hfn[l]*conj(fac);
      }
    }
  }

  #pragma acc parallel loop gang collapse(2) private(gij)
  for(int ir1=0; ir1<anumLIZ; ir1++)
  for(int ir2_0=0; ir2_0<anumLIZ-1; ir2_0++)
  {

    // Skip the diagonals
    int ir2 = ir2_0 + (ir2_0>=ir1?1:0);

    #pragma acc cache(gij)

    #pragma acc loop vector collapse(2)
    for(int j=0; j<kkr2; j++)
    for(int i=0; i<kkr1; i++)
    {
      int l1=lofk[j];
      int m1=mofk[j];
      int l2=lofk[i];
      int m2=mofk[i];
      int m3=m2-m1;

      int llow=(abs(m3)>abs(l1-l2)?abs(m3):abs(l1-l2));

      if (prel0) llow=l1+l2;

      gij[i+j*kkr1] = 0.0;

      for(int l3=l1+l2; l3>=llow;l3-=2)
        gij[i+j*kkr1]+=
          cgnt[l3/2+j*(lmax+1)+i*(lmax+1)*kkrsz]*dlm[MPdlm(l3*(l3+1)+m3,ir1,ir2)];

      gij[i+j*kkr1]=pi4*illp[i+kkrsz*j]*gij[i+j*kkr1];
    }

    #pragma acc loop vector tile(kkr1,kkr2)
    for(int j=0; j<kkr2_ns; j++)
    for(int i=0; i<kkr1_ns; i++)
    {
      double complex mij = 0.0;

      int kk=i+szBlk_tmat*ir1;
      int jj=j;

      if (2*j-kkr2_ns+1>0)
      {
        kk+=kkr1*kkr1_ns;
        jj-=kkr2;
      }
      jj=jj*kkr1;

      for(int k=0; k<kkr1; k++)
      {
        mij -= tmat_n[kk]*gij[jj];
        kk +=kkr1_ns;
        jj++;
      }

      m[kkr1_ns*ir1+i+(kkr2_ns*ir2+j)*kkr1_ns*anumLIZ]=mij;
    }

  }

  }  // Data region end

}
