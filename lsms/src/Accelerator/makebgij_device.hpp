
#include<math.h>
#include "cudaDoubleComplex.hpp"


#define FIDX(i,j,ldim) (j*ldim+i) 
#include "plglmax_device.hpp"  

__device__
__inline__ void makebgij_device(int lmaxi,int kkri,int lmaxj,int kkrj,
                    int lmax,int kkrsz,int ndlj,int ndlm, int nspin,
                    cudaDoubleComplex prel,double rij[3],double *sinmp,double *cosmp,
                    double *clm,double *plm,DeviceArray3d<double> &cgnt,int lmax_cg,int *lofk,int *mofk,
                    cudaDoubleComplex *ilp1,DeviceMatrix<cudaDoubleComplex> &illp,
                    cudaDoubleComplex *hfn,cudaDoubleComplex *dlm, cudaDoubleComplex *bgij,
                    double pi4) {
  const cudaDoubleComplex sqrtm1=cudaDoubleComplex(0.0,1.0);
  const double ptol=1.0e-6;
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
  double rmag=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);

#if 0
  if(threadIdx.x==0) {
  printf("ilp1 (lend: %d\n", lend);
  for(int i=0;i<=lend;i++)
    printf("%d, (%lg, %lg)\n",i,ilp1[i].x,ilp1[i].y);
  }
#endif
  /**************************SERIAL************************************/
  if(threadIdx.x==0) {
#if 0
    printf("prel: (%lg, %lg)\n", prel.x, prel.y);
    printf("sqrtm1: (%lg, %lg)\n", sqrtm1.x,sqrtm1.y);
    printf("rmag: %lg\n",rmag);
#endif
    //      if(rmag.lt.ptol) stop 'Error in makegij: rmag = 0.'
    if(prel.x==0.0 && prel.y==0.0)
    {
      // if prel is zero then we calculate Gij for multipole Coulomb interaction.
      hfn[0]=cudaDoubleComplex(1.0/rmag,0);
      for(int l=1; l<=lend; l++)
      {
        hfn[l]=((2.0*l-1))/rmag*sqrtm1*hfn[l-1];  //Why does this crash?
      }
    } else {
      cudaDoubleComplex z=prel*rmag;
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
      z=exp(sqrtm1*z)/rmag;
      //printf("z: (%lg,%lg)\n",z.x,z.y);
      for(int l=0; l<=lend;l++)  //TODO, this could be parallel, but probably not worth it...
      {
        hfn[l]=-hfn[l]*z*ilp1[l];     
      }
    }
#if 0
    if ( blockIdx.x == 0 && blockIdx.y == 1 )
    {
      printf("i, hfn(i)\n");
      for(int i=0;i<(2*lmax+1);i++) {
        printf("%d, (%lg, %lg)\n",i,hfn[i].real(), hfn[i].imag());
      }
    }
#endif
  }
  /*********************END  SERIAL************************************/
  //     calculate p(l,m)'s...............................................
  double costheta=rij[2]/rmag;
  //      plglmax_(&lend,&costheta,plm);


#if 0
  if(threadIdx.x==0) {
  printf("lend: %d, costheta: %lg\n",lend,costheta);
  }
#endif
  //plglmax_device(lend,costheta,plm); //JLARKIN replaced with associatedLegendreFunctionNormalized
  //JLARKIN - DEBUG
  //for(int i = threadIdx.x; i <= lend; i += blockDim.x) plm[i] = 0.0;
  associatedLegendreFunctionNormalizedDevice(costheta,lend,plm);
  __syncthreads(); //wait for plm to be updated
#if 0
  if(threadIdx.x==0 && blockIdx.x == 0 && blockIdx.y == 1) {
  //printf("Rij: %lg, %lg, %lg\n",rij[0],rij[1],rij[2]);
  //printf("cos(theata): %lg\n", costheta);
  printf("i, plm(i)\n");
  for(int i=0;i<ndlm;i++) {
    printf("%d, %lg\n",i,plm[i]);
  }
  }
 __syncthreads(); //wait for plm to be output
 return;
#endif
  //     multiply be the normalization constant...........................
  int ndlm_local=(lend+1)*(lend+2)/2;

  //for(int t=0;t<ndlm_local;t++) {  //SERIAL
  for(int t=threadIdx.x;t<ndlm_local;t+=blockDim.x) { //PARALLEL
    plm[t]=clm[t]*plm[t];
  }
  __syncthreads(); //wait for plm to be updated

  //     =================================================================
  //     calculate cos(phi) and sin(phi) .................................
  /**************************SERIAL************************************/
  if(threadIdx.x==0) { 
    double pmag=sqrt(rij[0]*rij[0]+rij[1]*rij[1]);
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
#if 0
  if(threadIdx.x==0) {
  printf("i, cosmp(i), sinmp(i)\n");
  for(int i=0;i<2*lmax+1;i++) {
    printf("%d, %lg, %lg\n",i, cosmp[i], sinmp[i]);
  }
  }
#endif
  
    int j=0;
    for(int l=0; l<=lend; l++)
    {
      int ll=l*(l+1);
      j=ll+1;
      ll=ll/2+1;
      double m1m=1.0;
      dlm[j-1]= hfn[l]*plm[ll-1];
      for(int m=1; m<=l; m++)
      {
        m1m=-m1m;
        cudaDoubleComplex fac=plm[ll+m-1]*cudaDoubleComplex(cosmp[m],sinmp[m]);
        dlm[j-m-1]= hfn[l]*m1m*fac;
        dlm[j+m-1]= hfn[l]*conj(fac);
      }
    }
  }
  /*********************END  SERIAL************************************/
  __syncthreads(); //wait for dlm to be updated
#if 0
  if(threadIdx.x==0) {
  printf("i, dlm(i)\n");
  for(int i=0;i<ndlj;i++) {
    printf("%d, (%lg, %lg)\n",i,dlm[i].real(), dlm[i].imag());
  }
  }
#endif
  

  //for(int t=0;t<kkrj*kkri;t++) //SERIAL
  for(int t=threadIdx.x;t<kkrj*kkri;t+=blockDim.x) //PARALLEL
  { 
    int lm1=t/kkri;
    int lm2=t%kkri;
    int l1=lofk[lm1];
    int m1=mofk[lm1];
    int l2=lofk[lm2];
    int m2=mofk[lm2];
    int m3=m2-m1;

    int llow=max(abs(m3),abs(l1-l2));
    if(prel.x==0.0 && prel.y==0.0) llow=l1+l2;
    cudaDoubleComplex sum=cudaDoubleComplex(0.,0.);
    for(int l3=l1+l2; l3>=llow;l3-=2)
    {
      int j=l3*(l3+1)+m3;
      sum = sum + cgnt(l3/2,lm1,lm2)*dlm[j];
    }
    sum = sum * pi4*illp(lm2,lm1);

    for(int x=0;x<nspin;x++) {
      int i=lm2;
      int j=lm1;
      int bi=x*kkri+i;//TODO  verify i and j are correct and shouldn't be swapped
      int bj=x*kkrj+j;//TODO
      int ldim=kkri*nspin; //TODO

      bgij[FIDX(bi,bj,ldim)]=sum;
    }
  }
}
