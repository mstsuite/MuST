#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <complex.h>
#include "cuComplex.h"
#include "acclib.hpp"
#include "accmath.hpp"

__constant__ int lmax_kkr_max = 8;
__constant__ int lmax_max = 16; // = 2*lmax_kkr_max
__constant__ int jmax_max = 153; // = (lmax_max+1)*(lmax_max+2)/2
__constant__ int kmax_max = 289; // = (lmax_max+1)**2

__device__ void return_posi_index(int liz_max, int idx, int *ida, int *ip, int *jp) {
   int num_pairs = liz_max*liz_max;
   *ida = idx/num_pairs;
   int idp = idx%num_pairs;
   *jp = idp/liz_max;
   *ip = idp%liz_max;
}

// The CUDA kernel
__global__ void computeSphericalHarmonicsKernel(int liz_max, 
                                                int lmax, 
                                                int N_size, 
                                                int *p_num_nbs_d,
                                                double *p_posi_d,
                                                double *p_Clm_d,
                                                cuDoubleComplex *p_Ylm_d) {
   double Plm[153]; //jmax=153
// int jmax = (lmax+1)*(lmax+2)/2;
// double *Plm = (double*)malloc(jmax * sizeof(double));
   cuDoubleComplex czero = make_cuDoubleComplex(0.0, 0.0);

   int idx = blockIdx.x * blockDim.x + threadIdx.x;

   if (idx < N_size) {
      int ida = 0;
      int ip = 0;
      int jp = 0;

      return_posi_index(liz_max, idx, &ida, &ip, &jp);

      int n = ida*(lmax+1)*(lmax+1)*liz_max*liz_max+jp*(lmax+1)*(lmax+1)*liz_max
                                                   +ip*(lmax+1)*(lmax+1);

      if (ip != jp && ip <= p_num_nbs_d[ida] && jp <= p_num_nbs_d[ida]) {
         int j3 = ida*liz_max*3+jp*3;
         int i3 = ida*liz_max*3+ip*3;
         double x = p_posi_d[j3]   - p_posi_d[i3];
         double y = p_posi_d[j3+1] - p_posi_d[i3+1];
         double z = p_posi_d[j3+2] - p_posi_d[i3+2];
         double cos_theta = z/sqrt(x*x + y*y + z*z);
         double phi = atan2(y,x);

         computeLegendre(lmax, cos_theta, Plm);
        
         cuDoubleComplex cone =  make_cuDoubleComplex(1.0,0.0);

         for (int l = 0; l <= lmax; ++l) {
            int jl = (l+1)*(l+2)/2-l;
            int kl = (l+1)*(l+1)-l;
            int nl = n + kl;
            double sign_factor = 1.0;
            for (int m = 0; m <= l; ++m) {
               double P_lm = Plm[jl-1];
               double C_lm = p_Clm_d[jl-1];
                   
               double m_phi = m * phi;
               double cr = C_lm*P_lm*cos(m_phi);
               double ci = C_lm*P_lm*sin(m_phi);
             //cuDoubleComplex e_im_phi = make_cuDoubleComplex(cos(m_phi), sin(m_phi));
             //cuDoubleComplex Ylm = C_lm * P_lm * e_im_phi;
             //cuDoubleComplex Ylm_conj = C_lm * P_lm / e_im_phi;
               cuDoubleComplex Ylm = make_cuDoubleComplex(cr, ci);
               cuDoubleComplex Ylm_nm = make_cuDoubleComplex(sign_factor*cr,-sign_factor*ci);
               p_Ylm_d[nl+m-1] = Ylm;
               p_Ylm_d[nl-m-1] = Ylm_nm;
               jl++;
               sign_factor = -sign_factor;
            }
         }
      }
      else {
         int kmax = (lmax+1)*(lmax+1);
         p_Ylm_d[n] = make_cuDoubleComplex(p_Clm_d[0],0.0);  // In the original CPU code, it is ylm(1) = clm(1)
         for (int kl = 2; kl <= kmax; kl++) {
            p_Ylm_d[n+kl-1] = czero;
         }
      }
      __syncthreads();
   }
// free(Plm);
}

__global__ void computeGijKernel(cuDoubleComplex kappa_d, int liz_max,
                                 int aid, int liz_size, int dsize, int lmax, int cant,
                                 int MaxJ3, int cgnt_kmax,
                                 double *p_posi_d, 
                                 int *p_lofk_d, 
                                 int *p_nj3_d, 
                                 int *p_kj3_d,
                                 double *p_cgnt_d,
                                 cuDoubleComplex *p_Ylm_d,
                                 cuDoubleComplex *p_gij_d) {
   const double pi4=3.14159265358979*4.0;
   cuDoubleComplex ylmcc[289]; // (lmax_max+1)**2=289
   cuDoubleComplex dlm[289];
   cuDoubleComplex hfn[17]; // lmax_max+1=17
   cuDoubleComplex i2l[17]; // lmax_max+1=17
   cuDoubleComplex pi4c = make_cuDoubleComplex(pi4, 0.0);
   cuDoubleComplex czero = make_cuDoubleComplex(0.0, 0.0);
   cuDoubleComplex cone = make_cuDoubleComplex(1.0, 0.0);
   cuDoubleComplex sqrtm1 = make_cuDoubleComplex(0.0, 1.0);
   cuDoubleComplex neg_sqrtm1 = make_cuDoubleComplex(0.0, -1.0);

   i2l[0] = cone;
   for (int l=1; l<=lmax_max; l++) {
      i2l[l] = cuCmul(sqrtm1,i2l[l-1]);
   }

   int idx = blockIdx.x * blockDim.x + threadIdx.x;

   // dsize = big matrix size
   if (idx < dsize*dsize) {
      int kmax = (lmax+1)*(lmax+1);
      int NK = dsize*kmax;
      int jp = idx/NK;         // jp = column-wise atom index (starts from 0)
      int m =  idx%NK;
      int kj = m/dsize;        // kj = column-wise kl index - 1
      int ip = (m%dsize)/kmax; // ip = row-wise atom index (starts from 0)
      int ki = (m%dsize)%kmax; // ki = row-wise kl index - 1

      if (ip != jp && ip < liz_size && jp < liz_size) {
         int ni = (aid-1)*liz_max*3 + 3*ip;
         int nj = (aid-1)*liz_max*3 + 3*jp;
         double x = p_posi_d[nj]   - p_posi_d[ni];
         double y = p_posi_d[nj+1] - p_posi_d[ni+1];
         double z = p_posi_d[nj+2] - p_posi_d[ni+2];
         cuDoubleComplex rmag = make_cuDoubleComplex(sqrt(x*x+y*y+z*z),0.0);

         cuDoubleComplex kr = cuCmul(kappa_d,rmag);
         int lmax_dlm = 2*lmax;
         int kmax_dlm = (lmax_dlm+1)*(lmax_dlm+1);

         int num_pairs = liz_max*liz_max;
         int n = (aid-1)*num_pairs*kmax_dlm + jp*liz_max*kmax_dlm + ip*kmax_dlm;
         for (int kl=0; kl<kmax_dlm; kl++) {
            ylmcc[kl] = cuConj(p_Ylm_d[n+kl]);
         }

         // z=kappa*rmag
         // hfn(0)=-sqrtm1
         // if ( lmax_dlm > 0 ) then
         //    hfn(1)=-(cone+sqrtm1/z)
         //    do l=2,lmax_dlm
         //       hfn(l)=(2*l-1)*hfn(l-1)/z - hfn(l-2)
         //    enddo
         // endif
         hfn[0]=neg_sqrtm1;
         if ( lmax_dlm > 0 ) {
            hfn[1] = cuCsub(czero,cuCadd(cone,cuCdiv(sqrtm1,kr)));
            for (int l=2; l <= lmax_dlm; l++) {
               cuDoubleComplex cfac = make_cuDoubleComplex(2.0*l-1.0,0.0);
               hfn[l] = cuCsub(cuCdiv(cuCmul(cfac,hfn[l-1]),kr),hfn[l-2]);
            }
         }
         // ==========================================================
         // generate the KKR real space lattice structure matrix for
         // the energy and store result in gij
         //
         //            l+1
         //    fac = -i   *h (k*R  )*sqrt(E)
         //                 l    ij
         // ==========================================================
         //   z = exp(sqrtm1*z)/rmag
         //   do kl = 1,kmax_dlm
         //      l =lofk(kl)
         //      fac =  hfn(l)*z/ilp1(l)
         //      dlm(kl) = fac*ylmcc(kl)
         //   enddo
         // ==========================================================
         cuDoubleComplex cz =cuCdiv(cuCexp(cuCmul(sqrtm1,kr)),rmag);
         for (int kl=0; kl<kmax_dlm; kl++) {
            int l = p_lofk_d[kl];
            cuDoubleComplex ilp1 = cuCmul(i2l[l],sqrtm1);
            cuDoubleComplex cfac = cuCdiv(cuCmul(hfn[l],cz),ilp1);
            dlm[kl] = cuCmul(cfac,ylmcc[kl]);
         }

         // ==========================================================
         // loop over klp.............................................
         // do klp=1,kmaxj
         //    lp=lofk(klp)
         //    =======================================================
         //    loop over kl...........................................
         //    =======================================================
         //    do kl=1,kmaxi
         //       l=lofk(kl)
         //       ====================================================
         //                     l-lp
         //       illp(l,lp) = i
         //
         //       perform sum over j with gaunt # ....................
         //       ====================================================
         //       nnj3 = nj3(klp,kl)
         //       pkj3 => kj3(1:nnj3,klp,kl)
         //       pcgnt=> cgnt(1:nnj3,klp,kl)
         //       gij_llp = CZERO
         //       do j = 1,nnj3
         //          gij_llp = gij_llp+pcgnt(j)*dlm(pkj3(j))
         //       enddo
         //       gij(kl,klp)=pi4*illp(kl,klp)*gij_llp
         //    enddo
         // enddo
         // ==========================================================
         n = ki*cgnt_kmax+kj;
         int nnj3 = p_nj3_d[n];
         int nc = ki*MaxJ3*cgnt_kmax+kj*MaxJ3;
         cuDoubleComplex gij_llp = czero;
         for (int j=0; j<nnj3; j++) {
            int pkj3 = p_kj3_d[nc+j]-1;
            // gij_llp = gij_llp+pcgnt(j)*dlm(pkj3(j))
            cuDoubleComplex c = make_cuDoubleComplex(p_cgnt_d[nc+j],0.0);
            gij_llp = cuCadd(gij_llp,cuCmul(c,dlm[pkj3]));
         }
         int lp = p_lofk_d[kj];
         int l  = p_lofk_d[ki];
         cuDoubleComplex c = cuCdiv(i2l[l],i2l[lp]);
         cuDoubleComplex gij = cuCmul(cuCmul(pi4c,c),gij_llp);
         if (cant == 1) {
            p_gij_d[idx] = gij;
         }
         else {
            n = jp*dsize*cant*kmax*cant+kj*dsize*cant+ip*kmax*cant+ki;
            p_gij_d[n] = gij;
         // p_gij_d[n+kmax] = czero;
         // p_gij_d[n+kmax*dsize*cant] = czero;
            p_gij_d[n+kmax*dsize*cant+kmax] = gij;
         }
      }
/*
      else {
         if (cant == 1) {
            p_gij_d[idx] = czero;
         }
         else {
            int n = jp*dsize*cant*kmax*cant+kj*dsize*cant+ip*kmax*cant+ki;
            p_gij_d[n] = czero;
            p_gij_d[n+kmax] = czero;
            p_gij_d[n+kmax*dsize*cant] = czero;
            p_gij_d[n+kmax*dsize*cant+kmax] = czero;
         }
      }
*/
   }
   __syncthreads();
}

__global__ void calculateGijKernel(cuDoubleComplex kappa_d, int liz_max,
                                   int aid, int liz_size, int dsize, int lmax, int cant,
                                   int MaxJ3, int cgnt_kmax,
                                   double *p_posi_d, 
                                   int *p_lofk_d, 
                                   int *p_nj3_d, 
                                   int *p_kj3_d,
                                   double *p_cgnt_d,
                                   cuDoubleComplex *p_Ylm_d,
                                   cuDoubleComplex *p_gij_d) {
   const double pi4=3.14159265358979*4.0;
   cuDoubleComplex ylmcc[289]; // (lmax_max+1)**2=289
   cuDoubleComplex dlm[289];
   cuDoubleComplex hfn[17]; // lmax_max+1=17
   cuDoubleComplex i2l[17]; // lmax_max+1=17
   cuDoubleComplex pi4c = make_cuDoubleComplex(pi4, 0.0);
   cuDoubleComplex czero = make_cuDoubleComplex(0.0, 0.0);
   cuDoubleComplex cone = make_cuDoubleComplex(1.0, 0.0);
   cuDoubleComplex sqrtm1 = make_cuDoubleComplex(0.0, 1.0);
   cuDoubleComplex neg_sqrtm1 = make_cuDoubleComplex(0.0, -1.0);

   i2l[0] = cone;
   for (int l=1; l<=lmax_max; l++) {
      i2l[l] = cuCmul(sqrtm1,i2l[l-1]);
   }

   int idx = blockIdx.x * blockDim.x + threadIdx.x;

   // dsize = big matrix size
   if (idx < dsize*dsize) {
      int kmax    = (lmax+1)*(lmax+1);
      int kmax_ns = kmax*cant;
      int NK      = dsize*kmax_ns;
      int jp      = idx/NK;             // jp = column-wise atom index (starts from 0)
      int MK      = idx%NK;
      int m       = MK/dsize;           // kj = column-wise kl index - 1
      int jcant   = m/kmax;             // jcant = 0, 1
      int kj      = m%kmax;             // kj = column-wise kl index - 1
      int ip      = (MK%dsize)/kmax_ns; // ip = row-wise atom index (starts from 0)
      int n       = (MK%dsize)%kmax_ns; //
      int icant   = n/kmax;             // icant = 0, 1
      int ki      = n%kmax;             // ki = row-wise kl index - 1

      if (ip != jp && icant == jcant && ip < liz_size && jp < liz_size) {
         int ni = (aid-1)*liz_max*3 + 3*ip;
         int nj = (aid-1)*liz_max*3 + 3*jp;
         double x = p_posi_d[nj]   - p_posi_d[ni];
         double y = p_posi_d[nj+1] - p_posi_d[ni+1];
         double z = p_posi_d[nj+2] - p_posi_d[ni+2];
         cuDoubleComplex rmag = make_cuDoubleComplex(sqrt(x*x+y*y+z*z),0.0);

         cuDoubleComplex kr = cuCmul(kappa_d,rmag);
         int lmax_dlm = 2*lmax;
         int kmax_dlm = (lmax_dlm+1)*(lmax_dlm+1);

         int num_pairs = liz_max*liz_max;
         n = (aid-1)*num_pairs*kmax_dlm + jp*liz_max*kmax_dlm + ip*kmax_dlm;
         for (int kl=0; kl<kmax_dlm; kl++) {
            ylmcc[kl] = cuConj(p_Ylm_d[n+kl]);
         }

         // z=kappa*rmag
         // hfn(0)=-sqrtm1
         // if ( lmax_dlm > 0 ) then
         //    hfn(1)=-(cone+sqrtm1/z)
         //    do l=2,lmax_dlm
         //       hfn(l)=(2*l-1)*hfn(l-1)/z - hfn(l-2)
         //    enddo
         // endif
         hfn[0]=neg_sqrtm1;
         if ( lmax_dlm > 0 ) {
            hfn[1] = cuCsub(czero,cuCadd(cone,cuCdiv(sqrtm1,kr)));
            for (int l=2; l <= lmax_dlm; l++) {
               cuDoubleComplex cfac = make_cuDoubleComplex(2.0*l-1.0,0.0);
               hfn[l] = cuCsub(cuCdiv(cuCmul(cfac,hfn[l-1]),kr),hfn[l-2]);
            }
         }
         // ==========================================================
         // generate the KKR real space lattice structure matrix for
         // the energy and store result in gij
         //
         //            l+1
         //    fac = -i   *h (k*R  )*sqrt(E)
         //                 l    ij
         // ==========================================================
         //   z = exp(sqrtm1*z)/rmag
         //   do kl = 1,kmax_dlm
         //      l =lofk(kl)
         //      fac =  hfn(l)*z/ilp1(l)
         //      dlm(kl) = fac*ylmcc(kl)
         //   enddo
         // ==========================================================
         cuDoubleComplex cz =cuCdiv(cuCexp(cuCmul(sqrtm1,kr)),rmag);
         for (int kl=0; kl<kmax_dlm; kl++) {
            int l = p_lofk_d[kl];
            cuDoubleComplex ilp1 = cuCmul(i2l[l],sqrtm1);
            cuDoubleComplex cfac = cuCdiv(cuCmul(hfn[l],cz),ilp1);
            dlm[kl] = cuCmul(cfac,ylmcc[kl]);
         }

         // ==========================================================
         // loop over klp.............................................
         // do klp=1,kmaxj
         //    lp=lofk(klp)
         //    =======================================================
         //    loop over kl...........................................
         //    =======================================================
         //    do kl=1,kmaxi
         //       l=lofk(kl)
         //       ====================================================
         //                     l-lp
         //       illp(l,lp) = i
         //
         //       perform sum over j with gaunt # ....................
         //       ====================================================
         //       nnj3 = nj3(klp,kl)
         //       pkj3 => kj3(1:nnj3,klp,kl)
         //       pcgnt=> cgnt(1:nnj3,klp,kl)
         //       gij_llp = CZERO
         //       do j = 1,nnj3
         //          gij_llp = gij_llp+pcgnt(j)*dlm(pkj3(j))
         //       enddo
         //       gij(kl,klp)=pi4*illp(kl,klp)*gij_llp
         //    enddo
         // enddo
         // ==========================================================
         n = ki*cgnt_kmax+kj;
         int nnj3 = p_nj3_d[n];
         int nc = ki*MaxJ3*cgnt_kmax+kj*MaxJ3;
         cuDoubleComplex gij_llp = czero;
         for (int j=0; j<nnj3; j++) {
            int pkj3 = p_kj3_d[nc+j]-1;
            cuDoubleComplex c = make_cuDoubleComplex(p_cgnt_d[nc+j],0.0);
            gij_llp = cuCadd(gij_llp,cuCmul(c,dlm[pkj3]));
         }
         int lp = p_lofk_d[kj];
         int l  = p_lofk_d[ki];
         cuDoubleComplex c = cuCdiv(i2l[l],i2l[lp]);
         cuDoubleComplex gij = cuCmul(cuCmul(pi4c,c),gij_llp);
         p_gij_d[idx] = gij;
      }
   }
   __syncthreads();
}

bool Ylm_allocated = false;
bool param_pushed = false;
bool initialized = false;

int kj3_MaxJ3 = 0;
int kj3_kmax = 0;
int *nj3_d;
int *kj3_d;
int *lofk_d;
double *cgnt_d;

double *Clm_d;
cuDoubleComplex *Ylm_d;

int *num_nbs_d;
double *posi_d;
cuDoubleComplex  *gij_d;

void runYlmOverLIZ(int na, int liz_max, int l_max) {

      // Define kernel launch parameters
      int N_size = na*liz_max*liz_max;
      int threadsPerBlock = 256;
      int numBlocks = (N_size + threadsPerBlock - 1) / threadsPerBlock;

//    testKernel<<<numBlocks, threadsPerBlock>>>(N_size,lofk_d,posi_d);
//    checkCudaErrors(cudaDeviceSynchronize());
  
      computeSphericalHarmonicsKernel<<<numBlocks, threadsPerBlock>>>(liz_max, 
                                                                      l_max, 
                                                                      N_size, 
                                                                      num_nbs_d,
                                                                      posi_d, 
                                                                      Clm_d, 
                                                                      Ylm_d);
      checkCudaErrors(cudaDeviceSynchronize());
      checkCudaErrors(cudaGetLastError()); // Check for launch errors
}

int tau_size = 0;
int mmat_size = 0;

double _Complex  *sine_h;  // mat_id = 1
double _Complex  *jinv_h;  // mat_id = 2
double _Complex  *gij_h;   // mat_id = 3

cuDoubleComplex  *sine_d;
cuDoubleComplex  *jinv_d;
cuDoubleComplex  *BigMat_d;
cuDoubleComplex  *BigMatInv_d;
cuDoubleComplex  *block_d;
int *pivotArray;
int *infoArray;

extern "C"
void calculate_gij_gpu_(double _Complex *kappa, int *n_spin_cant, int *numnb_max, 
                        int *ia, int *num_nbs, int *lmax_kkr) {
   if (!initialized) {
      fprintf(stderr, "\nError in calculate_gij_gpu: Needs to call init_lsms_gpu first.\n");
      exit(EXIT_FAILURE);
   }
   
   // double _Complex e2 = kappa*kappa;
   // std::complex<double> e2 = kappa*kappa;
   // printf("energy = %.5f + %.5fi\n", e2.real(), e2.imag());

   // total number of rij pairs
   int liz_max = *numnb_max+1;
   int aid = *ia;
   int liz_size = *num_nbs+1;
   int lmax = *lmax_kkr;
   int cant = *n_spin_cant;

   cuDoubleComplex kappa_d = make_cuDoubleComplex(creal(*kappa),cimag(*kappa));

   int dsize = liz_size*(lmax+1)*(lmax+1); // = big matrix rank/cant
   // Define kernel launch parameters
   int threadsPerBlock = 256;
   int numBlocks = (dsize*dsize + threadsPerBlock - 1) / threadsPerBlock;
   // Launch the kernel
   computeGijKernel<<<numBlocks,threadsPerBlock>>>(kappa_d,liz_max,aid,liz_size,dsize,lmax,cant,kj3_MaxJ3,kj3_kmax,
                                                   posi_d,lofk_d,nj3_d,kj3_d,cgnt_d,Ylm_d,gij_d);

/*
   int dsize = liz_size*(lmax+1)*(lmax+1)*cant; // = big matrix rank
   // Define kernel launch parameters
   int threadsPerBlock = 256;
   int numBlocks = (dsize*dsize + threadsPerBlock - 1) / threadsPerBlock;
   // Launch the kernel
   calculateGijKernel<<<numBlocks,threadsPerBlock>>>(kappa_d,liz_max,aid,liz_size,dsize,lmax,cant,kj3_MaxJ3,kj3_kmax,
                                                     posi_d,lofk_d,nj3_d,kj3_d,cgnt_d,Ylm_d,gij_d);
*/

   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaGetLastError()); // Check for launch errors
}

extern "C"
void get_gij_from_gpu_(int *ip, int *jp, int *lmax, double _Complex *gij) {
   int kmax = (*lmax+1)*(*lmax+1);
   if (initialized) {
      // --------------------------------------------
      // Aggregate the data into block_d array and then make one cudaMemcpy call
      // --------------------------------------------
      // Define grid and block dimensions for the kernel
      dim3 threadsPerBlock(16, 16); // Using a 32x32 block
      dim3 blocksPerGrid((kmax + threadsPerBlock.x - 1) / threadsPerBlock.x,
                         (kmax + threadsPerBlock.y - 1) / threadsPerBlock.y);

      // Launch the kernel and use block_d to store gij
      // ============================================
      copySubBlockKernel<<<blocksPerGrid, threadsPerBlock>>>(gij_d, mmat_size, *ip, *jp, block_d, kmax);
      checkCudaErrors(cudaPeekAtLastError());
      checkCudaErrors(cudaMemcpy(gij, block_d, sizeof(cuDoubleComplex)*kmax*kmax, 
                                 cudaMemcpyDeviceToHost));
   }
   else {
      printf("\nIt needs to call init_lsms_gpu first.\n");
   }
}

extern "C"
void init_lsms_gpu_(int *lmax, int *dsize, int *block_size, int *my_pe) {
   if (!initialized) {
      mmat_size = *dsize;
      tau_size = *block_size;

      // This code assumes that each atom in the LIZ has the same scattering matrix size
      // ==============================================================

      if (*my_pe == 0) printf("CUDA memory assigned \n");

      // cudaError_t error;
      size_t size_block = tau_size * tau_size * sizeof(cuDoubleComplex);
      size_t size_bigmat = mmat_size * mmat_size * sizeof(cuDoubleComplex);

      checkCudaErrors(cudaMalloc((void**)&BigMat_d, size_bigmat));

      checkCudaErrors(cudaMalloc((void**)&pivotArray, sizeof(int)*mmat_size));

      checkCudaErrors(cudaMalloc((void**)&infoArray, sizeof(int)));

      checkCudaErrors(cudaMalloc((void**)&BigMatInv_d, size_bigmat));

      checkCudaErrors(cudaMalloc((void**)&block_d, size_block));

      checkCudaErrors(cudaMalloc((void**)&sine_d, size_bigmat));

      checkCudaErrors(cudaMalloc((void**)&jinv_d, size_bigmat));

      sine_h = (double _Complex *) malloc(size_bigmat);
      memset(sine_h, 0, size_bigmat);    // set the host array to 0

      jinv_h = (double _Complex *) malloc(size_bigmat);
      memset(jinv_h, 0, size_bigmat);    // set the host array to 0

      gij_h = (double _Complex *) malloc(size_bigmat);
      memset(gij_h, 0, size_bigmat);     // set the host array to 0

      checkCudaErrors(cudaMalloc((void**)&gij_d, size_bigmat));
      cudaMemset(gij_d, 0, size_bigmat); // set the device array to 0

      initialized = true;
   }
}
extern "C"
void compute_ylm_gpu_(int *lmax, int *local_atoms, int *numnb_max, 
                     int *num_nbs, double *posi, int *my_pe) {
   if (!Ylm_allocated) {
      int lmax_kkr = *lmax;
      int lmax2 = 2*lmax_kkr;
      int liz_max = *numnb_max+1;
      int kmax2 = (lmax2+1)*(lmax2+1);

//    printf("\nNum of local atoms = %d\n",*local_atoms);
//    printf("\nNum of posi_d elements = %d\n",(*local_atoms)*liz_max*3);

      size_t size_nbs = (*local_atoms) * sizeof(int);
      checkCudaErrors(cudaMalloc((void**)&num_nbs_d, size_nbs));
      checkCudaErrors(cudaMemcpy(num_nbs_d, num_nbs, size_nbs, cudaMemcpyHostToDevice));

      size_t size_Ylm = (*local_atoms) * liz_max * liz_max * kmax2 * sizeof(cuDoubleComplex);
      checkCudaErrors(cudaMalloc((void**)&Ylm_d, size_Ylm));
      cudaMemset(Ylm_d, 0, size_Ylm); // set the device array to 0

      size_t size_posi = (*local_atoms) * liz_max * 3 * sizeof(double);
      checkCudaErrors(cudaMalloc((void**)&posi_d, size_posi));
      checkCudaErrors(cudaMemcpy(posi_d, posi, size_posi, cudaMemcpyHostToDevice));

//    printf("\nYlm_d and posi_d are allocated\n");
      Ylm_allocated = true;

      int na = *local_atoms;
      runYlmOverLIZ(na, liz_max, lmax2);
   }
}

extern "C"
void push_parameters_gpu_(int *lmax, int *lofk, int *kj3_size_1, int *kj3_size_2,
                          int *nj3, int *kj3, double *cgnt, double *Clm) {
   if (!param_pushed) {
      int jmax = (*lmax+1)*(*lmax+2)/2;
      int kmax = (*lmax+1)*(*lmax+1);
      size_t size_lofk_d = kmax * sizeof(int);
      // checkCudaErrors(cudaMalloc((void**)&lofk_d, size_lofk_d));
      checkCudaErrors(cudaMalloc((void**)&lofk_d, size_lofk_d));
      checkCudaErrors(cudaMemcpy(lofk_d, lofk, size_lofk_d, cudaMemcpyHostToDevice));
      // printf("\nLocation 1: lofk_d = %d,%d,%d,%d\n",lofk_d[0],lofk_d[1],lofk_d[2],lofk_d[3]);

      size_t size_Clm_d = jmax * sizeof(double);
      // checkCudaErrors(cudaMalloc((void**)&Clm_d, size_Clm_d));
      checkCudaErrors(cudaMalloc((void**)&Clm_d, size_Clm_d));
      checkCudaErrors(cudaMemcpy(Clm_d, Clm, size_Clm_d, cudaMemcpyHostToDevice));

      kj3_MaxJ3 = *kj3_size_1;
      kj3_kmax = *kj3_size_2;
      int nj3_size = kj3_kmax*kj3_kmax;
      int kj3_size = kj3_MaxJ3 * nj3_size;

      size_t size_nj3_d = nj3_size * sizeof(int);
      checkCudaErrors(cudaMalloc((void**)&nj3_d, size_nj3_d));
      checkCudaErrors(cudaMemcpy(nj3_d, nj3, size_nj3_d, cudaMemcpyHostToDevice));

      size_t size_kj3_d = kj3_size * sizeof(int);
      checkCudaErrors(cudaMalloc((void**)&kj3_d, size_kj3_d));
      checkCudaErrors(cudaMemcpy(kj3_d, kj3, size_kj3_d, cudaMemcpyHostToDevice));

      size_t size_cgnt_d = kj3_size * sizeof(double);
      checkCudaErrors(cudaMalloc((void**)&cgnt_d, size_cgnt_d));
      checkCudaErrors(cudaMemcpy(cgnt_d, cgnt, size_cgnt_d, cudaMemcpyHostToDevice));

      param_pushed = true;
   }
}

extern "C"
void finalize_lsms_gpu_() {
   if (initialized) {
      free(sine_h);
      free(jinv_h);
      free(gij_h);
      cudaFree(BigMat_d);
      cudaFree(pivotArray);
      cudaFree(infoArray);
      cudaFree(BigMatInv_d);
      cudaFree(block_d);
      cudaFree(sine_d);
      cudaFree(jinv_d);
      cudaFree(gij_d);
   }
   initialized = false;
}

extern "C"
void free_ylm_gpu_() {
   if (Ylm_allocated) {
      cudaFree(posi_d);
      posi_d = NULL;
      cudaFree(Ylm_d);
      Ylm_d = NULL;
   }
   Ylm_allocated = false;
}

extern "C"
void free_parameters_gpu_() {
   if (param_pushed) {
      cudaFree(lofk_d);
      cudaFree(nj3_d); 
      cudaFree(kj3_d); 
      cudaFree(cgnt_d); 
      cudaFree(Clm_d);
   }
   param_pushed = false;
}

extern "C"
void init_bigmatrix_gpu_(int *b_size) {
   if (*b_size != mmat_size) {
      fprintf(stderr,"\nError: b_size <> mmat_size, %d,%d\n",*b_size,mmat_size);
      exit(EXIT_FAILURE);
   }
   // Define kernel launch parameters and set BigMat_d to be a unit matrix
   // ========================================
   int threads_per_block = 512;
   int num_blocks = (mmat_size*mmat_size + threads_per_block - 1) / threads_per_block;
   createUnitMatrixKernel<<<num_blocks, threads_per_block>>>(BigMat_d, mmat_size);
   checkCudaErrors(cudaPeekAtLastError());
}

extern "C"
void push_submatrix_gpu_(int *mat_id, int *row, int *col, double _Complex *mat, int *msize ) {
   if (*msize != tau_size) {
      fprintf(stderr,"\nError: msize <> tau_size, %d,%d\n",*msize,tau_size);
      exit(EXIT_FAILURE);
   }
   else if (*mat_id == 1) {
      copySubBlockToMatrix(mat,msize,row,col,sine_h,&mmat_size);
   }
   else if (*mat_id == 2) {
      copySubBlockToMatrix(mat,msize,row,col,jinv_h,&mmat_size);
   }
   else if (*mat_id == 3) {
      copySubBlockToMatrix(mat,msize,row,col,gij_h,&mmat_size);
   }
   else {
      fprintf(stderr,"\nError: invalid matrix ID, %d\n",*mat_id);
      exit(EXIT_FAILURE);
   }
}

extern "C"
void push_gij_matrix_gpu_(int *cant, int *row, int *col, double _Complex *gij, int *kkri) {
   if (*cant != 1 && *cant != 2) {
      fprintf(stderr,"\nError: cant <> 1 and 2, %d\n",*cant);
      exit(EXIT_FAILURE);
   } 
   else if (*kkri * *cant != tau_size) {
      fprintf(stderr,"\nError: kkri*cant <> tau_size, %d,%d\n",*kkri * *cant,tau_size);
      exit(EXIT_FAILURE);
   }

   if (*cant == 1) {
      copySubBlockToMatrix(gij,kkri,row,col,gij_h,&mmat_size);
   }
   else {
      double _Complex *dij_h;
      dij_h = (double _Complex *) malloc(sizeof(double _Complex)*tau_size*tau_size);
      memset(dij_h, 0, sizeof(double _Complex)*tau_size*tau_size);
      int diag = 0;
      copySubBlockToMatrix(gij,kkri,&diag,&diag,dij_h,&tau_size);
      diag = 1;
      copySubBlockToMatrix(gij,kkri,&diag,&diag,dij_h,&tau_size);
      copySubBlockToMatrix(dij_h,&tau_size,row,col,gij_h,&mmat_size);
      free(dij_h);
   }
}

extern "C"
void commit_to_gpu_(int *mat_id) {
   size_t size = sizeof(cuDoubleComplex)*mmat_size*mmat_size;
   if (*mat_id == 1) {
      checkCudaErrors(cudaMemcpy(sine_d, sine_h, size, cudaMemcpyHostToDevice));
   }
   else if (*mat_id == 2) {
      checkCudaErrors(cudaMemcpy(jinv_d, jinv_h, size, cudaMemcpyHostToDevice));
   }
   else if (*mat_id == 3) {
      checkCudaErrors(cudaMemcpy(gij_d, gij_h, size, cudaMemcpyHostToDevice));
   }
   else {
      fprintf(stderr,"\nError: invalid matrix ID, %d\n",*mat_id);
      exit(EXIT_FAILURE);
   }
}

extern "C"
void construct_bigmatrix_gpu_(double _Complex *kappa) {
    cuDoubleComplex *jig_d;
    cuDoubleComplex minus_one = make_cuDoubleComplex(-1.0, 0.0);
    cuDoubleComplex one = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex zero = make_cuDoubleComplex(0.0, 0.0);
    double _Complex neg_kappa_inv = -1.0/(*kappa);
    cuDoubleComplex alpha = make_cuDoubleComplex(creal(neg_kappa_inv),cimag(neg_kappa_inv));
    // std::complex<double> neg_kappa_inv = -1.0/(*kappa);
    // cuDoubleComplex alpha = make_cuDoubleComplex(neg_kappa_inv.real(),neg_kappa_inv.imag());

    size_t size = mmat_size * mmat_size * sizeof(cuDoubleComplex);

    // Allocate device memory
    checkCudaErrors(cudaMalloc((void**)&jig_d, size));

    // Create cuBLAS handle
    cublasHandle_t handle;
    checkCublasErrors(cublasCreate(&handle));

    // Compute jig = jinv * gij
    checkCublasErrors(cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                  mmat_size, mmat_size, mmat_size, &one,
                                  jinv_d, mmat_size,
                                  gij_d, mmat_size,
                                  &zero, jig_d, mmat_size));

    // Compute BigMat = 1 - jig * sine / kappa
    checkCublasErrors(cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                  mmat_size, mmat_size, mmat_size, &alpha,
                                  jig_d, mmat_size,
                                  sine_d, mmat_size,
                                  &one, BigMat_d, mmat_size));

    // Cleanup
    cudaFree(jig_d);
    cublasDestroy(handle);
}

extern "C"
void invert_bigmatrix_gpu_(double _Complex *block, int *block_size) {
   // static cudaError_t error;
   static int Lwork;

   if (tau_size != *block_size) {
      fprintf(stderr,"\nError: tau_size <> block_size, %d,%d\n",tau_size,*block_size);
      exit(EXIT_FAILURE);
   }

   cusolverDnHandle_t cusolverHandle;
   cusolverStatus_t cusolverStatus;
   cusolverDnCreate(&cusolverHandle);
   cusolverDnZgetrf_bufferSize(cusolverHandle, mmat_size, mmat_size, BigMat_d, mmat_size, &Lwork);

   // Create a working array on the device
   // ========================================
   cuDoubleComplex  *workArray;
   checkCudaErrors(cudaMalloc((void**)&workArray, Lwork*sizeof(cuDoubleComplex)));

   // Define kernel launch parameters and create a unit matrix on device
   // ========================================
   int threads_per_block = 512;
   int num_blocks = (mmat_size*mmat_size + threads_per_block - 1) / threads_per_block;
   createUnitMatrixKernel<<<num_blocks, threads_per_block>>>(BigMatInv_d, mmat_size);
   checkCudaErrors(cudaPeekAtLastError());

   // Perform matrix inverse on the device
   // ============================================
   cusolverStatus = cusolverDnZgetrf(cusolverHandle, mmat_size, mmat_size, BigMat_d, mmat_size, 
                                     workArray, pivotArray, infoArray);
   if (cusolverStatus != CUSOLVER_STATUS_SUCCESS) {
      fprintf(stderr,"\nError: cuSOLVER ZGETRF UNSUCCESSFUL! \n");
   }
   cusolverStatus = cusolverDnZgetrs(cusolverHandle, CUBLAS_OP_N, mmat_size, mmat_size, BigMat_d, mmat_size, 
                                     pivotArray, BigMatInv_d, mmat_size, infoArray); 
   if (cusolverStatus != CUSOLVER_STATUS_SUCCESS) {
      fprintf(stderr,"\nError: cuSOLVER ZGETRS UNSUCCESSFUL! \n");
   }

   // --------------------------------------------
   // Aggregate the data into block_d array and then make one cudaMemcpy call
   // --------------------------------------------
   // Define grid and block dimensions for the kernel
   dim3 threadsPerBlock(32, 32); // Using a 32x32 block
   dim3 blocksPerGrid((tau_size + threadsPerBlock.x - 1) / threadsPerBlock.x,
                      (tau_size + threadsPerBlock.y - 1) / threadsPerBlock.y);

   // Launch the kernel
   // ============================================
   copyBlockKernel<<<blocksPerGrid, threadsPerBlock>>>(BigMatInv_d, mmat_size, block_d, tau_size);
   checkCudaErrors(cudaPeekAtLastError());
   checkCudaErrors(cudaMemcpy(block, block_d, sizeof(cuDoubleComplex)*tau_size*tau_size, 
                              cudaMemcpyDeviceToHost));

   // clean up
   // ============================================
   cusolverDnDestroy(cusolverHandle);
   cudaFree(workArray);
}
