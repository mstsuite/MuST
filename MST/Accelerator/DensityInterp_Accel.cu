/* *********************************************************************
 *  DensityInterp_Accel.cu
 *
 *  GPU acceleration of the radial-interpolation + spherical-harmonic
 *  evaluation of an L-expanded density (or any L-expanded radial
 *  function) on a structured (radial node) x (angular direction) grid.
 *
 *  This replaces the per-point CPU work performed by
 *      ChargeDensityModule::getChargeDensityAtPoint
 *  when it is called in the inner loop of
 *      PotentialGenerationModule::calExchangeJl
 *  where the evaluation points are the direct product
 *      r_i * u_g ,   i = 1..n_r ,  g = 1..n_g
 *  On a typical full-potential run n_r ~ 1500 and n_g = n_theta*n_phi
 *  = 50*80 = 4000, i.e. ~6x10^6 point evaluations per (atom, species),
 *  each of which performed a bisection search (hunt) plus a 5-point
 *  Neville interpolation for every jl component.
 *
 *  Mathematical content (identical to the CPU path):
 *
 *    rho(r_i, u_g) = sum_jl  fa2(jl) * Re[ A(i,jl) * Y_jl(u_g) ]
 *
 *    A(i,jl)       = sum_k    w_k(r_i) * rhoL( irp(i)+k-1, jl )
 *
 *  where w_k are the Lagrange weights of the n_inter-point polynomial
 *  through the mesh nodes irp(i) .. irp(i)+n_inter-1 and fa2(jl) = 1
 *  for m = 0 and 2 for m /= 0.  The Lagrange form is algebraically the
 *  same polynomial as the Neville recursion used by polint_inline, and
 *  it is exact (bitwise, since one numerator factor is exactly zero)
 *  when the target radius coincides with a mesh node -- which is the
 *  case for every call made from calExchangeJl.
 *
 *  Because the interpolation is a *linear* operator on rhoL, the
 *  angular sum collapses to two real matrix-matrix products:
 *
 *    rho(n_r x n_g) = Ar * Wr - Ai * Wi
 *
 *    Wr(jl,g) = fa2(jl)*Re[Y_jl(u_g)] ,  Wi(jl,g) = fa2(jl)*Im[Y_jl(u_g)]
 *
 *  Wr and Wi depend only on the (fixed) angular grid and are built once
 *  at initialisation.  The GEMMs are performed by cuBLAS.
 *
 *  Fortran interface (implicit interface, trailing underscore, every
 *  argument passed by reference -- same convention as LSMS_Accel.cu):
 *
 *    call init_density_interp_gpu(nr_max, jmax_max, ng, kmax, n_inter,
 *                                n_fields, my_pe)
 *    call push_angular_ylm_gpu(ylm, ng, kmax, kofj, fa2, jmax)
 *    call push_radial_mesh_gpu(r_mesh, nr)
 *    call push_density_l_gpu(rhoL, nr, jmax, ifield)
 *    call eval_density_sphere_gpu(r_target, n_target, jmax, ifield, rho)
 *    call finalize_density_interp_gpu()
 *
 *  Book-keeping:      date                       note
 *                   08-24-26     Initial version (B1 hotspot offload)
 * ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <complex.h>
#include "cuComplex.h"
#include "acclib.hpp"

/*  Maximum interpolation order kept in a register array.  The CPU side
 *  uses n_inter = 5 (ChargeDensityModule); 10 matches the nmax of
 *  polint_inline.                                                    */
#define DI_MAX_N_INTER 10
#define DI_MAX_FIELDS   4

/* -------------------------------------------------------------------
 *  Persistent device state
 * ------------------------------------------------------------------- */
static bool   di_initialized = false;
static int    di_nr_max   = 0;   /* max number of radial mesh points   */
static int    di_jmax_max = 0;   /* max number of (l,m>=0) components  */
static int    di_ng       = 0;   /* number of angular directions       */
static int    di_kmax     = 0;   /* (lmax+1)^2                         */
static int    di_n_inter  = 5;   /* interpolation order                */
static int    di_n_fields = 2;   /* e.g. 1 = charge, 2 = moment        */
static int    di_my_rank  = -1;

static bool   di_ylm_pushed  = false;
static bool   di_mesh_pushed = false;
static int    di_nr_pushed   = 0;

static double         *d_r_mesh = NULL;              /* (nr_max)            */
static double         *d_Wr     = NULL;              /* (jmax_max x ng)     */
static double         *d_Wi     = NULL;              /* (jmax_max x ng)     */
static cuDoubleComplex *d_ylm   = NULL;              /* (ng x kmax)         */
static int            *d_kofj   = NULL;              /* (jmax_max)          */
static double         *d_fa2    = NULL;              /* (jmax_max)          */
static cuDoubleComplex *d_rhoL[DI_MAX_FIELDS];       /* (nr_max x jmax_max) */
static double         *d_Ar     = NULL;              /* (nr_max x jmax_max) */
static double         *d_Ai     = NULL;              /* (nr_max x jmax_max) */
static double         *d_rho    = NULL;              /* (nr_max x ng)       */
static double         *d_rt     = NULL;              /* (nr_max)            */
static int            *d_irp    = NULL;              /* (nr_max)            */
static double         *d_wlag   = NULL;              /* (nr_max x n_inter)  */

/* ---- gradient (GGA) support ------------------------------------- *
 *  grad_ylm as returned by calYlm depends on the radius only through
 *  an overall 1/r factor (rfac = clm(jl)/r in SphericalHarmonics4), so
 *  its angular part can be tabulated once on the unit sphere exactly
 *  like Wr/Wi.  d_Gr[c]/d_Gi[c] hold fa2(jl)*Re/Im of that angular
 *  part for Cartesian component c.                                   */
static double         *d_Gr[3]  = {NULL,NULL,NULL};  /* (jmax_max x ng)     */
static double         *d_Gi[3]  = {NULL,NULL,NULL};  /* (jmax_max x ng)     */
static cuDoubleComplex *d_grady = NULL;              /* (ng x kmax x 3)     */
static double         *d_upos   = NULL;              /* (3 x ng)            */
static double         *d_D      = NULL;              /* (nr_max x ng)       */
static double         *d_T      = NULL;              /* (nr_max x ng)       */
static bool            di_grad_ready = false;

static cudaStream_t    di_stream;
static cublasHandle_t  di_cublas;

/* -------------------------------------------------------------------
 *  Device helpers
 * ------------------------------------------------------------------- */

/*  Bracket search reproducing the semantics of MST/lib/hunt.F90 for a
 *  monotonically increasing mesh: returns jlo (1-based) such that
 *  x[jlo-1] <= xt <= x[jlo].  Out-of-range targets are clamped to the
 *  end intervals, matching hunt's x==xx(1) -> 1, x==xx(n) -> n-1
 *  special cases.                                                    */
__device__ __forceinline__
int di_bracket(const double * __restrict__ x, int n, double xt)
{
   if (n < 2) return 1;
   if (xt <= x[0])   return 1;
   if (xt >= x[n-1]) return n-1;

   int lo = 1, hi = n;                 /* 1-based, x[lo-1] <= xt < x[hi-1] */
   while (hi - lo > 1) {
      int mid = (hi + lo) >> 1;
      if (xt >= x[mid-1]) lo = mid;
      else                hi = mid;
   }
   return lo;
}

/* -------------------------------------------------------------------
 *  Kernel 1: interpolation stencil offset + Lagrange weights
 *
 *  Reproduces exactly the CPU clamping logic of
 *  getChargeDensityAtPoint:
 *      if (ir > iend-(n_inter-1)/2)  irp = iend-n_inter+1
 *      else if (2*ir+1 > n_inter)    irp = ir-(n_inter-1)/2
 *      else                          irp = 1
 * ------------------------------------------------------------------- */
__global__
void di_buildLagrangeWeightsKernel(const double * __restrict__ r_mesh, int nr,
                                   const double * __restrict__ rt,     int nt,
                                   int n_inter,
                                   int    * __restrict__ irp_out,
                                   double * __restrict__ w_out)
{
   for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < nt;
            i += blockDim.x*gridDim.x) {

      const double x    = rt[i];
      const int    ir   = di_bracket(r_mesh, nr, x);
      const int    half = (n_inter-1)/2;

      int irp;
      if      (ir > nr - half)      irp = nr - n_inter + 1;
      else if (2*ir + 1 > n_inter)  irp = ir - half;
      else                          irp = 1;

      if (irp < 1)                irp = 1;
      if (irp > nr - n_inter + 1) irp = nr - n_inter + 1;
      irp_out[i] = irp;

      double xk[DI_MAX_N_INTER];
      for (int k = 0; k < n_inter; ++k) xk[k] = r_mesh[irp-1+k];

      /*  L_k(x) = prod_{m/=k} (x-x_m)/(x_k-x_m).  When x == x_j exactly
       *  the numerator of every k/=j contains an exact 0.0 and the j-th
       *  ratio is a quotient of two identical products, i.e. exactly
       *  1.0 -- so on-node targets return the node value bitwise.     */
      for (int k = 0; k < n_inter; ++k) {
         double num = 1.0, den = 1.0;
         for (int m = 0; m < n_inter; ++m) {
            if (m == k) continue;
            num *= (x     - xk[m]);
            den *= (xk[k] - xk[m]);
         }
         w_out[i + (size_t)k*nt] = num/den;
      }
   }
}

/* -------------------------------------------------------------------
 *  Kernel 2: radial interpolation of every L component,
 *            split into real and imaginary parts for the GEMMs.
 *
 *      A(i,jl) = sum_k w_k(i) * rhoL(irp(i)+k-1, jl)
 * ------------------------------------------------------------------- */
__global__
void di_interpolateDensityLKernel(const cuDoubleComplex * __restrict__ rhoL,
                                  int ld_rho, int nt, int jmax, int n_inter,
                                  const int    * __restrict__ irp,
                                  const double * __restrict__ w,
                                  double * __restrict__ Ar,
                                  double * __restrict__ Ai,
                                  int ld_A)
{
   const int total = nt*jmax;
   for (int idx = blockIdx.x*blockDim.x + threadIdx.x; idx < total;
            idx += blockDim.x*gridDim.x) {

      const int i    = idx % nt;
      const int jl   = idx / nt;
      const int base = irp[i] - 1;

      double sr = 0.0, si = 0.0;
      for (int k = 0; k < n_inter; ++k) {
         const double          wk = w[i + (size_t)k*nt];
         const cuDoubleComplex v  = rhoL[(base + k) + (size_t)jl*ld_rho];
         sr = fma(wk, cuCreal(v), sr);
         si = fma(wk, cuCimag(v), si);
      }
      Ar[i + (size_t)jl*ld_A] = sr;
      Ai[i + (size_t)jl*ld_A] = si;
   }
}

/* -------------------------------------------------------------------
 *  Kernel 3: pack the fixed angular weight tables
 *
 *      Wr(jl,g) = fa2(jl)*Re[Y_{kofj(jl)}(u_g)]
 *      Wi(jl,g) = fa2(jl)*Im[Y_{kofj(jl)}(u_g)]
 *
 *  ylm arrives in the Fortran layout ylm(ng,kmax) (column major).
 * ------------------------------------------------------------------- */
__global__
void di_packYlmWeightsKernel(const cuDoubleComplex * __restrict__ ylm,
                             int ng, int kmax, int jmax,
                             const int    * __restrict__ kofj,
                             const double * __restrict__ fa2,
                             double * __restrict__ Wr,
                             double * __restrict__ Wi)
{
   const int total = jmax*ng;
   for (int idx = blockIdx.x*blockDim.x + threadIdx.x; idx < total;
            idx += blockDim.x*gridDim.x) {

      const int jl = idx % jmax;
      const int g  = idx / jmax;
      const int kl = kofj[jl] - 1;            /* Fortran 1-based -> C */

      if (kl < 0 || kl >= kmax) {
         Wr[jl + (size_t)g*jmax] = 0.0;
         Wi[jl + (size_t)g*jmax] = 0.0;
         continue;
      }
      const cuDoubleComplex y = ylm[g + (size_t)kl*ng];
      const double          f = fa2[jl];
      Wr[jl + (size_t)g*jmax] = f*cuCreal(y);
      Wi[jl + (size_t)g*jmax] = f*cuCimag(y);
   }
}

/* -------------------------------------------------------------------
 *  Kernel 4: pack the fixed angular gradient tables
 *
 *      Gr(jl,g,c) = fa2(jl)*Re[ r * d/dx_c Y_{kofj(jl)}(u_g) ]
 *      Gi(jl,g,c) = fa2(jl)*Im[ ... ]
 *
 *  The host evaluates calYlm on the UNIT sphere, so the 1/r factor
 *  present in grad_ylm equals 1 there and the stored table is the pure
 *  angular part; the 1/r_i scaling is reapplied per radius in
 *  di_assembleGradKernel below.
 *
 *  grady arrives in the Fortran layout grady(ng,kmax,3) (column major).
 * ------------------------------------------------------------------- */
__global__
void di_packGradYlmWeightsKernel(const cuDoubleComplex * __restrict__ grady,
                                 int ng, int kmax, int jmax, int c,
                                 const int    * __restrict__ kofj,
                                 const double * __restrict__ fa2,
                                 double * __restrict__ Gr,
                                 double * __restrict__ Gi)
{
   const int total = jmax*ng;
   const size_t coff = (size_t)c*(size_t)ng*(size_t)kmax;

   for (int idx = blockIdx.x*blockDim.x + threadIdx.x; idx < total;
            idx += blockDim.x*gridDim.x) {

      const int jl = idx % jmax;
      const int g  = idx / jmax;
      const int kl = kofj[jl] - 1;

      if (kl < 0 || kl >= kmax) {
         Gr[jl + (size_t)g*jmax] = 0.0;
         Gi[jl + (size_t)g*jmax] = 0.0;
         continue;
      }
      const cuDoubleComplex y = grady[coff + g + (size_t)kl*ng];
      const double          f = fa2[jl];
      Gr[jl + (size_t)g*jmax] = f*cuCreal(y);
      Gi[jl + (size_t)g*jmax] = f*cuCimag(y);
   }
}

/* -------------------------------------------------------------------
 *  Kernel 5: assemble one Cartesian component of the gradient
 *
 *      grad_c(i,g) = D(i,g)*u_c(g) + T_c(i,g)/r_i
 *
 *  D is the radial-derivative contraction (same weight tables as the
 *  value) and T_c the angular one.  Reproduces exactly the CPU line
 *
 *      grad(i) += fa2(jl)*Re[ der_rho_in*ylm(kl)*er(i)
 *                           + rho_in(jl)*grad_ylm(kl,i) ]
 *
 *  of getChargeDensityAtPoint, with er = u_g.  The result overwrites T.
 * ------------------------------------------------------------------- */
__global__
void di_assembleGradKernel(const double * __restrict__ D,
                           double * __restrict__ T,
                           const double * __restrict__ upos,
                           const double * __restrict__ r_mesh,
                           int nt, int ng, int ld, int c)
{
   const int total = nt*ng;
   for (int idx = blockIdx.x*blockDim.x + threadIdx.x; idx < total;
            idx += blockDim.x*gridDim.x) {

      const int i = idx % nt;
      const int g = idx / nt;

      const double uc   = upos[c + 3*(size_t)g];
      const double rinv = 1.0/r_mesh[i];
      const size_t o    = (size_t)i + (size_t)g*ld;

      T[o] = D[o]*uc + T[o]*rinv;
   }
}

/* -------------------------------------------------------------------
 *  Host-side helpers
 * ------------------------------------------------------------------- */
static inline int di_blocks(int n, int tpb)
{
   int b = (n + tpb - 1)/tpb;
   if (b < 1)     b = 1;
   if (b > 65535) b = 65535;      /* grid-stride loops cover the rest  */
   return b;
}

/* ===================================================================
 *  init_density_interp_gpu(nr_max, jmax_max, ng, kmax, n_inter,
 *                          n_fields, my_pe)
 * =================================================================== */
/* ===================================================================
 *  query_density_interp_gpu(bytes_needed, ok)
 *
 *  Runtime gate.  This module previously had NONE: DensityOnGridModule
 *  set useGPU from #ifdef ACCEL alone, so an ACCEL binary died at the
 *  first cudaMalloc on a node with no usable GPU, and oversubscribing
 *  one GPU with several ranks aborted instead of falling back.  The
 *  host DGEMM path was always there; only the decision to take it was
 *  missing.
 *
 *  ok = 1 only if a device exists, can be selected, and has
 *  bytes_needed free.  Every failure clears the CUDA error state so a
 *  failed probe cannot poison later CUDA use.  MUST_DENSITY_GPU=0
 *  forces ok = 0 without touching the driver; the -cpu/--cpu-only
 *  command-line flag takes precedence over that and is handled on the
 *  Fortran side, which simply does not call this.
 * =================================================================== */
extern "C"
void query_density_interp_gpu_(long *bytes_needed, int *ok)
{
   *ok = 0;

   const char *env = getenv("MUST_DENSITY_GPU");
   if (env != NULL && (env[0] == '0' || env[0] == 'n' || env[0] == 'N' ||
                       env[0] == 'f' || env[0] == 'F')) {
      return;
   }

   int ndev = 0;
   if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev < 1) {
      cudaGetLastError();  return;
   }
   int dev = 0;
   if (cudaGetDevice(&dev) != cudaSuccess) { cudaGetLastError(); return; }
   if (cudaSetDevice(dev) != cudaSuccess)  { cudaGetLastError(); return; }

   size_t freeb = 0, totb = 0;
   if (cudaMemGetInfo(&freeb, &totb) != cudaSuccess) {
      cudaGetLastError();  return;
   }
   /*  Head room on top of this module's own buffers, for the cuBLAS
    *  work space and context.                                        */
   const long slack = 256L*1024L*1024L;
   if (*bytes_needed > 0 && (size_t)(*bytes_needed + slack) > freeb) return;

   *ok = 1;
}

extern "C"
void init_density_interp_gpu_(int *nr_max, int *jmax_max, int *ng, int *kmax,
                              int *n_inter, int *n_fields, int *my_pe)
{
   if (di_initialized) return;

   di_nr_max   = *nr_max;
   di_jmax_max = *jmax_max;
   di_ng       = *ng;
   di_kmax     = *kmax;
   di_n_inter  = *n_inter;
   di_n_fields = *n_fields;
   di_my_rank  = *my_pe;

   if (di_nr_max < 2 || di_jmax_max < 1 || di_ng < 1 || di_kmax < 1) {
      fprintf(stderr, "\nError in init_density_interp_gpu: bad dimensions "
                      "nr_max=%d jmax_max=%d ng=%d kmax=%d\n",
                      di_nr_max, di_jmax_max, di_ng, di_kmax);
      exit(EXIT_FAILURE);
   }
   if (di_n_inter < 2 || di_n_inter > DI_MAX_N_INTER) {
      fprintf(stderr, "\nError in init_density_interp_gpu: n_inter=%d out of "
                      "range [2,%d]\n", di_n_inter, DI_MAX_N_INTER);
      exit(EXIT_FAILURE);
   }
   if (di_n_fields < 1 || di_n_fields > DI_MAX_FIELDS) {
      fprintf(stderr, "\nError in init_density_interp_gpu: n_fields=%d out of "
                      "range [1,%d]\n", di_n_fields, DI_MAX_FIELDS);
      exit(EXIT_FAILURE);
   }

   /*  Assign a device per MPI rank, exactly as init_lsms_gpu_ does, so
    *  that the two accelerator modules share one context per rank.    */
   int ngpus = 0;
   checkCudaErrors(cudaGetDeviceCount(&ngpus));
   if (ngpus < 1) {
      fprintf(stderr, "\nError in init_density_interp_gpu: no GPU found\n");
      exit(EXIT_FAILURE);
   }
   if (di_my_rank >= 0) {
      checkCudaErrors(cudaSetDevice(di_my_rank % ngpus));
   }

   checkCudaErrors(cudaStreamCreate(&di_stream));
   checkCublasErrors(cublasCreate(&di_cublas));
   checkCublasErrors(cublasSetStream(di_cublas, di_stream));

   const size_t sz_r    = (size_t)di_nr_max*sizeof(double);
   const size_t sz_W    = (size_t)di_jmax_max*(size_t)di_ng*sizeof(double);
   const size_t sz_ylm  = (size_t)di_ng*(size_t)di_kmax*sizeof(cuDoubleComplex);
   const size_t sz_rhoL = (size_t)di_nr_max*(size_t)di_jmax_max*sizeof(cuDoubleComplex);
   const size_t sz_A    = (size_t)di_nr_max*(size_t)di_jmax_max*sizeof(double);
   const size_t sz_rho  = (size_t)di_nr_max*(size_t)di_ng*sizeof(double);
   const size_t sz_w    = (size_t)di_nr_max*(size_t)di_n_inter*sizeof(double);

   checkCudaErrors(cudaMalloc((void**)&d_r_mesh, sz_r));
   checkCudaErrors(cudaMalloc((void**)&d_rt,     sz_r));
   checkCudaErrors(cudaMalloc((void**)&d_Wr,     sz_W));
   checkCudaErrors(cudaMalloc((void**)&d_Wi,     sz_W));
   checkCudaErrors(cudaMalloc((void**)&d_ylm,    sz_ylm));
   checkCudaErrors(cudaMalloc((void**)&d_kofj,   (size_t)di_jmax_max*sizeof(int)));
   checkCudaErrors(cudaMalloc((void**)&d_fa2,    (size_t)di_jmax_max*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_Ar,     sz_A));
   checkCudaErrors(cudaMalloc((void**)&d_Ai,     sz_A));
   checkCudaErrors(cudaMalloc((void**)&d_rho,    sz_rho));
   checkCudaErrors(cudaMalloc((void**)&d_irp,    (size_t)di_nr_max*sizeof(int)));
   checkCudaErrors(cudaMalloc((void**)&d_wlag,   sz_w));

   for (int f = 0; f < DI_MAX_FIELDS; ++f) d_rhoL[f] = NULL;
   for (int f = 0; f < di_n_fields; ++f) {
      checkCudaErrors(cudaMalloc((void**)&d_rhoL[f], sz_rhoL));
      checkCudaErrors(cudaMemset(d_rhoL[f], 0, sz_rhoL));
   }

   if (di_my_rank == 0) {
      double mb = (double)(sz_r*2 + sz_W*2 + sz_ylm + sz_A*2 + sz_rho + sz_w
                           + sz_rhoL*di_n_fields)/(1024.0*1024.0);
      printf("DensityInterp GPU initialized: nr_max=%d jmax_max=%d ng=%d "
             "kmax=%d n_inter=%d, device memory = %.2f MB\n",
             di_nr_max, di_jmax_max, di_ng, di_kmax, di_n_inter, mb);
   }

   di_ylm_pushed  = false;
   di_mesh_pushed = false;
   di_grad_ready  = false;
   di_initialized = true;
}

/* ===================================================================
 *  push_angular_gradylm_gpu(grady, upos, ng, kmax, kofj, fa2, jmax)
 *
 *  Optional companion to push_angular_ylm_gpu, required only for GGA.
 *  grady(ng,kmax,3) must be evaluated on the UNIT sphere (|u_g| = 1) so
 *  that the 1/r factor carried by grad_ylm is unity; the radius scaling
 *  is applied per evaluation radius on the device.
 *
 *  Allocates the gradient work space on first call, so a run that never
 *  uses GGA pays neither the memory nor the transfer.
 * =================================================================== */
extern "C"
void push_angular_gradylm_gpu_(double _Complex *grady, double *upos,
                               int *ng, int *kmax, int *kofj, double *fa2,
                               int *jmax)
{
   if (!di_initialized) {
      fprintf(stderr, "\nError in push_angular_gradylm_gpu: "
                      "init_density_interp_gpu must be called first\n");
      exit(EXIT_FAILURE);
   }
   if (*ng != di_ng || *kmax != di_kmax || *jmax > di_jmax_max) {
      fprintf(stderr, "\nError in push_angular_gradylm_gpu: dimension "
                      "mismatch (ng %d/%d, kmax %d/%d, jmax %d/%d)\n",
                      *ng, di_ng, *kmax, di_kmax, *jmax, di_jmax_max);
      exit(EXIT_FAILURE);
   }

   const size_t sz_W     = (size_t)di_jmax_max*(size_t)di_ng*sizeof(double);
   const size_t sz_grady = (size_t)di_ng*(size_t)di_kmax*3*sizeof(cuDoubleComplex);
   const size_t sz_upos  = 3*(size_t)di_ng*sizeof(double);
   const size_t sz_rho   = (size_t)di_nr_max*(size_t)di_ng*sizeof(double);

   if (!di_grad_ready) {
      for (int c = 0; c < 3; ++c) {
         checkCudaErrors(cudaMalloc((void**)&d_Gr[c], sz_W));
         checkCudaErrors(cudaMalloc((void**)&d_Gi[c], sz_W));
      }
      checkCudaErrors(cudaMalloc((void**)&d_grady, sz_grady));
      checkCudaErrors(cudaMalloc((void**)&d_upos,  sz_upos));
      checkCudaErrors(cudaMalloc((void**)&d_D,     sz_rho));
      checkCudaErrors(cudaMalloc((void**)&d_T,     sz_rho));

      if (di_my_rank == 0) {
         double mb = (double)(sz_W*6 + sz_grady + sz_upos + sz_rho*2)
                     /(1024.0*1024.0);
         printf("DensityInterp GPU: gradient (GGA) tables allocated, "
                "additional device memory = %.2f MB\n", mb);
      }
   }

   checkCudaErrors(cudaMemcpyAsync(d_grady, grady, sz_grady,
                   cudaMemcpyHostToDevice, di_stream));
   checkCudaErrors(cudaMemcpyAsync(d_upos, upos, sz_upos,
                   cudaMemcpyHostToDevice, di_stream));

   const int tpb = 256;
   const int n   = (*jmax)*(*ng);
   for (int c = 0; c < 3; ++c) {
      di_packGradYlmWeightsKernel<<<di_blocks(n,tpb), tpb, 0, di_stream>>>
                        (d_grady, *ng, *kmax, *jmax, c, d_kofj, d_fa2,
                         d_Gr[c], d_Gi[c]);
      checkCudaErrors(cudaPeekAtLastError());
   }
   checkCudaErrors(cudaStreamSynchronize(di_stream));

   di_grad_ready = true;
}

/* ===================================================================
 *  push_angular_ylm_gpu(ylm, ng, kmax, kofj, fa2, jmax)
 *
 *  The angular grid is fixed for the whole run (AngularIntegrationModule
 *  builds it once in setAngularData), so this is called only once and
 *  the packed weight tables are reused by every evaluation.
 * =================================================================== */
extern "C"
void push_angular_ylm_gpu_(double _Complex *ylm, int *ng, int *kmax,
                           int *kofj, double *fa2, int *jmax)
{
   if (!di_initialized) {
      fprintf(stderr, "\nError in push_angular_ylm_gpu: "
                      "init_density_interp_gpu must be called first\n");
      exit(EXIT_FAILURE);
   }
   if (*ng != di_ng || *kmax != di_kmax || *jmax > di_jmax_max) {
      fprintf(stderr, "\nError in push_angular_ylm_gpu: dimension mismatch "
                      "(ng %d/%d, kmax %d/%d, jmax %d/%d)\n",
                      *ng, di_ng, *kmax, di_kmax, *jmax, di_jmax_max);
      exit(EXIT_FAILURE);
   }

   checkCudaErrors(cudaMemcpyAsync(d_ylm, ylm,
                   (size_t)(*ng)*(size_t)(*kmax)*sizeof(cuDoubleComplex),
                   cudaMemcpyHostToDevice, di_stream));
   checkCudaErrors(cudaMemcpyAsync(d_kofj, kofj,
                   (size_t)(*jmax)*sizeof(int),
                   cudaMemcpyHostToDevice, di_stream));
   checkCudaErrors(cudaMemcpyAsync(d_fa2, fa2,
                   (size_t)(*jmax)*sizeof(double),
                   cudaMemcpyHostToDevice, di_stream));

   const int tpb = 256;
   const int n   = (*jmax)*(*ng);
   di_packYlmWeightsKernel<<<di_blocks(n,tpb), tpb, 0, di_stream>>>
                          (d_ylm, *ng, *kmax, *jmax, d_kofj, d_fa2, d_Wr, d_Wi);
   checkCudaErrors(cudaPeekAtLastError());
   checkCudaErrors(cudaStreamSynchronize(di_stream));

   di_ylm_pushed = true;
}

/* ===================================================================
 *  push_radial_mesh_gpu(r_mesh, nr)
 * =================================================================== */
extern "C"
void push_radial_mesh_gpu_(double *r_mesh, int *nr)
{
   if (!di_initialized) {
      fprintf(stderr, "\nError in push_radial_mesh_gpu: "
                      "init_density_interp_gpu must be called first\n");
      exit(EXIT_FAILURE);
   }
   if (*nr < di_n_inter || *nr > di_nr_max) {
      fprintf(stderr, "\nError in push_radial_mesh_gpu: nr=%d outside "
                      "[n_inter=%d, nr_max=%d]\n", *nr, di_n_inter, di_nr_max);
      exit(EXIT_FAILURE);
   }
   checkCudaErrors(cudaMemcpyAsync(d_r_mesh, r_mesh,
                   (size_t)(*nr)*sizeof(double),
                   cudaMemcpyHostToDevice, di_stream));
   di_nr_pushed   = *nr;
   di_mesh_pushed = true;
}

/* ===================================================================
 *  push_density_l_gpu(rhoL, ld_rho, nr, jmax, ifield)
 *
 *  rhoL is the Fortran array rhoL(ld_rho, jmax) (column major); only
 *  its leading nr rows are transferred.  ifield is 1-based:
 *  1 = charge density, 2 = moment density, ...
 * =================================================================== */
extern "C"
void push_density_l_gpu_(double _Complex *rhoL, int *ld_rho, int *nr,
                         int *jmax, int *ifield)
{
   if (!di_initialized) {
      fprintf(stderr, "\nError in push_density_l_gpu: "
                      "init_density_interp_gpu must be called first\n");
      exit(EXIT_FAILURE);
   }
   const int f = *ifield - 1;
   if (f < 0 || f >= di_n_fields) {
      fprintf(stderr, "\nError in push_density_l_gpu: invalid ifield=%d\n",
                      *ifield);
      exit(EXIT_FAILURE);
   }
   if (*nr > di_nr_max || *jmax > di_jmax_max || *ld_rho < *nr) {
      fprintf(stderr, "\nError in push_density_l_gpu: dimension problem "
                      "(nr %d/%d, jmax %d/%d, ld_rho %d)\n",
                      *nr, di_nr_max, *jmax, di_jmax_max, *ld_rho);
      exit(EXIT_FAILURE);
   }

   /*  The host array has leading dimension ld_rho while the device
    *  buffer has leading dimension nr_max, so this is a strided (2-D)
    *  copy of nr rows x jmax columns.                                */
   checkCudaErrors(cudaMemcpy2DAsync(
                      d_rhoL[f], (size_t)di_nr_max*sizeof(cuDoubleComplex),
                      rhoL,      (size_t)(*ld_rho)*sizeof(cuDoubleComplex),
                      (size_t)(*nr)*sizeof(cuDoubleComplex),
                      (size_t)(*jmax),
                      cudaMemcpyHostToDevice, di_stream));
}

/* ===================================================================
 *  eval_density_sphere_gpu(n_target, jmax, ifield, rho, ld_rho_out)
 *
 *  The evaluation targets are the first n_target nodes of the radial
 *  mesh already resident on the device (the caller's grid points), so
 *  no target list has to be transferred.
 *
 *  Returns rho(ld_rho_out, ng) in the Fortran (column-major) layout:
 *      rho(i,g) = sum_jl fa2(jl)*Re[ A(i,jl)*Y_jl(u_g) ]
 * =================================================================== */
extern "C"
void eval_density_sphere_gpu_(int *n_target, int *jmax, int *ifield,
                              double *rho, int *ld_rho_out)
{
   if (!di_initialized || !di_ylm_pushed || !di_mesh_pushed) {
      fprintf(stderr, "\nError in eval_density_sphere_gpu: module state "
                      "incomplete (init=%d ylm=%d mesh=%d)\n",
                      (int)di_initialized, (int)di_ylm_pushed,
                      (int)di_mesh_pushed);
      exit(EXIT_FAILURE);
   }
   const int f   = *ifield - 1;
   const int nt  = *n_target;
   const int jm  = *jmax;
   const int ldo = *ld_rho_out;

   if (f < 0 || f >= di_n_fields || d_rhoL[f] == NULL) {
      fprintf(stderr, "\nError in eval_density_sphere_gpu: invalid ifield=%d\n",
                      *ifield);
      exit(EXIT_FAILURE);
   }
   if (nt < 1 || nt > di_nr_pushed || jm < 1 || jm > di_jmax_max || ldo < nt) {
      fprintf(stderr, "\nError in eval_density_sphere_gpu: bad sizes "
                      "(n_target %d/%d, jmax %d/%d, ld_out %d)\n",
                      nt, di_nr_pushed, jm, di_jmax_max, ldo);
      exit(EXIT_FAILURE);
   }

   const int tpb = 256;

   /* --- stencil offsets and Lagrange weights ----------------------
    *  Targets are the mesh nodes themselves (d_r_mesh), which makes
    *  the weights exactly Kronecker deltas; the general code path is
    *  retained so that the same kernel serves off-node targets.      */
   di_buildLagrangeWeightsKernel<<<di_blocks(nt,tpb), tpb, 0, di_stream>>>
                    (d_r_mesh, di_nr_pushed, d_r_mesh, nt, di_n_inter,
                     d_irp, d_wlag);
   checkCudaErrors(cudaPeekAtLastError());

   /* --- radial interpolation of every L component ----------------- */
   const int n2 = nt*jm;
   di_interpolateDensityLKernel<<<di_blocks(n2,tpb), tpb, 0, di_stream>>>
                    (d_rhoL[f], di_nr_max, nt, jm, di_n_inter,
                     d_irp, d_wlag, d_Ar, d_Ai, di_nr_max);
   checkCudaErrors(cudaPeekAtLastError());

   /* --- angular sum:  rho = Ar*Wr - Ai*Wi -------------------------
    *  Ar, Ai are (nt x jm) with leading dimension nr_max on the
    *  device; Wr, Wi are (jm x ng) with leading dimension jmax_max.  */
   const double one = 1.0, mone = -1.0, zero = 0.0;
   checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                 nt, di_ng, jm,
                                 &one,  d_Ar, di_nr_max,
                                        d_Wr, di_jmax_max,
                                 &zero, d_rho, di_nr_max));
   checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                 nt, di_ng, jm,
                                 &mone, d_Ai, di_nr_max,
                                        d_Wi, di_jmax_max,
                                 &one,  d_rho, di_nr_max));

   /*  Device buffer is (nr_max x ng); the host array is (ldo x ng).  */
   checkCudaErrors(cudaMemcpy2DAsync(
                      rho,   (size_t)ldo*sizeof(double),
                      d_rho, (size_t)di_nr_max*sizeof(double),
                      (size_t)nt*sizeof(double), (size_t)di_ng,
                      cudaMemcpyDeviceToHost, di_stream));
   checkCudaErrors(cudaStreamSynchronize(di_stream));
}

/* ===================================================================
 *  eval_density_grad_sphere_gpu(n_target, jmax, ifield_val, ifield_der,
 *                               rho, grad, ld_out)
 *
 *  GGA companion of eval_density_sphere_gpu.  Returns both the density
 *  and its Cartesian gradient on the (mesh node) x (direction) grid:
 *
 *      rho (ld_out, ng)      value
 *      grad(ld_out, ng, 3)   d(rho)/dx_c
 *
 *  ifield_val selects the buffer holding the L coefficients and
 *  ifield_der the buffer holding their radial derivatives; both must
 *  have been pushed with push_density_l_gpu.
 * =================================================================== */
extern "C"
void eval_density_grad_sphere_gpu_(int *n_target, int *jmax, int *ifield_val,
                                   int *ifield_der, double *rho, double *grad,
                                   int *ld_out)
{
   if (!di_initialized || !di_ylm_pushed || !di_mesh_pushed) {
      fprintf(stderr, "\nError in eval_density_grad_sphere_gpu: module state "
                      "incomplete (init=%d ylm=%d mesh=%d)\n",
                      (int)di_initialized, (int)di_ylm_pushed,
                      (int)di_mesh_pushed);
      exit(EXIT_FAILURE);
   }
   if (!di_grad_ready) {
      fprintf(stderr, "\nError in eval_density_grad_sphere_gpu: "
                      "push_angular_gradylm_gpu must be called first\n");
      exit(EXIT_FAILURE);
   }

   const int fv  = *ifield_val - 1;
   const int fd  = *ifield_der - 1;
   const int nt  = *n_target;
   const int jm  = *jmax;
   const int ldo = *ld_out;

   if (fv < 0 || fv >= di_n_fields || d_rhoL[fv] == NULL ||
       fd < 0 || fd >= di_n_fields || d_rhoL[fd] == NULL) {
      fprintf(stderr, "\nError in eval_density_grad_sphere_gpu: invalid "
                      "field indices (%d,%d)\n", *ifield_val, *ifield_der);
      exit(EXIT_FAILURE);
   }
   if (nt < 1 || nt > di_nr_pushed || jm < 1 || jm > di_jmax_max || ldo < nt) {
      fprintf(stderr, "\nError in eval_density_grad_sphere_gpu: bad sizes "
                      "(n_target %d/%d, jmax %d/%d, ld_out %d)\n",
                      nt, di_nr_pushed, jm, di_jmax_max, ldo);
      exit(EXIT_FAILURE);
   }

   const int    tpb  = 256;
   const double one = 1.0, mone = -1.0, zero = 0.0;
   const int    n2   = nt*jm;
   const int    nrg  = nt*di_ng;

   /* --- shared interpolation stencil (targets are the mesh nodes) -- */
   di_buildLagrangeWeightsKernel<<<di_blocks(nt,tpb), tpb, 0, di_stream>>>
                    (d_r_mesh, di_nr_pushed, d_r_mesh, nt, di_n_inter,
                     d_irp, d_wlag);
   checkCudaErrors(cudaPeekAtLastError());

   /* ============================================================== *
    *  (1) radial term first: D = sum_jl fa2*Re[ d(rho_jl)/dr * Y_jl ]
    *      from the DERIVATIVE coefficients, parked in d_D.  Doing this
    *      before the value pass means d_Ar/d_Ai can then be reused for
    *      the value without any device->host->device round trip.
    * ============================================================== */
   di_interpolateDensityLKernel<<<di_blocks(n2,tpb), tpb, 0, di_stream>>>
                    (d_rhoL[fd], di_nr_max, nt, jm, di_n_inter,
                     d_irp, d_wlag, d_Ar, d_Ai, di_nr_max);
   checkCudaErrors(cudaPeekAtLastError());

   checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                 nt, di_ng, jm, &one,  d_Ar, di_nr_max,
                                 d_Wr, di_jmax_max, &zero, d_D, di_nr_max));
   checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                 nt, di_ng, jm, &mone, d_Ai, di_nr_max,
                                 d_Wi, di_jmax_max, &one,  d_D, di_nr_max));

   /* ============================================================== *
    *  (2) value coefficients: rho = Ar*Wr - Ai*Wi
    * ============================================================== */
   di_interpolateDensityLKernel<<<di_blocks(n2,tpb), tpb, 0, di_stream>>>
                    (d_rhoL[fv], di_nr_max, nt, jm, di_n_inter,
                     d_irp, d_wlag, d_Ar, d_Ai, di_nr_max);
   checkCudaErrors(cudaPeekAtLastError());

   checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                 nt, di_ng, jm, &one,  d_Ar, di_nr_max,
                                 d_Wr, di_jmax_max, &zero, d_rho, di_nr_max));
   checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                 nt, di_ng, jm, &mone, d_Ai, di_nr_max,
                                 d_Wi, di_jmax_max, &one,  d_rho, di_nr_max));

   /* ============================================================== *
    *  (3) per component: angular term T_c from the same (value)
    *      coefficients against the gradient tables, then assemble
    *          grad_c = D*u_c + T_c/r
    *      entirely on the device and copy the finished component out.
    * ============================================================== */
   for (int c = 0; c < 3; ++c) {
      checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                    nt, di_ng, jm, &one,  d_Ar, di_nr_max,
                                    d_Gr[c], di_jmax_max, &zero,
                                    d_T, di_nr_max));
      checkCublasErrors(cublasDgemm(di_cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                    nt, di_ng, jm, &mone, d_Ai, di_nr_max,
                                    d_Gi[c], di_jmax_max, &one,
                                    d_T, di_nr_max));

      di_assembleGradKernel<<<di_blocks(nrg,tpb), tpb, 0, di_stream>>>
                        (d_D, d_T, d_upos, d_r_mesh, nt, di_ng,
                         di_nr_max, c);
      checkCudaErrors(cudaPeekAtLastError());

      checkCudaErrors(cudaMemcpy2DAsync(
                         grad + (size_t)c*(size_t)ldo*(size_t)di_ng,
                         (size_t)ldo*sizeof(double),
                         d_T, (size_t)di_nr_max*sizeof(double),
                         (size_t)nt*sizeof(double), (size_t)di_ng,
                         cudaMemcpyDeviceToHost, di_stream));
   }

   checkCudaErrors(cudaMemcpy2DAsync(
                      rho,   (size_t)ldo*sizeof(double),
                      d_rho, (size_t)di_nr_max*sizeof(double),
                      (size_t)nt*sizeof(double), (size_t)di_ng,
                      cudaMemcpyDeviceToHost, di_stream));
   checkCudaErrors(cudaStreamSynchronize(di_stream));
}

/* ===================================================================
 *  finalize_density_interp_gpu()
 * =================================================================== */
extern "C"
void finalize_density_interp_gpu_()
{
   if (!di_initialized) return;

   for (int f = 0; f < DI_MAX_FIELDS; ++f) {
      if (d_rhoL[f] != NULL) {
         checkCudaErrors(cudaFree(d_rhoL[f]));
         d_rhoL[f] = NULL;
      }
   }
   if (d_r_mesh) { checkCudaErrors(cudaFree(d_r_mesh)); d_r_mesh = NULL; }
   if (d_rt)     { checkCudaErrors(cudaFree(d_rt));     d_rt     = NULL; }
   if (d_Wr)     { checkCudaErrors(cudaFree(d_Wr));     d_Wr     = NULL; }
   if (d_Wi)     { checkCudaErrors(cudaFree(d_Wi));     d_Wi     = NULL; }
   if (d_ylm)    { checkCudaErrors(cudaFree(d_ylm));    d_ylm    = NULL; }
   if (d_kofj)   { checkCudaErrors(cudaFree(d_kofj));   d_kofj   = NULL; }
   if (d_fa2)    { checkCudaErrors(cudaFree(d_fa2));    d_fa2    = NULL; }
   if (d_Ar)     { checkCudaErrors(cudaFree(d_Ar));     d_Ar     = NULL; }
   if (d_Ai)     { checkCudaErrors(cudaFree(d_Ai));     d_Ai     = NULL; }
   if (d_rho)    { checkCudaErrors(cudaFree(d_rho));    d_rho    = NULL; }
   if (d_irp)    { checkCudaErrors(cudaFree(d_irp));    d_irp    = NULL; }
   if (d_wlag)   { checkCudaErrors(cudaFree(d_wlag));   d_wlag   = NULL; }

   if (di_grad_ready) {
      for (int c = 0; c < 3; ++c) {
         if (d_Gr[c]) { checkCudaErrors(cudaFree(d_Gr[c])); d_Gr[c] = NULL; }
         if (d_Gi[c]) { checkCudaErrors(cudaFree(d_Gi[c])); d_Gi[c] = NULL; }
      }
      if (d_grady) { checkCudaErrors(cudaFree(d_grady)); d_grady = NULL; }
      if (d_upos)  { checkCudaErrors(cudaFree(d_upos));  d_upos  = NULL; }
      if (d_D)     { checkCudaErrors(cudaFree(d_D));     d_D     = NULL; }
      if (d_T)     { checkCudaErrors(cudaFree(d_T));     d_T     = NULL; }
      di_grad_ready = false;
   }

   checkCublasErrors(cublasDestroy(di_cublas));
   checkCudaErrors(cudaStreamSynchronize(di_stream));
   checkCudaErrors(cudaStreamDestroy(di_stream));

   di_initialized = false;
   di_ylm_pushed  = false;
   di_mesh_pushed = false;
   di_nr_pushed   = 0;
   di_my_rank     = -1;
}
