/* *********************************************************************
 *  PseudoPotBackProj_Accel.cu
 *
 *  GPU acceleration of the "BackProjection" stage of
 *      PotentialGenerationModule::calFFTPseudoPot
 *  i.e. of the k-space summation loop inside
 *      PotentialGenerationModule::calRadialInterpolation
 *
 *  Despite the name of the enclosing routine, calFFTPseudoPot spends
 *  essentially none of its time in the Fourier transform (2 ms of
 *  69 s on the Fe3Ni full-potential reference case).  96 % of it is
 *  this back-projection: for every local atom, a sum over all local
 *  reciprocal-lattice points of a spherical Bessel function times a
 *  conjugated spherical harmonic.
 *
 *  MATHEMATICAL CONTENT (identical to the CPU path)
 *  ------------------------------------------------
 *  The CPU loop accumulates, for atom a and for every radial node
 *  r_n (n = 1..nnr) and every (l,m>=0) component jl,
 *
 *     pv(n,jl) = i2l(l) / r_n^l  *
 *                SUM_k  j_l(|k| r_n) * C_a(k,jl)
 *
 *  with
 *
 *     C_a(k,jl) = Re[kf_a(k)] * conjg(Y_{l m}(k))          l even
 *               = i * Im[kf_a(k)] * conjg(Y_{l m}(k))      l odd
 *
 *     kf_a(k)   = exp(i k.R_a) * rho~(k) * |k|^kpow
 *
 *     i2l(l)    = dummy * i^l ,
 *     dummy     = 2 (4 pi)^2   for kpow = -2   (Coulomb, 1/k^2)
 *               =    4 pi      for kpow =  0
 *
 *  The l-parity split is exactly the CPU branch
 *      even l :  pv += kfact_r * Ylm(kl) * Bj_l
 *      odd  l :  pv += i * kfact_i * Ylm(kl) * Bj_l
 *  with kfact_r = Re[kfact] and kfact_i = Im[kfact] (the CPU writes the
 *  latter as -Re[i*kfact], which is the same number).
 *
 *  Because j_l does not depend on m, and because the packed index
 *      jl(l,m) = (l+1)(l+2)/2 - l + m ,  m = 0..l
 *  is CONTIGUOUS in m for fixed l, the k-sum is a matrix-matrix
 *  product for each l separately:
 *
 *     pv(1:nnr, jl0(l) : jl0(l)+l)  =  B_l(nnr x numk) * C(numk, ...)
 *
 *     B_l(n,k) = j_l(|k| r_n)                        REAL
 *
 *  B_l is real and C is complex, so the product is issued as two real
 *  DGEMMs, one into the real accumulator and one into the imaginary
 *  accumulator.  Splitting C = Cr + i*Ci gives, per l,
 *
 *     even l :  Cr = Re[kf]*Re[Y] ,  Ci =  Re[kf]*Im[Y]
 *     odd  l :  Cr = -Im[kf]*Im[Y],  Ci =  Im[kf]*Re[Y]
 *
 *  which is what ppbp_coeffKernel writes.  This is the same
 *  real/imaginary DGEMM split that DensityInterp_Accel.cu already uses
 *  for the density interpolation.
 *
 *  The sum over k is tiled (tile width KT) so that the Bessel table
 *  B(nnr, KT, 0:lmax) stays bounded independently of the FFT grid.
 *  Tiles accumulate with DGEMM beta = 1.
 *
 *  KERNELS
 *     ppbp_kmagKernel     |k| and |k|^2 for one k-tile          (per tile)
 *     ppbp_coeffKernel    Cr, Ci  from R_a, rho~(k), Y*(k)      (per atom)
 *     ppbp_besselKernel   B(n,k,l) = j_l(|k| r_n)               (per group)
 *     ppbp_finalKernel    i2l(l)/r_n^l scaling, real -> complex (per atom)
 *
 *  The Bessel evaluation is a faithful port of
 *  MST/lib/BesselModule.F90 :: SphericalBesselReal0_v1 -- both
 *  recursion branches, the same branch-selection bit
 *  k = ishft(int(0.75*|x|)+2,-2), the same fixed-trip continued
 *  fraction of length ishft(lmax,4), the same |x| < tol early exit and
 *  the same final rescaling.  The Fortran routine itself is NOT
 *  modified.
 *
 *  Fortran interface (implicit interface, trailing underscore, every
 *  argument by reference -- same convention as DensityInterp_Accel.cu):
 *
 *    call query_pseudopot_gpu(bytes_needed, iok)
 *    call init_pseudopot_backproj_gpu(kt, kmax, nnr_max, jmax_max,
 *                                    lmax_max, na_max, my_pe)
 *    call ppbp_begin_atoms_gpu(na, lmax_a, nnr_a, jmax_a, rint_a, ld_rint)
 *    call ppbp_push_tile_gpu(kvec, fftc, ylm, nk, kmax, ldk)
 *    call ppbp_process_atom_gpu(ia, posi, nk, kpow, rebuild_bessel)
 *    call ppbp_finalize_atom_gpu(ia, dummy, pv_out)
 *    call finalize_pseudopot_backproj_gpu()
 *
 *  Book-keeping:      date                       note
 *                   09-09-26     Initial version (BackProjection offload)
 * ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <complex.h>
#include "cuComplex.h"
#include "acclib.hpp"

/*  Largest lmax for which the Bessel scratch array is kept in
 *  registers by a templated instantiation.  Above this the dynamic
 *  (local-memory) fallback is used -- correct, just slower.          */
#define PPBP_LMAX_TMPL 12

/*  Hard cap of the dynamic fallback's scratch array.                 */
#define PPBP_LMAX_CAP  32

/*  tol of MST/lib/BesselModule.F90 (tol = ten2m12).                  */
#define PPBP_BESSEL_TOL 1.0e-12

/* -------------------------------------------------------------------
 *  Persistent device state
 * ------------------------------------------------------------------- */
static bool ppbp_initialized = false;

static int  ppbp_KT       = 0;   /* k-tile width == leading dim of C/ylm */
static int  ppbp_kmax     = 0;   /* (lmax_max+1)^2                       */
static int  ppbp_nnr_max  = 0;   /* max radial nodes over atoms           */
static int  ppbp_jmax_max = 0;   /* max (l,m>=0) components               */
static int  ppbp_lmax_max = 0;
static int  ppbp_na_max   = 0;
static int  ppbp_my_rank  = -1;

/*  Per-atom description pushed once per calRadialInterpolation call.  */
static int  ppbp_na       = 0;
static int *ppbp_lmax_a   = NULL;   /* host copies, na_max              */
static int *ppbp_nnr_a    = NULL;
static int *ppbp_jmax_a   = NULL;

/*  Identity of the resident Bessel table.  ppbp_B_owner is the atom
 *  whose r_interp built it -- (nnr,lmax) alone is NOT enough to
 *  identify it, since two atoms can share both and still have
 *  different radial grids.  It is invalidated on every new k-tile,
 *  because B is a per-(tile,species) quantity.                       */
static int  ppbp_B_owner  = -1;
static int  ppbp_B_nnr    = -1;
static int  ppbp_B_lmax   = -1;

static double         *d_kvec = NULL;  /* (3 x KT)                      */
static cuDoubleComplex *d_fftc = NULL; /* (KT)                          */
static cuDoubleComplex *d_ylm  = NULL; /* (KT x kmax)  k-major          */
static double         *d_kmag = NULL;  /* (KT)                          */
static double         *d_k2   = NULL;  /* (KT)                          */
static double         *d_Cr   = NULL;  /* (KT x jmax_max)               */
static double         *d_Ci   = NULL;  /* (KT x jmax_max)               */
static double         *d_B    = NULL;  /* (nnr_max x KT x (lmax_max+1)) */
static double         *d_pvr  = NULL;  /* (nnr_max x jmax_max x na_max) */
static double         *d_pvi  = NULL;  /* (nnr_max x jmax_max x na_max) */
static double         *d_rint = NULL;  /* (nnr_max x na_max)            */
static int            *d_lofj = NULL;  /* (jmax_max)                    */
static int            *d_kofj = NULL;  /* (jmax_max)                    */
static cuDoubleComplex *d_pv  = NULL;  /* (nnr_max x jmax_max) staging  */

static cudaStream_t   ppbp_stream;
static cublasHandle_t ppbp_cublas;

static __forceinline__ int ppbp_blocks(long n, int tpb)
{
   long b = (n + tpb - 1)/tpb;
   if (b < 1)      b = 1;
   if (b > 1048576) b = 1048576;   /* grid-stride loops cover the rest */
   return (int)b;
}

/* -------------------------------------------------------------------
 *  Device spherical Bessel j_l(x), l = 0..LMAX.
 *
 *  Faithful port of SphericalBesselReal0_v1.  bj[l*stride] receives
 *  j_l(x).
 *
 *  The templated form gives the scratch array compile-time bounds and
 *  every loop a compile-time trip count, with the runtime branch point
 *  k applied as predication instead of as a loop bound.  That matters:
 *  the Fortran routine reads bscratr(k) with a RUNTIME k twice, and a
 *  dynamically indexed local array cannot live in registers -- it
 *  spills to local memory and the kernel loses roughly an order of
 *  magnitude.  Both of those reads are therefore replaced by scalars
 *  captured as the recursions run:
 *
 *    a_fwd  == bscratr(k) as left by the FORWARD recursion  (the last
 *              value it wrote, or sin(x) when k = 0)
 *    bk_bwd == bscratr(k) as left by the BACKWARD recursion (the value
 *              written by its last executed iteration, l = k+1, or the
 *              pre-loop bscr[LMAX-1] when that loop is empty)
 *
 *  so that every bscr[] subscript is a compile-time constant.
 * ------------------------------------------------------------------- */
template<int LMAX>
__device__ __forceinline__
void ppbp_bessel_t(double x, double * __restrict__ bj, long stride)
{
   double bscr[LMAX+1];

   if (fabs(x) < PPBP_BESSEL_TOL) {
      bj[0] = 1.0;
#pragma unroll
      for (int l = 1; l <= LMAX; ++l) bj[l*stride] = 0.0;
      return;
   }

   double as, ac;
   sincos(x, &as, &ac);
   const double x1 = 1.0/x;

   bj[0] = as*x1;

   /* ----------------------------------------------------------------
    *  Forward recursion for small l:  bscr(l) = (2l-1) bj(l-1) - bscr(l-2)
    * ---------------------------------------------------------------- */
   int k = ((int)(0.75*fabs(x)) + 2) >> 2;
   bscr[0] = as;
   double a_fwd = as;                 /* == bscr[k] for k = 0 */

   bool done = false;
   if (k >= 1) {
      if (k > LMAX) k = LMAX;
      bscr[1]    = bj[0] - ac;
      bj[stride] = bscr[1]*x1;
      a_fwd      = bscr[1];
#pragma unroll
      for (int l = 2; l <= LMAX; ++l) {
         if (l <= k) {
            bscr[l]      = (double)(2*l-1)*bj[(l-1)*stride] - bscr[l-2];
            bj[l*stride] = bscr[l]*x1;
            a_fwd        = bscr[l];
         }
      }
      if (k == LMAX) done = true;
   }

   if (!done) {
      /* -------------------------------------------------------------
       *  Backward recursion from very large l down to l = k, seeded by
       *  the continued fraction  aj = (2l+1) - x^2/aj  run from
       *  nm = ishft(lmax,4) down to lmax+2.
       * ------------------------------------------------------------- */
      const double a  = a_fwd;
      const int    nm = LMAX << 4;
      double aj = (double)(2*nm + 3);
      const double x2 = x*x;
      for (int l = nm; l >= LMAX+2; --l) {
         aj = (double)(2*l + 1) - x2/aj;
      }
      bscr[LMAX]      = (double)(2*LMAX+3)*aj*x1 - x;
      bj[LMAX*stride] = bscr[LMAX]*x1;
      bscr[LMAX-1]    = (double)(2*LMAX+1)*bj[LMAX*stride] - aj;
      double bk_bwd   = bscr[LMAX-1];     /* k == LMAX-1: loop is empty */
#pragma unroll
      for (int l = LMAX-1; l >= 1; --l) {
         if (l >= k+1) {
            bj[l*stride] = bscr[l]*x1;
            bscr[l-1]    = (double)(2*l+1)*bj[l*stride] - bscr[l+1];
            bk_bwd       = bscr[l-1];
         }
      }
      /* rescale so that bj(k) is the forward-recursion value */
      aj = a/bk_bwd;
#pragma unroll
      for (int l = 1; l <= LMAX; ++l) {
         if (l >= k+1) bj[l*stride] *= aj;
      }
   }
}

/*  Dynamic fallback: same algorithm, runtime lmax, scratch array in
 *  local memory.  Used for lmax = 0, 1 and lmax > PPBP_LMAX_TMPL.    */
__device__
void ppbp_bessel_dyn(int lmax, double x, double * __restrict__ bj, long stride)
{
   double bscr[PPBP_LMAX_CAP+1];

   if (fabs(x) < PPBP_BESSEL_TOL) {
      bj[0] = 1.0;
      for (int l = 1; l <= lmax; ++l) bj[l*stride] = 0.0;
      return;
   }

   double as, ac;
   sincos(x, &as, &ac);
   const double x1 = 1.0/x;

   bj[0] = as*x1;
   if (lmax == 0) return;
   if (lmax == 1) { bj[stride] = (bj[0] - ac)*x1; return; }

   int k = ((int)(0.75*fabs(x)) + 2) >> 2;
   bscr[0] = as;
   if (k >= 1) {
      if (k > lmax) k = lmax;
      bscr[1]    = bj[0] - ac;
      bj[stride] = bscr[1]*x1;
      for (int l = 2; l <= k; ++l) {
         bscr[l]      = (double)(2*l-1)*bj[(l-1)*stride] - bscr[l-2];
         bj[l*stride] = bscr[l]*x1;
      }
      if (k == lmax) return;
   }

   const double a  = bscr[k];
   const int    nm = lmax << 4;
   double aj = (double)(2*nm + 3);
   const double x2 = x*x;
   for (int l = nm; l >= lmax+2; --l) {
      aj = (double)(2*l + 1) - x2/aj;
   }
   bscr[lmax]      = (double)(2*lmax+3)*aj*x1 - x;
   bj[lmax*stride] = bscr[lmax]*x1;
   bscr[lmax-1]    = (double)(2*lmax+1)*bj[lmax*stride] - aj;
   for (int l = lmax-1; l >= k+1; --l) {
      bj[l*stride] = bscr[l]*x1;
      bscr[l-1]    = (double)(2*l+1)*bj[l*stride] - bscr[l+1];
   }
   aj = a/bscr[k];
   for (int l = k+1; l <= lmax; ++l) bj[l*stride] *= aj;
}

/* -------------------------------------------------------------------
 *  Kernel: |k| and |k|^2 for one k-tile
 * ------------------------------------------------------------------- */
__global__
void ppbp_kmagKernel(const double * __restrict__ kvec, int nk,
                     double * __restrict__ kmag, double * __restrict__ k2)
{
   for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < nk;
            i += blockDim.x*gridDim.x) {
      const double kx = kvec[3*(long)i    ];
      const double ky = kvec[3*(long)i + 1];
      const double kz = kvec[3*(long)i + 2];
      const double q  = kx*kx + ky*ky + kz*kz;
      k2[i]   = q;
      kmag[i] = sqrt(q);
   }
}

/* -------------------------------------------------------------------
 *  Kernel: build the complex coefficient matrix C(k,jl) for one atom
 *
 *     kf   = exp(i k.R) * rho~(k) * w ,   w = 1/k^2 (kpow=-2) or 1
 *     even l :  Cr = Re[kf]*Re[Y*] ,  Ci = Re[kf]*Im[Y*]
 *     odd  l :  Cr = -Im[kf]*Im[Y*],  Ci = Im[kf]*Re[Y*]
 *
 *  One thread per k-point, looping over jl, so exp(i k.R) is evaluated
 *  once per k instead of once per (k,jl) -- and the stores are
 *  coalesced because k is the fast index of C.
 * ------------------------------------------------------------------- */
__global__
void ppbp_coeffKernel(int nk, long ldk, int jmax,
                      double px, double py, double pz, int kpow,
                      const double * __restrict__ kvec,
                      const cuDoubleComplex * __restrict__ fftc,
                      const cuDoubleComplex * __restrict__ ylm,
                      const double * __restrict__ k2,
                      const int * __restrict__ lofj,
                      const int * __restrict__ kofj,
                      double * __restrict__ Cr, double * __restrict__ Ci)
{
   for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < nk;
            i += blockDim.x*gridDim.x) {

      const double kdr = kvec[3*(long)i    ]*px
                       + kvec[3*(long)i + 1]*py
                       + kvec[3*(long)i + 2]*pz;
      double s, c;
      sincos(kdr, &s, &c);

      const cuDoubleComplex f = fftc[i];
      /* (c + i s)*(fr + i fi) */
      double ar = c*cuCreal(f) - s*cuCimag(f);
      double ai = c*cuCimag(f) + s*cuCreal(f);

      if (kpow == -2) {
         const double w = 1.0/k2[i];
         ar *= w;
         ai *= w;
      }
      /* kpow == 0 : no weight.  Other kpow are rejected on the host. */

      for (int jl = 0; jl < jmax; ++jl) {
         const int l  = lofj[jl];
         const int kl = kofj[jl] - 1;                 /* to 0-based */
         const cuDoubleComplex y = ylm[(long)i + ldk*(long)kl];
         double cr, ci;
         if ((l & 1) == 0) {
            cr = ar*cuCreal(y);
            ci = ar*cuCimag(y);
         } else {
            cr = -ai*cuCimag(y);
            ci =  ai*cuCreal(y);
         }
         Cr[(long)i + ldk*(long)jl] = cr;
         Ci[(long)i + ldk*(long)jl] = ci;
      }
   }
}

/* -------------------------------------------------------------------
 *  Kernel: Bessel table  B(n, k, l) = j_l(|k| r_n)
 *
 *  Layout  B[n + nnr*kt + nnr*ldk*l]  -- n fastest, so the stores are
 *  coalesced AND each B_l is a contiguous column-major (nnr x ldk)
 *  block with lda = nnr, exactly what cublasDgemm wants.  No transpose
 *  and no packing anywhere.
 * ------------------------------------------------------------------- */
template<int LMAX>
__global__
void ppbp_besselKernel(int nnr, int nk, long ldk,
                       const double * __restrict__ kmag,
                       const double * __restrict__ rint,
                       double * __restrict__ B)
{
   const long ntot = (long)nnr*(long)nk;
   for (long idx = (long)blockIdx.x*blockDim.x + threadIdx.x; idx < ntot;
             idx += (long)blockDim.x*gridDim.x) {
      const int n  = (int)(idx % (long)nnr);
      const int kt = (int)(idx / (long)nnr);
      ppbp_bessel_t<LMAX>( kmag[kt]*rint[n],
                           B + (long)n + (long)nnr*(long)kt,
                           (long)nnr*ldk );
   }
}

__global__
void ppbp_besselKernelDyn(int lmax, int nnr, int nk, long ldk,
                          const double * __restrict__ kmag,
                          const double * __restrict__ rint,
                          double * __restrict__ B)
{
   const long ntot = (long)nnr*(long)nk;
   for (long idx = (long)blockIdx.x*blockDim.x + threadIdx.x; idx < ntot;
             idx += (long)blockDim.x*gridDim.x) {
      const int n  = (int)(idx % (long)nnr);
      const int kt = (int)(idx / (long)nnr);
      ppbp_bessel_dyn( lmax, kmag[kt]*rint[n],
                       B + (long)n + (long)nnr*(long)kt,
                       (long)nnr*ldk );
   }
}

/* -------------------------------------------------------------------
 *  Kernel: final scaling and real -> complex packing
 *
 *     pv(n,jl) = i2l(l) * ( pvr(n,jl) + i pvi(n,jl) ) / r_n^l
 *     i2l(l)   = dummy * i^l
 *
 *  The multiply-by-i2l-then-divide-by-r^l order, and the repeated
 *  multiplication used to form r^l, reproduce the CPU expression
 *      pv_interp(i,jl) = i2l*pv_interp(i,jl)/(r_interp(i)**l)
 *  operation for operation.
 * ------------------------------------------------------------------- */
__global__
void ppbp_finalKernel(int nnr, int jmax, double dummy,
                      const double * __restrict__ rint,
                      const int * __restrict__ lofj,
                      const double * __restrict__ pvr,
                      const double * __restrict__ pvi,
                      cuDoubleComplex * __restrict__ pv)
{
   const long ntot = (long)nnr*(long)jmax;
   for (long idx = (long)blockIdx.x*blockDim.x + threadIdx.x; idx < ntot;
             idx += (long)blockDim.x*gridDim.x) {
      const int n  = (int)(idx % (long)nnr);
      const int jl = (int)(idx / (long)nnr);
      const int l  = lofj[jl];

      const double a = pvr[idx];
      const double b = pvi[idx];

      /* i2l * (a + i b) with i2l = dummy * i^l */
      double tr, ti;
      switch (l & 3) {
         case 0:  tr =  dummy*a;  ti =  dummy*b;  break;
         case 1:  tr = -dummy*b;  ti =  dummy*a;  break;
         case 2:  tr = -dummy*a;  ti = -dummy*b;  break;
         default: tr =  dummy*b;  ti = -dummy*a;  break;
      }

      double rl = 1.0;
      const double r = rint[n];
      for (int j = 0; j < l; ++j) rl *= r;

      pv[idx] = make_cuDoubleComplex(tr/rl, ti/rl);
   }
}

/* ===================================================================
 *  query_pseudopot_gpu(bytes_needed, ok)
 *
 *  Runtime gate.  Returns ok = 1 only if a device exists, can be
 *  selected, has at least bytes_needed free, and a cuBLAS handle can be
 *  created.  Any failure clears the CUDA error state and returns 0, so
 *  a failed probe cannot poison later CUDA use.  MUST_PSEUDOPOT_GPU=0
 *  forces ok = 0 without touching the driver at all.
 * =================================================================== */
extern "C"
void query_pseudopot_gpu_(long *bytes_needed, int *ok)
{
   *ok = 0;

   const char *env = getenv("MUST_PSEUDOPOT_GPU");
   if (env != NULL && (env[0] == '0' || env[0] == 'n' || env[0] == 'N' ||
                       env[0] == 'f' || env[0] == 'F')) {
      return;
   }

   int ndev = 0;
   if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev < 1) {
      cudaGetLastError();
      return;
   }

   int dev = 0;
   if (cudaGetDevice(&dev) != cudaSuccess) {
      cudaGetLastError();
      return;
   }
   if (cudaSetDevice(dev) != cudaSuccess) {
      cudaGetLastError();
      return;
   }

   size_t freeb = 0, totb = 0;
   if (cudaMemGetInfo(&freeb, &totb) != cudaSuccess) {
      cudaGetLastError();
      return;
   }
   if (*bytes_needed > 0 && (size_t)(*bytes_needed) > freeb) {
      return;
   }

   cublasHandle_t h;
   if (cublasCreate(&h) != CUBLAS_STATUS_SUCCESS) {
      cudaGetLastError();
      return;
   }
   cublasDestroy(h);
   cudaGetLastError();

   *ok = 1;
}

/* ===================================================================
 *  init_pseudopot_backproj_gpu(kt, kmax, nnr_max, jmax_max, lmax_max,
 *                              na_max, my_pe)
 * =================================================================== */
extern "C"
void init_pseudopot_backproj_gpu_(int *kt, int *kmax, int *nnr_max,
                                  int *jmax_max, int *lmax_max,
                                  int *na_max, int *my_pe)
{
   if (ppbp_initialized) return;

   ppbp_KT       = *kt;
   ppbp_kmax     = *kmax;
   ppbp_nnr_max  = *nnr_max;
   ppbp_jmax_max = *jmax_max;
   ppbp_lmax_max = *lmax_max;
   ppbp_na_max   = *na_max;
   ppbp_my_rank  = *my_pe;

   if (ppbp_KT < 1 || ppbp_kmax < 1 || ppbp_nnr_max < 1 ||
       ppbp_jmax_max < 1 || ppbp_lmax_max < 0 || ppbp_na_max < 1) {
      fprintf(stderr, "\nError in init_pseudopot_backproj_gpu: bad "
                      "dimensions (KT %d, kmax %d, nnr %d, jmax %d, "
                      "lmax %d, na %d)\n", ppbp_KT, ppbp_kmax,
                      ppbp_nnr_max, ppbp_jmax_max, ppbp_lmax_max,
                      ppbp_na_max);
      exit(EXIT_FAILURE);
   }
   if (ppbp_lmax_max > PPBP_LMAX_CAP) {
      fprintf(stderr, "\nError in init_pseudopot_backproj_gpu: lmax %d "
                      "exceeds PPBP_LMAX_CAP %d\n",
                      ppbp_lmax_max, PPBP_LMAX_CAP);
      exit(EXIT_FAILURE);
   }

   const long KT   = ppbp_KT;
   const long nnr  = ppbp_nnr_max;
   const long jmx  = ppbp_jmax_max;
   const long nl   = ppbp_lmax_max + 1;
   const long na   = ppbp_na_max;

   checkCudaErrors(cudaStreamCreate(&ppbp_stream));
   checkCublasErrors(cublasCreate(&ppbp_cublas));
   checkCublasErrors(cublasSetStream(ppbp_cublas, ppbp_stream));

   checkCudaErrors(cudaMalloc((void**)&d_kvec, 3*KT*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_fftc, KT*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_ylm,
                              KT*(long)ppbp_kmax*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_kmag, KT*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_k2,   KT*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_Cr,   KT*jmx*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_Ci,   KT*jmx*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_B,    nnr*KT*nl*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_pvr,  nnr*jmx*na*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_pvi,  nnr*jmx*na*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_rint, nnr*na*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_lofj, jmx*sizeof(int)));
   checkCudaErrors(cudaMalloc((void**)&d_kofj, jmx*sizeof(int)));
   checkCudaErrors(cudaMalloc((void**)&d_pv,
                              nnr*jmx*sizeof(cuDoubleComplex)));

   ppbp_lmax_a = (int*)malloc(na*sizeof(int));
   ppbp_nnr_a  = (int*)malloc(na*sizeof(int));
   ppbp_jmax_a = (int*)malloc(na*sizeof(int));
   if (!ppbp_lmax_a || !ppbp_nnr_a || !ppbp_jmax_a) {
      fprintf(stderr, "\nError in init_pseudopot_backproj_gpu: host "
                      "allocation failed\n");
      exit(EXIT_FAILURE);
   }

   /*  lofj / kofj for the packed (l, m>=0) index
    *     jl(l,m) = (l+1)(l+2)/2 - l + m ,  kl(l,m) = (l+1)^2 - l + m   */
   int *h_lofj = (int*)malloc(jmx*sizeof(int));
   int *h_kofj = (int*)malloc(jmx*sizeof(int));
   for (long j = 0; j < jmx; ++j) { h_lofj[j] = 0; h_kofj[j] = 1; }
   for (int l = 0; l <= ppbp_lmax_max; ++l) {
      for (int m = 0; m <= l; ++m) {
         const int jl = ((l+1)*(l+2))/2 - l + m;      /* 1-based */
         const int kl =  (l+1)*(l+1)    - l + m;      /* 1-based */
         if (jl >= 1 && jl <= (int)jmx) {
            h_lofj[jl-1] = l;
            h_kofj[jl-1] = kl;
         }
      }
   }
   checkCudaErrors(cudaMemcpy(d_lofj, h_lofj, jmx*sizeof(int),
                              cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(d_kofj, h_kofj, jmx*sizeof(int),
                              cudaMemcpyHostToDevice));
   free(h_lofj);
   free(h_kofj);

   ppbp_B_nnr  = -1;
   ppbp_B_lmax = -1;
   ppbp_initialized = true;
}

/* ===================================================================
 *  ppbp_begin_atoms_gpu(na, lmax_a, nnr_a, jmax_a, rint_a, ld_rint)
 *
 *  Starts one calRadialInterpolation call: records the per-atom
 *  geometry, uploads the radial interpolation nodes and zeroes the
 *  accumulators.
 * =================================================================== */
extern "C"
void ppbp_begin_atoms_gpu_(int *na, int *lmax_a, int *nnr_a, int *jmax_a,
                           double *rint_a, int *ld_rint)
{
   if (!ppbp_initialized) {
      fprintf(stderr, "\nError in ppbp_begin_atoms_gpu: "
                      "init_pseudopot_backproj_gpu must be called first\n");
      exit(EXIT_FAILURE);
   }
   if (*na < 1 || *na > ppbp_na_max) {
      fprintf(stderr, "\nError in ppbp_begin_atoms_gpu: na %d outside "
                      "[1, %d]\n", *na, ppbp_na_max);
      exit(EXIT_FAILURE);
   }

   ppbp_na = *na;
   for (int i = 0; i < ppbp_na; ++i) {
      ppbp_lmax_a[i] = lmax_a[i];
      ppbp_nnr_a[i]  = nnr_a[i];
      ppbp_jmax_a[i] = jmax_a[i];
      if (nnr_a[i] < 1 || nnr_a[i] > ppbp_nnr_max ||
          jmax_a[i] < 1 || jmax_a[i] > ppbp_jmax_max ||
          lmax_a[i] < 0 || lmax_a[i] > ppbp_lmax_max) {
         fprintf(stderr, "\nError in ppbp_begin_atoms_gpu: atom %d has "
                         "nnr %d/%d, jmax %d/%d, lmax %d/%d\n", i+1,
                         nnr_a[i], ppbp_nnr_max, jmax_a[i],
                         ppbp_jmax_max, lmax_a[i], ppbp_lmax_max);
         exit(EXIT_FAILURE);
      }
   }

   /*  rint_a is the Fortran array rint_a(ld_rint, na); only the leading
    *  nnr rows of each column are meaningful.                        */
   checkCudaErrors(cudaMemcpy2DAsync(d_rint, ppbp_nnr_max*sizeof(double),
                   rint_a, (*ld_rint)*sizeof(double),
                   ppbp_nnr_max*sizeof(double), (size_t)ppbp_na,
                   cudaMemcpyHostToDevice, ppbp_stream));

   const long nz = (long)ppbp_nnr_max*(long)ppbp_jmax_max*(long)ppbp_na;
   checkCudaErrors(cudaMemsetAsync(d_pvr, 0, nz*sizeof(double), ppbp_stream));
   checkCudaErrors(cudaMemsetAsync(d_pvi, 0, nz*sizeof(double), ppbp_stream));
   checkCudaErrors(cudaStreamSynchronize(ppbp_stream));

   ppbp_B_owner = -1;
   ppbp_B_nnr   = -1;
   ppbp_B_lmax  = -1;
}

/* ===================================================================
 *  ppbp_push_tile_gpu(kvec, fftc, ylm, nk, kmax, ldk)
 *
 *  Uploads one k-tile.  ylm is the Fortran array ylm(ldk, kmax) holding
 *  conjg(Y_lm(k_hat)) with the k index leading, so that the coefficient
 *  kernel's loads are coalesced.  Only the first nk rows are copied.
 * =================================================================== */
extern "C"
void ppbp_push_tile_gpu_(double *kvec, double _Complex *fftc,
                         double _Complex *ylm, int *nk, int *kmax, int *ldk)
{
   if (!ppbp_initialized) {
      fprintf(stderr, "\nError in ppbp_push_tile_gpu: not initialised\n");
      exit(EXIT_FAILURE);
   }
   if (*nk < 1 || *nk > ppbp_KT || *ldk != ppbp_KT || *kmax > ppbp_kmax) {
      fprintf(stderr, "\nError in ppbp_push_tile_gpu: nk %d, ldk %d/%d, "
                      "kmax %d/%d\n", *nk, *ldk, ppbp_KT, *kmax, ppbp_kmax);
      exit(EXIT_FAILURE);
   }

   const int nk_l = *nk;

   checkCudaErrors(cudaMemcpyAsync(d_kvec, kvec,
                   3*(size_t)nk_l*sizeof(double),
                   cudaMemcpyHostToDevice, ppbp_stream));
   checkCudaErrors(cudaMemcpyAsync(d_fftc, fftc,
                   (size_t)nk_l*sizeof(cuDoubleComplex),
                   cudaMemcpyHostToDevice, ppbp_stream));
   checkCudaErrors(cudaMemcpy2DAsync(d_ylm,
                   (size_t)ppbp_KT*sizeof(cuDoubleComplex),
                   ylm, (size_t)(*ldk)*sizeof(cuDoubleComplex),
                   (size_t)nk_l*sizeof(cuDoubleComplex), (size_t)(*kmax),
                   cudaMemcpyHostToDevice, ppbp_stream));

   const int tpb = 256;
   ppbp_kmagKernel<<<ppbp_blocks(nk_l,tpb), tpb, 0, ppbp_stream>>>
                  (d_kvec, nk_l, d_kmag, d_k2);
   checkCudaErrors(cudaPeekAtLastError());
   checkCudaErrors(cudaStreamSynchronize(ppbp_stream));

   /*  B depends on |k|, so the table from the previous tile is stale.
    *  Invalidating here means a caller that forgets to rebuild it on a
    *  new tile gets a hard error instead of a silently wrong answer.  */
   ppbp_B_owner = -1;
   ppbp_B_nnr   = -1;
   ppbp_B_lmax  = -1;
}

/* -------------------------------------------------------------------
 *  Bessel-table launcher: dispatches on lmax so that the scratch array
 *  of the templated device routine has compile-time bounds.
 *
 *  ia_owner is the atom whose radial grid the table is built from.  It
 *  MUST be applied to d_rint, which is an (nnr_max x na) array whose
 *  column ia holds that atom's r_interp: without the offset every
 *  table is built from atom 1's grid, while ppbp_finalize_atom_gpu
 *  divides by the atom's OWN r_n^l.  The first species group happens
 *  to be right (its representative is atom 1) and every later group is
 *  silently wrong, the more so the higher l, since the 1/r_n^l
 *  division amplifies a grid mismatch at the innermost nodes by up to
 *  r^-lmax.
 * ------------------------------------------------------------------- */
static void ppbp_launch_bessel(int ia_owner, int lmax, int nnr,
                               int nk, long ldk)
{
   if (ia_owner < 0 || ia_owner >= ppbp_na) {
      fprintf(stderr, "\nError in ppbp_launch_bessel: ia_owner %d outside "
                      "[0, %d)\n", ia_owner, ppbp_na);
      exit(EXIT_FAILURE);
   }
   if (ppbp_nnr_a[ia_owner] != nnr || ppbp_lmax_a[ia_owner] != lmax) {
      fprintf(stderr, "\nError in ppbp_launch_bessel: table for atom %d "
                      "requested with nnr %d, lmax %d, but that atom has "
                      "nnr %d, lmax %d\n", ia_owner + 1, nnr, lmax,
                      ppbp_nnr_a[ia_owner], ppbp_lmax_a[ia_owner]);
      exit(EXIT_FAILURE);
   }

   const double * const rint_a = d_rint
                               + (long)ppbp_nnr_max*(long)ia_owner;

   const int  tpb  = 128;
   const long ntot = (long)nnr*(long)nk;
   const int  nb   = ppbp_blocks(ntot, tpb);

#define PPBP_CASE(L)                                                   \
   case L: ppbp_besselKernel<L><<<nb,tpb,0,ppbp_stream>>>              \
                            (nnr, nk, ldk, d_kmag, rint_a, d_B); break;

   switch (lmax) {
      PPBP_CASE(2)  PPBP_CASE(3)  PPBP_CASE(4)  PPBP_CASE(5)
      PPBP_CASE(6)  PPBP_CASE(7)  PPBP_CASE(8)  PPBP_CASE(9)
      PPBP_CASE(10) PPBP_CASE(11) PPBP_CASE(12)
      default:
         ppbp_besselKernelDyn<<<nb,tpb,0,ppbp_stream>>>
                             (lmax, nnr, nk, ldk, d_kmag, rint_a, d_B);
         break;
   }
#undef PPBP_CASE

   checkCudaErrors(cudaPeekAtLastError());
}

/* ===================================================================
 *  ppbp_process_atom_gpu(ia, posi, nk, kpow, ia_bessel)
 *
 *  Adds the contribution of the resident k-tile to atom ia.
 *
 *  ia_bessel names the atom whose (r_interp, nnr, lmax) the Bessel
 *  table is to be built from.  ia_bessel == ia rebuilds it; otherwise
 *  the resident table is reused and MUST already belong to ia_bessel.
 *  Atoms sharing a radial grid -- which is per SPECIES, not per atom,
 *  since r_interp is generated from the radial mesh ends, the inscribed
 *  sphere radius and d_ir -- therefore share the table, so the
 *  expensive Bessel kernel runs once per (tile, species) rather than
 *  once per (tile, atom).
 *
 *  The reuse check compares the OWNING ATOM, not just (nnr, lmax):
 *  two atoms can agree on both and still have different radial grids,
 *  in which case a (nnr, lmax) test would pass and the result would be
 *  silently wrong.
 * =================================================================== */
extern "C"
void ppbp_process_atom_gpu_(int *ia, double *posi, int *nk, int *kpow,
                            int *ia_bessel)
{
   if (!ppbp_initialized || ppbp_na < 1) {
      fprintf(stderr, "\nError in ppbp_process_atom_gpu: "
                      "ppbp_begin_atoms_gpu must be called first\n");
      exit(EXIT_FAILURE);
   }
   const int a = *ia - 1;
   if (a < 0 || a >= ppbp_na) {
      fprintf(stderr, "\nError in ppbp_process_atom_gpu: ia %d outside "
                      "[1, %d]\n", *ia, ppbp_na);
      exit(EXIT_FAILURE);
   }
   if (*kpow != -2 && *kpow != 0) {
      fprintf(stderr, "\nError in ppbp_process_atom_gpu: unsupported "
                      "kpow %d\n", *kpow);
      exit(EXIT_FAILURE);
   }

   const int  nk_l = *nk;
   const long ldk  = ppbp_KT;
   const int  nnr  = ppbp_nnr_a[a];
   const int  lmax = ppbp_lmax_a[a];
   const int  jmax = ppbp_jmax_a[a];

   /* ---- C(k,jl) for this atom ---------------------------------- */
   {
      const int tpb = 256;
      ppbp_coeffKernel<<<ppbp_blocks(nk_l,tpb), tpb, 0, ppbp_stream>>>
                      (nk_l, ldk, jmax, posi[0], posi[1], posi[2], *kpow,
                       d_kvec, d_fftc, d_ylm, d_k2, d_lofj, d_kofj,
                       d_Cr, d_Ci);
      checkCudaErrors(cudaPeekAtLastError());
   }

   /* ---- B(n,k,l) for this atom's radial grid -------------------- */
   const int ab = *ia_bessel - 1;
   if (ab < 0 || ab >= ppbp_na) {
      fprintf(stderr, "\nError in ppbp_process_atom_gpu: ia_bessel %d "
                      "outside [1, %d]\n", *ia_bessel, ppbp_na);
      exit(EXIT_FAILURE);
   }
   if (ab == a) {
      ppbp_launch_bessel(ab, lmax, nnr, nk_l, ldk);
      ppbp_B_owner = a;
      ppbp_B_nnr   = nnr;
      ppbp_B_lmax  = lmax;
   } else if (ppbp_B_owner != ab || ppbp_B_nnr != nnr ||
              ppbp_B_lmax != lmax) {
      fprintf(stderr, "\nError in ppbp_process_atom_gpu: atom %d asked to "
                      "reuse the Bessel table of atom %d, but the resident "
                      "table belongs to atom %d (nnr %d/%d, lmax %d/%d). "
                      "The table is invalidated on every k-tile, so the "
                      "owner must be processed first within each tile.\n",
                      *ia, *ia_bessel, ppbp_B_owner + 1,
                      ppbp_B_nnr, nnr, ppbp_B_lmax, lmax);
      exit(EXIT_FAILURE);
   }

   /* ---- accumulate:  pv_l += B_l * C_l ------------------------- *
    *  For fixed l the packed index jl is contiguous in m, so the
    *  columns jl0(l) .. jl0(l)+l of C are already the right block --
    *  no gather, no packing.                                        */
   {
      const double one = 1.0;
      const long   off = (long)ppbp_nnr_max*(long)ppbp_jmax_max*(long)a;
      for (int l = 0; l <= lmax; ++l) {
         const int  n    = l + 1;
         const long col0 = ((long)(l+1)*(l+2))/2 - l - 1;   /* 0-based */
         if (col0 + n > jmax) break;
         checkCublasErrors(cublasDgemm(ppbp_cublas,
                    CUBLAS_OP_N, CUBLAS_OP_N,
                    nnr, n, nk_l,
                    &one, d_B  + (long)nnr*ldk*(long)l, nnr,
                          d_Cr + ldk*col0,              (int)ldk,
                    &one, d_pvr + off + (long)nnr*col0, nnr));
         checkCublasErrors(cublasDgemm(ppbp_cublas,
                    CUBLAS_OP_N, CUBLAS_OP_N,
                    nnr, n, nk_l,
                    &one, d_B  + (long)nnr*ldk*(long)l, nnr,
                          d_Ci + ldk*col0,              (int)ldk,
                    &one, d_pvi + off + (long)nnr*col0, nnr));
      }
   }

   checkCudaErrors(cudaStreamSynchronize(ppbp_stream));
}

/* ===================================================================
 *  ppbp_finalize_atom_gpu(ia, dummy, pv_out)
 *
 *  Applies i2l(l)/r^l and copies the (nnr x jmax) complex block back
 *  to the host, packed with leading dimension nnr -- which is exactly
 *  the layout that calRadialProjection expects from w_interp(:,ia).
 * =================================================================== */
extern "C"
void ppbp_finalize_atom_gpu_(int *ia, double *dummy, double _Complex *pv_out)
{
   if (!ppbp_initialized || ppbp_na < 1) {
      fprintf(stderr, "\nError in ppbp_finalize_atom_gpu: not started\n");
      exit(EXIT_FAILURE);
   }
   const int a = *ia - 1;
   if (a < 0 || a >= ppbp_na) {
      fprintf(stderr, "\nError in ppbp_finalize_atom_gpu: ia %d outside "
                      "[1, %d]\n", *ia, ppbp_na);
      exit(EXIT_FAILURE);
   }

   const int  nnr  = ppbp_nnr_a[a];
   const int  jmax = ppbp_jmax_a[a];
   const long off  = (long)ppbp_nnr_max*(long)ppbp_jmax_max*(long)a;

   const int  tpb  = 256;
   const long ntot = (long)nnr*(long)jmax;
   ppbp_finalKernel<<<ppbp_blocks(ntot,tpb), tpb, 0, ppbp_stream>>>
                   (nnr, jmax, *dummy,
                    d_rint + (long)ppbp_nnr_max*(long)a, d_lofj,
                    d_pvr + off, d_pvi + off, d_pv);
   checkCudaErrors(cudaPeekAtLastError());

   checkCudaErrors(cudaMemcpyAsync(pv_out, d_pv,
                   (size_t)ntot*sizeof(cuDoubleComplex),
                   cudaMemcpyDeviceToHost, ppbp_stream));
   checkCudaErrors(cudaStreamSynchronize(ppbp_stream));
}

/* ===================================================================
 *  finalize_pseudopot_backproj_gpu()
 * =================================================================== */
extern "C"
void finalize_pseudopot_backproj_gpu_(void)
{
   if (!ppbp_initialized) return;

   cudaFree(d_kvec);  d_kvec = NULL;
   cudaFree(d_fftc);  d_fftc = NULL;
   cudaFree(d_ylm);   d_ylm  = NULL;
   cudaFree(d_kmag);  d_kmag = NULL;
   cudaFree(d_k2);    d_k2   = NULL;
   cudaFree(d_Cr);    d_Cr   = NULL;
   cudaFree(d_Ci);    d_Ci   = NULL;
   cudaFree(d_B);     d_B    = NULL;
   cudaFree(d_pvr);   d_pvr  = NULL;
   cudaFree(d_pvi);   d_pvi  = NULL;
   cudaFree(d_rint);  d_rint = NULL;
   cudaFree(d_lofj);  d_lofj = NULL;
   cudaFree(d_kofj);  d_kofj = NULL;
   cudaFree(d_pv);    d_pv   = NULL;

   free(ppbp_lmax_a); ppbp_lmax_a = NULL;
   free(ppbp_nnr_a);  ppbp_nnr_a  = NULL;
   free(ppbp_jmax_a); ppbp_jmax_a = NULL;

   cublasDestroy(ppbp_cublas);
   cudaStreamDestroy(ppbp_stream);

   ppbp_na = 0;
   ppbp_initialized = false;
}

/* ===================================================================
 *  ppbp_bessel_ref_gpu(lmax, nx, x, bj)
 *
 *  Validation entry point (test V1): evaluates the DEVICE Bessel
 *  routine on an arbitrary argument list so that it can be compared
 *  against BesselModule::SphericalBessel on the host.  bj is the
 *  Fortran array bj(nx, 0:lmax).
 * =================================================================== */
__global__
void ppbp_besselRefKernel(int lmax, int nx, const double * __restrict__ x,
                          double * __restrict__ bj)
{
   for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < nx;
            i += blockDim.x*gridDim.x) {
      ppbp_bessel_dyn(lmax, x[i], bj + i, (long)nx);
   }
}

template<int LMAX>
__global__
void ppbp_besselRefKernelT(int nx, const double * __restrict__ x,
                           double * __restrict__ bj)
{
   for (int i = blockIdx.x*blockDim.x + threadIdx.x; i < nx;
            i += blockDim.x*gridDim.x) {
      ppbp_bessel_t<LMAX>(x[i], bj + i, (long)nx);
   }
}

extern "C"
void ppbp_bessel_ref_gpu_(int *lmax, int *nx, double *x, double *bj,
                          int *use_template)
{
   const int n = *nx;
   const int L = *lmax;
   double *dx = NULL, *dbj = NULL;
   checkCudaErrors(cudaMalloc((void**)&dx,  (size_t)n*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&dbj,
                              (size_t)n*(size_t)(L+1)*sizeof(double)));
   checkCudaErrors(cudaMemcpy(dx, x, (size_t)n*sizeof(double),
                              cudaMemcpyHostToDevice));

   const int tpb = 128;
   const int nb  = ppbp_blocks(n, tpb);

#define PPBP_REFCASE(LL)                                               \
   case LL: ppbp_besselRefKernelT<LL><<<nb,tpb>>>(n, dx, dbj); break;

   if (*use_template != 0) {
      switch (L) {
         PPBP_REFCASE(2)  PPBP_REFCASE(3)  PPBP_REFCASE(4)
         PPBP_REFCASE(5)  PPBP_REFCASE(6)  PPBP_REFCASE(7)
         PPBP_REFCASE(8)  PPBP_REFCASE(9)  PPBP_REFCASE(10)
         PPBP_REFCASE(11) PPBP_REFCASE(12)
         default: ppbp_besselRefKernel<<<nb,tpb>>>(L, n, dx, dbj); break;
      }
   } else {
      ppbp_besselRefKernel<<<nb,tpb>>>(L, n, dx, dbj);
   }
#undef PPBP_REFCASE

   checkCudaErrors(cudaPeekAtLastError());
   checkCudaErrors(cudaDeviceSynchronize());
   checkCudaErrors(cudaMemcpy(bj, dbj,
                   (size_t)n*(size_t)(L+1)*sizeof(double),
                   cudaMemcpyDeviceToHost));
   cudaFree(dx);
   cudaFree(dbj);
}
