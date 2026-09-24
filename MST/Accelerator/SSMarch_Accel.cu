/* *********************************************************************
 *  SSMarch_Accel.cu
 *
 *  GPU acceleration of the single-site radial march
 *      SSSolverModule::solveSCr
 *  which is the Adams-Bashforth predictor / Adams-Moulton corrector
 *  integration at the heart of the full-potential single-site solver.
 *
 *  WHY THIS ROUTINE
 *  ----------------
 *  On the Fe3Ni full-potential reference case (8 atoms, 1 MPI rank,
 *  LIZ 135, lmax_kkr = lmax_phi = 4), SCF iteration 2 spends
 *
 *      calValenceStates                        2,636 s
 *        solveSingleScattering (6,234 calls)   1,696 s   64 %
 *          -> calPhiLr -> solveSCr             <== here
 *
 *  at roughly 3.5 GFLOP/s, about 7 % of ONE Grace core, with the GH200
 *  idle.  The KKR matrix inversion that one would expect to dominate an
 *  LSMS code is already on the device and is inside the 2.6 % remainder.
 *
 *  WHAT IS AND IS NOT PARALLEL
 *  ---------------------------
 *  The march itself is strictly sequential: sx(ir) depends on sx(ir-1)
 *  and on dsx at the previous four points, and the five corrector
 *  iterations are a dependency chain.  Nothing here tries to break
 *  that, and nothing should -- parallel-in-time on a complex-energy
 *  Adams-Moulton integrator is a research problem, not an engineering
 *  one.
 *
 *  What IS parallel is the angular-momentum index.  One march couples
 *  kmax_phi components through a single kmax_phi x kmax_phi matvec per
 *  corrector, so this kernel runs ONE BLOCK PER MARCH with one thread
 *  per component.  Thread klp owns output component klp and reads
 *  COLUMN klp of V, whose elements are contiguous, so the matvec needs
 *  no cross-thread reduction at all -- only a __syncthreads() between
 *  the pwave update and the matvec.
 *
 *  THE CONTRACTION, AND WHY V IS TRANSPOSED ON UPLOAD
 *  --------------------------------------------------
 *  The Fortran call is
 *      zgemv('t', kmax_phi, kmax_phi, CONE, V, ldV, pwave, 1, ...)
 *  i.e.  vpsum(klp) = sum_klpp V(klpp,klp) * pwave(klpp)
 *  summing over the FIRST index.  A silent transpose here yields
 *  results that look physical and are wrong.
 *
 *  In Fortran V is column-major with klpp contiguous, so column klp is
 *  V[klpp + ldV*klp].  The first version of this kernel had thread klp
 *  walk that column with unit stride -- contiguous PER THREAD, which is
 *  exactly the wrong thing on a GPU.  Across a warp, at summation step
 *  s the threads then touch addresses 16*ldV apart: one cache line per
 *  lane per load, 25 loads per corrector, ~625 serialised transactions.
 *  Measured on the Si case that cost 11.58 ms per march and 822 s of
 *  device time over three SCF iterations -- 99.5 % of all GPU time, and
 *  slower than the CPU path it replaced.
 *
 *  V is therefore TRANSPOSED once per upload into Vt(klp,klpp,ir), so
 *  element (t,s) sits at 2*(t + ldV*s).  At each summation step the
 *  warp reads consecutive complex values: coalesced, and the 10 KB
 *  radial slab stays in L1 across the five correctors.  This is a data
 *  layout change only; the contraction still sums over the first index
 *  of the original V.
 *
 *  THE DIAGONAL SHIFT
 *  ------------------
 *  gaunt_pot = V(:,:,irmn0) + vshift*I, where
 *      vshift = llp1*(1-cm0(ir))/r(ir)^2 + e^2*c2inv + PotShift
 *  carries all of the l and energy dependence (see
 *  SSSolverModule::buildVTables).  It is applied to vpsum, never to the
 *  matrix, so the matrix is read straight from the resident table.
 *
 *  RESIDENCY
 *  ---------
 *  V is (ldV x ldV x nr) and is 10.4 MB for the reference case.  It
 *  depends on (site, spin, potential) only, so it is uploaded when that
 *  key changes rather than per march.
 *
 *  Transfer is NOT the limiting cost here, contrary to what an earlier
 *  version of this file argued from arithmetic.  Nsight on the Si case
 *  reports 60 GB device-to-host in 0.71 s and 10.6 GB host-to-device in
 *  0.60 s across the WHOLE three-iteration run, against 822 s of
 *  kernel.  The keyed upload is still worth having, but do not
 *  re-derive a transfer ceiling from bandwidth estimates -- measure it.
 *
 *  Fortran interface (implicit interface, trailing underscore, every
 *  argument by reference -- same convention as the other .cu files in
 *  this directory):
 *
 *    call query_ssmarch_gpu(bytes_needed, iok)
 *    call init_ssmarch_gpu(nr_max, kmax_phi_max, lmax_phi_max, my_pe)
 *    call ssm_push_vtable_gpu(V, ldv, nr, iregion)
 *    call ssm_push_energy_gpu(r_mesh, bjl, bnl, cm0, lofk, iend, nlj,
 *                             kmax_phi)
 *    call ssm_march_gpu(...)
 *    call finalize_ssmarch_gpu()
 *
 *  Book-keeping:      date                       note
 *                   09-09-26     Initial version (solveSCr offload)
 * ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <complex.h>
#include "cuComplex.h"
#include "acclib.hpp"

/*  Maximum angular-momentum components carried in registers/shared by
 *  the march kernel.  kmax_phi = (lmax_phi+1)^2, so 100 covers
 *  lmax_phi <= 9.                                                    */
#define SSM_KMAX_CAP 100

/*  Adams-Moulton corrector count; matches icmax in solveSCr.         */
#define SSM_ICMAX_CAP 8

/*  Reduction lanes per output component.  One warp per march leaves the
 *  SM with nothing to interleave: measured at ~2,040 cycles per
 *  corrector for ~500 instructions, i.e. ~4 cycles each, which is a
 *  stalled single warp rather than any real work.  Splitting the
 *  summation over SSM_NRED lanes gives the scheduler that many warps.
 *  The state update that follows the reduction is inherently one warp
 *  wide -- only batching whole marches into a block fixes that -- so
 *  this is a partial remedy, not the end of the story.               */
#define SSM_NRED 8

/* -------------------------------------------------------------------
 *  Persistent device state
 * ------------------------------------------------------------------- */
static bool ssm_initialized = false;

static int ssm_nr_max   = 0;   /* max radial points (iend)            */
static int ssm_kmax_max = 0;   /* max kmax_phi                        */
static int ssm_nlj_max  = 0;   /* lmax_phi+2 : second dim of bjl/bnl  */
static int ssm_my_rank  = -1;

static int ssm_ldV      = 0;   /* leading dim of the resident V table */
static int ssm_nr_in    = 0;   /* radial extent of the N0 == 0 table  */
static int ssm_nr_out   = 0;

static double         *d_Vin  = NULL;  /* Vt(klp,klpp,ir), nr_in slabs  */
static double         *d_Vout = NULL;  /* Vt(klp,klpp,ir), nr_out slabs */
static double         *d_Vstg = NULL;  /* upload staging, pre-transpose */
static double         *d_r    = NULL;  /* (nr_max)          real      */
static cuDoubleComplex *d_bjl = NULL;  /* (nr_max x nlj)              */
static cuDoubleComplex *d_bnl = NULL;
static cuDoubleComplex *d_cm0 = NULL;  /* (nr_max)                    */
static int            *d_lofk = NULL;  /* (kmax_max)                  */
static cuDoubleComplex *d_sx  = NULL;  /* (nr_max x kmax_max)         */
static cuDoubleComplex *d_cx  = NULL;
static cuDoubleComplex *d_ds  = NULL;  /* (kmax_max x 4)              */
static cuDoubleComplex *d_dc  = NULL;

static cudaStream_t ssm_stream;

/* ---- small complex helpers, kept explicit for readability -------- */
static __device__ __forceinline__
cuDoubleComplex cmul(cuDoubleComplex a, cuDoubleComplex b)
{ return make_cuDoubleComplex(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }

static __device__ __forceinline__
cuDoubleComplex cadd(cuDoubleComplex a, cuDoubleComplex b)
{ return make_cuDoubleComplex(a.x+b.x, a.y+b.y); }

static __device__ __forceinline__
cuDoubleComplex csub(cuDoubleComplex a, cuDoubleComplex b)
{ return make_cuDoubleComplex(a.x-b.x, a.y-b.y); }

static __device__ __forceinline__
cuDoubleComplex cscale(double s, cuDoubleComplex a)
{ return make_cuDoubleComplex(s*a.x, s*a.y); }

static __device__ __forceinline__
cuDoubleComplex cdiv(cuDoubleComplex a, cuDoubleComplex b)
{
   const double d = b.x*b.x + b.y*b.y;
   return make_cuDoubleComplex((a.x*b.x + a.y*b.y)/d, (a.y*b.x - a.x*b.y)/d);
}

/* -------------------------------------------------------------------
 *  Kernel: transpose each radial slab of V on upload.
 *
 *  src is the Fortran V(klpp,klp,ir), i.e. src[a + ldV*b] with a=klpp.
 *  dst is Vt(klp,klpp,ir) = dst[b + ldV*a], so that the march kernel's
 *  thread b reads consecutive addresses as a advances.  Both reads and
 *  writes here are coalesced along one of the two indices, and this
 *  runs once per V upload (16 per SCF iteration on the real-axis
 *  path), so its own efficiency is irrelevant next to the march.
 * ------------------------------------------------------------------- */
__global__
void ssm_transposeKernel(const double * __restrict__ src,
                         double * __restrict__ dst, int ldV, int nr)
{
   const long ntot = (long)ldV*ldV*nr;
   for (long idx = (long)blockIdx.x*blockDim.x + threadIdx.x; idx < ntot;
             idx += (long)blockDim.x*gridDim.x) {
      const int slab = (int)(idx / ((long)ldV*ldV));
      const int rem  = (int)(idx % ((long)ldV*ldV));
      const int a    = rem % ldV;          /* klpp */
      const int b    = rem / ldV;          /* klp  */
      const long so  = 2*((long)ldV*ldV*slab + a + (long)ldV*b);
      const long dofs= 2*((long)ldV*ldV*slab + b + (long)ldV*a);
      dst[dofs    ] = src[so    ];
      dst[dofs + 1] = src[so + 1];
   }
}

/* -------------------------------------------------------------------
 *  Kernel: one block per march.
 *
 *  Reproduces the do-ir loop of solveSCr operation for operation:
 *  the 4-step Adams-Bashforth predictor, then icmax Adams-Moulton
 *  corrector iterations, with the same j/jm1/jm2/jm3 rotation of the
 *  derivative history.
 *
 *  Layout note: sx and cx are the Fortran arrays sx(iend,kmax_phi),
 *  column-major, so element (ir,klp) is at ir + nr*klp.
 * ------------------------------------------------------------------- */
__global__
void ssm_marchKernel(const double * __restrict__ V, int ldV,
                     const double * __restrict__ rmesh,
                     const cuDoubleComplex * __restrict__ bjl,
                     const cuDoubleComplex * __restrict__ bnl,
                     const cuDoubleComplex * __restrict__ cm0,
                     const int * __restrict__ lofk,
                     int nr, int nlj, int kmax,
                     int N1, int N2, int N0, int nstep, int icmax,
                     double hfac, int llp1, int has_l0,
                     cuDoubleComplex kappa, cuDoubleComplex e2oc2,
                     cuDoubleComplex potshift,
                     cuDoubleComplex * __restrict__ sx,
                     cuDoubleComplex * __restrict__ cx,
                     cuDoubleComplex * __restrict__ dsx,
                     cuDoubleComplex * __restrict__ dcx)
{
   __shared__ cuDoubleComplex s_pwave[SSM_KMAX_CAP];
   __shared__ cuDoubleComplex s_sxir[SSM_KMAX_CAP];
   __shared__ cuDoubleComplex s_cxir[SSM_KMAX_CAP];
   __shared__ cuDoubleComplex s_sxm1[SSM_KMAX_CAP];
   __shared__ cuDoubleComplex s_cxm1[SSM_KMAX_CAP];
   __shared__ cuDoubleComplex s_ds[SSM_KMAX_CAP][4];
   __shared__ cuDoubleComplex s_dc[SSM_KMAX_CAP][4];
   __shared__ cuDoubleComplex s_red[SSM_NRED][SSM_KMAX_CAP];

   const int t    = threadIdx.x;
   const int ty   = threadIdx.y;          /* reduction lane            */
   const int nred = blockDim.y;
   const bool lead = (ty == 0);           /* owns the march state      */
   const bool act  = (t < kmax);
   const int lp = act ? lofk[t] : 0;

   /*  pwave is divided by kappa five times per radial step.  An FP64
    *  complex divide is two FP64 divisions at ~40-80 cycles each, so
    *  the reciprocal is formed once and multiplied thereafter.  This
    *  rounds differently from a divide in the last bit; the CPU path
    *  uses Fortran's complex divide and neither matches the other
    *  exactly in any case.                                           */
   const double kden = kappa.x*kappa.x + kappa.y*kappa.y;
   const cuDoubleComplex inv_kappa =
         make_cuDoubleComplex(kappa.x/kden, -kappa.y/kden);

   /* ---- load the incoming state ---------------------------------- */
   if (lead && act) {
      const int irm1 = N1 - nstep;          /* 0-based below */
      s_sxm1[t] = sx[(irm1-1) + (long)nr*t];
      s_cxm1[t] = cx[(irm1-1) + (long)nr*t];
#pragma unroll
      for (int q = 0; q < 4; ++q) {
         s_ds[t][q] = dsx[t + (long)kmax*q];
         s_dc[t][q] = dcx[t + (long)kmax*q];
      }
   }
   __syncthreads();

   int j = 3, jm1 = 2, jm2 = 1, jm3 = 0;    /* 0-based j=4,3,2,1 */

   for (int ir = N1; nstep > 0 ? ir <= N2 : ir >= N2; ir += nstep) {
      const int ir0   = ir - 1;             /* 0-based radial index   */
      const int irmn0 = ir - N0 - 1;        /* 0-based index into V   */
      const double rr = rmesh[ir0];

      /* ---- predictor --------------------------------------------- */
      if (lead && act) {
         cuDoubleComplex a = cscale(55.0, s_ds[t][j]);
         a = csub(a, cscale(59.0, s_ds[t][jm1]));
         a = cadd(a, cscale(37.0, s_ds[t][jm2]));
         a = csub(a, cscale( 9.0, s_ds[t][jm3]));
         s_sxir[t] = cadd(s_sxm1[t], cscale(hfac, a));

         cuDoubleComplex b = cscale(55.0, s_dc[t][j]);
         b = csub(b, cscale(59.0, s_dc[t][jm1]));
         b = cadd(b, cscale(37.0, s_dc[t][jm2]));
         b = csub(b, cscale( 9.0, s_dc[t][jm3]));
         s_cxir[t] = cadd(s_cxm1[t], cscale(hfac, b));
      }

      /* ---- diagonal shift: llp1*(1-cm0)/r^2 + e2oc2 + PotShift ---- */
      cuDoubleComplex vshift = make_cuDoubleComplex(0.0, 0.0);
      if (has_l0) {
         const cuDoubleComplex one_m_cm0 =
            make_cuDoubleComplex(1.0 - cm0[ir0].x, -cm0[ir0].y);
         vshift = cscale((double)llp1/(rr*rr), one_m_cm0);
         vshift = cadd(vshift, e2oc2);
         vshift = cadd(vshift, potshift);
      }

      /* ---- rotate the history exactly as the Fortran does --------- */
      /*  Fortran, 1-based:  j = jm3
       *                       jm1 = mod(j+2,4)+1
       *                       jm2 = mod(j+1,4)+1
       *                       jm3 = mod(j  ,4)+1
       *  Subtracting one throughout gives the 0-based form below.
       *  Starting from j=3,jm1=2,jm2=1,jm3=0 the first rotation must
       *  yield j=0,jm1=3,jm2=2,jm3=1, which it does.                 */
      j   = jm3;
      jm1 = (j + 3) % 4;
      jm2 = (j + 2) % 4;
      jm3 = (j + 1) % 4;

      /*  V holds complex as adjacent double pairs, so one radial slab
       *  is 2*ldV*ldV doubles, not ldV*ldV.                          */
      const double * __restrict__ Vir = V + (long)2*ldV*ldV*irmn0;

      /*  bjl and bnl depend only on (ir, lp), so they are loaded once
       *  per radial step rather than twice per corrector.            */
      const cuDoubleComplex bj = act ? bjl[ir0 + (long)nr*lp]
                                     : make_cuDoubleComplex(0.0,0.0);
      const cuDoubleComplex bn = act ? bnl[ir0 + (long)nr*lp]
                                     : make_cuDoubleComplex(0.0,0.0);

      for (int icor = 0; icor < icmax; ++icor) {
         /*  pwave is march state, so only the lead lane forms it. */
         if (lead && act) {
            s_pwave[t] = cmul(csub(cmul(s_sxir[t], bn),
                                   cmul(s_cxir[t], bj)), inv_kappa);
         }
         __syncthreads();

         /*  vpsum(t) = sum_s V(s,t)*pwave(s), split over nred lanes.
          *  V is held transposed as Vt(t,s) at 2*(t + ldV*s), so at
          *  every s the warp reads consecutive complex values.  Two
          *  accumulators break the FP64 dependency chain within a
          *  lane; the lanes themselves are independent.            */
         if (act) {
            double vr0=0.0, vr1=0.0, vi0=0.0, vi1=0.0;
            int s = ty;
            for (; s + nred < kmax; s += 2*nred) {
               const long o0 = (long)2*(t + (long)ldV*s);
               const long o1 = (long)2*(t + (long)ldV*(s+nred));
               const double a0r=Vir[o0], a0i=Vir[o0+1];
               const double a1r=Vir[o1], a1i=Vir[o1+1];
               const double b0r=s_pwave[s].x,      b0i=s_pwave[s].y;
               const double b1r=s_pwave[s+nred].x, b1i=s_pwave[s+nred].y;
               vr0 += a0r*b0r - a0i*b0i;  vi0 += a0r*b0i + a0i*b0r;
               vr1 += a1r*b1r - a1i*b1i;  vi1 += a1r*b1i + a1i*b1r;
            }
            for (; s < kmax; s += nred) {
               const long o = (long)2*(t + (long)ldV*s);
               const double ar=Vir[o], ai=Vir[o+1];
               const double br=s_pwave[s].x, bi=s_pwave[s].y;
               vr0 += ar*br - ai*bi;  vi0 += ar*bi + ai*br;
            }
            s_red[ty][t] = make_cuDoubleComplex(vr0+vr1, vi0+vi1);
         }
         __syncthreads();

         if (lead && act) {
            cuDoubleComplex vps = s_red[0][t];
            for (int q = 1; q < nred; ++q) vps = cadd(vps, s_red[q][t]);
            vps = cadd(vps, cmul(vshift, s_pwave[t]));

            const cuDoubleComplex rv = cscale((double)nstep*rr, vps);
            s_ds[t][j] = cmul(rv, bj);
            s_dc[t][j] = cmul(rv, bn);

            cuDoubleComplex a = cscale( 9.0, s_ds[t][j]);
            a = cadd(a, cscale(19.0, s_ds[t][jm1]));
            a = csub(a, cscale( 5.0, s_ds[t][jm2]));
            a = cadd(a,             s_ds[t][jm3]);
            s_sxir[t] = cadd(s_sxm1[t], cscale(hfac, a));

            cuDoubleComplex b = cscale( 9.0, s_dc[t][j]);
            b = cadd(b, cscale(19.0, s_dc[t][jm1]));
            b = csub(b, cscale( 5.0, s_dc[t][jm2]));
            b = cadd(b,             s_dc[t][jm3]);
            s_cxir[t] = cadd(s_cxm1[t], cscale(hfac, b));
         }
         __syncthreads();
      }

      if (lead && act) {
         sx[ir0 + (long)nr*t] = s_sxir[t];
         cx[ir0 + (long)nr*t] = s_cxir[t];
         s_sxm1[t] = s_sxir[t];
         s_cxm1[t] = s_cxir[t];
      }
   }

   /* ---- hand the derivative history back ------------------------- */
   __syncthreads();
   if (lead && act) {
#pragma unroll
      for (int q = 0; q < 4; ++q) {
         dsx[t + (long)kmax*q] = s_ds[t][q];
         dcx[t + (long)kmax*q] = s_dc[t][q];
      }
   }
}

/* ===================================================================
 *  query_ssmarch_gpu(bytes_needed, ok)
 *
 *  Runtime gate, and it is OPT-IN: the device march is OFF unless
 *  MUST_SSMARCH_GPU is set to 1/y/t.
 *
 *  WHY OFF BY DEFAULT.  Measured on the Fe3Ni full-potential
 *  spin-canted case, same binary, same GH200 node, only this variable
 *  changed:
 *
 *      calValenceStates, steady-state iteration
 *         march on the GPU   2013.5 s
 *         march on the CPU   1430.9 s     <-- 582.6 s FASTER
 *
 *      inner march, 190,500 of them
 *         GPU   5.16 ms each   (983 s)
 *         CPU   2.10 ms each   (401 s)    <-- CPU is 2.45x faster
 *
 *  computeNewPotential was 42.6 s against 42.8 s across the pair, so
 *  that is a clean single-variable comparison.
 *
 *  This is precisely the configuration the plan's Gate 3 warned about:
 *  "a single march on the GPU is slower than on the CPU; set the
 *  threshold at ~64 marches".  It meant 64 marches PER LAUNCH, and
 *  this kernel runs one march per launch -- one block, a handful of
 *  warps, on one SM of 132, driving a 5,000-long dependency chain with
 *  barriers at every corrector.  Meanwhile the CPU march, once the
 *  Gaunt contraction is hoisted, is five small zgemv calls per radial
 *  step on data that stays in L1 of a 3.1 GHz core.  The CPU simply
 *  wins that shape of problem.
 *
 *  The code is kept, and kept correct, for two reasons: batching whole
 *  marches into a block would change the verdict, and so might a
 *  larger kmax_phi or different hardware.  Turn it on with
 *  MUST_SSMARCH_GPU=1 and compare -- both paths now report their own
 *  cost, so one run is enough.
 *
 *  Any probe failure clears the CUDA error state so it cannot poison
 *  later CUDA use.
 * =================================================================== */
extern "C"
void query_ssmarch_gpu_(long *bytes_needed, int *ok)
{
   *ok = 0;

   const char *env = getenv("MUST_SSMARCH_GPU");
   if (env == NULL) {
      return;                       /* off unless explicitly requested */
   }
   if (!(env[0] == '1' || env[0] == 'y' || env[0] == 'Y' ||
         env[0] == 't' || env[0] == 'T')) {
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
   if (*bytes_needed > 0 && (size_t)(*bytes_needed) > freeb) return;

   *ok = 1;
}

/* ===================================================================
 *  init_ssmarch_gpu(nr_max, kmax_max, nlj_max, ldV, nr_in, nr_out, my_pe)
 * =================================================================== */
extern "C"
void init_ssmarch_gpu_(int *nr_max, int *kmax_max, int *nlj_max,
                       int *ldV, int *nr_in, int *nr_out, int *my_pe)
{
   if (ssm_initialized) return;

   ssm_nr_max   = *nr_max;
   ssm_kmax_max = *kmax_max;
   ssm_nlj_max  = *nlj_max;
   ssm_ldV      = *ldV;
   ssm_nr_in    = *nr_in;
   ssm_nr_out   = *nr_out;
   ssm_my_rank  = *my_pe;

   if (ssm_kmax_max > SSM_KMAX_CAP) {
      fprintf(stderr, "\nError in init_ssmarch_gpu: kmax_phi %d exceeds "
                      "SSM_KMAX_CAP %d\n", ssm_kmax_max, SSM_KMAX_CAP);
      exit(EXIT_FAILURE);
   }
   if (ssm_nr_max < 1 || ssm_kmax_max < 1 || ssm_ldV < 1) {
      fprintf(stderr, "\nError in init_ssmarch_gpu: bad dimensions "
                      "(nr %d, kmax %d, ldV %d)\n",
                      ssm_nr_max, ssm_kmax_max, ssm_ldV);
      exit(EXIT_FAILURE);
   }

   const long nr  = ssm_nr_max;
   const long km  = ssm_kmax_max;
   const long ldv = ssm_ldV;

   checkCudaErrors(cudaStreamCreate(&ssm_stream));
   checkCudaErrors(cudaMalloc((void**)&d_Vin,
                   2*ldv*ldv*(long)ssm_nr_in*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_Vout,
                   2*ldv*ldv*(long)ssm_nr_out*sizeof(double)));
   {  const long nstg = (ssm_nr_in > ssm_nr_out) ? ssm_nr_in : ssm_nr_out;
      checkCudaErrors(cudaMalloc((void**)&d_Vstg,
                      2*ldv*ldv*nstg*sizeof(double)));
   }
   checkCudaErrors(cudaMalloc((void**)&d_r,    nr*sizeof(double)));
   checkCudaErrors(cudaMalloc((void**)&d_bjl,
                   nr*(long)ssm_nlj_max*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_bnl,
                   nr*(long)ssm_nlj_max*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_cm0,  nr*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_lofk, km*sizeof(int)));
   checkCudaErrors(cudaMalloc((void**)&d_sx,   nr*km*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_cx,   nr*km*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_ds,   km*4*sizeof(cuDoubleComplex)));
   checkCudaErrors(cudaMalloc((void**)&d_dc,   km*4*sizeof(cuDoubleComplex)));

   ssm_initialized = true;
}

/* ===================================================================
 *  ssm_push_vtable_gpu(V, ldv, nr, iregion)
 *
 *  iregion = 0 : the N0 == 0 table;  1 : the truncated-region table.
 *  Called only when the (site, spin, SCF iteration) key changes -- see
 *  SSMarchModule.  Pushing per march would move ~129 GB per SCF
 *  iteration and erase the whole benefit.
 * =================================================================== */
extern "C"
void ssm_push_vtable_gpu_(double _Complex *V, int *ldv, int *nr, int *iregion)
{
   if (!ssm_initialized) {
      fprintf(stderr, "\nError in ssm_push_vtable_gpu: not initialised\n");
      exit(EXIT_FAILURE);
   }
   if (*ldv != ssm_ldV) {
      fprintf(stderr, "\nError in ssm_push_vtable_gpu: ldv %d/%d\n",
                      *ldv, ssm_ldV);
      exit(EXIT_FAILURE);
   }
   const int nrmax = (*iregion == 0) ? ssm_nr_in : ssm_nr_out;
   if (*nr < 1 || *nr > nrmax) {
      fprintf(stderr, "\nError in ssm_push_vtable_gpu: nr %d/%d "
                      "(region %d)\n", *nr, nrmax, *iregion);
      exit(EXIT_FAILURE);
   }
   double *dst = (*iregion == 0) ? d_Vin : d_Vout;
   checkCudaErrors(cudaMemcpyAsync(d_Vstg, V,
                   2*(size_t)ssm_ldV*ssm_ldV*(*nr)*sizeof(double),
                   cudaMemcpyHostToDevice, ssm_stream));
   {
      const long ntot = (long)ssm_ldV*ssm_ldV*(*nr);
      int tpb = 256;
      long nb  = (ntot + tpb - 1)/tpb;
      if (nb > 65535) nb = 65535;
      if (nb < 1)     nb = 1;
      ssm_transposeKernel<<<(int)nb, tpb, 0, ssm_stream>>>
                         (d_Vstg, dst, ssm_ldV, *nr);
      checkCudaErrors(cudaPeekAtLastError());
   }
   checkCudaErrors(cudaStreamSynchronize(ssm_stream));
}

/* ===================================================================
 *  ssm_push_energy_gpu(rmesh, bjl, bnl, cm0, lofk, iend, nlj, kmax)
 *
 *  The per-(site,energy) tables.  Small: ~0.5 MB for the reference
 *  case, against 10.4 MB for V.
 * =================================================================== */
extern "C"
void ssm_push_energy_gpu_(double *rmesh, double _Complex *bjl,
                          double _Complex *bnl, double _Complex *cm0,
                          int *lofk, int *iend, int *nlj, int *kmax)
{
   if (!ssm_initialized) {
      fprintf(stderr, "\nError in ssm_push_energy_gpu: not initialised\n");
      exit(EXIT_FAILURE);
   }
   if (*iend > ssm_nr_max || *nlj > ssm_nlj_max || *kmax > ssm_kmax_max) {
      fprintf(stderr, "\nError in ssm_push_energy_gpu: iend %d/%d, nlj "
                      "%d/%d, kmax %d/%d\n", *iend, ssm_nr_max, *nlj,
                      ssm_nlj_max, *kmax, ssm_kmax_max);
      exit(EXIT_FAILURE);
   }
   const size_t n = (size_t)(*iend);
   checkCudaErrors(cudaMemcpyAsync(d_r, rmesh, n*sizeof(double),
                   cudaMemcpyHostToDevice, ssm_stream));
   /*  bjl and bnl are the module arrays bjl(iend_max, 0:lmax_phi+1):
    *  their leading dimension is iend_max, NOT the current site's
    *  iend, so the source pitch is ssm_nr_max.  Using iend here would
    *  silently shear the table for any site with iend < iend_max.   */
   checkCudaErrors(cudaMemcpy2DAsync(d_bjl, (size_t)ssm_nr_max*sizeof(cuDoubleComplex),
                   bjl, (size_t)ssm_nr_max*sizeof(cuDoubleComplex),
                   n*sizeof(cuDoubleComplex), (size_t)(*nlj),
                   cudaMemcpyHostToDevice, ssm_stream));
   checkCudaErrors(cudaMemcpy2DAsync(d_bnl, (size_t)ssm_nr_max*sizeof(cuDoubleComplex),
                   bnl, (size_t)ssm_nr_max*sizeof(cuDoubleComplex),
                   n*sizeof(cuDoubleComplex), (size_t)(*nlj),
                   cudaMemcpyHostToDevice, ssm_stream));
   checkCudaErrors(cudaMemcpyAsync(d_cm0, cm0, n*sizeof(cuDoubleComplex),
                   cudaMemcpyHostToDevice, ssm_stream));
   checkCudaErrors(cudaMemcpyAsync(d_lofk, lofk, (size_t)(*kmax)*sizeof(int),
                   cudaMemcpyHostToDevice, ssm_stream));
   checkCudaErrors(cudaStreamSynchronize(ssm_stream));
}

/* ===================================================================
 *  ssm_march_gpu(...)
 *
 *  One march.  sx and cx are the host arrays sx(nr,kmax); only the
 *  row N1-nstep is uploaded and only rows [min(N1,N2), max(N1,N2)] are
 *  brought back, so the transfer is the march's own footprint and not
 *  the whole array.
 * =================================================================== */
extern "C"
void ssm_march_gpu_(int *N1, int *N2, int *N0, int *nstep, int *icmax,
                    int *nr, int *kmax, int *nlj, int *llp1, int *has_l0,
                    double *hfac, double _Complex *kappa,
                    double _Complex *e2oc2, double _Complex *potshift,
                    double _Complex *sx, double _Complex *cx,
                    double _Complex *dsx, double _Complex *dcx)
{
   if (!ssm_initialized) {
      fprintf(stderr, "\nError in ssm_march_gpu: not initialised\n");
      exit(EXIT_FAILURE);
   }
   const int nrl = *nr, km = *kmax;
   const int irm1 = *N1 - *nstep;
   if (irm1 < 1 || irm1 > nrl) {
      fprintf(stderr, "\nError in ssm_march_gpu: seed row %d outside "
                      "[1,%d]\n", irm1, nrl);
      exit(EXIT_FAILURE);
   }

   const int lo = (*N1 < *N2) ? *N1 : *N2;
   const int hi = (*N1 < *N2) ? *N2 : *N1;
   const size_t esz = sizeof(cuDoubleComplex);

   /* seed row, plus the derivative history */
   checkCudaErrors(cudaMemcpy2DAsync(d_sx + (irm1-1), (size_t)ssm_nr_max*esz,
                   (char*)sx + (size_t)(irm1-1)*esz, (size_t)nrl*esz,
                   esz, (size_t)km, cudaMemcpyHostToDevice, ssm_stream));
   checkCudaErrors(cudaMemcpy2DAsync(d_cx + (irm1-1), (size_t)ssm_nr_max*esz,
                   (char*)cx + (size_t)(irm1-1)*esz, (size_t)nrl*esz,
                   esz, (size_t)km, cudaMemcpyHostToDevice, ssm_stream));
   checkCudaErrors(cudaMemcpyAsync(d_ds, dsx, (size_t)km*4*esz,
                   cudaMemcpyHostToDevice, ssm_stream));
   checkCudaErrors(cudaMemcpyAsync(d_dc, dcx, (size_t)km*4*esz,
                   cudaMemcpyHostToDevice, ssm_stream));

   const double *V = (*N0 == 0) ? d_Vin : d_Vout;

   /*  Fortran complex(kind=CmplxKind) and cuDoubleComplex are both two
    *  contiguous doubles, real part first, so the scalars are
    *  reinterpreted rather than decomposed: <complex.h>'s creal/cimag
    *  are not available to the C++ compiler nvcc uses for this file. */

   /*  blockDim.x covers the angular components, blockDim.y the
    *  reduction lanes.  For kmax_phi = 25 that is 32 x 8 = 256
    *  threads, i.e. eight warps for the scheduler to interleave
    *  instead of the single warp the first version ran.            */
   int tpx = 32;
   while (tpx < km) tpx += 32;
   const dim3 blk(tpx, SSM_NRED, 1);

   ssm_marchKernel<<<1, blk, 0, ssm_stream>>>(
        V, ssm_ldV, d_r, d_bjl, d_bnl, d_cm0, d_lofk,
        ssm_nr_max, *nlj, km, *N1, *N2, *N0, *nstep, *icmax,
        *hfac, *llp1, *has_l0,
        *reinterpret_cast<const cuDoubleComplex*>(kappa),
        *reinterpret_cast<const cuDoubleComplex*>(e2oc2),
        *reinterpret_cast<const cuDoubleComplex*>(potshift),
        d_sx, d_cx, d_ds, d_dc);
   checkCudaErrors(cudaPeekAtLastError());

   /* bring back only the rows this march wrote */
   const size_t wrows = (size_t)(hi - lo + 1);
   checkCudaErrors(cudaMemcpy2DAsync((char*)sx + (size_t)(lo-1)*esz,
                   (size_t)nrl*esz,
                   d_sx + (lo-1), (size_t)ssm_nr_max*esz,
                   wrows*esz, (size_t)km, cudaMemcpyDeviceToHost, ssm_stream));
   checkCudaErrors(cudaMemcpy2DAsync((char*)cx + (size_t)(lo-1)*esz,
                   (size_t)nrl*esz,
                   d_cx + (lo-1), (size_t)ssm_nr_max*esz,
                   wrows*esz, (size_t)km, cudaMemcpyDeviceToHost, ssm_stream));
   checkCudaErrors(cudaMemcpyAsync(dsx, d_ds, (size_t)km*4*esz,
                   cudaMemcpyDeviceToHost, ssm_stream));
   checkCudaErrors(cudaMemcpyAsync(dcx, d_dc, (size_t)km*4*esz,
                   cudaMemcpyDeviceToHost, ssm_stream));
   checkCudaErrors(cudaStreamSynchronize(ssm_stream));
}

/* ===================================================================
 *  finalize_ssmarch_gpu()
 * =================================================================== */
extern "C"
void finalize_ssmarch_gpu_(void)
{
   if (!ssm_initialized) return;
   cudaFree(d_Vin);  d_Vin = NULL;
   cudaFree(d_Vout); d_Vout = NULL;
   cudaFree(d_Vstg); d_Vstg = NULL;
   cudaFree(d_r);    d_r   = NULL;
   cudaFree(d_bjl);  d_bjl = NULL;
   cudaFree(d_bnl);  d_bnl = NULL;
   cudaFree(d_cm0);  d_cm0 = NULL;
   cudaFree(d_lofk); d_lofk= NULL;
   cudaFree(d_sx);   d_sx  = NULL;
   cudaFree(d_cx);   d_cx  = NULL;
   cudaFree(d_ds);   d_ds  = NULL;
   cudaFree(d_dc);   d_dc  = NULL;
   cudaStreamDestroy(ssm_stream);
   ssm_initialized = false;
}
