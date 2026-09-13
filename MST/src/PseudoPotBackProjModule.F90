!  *******************************************************************
!  *                                                                 *
!  *  PseudoPotBackProjModule                                        *
!  *                                                                 *
!  *  GPU offload of the "BackProjection" stage of                    *
!  *     PotentialGenerationModule::calFFTPseudoPot                   *
!  *  which is the k-space summation loop of                          *
!  *     PotentialGenerationModule::calRadialInterpolation            *
!  *                                                                 *
!  *  Motivation                                                     *
!  *  ----------                                                     *
!  *  calFFTPseudoPot is named after its Fourier transform but does   *
!  *  not spend its time there.  Measured on the Fe3Ni full-potential  *
!  *  reference case (8 atoms, 1 MPI rank, 64^3 FFT grid), per SCF     *
!  *  iteration:                                                      *
!  *                                                                 *
!  *     Time in UniformGrid            2.707 s     3.92 %            *
!  *     Time in performTransformR2C    0.002 s     0.003 %           *
!  *     Time in Dipole                 0.058 s     0.08 %            *
!  *     Time in BackProjection        66.294 s    96.00 %  <== here  *
!  *                                                                 *
!  *  The 66 s is a single-threaded, scalar triple loop over          *
!  *  (atom, k-point, radial node) that, for every one of the         *
!  *  na*numk = 8*135168 (atom, k) pairs, evaluates a conjugated      *
!  *  spherical harmonic and a spherical Bessel recursion over nnr    *
!  *  radial nodes.  It runs at ~0.6 GFLOP/s -- about 2 % of ONE      *
!  *  Grace core -- with the GH200 idle.                              *
!  *                                                                 *
!  *  Method                                                         *
!  *  ------                                                         *
!  *  The CPU loop accumulates, for atom a,                          *
!  *                                                                 *
!  *     pv(n,jl) = i2l(l)/r_n^l * SUM_k j_l(|k| r_n) * C_a(k,jl)    *
!  *                                                                 *
!  *     C_a(k,jl) = Re[kf_a(k)] * conjg(Y_lm(k))        l even       *
!  *               = i*Im[kf_a(k)] * conjg(Y_lm(k))      l odd        *
!  *     kf_a(k)   = exp(i k.R_a) * rho~(k) * |k|^kpow                *
!  *     i2l(l)    = dummy * i^l                                      *
!  *                                                                 *
!  *  j_l does not depend on m, and the packed index                  *
!  *     jl(l,m) = (l+1)(l+2)/2 - l + m ,  m = 0..l                   *
!  *  is CONTIGUOUS in m for fixed l, so the k-sum is a matrix-matrix *
!  *  product for each l separately:                                  *
!  *                                                                 *
!  *     pv(1:nnr, jl0(l):jl0(l)+l) = B_l(nnr x numk) * C(numk, ...)  *
!  *     B_l(n,k) = j_l(|k| r_n)                       REAL           *
!  *                                                                 *
!  *  B_l is real and C complex, so each l is issued as TWO real      *
!  *  DGEMMs -- the same real/imaginary split that                    *
!  *  DensityOnGridModule already uses for the density interpolation. *
!  *  The k-sum is tiled so the Bessel table stays bounded            *
!  *  independently of the FFT grid; tiles accumulate with beta = 1.  *
!  *                                                                 *
!  *  Two structural savings come for free with this restructuring    *
!  *  and are the reason the host side is worth touching at all:      *
!  *                                                                 *
!  *    - conjg(Y_lm(k)) depends only on k, not on the atom, yet the  *
!  *      CPU loop calls calYlmConjg inside the atom loop -- na times *
!  *      per k-point.  Here it is evaluated once per k, on the host, *
!  *      with the existing trusted routine, so calYlmConjg never has *
!  *      to be ported to CUDA.                                       *
!  *    - r_interp is generated from the radial mesh ends, the        *
!  *      inscribed sphere radius and d_ir, so it is a per-SPECIES    *
!  *      quantity, not a per-atom one.  Atoms sharing                *
!  *      (lmax, nnr, r_interp) share the Bessel table, so the        *
!  *      expensive kernel runs once per (tile, species).  Fe3Ni's 8  *
!  *      atoms have only 2 distinct r_interp.                        *
!  *                                                                 *
!  *  The CPU path in calRadialInterpolation is NOT removed.  It      *
!  *  stays reachable at runtime, permanently, because it is the      *
!  *  correctness reference, it is what a CPU-architecture build      *
!  *  runs, and it is what keeps the baseline measurement above       *
!  *  reproducible.  MUST_PSEUDOPOT_GPU=0 selects it without a        *
!  *  rebuild.                                                        *
!  *                                                                 *
!  *  Gating                                                         *
!  *  ------                                                         *
!  *  Two independent gates, unlike DensityOnGridModule which decides *
!  *  at compile time only:                                          *
!  *     compile time : #ifdef ACCEL                                 *
!  *     run time     : query_pseudopot_gpu() -- device present, can  *
!  *                    be selected, enough free memory for THIS      *
!  *                    problem size, cuBLAS handle creatable, and    *
!  *                    MUST_PSEUDOPOT_GPU not set to 0.             *
!  *  A rank that fails either gate falls back to the CPU path        *
!  *  instead of aborting, which is what makes an ACCEL binary safe   *
!  *  to launch on a CPU-only node and safe to oversubscribe on a     *
!  *  one-GPU-per-node machine such as Vista.                        *
!  *                                                                 *
!  *  version 1.0, Sep 9, 2026                                       *
!  *                                                                 *
!  *******************************************************************
module PseudoPotBackProjModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : ZERO, ONE, CZERO
   use TimerModule, only : getTime
   use CmdLineOptionModule, only : getCmdLineOption
!
   implicit none
!
public :: initPseudoPotBackProj,     &
          endPseudoPotBackProj,      &
          isPseudoPotBackProjGPU,    &
          calRadialInterpGPU,        &
          getPseudoPotBackProjTime,  &
          printPseudoPotBackProjInfo
!
private
!
!  ===================================================================
!  Bessel-table memory budget, in MB.  The k-tile width KT is chosen so
!  that B(nnr_max, KT, 0:lmax_max) fits in it, which is what bounds the
!  device footprint independently of the FFT grid: numk grows as the
!  cube of the grid dimension (135168 at 64^3, 1064960 at 128^3) but KT
!  does not move.
!  ===================================================================
   integer (kind=IntKind), parameter :: BesselBudgetMB = 512
   integer (kind=IntKind), parameter :: KT_min = 2048
   integer (kind=IntKind), parameter :: KT_max = 65536
!
!  Device memory head room demanded of the runtime gate, in MB, on top
!  of the buffers this module allocates (cuBLAS work space, context).
   integer (kind=IntKind), parameter :: SlackMB = 256
!
   logical :: Initialized = .false.
   logical :: useGPU      = .false.
!
   integer (kind=IntKind) :: KT       = 0    ! k-tile width
   integer (kind=IntKind) :: numk_max = 0
   integer (kind=IntKind) :: kmax_max = 0
   integer (kind=IntKind) :: nnr_max  = 0
   integer (kind=IntKind) :: jmax_max = 0
   integer (kind=IntKind) :: lmax_max = 0
   integer (kind=IntKind) :: na_max   = 0
   integer (kind=IntKind) :: my_pe    = 0
   integer (kind=IntKind) :: print_level = -1
!
   integer (kind=IntKind) :: num_calls = 0
!
   real (kind=RealKind) :: t_backproj = ZERO
   real (kind=RealKind) :: dev_MB     = ZERO
!
!  ===================================================================
!  Host staging for one k-tile.  ylm_h is (KT, kmax_max) -- k leading --
!  so that the device coefficient kernel, which runs one thread per
!  k-point, loads consecutive addresses across a warp.
!  ===================================================================
   real (kind=RealKind), allocatable :: kvec_h(:,:)
   complex (kind=CmplxKind), allocatable :: fftc_h(:)
   complex (kind=CmplxKind), allocatable :: ylm_h(:,:)
   complex (kind=CmplxKind), allocatable :: ylm_1(:)
!
!  Species grouping: atoms sharing (lmax, nnr, r_interp) share a Bessel
!  table.  a_order is the processing order, with the members of a group
!  adjacent and its representative first; a_group(ia) is that
!  representative, and is what the device is told the resident table
!  should belong to.
   integer (kind=IntKind), allocatable :: a_order(:)
   integer (kind=IntKind), allocatable :: a_group(:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPseudoPotBackProj(numk, lmaxmax, nnrmax, jmaxmax,     &
                                    namax, mype, iprint)
!  ===================================================================
!  numk    : number of local reciprocal grid points that will be summed
!  lmaxmax : maximum potential lmax over the atoms involved
!  nnrmax  : maximum number of radial interpolation nodes (n_interp*nr_int)
!  jmaxmax : maximum number of (l, m>=0) components
!  namax   : maximum number of atoms handled in one call (LocalNumAtoms
!            for serial FFT, GlobalNumAtoms for parallel FFT)
!  mype    : MPI rank
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: numk, lmaxmax, nnrmax
   integer (kind=IntKind), intent(in) :: jmaxmax, namax, mype
   integer (kind=IntKind), intent(in), optional :: iprint
!
   integer (kind=IntKind) :: kt_try, iok
!
   logical :: cpu_only
   integer (kind=IntKind) :: kt_l, km_l, nr_l, jm_l, lm_l, na_l, mp_l
!
   integer (kind=8) :: bytes_needed
!
   if (Initialized) then
      return
   endif
!
   if (present(iprint)) then
      print_level = iprint
   else
      print_level = -1
   endif
!
   if (numk < 1) then
      call ErrorHandler('initPseudoPotBackProj','invalid numk',numk)
   else if (nnrmax < 1) then
      call ErrorHandler('initPseudoPotBackProj','invalid nnr_max',nnrmax)
   else if (jmaxmax < 1) then
      call ErrorHandler('initPseudoPotBackProj','invalid jmax_max',jmaxmax)
   else if (lmaxmax < 0) then
      call ErrorHandler('initPseudoPotBackProj','invalid lmax_max',lmaxmax)
   else if (namax < 1) then
      call ErrorHandler('initPseudoPotBackProj','invalid na_max',namax)
   endif
!
   numk_max = numk
   lmax_max = lmaxmax
   nnr_max  = nnrmax
   jmax_max = jmaxmax
   na_max   = namax
   my_pe    = mype
   kmax_max = (lmax_max+1)*(lmax_max+1)
!
!  ===================================================================
!  k-tile width: largest KT for which B(nnr_max, KT, 0:lmax_max) stays
!  inside the budget, clamped and never larger than the k range itself.
!  ===================================================================
   kt_try = (BesselBudgetMB*1024*1024)/(8*nnr_max*(lmax_max+1))
   KT = min(kt_try, KT_max)
   KT = max(KT, KT_min)
   KT = min(KT, numk_max)
   if (KT < 1) then
      KT = 1
   endif
!
!  ===================================================================
!  Device footprint of the buffers allocated by
!  init_pseudopot_backproj_gpu, used by the runtime memory gate.
!  ===================================================================
   bytes_needed = 8_8*3_8*int(KT,8)                                   & ! kvec
                + 16_8*int(KT,8)                                      & ! fftc
                + 16_8*int(KT,8)*int(kmax_max,8)                      & ! ylm
                + 8_8*2_8*int(KT,8)                                   & ! kmag,k2
                + 8_8*2_8*int(KT,8)*int(jmax_max,8)                   & ! Cr,Ci
                + 8_8*int(nnr_max,8)*int(KT,8)*int(lmax_max+1,8)      & ! B
                + 8_8*2_8*int(nnr_max,8)*int(jmax_max,8)*int(na_max,8)& ! pvr,pvi
                + 8_8*int(nnr_max,8)*int(na_max,8)                    & ! rint
                + 4_8*2_8*int(jmax_max,8)                             & ! lofj,kofj
                + 16_8*int(nnr_max,8)*int(jmax_max,8)                   ! pv
!
   dev_MB = real(bytes_needed,kind=RealKind)/(1024.0d0*1024.0d0)
   bytes_needed = bytes_needed + int(SlackMB,8)*1024_8*1024_8
!
   cpu_only = ( getCmdLineOption('Run on CPU without Acceleration') == 0 )
!
   useGPU = .false.
#ifdef ACCEL
!  ===================================================================
!  Precedence, highest first:
!     1. the -cpu / --cpu-only / --cpu_only command-line flag.  This is
!        the project-wide switch honoured by ClusterMatrixModule, and it
!        means "use no GPU at all".  It OVERRIDES the environment
!        variable below, so a user who asks for a CPU run gets one even
!        if a per-module override is set in their environment.
!     2. the per-module environment variable, checked inside the device
!        probe.  This exists for A/B measurement: --cpu-only is
!        all-or-nothing, whereas comparing one offload at a time is how
!        the radial march was found to be slower than the CPU.
!     3. the device probe itself -- is there a usable GPU with room.
!  ===================================================================
!  A failure of the probe is not fatal either: the caller keeps the CPU
!  reference path.
!  ===================================================================
   iok = 0
   if (.not.cpu_only) then
!     ----------------------------------------------------------------
      call query_pseudopot_gpu(bytes_needed, iok)
!     ----------------------------------------------------------------
   endif
   if (iok == 1) then
      kt_l = KT;       km_l = kmax_max; nr_l = nnr_max
      jm_l = jmax_max; lm_l = lmax_max; na_l = na_max
      mp_l = my_pe
!     ----------------------------------------------------------------
      call init_pseudopot_backproj_gpu(kt_l, km_l, nr_l, jm_l, lm_l,   &
                                       na_l, mp_l)
!     ----------------------------------------------------------------
      useGPU = .true.
   else if (print_level >= 0) then
      if (cpu_only) then
!        -------------------------------------------------------------
         call WarningHandler('initPseudoPotBackProj',                  &
                 '-cpu/--cpu-only given; the CPU path is used')
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call WarningHandler('initPseudoPotBackProj',                  &
                 'GPU probe failed or was disabled; using the CPU path')
!        -------------------------------------------------------------
      endif
   endif
#endif
!
!  ===================================================================
!  Host staging is allocated only when the device path was actually
!  taken.  ylm_h alone is 42 MB at KT = 32768 and kmax = 81, and a rank
!  that fell back to the CPU path has no use for any of it.
!  ===================================================================
   if (useGPU) then
      allocate( kvec_h(3,KT), fftc_h(KT), ylm_h(KT,kmax_max) )
      allocate( ylm_1(kmax_max) )
      allocate( a_order(na_max), a_group(na_max) )
   endif
!
   t_backproj = ZERO
   num_calls  = 0
   Initialized = .true.
!
   if (print_level >= 0) then
!     ----------------------------------------------------------------
      call printPseudoPotBackProjInfo()
!     ----------------------------------------------------------------
   endif
!
   end subroutine initPseudoPotBackProj
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endPseudoPotBackProj()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
      return
   endif
!
#ifdef ACCEL
   if (useGPU) then
!     ----------------------------------------------------------------
      call finalize_pseudopot_backproj_gpu()
!     ----------------------------------------------------------------
   endif
#endif
!
   if (allocated(kvec_h))  deallocate( kvec_h )
   if (allocated(fftc_h))  deallocate( fftc_h )
   if (allocated(ylm_h))   deallocate( ylm_h )
   if (allocated(ylm_1))   deallocate( ylm_1 )
   if (allocated(a_order)) deallocate( a_order )
   if (allocated(a_group)) deallocate( a_group )
!
   useGPU = .false.
   Initialized = .false.
!
   end subroutine endPseudoPotBackProj
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isPseudoPotBackProjGPU() result(y)
!  ===================================================================
   implicit none
   logical :: y
!
   y = ( Initialized .and. useGPU )
!
   end function isPseudoPotBackProjGPU
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPseudoPotBackProjTime() result(t)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: t
!
   t = t_backproj
!
   end function getPseudoPotBackProjTime
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPseudoPotBackProjInfo()
!  ===================================================================
   implicit none
!
   write(6,'(/,a)') 'PseudoPotential back-projection: batched k-space summation'
   write(6,'(a,i12)')   '   local k-points      numk       = ', numk_max
   write(6,'(a,i12)')   '   radial nodes        nnr_max    = ', nnr_max
   write(6,'(a,i12)')   '   L components        jmax_max   = ', jmax_max
   write(6,'(a,i12)')   '   potential lmax      lmax_max   = ', lmax_max
   write(6,'(a,i12)')   '   atoms per call      na_max     = ', na_max
   write(6,'(a,i12)')   '   k-tile              KT         = ', KT
   write(6,'(a,i12)')   '   k-tiles per call               = ',        &
                        (numk_max+KT-1)/KT
   write(6,'(a,f12.2)') '   device memory       MB         = ', dev_MB
   if (useGPU) then
      write(6,'(a)')    '   evaluation path                = GPU (CUDA + cuBLAS)'
   else
      write(6,'(a)')    '   evaluation path                = CPU (reference)'
   endif
   write(6,'(a)') ' '
!
   end subroutine printPseudoPotBackProjInfo
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRadialInterpGPU(p_fft_c, k_first, k_last, kpow, dummy, &
                                 na, posi_a, lmax_a, nnr_a, jmax_a,     &
                                 r_interp_a, ld_rint, w_out, ld_w)
!  ===================================================================
!  Device evaluation of
!
!     w_out(1:nnr*jmax, ia)  <-  pv_interp(1:nnr, 1:jmax)
!
!  for every atom ia = 1..na, summed over k = k_first .. k_last.
!
!  p_fft_c    : rho~(k), indexed over the FULL local k range
!  k_first    : first k index to sum (idk0+1; the k = 0 point is
!               excluded for kpow < 0)
!  k_last     : last  k index to sum (numk_local)
!  kpow       : -2 (Coulomb) or 0
!  dummy      : the real prefactor, 2*(4 pi)^2 or 4 pi
!  posi_a     : atom position relative to the FFT grid origin
!  lmax_a/nnr_a/jmax_a : per-atom dimensions
!  r_interp_a : (ld_rint, na) radial interpolation nodes
!  w_out      : (ld_w, na) destination, written packed with leading
!               dimension nnr_a(ia) -- the layout calRadialProjection
!               expects.  Untouched beyond nnr*jmax, so a caller that
!               pre-zeroed it keeps its zeros.
!  ===================================================================
   use ParallelFFTModule, only : getGridPointCoord
!
   use SphericalHarmonicsModule, only : calYlmConjg
!
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: k_first, k_last, kpow, na
   integer (kind=IntKind), intent(in) :: ld_rint, ld_w
   integer (kind=IntKind), intent(in) :: lmax_a(na), nnr_a(na), jmax_a(na)
!
   real (kind=RealKind), intent(in) :: dummy
   real (kind=RealKind), intent(in) :: posi_a(3,na)
   real (kind=RealKind), intent(in) :: r_interp_a(ld_rint,na)
!
   complex (kind=CmplxKind), intent(in) :: p_fft_c(*)
   complex (kind=CmplxKind), intent(inout) :: w_out(ld_w,*)
!
   integer (kind=IntKind) :: ia, ib, io, i, ng, nk, k0, kl
   integer (kind=IntKind) :: lmax_l, kmax_l, nk_total, ngroup
   integer (kind=IntKind) :: na_l, kp_l, nk_l, km_l, kt_l, ld_l
   integer (kind=IntKind) :: ia_l, ib_l
!
   real (kind=RealKind) :: t0, dmy_l
   real (kind=RealKind) :: pos_l(3)
!
   logical :: same
!
   if (.not.Initialized) then
      call ErrorHandler('calRadialInterpGPU','module is not initialized')
   else if (.not.useGPU) then
      call ErrorHandler('calRadialInterpGPU','the GPU path is not available')
   endif
!
   nk_total = k_last - k_first + 1
   if (nk_total < 1) then
      return
   endif
   if (na < 1 .or. na > na_max) then
      call ErrorHandler('calRadialInterpGPU','na out of range',na,na_max)
   else if (ld_rint < nnr_max) then
      call ErrorHandler('calRadialInterpGPU',                          &
                        'ld_rint < nnr_max',ld_rint,nnr_max)
   else if (kpow /= -2 .and. kpow /= 0) then
      call ErrorHandler('calRadialInterpGPU','unsupported kpow',kpow)
   endif
   do ia = 1, na
      if (nnr_a(ia) < 1 .or. nnr_a(ia) > nnr_max) then
         call ErrorHandler('calRadialInterpGPU','nnr out of range',     &
                           nnr_a(ia),nnr_max)
      else if (jmax_a(ia) < 1 .or. jmax_a(ia) > jmax_max) then
         call ErrorHandler('calRadialInterpGPU','jmax out of range',    &
                           jmax_a(ia),jmax_max)
      else if (lmax_a(ia) < 0 .or. lmax_a(ia) > lmax_max) then
         call ErrorHandler('calRadialInterpGPU','lmax out of range',    &
                           lmax_a(ia),lmax_max)
      else if (nnr_a(ia)*jmax_a(ia) > ld_w) then
         call ErrorHandler('calRadialInterpGPU','w_out is too small',   &
                           nnr_a(ia)*jmax_a(ia),ld_w)
      endif
   enddo
!
   t0 = getTime()
!
   lmax_l = 0
   do ia = 1, na
      lmax_l = max(lmax_l,lmax_a(ia))
   enddo
   kmax_l = (lmax_l+1)*(lmax_l+1)
!
!  ===================================================================
!  Group the atoms by (lmax, nnr, r_interp).  r_interp is generated
!  from rmesh(1), rmesh(n_rmesh), the inscribed sphere radius and d_ir,
!  so it is a per-SPECIES quantity: grouping lets the Bessel kernel run
!  once per (k-tile, species) instead of once per (k-tile, atom).  The
!  test is exact equality, so a failure to match only costs speed.
!  ===================================================================
   ngroup = 0
   do ia = 1, na
      a_group(ia) = 0
      do ib = 1, ia-1
         if (a_group(ib) /= ib) then
            cycle                        ! ib is not a group representative
         else if (lmax_a(ib) /= lmax_a(ia) .or. nnr_a(ib) /= nnr_a(ia)) then
            cycle
         endif
         same = .true.
         do i = 1, nnr_a(ia)
            if (r_interp_a(i,ib) /= r_interp_a(i,ia)) then
               same = .false.
               exit
            endif
         enddo
         if (same) then
            a_group(ia) = ib
            exit
         endif
      enddo
      if (a_group(ia) == 0) then
         a_group(ia) = ia
         ngroup = ngroup + 1
      endif
   enddo
!
!  Processing order: group members adjacent, representative FIRST.  The
!  device invalidates the Bessel table on every k-tile and checks the
!  owning atom on reuse, so this ordering is not merely an optimisation
!  -- getting it wrong is a hard error rather than a wrong answer.
   io = 0
   do ia = 1, na
      if (a_group(ia) /= ia) then
         cycle
      endif
      io = io + 1
      a_order(io) = ia
      do ib = ia+1, na
         if (a_group(ib) == ia) then
            io = io + 1
            a_order(io) = ib
         endif
      enddo
   enddo
!
!  ===================================================================
!  Hand the per-atom geometry to the device and zero the accumulators.
!  ===================================================================
   na_l = na; ld_l = ld_rint
!  -------------------------------------------------------------------
   call ppbp_begin_atoms_gpu(na_l, lmax_a, nnr_a, jmax_a, r_interp_a, ld_l)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  k-tile loop.  Each tile is built once on the host -- including the
!  conjugated spherical harmonics, which are atom independent -- then
!  reused by every atom.
!  ===================================================================
   k0 = k_first
   do while (k0 <= k_last)
      nk = min(KT, k_last-k0+1)
!
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
      do i = 1, nk
         kvec_h(1:3,i) = getGridPointCoord('K',k0+i-1)
         fftc_h(i)     = p_fft_c(k0+i-1)
!        -------------------------------------------------------------
         call calYlmConjg( kvec_h(1:3,i), lmax_l, ylm_1 )
!        -------------------------------------------------------------
         do kl = 1, kmax_l
            ylm_h(i,kl) = ylm_1(kl)
         enddo
      enddo
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
      nk_l = nk; km_l = kmax_l; kt_l = KT
!     ----------------------------------------------------------------
      call ppbp_push_tile_gpu(kvec_h, fftc_h, ylm_h, nk_l, km_l, kt_l)
!     ----------------------------------------------------------------
!
      do io = 1, na
         ia    = a_order(io)
         ia_l  = ia
         ib_l  = a_group(ia)
         kp_l  = kpow
         nk_l  = nk
         pos_l(1:3) = posi_a(1:3,ia)
!        -------------------------------------------------------------
!        The shim synchronizes its private stream before returning, so
!        this range reports device time rather than launch time.
!        -------------------------------------------------------------
         call ppbp_process_atom_gpu(ia_l, pos_l, nk_l, kp_l, ib_l)
!        -------------------------------------------------------------
      enddo
!
      k0 = k0 + nk
   enddo
!
!  ===================================================================
!  i2l(l)/r^l scaling and copy back.
!  ===================================================================
!  -------------------------------------------------------------------
!  -------------------------------------------------------------------
   dmy_l = dummy
   do ia = 1, na
      ia_l = ia
!     ----------------------------------------------------------------
      call ppbp_finalize_atom_gpu(ia_l, dmy_l, w_out(1,ia))
!     ----------------------------------------------------------------
   enddo
!  -------------------------------------------------------------------
!  -------------------------------------------------------------------
!
   t_backproj = t_backproj + getTime() - t0
   num_calls  = num_calls + 1
!
   if (print_level >= 1) then
      write(6,'(a,i5,a,i5,a,i8,a,i5)')                                  &
         'calRadialInterpGPU:: atoms = ',na,', species groups = ',ngroup,&
         ', k-points = ',nk_total,', tiles = ',(nk_total+KT-1)/KT
   endif
!
   end subroutine calRadialInterpGPU
!  ===================================================================
!
end module PseudoPotBackProjModule
