!  *******************************************************************
!  *                                                                 *
!  *  SSMarchModule                                                  *
!  *                                                                 *
!  *  GPU offload of the single-site radial march                    *
!  *     SSSolverModule::solveSCr                                    *
!  *  the Adams-Bashforth / Adams-Moulton integration at the heart of *
!  *  the full-potential single-site solver.                          *
!  *                                                                 *
!  *  Motivation                                                     *
!  *  ----------                                                     *
!  *  On the Fe3Ni full-potential reference case (8 atoms, 1 MPI      *
!  *  rank, LIZ 135, lmax_kkr = lmax_phi = 4), SCF iteration 2 spends *
!  *                                                                 *
!  *     calValenceStates                       2,636 s              *
!  *       solveSingleScattering, 6,234 calls   1,696 s   64 %        *
!  *         -> calPhiLr -> solveSCr            <== this module       *
!  *                                                                 *
!  *  at about 3.5 GFLOP/s, roughly 7 % of ONE Grace core, on a       *
!  *  72-core socket with the GH200 idle.  The KKR matrix inversion   *
!  *  one would expect to dominate an LSMS code is already on the     *
!  *  device and sits inside the 2.6 % remainder.                     *
!  *                                                                 *
!  *  What is offloaded, and what is not                             *
!  *  ----------------------------------                             *
!  *  The march is strictly sequential in the radial index and in the *
!  *  five corrector iterations, and nothing here tries to change     *
!  *  that.  What is parallel is the angular-momentum index: one      *
!  *  march couples kmax_phi components through a single              *
!  *  kmax_phi x kmax_phi matvec per corrector.  The kernel therefore *
!  *  runs one block per march with one thread per component, and     *
!  *  because thread klp owns output klp and reads the CONTIGUOUS     *
!  *  column klp of V, the matvec needs no cross-thread reduction.    *
!  *                                                                 *
!  *  This is deliberately a drop-in replacement for one leaf         *
!  *  routine.  The surrounding control flow of calPhiLr and          *
!  *  solveSingleScattering is untouched, which is what makes the     *
!  *  change testable one march at a time against the CPU path.       *
!  *                                                                 *
!  *  It rests on the hoist done by SSSolverModule::buildVTables:     *
!  *  without a resident V table the kernel would have to rebuild the *
!  *  Gaunt contraction of the potential at every radial step, which  *
!  *  is the very thing that made the CPU loop slow.                  *
!  *                                                                 *
!  *  Transfer is NOT the ceiling -- an earlier version of this       *
!  *  comment said it was, from arithmetic, and was wrong             *
!  *  -----------------------------------------------------------    *
!  *  Nsight on the Si case reports 60 GB device-to-host in 0.71 s    *
!  *  and 10.6 GB host-to-device in 0.60 s across a WHOLE three-      *
!  *  iteration run, against 822 s of march kernel.  The earlier      *
!  *  estimate of ~320 s per iteration was wrong twice over: it used  *
!  *  the Fe3Ni site and spin counts rather than this case's, and it  *
!  *  assumed a C2C rate far below what these copies achieve.         *
!  *                                                                 *
!  *  The cost is the kernel, and in the first version it was the     *
!  *  memory ACCESS PATTERN inside it, not the arithmetic: V was read *
!  *  so that each thread walked its own contiguous column, which     *
!  *  across a warp is one cache line per lane per load.  V is now    *
!  *  transposed on upload so the warp reads consecutive addresses.   *
!  *  Do not re-derive a transfer ceiling from bandwidth estimates.   *
!  *                                                                 *
!  *  V, by contrast, is cached on the device against a               *
!  *  (site, atom, spin, SCF iteration) key.  It is 10.4 MB, so       *
!  *  pushing it per march would add another ~129 GB per iteration    *
!  *  and there would be no win at all; keyed, it is 16 uploads per   *
!  *  iteration on the real-axis path, which is 81 % of the work.     *
!  *                                                                 *
!  *  Gating                                                         *
!  *  ------                                                         *
!  *     compile time : #ifdef ACCEL                                 *
!  *     run time     : query_ssmarch_gpu -- a device exists, can be  *
!  *                    selected, and has room for this problem;      *
!  *                    MUST_SSMARCH_GPU=0 forces the CPU path        *
!  *     per march    : a minimum step count, because a short march   *
!  *                    is not worth a launch.  MUST_SSMARCH_MINSTEPS *
!  *                    overrides the default.                        *
!  *                                                                 *
!  *  A rank that fails any gate uses the CPU path rather than        *
!  *  aborting, so an ACCEL binary stays safe on a CPU-only node.     *
!  *  The CPU march in solveSCr is preserved verbatim and stays       *
!  *  reachable at run time permanently: it is the correctness        *
!  *  reference, it is what a CPU-architecture build runs, and it is  *
!  *  what keeps the measurements above reproducible.                 *
!  *                                                                 *
!  *  version 1.0, Sep 9, 2026                                       *
!  *                                                                 *
!  *******************************************************************
module SSMarchModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : ZERO
   use TimerModule, only : getTime
   use CmdLineOptionModule, only : getCmdLineOption
!
   implicit none
!
public :: initSSMarch,           &
          endSSMarch,            &
          isSSMarchGPU,          &
          pushSSMarchVTable,     &
          pushSSMarchEnergy,     &
          marchSCrGPU,           &
          addSSMarchCpuTime,     &
          getSSMarchTime,        &
          getSSMarchCount,       &
          printSSMarchInfo
!
private
!
!  ===================================================================
!  A march shorter than this many radial steps stays on the CPU: the
!  launch plus the seed/result transfers cost more than the steps save.
!  The outer march of the reference case is only ~35 steps, the inner
!  one ~997, so this cleanly separates the two.
!  ===================================================================
   integer (kind=IntKind), parameter :: MinStepsDefault = 128
!
!  ===================================================================
!  Mirrors SSM_NRED in SSMarch_Accel.cu, for reporting only.  Note the
!  march is NOT batched across marches: the kernel is launched with a
!  grid of one block per call, so the 25 angular-momentum marches of a
!  solve go out as 25 sequential launches.  All the parallelism is
!  INSIDE one march -- kmax_phi components across threads, times this
!  many reduction lanes -- which is why it occupies a single SM.
!  ===================================================================
   integer (kind=IntKind), parameter :: SSM_NRED_INFO = 8
!
!  Device head room demanded of the runtime gate, in MB, on top of the
!  buffers this module allocates.
   integer (kind=IntKind), parameter :: SlackMB = 256
!
   logical :: Initialized = .false.
   logical :: useGPU      = .false.
!
   integer (kind=IntKind) :: nr_max    = 0
   integer (kind=IntKind) :: kmax_max  = 0
   integer (kind=IntKind) :: nlj_max   = 0
   integer (kind=IntKind) :: ldV_dev   = 0
   integer (kind=IntKind) :: nrV_in    = 0
   integer (kind=IntKind) :: nrV_out   = 0
   integer (kind=IntKind) :: min_steps = MinStepsDefault
   integer (kind=IntKind) :: my_pe     = 0
   integer (kind=IntKind) :: print_level = -1
!
!  Identity of the resident V table.  All four components matter: the
!  potential changes between SCF iterations, and site/atom/spin select
!  which potential.  A stale table would be a wrong potential, not a
!  crash, so the key is compared in full.
   integer (kind=IntKind) :: v_site = -1
   integer (kind=IntKind) :: v_atom = -1
   integer (kind=IntKind) :: v_spin = -1
   integer (kind=IntKind) :: v_iter = -1
!
   integer (kind=IntKind) :: n_march_gpu = 0
   integer (kind=IntKind) :: n_march_cpu = 0
   integer (kind=IntKind) :: n_vpush     = 0
!
   real (kind=RealKind) :: t_march     = ZERO
   real (kind=RealKind) :: t_march_cpu = ZERO
   real (kind=RealKind) :: dev_MB  = ZERO
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSSMarch(nrmax, kmaxmax, nljmax, ldv, nrin, nrout,     &
                          mype, iprint)
!  ===================================================================
!  nrmax   : maximum number of radial points (iend) over the sites
!  kmaxmax : maximum kmax_phi over the sites
!  nljmax  : second dimension of bjl/bnl, i.e. lmax_phi+2
!  ldv     : leading dimension of the V tables built by SSSolverModule
!  nrin    : radial extent of the N0 == 0 table  (V_in)
!  nrout   : radial extent of the truncated-region table (V_out)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nrmax, kmaxmax, nljmax
   integer (kind=IntKind), intent(in) :: ldv, nrin, nrout, mype
   integer (kind=IntKind), intent(in), optional :: iprint
!
   integer (kind=IntKind) :: iok
   integer (kind=IntKind) :: a1, a2, a3, a4, a5, a6, a7
!
   logical :: cpu_only
!
   integer (kind=8) :: bytes_needed
!
   character (len=32) :: envval
   integer (kind=IntKind) :: estat, ival
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
   if (nrmax < 1) then
      call ErrorHandler('initSSMarch','invalid nr_max',nrmax)
   else if (kmaxmax < 1) then
      call ErrorHandler('initSSMarch','invalid kmax_max',kmaxmax)
   else if (nljmax < 1) then
      call ErrorHandler('initSSMarch','invalid nlj_max',nljmax)
   else if (ldv < kmaxmax) then
      call ErrorHandler('initSSMarch','ldV < kmax_max',ldv,kmaxmax)
   endif
!
   nr_max   = nrmax
   kmax_max = kmaxmax
   nlj_max  = nljmax
   ldV_dev  = ldv
   nrV_in   = max(nrin,1)
   nrV_out  = max(nrout,1)
   my_pe    = mype
!
!  ===================================================================
!  Device footprint of the buffers allocated by init_ssmarch_gpu.
!  V dominates: it scales as ldV^2 * numrs, so an lmax_phi = 6 case is
!  four times an lmax_phi = 4 one.  Size from the parameters, not from
!  constants.
!  ===================================================================
   bytes_needed = 16_8*int(ldV_dev,8)*int(ldV_dev,8)*int(nrV_in,8)     &
                + 16_8*int(ldV_dev,8)*int(ldV_dev,8)*int(nrV_out,8)    &
                +  8_8*int(nr_max,8)                                   &
                + 32_8*int(nr_max,8)*int(nlj_max,8)                    &
                + 16_8*int(nr_max,8)                                   &
                +  4_8*int(kmax_max,8)                                 &
                + 32_8*int(nr_max,8)*int(kmax_max,8)                   &
                + 32_8*int(kmax_max,8)*4_8
!
   dev_MB = real(bytes_needed,kind=RealKind)/(1024.0d0*1024.0d0)
   bytes_needed = bytes_needed + int(SlackMB,8)*1024_8*1024_8
!
!  ===================================================================
!  Per-march step threshold, overridable without a rebuild.
!  ===================================================================
   min_steps = MinStepsDefault
   call get_environment_variable('MUST_SSMARCH_MINSTEPS',envval,status=estat)
   if (estat == 0) then
      read(envval,*,iostat=estat) ival
      if (estat == 0 .and. ival >= 0) then
         min_steps = ival
      endif
   endif
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
!  march, which is what makes an ACCEL binary safe to launch on a
!  CPU-only node and safe to oversubscribe on a one-GPU-per-node
!  machine such as Vista.
!  ===================================================================
   iok = 0
   if (.not.cpu_only) then
!     ----------------------------------------------------------------
      call query_ssmarch_gpu(bytes_needed, iok)
!     ----------------------------------------------------------------
   endif
   if (iok == 1) then
      a1 = nr_max;   a2 = kmax_max; a3 = nlj_max; a4 = ldV_dev
      a5 = nrV_in;   a6 = nrV_out;  a7 = my_pe
!     ----------------------------------------------------------------
      call init_ssmarch_gpu(a1, a2, a3, a4, a5, a6, a7)
!     ----------------------------------------------------------------
      useGPU = .true.
   else if (print_level >= 0) then
      if (cpu_only) then
!        -------------------------------------------------------------
         call WarningHandler('initSSMarch',                             &
                 '-cpu/--cpu-only given; the CPU march is used')
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call WarningHandler('initSSMarch',                             &
                 'GPU probe failed or was disabled; using the CPU march')
!        -------------------------------------------------------------
      endif
   endif
#endif
!
   v_site = -1; v_atom = -1; v_spin = -1; v_iter = -1
   n_march_gpu = 0; n_march_cpu = 0; n_vpush = 0
   t_march     = ZERO
   t_march_cpu = ZERO
   Initialized = .true.
!
   if (print_level >= 0) then
!     ----------------------------------------------------------------
      call printSSMarchInfo()
!     ----------------------------------------------------------------
   endif
!
   end subroutine initSSMarch
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSSMarch()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
      return
   endif
!
!  ===================================================================
!  Report what the march actually cost.  t_march is host wall time
!  around ssm_march_gpu, which synchronises its stream, so it is real
!  device plus transfer time and not launch time.  Without this the
!  only way to know the kernel's share was to run under Nsight and
!  subtract, which is how the first version's 822 s went unnoticed
!  until it was profiled.
!  ===================================================================
   if (print_level >= 0 .and. (n_march_gpu > 0 .or. n_march_cpu > 0)) then
      write(6,'(/,a)') 'Single-site radial march: cost'
      write(6,'(a,i12)')   '   marches on the GPU             = ', n_march_gpu
      write(6,'(a,i12)')   '   marches on the CPU             = ', n_march_cpu
      write(6,'(a,i12)')   '   V-table uploads                = ', n_vpush
      write(6,'(a,f12.3)') '   total GPU march time    (sec)  = ', t_march
      if (n_march_gpu > 0) then
         write(6,'(a,f12.5)') '   per GPU march           (msec) = ',      &
                              1000.0d0*t_march/real(n_march_gpu,kind=RealKind)
      endif
      write(6,'(a,f12.3)') '   total CPU march time    (sec)  = ', t_march_cpu
      if (n_march_cpu > 0) then
         write(6,'(a,f12.5)') '   per CPU march           (msec) = ',      &
                          1000.0d0*t_march_cpu/real(n_march_cpu,kind=RealKind)
      endif
!     ================================================================
!     Compare the two directly.  The per-march numbers are not
!     like-for-like on their own -- the CPU count includes the short
!     outer marches that never go to the device -- so read the totals
!     alongside the counts, and use MUST_SSMARCH_GPU=1 versus unset for
!     the decisive comparison.
!     ================================================================
      write(6,'(a)') ' '
   endif
!
#ifdef ACCEL
   if (useGPU) then
!     ----------------------------------------------------------------
      call finalize_ssmarch_gpu()
!     ----------------------------------------------------------------
   endif
#endif
!
   useGPU = .false.
   Initialized = .false.
!
   end subroutine endSSMarch
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSSMarchGPU() result(y)
!  ===================================================================
   implicit none
   logical :: y
!
   y = ( Initialized .and. useGPU )
!
   end function isSSMarchGPU
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine addSSMarchCpuTime(t)
!  ===================================================================
!  Accumulate the cost of one CPU reference march.  solveSCr calls this
!  whenever marchSCrGPU declines, so a single run reports what BOTH
!  paths cost and the question "is the offload worth it on this system"
!  no longer needs a two-run A/B.  That question turned out to matter:
!  on the Fe3Ni full-potential case the CPU march is 2.45x faster than
!  the device one, which is why the device path is now opt-in.
!  ===================================================================
   implicit none
   real (kind=RealKind), intent(in) :: t
   t_march_cpu = t_march_cpu + t
   end subroutine addSSMarchCpuTime
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSSMarchTime() result(t)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: t
   t = t_march
   end function getSSMarchTime
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getSSMarchCount(ngpu,ncpu,nvpush)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(out) :: ngpu, ncpu, nvpush
   ngpu = n_march_gpu; ncpu = n_march_cpu; nvpush = n_vpush
   end subroutine getSSMarchCount
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSSMarchInfo()
!  ===================================================================
   implicit none
!
   write(6,'(/,a)') 'Single-site radial march: Adams-Moulton, one march per launch'
   write(6,'(a,i12)')   '   radial points       nr_max     = ', nr_max
   write(6,'(a,i12)')   '   angular channels    kmax_phi   = ', kmax_max
   write(6,'(a,i12)')   '   V leading dim       ldV        = ', ldV_dev
   write(6,'(a,i12)')   '   V radial extent     inner      = ', nrV_in
   write(6,'(a,i12)')   '   V radial extent     outer      = ', nrV_out
   write(6,'(a,i12)')   '   min steps for GPU              = ', min_steps
   write(6,'(a,i12)')   '   marches per kernel launch      = ', 1
   write(6,'(a,i12)')   '   reduction lanes per component  = ', SSM_NRED_INFO
   write(6,'(a,f12.2)') '   device memory       MB         = ', dev_MB
   if (useGPU) then
      write(6,'(a)')    '   evaluation path                = GPU (CUDA)'
      write(6,'(a)')    '      enabled by MUST_SSMARCH_GPU; unset it to use the CPU path'
   else
      write(6,'(a)')    '   evaluation path                = CPU (reference)'
      write(6,'(a)')    '      The device march is OFF BY DEFAULT.  It was measured'
      write(6,'(a)')    '      2.45x SLOWER than this CPU path on GH200 (Fe3Ni full'
      write(6,'(a)')    '      potential: 5.16 ms vs 2.10 ms per march), because one'
      write(6,'(a)')    '      march is a long sequential chain that occupies a single'
      write(6,'(a)')    '      warp on one SM.  Set MUST_SSMARCH_GPU=1 to enable it and'
      write(6,'(a)')    '      compare -- both paths report their cost at end of run.'
      write(6,'(a)')    '      MUST_SSMARCH_MINSTEPS=<n> changes the per-march threshold.'
      write(6,'(a)')    '      -cpu/--cpu-only overrides both and disables every offload.'
   endif
   write(6,'(a)') ' '
!
   end subroutine printSSMarchInfo
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pushSSMarchVTable(site, atom, spin, iter, V_in, V_out,    &
                                ldv, nrin, nrout, haveOut, pushed)
!  ===================================================================
!  Upload the hoisted Gaunt contraction of the potential, but only when
!  the (site, atom, spin, SCF iteration) key has changed.  On the
!  real-axis path a single (site,spin) is held across hundreds of
!  energies, so this collapses thousands of 10.4 MB uploads into one.
!
!  pushed reports whether an upload actually happened, so the caller
!  can account for it.
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, spin, iter
   integer (kind=IntKind), intent(in) :: ldv, nrin, nrout
   logical, intent(in) :: haveOut
   logical, intent(out) :: pushed
!
   complex (kind=CmplxKind), intent(in) :: V_in(*)
   complex (kind=CmplxKind), intent(in) :: V_out(*)
!
   integer (kind=IntKind) :: a1, a2, a3
!
   pushed = .false.
   if (.not.isSSMarchGPU()) then
      return
   endif
!
   if (site == v_site .and. atom == v_atom .and. spin == v_spin        &
                      .and. iter == v_iter) then
      return
   endif
!
#ifdef ACCEL
   a1 = ldv; a2 = nrin; a3 = 0
!  -------------------------------------------------------------------
   call ssm_push_vtable_gpu(V_in, a1, a2, a3)
!  -------------------------------------------------------------------
   if (haveOut) then
      a2 = nrout; a3 = 1
!     ----------------------------------------------------------------
      call ssm_push_vtable_gpu(V_out, a1, a2, a3)
!     ----------------------------------------------------------------
   endif
#endif
!
   v_site = site; v_atom = atom; v_spin = spin; v_iter = iter
   n_vpush = n_vpush + 1
   pushed = .true.
!
   end subroutine pushSSMarchVTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pushSSMarchEnergy(r_mesh, bjl, bnl, cm0, lofk_a, iend,    &
                                nlj, kmax)
!  ===================================================================
!  The per-(site,energy) tables.  Small next to V: about 0.5 MB against
!  10.4 MB for the reference case, so these are pushed unconditionally
!  once per solveSingleScattering rather than keyed.
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: iend, nlj, kmax
   integer (kind=IntKind), intent(in) :: lofk_a(kmax)
!
   real (kind=RealKind), intent(in) :: r_mesh(iend)
!
   complex (kind=CmplxKind), intent(in) :: bjl(*), bnl(*), cm0(iend)
!
   integer (kind=IntKind) :: a1, a2, a3
!
   if (.not.isSSMarchGPU()) then
      return
   endif
!
   if (iend > nr_max) then
      call ErrorHandler('pushSSMarchEnergy','iend > nr_max',iend,nr_max)
   else if (nlj > nlj_max) then
      call ErrorHandler('pushSSMarchEnergy','nlj > nlj_max',nlj,nlj_max)
   else if (kmax > kmax_max) then
      call ErrorHandler('pushSSMarchEnergy','kmax > kmax_max',kmax,kmax_max)
   endif
!
#ifdef ACCEL
   a1 = iend; a2 = nlj; a3 = kmax
!  -------------------------------------------------------------------
   call ssm_push_energy_gpu(r_mesh, bjl, bnl, cm0, lofk_a, a1, a2, a3)
!  -------------------------------------------------------------------
#endif
!
   end subroutine pushSSMarchEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function marchSCrGPU(N1, N2, N0, nstep, icmax, nr, kmax, nlj,        &
                        llp1, has_l0, hfac, kappa, e2oc2, potshift,     &
                        sx, cx, dsx, dcx) result(done)
!  ===================================================================
!  Run one march on the device.  Returns .false. WITHOUT touching any
!  argument if the device path is unavailable or the march is too short
!  to be worth a launch, in which case the caller runs the CPU march.
!
!  sx and cx are the caller's sx(nr,kmax) / cx(nr,kmax); only the seed
!  row N1-nstep is uploaded and only the rows this march writes are
!  brought back.  dsx/dcx are the (kmax,4) derivative history and are
!  updated in place, because the next march of the same calPhiLr call
!  continues from them.
!  ===================================================================
   implicit none
!
   logical :: done
!
   integer (kind=IntKind), intent(in) :: N1, N2, N0, nstep, icmax
   integer (kind=IntKind), intent(in) :: nr, kmax, nlj, llp1
   logical, intent(in) :: has_l0
!
   real (kind=RealKind), intent(in) :: hfac
!
   complex (kind=CmplxKind), intent(in) :: kappa, e2oc2, potshift
   complex (kind=CmplxKind), intent(inout) :: sx(nr,kmax), cx(nr,kmax)
   complex (kind=CmplxKind), intent(inout) :: dsx(kmax,4), dcx(kmax,4)
!
   integer (kind=IntKind) :: nsteps, ihas, irm1
   integer (kind=IntKind) :: b1,b2,b3,b4,b5,b6,b7,b8,b9
   real (kind=RealKind) :: t0, hl
!
   done = .false.
   if (.not.isSSMarchGPU()) then
      n_march_cpu = n_march_cpu + 1
      return
   endif
!
   nsteps = abs(N2-N1) + 1
   if (nsteps < min_steps) then
      n_march_cpu = n_march_cpu + 1
      return
   endif
!
   irm1 = N1 - nstep
   if (irm1 < 1 .or. irm1 > nr) then
!     ================================================================
!     The seed row lies outside the array; this march cannot be run on
!     the device.  Fall back rather than reading out of bounds.
!     ================================================================
      n_march_cpu = n_march_cpu + 1
      return
   endif
!
   t0 = getTime()
!
   if (has_l0) then
      ihas = 1
   else
      ihas = 0
   endif
   b1=N1; b2=N2; b3=N0; b4=nstep; b5=icmax
   b6=nr; b7=kmax; b8=nlj; b9=llp1
   hl = hfac
!
#ifdef ACCEL
!  -------------------------------------------------------------------
   call ssm_march_gpu(b1,b2,b3,b4,b5,b6,b7,b8,b9,ihas,hl,               &
                      kappa,e2oc2,potshift,sx,cx,dsx,dcx)
!  -------------------------------------------------------------------
#endif
!
   t_march = t_march + getTime() - t0
   n_march_gpu = n_march_gpu + 1
   done = .true.
!
   end function marchSCrGPU
!  ===================================================================
!
end module SSMarchModule
