!  *******************************************************************
!  *                                                                 *
!  *  DensityOnGridModule                                            *
!  *                                                                 *
!  *  Batched evaluation of an L-expanded density on the structured   *
!  *  grid  (radial mesh node) x (angular direction), i.e. on the     *
!  *  direct product of the radial mesh of an atom and the fixed      *
!  *  spherical grid held by AngularIntegrationModule.                *
!  *                                                                 *
!  *  Motivation                                                     *
!  *  ----------                                                     *
!  *  PotentialGenerationModule::calExchangeJl used to evaluate the   *
!  *  charge (and moment) density one point at a time through         *
!  *  ChargeDensityModule::getChargeDensityAtPoint.  Each such call   *
!  *  performs a bisection search (hunt) over the radial mesh plus an *
!  *  n_inter-point Neville interpolation for every jl component,     *
!  *  even though                                                    *
!  *     (a) the evaluation points are r_i*u_g with r_i EXACTLY on    *
!  *         radial mesh nodes, and                                  *
!  *     (b) the angular directions u_g are fixed for the whole run.  *
!  *  With n_r ~ 1500 radial points and n_g = 50*80 = 4000 angular    *
!  *  directions that is ~6x10^6 interpolations per (atom, species)   *
!  *  per SCF iteration -- the dominant full-potential cost.          *
!  *                                                                 *
!  *  Method                                                         *
!  *  ------                                                         *
!  *  The radial interpolation is a LINEAR operator on the stored     *
!  *  coefficients, so it can be separated from the angular sum:      *
!  *                                                                 *
!  *     A(i,jl)      = sum_k w_k(r_i)*den_l(irp(i)+k-1, jl)         *
!  *     rho(r_i,u_g) = sum_jl fa2(jl)*Re[ A(i,jl)*Y_jl(u_g) ]       *
!  *                  = [Ar*Wr - Ai*Wi](i,g)                         *
!  *                                                                 *
!  *  with Wr(jl,g) = fa2(jl)*Re[Y_jl(u_g)] and Wi likewise.  The     *
!  *  angular sum is therefore two real GEMMs, and the weight tables  *
!  *  Wr, Wi are built once.  w_k are the Lagrange weights of the     *
!  *  same n_inter-point polynomial that polint_inline evaluates by   *
!  *  the Neville recursion, so the result is numerically identical   *
!  *  and is bitwise exact for on-node targets (one numerator factor  *
!  *  is then exactly zero and the surviving ratio is exactly one).   *
!  *                                                                 *
!  *  Two execution paths are provided and are numerically equivalent:*
!  *     - CPU: BLAS DGEMM (always available)                         *
!  *     - GPU: Accelerator/DensityInterp_Accel.cu, under #ifdef ACCEL*
!  *                                                                 *
!  *  All array dummies that are handed to BLAS or to the CUDA layer  *
!  *  are declared explicit-shape / assumed-size so that they are     *
!  *  sequence associated (no descriptors, no copy-in/copy-out).      *
!  *                                                                 *
!  *  version 1.0, Aug 24, 2026                                      *
!  *                                                                 *
!  *******************************************************************
module DensityOnGridModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : ZERO, ONE, TWO
   use IntegerFactorsModule, only : mofj, kofj
   use TimerModule, only : getTime
   use CmdLineOptionModule, only : getCmdLineOption
!
   implicit none
!
public :: initDensityOnGrid,       &
          endDensityOnGrid,        &
          calDensityOnAngularGrid, &
          calDensityGradOnAngularGrid, &
          isDensityOnGridGPU,      &
          isDensityGradAvailable,  &
          getDensityOnGridTime,    &
          printDensityOnGridInfo,  &
          FIELD_CHARGE,            &
          FIELD_MOMENT,            &
          FIELD_DER_CHARGE,        &
          FIELD_DER_MOMENT
!
private
!
!  ===================================================================
!  Field identifiers selecting the device-side buffer.  They must match
!  the ifield convention of DensityInterp_Accel.cu (1-based).
!  ===================================================================
   integer (kind=IntKind), parameter :: FIELD_CHARGE     = 1
   integer (kind=IntKind), parameter :: FIELD_MOMENT     = 2
   integer (kind=IntKind), parameter :: FIELD_DER_CHARGE = 3
   integer (kind=IntKind), parameter :: FIELD_DER_MOMENT = 4
   integer (kind=IntKind), parameter :: NumFields        = 4
!
!  ===================================================================
!  n_inter must match ChargeDensityModule's interpolation order so that
!  this module reproduces getChargeDensityAtPoint exactly.
!  ===================================================================
   integer (kind=IntKind), parameter :: n_inter = 5
!
   logical :: Initialized = .false.
   logical :: useGPU = .false.
   logical :: gradReady = .false.
!
   integer (kind=IntKind) :: nr_max   = 0
   integer (kind=IntKind) :: jmax_max = 0
   integer (kind=IntKind) :: kmax_max = 0
   integer (kind=IntKind) :: num_ang  = 0
   integer (kind=IntKind) :: print_level = -1
!
!  Fixed angular weight tables:  W(jl,g) = fa2(jl)*Y_jl(u_g)
   real (kind=RealKind), allocatable :: Wr(:,:), Wi(:,:)
   real (kind=RealKind), allocatable :: fa2(:)
   integer (kind=IntKind), allocatable :: kofj_tab(:)
!
!  Fixed angular GRADIENT tables (GGA only), allocated on first use:
!     G(jl,g,c) = fa2(jl) * [ r * d/dx_c Y_jl(u_g) ]  evaluated at r = 1
!  grad_ylm carries an overall 1/r (rfac = clm(jl)/r in
!  SphericalHarmonics4), so the unit-sphere value is the pure angular
!  part and the 1/r_i scaling is reapplied per radius.
   real (kind=RealKind), allocatable :: Gr(:,:,:), Gi(:,:,:)
   real (kind=RealKind), allocatable :: upos_tab(:,:)
!
!  Work space
   real (kind=RealKind), allocatable :: Ar(:,:), Ai(:,:)
   real (kind=RealKind), allocatable :: Dr(:,:), Tc(:,:)
   real (kind=RealKind), allocatable :: wlag(:,:)
   integer (kind=IntKind), allocatable :: irp_tab(:)
!
!  Identifier of the radial mesh currently resident on the device
   integer (kind=IntKind) :: mesh_id_on_gpu = -1
!
   real (kind=RealKind) :: t_interp = ZERO
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initDensityOnGrid(nrmax, jmaxmax, lmax, mype, iprint, needGrad)
!  ===================================================================
!  nrmax   : maximum number of radial mesh points over the local atoms
!  jmaxmax : maximum number of (l,m>=0) density components
!  lmax    : lmax of the Ylm table held by AngularIntegrationModule
!  mype    : MPI rank, used for the GPU device assignment (rank modulo
!            the number of visible devices, as in init_lsms_gpu)
!  ===================================================================
   use AngularIntegrationModule, only : getNumSphericalGridPoints,     &
                                        getYlmTable, getUnitVecTable
!
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nrmax, jmaxmax, lmax, mype
   integer (kind=IntKind), intent(in), optional :: iprint
   logical, intent(in), optional :: needGrad
!
   integer (kind=IntKind) :: jl, kl, g, c
   integer (kind=IntKind) :: nf, ni, ng_l, km_l, nr_l, jm_l, mp_l
   integer (kind=IntKind) :: iok
!
   integer (kind=8) :: bytes_needed
!
   logical :: want_grad
   logical :: cpu_only
!
   real (kind=RealKind), pointer :: p_upos(:,:)
!
   complex (kind=CmplxKind), pointer :: p_ylm(:,:)
   complex (kind=CmplxKind), allocatable :: ylm_u(:), grady_u(:,:)
   complex (kind=CmplxKind), allocatable :: grady_flat(:,:,:)
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
   nr_max   = nrmax
   jmax_max = jmaxmax
   kmax_max = (lmax+1)*(lmax+1)
   num_ang  = getNumSphericalGridPoints()
!
   if (nr_max < n_inter) then
      call ErrorHandler('initDensityOnGrid','nr_max < n_inter',nr_max,n_inter)
   else if (jmax_max < 1) then
      call ErrorHandler('initDensityOnGrid','invalid jmax_max',jmax_max)
   else if (num_ang < 1) then
      call ErrorHandler('initDensityOnGrid',                             &
                        'invalid number of angular points',num_ang)
   endif
!
   if (present(needGrad)) then
      want_grad = needGrad
   else
      want_grad = .false.
   endif
!
   allocate( fa2(jmax_max), kofj_tab(jmax_max) )
   allocate( Wr(jmax_max,num_ang), Wi(jmax_max,num_ang) )
   allocate( Ar(nr_max,jmax_max), Ai(nr_max,jmax_max) )
   allocate( wlag(nr_max,n_inter), irp_tab(nr_max) )
!
!  ===================================================================
!  fa2 folds the m < 0 components onto m > 0: the stored expansion runs
!  over jl = (l,|m|) only, and for a real function
!     den_(l,-m) = (-1)^m conjg(den_(l,m)),
!  so each m /= 0 pair contributes 2*Re[den_(l,m) Y_(l,m)].  This is the
!  same fa2 array that getChargeDensityAtPoint builds.
!  ===================================================================
   do jl = 1, jmax_max
      if ( mofj(jl) == 0 ) then
         fa2(jl) = ONE
      else
         fa2(jl) = TWO
      endif
      kofj_tab(jl) = kofj(jl)
   enddo
!
!  ===================================================================
!  Pack the fixed angular weight tables from the Ylm table that
!  AngularIntegrationModule already built on its spherical grid.
!  ===================================================================
   p_ylm => getYlmTable()
   if ( size(p_ylm,1) /= num_ang ) then
      call ErrorHandler('initDensityOnGrid',                             &
                        'Ylm table angular size mismatch',                &
                        size(p_ylm,1),num_ang)
   else if ( size(p_ylm,2) < kmax_max ) then
      call ErrorHandler('initDensityOnGrid','Ylm table kmax too small',   &
                        size(p_ylm,2),kmax_max)
   endif
!
   do g = 1, num_ang
      do jl = 1, jmax_max
         kl = kofj_tab(jl)
         if ( kl >= 1 .and. kl <= kmax_max ) then
            Wr(jl,g) = fa2(jl)*real(p_ylm(g,kl),kind=RealKind)
            Wi(jl,g) = fa2(jl)*aimag(p_ylm(g,kl))
         else
            Wr(jl,g) = ZERO
            Wi(jl,g) = ZERO
         endif
      enddo
   enddo
!
!  ===================================================================
!  Fixed angular GRADIENT tables (GGA only).
!
!  calYlm is evaluated on the UNIT sphere, where the 1/r factor carried
!  by grad_ylm equals one, so grady_u is the pure angular part.  The
!  per-radius 1/r_i scaling is applied at evaluation time.  This mirrors
!  the CPU line in getChargeDensityAtPoint,
!     grad(i) += fa2(jl)*Re[ der_rho_in*ylm(kl)*er(i)
!                          + rho_in(jl)*grad_ylm(kl,i) ]
!  with er = u_g.
!  ===================================================================
   gradReady = .false.
   if ( want_grad ) then
      allocate( Gr(jmax_max,num_ang,3), Gi(jmax_max,num_ang,3) )
      allocate( upos_tab(3,num_ang) )
      allocate( Dr(nr_max,num_ang), Tc(nr_max,num_ang) )
      allocate( ylm_u(kmax_max), grady_u(kmax_max,3) )
      allocate( grady_flat(num_ang,kmax_max,3) )
!
      p_upos => getUnitVecTable()
      do g = 1, num_ang
         upos_tab(1:3,g) = p_upos(1:3,g)
!        -------------------------------------------------------------
         call calYlm(p_upos(1:3,g), lmax, ylm_u, grady_u)
!        -------------------------------------------------------------
         do c = 1, 3
            do kl = 1, kmax_max
               grady_flat(g,kl,c) = grady_u(kl,c)
            enddo
            do jl = 1, jmax_max
               kl = kofj_tab(jl)
               if ( kl >= 1 .and. kl <= kmax_max ) then
                  Gr(jl,g,c) = fa2(jl)*real(grady_u(kl,c),kind=RealKind)
                  Gi(jl,g,c) = fa2(jl)*aimag(grady_u(kl,c))
               else
                  Gr(jl,g,c) = ZERO
                  Gi(jl,g,c) = ZERO
               endif
            enddo
         enddo
      enddo
      nullify( p_upos )
      gradReady = .true.
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
!  Until now this module had NONE of the three: it set useGPU = .true.
!  from #ifdef ACCEL alone, with no flag, no environment override and no
!  device query, so it ignored -cpu/--cpu-only and died at its first
!  cudaMalloc on a node without a usable GPU instead of falling back.
!  The CPU path below (host DGEMM) was already there and reachable; only
!  the decision to take it was missing.
!  ===================================================================
   iok = 0
   if (.not.cpu_only) then
      bytes_needed = 16_8*int(nr_max,8)*int(jmax_max,8)                 &
                   + 16_8*int(nr_max,8)*int(num_ang,8)*int(NumFields,8) &
                   + 16_8*int(num_ang,8)*int(kmax_max,8)                &
                   +  8_8*int(nr_max,8)*int(n_inter,8)
!     ----------------------------------------------------------------
      call query_density_interp_gpu(bytes_needed, iok)
!     ----------------------------------------------------------------
   endif
   if (iok == 1) then
      nr_l = nr_max; jm_l = jmax_max; ng_l = num_ang
      km_l = kmax_max; ni = n_inter; nf = NumFields; mp_l = mype
!     ----------------------------------------------------------------
      call init_density_interp_gpu(nr_l, jm_l, ng_l, km_l, ni, nf, mp_l)
!     ----------------------------------------------------------------
      call push_angular_ylm_gpu(p_ylm, ng_l, km_l, kofj_tab, fa2, jm_l)
!     ----------------------------------------------------------------
      if ( gradReady ) then
!        -------------------------------------------------------------
         call push_angular_gradylm_gpu(grady_flat, upos_tab, ng_l, km_l,&
                                       kofj_tab, fa2, jm_l)
!        -------------------------------------------------------------
      endif
      useGPU = .true.
   else if (print_level >= 0) then
      if (cpu_only) then
!        -------------------------------------------------------------
         call WarningHandler('initDensityOnGrid',                       &
                 '-cpu/--cpu-only given; the host DGEMM path is used')
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call WarningHandler('initDensityOnGrid',                       &
                 'GPU probe failed or was disabled; using host DGEMM')
!        -------------------------------------------------------------
      endif
   endif
#endif
!
   if ( want_grad ) then
      deallocate( ylm_u, grady_u, grady_flat )
   endif
   nullify( p_ylm )
   mesh_id_on_gpu = -1
   t_interp = ZERO
   Initialized = .true.
!
   if (print_level >= 0) then
!     ----------------------------------------------------------------
      call printDensityOnGridInfo()
!     ----------------------------------------------------------------
   endif
!
   end subroutine initDensityOnGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endDensityOnGrid()
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
      call finalize_density_interp_gpu()
!     ----------------------------------------------------------------
   endif
#endif
!
   deallocate( fa2, kofj_tab, Wr, Wi, Ar, Ai, wlag, irp_tab )
   if ( gradReady ) then
      deallocate( Gr, Gi, upos_tab, Dr, Tc )
      gradReady = .false.
   endif
!
   nr_max = 0; jmax_max = 0; kmax_max = 0; num_ang = 0
   mesh_id_on_gpu = -1
   useGPU = .false.
   Initialized = .false.
!
   end subroutine endDensityOnGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isDensityOnGridGPU() result(y)
!  ===================================================================
   implicit none
   logical :: y
   y = useGPU
   end function isDensityOnGridGPU
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isDensityGradAvailable() result(y)
!  ===================================================================
!  True when the module was initialised with needGrad = .true., i.e.
!  when the angular gradient tables exist and
!  calDensityGradOnAngularGrid may be called.
!  ===================================================================
   implicit none
   logical :: y
   y = gradReady
   end function isDensityGradAvailable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDensityOnGridTime() result(t)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: t
   t = t_interp
   end function getDensityOnGridTime
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDensityOnGridInfo()
!  ===================================================================
   implicit none
!
   write(6,'(/,80(''=''))')
   write(6,'(a)')    'DensityOnGridModule: batched density evaluation'
   write(6,'(a,i8)') '   max radial points   nr_max   = ',nr_max
   write(6,'(a,i8)') '   max L components    jmax_max = ',jmax_max
   write(6,'(a,i8)') '   angular directions  num_ang  = ',num_ang
   write(6,'(a,i8)') '   interpolation order n_inter  = ',n_inter
   if (useGPU) then
      write(6,'(a)') '   evaluation path              = GPU (CUDA + cuBLAS)'
   else
      write(6,'(a)') '   evaluation path              = CPU (BLAS DGEMM)'
   endif
   write(6,'(80(''=''),/)')
!
   end subroutine printDensityOnGridInfo
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDensityOnAngularGrid(den_l, ld_den, r_mesh, nr, jmax,  &
                                      ifield, den_grid, ld_grid, mesh_id)
!  ===================================================================
!  Evaluates
!      den_grid(i,g) = sum_jl fa2(jl)*Re[ A(i,jl)*Y_jl(u_g) ]
!  for i = 1..nr (the targets are the radial mesh nodes themselves) and
!  g = 1..num_ang.
!
!  den_l    : (ld_den,jmax) complex L-expansion coefficients, e.g. the
!             array returned by getChargeDensity("TotalNew",id,ia)
!  ld_den   : leading dimension of den_l (= size(den_l,1))
!  r_mesh   : radial mesh, at least nr points
!  nr       : number of radial nodes to evaluate on
!  jmax     : number of L components to include
!  ifield   : FIELD_CHARGE or FIELD_MOMENT (selects the device buffer)
!  den_grid : (ld_grid,num_ang) real output
!  ld_grid  : leading dimension of den_grid
!  mesh_id  : optional identifier of the radial mesh; when unchanged
!             from the previous call the device-side mesh copy is reused
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ld_den, nr, jmax, ifield, ld_grid
   integer (kind=IntKind), intent(in), optional :: mesh_id
!
   real (kind=RealKind), intent(in) :: r_mesh(*)
   real (kind=RealKind), intent(out) :: den_grid(ld_grid,*)
!
   complex (kind=CmplxKind), intent(in) :: den_l(ld_den,*)
!
   integer (kind=IntKind) :: i, jl, k, ir, irp, half, mid, ip
   integer (kind=IntKind) :: nr_l, jm_l, fld, ldd
!
   real (kind=RealKind) :: x, num, den, t0
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('calDensityOnAngularGrid','module not initialized')
   else if (nr < n_inter .or. nr > nr_max) then
      call ErrorHandler('calDensityOnAngularGrid','nr out of range',nr,nr_max)
   else if (jmax < 1 .or. jmax > jmax_max) then
      call ErrorHandler('calDensityOnAngularGrid','jmax out of range',      &
                        jmax,jmax_max)
   else if (ifield < 1 .or. ifield > NumFields) then
      call ErrorHandler('calDensityOnAngularGrid','invalid ifield',ifield)
   else if (ld_grid < nr) then
      call ErrorHandler('calDensityOnAngularGrid','ld_grid < nr',ld_grid,nr)
   else if (ld_den < nr) then
      call ErrorHandler('calDensityOnAngularGrid','ld_den < nr',ld_den,nr)
   endif
!
   t0 = getTime()
!
   nr_l = nr; jm_l = jmax; fld = ifield; ldd = ld_den
!
#ifdef ACCEL
   if (useGPU) then
      if ( .not.present(mesh_id) ) then
!        -------------------------------------------------------------
         call push_radial_mesh_gpu(r_mesh, nr_l)
!        -------------------------------------------------------------
      else if ( mesh_id /= mesh_id_on_gpu ) then
!        -------------------------------------------------------------
         call push_radial_mesh_gpu(r_mesh, nr_l)
!        -------------------------------------------------------------
         mesh_id_on_gpu = mesh_id
      endif
!     ----------------------------------------------------------------
      call push_density_l_gpu(den_l, ldd, nr_l, jm_l, fld)
!     ----------------------------------------------------------------
!     The evaluation targets are the mesh nodes, so the device reuses
!     the radial mesh it already holds as the target list.
!     ----------------------------------------------------------------
      call eval_density_sphere_gpu(nr_l, jm_l, fld, den_grid, ld_grid)
!     ----------------------------------------------------------------
      t_interp = t_interp + (getTime()-t0)
      return
   endif
#endif
!
!  ===================================================================
!  CPU path.  Identical algebra to the GPU path.
!
!  Step 1: interpolation stencil offset and Lagrange weights, using the
!          same clamping logic as getChargeDensityAtPoint.
!  ===================================================================
   half = (n_inter-1)/2
   ir = 1
   do i = 1, nr
      x = r_mesh(i)
!     ----------------------------------------------------------------
      call hunt(nr, r_mesh(1:nr), x, ir)
!     ----------------------------------------------------------------
      if ( ir > nr-half ) then
         irp = nr-n_inter+1
      else if ( 2*ir+1 > n_inter ) then
         irp = ir-half
      else
         irp = 1
      endif
      irp = max(1, min(irp, nr-n_inter+1))
      irp_tab(i) = irp
!
      do k = 1, n_inter
         num = ONE
         den = ONE
         do mid = 1, n_inter
            if ( mid == k ) then
               cycle
            endif
            num = num*(x - r_mesh(irp+mid-1))
            den = den*(r_mesh(irp+k-1) - r_mesh(irp+mid-1))
         enddo
         wlag(i,k) = num/den
      enddo
   enddo
!
!  ===================================================================
!  Step 2: radial interpolation of every L component, split into real
!          and imaginary parts.
!  ===================================================================
   do jl = 1, jmax
      do i = 1, nr
         Ar(i,jl) = ZERO
         Ai(i,jl) = ZERO
      enddo
      do k = 1, n_inter
         do i = 1, nr
            ip = irp_tab(i) + k - 1
            Ar(i,jl) = Ar(i,jl) + wlag(i,k)*real(den_l(ip,jl),kind=RealKind)
            Ai(i,jl) = Ai(i,jl) + wlag(i,k)*aimag(den_l(ip,jl))
         enddo
      enddo
   enddo
!
!  ===================================================================
!  Step 3: angular sum as two real GEMMs,
!             den_grid = Ar*Wr - Ai*Wi
!  ===================================================================
!  -------------------------------------------------------------------
   call dgemm('n','n', nr, num_ang, jmax,  ONE, Ar, nr_max,             &
              Wr, jmax_max, ZERO, den_grid, ld_grid)
!  -------------------------------------------------------------------
   call dgemm('n','n', nr, num_ang, jmax, -ONE, Ai, nr_max,             &
              Wi, jmax_max,  ONE, den_grid, ld_grid)
!  -------------------------------------------------------------------
!
   t_interp = t_interp + (getTime()-t0)
!
   end subroutine calDensityOnAngularGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDensityGradOnAngularGrid(den_l, der_den_l, ld_den,     &
                                          r_mesh, nr, jmax, ifield_val,  &
                                          ifield_der, den_grid,          &
                                          grad_grid, ld_grid, mesh_id)
!  ===================================================================
!  GGA companion of calDensityOnAngularGrid: returns both the density
!  and its Cartesian gradient on the (radial node) x (angular direction)
!  grid.
!
!      den_grid(i,g)    = sum_jl fa2(jl)*Re[ A(i,jl)*Y_jl(u_g) ]
!      grad_grid(i,g,c) = D(i,g)*u_c(g) + T_c(i,g)/r_i
!
!  with
!      D(i,g)   = sum_jl fa2(jl)*Re[ Adot(i,jl)*Y_jl(u_g) ]
!      T_c(i,g) = sum_jl fa2(jl)*Re[ A(i,jl)*G_jl,c(u_g) ]
!
!  where A and Adot are the Lagrange-interpolated value and radial
!  derivative coefficients.  This reproduces exactly the CPU expression
!  in getChargeDensityAtPoint,
!      grad(i) = sum_jl fa2(jl)*Re[ der_rho_in*ylm(kl)*er(i)
!                                 + rho_in(jl)*grad_ylm(kl,i) ]
!  with er = u_g and grad_ylm = G/r.
!
!  der_den_l : (ld_den,jmax) radial derivatives of den_l, e.g. the array
!              returned as the optional 4th argument of
!              getChargeDensity("TotalNew",id,ia,p_grad_den)
!  grad_grid : (ld_grid,num_ang,3)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ld_den, nr, jmax, ld_grid
   integer (kind=IntKind), intent(in) :: ifield_val, ifield_der
   integer (kind=IntKind), intent(in), optional :: mesh_id
!
   real (kind=RealKind), intent(in) :: r_mesh(*)
   real (kind=RealKind), intent(out) :: den_grid(ld_grid,*)
   real (kind=RealKind), intent(out) :: grad_grid(ld_grid,*)
!
   complex (kind=CmplxKind), intent(in) :: den_l(ld_den,*)
   complex (kind=CmplxKind), intent(in) :: der_den_l(ld_den,*)
!
   integer (kind=IntKind) :: i, jl, k, ir, irp, half, mid, ip, g, c, off
   integer (kind=IntKind) :: nr_l, jm_l, fv, fd, ldd
!
   real (kind=RealKind) :: x, num, den, t0, rinv
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('calDensityGradOnAngularGrid','module not initialized')
   else if (.not.gradReady) then
      call ErrorHandler('calDensityGradOnAngularGrid',                   &
                        'gradient tables not built; initialize with needGrad=.true.')
   else if (nr < n_inter .or. nr > nr_max) then
      call ErrorHandler('calDensityGradOnAngularGrid','nr out of range',nr,nr_max)
   else if (jmax < 1 .or. jmax > jmax_max) then
      call ErrorHandler('calDensityGradOnAngularGrid','jmax out of range', &
                        jmax,jmax_max)
   else if (ld_grid < nr) then
      call ErrorHandler('calDensityGradOnAngularGrid','ld_grid < nr',ld_grid,nr)
   else if (ld_den < nr) then
      call ErrorHandler('calDensityGradOnAngularGrid','ld_den < nr',ld_den,nr)
   endif
!
   t0 = getTime()
!
   nr_l = nr; jm_l = jmax; fv = ifield_val; fd = ifield_der; ldd = ld_den
!
#ifdef ACCEL
   if (useGPU) then
      if ( .not.present(mesh_id) ) then
!        -------------------------------------------------------------
         call push_radial_mesh_gpu(r_mesh, nr_l)
!        -------------------------------------------------------------
      else if ( mesh_id /= mesh_id_on_gpu ) then
!        -------------------------------------------------------------
         call push_radial_mesh_gpu(r_mesh, nr_l)
!        -------------------------------------------------------------
         mesh_id_on_gpu = mesh_id
      endif
!     ----------------------------------------------------------------
      call push_density_l_gpu(den_l,     ldd, nr_l, jm_l, fv)
      call push_density_l_gpu(der_den_l, ldd, nr_l, jm_l, fd)
!     ----------------------------------------------------------------
      call eval_density_grad_sphere_gpu(nr_l, jm_l, fv, fd, den_grid,   &
                                        grad_grid, ld_grid)
!     ----------------------------------------------------------------
      t_interp = t_interp + (getTime()-t0)
      return
   endif
#endif
!
!  ===================================================================
!  CPU path.  Identical algebra to the GPU path.
!
!  Step 1: shared interpolation stencil and Lagrange weights.
!  ===================================================================
   half = (n_inter-1)/2
   ir = 1
   do i = 1, nr
      x = r_mesh(i)
!     ----------------------------------------------------------------
      call hunt(nr, r_mesh(1:nr), x, ir)
!     ----------------------------------------------------------------
      if ( ir > nr-half ) then
         irp = nr-n_inter+1
      else if ( 2*ir+1 > n_inter ) then
         irp = ir-half
      else
         irp = 1
      endif
      irp = max(1, min(irp, nr-n_inter+1))
      irp_tab(i) = irp
!
      do k = 1, n_inter
         num = ONE
         den = ONE
         do mid = 1, n_inter
            if ( mid == k ) then
               cycle
            endif
            num = num*(x - r_mesh(irp+mid-1))
            den = den*(r_mesh(irp+k-1) - r_mesh(irp+mid-1))
         enddo
         wlag(i,k) = num/den
      enddo
   enddo
!
!  ===================================================================
!  Step 2: radial term D from the DERIVATIVE coefficients.
!  ===================================================================
   do jl = 1, jmax
      do i = 1, nr
         Ar(i,jl) = ZERO
         Ai(i,jl) = ZERO
      enddo
      do k = 1, n_inter
         do i = 1, nr
            ip = irp_tab(i) + k - 1
            Ar(i,jl) = Ar(i,jl) + wlag(i,k)*real(der_den_l(ip,jl),kind=RealKind)
            Ai(i,jl) = Ai(i,jl) + wlag(i,k)*aimag(der_den_l(ip,jl))
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call dgemm('n','n', nr, num_ang, jmax,  ONE, Ar, nr_max,             &
              Wr, jmax_max, ZERO, Dr, nr_max)
   call dgemm('n','n', nr, num_ang, jmax, -ONE, Ai, nr_max,             &
              Wi, jmax_max,  ONE, Dr, nr_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Step 3: value coefficients -> den_grid, and the angular terms T_c.
!  ===================================================================
   do jl = 1, jmax
      do i = 1, nr
         Ar(i,jl) = ZERO
         Ai(i,jl) = ZERO
      enddo
      do k = 1, n_inter
         do i = 1, nr
            ip = irp_tab(i) + k - 1
            Ar(i,jl) = Ar(i,jl) + wlag(i,k)*real(den_l(ip,jl),kind=RealKind)
            Ai(i,jl) = Ai(i,jl) + wlag(i,k)*aimag(den_l(ip,jl))
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call dgemm('n','n', nr, num_ang, jmax,  ONE, Ar, nr_max,             &
              Wr, jmax_max, ZERO, den_grid, ld_grid)
   call dgemm('n','n', nr, num_ang, jmax, -ONE, Ai, nr_max,             &
              Wi, jmax_max,  ONE, den_grid, ld_grid)
!  -------------------------------------------------------------------
!
   do c = 1, 3
!     ----------------------------------------------------------------
      call dgemm('n','n', nr, num_ang, jmax,  ONE, Ar, nr_max,          &
                 Gr(1,1,c), jmax_max, ZERO, Tc, nr_max)
      call dgemm('n','n', nr, num_ang, jmax, -ONE, Ai, nr_max,          &
                 Gi(1,1,c), jmax_max,  ONE, Tc, nr_max)
!     ----------------------------------------------------------------
      off = (c-1)*num_ang
      do g = 1, num_ang
         do i = 1, nr
            rinv = ONE/r_mesh(i)
            grad_grid(i,off+g) = Dr(i,g)*upos_tab(c,g) + Tc(i,g)*rinv
         enddo
      enddo
   enddo
!
   t_interp = t_interp + (getTime()-t0)
!
   end subroutine calDensityGradOnAngularGrid
!  ===================================================================
end module DensityOnGridModule
