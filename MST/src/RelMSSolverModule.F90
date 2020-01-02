!**************************************************************************
! Name: RelMSSolverModule                                                 *
!                                                                         *
! Version 1.0: Jun. 2016 by Xianglin Liu                                  *
!              Apr. 2017 by Xianglin Liu                                  *
!                        modified to add spin-polaried calculation.       *
! Description:                                                            *
!    Solve for the relativistic tau matrix and construct the              *
!    MST Green function with solutions of the Dirac equation from         *
!    RelSSSolverModule.                                                   *
!**************************************************************************
!   subroutine initRelMSSolver(num_latoms,index,lmaxkkr,lmaxphi,&          *
!              lmaxgreen,local_posi,isBxyz, istop, iprint)                 *
!    Initializes variables and allocates the memory needed for the        *
!    routines in this module. Must be called first.                       *
!    input: num_latoms - integer number of atoms on local processor       *
!           index - an integer array of size num_latoms. index(i)         *
!                   is the global index of local atom (i), where i =      *
!                   1,2,..., num_latoms.                                  *
!           lmaxkkr - an integer array of size num_latoms; is the         *
!                     angular momentum cutoff of KKR matrices.            *
!           lmaxphi - always the same as lmaxkkr, artifact from MSSolver  *
!           lmaxgreen - should be 2*lmaxkkr; is the angular momentum      *
!                       cutoff of the spherical harmonic expanded MST     *
!                       Green function.                                   *
!           local_posi - a real array of size 3*num_latoms; posi(1:3,i)   *
!                       is 3D space coordinates of atom i.                *
!           isBxyz - a logical variable for whether use B field           *
!           istop - routine name to stop (character string)               *
!           iprint - iprint = print instruction parameter (integer)       *
!**************************************************************************
module RelMSSolverModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : ZERO, CZERO, CONE, TEN2m8, HALF, SQRTm1, PI2, PI
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : GridStruct
!   use PublicTypeDefinitionsModule, only : NeighborStruct
!   use NeighborModule, only : getNeighbor, sortNeighbors
!#ifdef TIMING
   use TimerModule, only : getTime
   use RelSSSolverModule, only : getRelSineMatrix,getRelCosineMatrix
   use RelSSSolverModule, only : computeRelSingleSiteDOS, getRelOmegaHatMatrix
   use RelSSSolverModule, only : getRelkyk, getRelkyk_Bx, getRelkyk_By, getRelkyk_Bz
   use RelSSSolverModule, only : calNum_kyk
   use PhysParamModule, only : LightSpeed
   use BesselModule,only: SphericalBessel, SphericalNeumann
!#endif

public :: initRelMSSolver, &
          endRelMSSolver,  &
          computeRelMST,   &
          getRelMSDOS,      &
          getRelMSGreenfunction

private
   logical :: Initialized = .false.
   logical :: isRealSpace = .false.
   logical :: isScreenKKR = .false.
   logical :: InitializedFactors = .false.

   integer (kind=IntKind), allocatable :: print_instruction(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: kmax_kkr(:)
   integer (kind=IntKind), allocatable :: kmax_phi(:)
   integer (kind=IntKind), allocatable :: iend_list(:)
   integer (kind=IntKind) :: lmax_phi_max, kmax_kkr_max, lmax_green_max
   integer (kind=IntKind) :: kmax_phi_max, kmax_green_max, iend_max
   integer (kind=IntKind) :: MaxPrintLevel

   integer (kind=IntKind) :: LocalNumAtoms
   character (len=50) :: stop_routine
   logical :: is_Bxyz=.false.

   real (kind=RealKind) :: Kvec(3)
   real (kind=RealKind), parameter :: Me = 0.5d0

   type MSTStruct
      integer :: lmax !this is lmax_green
      integer :: lamax !this is 2*(lmax+1)^2
      integer :: iend
      real (kind=RealKind), pointer :: dos(:)
      complex (kind=CmplxKind), pointer :: dos_rj(:,:,:)
      complex (kind=CmplxKind), pointer :: green(:,:,:)  
      ! pointer used as allocatable, target
      ! Stores the multiple scattering component of the Green fucntion (Expanded by Y_lm), 
      ! the third dimension is for density+Magnetization
   end type MSTStruct

   type (MSTStruct), allocatable :: mst(:)
   type(GridStruct), pointer :: Grid
!   complex (kind=CmplxKind), allocatable, target :: wspace(:),  gspace(:)
!   complex (kind=CmplxKind), allocatable, target :: wspacep(:), gspacep(:)

!   complex (kind=CmplxKind), allocatable, target :: gfws_comp(:)
!   complex (kind=CmplxKind), allocatable, target :: dosws_comp(:)
contains
   include '../lib/arrayTools.F90'

!======================================================================
   subroutine initRelMSSolver(num_latoms,index,lmaxkkr,lmaxphi, lmaxgreen,&
                              local_posi,isBxyz, istop, iprint)
!======================================================================
   !num_latoms is used for LSMS method
   use SystemModule, only  : getNumAtoms, getAtomPosition, getLmaxMax!, getBravaisLattice !xianglin
   !
   use ScfDataModule, only : isKKR, isScreenKKR, isLSMS
   !
!   use RSpaceStrConstModule, only : initRSpaceStrConst
   !
!   use StrConstModule, only : initStrConst !xianglin
   !
!   use ClusterMatrixModule, only : initClusterMatrix
   !
   use CrystalMatrixModule, only : initCrystalMatrix

   implicit none

   character (len=*), intent(in) :: istop
   integer (kind=IntKind), intent(in) :: num_latoms
   integer (kind=IntKind), intent(in) :: index(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxkkr(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxphi(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxgreen(num_latoms)
   logical, intent(in) :: isBxyz
   integer (kind=IntKind), intent(in) :: iprint(num_latoms)

   real (kind=RealKind), intent(in) :: local_posi(3,num_latoms)

   integer (kind=IntKind) :: i!, NumSysAtoms !commented by xianglin
!   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: global_posi(:,:)

   integer (kind=IntKind) :: rel=2, cant=2 !Default parameters for relativistic calculation

!  -------------------------------------------------------------------
   call initParameters(num_latoms,lmaxkkr,lmaxphi,lmaxgreen,istop,iprint)
!  -------------------------------------------------------------------
   is_Bxyz = isBxyz

   if ( isKKR()) then !initStrConst commented because it has been included in initCrystalMatrix in MST2.1
!      NumSysAtoms = getNumAtoms()
!      allocate( global_posi(1:3,1:NumSysAtoms) )
!      do i = 1, NumSysAtoms
!         global_posi(1:3,i) = getAtomPosition(i)
!      enddo
!      bravais(1:3,1:3) = getBravaisLattice()
!      call initStrConst(lmax_phi_max, NumSysAtoms, global_posi, bravais, istop, MaxPrintLevel)
      call initCrystalMatrix( LocalNumAtoms, cant, lmaxkkr, rel, istop, iprint)
      isRealSpace=.false.
   else if (isLSMS()) then
!      call initRSpaceStrConst(lmax_phi_max,istop,-1)
!      call initClusterMatrix(num_latoms,index,lmaxkkr,lmaxphi,local_posi,cant,rel,istop,iprint)
      isRealSpace=.true.
   else
      call ErrorHandler('initRelMSSolver','unimplemented Non-KKR method is called')
   endif

   Initialized = .true.
   end subroutine initRelMSSolver

!  ===================================================================
   subroutine initParameters(num_atoms, lmaxkkr, lmaxphi, lmaxgreen,  &
                              istop, iprint)
!  ===================================================================
   use RadialGridModule, only : getNumRmesh, getMaxNumRmesh
   implicit none
!
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: lmaxkkr(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxphi(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxgreen(num_atoms)
   integer (kind=IntKind), intent(in) :: iprint(num_atoms)
   integer (kind=LongIntKind) :: wspace_size, gspace_size
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind) :: i, lmax_max, jmax, iend, kmax
!
   if (Initialized) then
!     ----------------------------------------------------------------
      call WarningHandler('initRelMSSolver',                   &
                  'RelMSSolverModule has already been initialized')
!     ----------------------------------------------------------------
      return
   else if (num_atoms < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initRelMSSolver','num_atoms < 1',num_atoms)
!     ----------------------------------------------------------------
   endif
!
   LocalNumAtoms = num_atoms
   stop_routine = istop
!
   allocate( print_instruction(LocalNumAtoms) )
   allocate( iend_list(LocalNumAtoms) )
   allocate( lmax_kkr(LocalNumAtoms) )
   allocate( lmax_phi(LocalNumAtoms) )
   allocate( kmax_kkr(LocalNumAtoms) )
   allocate( kmax_phi(LocalNumAtoms) )
!
   kmax_kkr_max = 1
   kmax_phi_max = 1
   kmax_green_max = 1
   lmax_phi_max = 0
   lmax_green_max = 0
   lmax_max = 0
   do i=1, LocalNumAtoms
      if (lmaxkkr(i) .ne. lmaxphi(i)) then
         print*,"lmaxkkr(i)=", lmaxkkr(i), " lmaxphi(i)=", lmaxphi(i)
         call ErrorHandler('initRelMSSolver','lmax_kkr .ne. lmax_phi, unimplemented')
      endif
      lmax_kkr(i) = lmaxkkr(i)
      kmax_kkr(i) = (lmaxkkr(i)+1)**2
      lmax_phi(i) = lmaxphi(i)
      kmax_phi(i) = (lmaxphi(i)+1)**2
      print_instruction(i) = iprint(i)
      kmax_kkr_max = max(kmax_kkr_max, (lmaxkkr(i)+1)**2)
      kmax_phi_max = max(kmax_phi_max, (lmaxphi(i)+1)**2)
      kmax_green_max = max(kmax_green_max, (lmaxgreen(i)+1)**2)
      lmax_phi_max = max(lmax_phi_max, lmaxphi(i))
      lmax_green_max = max(lmax_green_max, lmaxgreen(i))
      lmax_max = max(lmax_max, lmaxgreen(i), lmaxkkr(i), lmaxphi(i))
   enddo
!  -------------------------------------------------------------------
   MaxPrintLevel = maxval(print_instruction(1:LocalNumAtoms))
!
   allocate( mst(LocalNumAtoms) )
   do i=1, LocalNumAtoms
      kmax = (lmaxgreen(i)+1)**2
      jmax = (lmaxgreen(i)+1)*(lmaxgreen(i)+2)/2
      iend = getNumRmesh(i)
      iend_list(i) = iend
      mst(i)%lmax = lmaxgreen(i)
      mst(i)%lamax = 2*(lmaxkkr(i)+1)**2 ! no distinction between lmaxkkr and lmaxphi
      mst(i)%iend = iend
      allocate( mst(i)%dos(4) )
      allocate( mst(i)%dos_rj(iend,jmax,4) )
      allocate( mst(i)%green(iend,kmax,4) )
      mst(i)%dos=ZERO
      mst(i)%dos_rj=CZERO
      mst(i)%green=CZERO
   enddo
   iend_max = getMaxNumRmesh()
!   wspace_size = 4*iend_max*kmax_phi_max*kmax_kkr_max
!   gspace_size = 4*iend_max*kmax_green_max*kmax_phi_max
!   allocate( wspace(wspace_size) )
!   allocate( gspace(gspace_size) )
!
   end subroutine initParameters
!  ===================================================================

   !=================================================================
   subroutine endRelMSSolver
   !=================================================================
!   use StrConstModule, only :  endStrConst
!   use RSpaceStrConstModule, only : endRSpaceStrConst
!   use ClusterMatrixModule, only : endClusterMatrix
   use CrystalMatrixModule, only : endCrystalMatrix

   implicit none
!
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endRelMSSolver',                      &
                        'RelMSSolverModule has not been initialized')
!     ----------------------------------------------------------------
   endif
!
   if (isRealSpace) then
!     ----------------------------------------------------------------
!      call endRSpaceStrConst()
      call ErrorHandler('endRelMSSolver',                      &
                        'Relativistic LSMS has not been implemented yet')
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
!      call endStrConst() !xianglin
      call endCrystalMatrix()
!     ----------------------------------------------------------------
   endif
!
   deallocate( print_instruction )
   deallocate( iend_list)
   deallocate( lmax_kkr )
   deallocate( lmax_phi )
   deallocate( kmax_kkr )
   deallocate( kmax_phi )
   do i=1, LocalNumAtoms
      deallocate( mst(i)%dos )
      deallocate( mst(i)%dos_rj)
      deallocate( mst(i)%green )
   enddo
   deallocate(mst)
   Initialized = .false.
   isRealSpace = .false.

   end subroutine endRelMSSolver

   !=================================================================
   subroutine computeRelMST(Ek)
   !=================================================================
   use CrystalMatrixModule, only : calCrystalMatrix,getTau
   use RelSSSolverModule, only : SingleDiracScattering, getRelScatteringMatrix

   implicit none

   complex (kind=CmplxKind), intent(in) :: Ek
   logical :: tau_needed

   integer (kind=IntKind) :: id
   do id = 1, LocalNumAtoms
      call SingleDiracScattering (id,Ek)
   enddo
   call calCrystalMatrix(Ek, getRelScatteringMatrix, tau_needed=.true.)
   !open (102,file="tau_matrix",action="write")
   !write (102,*) getTau(1,1) 
   !close (102)

   do id = 1, LocalNumAtoms
      call computeRelMSGF_part(id,Ek)
      call computeRelDOS_part(id,Ek)
      call computeRelVolumeIDOS(id)
      if (is_Bxyz) then
         call computeRelVolumeMMIDOS(id)
      endif
   enddo
   !open (102,file="MST_dos_rj",action="write")
   !write (102,*) mst(1)%dos_rj(:,:,1)
   !close (102)
   !print*,"mst dos=", mst(1)%dos(1)
   end subroutine computeRelMST

   !=================================================================
   subroutine computeRelMSGF_part(ia,Ek)
   !=================================================================
!   use ClusterMatrixModule, only : getClusterKau => getKau
   use CrystalMatrixModule, only : getCrystalKau => getKau
   use RadialGridModule, only : getGrid
   use RelSSSolverModule, only : getRelRegSolutionSine, getRelRegSolutionCosine
   use RelSSSolverModule, only : getldex,getldex_bar

   implicit none

   integer (kind=IntKind), intent(in) :: ia
   complex (kind=CmplxKind), intent(in) :: Ek
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: sz, loc1, loc2, nr
   integer (kind=IntKind) :: nofr, km
   integer (kind=IntKind) :: kap_lp, my_lp, l_lp, kapbar_lp, lambar_lp, lbar_lp
   integer (kind=IntKind) :: l_lpp, lbar_lpp, lmax_p2
   integer (kind=IntKind) :: i, i3, ir, lp, lpp, lamax, k_1, k_2, kl1, lmax, lmax_green
!   integer (kind=IntKind) :: nj3_i3
!   integer (kind=IntKind), pointer :: kj3_i3(:)
!   real (kind=RealKind), pointer :: cgnt_i3(:)
   complex (kind=CmplxKind), pointer :: ss_mat_r(:,:,:)
   complex (kind=CmplxKind), pointer :: cc_mat_r(:,:,:)
   complex (kind=CmplxKind), pointer :: green(:,:), green_Bxyz(:,:,:)
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind) :: energy !Total energy
   complex (kind=CmplxKind) :: bjtmp( 0:lmax_kkr(ia)+2 )
   complex (kind=CmplxKind) :: bntmp( 0:lmax_kkr(ia)+2 )
   complex (kind=CmplxKind) :: mid_matrix( mst(ia)%lamax,mst(ia)%lamax,4)
   complex (kind=CmplxKind) :: prefactor,prefactor_bar
   !complex (kind=CmplxKind), pointer :: t_matrix(:,:)
   !complex (kind=CmplxKind), pointer :: s_matrix(:,:)
   !complex (kind=CmplxKind), pointer :: c_matrix(:,:)
!   complex (kind=CmplxKind), pointer :: OmegaHat(:,:)
!   complex (kind=CmplxKind), pointer :: d_s_matrix(:,:)
!   complex (kind=CmplxKind), pointer :: d_c_matrix(:,:)
   complex (kind=CmplxKind), pointer :: d_ss_mat_r(:,:,:)
   complex (kind=CmplxKind), pointer :: d_cc_mat_r(:,:,:)
   real (kind=RealKind), pointer :: kyk_loc(:,:,:,:)
   real (kind=RealKind), pointer :: kyk_loc_Bx(:,:,:,:), kyk_loc_By(:,:,:,:), kyk_loc_Bz(:,:,:,:)
   integer (kind=IntKind) :: n_kyk(mst(ia)%lamax,mst(ia)%lamax)
   integer (kind=IntKind) :: ind_kyk((mst(ia)%lmax+1)**2,mst(ia)%lamax,mst(ia)%lamax)
   integer (kind=IntKind) :: n_kyk_bar(mst(ia)%lamax,mst(ia)%lamax)
   integer (kind=IntKind) :: ind_kyk_bar((mst(ia)%lmax+1)**2,&
                             mst(ia)%lamax,mst(ia)%lamax)
   integer (kind=IntKind), allocatable :: n_kyk_Bxyz(:,:,:)
   integer (kind=IntKind), allocatable :: ind_kyk_Bxyz(:,:,:,:)
   integer (kind=IntKind), allocatable :: n_kyk_bar_Bxyz(:,:,:)
   integer (kind=IntKind), allocatable :: ind_kyk_bar_Bxyz(:,:,:,:)

   complex (kind=CmplxKind) :: tmpI
   complex (kind=CmplxKind), allocatable :: bjl(:,:), bnl(:,:)
   complex (kind=CmplxKind), pointer :: kau00(:,:,:),OmegaHat(:,:) !for relativistic, the third component is 1d
   integer (kind=IntKind), pointer :: ldex(:),ldex_bar(:)

   if (is_Bxyz) then
      allocate( n_kyk_Bxyz(mst(ia)%lamax,mst(ia)%lamax,3) )
      allocate( ind_kyk_Bxyz((mst(ia)%lmax+1)**2,mst(ia)%lamax,mst(ia)%lamax,3) )
      allocate( n_kyk_bar_Bxyz(mst(ia)%lamax,mst(ia)%lamax,3) )
      allocate( ind_kyk_bar_Bxyz((mst(ia)%lmax+1)**2,mst(ia)%lamax,mst(ia)%lamax,3) )
   endif

   if (isRealSpace) then
      call ErrorHandler('computeRelMSGF_part',                      &
                        'Relativistic LSMS not implemented yet')
!      kau00 => getClusterKau(ia)   ! Kau00 = kappa**2 * S^{-1} * [Tau00 - t_matrix] * S^{-1*}
   else
      kau00 => getCrystalKau(ia)   ! Kau00 = kappa**2 * S^{-1} * [Tau00 - t_matrix] * S^{-1*}
   endif

!   OmegaHat => getRelOmegaHatMatrix(ia) ! to test single-site scattering


!   open (102,file="kau00",action="write")
!   write (102,*) kau00
!   close (102)
   !s_matrix => Scatter(ia)%sin_mat
   !c_matrix => Scatter(ia)%cos_mat
   ss_mat_r => getRelRegSolutionSine(ia)
   cc_mat_r => getRelRegSolutionCosine(ia)

   !Right now no distinction between dot and undot sine matrix, to be generalized
   !d_s_matrix => Scatter(ia)%sin_mat
   !d_c_matrix => Scatter(ia)%cos_mat
   d_ss_mat_r => ss_mat_r
   d_cc_mat_r => cc_mat_r

   kappa = sqrt(2.d0*Me*Ek + Ek**2/LightSpeed**2)
   energy = Ek + Me*LightSpeed**2

   lamax = mst(ia)%lamax
   lmax_green = mst(ia)%lmax !note lmax here is actually lmax_green
   !Note lmax_p

   kyk_loc => getRelkyk(ia)

 !  open (101,file="kyk_mst",action="write")
 !  write (101,*) kyk_loc(:,1,:,2)
 !  close(101)
!   print*, "mst%lmax=", mst(ia)%lmax
   call calNum_kyk(lmax_green,lamax,kyk_loc(:,1,:,:),n_kyk,ind_kyk)
   call calNum_kyk(lmax_green,lamax,kyk_loc(:,2,:,:),n_kyk_bar,ind_kyk_bar)
   if (is_Bxyz) then
      kyk_loc_Bx => getRelkyk_Bx(ia)
      kyk_loc_By => getRelkyk_By(ia)
      kyk_loc_Bz => getRelkyk_Bz(ia)
      call calNum_kyk(lmax_green,lamax,kyk_loc_Bx(:,1,:,:),n_kyk_Bxyz(:,:,1),ind_kyk_Bxyz(:,:,:,1))
      call calNum_kyk(lmax_green,lamax,kyk_loc_Bx(:,2,:,:),n_kyk_bar_Bxyz(:,:,1),ind_kyk_bar_Bxyz(:,:,:,1))
      call calNum_kyk(lmax_green,lamax,kyk_loc_By(:,1,:,:),n_kyk_Bxyz(:,:,2),ind_kyk_Bxyz(:,:,:,2))
      call calNum_kyk(lmax_green,lamax,kyk_loc_By(:,2,:,:),n_kyk_bar_Bxyz(:,:,2),ind_kyk_bar_Bxyz(:,:,:,2))
      call calNum_kyk(lmax_green,lamax,kyk_loc_Bz(:,1,:,:),n_kyk_Bxyz(:,:,3),ind_kyk_Bxyz(:,:,:,3))
      call calNum_kyk(lmax_green,lamax,kyk_loc_Bz(:,2,:,:),n_kyk_bar_Bxyz(:,:,3),ind_kyk_bar_Bxyz(:,:,:,3))
      green_Bxyz => mst(ia)%green(:,:,2:4)
   endif
!      print*,"MST n_kyk="
!      print*, n_kyk(2,2)
!      print*,"ind_kyk="
!      print*, ind_kyk(:,2,2)

   nr = iend_list(ia)
   Grid => getGrid(ia)
   mst(ia)%green = CZERO
   green => mst(ia)%green(:,:,1)
   lmax_p2=lmax_kkr(ia)+2
   allocate( bjl(0:lmax_p2,1:nr), bnl(0:lmax_p2,1:nr) )
   do ir=1,nr
      call SphericalBessel(lmax_p2,kappa*Grid%r_mesh(ir),bjl(:,ir))
      call SphericalNeumann(lmax_p2,kappa*Grid%r_mesh(ir),bnl(:,ir))
   enddo
   ldex => getldex()
   ldex_bar => getldex_bar()

!   green=CZERO
   do ir=1, nr
   !note divided by pi in calculating DOS or density is not included in prefactor
      prefactor=kappa*(energy+Me*LightSpeed**2)/LightSpeed**2
      !*kappa removed since it has been included in Kau00
      prefactor_bar=kappa*(energy+Me*LightSpeed**2)/LightSpeed**2*&
                    (kappa**2*LightSpeed**2)/(energy+Me*LightSpeed**2)**2
      bjtmp=bjl(:,ir)
      bntmp=bnl(:,ir)
      call get_midM(lamax,ss_mat_r(:,:,ir),&
           cc_mat_r(:,:,ir),kau00(:,:,1)/kappa,mid_matrix,d_ss_mat_r(:,:,ir),&
           d_cc_mat_r(:,:,ir)) !mid_matrix has been initialized in get_midM
      do lp=1,lamax
   !      kap_lp=kapdex(lp)
   !      my_lp=mydex(lp)
         l_lp=ldex(lp)
         lbar_lp=ldex_bar(lp)
         do lpp=1,lamax
   !         kap_lpp=kapdex(lpp)
   !         my_lpp=mydex(lpp)
            l_lpp=ldex(lpp)
            lbar_lpp=ldex_bar(lpp)

            tmpI=bntmp(l_lpp)*mid_matrix(lpp,lp,1)*bntmp(l_lp)
!  if(ir==nr) then
!  print*,"tmpI 1", tmpI, kyk_loc(kl1,1,lpp,lp), n_kyk(lpp,lp)
!  endif
            tmpI=tmpI - bntmp(l_lpp)*mid_matrix(lpp,lp,2)*bjtmp(l_lp)
!  if(ir==nr) then
!  print*,"tmpI 2", tmpI
!  endif
            tmpI=tmpI - bjtmp(l_lpp)*mid_matrix(lpp,lp,3)*bntmp(l_lp)
!  if(ir==nr) then
!  print*,"tmpI 3", tmpI
!  endif
            tmpI=tmpI + bjtmp(l_lpp)*mid_matrix(lpp,lp,4)*bjtmp(l_lp)
!  if(ir==nr) then
!  print*,"tmpI 4", tmpI
!  endif
            tmpI=tmpI*prefactor
            do i=1,n_kyk(lpp,lp)
               kl1=ind_kyk(i,lpp,lp)
               green(ir,kl1)=green(ir,kl1) + tmpI*kyk_loc(kl1,1,lpp,lp)
            enddo
            if (is_Bxyz) then
               do i=1,n_kyk_Bxyz(lpp,lp,1)
                  kl1=ind_kyk_Bxyz(i,lpp,lp,1)
                  green_Bxyz(ir,kl1,1)=green_Bxyz(ir,kl1,1) + tmpI*kyk_loc_Bx(kl1,1,lpp,lp)
               enddo
               do i=1,n_kyk_Bxyz(lpp,lp,2)
                  kl1=ind_kyk_Bxyz(i,lpp,lp,2)
                  green_Bxyz(ir,kl1,2)=green_Bxyz(ir,kl1,2) + tmpI*kyk_loc_By(kl1,1,lpp,lp)
               enddo
               do i=1,n_kyk_Bxyz(lpp,lp,3)
                  kl1=ind_kyk_Bxyz(i,lpp,lp,3)
                  green_Bxyz(ir,kl1,3)=green_Bxyz(ir,kl1,3) + tmpI*kyk_loc_Bz(kl1,1,lpp,lp)
               enddo
            endif
!  if(ir==nr) then
!  print*,"green(1,1)", green(1,1)
!  endif
   ! the lower "bar" part
            tmpI=bntmp(lbar_lpp)*mid_matrix(lpp,lp,1)*bntmp(lbar_lp)
            tmpI=tmpI - bntmp(lbar_lpp)*mid_matrix(lpp,lp,2)*bjtmp(lbar_lp)
            tmpI=tmpI - bjtmp(lbar_lpp)*mid_matrix(lpp,lp,3)*bntmp(lbar_lp)
            tmpI=tmpI + bjtmp(lbar_lpp)*mid_matrix(lpp,lp,4)*bjtmp(lbar_lp)
            tmpI=tmpI*prefactor_bar
            do i=1,n_kyk_bar(lpp,lp)
               kl1=ind_kyk_bar(i,lpp,lp)
               green(ir,kl1)=green(ir,kl1) + tmpI*kyk_loc(kl1,2,lpp,lp)
            enddo
            if (is_Bxyz) then
               do i=1,n_kyk_Bxyz(lpp,lp,1)
                  kl1=ind_kyk_Bxyz(i,lpp,lp,1)
                  green_Bxyz(ir,kl1,1)=green_Bxyz(ir,kl1,1) + tmpI*kyk_loc_Bx(kl1,2,lpp,lp)
               enddo
               do i=1,n_kyk_Bxyz(lpp,lp,2)
                  kl1=ind_kyk_Bxyz(i,lpp,lp,2)
                  green_Bxyz(ir,kl1,2)=green_Bxyz(ir,kl1,2) + tmpI*kyk_loc_By(kl1,2,lpp,lp)
               enddo
               do i=1,n_kyk_Bxyz(lpp,lp,3)
                  kl1=ind_kyk_Bxyz(i,lpp,lp,3)
                  green_Bxyz(ir,kl1,3)=green_Bxyz(ir,kl1,3) + tmpI*kyk_loc_Bz(kl1,2,lpp,lp)
               enddo
            endif
         enddo
      enddo
      green(ir,:) = green(ir,:)*Grid%r_mesh(ir)**2 !fit the GF notation of YangWang's code
      if (is_Bxyz) then
         green_Bxyz(ir,:,:) = green_Bxyz(ir,:,:)*Grid%r_mesh(ir)**2
      endif
   enddo
   deallocate(bjl,bnl)
 !  open (101,file="mid_M_mst",action="write")
 !  write (101,*) green(:,1)!mid_matrix(:,:,1)
 !  close(101)

   if (is_Bxyz) then
      deallocate( n_kyk_Bxyz )
      deallocate( ind_kyk_Bxyz )
      deallocate( n_kyk_bar_Bxyz )
      deallocate( ind_kyk_bar_Bxyz )
   endif
   end subroutine computeRelMSGF_part

   !=================================================================
   subroutine computeRelDOS_part(atom,Ek)
   !=================================================================
!   use GroupCommModule, only : GlobalSumInGroup
!   use RadialGridModule, only : getGrid
   use IntegerFactorsModule, only : lofk, mofk, &
       jofk, lofj, mofj, kofj, m1m
!
   implicit none

!   logical, intent(in), optional :: add_highl_fec   ! Add high L contribution by free electron
   integer (kind=IntKind), intent(in) :: atom
   complex (kind=CmplxKind), intent(in) :: Ek
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: ir, jl, kl, klc, l, m, kl1, kl1p, kl2, kl2c, kl2p, &
         kl2pc, m1, m2, m2p, ma, i, l2
   integer (kind=IntKind) :: js, ia, nr, lm, jm, sz, loc1, loc2
   integer (kind=IntKind) :: lmax_dos, jmax_dos, kmax_kkr_loc
   complex (kind=CmplxKind), pointer :: dos(:,:), dos_Bxyz(:,:,:)
   complex (kind=CmplxKind) :: cfac, kappa, energy
   complex (kind=CmplxKind), pointer :: green(:,:), green_Bxyz(:,:,:)

   id = atom
   kappa=sqrt(2.d0*Me*Ek + Ek**2/LightSpeed**2)

!   if (.not. allocated(wks_dos)) then
!      sz = 0
!      do ia = 1, LocalNumAtoms
!         nr = Scatter(ia)%numrs_cs
!         lm = Scatter(ia)%lmax_green
!         jm = (lm+1)*(lm+2)/2
!         sz = sz + nr*jm
!      enddo
!      allocate( wks_dos(sz) )
!      loc1 = 0
!      do ia = 1, LocalNumAtoms
!         nr = Scatter(ia)%numrs_cs
!         lm = Scatter(ia)%lmax_green
!         jm = (lm+1)*(lm+2)/2
!         loc2 = loc1 + nr*jm
!         Scatter(ia)%dos => aliasArray2_c(wks_dos(loc1+1:loc2),nr,jm)
!         loc1 = loc2
!      enddo
!   endif

   lmax_dos = mst(id)%lmax
   jmax_dos = (lmax_dos+1)*(lmax_dos+2)/2
   green => mst(id)%green(:,:,1)
   mst(id)%dos_rj=CZERO
   dos => mst(id)%dos_rj(:,:,1)
   nr = iend_list(id)
   dos=CZERO
   do jl = 1, jmax_dos
      kl = kofj(jl)
      m = mofj(jl)
      klc = kl-2*m
      do ir = 1, nr
         dos(ir,jl) = green(ir,kl) - m1m(m)*conjg(green(ir,klc))
      enddo
   enddo
   cfac = sqrtm1/PI2 !divid by 2*Pi, the 2 is from 1/2*(green-conjg(green))
!
   dos = cfac*dos

   if (is_Bxyz) then
      green_Bxyz => mst(id)%green(:,:,2:4)
      dos_Bxyz => mst(id)%dos_rj(:,:,2:4)
      do jl = 1, jmax_dos
         kl = kofj(jl)
         m = mofj(jl)
         klc = kl-2*m
         do i = 1, 3
            do ir = 1, nr   
               dos_Bxyz(ir,jl,i) = green_Bxyz(ir,kl,i) - m1m(m)*conjg(green_Bxyz(ir,klc,i))
            enddo
         enddo
      enddo
      cfac = sqrtm1/PI2 !divid by 2*Pi, the 2 is from 1/2*(green-conjg(green))
   !
      dos_Bxyz = cfac*dos_Bxyz
   endif
!
!  add high l free space part. Already included in RelSSSolverModule
!    Grid => getGrid(id)
!   if ( present(add_highl_fec) ) then
!      if (add_highl_fec) then
!         TmpSpace = CZERO
!         do l = lofk(kmax_kkr_loc), 0, -1
!            do ir = 1, nr
!               TmpSpace(ir) = TmpSpace(ir) + (2*l+1)*bjl(ir,l)**2/kappa**2
!            enddo
!         enddo
!         cfac = kappa*sqrt(PI)/(2.0d0*PI**2)*(4.d0*energy/LightSpeed**2)
!         ! =(4*W/c^2)*kappa/(4.0d0*PI**2)/Y00
!         do ir = 1, nr
!            dos(ir,1) = dos(ir,1) + cfac*(Grid%r_mesh(ir)**2-TmpSpace(ir))
!         enddo
!      endif
!   endif

   end subroutine computeRelDOS_part

!  ===================================================================
   subroutine computeRelVolumeIDOS(id)
!  ===================================================================
!  
!  This function returns single site DOS for a given energy on the real
!  energy axis.
!  
!  ===================================================================
   use IntegrationModule, only : calIntegration
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration


   implicit none
!  
   integer (kind=IntKind), intent(in) :: id
!  
   real (kind=RealKind) :: dos, dos_mt
   integer (kind=IntKind) :: ir, nr, jmax_dos
   real (kind=RealKind) :: dos_r(iend_list(id)), dos_int(iend_list(id))
   complex (kind=CmplxKind) :: energy !Total energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   logical :: SphericalInt=.false., truncated !to test getVolumeIntegration, .false. by default

   nr=iend_list(id)


   Grid => getGrid(id)
   dos_r_jl => mst(id)%dos_rj(:,:,1)

!   mst(id)%dos = 0.d0
   if (SphericalInt) then
      do ir=1,nr
         dos_r(ir)=Grid%r_mesh(ir)*Grid%r_mesh(ir)*Real(dos_r_jl(ir,1))*Sqrt(4.d0*PI)
         !r^2*dos*Y_0/Y_0^2, note scatter%dos are basically GF expanded on spherical harmonics
      enddo
      dos_int=0.d0
      call calIntegration(0, nr, Grid%r_mesh, dos_r(:), dos_int(:))
      dos=dos_int(nr)
      mst(id)%dos(1) = dos
   else
      jmax_dos = size(dos_r_jl,2)
      mst(id)%dos(1) = getVolumeIntegration( id, nr, Grid%r_mesh,   &
                           jmax_dos, 2, dos_r_jl, dos_mt, truncated=.true. )
   endif

   end subroutine computeRelVolumeIDOS

!  ===================================================================
   subroutine computeRelVolumeMMIDOS(id)
!  ===================================================================
!  
!  This function returns single site DOS for a given energy on the real
!  energy axis.
!  
!  ===================================================================
   use IntegrationModule, only : calIntegration
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration


   implicit none
!  
   integer (kind=IntKind), intent(in) :: id
!  
   real (kind=RealKind) ::  dos_mt(3)
   integer (kind=IntKind) :: ir, nr, jmax_dos, i
   real (kind=RealKind) :: dos_r(iend_list(id)), dos_int(iend_list(id))
   complex (kind=CmplxKind) :: energy !Total energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:,:)
   logical :: SphericalInt=.false., truncated !to test getVolumeIntegration, .false. by default

   nr=iend_list(id)


   Grid => getGrid(id)
   dos_r_jl => mst(id)%dos_rj(:,:,2:4)

   jmax_dos = size(dos_r_jl,2)
   if (is_Bxyz) then
      do i=1,3
         mst(id)%dos(1+i) = getVolumeIntegration( id, nr, Grid%r_mesh,   &
                        jmax_dos, 2, dos_r_jl(:,:,i), dos_mt(i), truncated=.true. )
      enddo
   else
      call ErrorHandler ('computeRelVolumeMMIDOS',&
      'non-magnetic calculation, but magnetic moment density DOS is called')
   endif

   end subroutine computeRelVolumeMMIDOS

   !=================================================================
   subroutine get_midM(lamax,ss_mat_r,cc_mat_r,OmegaHat,mid_matrix,&
              d_ss_mat_r,d_cc_mat_r)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : CZERO, CONE, SQRTm1

   implicit none

   integer (kind=IntKind),intent(in) :: lamax
   integer (kind=IntKind) :: i
   complex (kind=CmplxKind) :: OmegaHat(lamax,lamax)
   complex (kind=CmplxKind) :: ss_mat_r(lamax,lamax),cc_mat_r(lamax,lamax),Omat(lamax,lamax)
   complex (kind=CmplxKind) :: d_ss_mat_r(lamax,lamax),d_cc_mat_r(lamax,lamax)
   complex (kind=CmplxKind) :: mid_matrix(lamax,lamax,4)
   complex (kind=CmplxKind) :: mat_tmp(lamax,lamax)

   mat_tmp=OmegaHat
   !mat_tmp=(S^t*(i*S-C))^-1

   Omat=(0.d0,0.d0)
   call zgemm('n','n',lamax,lamax,lamax,cone,&
   ss_mat_r, lamax, mat_tmp, lamax, czero, Omat, lamax)
   !Omat=ss_r*OmegaHat

   mid_matrix(:,:,1)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_ss_mat_r, lamax, czero, mid_matrix(:,:,1), lamax)
   !mid_matrix(1)=ss_r*OmegaHat*d_ss_r

   mid_matrix(:,:,2)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_cc_mat_r, lamax, czero, mid_matrix(:,:,2), lamax)
   !mid_matrix(2)=ss_r*OmegaHat*d_cc_r

   Omat=(0.d0,0.d0)
   call zgemm('n','n',lamax,lamax,lamax,cone,&
   cc_mat_r, lamax, mat_tmp, lamax, czero, Omat, lamax)
   !Omat=cc_r*OmegaHat

   mid_matrix(:,:,3)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_ss_mat_r, lamax, czero, mid_matrix(:,:,3), lamax)
   !mid_matrix(3)=cc_r*OmegaHat*d_ss_r

   mid_matrix(:,:,4)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_cc_mat_r, lamax, czero, mid_matrix(:,:,4), lamax)
   !mid_matrix(4)=cc_r*OmegaHat*d_cc_r

   end subroutine get_midM

   !=================================================================
   function getRelMSDOS(atom) result(dos_loc)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   real (kind=RealKind), pointer :: dos_loc(:)

   if (.not.Initialized) then
      call ErrorHandler('getRelDOS', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelDOS', 'invalid number of &
                         local atoms', LocalNumAtoms)
   endif


   dos_loc => mst(atom)%dos(:)

   end function getRelMSDOS

   !=================================================================
   function getRelMSGreenFunction(atom) result(green_loc)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: green_loc(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelMSGreenFunction', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelMSGreenFunction', 'invalid number of &
                         local atoms', LocalNumAtoms)
   endif


   green_loc => mst(atom)%green

   end function getRelMSGreenFunction
end module RelMSSolverModule

