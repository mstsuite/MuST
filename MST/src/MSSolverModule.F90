!  *******************************************************************
!  *                                                                 *
!  * NAME:    MSSolverModule                                         *
!  *                                                                 *
!  * VERSION: 1.0                                                    *
!  * DATE:    05/21/13                                               *
!  *          Notes :                                                *
!  *                                                                 *
!  * DESCRIPTION:                                                    *
!  * Module for solving Multiple Scattering Theory equations         *
!  *                                                                 *
!  * EXTERNAL MODULE DEPENDENCE                                      *
!  *    KindParamModule                                              *
!  *    MathParamModule                                              *
!  *    ErrorHandlerModule                                           *
!  *    WriteMatrixModule                                            *
!  *    GauntFactorsModule                                           *
!  *    RadialGridModule                                             *
!  *    NeighborModule                                               *
!  *    SpinRotationModule                                           *
!  *    SurfElementsModule                                           *
!  *    SSSolverModule                                               *
!  *                                                                 *
!  * USAGE:                                                          *
!  *    ============================================================ *
!  *    initMSSolver(num_latoms, z, index, lmaxkkr,                  *
!  *                           lmaxphi, lmaxgreen, posi, neighbor,   *
!  *                           pola, cant, rel, istop, iprint)       *
!  *    Purpose: initialize the module for solving cluster MST       *
!  *             equation                                            *
!  *    Input:   num_latoms = no. of atoms that need to be solved    *
!  *                          on local processor                     *
!  *             z          = an integer array of size num_latoms.   *
!  *                          z(i) is the atomic number of local     *
!  *                          atom i, where i = 1,2,..., num_latoms. *
!  *             index      = an integer array of size num_latoms.   *
!  *                          index(i) is the global index of local  *
!  *                          atom i, where i = 1,2,..., num_latoms, *
!  *                          and 1 <= index(i) <= no. of total atoms*
!  *             lmaxkkr    = an integer array of size num_latoms.   *
!  *                          lmaxkkr(i) is the lmax cut off for     *
!  *                          single-site wave function L-index for  *
!  *                          atom i, where i = 1,2,..., num_latoms. *
!  *             lmaxphi    = an integer array of size num_latoms.   *
!  *                          lmaxphi(i) is the lmax cut off for     *
!  *                          single-site wave function expansion for*
!  *                          atom i, where i = 1,2,..., num_latoms. *
!  *             lmaxgreen    = an integer array of size num_latoms. *
!  *                          lmaxgreen(i) is the lmax cut off for   *
!  *                          Green function expansion for atom i,   *
!  *                          where i = 1,2,..., num_latoms.         *
!  *             posi       = a real array of size 3*num_latoms.     *
!  *                          posi(1:3,i) is 3D space coordinates of *
!  *                          atom i, where i = 1,2,..., num_latoms. *
!  *             pola       = 1, if non-spin-polarized;              *
!  *                          2, if spin-polarized.                  *
!  *             cant       = 1, if non-spin-canted;                 *
!  *                          2, if spin-canted.                     *
!  *             Note: Both pola and cant are of integer type.       *
!  *                   If pola = 1 and cant = 2, an error message    *
!  *                   will occur.                                   *
!  *             rel        = 0, if non-relativistic                 *
!  *                          1, if semi-relativistic                *
!  *                          2, if fully-relativistic               *
!  *             istop = routine name to stop (character string)     *
!  *             iprint = print instruction parameter (integer)      *
!  *    Output:  none                                                *
!  *    ============================================================ *
!  *    endMSSolver()                                                *
!  *    Purpose: clean the memory allocated within the module.       *
!  *    Input:   none                                                *
!  *    Output:  none                                                *
!  *    ============================================================ *
!  *    solveMSTeqns(e,is)                                           *
!  *    Purpose: solve multiple scattering theory equations for (all *
!  *             the local atoms.                                    *
!  *    Input:   e      = energy (complex).                          *
!  *             is     = spin index (integer).                      *
!  *    Output:  N/A                                                 *
!  *******************************************************************
module MSSolverModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors
!#ifdef TIMING
   use TimerModule, only : getTime
!#endif
!
public :: initMSSolver,            &
          endMSSolver,             &
          computeMSGreenFunction,  &
          computeMSPDOS,           &
          getMSGreenFunction,      &        ! Returns Green function in Local frame
          getMSGreenMatrix,        &
          getMSCellDOS,            &
          getMSMTSphereDOS,        &
          getMSCellPDOS,           &
          getMSMTSpherePDOS,       &
          getMSGreenFunctionDerivative
!
   interface getDOS
      module procedure getDOS_is, getDOS_sc
   end interface ! getDOS
!
private
!
   logical :: Initialized = .false.
   logical :: InitializedFactors = .false.
!
   integer (kind=IntKind) :: Relativity
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
!
   character (len=50) :: stop_routine
!
   complex (kind=CmplxKind) :: Energy
!
   logical :: isPivoting = .false.
   complex (kind=CmplxKind) :: PivoteEnergy = CZERO
!
   integer (kind=IntKind), allocatable :: print_instruction(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: kmax_kkr(:)
   integer (kind=IntKind), allocatable :: kmax_phi(:)
   integer (kind=Intkind), allocatable :: lofk(:), mofk(:), jofk(:), m1m(:)
   integer (kind=Intkind), allocatable :: lofj(:), mofj(:)
   integer (kind=IntKind) :: lmax_phi_max, kmax_kkr_max, lmax_green_max
   integer (kind=IntKind) :: kmax_phi_max, kmax_green_max, iend_max
   integer (kind=IntKind) :: MaxPrintLevel
!
   real (kind=RealKind), allocatable :: Position(:,:)
!
   type (NeighborStruct), pointer :: Neighbor
!
   real (kind=RealKind) :: Kvec(3)
!
   type MSTStruct
      integer :: lmax
      integer :: iend
      complex (kind=CmplxKind), pointer :: dos(:,:)
      complex (kind=CmplxKind), pointer :: green(:,:,:,:)      ! Stores the multiple scattering component of the Green fucntion, and
                                                               ! the single site scattering term may be included.
      complex (kind=CmplxKind), pointer :: der_green(:,:,:,:)  ! Stores the multiple scattering component of the Green fucntion derivative,
   end type MSTStruct                                          ! and the single site scattering term may be included.
!
   type (MSTStruct), allocatable :: mst(:)
   complex (kind=CmplxKind), allocatable, target :: wspace(:), wspacep(:), gspace(:)
   complex (kind=CmplxKind), allocatable, target :: dwspace(:), dgspace(:)
   complex (kind=CmplxKind), allocatable, target :: gspacep(:)
!
   complex (kind=CmplxKind), pointer :: gaunt(:,:,:)
   complex (kind=CmplxKind), allocatable :: store_space(:)
!
   integer (kind=IntKind), parameter :: method = 0 ! The method for performing the multi-L summations
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, kGID
   integer (kind=IntKind) :: MSGF_Form   ! = 0, Green Function = Z*(Tau-t)*Z
                                         ! = 1, Green Function = Z*Tau*Z - Z*J
                                         ! = 2, Green Function = Z*Tau*Z
   integer (kind=IntKind) :: MSDOS_Form  ! = 0, Green Function = Z*(Tau-t)*Z
                                         ! = 1, Green Function = Z*Tau*Z - Z*J
!
   type IDOSStruct
      real (kind=RealKind), pointer :: msdos_cell(:,:)
      real (kind=RealKind), pointer :: msdos_mt(:,:)
      real (kind=RealKind), pointer :: mspdos_cell(:,:,:)
      real (kind=RealKind), pointer :: mspdos_mt(:,:,:)
   end type IDOSStruct
!
   type (IDOSStruct), allocatable :: SpaceIntegratedDOS(:)
!
   complex (kind=CmplxKind), allocatable, target :: gfws_comp(:)
   complex (kind=CmplxKind), allocatable, target :: dosws_comp(:)
!
   logical :: isDosSymmOn = .false.
   logical :: rad_deriv = .false.
!   
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMSSolver( num_latoms, index, lmaxkkr, lmaxphi, lmaxgreen, &
                            local_posi, pola, cant, rel, istop, iprint,     &
                            derivative )
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
   use GauntFactorsModule, only : getK3, getNumK3, getGauntFactor
!
   use ScfDataModule, only : isKKR, isScreenKKR, isLSMS, isKKRCPA, isEmbeddedCluster
   use ScfDataModule, only : isChargeSymm
   use ScfDataModule, only : retrieveEffectiveMediumParams
!
   use ClusterMatrixModule, only : initClusterMatrix
!
   use CrystalMatrixModule, only : initCrystalMatrix
!
   use CPAMediumModule, only : initCPAMedium
!
   implicit none
!
   character (len=22), parameter :: sname = 'initMSSolver'
!
   character (len=*), intent(in) :: istop
!
   logical, intent(in), optional :: derivative
!
   integer (kind=IntKind), intent(in) :: num_latoms
   integer (kind=IntKind), intent(in) :: index(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxkkr(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxphi(num_latoms)
   integer (kind=IntKind), intent(in) :: lmaxgreen(num_latoms)
   integer (kind=IntKind), intent(in) :: pola
   integer (kind=IntKind), intent(in) :: cant
   integer (kind=IntKind), intent(in) :: rel
   integer (kind=IntKind), intent(in) :: iprint(num_latoms)
!
   real (kind=RealKind), intent(in) :: local_posi(3,num_latoms)
   real (kind=RealKind) :: em_mix_0, em_mix_1, em_eswitch, em_tol
!
   integer (kind=IntKind) :: lmax, i
   integer (kind=IntKind) :: klp1, klp2, i3, klg
   integer (kind=IntKind) :: em_mix_type, em_max_iter
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!
   if (present(derivative)) then
      rad_deriv = derivative
   else
      rad_deriv = .false.
   endif
!
!  -------------------------------------------------------------------
   call initParameters(num_latoms,lmaxkkr,lmaxphi,lmaxgreen,pola,cant,rel,istop,iprint)
!  -------------------------------------------------------------------
!
   lmax = 0
   do i=1,num_latoms
      lmax = max( lmax, lmaxkkr(i), lmaxphi(i), lmaxgreen(i) )
   enddo
!  -------------------------------------------------------------------
   call initGauntFactors(lmax,istop,-1)
!  -------------------------------------------------------------------
!
   if ( isLSMS() ) then
!     ----------------------------------------------------------------
      call initClusterMatrix(num_latoms,index,lmaxkkr,lmaxphi,local_posi, &
                             cant,rel,istop,iprint)
!     ----------------------------------------------------------------
   else if ( isScreenKKR() ) then
      call ErrorHandler('initMSSolver','Screen KKR is not implemented yet.')
!     ----------------------------------------------------------------
!     call initTauScreenKKR( bravais, LocalNumAtoms, cant, pola, local_posi, &
!                            lmaxkkr, rel, iprint, istop )
!     ----------------------------------------------------------------
   else if ( isKKR() ) then
!     ----------------------------------------------------------------
      call initCrystalMatrix( LocalNumAtoms, cant, lmaxkkr, rel, istop, iprint)
!     ----------------------------------------------------------------
   else if ( isKKRCPA() ) then
      call retrieveEffectiveMediumParams(mix_type = em_mix_type,      &
                                         max_iter = em_max_iter,      &
                                         alpha_0  = em_mix_0,         &
                                         alpha_1  = em_mix_1,         &
                                         eSwitch  = em_eswitch, tol = em_tol)
!     ----------------------------------------------------------------
      call initCPAMedium(cant=cant, lmax_kkr=lmaxkkr, rel=rel,        &
                         cpa_mix_type=em_mix_type,                    &
                         cpa_max_iter=em_max_iter,                    &
                         cpa_mix_0=em_mix_0, cpa_mix_1=em_mix_1,      &
                         cpa_eswitch=em_eswitch, cpa_tol=em_tol,      &
                         istop=istop, iprint=iprint)
!     ----------------------------------------------------------------
   else if (isEmbeddedCluster()) then
!     ----------------------------------------------------------------
      call retrieveEffectiveMediumParams(mix_type = em_mix_type,      &
                                         max_iter = em_max_iter,      &
                                         alpha_0  = em_mix_0,         &
                                         alpha_1  = em_mix_1,         &
                                         eSwitch  = em_eswitch, tol = em_tol)
!     ----------------------------------------------------------------
      call initClusterMatrix(num_latoms,index,lmaxkkr,lmaxphi,local_posi, &
                             cant,rel,istop,iprint)
!     ----------------------------------------------------------------
      call initCPAMedium(cant=cant, lmax_kkr=lmaxkkr, rel=rel,        &
                         cpa_mix_type=em_mix_type,                    &
                         cpa_max_iter=em_max_iter,                    &
                         cpa_mix_0=em_mix_0, cpa_mix_1=em_mix_1,      &
                         cpa_eswitch=em_eswitch, cpa_tol=em_tol,      &
                         istop=istop, iprint=iprint)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initMSSolver','Unknown MST calculation category')
!     ----------------------------------------------------------------
   endif
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
   gaunt => aliasArray3_c(gspacep,kmax_phi_max,kmax_green_max,kmax_phi_max)
   gaunt = CZERO
   if (method == 0 .or. method == 2) then
      do klp2 = 1, kmax_phi_max
         do klg = 1, kmax_green_max
            do klp1 = 1, kmax_phi_max
               do i3 = 1, nj3(klp1,klg)
                  if (kj3(i3,klp1,klg) == klp2) then
                     gaunt(klp1,klg,klp2) = cgnt(i3,klp1,klg)
                  endif
               enddo
            enddo
         enddo
      enddo
   else  ! Store gaunt differently to help speeding up the data access...
      do klp1 = 1, kmax_phi_max
         do klg = 1, kmax_green_max
            do klp2 = 1, kmax_phi_max
               do i3 = 1, nj3(klp1,klg)
                  if (kj3(i3,klp1,klg) == klp2) then
                     gaunt(klp2,klg,klp1) = cgnt(i3,klp1,klg)
                  endif
               enddo
            enddo
         enddo
      enddo
   endif
!
   nullify(nj3, kj3, cgnt)
!
!  -------------------------------------------------------------------
   call endGauntFactors()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Use the existing "K-Mesh" MPI group to create a parallelization
!  over the kl loop in the single site solver...
!  -------------------------------------------------------------------
   kGID = getGroupID('K-Mesh')
   NumPEsInGroup = getNumPEsInGroup(kGID)
   MyPEinGroup = getMyPEinGroup(kGID)
!  -------------------------------------------------------------------
!
   Initialized = .true.
   Energy = -10.0d0
   MSGF_Form = 0
   MSDOS_Form = 0
!
   isDosSymmOn = isChargeSymm()
!
   end subroutine initMSSolver
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initParameters(num_atoms, lmaxkkr, lmaxphi, lmaxgreen,  &
                             pola, cant, rel, istop, iprint)
!  ===================================================================
   use RadialGridModule, only : getNumRmesh, getMaxNumRmesh
!
   use AtomModule, only : getLocalNumSpecies
   implicit none
!
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: lmaxkkr(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxphi(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxgreen(num_atoms)
   integer (kind=IntKind), intent(in) :: pola
   integer (kind=IntKind), intent(in) :: cant
   integer (kind=IntKind), intent(in) :: rel
   integer (kind=IntKind), intent(in) :: iprint(num_atoms)
   integer (kind=LongIntKind) :: wspace_size, gspace_size
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind) :: i, lmax_max, jmax, iend, kmax, NumSpecies
!
   if (Initialized) then
!     ----------------------------------------------------------------
      call WarningHandler('initMSSolver',                   &
                  'MSSolverModule has already been initialized')
!     ----------------------------------------------------------------
      return
   else if (num_atoms < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMSSolver','num_atoms < 1',num_atoms)
!     ----------------------------------------------------------------
   else if (pola < 1 .or. pola > 2) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMSSolver',                     &
                        'Invalid spin polarization index',pola)
!     ----------------------------------------------------------------
   else if (cant < 1 .or. cant > 2) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMSSolver',                     &
                        'Invalid spin canting index',cant)
!     ----------------------------------------------------------------
   else if (pola < cant ) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMSSolver',                     &
                        'Polarization = 1, and Canting = 2')
!     ----------------------------------------------------------------
   else if (rel < 0 .or. rel > 2) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMSSolver','rel < 0 or rel > 2', rel)
!     ----------------------------------------------------------------
   endif
!
   LocalNumAtoms = num_atoms
   n_spin_pola = pola
   n_spin_cant = cant
   Relativity = rel
   stop_routine = istop
!
   allocate( print_instruction(LocalNumAtoms) )
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
   call genFactors(lmax_max)
!  -------------------------------------------------------------------
   MaxPrintLevel = maxval(print_instruction(1:LocalNumAtoms))
!
   allocate( mst(LocalNumAtoms) )
   do i=1, LocalNumAtoms
      kmax = (lmaxgreen(i)+1)**2
      iend = getNumRmesh(i)
      mst(i)%lmax = lmaxgreen(i)
      mst(i)%iend = iend
      NumSpecies = getLocalNumSpecies(i)
      allocate( mst(i)%dos(n_spin_cant*n_spin_cant,NumSpecies) )
      allocate( mst(i)%green(iend,kmax,n_spin_cant*n_spin_cant,NumSpecies) )
      if (rad_deriv) then
         allocate( mst(i)%der_green(iend,kmax,n_spin_cant*n_spin_cant,NumSpecies) )
      endif
   enddo
   iend_max = getMaxNumRmesh()
   wspace_size = iend_max*kmax_phi_max*kmax_kkr_max
   gspace_size = iend_max*kmax_green_max*kmax_phi_max
   allocate( wspace(wspace_size), wspacep(kmax_phi_max*kmax_green_max) )
   allocate( gspace(gspace_size), gspacep(kmax_phi_max*kmax_green_max*kmax_phi_max) )
   if (rad_deriv) then
      allocate( dgspace(gspace_size), dwspace(wspace_size) )
   endif
!
   end subroutine initParameters
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMSSolver()
!  ===================================================================
   use ScfDataModule, only : isKKR, isScreenKKR, isLSMS, isKKRCPA,    &
                             isEmbeddedCluster
!
   use ClusterMatrixModule, only : endClusterMatrix
!
   use CrystalMatrixModule, only : endCrystalMatrix
!
   use CPAMediumModule, only : endCPAMedium
!
   implicit none
!
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endMSSolver',                      &
                        'MSSolverModule has not been initialized')
!     ----------------------------------------------------------------
   endif
!
   deallocate( print_instruction )
   deallocate( lmax_kkr )
   deallocate( lmax_phi )
   deallocate( kmax_kkr )
   deallocate( kmax_phi )
   do i=1, LocalNumAtoms
      deallocate( mst(i)%dos )
      deallocate( mst(i)%green )
      if (rad_deriv) then
         deallocate( mst(i)%der_green )
      endif
   enddo
   deallocate( mst, wspace, wspacep, gspace, gspacep )
   if (rad_deriv) then
      deallocate( dwspace, dgspace )
   endif
   nullify( Neighbor )
!
   if ( isLSMS() ) then
!     ----------------------------------------------------------------
      call endClusterMatrix()
!     ----------------------------------------------------------------
   else if (isScreenKKR()) then
!     ----------------------------------------------------------------
!     call endTauScreenKKR()
!     ----------------------------------------------------------------
   else if ( isKKR() ) then
!     ----------------------------------------------------------------
      call endCrystalMatrix()
!     ----------------------------------------------------------------
   else if (isKKRCPA()) then
!     ----------------------------------------------------------------
      call endCPAMedium()
!     ----------------------------------------------------------------
   else if (isEmbeddedCluster()) then
!     ----------------------------------------------------------------
      call endCPAMedium()
!     ----------------------------------------------------------------
      call endClusterMatrix()
!     ----------------------------------------------------------------
   endif
!
   if (allocated(store_space)) then
      deallocate(store_space)
   endif
!
   deallocate( lofk, mofk, jofk, m1m, lofj, mofj )
!
   if (allocated(SpaceIntegratedDOS) ) then
      do i=1, LocalNumAtoms
         deallocate( SpaceIntegratedDOS(i)%msdos_cell )
         deallocate( SpaceIntegratedDOS(i)%msdos_mt )
         deallocate( SpaceIntegratedDOS(i)%mspdos_cell )
         deallocate( SpaceIntegratedDOS(i)%mspdos_mt )
      enddo
      deallocate( gfws_comp, dosws_comp, SpaceIntegratedDOS )
   endif
!
   Initialized = .false.
   Energy = CZERO
   isDosSymmOn = .false.
!
   end subroutine endMSSolver
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDOS_is(is,id,ia,mst_term_only) result(dos)
!  ===================================================================
   use AtomModule, only : getLocalNumSpecies
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
!
   logical, intent(out), optional :: mst_term_only
!
   complex (kind=CmplxKind) :: dos
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDOS_is','Invalid site index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('getDOS_is','Invalid spin index',is)
   endif
   if (ia < 1 .or. ia > getLocalNumSpecies(id)) then
      call ErrorHandler('getDOS_is','Invalid species index',ia)
   endif
!
   dos = mst(id)%dos(is,ia)
   if (present(mst_term_only)) then
      if (MSDOS_Form == 0) then
         mst_term_only = .true.
      else
         mst_term_only = .false.
      endif
   endif
!
   end function getDOS_is
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDOS_sc(id,ia,mst_term_only) result(dos)
!  ===================================================================
   use AtomModule, only : getLocalNumSpecies
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   logical, intent(out), optional :: mst_term_only
!
   complex (kind=CmplxKind) :: dos(n_spin_cant*n_spin_cant)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDOS_is','Invalid site index',id)
   endif
   if (ia < 1 .or. ia > getLocalNumSpecies(id)) then
      call ErrorHandler('getDOS_is','Invalid species index',ia)
   endif
!
   dos = mst(id)%dos(:,ia)
   if (present(mst_term_only)) then
      if (MSDOS_Form == 0) then
         mst_term_only = .true.
      else
         mst_term_only = .false.
      endif
   endif
!
   end function getDOS_sc
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMSGreenFunction(id,gform) result(green)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out), optional :: gform
!
   complex (kind=CmplxKind), pointer :: green(:,:,:,:)
!
   green => mst(id)%green
   if (present(gform)) then
      gform = MSGF_Form
   endif
!
   end function getMSGreenFunction
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMSGreenFunctionDerivative(id,gform) result(der_green)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out), optional :: gform
!
   complex (kind=CmplxKind), pointer :: der_green(:,:,:,:)
!
   der_green => mst(id)%der_green
   if (present(gform)) then
      gform = MSGF_Form
   endif
!
   end function getMSGreenFunctionDerivative
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMSGreenMatrix(id) result(mat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   complex (kind=CmplxKind), pointer :: mat(:,:,:,:)
!
   nullify(mat)
!
   call ErrorHandler('getMSGreenMatrix','Density matrix calculation is not implemented yet')
!
   end function getMSGreenMatrix
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeMSTMatrix(is,e)
!  ===================================================================
   use ScfDataModule, only : isLSMS, isScreenKKR, isKKRCPA, isKKR, isEmbeddedCluster
!
   use SSSolverModule, only : getScatteringMatrix
!
   use ClusterMatrixModule, only : calClusterMatrix
!
   use CrystalMatrixModule, only : calCrystalMatrix
!
   use CPAMediumModule, only : computeCPAMedium
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   complex (kind=CmplxKind), intent(in) :: e
!
!  ===================================================================
!  call calClusterMatrix or calCrystalMatrix to calculate the TAU(0,0) matrix
!  NOTE: Needs to be checked for is = 2 and n_spin_cant = 2
!  ===================================================================
   if ( isLSMS() ) then
!     ----------------------------------------------------------------
      call calClusterMatrix(e,getScatteringMatrix)
!     ----------------------------------------------------------------
   else if (isScreenKKR()) then
!     ----------------------------------------------------------------
      call ErrorHandler('computeMSTMatrix','ScreenKKR is yet to be implemented')
!     call calScreenTauBZ( e, is )
!     ----------------------------------------------------------------
   else if (isKKR()) then
!     ----------------------------------------------------------------
      call calCrystalMatrix(e,getScatteringMatrix)
!     call calCrystalMatrix(e,getScatteringMatrix,tau_needed=.true.)
!     call calCrystalMatrix(e,getScatteringMatrix,use_tmat=.true.,tau_needed=.true.)
!     ----------------------------------------------------------------
   else if (isKKRCPA()) then
!     ----------------------------------------------------------------
      call computeCPAMedium(e)
!     ----------------------------------------------------------------
   else if (isEmbeddedCluster()) then  ! Needs further work.....
!     ----------------------------------------------------------------
      call calClusterMatrix(e,getScatteringMatrix)
!     ----------------------------------------------------------------
      call computeCPAMedium(e)
!     ----------------------------------------------------------------
   else
      call ErrorHandler('computeMSTMatrix','Undefined method')
   endif
!
   end subroutine computeMSTMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeMSGreenFunction(is, e, add_Ts, add_Gs, isSphSolver)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
!
   use GroupCommModule, only : GlobalSumInGroup
!
   use ScfDataModule, only : isKKR, isScreenKKR, isLSMS, isKKRCPA, isEmbeddedCluster
!
   use SystemSymmetryModule, only : getSymmetryFlags
!
   use AtomModule, only : getLocalNumSpecies
!
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
   use SSSolverModule, only : getRegSolutionDerivative
   use SSSolverModule, only : getSolutionFlags, getOmegaHatMatrix
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : computeGreenFunction, getGreenFunction
   use SSSolverModule, only : getGreenFunctionDerivative
!
   use ClusterMatrixModule, only : getClusterKau => getKau
!
   use CrystalMatrixModule, only : getCrystalKau => getKau
   use CrystalMatrixModule, only : getCrystalTau => getTau
!
   use CPAMediumModule, only : getImpurityMatrix, getCPAMatrix
!
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   logical, optional, intent(in) :: add_Ts, add_Gs, isSphSolver
   logical :: add_SingleSiteT, add_SingleSiteG
!
   integer (kind=IntKind), intent(in) :: is
!
   complex (kind=CmplxKind), intent(in) :: e
!
   integer (kind=IntKind) :: n, info, id, js1, js2, ns, kmaxk, kmaxp, kmaxg, irmax
   integer (kind=IntKind) :: klg, kl1, kl2, klp1, klp2, ir, kl2c, m2, np, ia
   integer (kind=IntKind), pointer :: green_flags(:)
!
   complex (kind=CmplxKind), pointer :: tfac(:,:), gfs(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr_right(:,:,:), PhiLr_left(:,:,:), kau00(:,:,:)
   complex (kind=CmplxKind), pointer :: der_PhiLr_right(:,:,:), der_PhiLr_left(:,:,:)
   complex (kind=CmplxKind), pointer :: gf(:,:), pp(:,:), ppr(:,:), ppg(:,:,:)
   complex (kind=CmplxKind), pointer :: dgf(:,:), dpp(:,:), dppr(:,:), dppg(:,:,:)
   complex (kind=CmplxKind), pointer :: pau00(:,:), OmegaHat(:,:), p_kau00(:,:)
   complex (kind=CmplxKind) :: cfac, kappa
!
   character (len=20), parameter :: sname = 'computeMSGreenFunction'
!
!#ifdef TIMING
   real (kind=RealKind) :: time
!
   time = getTime()
!#endif
   if (.not.Initialized) then
      call ErrorHandler(sname,'module not initialized')
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler(sname,'invalid spin index',is)
   endif
! 
   if (present(add_Gs)) then
      add_SingleSiteG = add_Gs
   else
      add_SingleSiteG = .false.
   endif
!
   if (present(add_Ts)) then
      add_SingleSiteT = add_Ts
   else
      add_SingleSiteT = .false.
   endif
!
   if (add_SingleSiteG) then
      MSGF_Form = 1
   else if (add_SingleSiteT) then
      MSGF_Form = 2
   else
      MSGF_Form = 0
   endif
!
   if (add_SingleSiteG .or. add_SingleSiteT) then
!     ================================================================
!     This is the case when [Kau00 + kappa*OmegaHat] is used in the 
!     calculation of the Green function and the DOS.
!     ================================================================
      if (.not.allocated(store_space)) then
         allocate(store_space(kmax_kkr_max*kmax_kkr_max))
      endif
   endif
!
   if (add_SingleSiteG) then
      do id = 1, LocalNumAtoms
         do js1 = 1, n_spin_cant
            ns = max(js1,is)
            if (present(isSphSolver)) then
!              -------------------------------------------------------
               call solveSingleScattering(spin=ns,site=id,e=e,        &
                                          vshift=CZERO,isSphSolver=isSphSolver,useIrrSol='H')
!              -------------------------------------------------------
            else
!              -------------------------------------------------------
               call solveSingleScattering(spin=ns,site=id,e=e,        &
                                          vshift=CZERO,useIrrSol='H')
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            call computeGreenFunction(ns,id)
!           ----------------------------------------------------------
         enddo
      enddo
   else
      do id = 1, LocalNumAtoms
         do js1 = 1, n_spin_cant
            ns = max(js1,is)
!           ----------------------------------------------------------
            call solveSingleScattering(spin=ns,site=id,e=e,vshift=CZERO)
!           ----------------------------------------------------------
         enddo
      enddo
   endif
!
!  -------------------------------------------------------------------
   call computeMSTMatrix(is,e)
!  -------------------------------------------------------------------
!
#ifdef TIMING
   if (MaxPrintLevel >= 0) then
      write(6,*)'MSSolver:: Total Time  : ', getTime()-time
   endif
#endif
!
   Energy = e
   kappa = sqrt(e)
!
   do id = 1, LocalNumAtoms
      kmaxk = kmax_kkr(id)
      kmaxp = kmax_phi(id)
      kmaxg = (mst(id)%lmax+1)**2
!     ================================================================
      irmax = getSolutionRmeshSize(id)
      if (irmax < mst(id)%iend) then
         call ErrorHandler(sname,'Source code error: irmax < iend',irmax,mst(id)%iend)
      endif
!     pp => aliasArray2_c(wspace,mst(id)%iend,kmaxp)
      pp => aliasArray2_c(wspace,irmax,kmaxp)
!     ppr => aliasArray2_c(wspace,mst(id)%iend*kmaxp,kmaxk)
      ppr => aliasArray2_c(wspace,irmax*kmaxp,kmaxk)
!     ppg => aliasArray3_c(gspace,mst(id)%iend,kmaxg,kmaxp)
      ppg => aliasArray3_c(gspace,irmax,kmaxg,kmaxp)
      tfac => aliasArray2_c(wspacep,kmaxp,kmaxg)
      if (rad_deriv) then
         dpp => aliasArray2_c(dwspace,irmax,kmaxp)
         dppr => aliasArray2_c(dwspace,irmax*kmaxp,kmaxk)
         dppg => aliasArray3_c(dgspace,irmax,kmaxg,kmaxp)
      endif
      if (add_SingleSiteT) then
         pau00 => aliasArray2_c(store_space,kmaxk,kmaxk)
      endif
      do ia = 1, getLocalNumSpecies(id)
         if (isLSMS()) then
            kau00 => getClusterKau(local_id=id) ! Kau00 = energy * S^{-1} * [Tau00 - t_matrix] * S^{-T*}
         else if (isKKR()) then
!kau00=>getCrystalTau(local_id=id)
!call writeMatrix('Tau',kau00(:,:,1),kmaxk,kmaxk,TEN2m6)
            kau00 => getCrystalKau(local_id=id) ! Kau00 = energy * S^{-1} * [Tau00 - t_matrix] * S^{-T*}
         else if (isKKRCPA() .or. isEmbeddedCluster()) then
!p_kau00 => getCPAMatrix('Tcpa',site=id)
!call writeMatrix('Tcpa',p_kau00,kmaxk,kmaxk,TEN2m6)
!kau00=>getImpurityMatrix('Tau_a',site=id,atom=0)
!call writeMatrix('Tau_a',kau00(:,:,1),kmaxk,kmaxk,TEN2m6)
            kau00 => getImpurityMatrix('Kau_a',site=id,atom=ia) ! Kau00 = energy * S^{-1} * [Tau_a - t_matrix] * S^{-T*}
         endif
!        -------------------------------------------------------------
!        call writeMatrix('Kau_a',kau00(:,:,1),kmaxk,kmaxk,TEN2m6)
!        -------------------------------------------------------------
         ns = 0
         do js2 = 1, n_spin_cant
!           ==========================================================
!           If needed, add kappa*Omega to Kau00, which is equivalent to add
!           t_mat to [tau00 - t_mat]
!           ==========================================================
            if (add_SingleSiteT) then
!              -------------------------------------------------------
               OmegaHat => getOmegaHatMatrix(js2,site=id,atom=ia)
!              -------------------------------------------------------
            endif
            PhiLr_right => getRegSolution(js2,site=id,atom=ia)
            if (rad_deriv) then
               der_PhiLr_right => getRegSolutionDerivative(js2,site=id,atom=ia)
            endif
            do js1 = 1, n_spin_cant
               PhiLr_left => getRegSolution(js1,site=id,atom=ia)
               if (rad_deriv) then
                  der_PhiLr_left => getRegSolutionDerivative(js1,site=id,atom=ia)
               endif
               ns = ns + 1
               gf => mst(id)%green(:,:,ns,ia)
               gf = CZERO
               if (rad_deriv) then
                  dgf => mst(id)%der_green(:,:,ns,ia)
                  dgf = CZERO
               endif
               p_kau00 => kau00(:,:,ns)
!              =======================================================
!              gf is the multiple scattering part of the Green function
!              multiplied by r^2:
!                  gf = Z_L*(Tau00-t_matrix)*Z_L^{*}*r^2
!                     = Phi_L*Kau00*Phi_L^{*}*r^2
!              Here implements three different methods for checking against each other
!
!              If needed, add kappa*Omega to Kau00, which is equivalent to add
!              t_mat to [tau00 - t_mat]
!              =======================================================
               if (add_SingleSiteT .and. js1 == js2) then
!                 ----------------------------------------------------
                  call zcopy(kmaxk*kmaxk,p_kau00,1,pau00,1)
                  call zaxpy(kmaxk*kmaxk,kappa,OmegaHat,1,pau00,1)
!                 ----------------------------------------------------
                  p_kau00 => pau00
               endif
!              =======================================================
               if (method == 0) then
!                 ====================================================
!                 ppr(ir,klp1,kl2) = sum_kl1 PhiLr_left(ir,klp1,kl1) * 
!                                            p_kau00(kl1,kl2)
!                 ----------------------------------------------------
                  call zgemm('n','n',irmax*kmaxp,kmaxk,kmaxk,CONE,PhiLr_left, &
                             irmax*kmaxp,p_kau00,kmaxk,CZERO,ppr,irmax*kmaxp)
!                 ----------------------------------------------------
                  if (rad_deriv) then
!                    -------------------------------------------------
                     call zgemm('n','n',irmax*kmaxp,kmaxk,kmaxk,CONE,der_PhiLr_left, &
                                 irmax*kmaxp,p_kau00,kmaxk,CZERO,dppr,irmax*kmaxp)
!                    -------------------------------------------------
                  endif
                  np = mod(kmaxk,NumPEsInGroup)
                  do kl2 = MyPEinGroup+1, kmaxk-np, NumPEsInGroup
                     m2 = mofk(kl2)
                     kl2c = kl2 -2*m2
                     cfac = m1m(m2)
!                    =================================================
!                    ppg(ir,klg,klp2;kl2) = sum_klp1 (-1)^m2 * ppr(ir,klp1,kl2c) *
!                                                    gaunt(klp1,klg,klp2)
!                    -------------------------------------------------
                     call zgemm('n','n',irmax,kmaxg*kmaxp,kmaxp,cfac,ppr(1,kl2c), &
                                 irmax,gaunt,kmaxp,CZERO,ppg,irmax)
!                    -------------------------------------------------
                     if (rad_deriv) then
!                       ----------------------------------------------
                        call zgemm('n','n',irmax,kmaxg*kmaxp,kmaxp,cfac, &
                                   dppr(1,kl2c),irmax,gaunt,kmaxp,CZERO,dppg,irmax)
!                       ----------------------------------------------
                     endif
!
!                    =================================================
!                    gf(ir,klg) = sum_{kl2,klp2} ppg(ir,klg,klp2;kl2) * 
!                                                PhiLr_right(ir,klp2,kl2)
!                    =================================================
                     do klp2 = 1, kmaxp
                        do klg = 1, kmaxg
                           do ir = 1, mst(id)%iend
                              gf(ir,klg) =  gf(ir,klg) + ppg(ir,klg,klp2)*PhiLr_right(ir,klp2,kl2)
                           enddo
                           if (rad_deriv) then
                              do ir = 1, mst(id)%iend
                                 dgf(ir,klg) =  dgf(ir,klg) + dppg(ir,klg,klp2)*PhiLr_right(ir,klp2,kl2)  &
                                                            + ppg(ir,klg,klp2)*der_PhiLr_right(ir,klp2,kl2)
                              enddo
                           endif
                        enddo
                     enddo
                  enddo ! kl2
                  if (NumPEsInGroup > 1) then
!                    -------------------------------------------------
                     call GlobalSumInGroup(kGID,gf,mst(id)%iend,kmaxg)
!                    -------------------------------------------------
                     if (rad_deriv) then
!                       ----------------------------------------------
                        call GlobalSumInGroup(kGID,dgf,mst(id)%iend,kmaxg)
!                       ----------------------------------------------
                     endif
                  endif
                  do kl2 = kmaxk-np+1,kmaxk
                     m2 = mofk(kl2)
                     kl2c = kl2 -2*m2
                     cfac = m1m(m2)
!                    =================================================
!                    ppg(ir,klg,klp2;kl2) = sum_klp1 (-1)^m2 * ppr(ir,klp1,kl2c) * gaunt(klp1,klg,klp2)
!                    -------------------------------------------------
                     call zgemm('n','n',irmax,kmaxg*kmaxp,kmaxp,cfac,ppr(1,kl2c),irmax,gaunt,kmaxp,CZERO,ppg,irmax)
!                    -------------------------------------------------
                     if (rad_deriv) then
!                       ----------------------------------------------
                        call zgemm('n','n',irmax,kmaxg*kmaxp,kmaxp,cfac,dppr(1,kl2c),irmax,gaunt,kmaxp,CZERO, &
                                   dppg,irmax)
!                       ----------------------------------------------
                     endif
!
!                    =================================================
!                    gf(ir,klg) = sum_{kl2,klp2} ppg(ir,klg,klp2;kl2) * PhiLr_right(ir,klp2,kl2)
!                    =================================================
                     do klp2 = 1, kmaxp
                        do klg = 1, kmaxg
                           do ir = 1, mst(id)%iend
                              gf(ir,klg) =  gf(ir,klg) + ppg(ir,klg,klp2)*PhiLr_right(ir,klp2,kl2)
                           enddo
                           if (rad_deriv) then
                              do ir = 1, mst(id)%iend
                                 dgf(ir,klg) =  dgf(ir,klg) + dppg(ir,klg,klp2)*PhiLr_right(ir,klp2,kl2)  &
                                                            + ppg(ir,klg,klp2)*der_PhiLr_right(ir,klp2,kl2)
                              enddo
                           endif
                        enddo
                     enddo
                  enddo
               else if (method == 1) then
                  do kl1 = kmaxk,1,-1
                     do klp1 = kmaxp,1,-1
                        do kl2 = kmaxk,1,-1
                           m2 = mofk(kl2)
                           kl2c = kl2 -2*m2
!                          ===========================================
!                          pp(ir,klp2;kl2,klp1,kl1) =
!                              PhiLr_left(ir,klp1,kl1) * PhiLr_right(ir,klp2,kl2)
!                          ===========================================
                           pp = CZERO
                           do klp2 = 1, kmaxp
                              do ir = 1, irmax
                                 pp(ir,klp2) = PhiLr_left(ir,klp1,kl1)*PhiLr_right(ir,klp2,kl2)
                              enddo
                           enddo
                           if (rad_deriv) then
                              dpp = CZERO
                              do klp2 = 1, kmaxp
                                 do ir = 1, irmax
                                    dpp(ir,klp2) = der_PhiLr_left(ir,klp1,kl1)*PhiLr_right(ir,klp2,kl2) &
                                                 + PhiLr_left(ir,klp1,kl1)*der_PhiLr_right(ir,klp2,kl2)
                                 enddo
                              enddo
                           endif
!                          ===========================================
!                          tfac(klp2,klg;kl2,klp1,kl1) = gaunt(klp2,klg,klp1)*p_kau00(kl1,kl2)
!                          ===========================================
                           do klg = 1, kmaxg
                              do klp2 = 1, kmaxp
                                 tfac(klp2,klg) = m1m(m2)*gaunt(klp2,klg,klp1)*p_kau00(kl1,kl2c)
                              enddo
                           enddo
!                          ===========================================
!                          gf(ir,klg) = sum_{kl1,klp1,kl2,klp2} pp(ir,klp2;kl2,klp1,kl1) * tfac(klp2,klg;kl2,klp1,kl1)
!                          -------------------------------------------
                           call zgemm('n','n',mst(id)%iend,kmaxg,kmaxp,CONE,pp,irmax,tfac,kmaxp,CONE,gf,mst(id)%iend)
!                          -------------------------------------------
                           if (rad_deriv) then
!                             ----------------------------------------
                              call zgemm('n','n',mst(id)%iend,kmaxg,kmaxp,CONE,dpp,irmax,tfac,kmaxp,CONE,dgf,mst(id)%iend)
!                             ----------------------------------------
                           endif
                        enddo ! kl2
                     enddo ! do klp1
                  enddo ! do kl1
               else
                  do kl2 = 1, kmaxk
                     m2 = mofk(kl2)
                     kl2c = kl2 -2*m2
                     cfac = m1m(m2)
!                    -------------------------------------------------
                     call zgemv('n',irmax*kmaxp,kmaxk,cfac,PhiLr_left,irmax*kmaxp,p_kau00(1,kl2c),1,CZERO,pp,1)
!                    -------------------------------------------------
                     if (rad_deriv) then
!                       ----------------------------------------------
                        call zgemv('n',irmax*kmaxp,kmaxk,cfac,der_PhiLr_left,irmax*kmaxp,p_kau00(1,kl2c),1,CZERO, &
                                   dpp,1)
!                       ----------------------------------------------
                     endif
                     do klp2 = 1, kmaxp
                        do klg = 1, kmaxg
                           do klp1 = 1, kmaxp
                              do ir = 1, mst(id)%iend
                                 gf(ir,klg) = gf(ir,klg) + gaunt(klp1,klg,klp2)*pp(ir,klp1)*PhiLr_right(ir,klp2,kl2)
                              enddo
                              if (rad_deriv) then
                                 do ir = 1, mst(id)%iend
                                    dgf(ir,klg) = dgf(ir,klg) +                                                   &
                                                  gaunt(klp1,klg,klp2)*( dpp(ir,klp1)*PhiLr_right(ir,klp2,kl2) +  &
                                                                         pp(ir,klp1)*der_PhiLr_right(ir,klp2,kl2) )
                                                             
                                 enddo
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               endif
               if (add_SingleSiteG .and. js1 == js2) then  ! This "js1==js2" logic
                                                           ! needs to be checked for spin-canted case
!                 ====================================================
!                 Add the singe site Green function to gf = Z*(tau-t)*Z, 
!                 so that gf = Z*tau*Z - Z*J
!                 ====================================================
                  gfs => getGreenFunction(spin=max(js1,is),site=id,atom=ia)
                  if (size(gfs,1) < size(gf,1)) then
                     call ErrorHandler('computeMSGreenFunction',                                  &
                                       '1st dim. of ss < 1st dim of ms green function arrays',    &
                                       size(gfs,1), size(gf,1))
                  else if (size(gfs,2) < size(gf,2)) then
                     call ErrorHandler('computeMSGreenFunction',                                  &
                                       '2nd dim. of ss < 2nd dim of ms green function arrays',    &
                                       size(gfs,2), size(gf,2))
                  endif
                  do klg = 1, kmaxg
                     do ir = 1, mst(id)%iend
                        gf(ir,klg) = gf(ir,klg) + gfs(ir,klg)
                     enddo
                  enddo
                  if (rad_deriv) then
                     gfs => getGreenFunctionDerivative(spin=max(js1,is),site=id,atom=ia)
                     do klg = 1, kmaxg
                        do ir = 1, mst(id)%iend
                           dgf(ir,klg) = dgf(ir,klg) + gfs(ir,klg)
                        enddo
                     enddo
                  endif
               endif
!
!              =======================================================
!              Symmetrizing the green function, if needed.  Added by Yang on 09-21-2018
!              =======================================================
               if (isDosSymmOn) then
                  green_flags => getSymmetryFlags(id)
                  do klg = 1, kmaxg
                     if (green_flags(jofk(klg)) == 0) then
                        gf(:,klg) = CZERO
                        if (rad_deriv) then
                           dgf(:,klg) = CZERO
                        endif
                     endif
                  enddo
               endif
!              =======================================================
            enddo ! do js1
         enddo ! do js2
      enddo ! do ia
   enddo ! do id
!
   nullify(kau00, p_kau00, pau00, OmegaHat, gf, dgf, pp, ppr, ppg, tfac)
   nullify(PhiLr_right, PhiLr_left, der_PhiLr_right, der_PhiLr_left)
!
   end subroutine computeMSGreenFunction
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeMSPDOS(is, e, add_Gs, isSphSolver)
!  ===================================================================
   use MathParamModule, only : PI2
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use MPPModule, only : MyPE, syncAllPEs
   use GroupCommModule, only : GlobalSumInGroup
!
   use ScfDataModule, only : isKKR, isScreenKKR, isLSMS, isKKRCPA, isEmbeddedCluster
!
   use SystemSymmetryModule, only : getSymmetryFlags
!
   use AtomModule, only : getLocalNumSpecies
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RadialGridModule, only : getGrid
!
   use ClusterMatrixModule, only : getClusterKau => getKau
!
   use CrystalMatrixModule, only : getCrystalKau => getKau
!
   use CPAMediumModule, only : getImpurityMatrix
!
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
   use SSSolverModule, only : getSolutionFlags, getOmegaHatMatrix
   use SSSolverModule, only : solveSingleScattering
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   complex (kind=CmplxKind), intent(in) :: e
!
   logical, optional, intent(in) :: add_Gs, isSphSolver
   logical :: add_SingleSite
!
   integer (kind=IntKind) :: id, ns, kmaxk, kmaxp, kmaxg, jmaxg, irmax, jmax_green_max
   integer (kind=IntKind) :: n, info, js1, js2, jlg, lg, mg, js, ia
   integer (kind=IntKind) :: klg, klgc, kl, klp1, klp2, ir, klc, m, np
   integer (kind=IntKind) :: num_species
   integer (kind=IntKind), pointer :: dos_flags(:)
!
   real (kind=RealKind) :: dos_buf(kmax_phi_max,2)
!
   complex (kind=CmplxKind), pointer :: PhiLr_right(:,:,:), PhiLr_left(:,:,:), kau00(:,:,:)
   complex (kind=CmplxKind), pointer :: gf(:,:), dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: pp(:,:), ppr(:,:), ppg(:,:,:)
   complex (kind=CmplxKind), pointer :: pau00(:,:), OmegaHat(:,:), p_kau00(:,:)
   complex (kind=CmplxKind) :: cfac, kappa
!
   type (GridStruct), pointer :: Grid
!
   character (len=20), parameter :: sname = 'computeMSPDOS'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'module not initialized')
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler(sname,'invalid spin index',is)
   else if (present(add_Gs)) then
      add_SingleSite = add_Gs
   else
      add_SingleSite = .false.
   endif
!
   if (add_SingleSite) then
      MSDOS_Form = 1
!     ================================================================
!     This is the case when [Kau00 + kappa*OmegaHat] is used in the 
!     calculation of the Green function and the DOS.
!     ================================================================
      if (.not.allocated(store_space)) then
         allocate(store_space(kmax_kkr_max*kmax_kkr_max))
      endif
!
      do id = 1, LocalNumAtoms
         do js1 = 1, n_spin_cant
            ns = max(js1,is)
            if (present(isSphSolver)) then
!              -------------------------------------------------------
               call solveSingleScattering(spin=ns,site=id,e=e,vshift=CZERO, &
                                          isSphSolver=isSphSolver,useIrrSol='H')
!              -------------------------------------------------------
            else
!              -------------------------------------------------------
               call solveSingleScattering(spin=ns,site=id,e=e,vshift=CZERO,useIrrSol='H')
!              -------------------------------------------------------
            endif
         enddo
      enddo
   else
      MSDOS_Form = 0
      do id = 1, LocalNumAtoms
         do js1 = 1, n_spin_cant
            ns = max(js1,is)
!           ----------------------------------------------------------
            call solveSingleScattering(spin=ns,site=id,e=e,vshift=CZERO)
!           ----------------------------------------------------------
         enddo
      enddo
   endif
!
!  -------------------------------------------------------------------
   call computeMSTMatrix(is,e)
!  -------------------------------------------------------------------
!
   Energy = e
   kappa = sqrt(Energy)
!
   if (.not.allocated(SpaceIntegratedDOS)) then
      allocate( SpaceIntegratedDOS(LocalNumAtoms) )
      do id = 1, LocalNumAtoms
         num_species = getLocalNumSpecies(id)
         allocate( SpaceIntegratedDOS(id)%msdos_cell(n_spin_cant*n_spin_pola,num_species) )
         allocate( SpaceIntegratedDOS(id)%msdos_mt(n_spin_cant*n_spin_pola,num_species) )
         allocate( SpaceIntegratedDOS(id)%mspdos_cell(kmax_phi_max,n_spin_cant*n_spin_pola,num_species) )
         allocate( SpaceIntegratedDOS(id)%mspdos_mt(kmax_phi_max,n_spin_cant*n_spin_pola,num_species) )
      enddo
      jmax_green_max = (lmax_green_max+1)*(lmax_green_max+2)/2
      allocate( gfws_comp(iend_max*kmax_green_max), dosws_comp(iend_max*jmax_green_max) )
   endif
!
   do id = 1, LocalNumAtoms
      kmaxk = kmax_kkr(id)
      kmaxp = kmax_phi(id)
      kmaxg = (mst(id)%lmax+1)**2
      jmaxg = (mst(id)%lmax+1)*(mst(id)%lmax+2)/2
!     ================================================================
      irmax = getSolutionRmeshSize(id)
      if (irmax < mst(id)%iend) then
         call ErrorHandler(sname,'Source code error: irmax < iend',irmax,mst(id)%iend)
      endif
      pp => aliasArray2_c(wspace,irmax,kmaxp)
      ppr => aliasArray2_c(wspace,irmax*kmaxp,kmaxk)
      ppg => aliasArray3_c(gspace,irmax,kmaxg,kmaxp)
      gf => aliasArray2_c(gfws_comp,irmax,kmaxg)
      dos_r_jl => aliasArray2_c(dosws_comp,irmax,jmaxg)
      Grid => getGrid(id)
!
      if (add_SingleSite) then
         pau00 => aliasArray2_c(store_space,kmaxk,kmaxk)
      endif
      do ia = 1, getLocalNumSpecies(id)
         if (isLSMS()) then
            kau00 => getClusterKau(local_id=id) ! Kau00 = energy * S^{-1} * [Tau00 - t_matrix] * S^{-1*}
         else if (isKKR()) then
            kau00 => getCrystalKau(local_id=id) ! Kau00 = energy * S^{-1} * [Tau00 - t_matrix] * S^{-1*}
         else if ( isKKRCPA() .or. isEmbeddedCluster()) then
            kau00 => getImpurityMatrix('Kau_a',site=id,atom=ia)
         endif
         js = 0
         do js2 = 1, n_spin_cant
!           ==========================================================
!           If needed, add kappa*Omega to Kau00, which is equivalent to add
!           t_mat to [tau00 - t_mat]
!           ==========================================================
            if (add_SingleSite) then
!              -------------------------------------------------------
               OmegaHat => getOmegaHatMatrix(js2,site=id,atom=ia)
!              -------------------------------------------------------
            endif
            PhiLr_right => getRegSolution(js2,site=id,atom=ia)
            do js1 = 1, n_spin_cant
               PhiLr_left => getRegSolution(js1,site=id,atom=ia)
               js = js + 1
               p_kau00 => kau00(:,:,js)
!              =======================================================
!              gf is the multiple scattering part of the Green function
!              multiplied by r^2:
!                  gf = Z_L*(Tau00-t_matrix)*Z_L^{*}*r^2
!                     = Phi_L*Kau00*Phi_L^{*}*r^2
!              Here implements the 1st method
!
!              If needed, add kappa*Omega to Kau00, which is equivalent to add
!              t_mat to [tau00 - t_mat]
!              =======================================================
               if (add_SingleSite .and. js1 == js2) then
!                 ----------------------------------------------------
                  call zcopy(kmaxk*kmaxk,p_kau00,1,pau00,1)
                  call zaxpy(kmaxk*kmaxk,kappa,OmegaHat,1,pau00,1)
!                 ----------------------------------------------------
                  p_kau00 => pau00
               endif
!
!              *******************************************************
!              =======================================================
!              ppr(ir,klp1,kl) = sum_kl1 PhiLr_left(ir,klp1,kl1) * 
!                                        p_kau00(kl1,kl)
!              -------------------------------------------------------
               call zgemm('n','n',irmax*kmaxp,kmaxk,kmaxk,CONE,PhiLr_left, &
                          irmax*kmaxp,p_kau00,kmaxk,CZERO,ppr,irmax*kmaxp)
!              -------------------------------------------------------
               np = mod(kmaxk,NumPEsInGroup)
               dos_buf = ZERO
               do kl = MyPEinGroup+1, kmaxk-np, NumPEsInGroup
                  m = mofk(kl)
                  klc = kl -2*m
                  cfac = m1m(m)
!                 ====================================================
!                 ppg(ir,klg,klp2;kl) = sum_klp1 (-1)^m * ppr(ir,klp1,klc) *
!                                                         gaunt(klp1,klg,klp2)
!                 ----------------------------------------------------
                  call zgemm('n','n',irmax,kmaxg*kmaxp,kmaxp,cfac,ppr(1,klc), &
                             irmax,gaunt,kmaxp,CZERO,ppg,irmax)
!                 ----------------------------------------------------
!
!                 ====================================================
!                 gf(ir,klg) = sum_{klp2} ppg(ir,klg,klp2;kl) * PhiLr_right(ir,klp2,kl)
!                 ====================================================
                  gf = CZERO
                  do klp2 = 1, kmaxp
                     do klg = 1, kmaxg
                        do ir = 1, mst(id)%iend
                           gf(ir,klg) =  gf(ir,klg) + ppg(ir,klg,klp2)*PhiLr_right(ir,klp2,kl)
                        enddo
                     enddo
                  enddo
!
                  cfac = SQRTm1/PI2
                  do jlg = 1, jmaxg
                     lg = lofj(jlg); mg = mofj(jlg); 
                     klg = (lg+1)*(lg+1)-lg+mg; klgc = (lg+1)*(lg+1)-lg-mg
                     do ir = 1, mst(id)%iend
                        dos_r_jl(ir,jlg) = cfac*(gf(ir,klg) - m1m(mg)*conjg(gf(ir,klgc)))
                     enddo
                  enddo
!
!                 ====================================================
!                 Symmetrizing the DOS, if needed.  Added by Yang on 09-21-2018
!                 ====================================================
                  if (isDosSymmOn) then
                     dos_flags => getSymmetryFlags(id)
                     do jlg = 1, jmaxg
                        if (dos_flags(jlg) == 0) then
                           dos_r_jl(:,jlg) = CZERO
                        endif
                     enddo
                  endif
!                 ====================================================
!
!                 ----------------------------------------------------
                  dos_buf(kl,1) = getVolumeIntegration( id, mst(id)%iend, Grid%r_mesh, &
                                                        jmaxg, 2, dos_r_jl, dos_buf(kl,2) )
!                 ----------------------------------------------------
               enddo ! kl
!
               if (NumPEsInGroup > 1) then
!                 ----------------------------------------------------
                  call GlobalSumInGroup(kGID,dos_buf,kmax_phi_max,2)
!                 ----------------------------------------------------
               endif
!
               do kl = kmaxk-np+1,kmaxk
                  m = mofk(kl)
                  klc = kl -2*m
                  cfac = m1m(m)
!                 ====================================================
!                 ppg(ir,klg,klp2;kl) = sum_klp1 (-1)^m * ppr(ir,klp1,klc) * gaunt(klp1,klg,klp2)
!                 ----------------------------------------------------
                  call zgemm('n','n',irmax,kmaxg*kmaxp,kmaxp,cfac,ppr(1,klc), &
                             irmax,gaunt,kmaxp,CZERO,ppg,irmax)
!                 ----------------------------------------------------
!
!                 ====================================================
!                 gf(ir,klg) = sum_{klp2} ppg(ir,klg,klp2;kl) * PhiLr_right(ir,klp2,kl)
!                 ====================================================
                  gf = CZERO
                  do klp2 = 1, kmaxp
                     do klg = 1, kmaxg
                        do ir = 1, mst(id)%iend
                           gf(ir,klg) =  gf(ir,klg) + ppg(ir,klg,klp2)*PhiLr_right(ir,klp2,kl)
                        enddo
                     enddo
                  enddo
!
                  cfac = SQRTm1/PI2
                  do jlg = 1, jmaxg
                     lg = lofj(jlg); mg = mofj(jlg); 
                     klg = (lg+1)*(lg+1)-lg+mg; klgc = (lg+1)*(lg+1)-lg-mg
                     do ir = 1, mst(id)%iend
                        dos_r_jl(ir,jlg) = cfac*(gf(ir,klg) - m1m(mg)*conjg(gf(ir,klgc)))
                     enddo
                  enddo
!
!                 ====================================================
!                 Symmetrizing the DOS, if needed.  Added by Yang on 09-21-2018
!                 ====================================================
                  if (isDosSymmOn) then
                     dos_flags => getSymmetryFlags(id)
                     do jlg = 1, jmaxg
                        if (dos_flags(jlg) == 0) then
                           dos_r_jl(:,jlg) = CZERO
                        endif
                     enddo
                  endif
!                 ====================================================
!
!                 ----------------------------------------------------
                  dos_buf(kl,1) = getVolumeIntegration( id, mst(id)%iend, Grid%r_mesh, &
                                                        jmaxg, 2, dos_r_jl, dos_buf(kl,2) )
!                 ----------------------------------------------------
               enddo
!
               ns = max(js,is)
               SpaceIntegratedDOS(id)%msdos_cell(ns,ia) = ZERO
               SpaceIntegratedDOS(id)%msdos_mt(ns,ia) = ZERO
               do kl = kmaxk, 1, -1
                  SpaceIntegratedDOS(id)%mspdos_cell(kl,ns,ia) = dos_buf(kl,1)
                  SpaceIntegratedDOS(id)%mspdos_mt(kl,ns,ia) = dos_buf(kl,2)
                  SpaceIntegratedDOS(id)%msdos_cell(ns,ia) =          &
                           SpaceIntegratedDOS(id)%msdos_cell(ns,ia) + &
                           SpaceIntegratedDOS(id)%mspdos_cell(kl,ns,ia)
                  SpaceIntegratedDOS(id)%msdos_mt(ns,ia) =            &
                           SpaceIntegratedDOS(id)%msdos_mt(ns,ia) +   &
                           SpaceIntegratedDOS(id)%mspdos_mt(kl,ns,ia)
               enddo
!              *******************************************************
            enddo ! do js1
         enddo ! do js2
      enddo ! do ia
   enddo ! do id
!
   nullify(p_kau00, pau00, OmegaHat)
!
   end subroutine computeMSPDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMSCellDOS(ks,id,ia) result(dos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ks, id, ia
!
   real (kind=RealKind) :: dos
!
   dos = SpaceIntegratedDOS(id)%msdos_cell(ks,ia)
!
   end function getMSCellDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMSMTSphereDOS(ks,id,ia) result(dos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ks, id, ia
!
   real (kind=RealKind) :: dos
!
   dos = SpaceIntegratedDOS(id)%msdos_mt(ks,ia)
!
   end function getMSMTSphereDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMSCellPDOS(ks,id,ia) result(pdos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ks, id, ia
!
   real (kind=RealKind), pointer :: pdos(:)
!
   pdos => SpaceIntegratedDOS(id)%mspdos_cell(:,ks,ia)
!
   end function getMSCellPDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMSMTSpherePDOS(ks,id,ia) result(pdos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ks, id, ia
!
   real (kind=RealKind), pointer :: pdos(:)
!
   pdos => SpaceIntegratedDOS(id)%mspdos_mt(:,ks,ia)
!
   end function getMSMTSpherePDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genFactors(lmax)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: kmax, jmax, l, m, kl, n, jl
!
   kmax=(lmax+1)*(lmax+1)
   jmax=(lmax+1)*(lmax+2)/2
!
   allocate( lofk(kmax), mofk(kmax), jofk(kmax), m1m(-lmax:lmax) )
   allocate( lofj(jmax), mofj(jmax) )
!
!  ===================================================================
!  calculate the factors: lofk, mofk, jofk, lofj, and mofj............
!  ===================================================================
   kl=0; jl = 0
   do l=0,lmax
      n=(l+1)*(l+2)/2-l
      do m=-l,l
         kl=kl+1
         lofk(kl)=l
         mofk(kl)=m
         jofk(kl)=n+abs(m)
         if (m >= 0) then
            jl = jl + 1
            lofj(jl) = l
            mofj(jl) = m
         endif
      enddo
   enddo
!
!  ===================================================================
!  calculate the factor (-1)**m and store in m1m(-lmax:lmax)..........
!  ===================================================================
   m1m(0)=1
   do m=1,lmax
      m1m(m)=-m1m(m-1)
   enddo
   do m=-1,-lmax,-1
      m1m(m)=-m1m(m+1)
   enddo
!
   end subroutine genFactors
!  ===================================================================
end module MSSolverModule
