module RelScattererModule

!*************************************************************************
! Name: RelScattererModule                                               *
!                                                                        *
! Version 1.0: June 2011 by wkb                                          *
!                                                                        *
! Description:                                                           *
!    Driver/interface for the single scattering routines in              *
!    DiracSolverModule                                                   *
!                                                                        *
! ---------------------------------------------------------------------- *
! subroutine initRelScatterer(num_latoms, l_max, atomic_num, &           *
!                                iprint, stop_routine)                   *
!    Initializes variables and allocates memory used by the routines in  *
!    this module.                                                        *
!  Input: num_latoms - Integer, number of local atoms                    *
!         l_max - 1D integer array of size num_latoms containing the     *
!                 lmax cutoff for each local atom                        *
!         atomic_num - 1D integer array containing the atomic number of  *
!                      each local atom                                   *
!         iprint - integer, output level                                 *
!         stop_routine - character string, name of the routine in which  *
!                        to stop processing                              *
! ---------------------------------------------------------------------- *
! subroutine endRelScatterer()                                           *
!    Cleans the memory used by this module. Must be the last routine     *
!    called in this module unless initRelScatterer is called again.      *
! ---------------------------------------------------------------------- *
! subroutine solveRelSST(e, id)                                          *
!    Run the single-site scattering routines to solve the Dirac equation *
!    and generate t-matrices.                                            *
!  Input: e - complex energy at which to solve the single-site equations *
!         id - (optional) integer, local atom index. If not provided,    *
!              this will solve all atoms on the local processor.         *
! ---------------------------------------------------------------------- *
! function getNumRadialPts(atom) result(nr)                              *
!    Returns the number of points on a given atom's radial mesh.         *
!  Input: atom - integer, local atom index                               *
!  Output: nr - integer, number of radial points on the mesh             *
! ---------------------------------------------------------------------- *
! function getRelTmat(atom) result(tmat)                                 *
!    Returns the single-site scattering matrix for the given atom in     *
!    the global frame of reference.                                      *
!  Input: atom - integer, local atom index                               *
!  Output: tmat - pointer to a 2D complex array containing the           *
!                 single-site t-matrix.                                  *
! ---------------------------------------------------------------------- *
! function getGZ(atom) result(gz)                                        *
!    Returns the large component of the regular solution to the Dirac    *
!    equation on the given atom's radial mesh.                           *
!  Input: atom - integer, local atom index                               *
!  Output: gz - pointer to 3D complex array containing the large         *
!               component of the regular solution                        *
! ---------------------------------------------------------------------- *
! function getFZ(atom) result(fz)                                        *
!    Returns the small component of the regular solution to the Dirac    *
!    equation on the given atom's radial mesh.                           *
!  Input: atom - integer, local atom index                               *
!  Output: fz - pointer to 3D complex array containing the small         *
!               component of the regular solution                        *
! ---------------------------------------------------------------------- *
! function getGJ(atom) result(gj)                                        *
!    Returns the large component of the irregular solution to the Dirac  *
!    equation on the given atom's radial mesh.                           *
!  Input: atom - integer, local atom index                               *
!  Output: gj - pointer to 3D complex array containing the large         *
!               component of the irregular solution                      *
! ---------------------------------------------------------------------- *
! function getFJ(atom) result(fj)                                        *
!    Returns the small component of the irregular solution to the Dirac  *
!    equation on the given atom's radial mesh.                           *
!  Input: atom - integer, local atom index                               *
!  Output: fj - pointer to 3D complex array containing the small         *
!               component of the irregular solution                      *
! ---------------------------------------------------------------------- *
! function get_Nuz(atom) result(nuz)                                     *
!    Returns the number of (kappa',mu') components for (kappa,mu).       *
!  Input: same as above                                                  *
!  Output: nuz - pointer to 1D integer array containing the number of    *
!                (kappa',mu') components for a given atom                *
! ---------------------------------------------------------------------- *
! function get_Indz(atom) result(indz)                                   *
!    Returns the (kappa',mu') indexes for (kappa,mu).                    *
!  Input: same as above                                                  *
!  Output: nuz - pointer to 2D integer array containing the              *
!                (kappa',mu') indexes for a given atom                   *
!*************************************************************************

use KindParamModule, only : IntKind, RealKind, CmplxKind
use MathParamModule, only : CZERO
use ErrorHandlerModule, only : ErrorHandler, WarningHandler

public :: initRelScatterer, &
	  endRelScatterer, &
	  solveRelSST, &
	  getNumRadialPts, &
	  getRelTmat, &
	  getGZ, &
	  getFZ, &
	  getGJ, &
	  getFJ, &
	  get_Nuz, &
	  get_Indz

integer(kind=IntKind), pointer, public :: pI_MemID_R, pQ_MemID_R
integer(kind=IntKind), pointer, public :: pJ_MemID_R, pP_MemID_R

private

type RelScatterCompStruct
	integer (kind=IntKind) :: NumRs
	integer (kind=IntKind), pointer :: nuz(:)
	integer (kind=IntKind), pointer :: indz(:,:)
	complex (kind=CmplxKind), pointer :: gz(:,:,:)
	complex (kind=CmplxKind), pointer :: fz(:,:,:)
	complex (kind=CmplxKind), pointer :: gj(:,:,:)
	complex (kind=CmplxKind), pointer :: fj(:,:,:)
	complex (kind=CmplxKind), pointer :: tmat(:,:)
end type RelScatterCompStruct

type (RelScatterCompStruct), allocatable, target :: RelSolutionComp(:)

logical :: Initialized = .false.
logical :: SSTFlag = .false.

integer (kind=IntKind) :: LocalNumAtoms, GroupID
integer (kind=IntKind), allocatable :: lmax(:), ipvt(:), kmax_kkr(:)
integer (kind=IntKind), parameter :: n_spin_pola = 2 ! Relativistic calculation
integer (kind=IntKind), parameter :: n_spin_cant = 2 ! must be spin canted
integer (kind=IntKind), pointer :: pindex(:,:)
integer (kind=IntKind), allocatable, target :: Q_MemID(:), I_MemID(:)
integer (kind=IntKind), allocatable, target :: P_MemID(:), J_MemID(:)


character (len=50) :: istop

complex (kind=CmplxKind) :: Energy
complex (kind=CmplxKind), pointer :: tmat_g(:,:)
complex (kind=CmplxKind), allocatable, target :: pmat_spinor(:,:)
complex (kind=CmplxKind), allocatable, target :: pmat_s2evec(:,:)

contains

include '../lib/arrayTools.F90'

!=========================================================================
   subroutine initRelScatterer(num_latoms, l_max, getAtomicNumber, &
                               iprint, stop_routine)
!=========================================================================

   use DiracSolverModule, only : initDiracSolver, getSizeOfRMesh
   use GroupCommModule, only : openLocalMemoryInGroup, getGroupID
   use DataServiceCenterModule, only : createDataStorage,        &
                                       setDataStorageLDA,        &
                                       getDataStorage,           &
                                       IntegerType, ComplexType, &
                                       IntegerMark, ComplexMark

   implicit none

   character (len=22), parameter :: sname = 'InitRelScatterer'
   character (len=*), intent(in) :: stop_routine

   integer (kind=IntKind), intent(in) :: iprint(:)
   integer (kind=IntKind), intent(in) :: num_latoms
   integer (kind=IntKind), intent(in) :: l_max(:)
   integer (kind=IntKind), allocatable :: atomic_num(:)
   integer (kind=IntKind) :: i, kmax_kkr_max
   integer (kind=IntKind) :: tsize, DataSize
!
   interface
      function getAtomicNumber(id,ic) result(z)
         use KindParamModule, only : IntKind
         implicit none
         integer (kind=IntKind), intent(in) :: id
         integer (kind=IntKind), intent(in), optional :: ic
         integer (kind=IntKind) :: z
      end function getAtomicNumber
   end interface

   if (Initialized) then
      call WarningHandler('initRelScatterer', 'relScatterer has already been initialized')
   else if (num_latoms < 1) then
      call ErrorHandler('initRelScatterer', 'num_latoms < 1', num_latoms)
   endif

#ifdef DEBUG
   write(6,*) "initRelScatterer: Begin"
   call FlushFile(6)
#endif

   LocalNumAtoms = num_latoms
   SSTFlag = .false.
   istop = stop_routine

   GroupID = getGroupID('Unit Cell')

   allocate( lmax(LocalNumAtoms) )
   allocate( kmax_kkr(LocalNumAtoms) )
   allocate( RelSolutionComp(LocalNumAtoms) )
   allocate( atomic_num(LocalNumAtoms) )

   kmax_kkr_max = 1
   do i=1, LocalNumAtoms
      lmax(i) = l_max(i)
      kmax_kkr(i) = (lmax(i)+1)*(lmax(i)+1)
      kmax_kkr_max = max(kmax_kkr_max, kmax_kkr(i))
      atomic_num(i) = getAtomicNumber(i,1) ! temporary fix!!!
   enddo
   tsize = kmax_kkr_max*kmax_kkr_max

   call initDiracSolver(LocalNumAtoms,lmax,atomic_num,iprint,istop)

   do i=1, LocalNumAtoms
      RelSolutionComp(i)%NumRs = getSizeOfRMesh(i)
   enddo

   allocate( Q_MemID(1), I_MemID(1) )
   pI_MemID_R => I_MemID(1)
   pQ_MemID_R => Q_MemID(1)

   DataSize = tsize*n_spin_cant*n_spin_cant*LocalNumAtoms
   call createDataStorage('Local T-Mat Buffer',DataSize,ComplexType)
   call setDataStorageLDA('Local T-Mat Buffer',tsize*n_spin_cant*n_spin_cant)

   DataSize = 2*LocalNumAtoms
   call createDataStorage('Local Index Buffer',DataSize,IntegerType)
   call setDataStorageLDA('Local Index Buffer',2)

   pindex => getDataStorage('Local Index Buffer',2,LocalNumAtoms,IntegerMark)
   tmat_g => getDataStorage('Local T-Mat Buffer',                        &
                            tsize*n_spin_cant*n_spin_cant,LocalNumAtoms, &
                            ComplexMark)

   tmat_g(:,:) = czero
   pindex(1,1) = lmax(1)
   pindex(2,1) = 1
   do i = 2,LocalNumAtoms
      pindex(1,i) = lmax(i)
      pindex(2,i) = (kmax_kkr(i)*n_spin_cant)**2+pindex(2,i-1)
   enddo

   allocate( pmat_spinor(tsize*n_spin_cant*n_spin_cant,LocalNumAtoms) )
   allocate( ipvt(kmax_kkr_max) )

   call openLocalMemoryInGroup(GroupID,pindex,2,LocalNumAtoms,pI_MemID_R)
   call openLocalMemoryInGroup(GroupID,tmat_g,n_spin_cant*n_spin_cant*tsize, &
                               LocalNumAtoms, pQ_MemID_R)

   Initialized = .true.

#ifdef DEBUG
   write(6,*) "initRelScatterer: End"
#endif
!
   deallocate(atomic_num)

   end subroutine initRelScatterer

!=========================================================================
   subroutine endRelScatterer()
!=========================================================================

   use DiracSolverModule, only : endDiracSolver
   use GroupCommModule, only : closeLocalMemoryInGroup
   use DataServiceCenterModule, only : deleteDataStorage

   implicit none

   integer (kind=IntKind) :: ia

   if (.not.Initialized) then
      call ErrorHandler('endRelScatterer', 'module has not been initialized')
   endif

#ifdef DEBUG
   write(6,*) "endRelScatterer: Begin"
   call FlushFile(6)
#endif

   Initialized = .false.

   do ia=1, LocalNumAtoms
      nullify( RelSolutionComp(ia)%tmat )
      nullify( RelSolutionComp(ia)%gz ) 
      nullify( RelSolutionComp(ia)%fz )
      nullify( RelSolutionComp(ia)%gj )
      nullify( RelSolutionComp(ia)%fj )
   enddo
   deallocate( RelSolutionComp )
   deallocate( lmax, kmax_kkr )

   call closeLocalMemoryInGroup(I_MemID(1))
   call closeLocalMemoryInGroup(Q_MemID(1))
   nullify(tmat_g, pindex)
   call deleteDataStorage('Local Index Buffer')
   call deleteDataStorage('Local T-Mat Buffer')

   nullify( pI_MemID_R, pQ_MemID_R )
   deallocate( Q_MemID, I_MemID )
   deallocate( pmat_spinor, ipvt )

   SSTFlag = .false.

   call endDiracSolver()

#ifdef DEBUG
   write(6,*) "endRelScatterer: End"
   FlushFile(6)
#endif

   end subroutine endRelScatterer

!=========================================================================
   subroutine solveRelSST(e, id)
!=========================================================================

   use DiracSolverModule, only : getRelTMatrix, getGReg, getFReg, &
                                 getGIrr, getFIrr, getnuz, getindz
   use GroupCommModule, only : syncLocalMemoryInGroup
   use RelativityToolsModule, only : repl
   use SpinRotationModule, only : rotateLtoG

   implicit none

   integer (kind=IntKind), intent(in), optional :: id
   integer (kind=IntKind) :: i, i0, i1, t0size, info
   integer (kind=IntKind) :: t0size_ns, kkri_ns

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: gmat(:,:), tmp1(:,:), tmp2(:,:), p1(:)

   if (.not.Initialized) then
      call ErrorHandler('solveRelSST', 'module is not initialized')
   endif

#ifdef DEBUG
   if ( present(id) ) then
      write(6,*) "solveRelSST: Begin Id=",id
   else
      write(6,*) "solveRelSST: Begin"
   endif
   call FlushFile(6)
#endif

   if (present(id)) then
      i0=id
      i1=id
   else
      i0=1
      i1=LocalNumAtoms
   endif

   do i=i0, i1
      RelSolutionComp(i)%tmat => getRelTMatrix(i, e)
      RelSolutionComp(i)%gz => getGReg(i, e)
      RelSolutionComp(i)%fz => getFReg(i, e)
      RelSolutionComp(i)%gj => getGIrr(i, e)
      RelSolutionComp(i)%fj => getFIrr(i, e)
      RelSolutionComp(i)%nuz => getnuz(i, e)
      RelSolutionComp(i)%indz => getindz(i, e)
   enddo

   SSTFlag = .true.

   do i = i0,i1
      t0size = 4*kmax_kkr(i)*kmax_kkr(i)
      kkri_ns =  kmax_kkr(i)*n_spin_cant
      p1 => tmat_g(1:t0size,i)
      gmat => aliasArray2_c(p1,kkri_ns,kkri_ns)

      call repl(gmat, RelSolutionComp(i)%tmat, kkri_ns, kkri_ns)

      nullify(gmat)
   enddo

   if ( .not.present(id) ) then
      call syncLocalMemoryInGroup(pI_MemID_R)
      call syncLocalMemoryInGroup(pQ_MemID_R)
   endif

#ifdef DEBUG
   write(6,*) "solveRelSST: End"
   call FlushFile(6)
#endif

   end subroutine solveRelSST

!===============================================================================
   function getNumRadialPts(atom) result (nr)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom
   integer (kind=IntKind) :: nr

   if (.not.Initialized) then
      call ErrorHandler('getNumRadialPts', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getNumRadialPts', 'invalid atom index', atom)
   endif

   nr = RelSolutionComp(atom)%NumRs

   end function getNumRadialPts

!===============================================================================
   function getRelTmat(atom) result(tmat)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: tmat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelTmat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelTmat', 'invalid number of local atoms', LocalNumAtoms)
   else if (.not.SSTFlag) then
      call ErrorHandler('getRelTmat', 'solveRelSST must be called first')
   endif

   tmat => RelSolutionComp(atom)%tmat

   end function getRelTmat

!===============================================================================
   function getGZ(atom) result(gz)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: gz(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getGZ', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getGZ', 'invalid number of local atoms', LocalNumAtoms)
   else if (.not.SSTFlag) then
      call ErrorHandler('getGZ', 'solveRelSST must be called first')
   endif

   gz => RelSolutionComp(atom)%gz

   end function getGZ

!===============================================================================
   function getFZ(atom) result(fz)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: fz(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getFZ', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getFZ', 'invalid number of local atoms', LocalNumAtoms)
   else if (.not.SSTFlag) then
      call ErrorHandler('getFZ', 'solveRelSST must be called first')
   endif

   fz => RelSolutionComp(atom)%fz

   end function getFZ


!===============================================================================
   function getGJ(atom) result(gj)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: gj(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getGJ', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getGJ', 'invalid number of local atoms', LocalNumAtoms)
   else if (.not.SSTFlag) then
      call ErrorHandler('getGJ', 'solveRelSST must be called first')
   endif

   gj => RelSolutionComp(atom)%gj

   end function getGJ

!===============================================================================
   function getFJ(atom) result(fj)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: fj(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getFJ', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getFJ', 'invalid number of local atoms', LocalNumAtoms)
   else if (.not.SSTFlag) then
      call ErrorHandler('getFJ', 'solveRelSST must be called first')
   endif

   fj => RelSolutionComp(atom)%fj

   end function getFJ

!===============================================================================
   function get_Nuz(atom) result(nuz)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   integer (kind=IntKind), pointer :: nuz(:)

   if (.not.Initialized) then
      call ErrorHandler('get_Nuz', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('get_Nuz', 'invalid number of local atoms', LocalNumAtoms)
   else if (.not.SSTFlag) then
      call ErrorHandler('get_Nuz', 'solveRelSST must be called first')
   endif

   nuz => RelSolutionComp(atom)%nuz

   end function get_Nuz

!===============================================================================
   function get_Indz(atom) result(indz)
!===============================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   integer (kind=IntKind), pointer :: indz(:,:)

   if (.not.Initialized) then
      call ErrorHandler('get_Indz', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('get_Indz', 'invalid number of local atoms', LocalNumAtoms)
   else if (.not.SSTFlag) then
      call ErrorHandler('get_Indz', 'solveRelSST must be called first')
   endif

   indz => RelSolutionComp(atom)%indz

   end function get_Indz

end module RelScattererModule
