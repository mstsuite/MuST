program testCoreStates
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use DataServiceCenterModule, only : isDataStorageExisting, &
                                       getDataStorage, RealMark
!
   use SystemModule, only : getNumAtoms, getAtomPosition
!
   use ScfDataModule, only : n_spin_pola, n_spin_cant
   use ScfDataModule, only : istop, EvBottom, isNonRelativisticCore
   use ScfDataModule, only : ngaussr, ngaussq
   use ScfDataModule, only : isFrozenCore
!
   use AtomModule, only : getPhiLmax, getStepFuncLmax, getPotLmax
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, SQRTm1, TEN2m8, PI4, PI
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use MadelungModule, only : initMadelung, endMadelung
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags
!
   use CoreStatesModule, only : initCoreStates, calCoreStates, endCoreStates
   use CoreStatesModule, only : readCoreStates, printCoreStates, printCoreDensity
   use CoreStatesModule, only : readCoreDensity, writeCoreDensity
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : getVolumeIntegration
!
   use OutputModule, only : getStandardOutputLevel
!
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential
!
   implicit   none
!
   logical :: FrozenCoreFileExist = .false.
!
   character (len=50) :: FrozenCoreFileName = 'FrozenCoreDensity.dat'
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: NumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: lmax_max, lmax_pot_max
   integer (kind=IntKind) :: id, ig
!
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:), lmax_step(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
!
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
!  -------------------------------------------------------------------
   call startProcess()
   NumAtoms = getNumAtoms()
   LocalNumAtoms=getLocalNumAtoms()
!  -------------------------------------------------------------------
!
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('testSineMatrixPole','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
   allocate(atom_print_level(1:LocalNumAtoms))
   allocate(lmax_step(LocalNumAtoms),lmax_pot(LocalNumAtoms))
   allocate(GlobalIndex(LocalNumAtoms), AtomPosition(3,NumAtoms))
!
   do ig = 1, NumAtoms
      AtomPosition(1:3,ig)=getAtomPosition(ig)
   enddo
!
   lmax_max = 0
   lmax_pot_max = 0
   do id = 1, LocalNumAtoms
      lmax_pot(id) = getPotLmax(id)
      lmax_step(id)  = getStepFuncLmax(id)
      lmax_max = max(lmax_max,2*getPhiLmax(id),lmax_step(id))
      lmax_pot_max = max(lmax_pot_max,lmax_pot(id))
      atom_print_level(id) = getStandardOutputLevel(id)
      GlobalIndex(id)=getGlobalIndex(id)
   enddo
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_max,istop,iprint)
!  -------------------------------------------------------------------
   call setupRadGridAndCell(LocalNumAtoms,lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize potential module
!  -------------------------------------------------------------------
   call initMadelung(LocalNumAtoms,NumAtoms,GlobalIndex,              &
                     lmax_pot_max,lmax_pot_max,bravais,AtomPosition,0)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize SystemSymmetry, Potential, and CoreStates
!  ===================================================================
   call initSystemSymmetry( NumAtoms, LocalNumAtoms, lmax_pot, lmax_step, atom_print_level )
!  -------------------------------------------------------------------
   call calSymmetryFlags()
!  -------------------------------------------------------------------
   call initPotential(LocalNumAtoms,lmax_pot,lmax_step,               &
                      n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
   call initCoreStates(LocalNumAtoms,EvBottom,n_spin_pola,         &
                       isNonRelativisticCore(),istop,1)
!  ===================================================================
!  read potential data
!  -------------------------------------------------------------------
   call readPotential()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  read core states and density data of local atoms from file
!  ===================================================================
   call readCoreStates()
!  -------------------------------------------------------------------
   if (isFrozenCore(fcf_name=FrozenCoreFileName,fcf_exist=FrozenCoreFileExist)) then
      if (FrozenCoreFileExist) then
         write(6,'(/,a)')'Frozen core density is read from a file'
!        -------------------------------------------------------------
         call readCoreDensity(FrozenCoreFileName)
!        -------------------------------------------------------------
      else
         write(6,'(/,a)')'Frozen core calculation'
!        -------------------------------------------------------------
         call calCoreStates()
!        -------------------------------------------------------------
         write(6,'(/,a)')'Frozen core density is written to a file'
!        -------------------------------------------------------------
         call writeCoreDensity(FrozenCoreFileName)
!        -------------------------------------------------------------
      endif
   else
      write(6,'(/,a)')'Calculate the core states'
!     ----------------------------------------------------------------
      call calCoreStates()
!     ----------------------------------------------------------------
   endif
!
   do id = 1,LocalNumAtoms
      if (atom_print_level(id) >= 0) then
!        -------------------------------------------------------------
         call printCoreStates(id)
         call printCoreDensity(id,derivative=.true.)
!        -------------------------------------------------------------
      endif
   enddo
!
   deallocate(GlobalIndex,AtomPosition)
   deallocate(atom_print_level,lmax_step,lmax_pot)
!
!  -------------------------------------------------------------------
   call endCoreStates()
   call endPotential()
   call endSystemSymmetry()
   call endMadelung()
   call endGauntFactors()
   call endSphericalHarmonics()
   call finishProcess()
!  -------------------------------------------------------------------
!
end program testCoreStates
