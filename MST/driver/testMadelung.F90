program  testMadelung
!  ********************************************************************
!  main to test the Madelung matrix code
!  ********************************************************************
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use TimerModule, only : initTimer, getTime
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       isDataStorageExisting, &
                                       getDataStorage, RealMark
!
   use MadelungModule, only : initMadelung, endMadelung
   use MadelungModule, only : getMadelungMatrix
   use MadelungModule, only : getDLMatrix
   use MadelungModule, only : getDLFactor
   use MadelungModule, only : printMadelungMatrix
!
   use InputModule, only : initInput, endInput, getKeyValue
!
   use SystemModule, only : initSystem, endSystem
   use SystemModule, only : getNumAtoms, getAtomPosition
!
   use ScfDataModule, only : initScfData, endScfData
   use ScfDataModule, only : getPotentialTypeParam
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType, &
                                   isFullPotential
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, CZERO, SQRTm1
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use MPPModule, only : initMPP, endMPP, syncAllPEs, MyPE
!
   use GroupCommModule, only : initGroupComm
!
   use ProcMappingModule, only : initProcMapping, endProcMapping,            &
                                 createParallelization
!
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms
!
   implicit   none
!
   character (len=4) :: istop = 'none'
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: lmax_kkr
   integer (kind=IntKind) :: lmax_phi
   integer (kind=IntKind) :: lmax_pot
   integer (kind=IntKind) :: lmax_rho
   integer (kind=IntKind) :: lmax_mad
   integer (kind=IntKind) :: fstatus
!  
   real (kind=RealKind), pointer :: Bravais(:,:)
!
   integer (kind=IntKind) :: i, kl, n, l, m
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
!
   real (kind=RealKind) :: t0
   real (kind=RealKind) :: a0
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
   real (kind=RealKind), pointer :: madmat(:)
!
   complex (kind=CmplxKind), pointer :: dlm(:,:)
!
   logical :: print_out
!
!  ===================================================================
!  initilize AtomInfo parameters
!  ===================================================================
!
!  -------------------------------------------------------------------
   call initMPP()
   call initGroupComm()
   call initDataServiceCenter()
   call initInput()
   call readInputs(def_id,info_id)
   call initScfData(def_id)
   call initPotentialType(getPotentialTypeParam())
   call initSystem(def_id)
!  -------------------------------------------------------------------
!
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      Bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else 
!     ----------------------------------------------------------------
      call ErrorHandler('testMadelung','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
   NumAtoms = getNumAtoms()
!
!  ===================================================================
!  Initialize the processes mapping module that determines how the
!  parallization will be performed
!  -------------------------------------------------------------------
   call initProcMapping(NumAtoms, 1, 1, isFullPotential(), istop, 0, NumAtoms)
!  -------------------------------------------------------------------
   call createParallelization()
!  -------------------------------------------------------------------
   call initAtom2Proc(NumAtoms, NumAtoms)
!  -------------------------------------------------------------------
!
   LocalNumAtoms=getLocalNumAtoms()
!
!  ===================================================================
!
   allocate(AtomPosition(1:3,1:NumAtoms))
   allocate(GlobalIndex(LocalNumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
   enddo
   do i = 1, LocalNumAtoms
      GlobalIndex(i)=getGlobalIndex(i)
   enddo
!
!  -------------------------------------------------------------------
   fstatus = getKeyValue(info_id,'Default Lmax-T matrix',lmax_kkr)
   fstatus = getKeyValue(info_id,'Default Lmax-Wave Func',lmax_phi)
   fstatus = getKeyValue(info_id,'Default Lmax-Charge Den',lmax_rho)
   fstatus = getKeyValue(info_id,'Default Lmax-Potential',lmax_pot)
!  -------------------------------------------------------------------
   lmax_mad = 2*max(lmax_pot,lmax_rho)
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_mad)
!  -------------------------------------------------------------------
!  
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
   call initGauntFactors(lmax_mad,istop,0)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call initTimer()
   t0 = getTime()
   call initMadelung(LocalNumAtoms,NumAtoms,GlobalIndex,              &
                     lmax_rho,lmax_pot,bravais,AtomPosition,0)
   write(6,'('' time for set up MadelungMatrix ='',1f10.5)')getTime()-t0
!  -------------------------------------------------------------------
   do n=1,LocalNumAtoms
      i = GlobalIndex(n)
      print_out = .false.
      if (MyPE == 0 .and. n < 101) then
         print_out = .true.
      endif
      if (print_out) then
         write(6,'(a,i5,a,3f15.8)')'Atom Index:',n,',  Position:',        &
                                   AtomPosition(1:3,i)
!        -------------------------------------------------------------
         madmat => getMadelungMatrix(n)
!        -------------------------------------------------------------
!        write(6,'(a)')'Madelung Matrix:'
!        write(6,'(i5,2x,d15.8)')(i,madmat(i), i=1,NumAtoms)
         call printMadelungMatrix(1)
         if (lmax_mad > 0 .and. iprint > 0) then
            dlm => getDLMatrix(n,a0)
            write(6,'(a)')'DL Matrix:'
            kl = 0
            do l = 0, lmax_mad
               do m = -l, l
                  kl = kl + 1
                  do i = 1, NumAtoms
                     if (abs(dlm(i,kl)) > 1.0d-6) then
                        write(6,'(3i5,2x,2d15.8)')i,l,m,dlm(i,kl)/a0**l
                     endif
                  enddo
               enddo
            enddo
         endif
      endif
   enddo
!  --------------------------------------------------------------------
   call endMadelung()
   call syncAllPEs()
   call endMPP()
!  --------------------------------------------------------------------
!
end program testMadelung
