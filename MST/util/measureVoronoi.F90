program measureVoronoi
!  ********************************************************************
!  main to measure the voronoi polyhedra distribution
!  ********************************************************************
   use KindParamModule, only : IntKind, RealKind
!
   use GroupCommModule, only : getGroupID
!
   use SystemModule, only : getNumAtoms, getAtomPosition
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use DataServiceCenterModule, only : getDataStorage, RealMark,      &
                                       isDataStorageExisting
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : printPolyhedron
   use PolyhedraModule, only : printPolyhedraTable
   use PolyhedraModule, only : getOutscrSphRadius,                    &
                               getInscrSphRadius, getWignerSeitzRadius
!
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms
!
   implicit none
!
   integer (kind=IntKind) :: i, ig, NumLocalAtoms, NumAtoms
!
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: Bravais(:,:)
   real (kind=RealKind) :: skewness
!
!  -------------------------------------------------------------------
   call StartUp()
!  -------------------------------------------------------------------
!
   NumAtoms = getNumAtoms()
!
   allocate(AtomPosition(1:3,1:NumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
   enddo
!  
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      Bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('testProcMapping','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
   NumLocalAtoms = getLocalNumAtoms()
!
!  -------------------------------------------------------------------
   call initPolyhedra(NumLocalAtoms,Bravais,'main',0)
!  -------------------------------------------------------------------
!
   do i=1,NumLocalAtoms
      ig = getGlobalIndex(i)
!     ----------------------------------------------------------------
      call genPolyhedron(i,ig,NumAtoms,AtomPosition)
!     ----------------------------------------------------------------
      if (NumLocalAtoms < 10000) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        -------------------------------------------------------------
         skewness = (getOutscrSphRadius(i)-getInscrSphRadius(i))/getWignerSeitzRadius(i)
         print *,'skewness = ',skewness
      endif
   enddo
!
!  -------------------------------------------------------------------
   call printPolyhedraTable(getGroupID('Unit Cell'),1,0)
!  -------------------------------------------------------------------
!
   call endPolyhedra()
!
   deallocate(AtomPosition)
!
!  -------------------------------------------------------------------
   call CleanUp()
!  -------------------------------------------------------------------
!
   stop 'OK'
!  
end program measureVoronoi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StartUp()
!  ===================================================================
   use KindParamModule, only : IntKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use DataServiceCenterModule, only : initDataServiceCenter
!
   use InputModule, only : initInput
!
   use SystemModule, only : initSystem, getNumAtoms
!
   use ScfDataModule, only : initScfData
   use ScfDataModule, only : getPotentialTypeParam
!
   use PotentialTypeModule, only : initPotentialType
   use PotentialTypeModule, only : isFullPotential
!
   use MPPModule, only : initMPP, endMPP
!
   use GroupCommModule, only : initGroupComm
!
   use ProcMappingModule, only : initProcMapping, createParallelization
!
   use Atom2ProcModule, only : initAtom2Proc
!
   implicit   none
!
   integer (kind=IntKind) :: def_id, info_id, n
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
   n = getNumAtoms()
!  -------------------------------------------------------------------
   call initProcMapping(isFullPotential(), 'none', 0, n)
!  -------------------------------------------------------------------
   call createParallelization()
!  -------------------------------------------------------------------
   call initAtom2Proc(n, n)
!  -------------------------------------------------------------------
!
   end subroutine StartUp
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine CleanUp()
!  ===================================================================
   use DataServiceCenterModule, only : endDataServiceCenter
!
   use InputModule, only : endInput
!
   use ScfDataModule, only : endScfData
!
   use SystemModule, only : endSystem
!
   use PotentialTypeModule, only : endPotentialType
!
   use MPPModule, only : endMPP
!
   use GroupCommModule, only : endGroupComm
!
   use ProcMappingModule, only : endProcMapping
!
   use Atom2ProcModule, only : endAtom2Proc
!
   implicit none
!
   call endAtom2Proc()
   call endProcMapping()
   call endSystem()
   call endPotentialType()
   call endScfData()
   call endInput()
   call endDataServiceCenter()
   call endGroupComm()
   call endMPP()
!
   end subroutine CleanUp
!  ===================================================================
