program testGroupComm
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE
   use MPPModule, only : MyPE, NumPEs, initMPP, endMPP, syncAllPEs
   use GroupCommModule
!
   use ErrorHandlerModule, only : ErrorHandler
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       isDataStorageExisting, &
                                       getDataStorage, RealMark
   use InputModule, only : initInput, endInput
   use SystemModule, only : initSystem, endSystem
   use ScfDataModule, only : initScfData, endScfData
   use ScfDataModule, only : getPotentialTypeParam
   use PotentialTypeModule, only : initPotentialType, endPotentialType
   use ProcMappingModule
!
   implicit none
!
   character (len=20):: istop
   logical :: ifp
   integer (kind=IntKind) :: iprint
!
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: proc_dim(3), box(3)
   integer (kind=IntKind) :: NumGrids, NumGroups, d, i, j, m, ival(1:3)
!
   real (kind=RealKind), pointer :: Bravais(:,:)
   real (kind=RealKind) :: rval(1:3)
!
   complex (kind=CmplxKind) :: cval(1:3)
!
!  -------------------------------------------------------------------
   call initMPP()
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
      call ErrorHandler('testGroupComm','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
   call initGroupComm()
!
   ifp = .false.
   istop = 'test'
   iprint = 0
!
   call initProcMapping(ifp,istop,iprint)
!
   call createParallelization()
!   call setupProcDim(NumPEs,proc_dim)
!
   box(1:3) = proc_dim(1:3)
!
   call createProcGrid(3,'3-D Proc Grid',proc_dim)
!
   NumGrids = getNumProcGrids()
   print *,'Number of Proc Grids = ',NumGrids
   print *,'Proc Grid ID: ',getProcGridID('3-D Proc Grid')
   do i = 1, NumGrids
      d = getProcGridDimention(i)
      print *,'Grid Dimension: ',d
      print *,'Grid Size: ',getProcGridSize(i)
      do j = 1, d
         print *,'For dim = ',j,', size = ',getProcGridSize(i,j)
         print *,'For dim = ',j,',my coordinate: ',getMyCoordInProcGrid(i,j)
      enddo
      call createBoxGroupFromGrid(i,box,'Unit Cell')
   enddo
!
   print *,'Group ID = ',getGroupID('3-D Proc Grid')
   print *,'Group ID = ',getGroupID('Unit Cell')
   NumGroups = getNumGroups()
   print *,'Number of Groups = ',NumGroups
   do i = 1, NumGroups
      print *,'For group title: ',getGroupLabel(i),', size = ',getNumPEsInGroup(i)
      print *,'MyPE = ',MyPE,', My rank in the group = ',getMyPEinGroup(i)
!
      if (getMyPEinGroup(i) == 0) then
         ival(1:3) = 10
         rval(1:3) = 5.0d0
         cval(1:3) = (1.2d0,2.0d0)
      endif
      call bcastMessageInGroup(i,ival,3,0)
      call bcastMessageInGroup(i,rval,3,0)
      call bcastMessageInGroup(i,cval,3,0)
      print *,'MyPE = ',MyPE,', bcast ival(1:3) = ',ival(1:3)
      print *,'MyPE = ',MyPE,', bcast rval(1:3) = ',rval(1:3)
      print *,'MyPE = ',MyPE,', bcast cval(1:3) = ',cval(1:3)
!
      call GlobalSumInGroup(i,ival,3)
      call GlobalSumInGroup(i,rval,3)
      call GlobalSumInGroup(i,cval,3)
      print *,'MyPE = ',MyPE,', sum ival(1:3) = ',ival(1:3)
      print *,'MyPE = ',MyPE,', sum rval(1:3) = ',rval(1:3)
      print *,'MyPE = ',MyPE,', sum cval(1:3) = ',cval(1:3)
   enddo
!
   call endProcMapping()
   call endGroupComm()
!  ===================================================================
!
!  -------------------------------------------------------------------
   call endSystem()
   call endPotentialType()
   call endScfData()
   call endInput()
   call endDataServiceCenter()
   call endMPP()
!  -------------------------------------------------------------------
   stop
!
end program testGroupComm
!  ===================================================================
