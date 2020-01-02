module SendRecvTmatModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct, &
                                           CommTableStruct
!
public :: initSendRecvTmat,        &
          endSendRecvTmat,         &
          getMaxNeighbors,         &
          getNumProcWiseNeighbors, &      ! returns number of unique neighbors of all the local atoms
          getNumNonLocalNeighbors, &      ! returns number of unique neighbors of all the local atoms which are on remote proc
          getProcWiseNeighbors,    &
          getNonLocalNeighbors,    &
          getNumProcWiseSends,     &
          getNumProcWiseReceives,  &
          getProcWiseTargetProcs,  &
          getProcWiseSourceProcs
!
private
   integer (kind=IntKind) :: MyPEinGroup, GroupID
!
   integer (kind=IntKind) :: NumProcWiseNeighbors = 0
   integer (kind=IntKind) :: NumNonLocalNeighbors = 0
   integer (kind=IntKind) :: MaxNumNeighbors = 0
   integer (kind=IntKind) :: NumProcWiseSends = 0
   integer (kind=IntKind) :: NumProcWiseReceives = 0
   integer (kind=IntKind) :: MaxNumTargets = 0
!
   integer (kind=IntKind), allocatable, target :: ProcWiseNeighbors(:)
   integer (kind=IntKind), allocatable, target :: NonLocalNeighbors(:)
   integer (kind=IntKind), allocatable, target :: ProcWiseTargetProcs(:)
   integer (kind=IntKind), allocatable, target :: ProcWiseSourceProcs(:)
!
   integer (kind=IntKind) :: LocalNumAtoms = 0
   integer (kind=IntKind) :: print_level = -1
!
   logical :: Initialized = .false.
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSendRecvTmat(na,iprint)
!  ===================================================================
   use PublicTypeDefinitionsModule, only : NeighborStruct, CommTableStruct
!
   use GroupCommModule, only : getGroupID, getMyPEinGroup
!
   use NeighborModule, only : getNumNeighbors
   use NeighborModule, only : getNeighbor, getSendTable, getRecvTable
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind) :: i, j, n, js, jr
!
   logical :: redundant
!
   type (NeighborStruct), pointer :: Neighbor
   type (CommTableStruct), pointer :: SendTable, RecvTable
!
   LocalNumAtoms = na
   print_level = iprint
!
   GroupID = getGroupID('Unit Cell')
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   MaxNumNeighbors = 0
   do i = 1, LocalNumAtoms
      MaxNumNeighbors = max(MaxNumNeighbors, getNumNeighbors(i))
   enddo
   if (MaxNumNeighbors > 0) then
      allocate( ProcWiseNeighbors(MaxNumNeighbors*LocalNumAtoms) )
      allocate( NonLocalNeighbors(MaxNumNeighbors*LocalNumAtoms) )
      allocate( ProcWiseSourceProcs(MaxNumNeighbors*LocalNumAtoms) )
   endif
!
   MaxNumTargets = 0
   do i = 1, LocalNumAtoms
      SendTable => getSendTable(i)
      MaxNumTargets = max(MaxNumTargets, SendTable%NumComms)
   enddo
   if (MaxNumTargets > 0) then
      allocate( ProcWiseTargetProcs(MaxNumTargets*LocalNumAtoms) )
   endif
!
   NumProcWiseNeighbors = 0
   NumProcWiseReceives = 0
   NumProcWiseSends = 0
   do i = 1, LocalNumAtoms
      Neighbor => getNeighbor(i)
      do n = 1, Neighbor%NumAtoms
         redundant = .false.
         LOOP_j: do j = 1, NumProcWiseNeighbors
            if (Neighbor%GlobalIndex(n) == ProcWiseNeighbors(j)) then
               redundant = .true.
               exit LOOP_j
            endif
         enddo LOOP_j
         if (.not.redundant) then
            NumProcWiseNeighbors = NumProcWiseNeighbors + 1
            ProcWiseNeighbors(NumProcWiseNeighbors) = Neighbor%GlobalIndex(n)
            if (Neighbor%ProcIndex(n) /= MyPEinGroup) then
               NumNonLocalNeighbors = NumNonLocalNeighbors + 1
               NonLocalNeighbors(NumNonLocalNeighbors) = Neighbor%GlobalIndex(n)
            endif
         endif
      enddo
!
      SendTable => getSendTable(i)
      do n = 1, SendTable%NumComms
         redundant = .false.
         LOOP_js: do js = 1, NumProcWiseSends
            if (SendTable%RemoteProcList(n) == ProcWiseTargetProcs(js)) then
               redundant = .true.
               exit LOOP_js
            endif
         enddo LOOP_js
         if (.not.redundant) then
            NumProcWiseSends = NumProcWiseSends + 1
            ProcWiseTargetProcs(NumProcWiseSends) = SendTable%RemoteProcList(n)
         endif
      enddo
!
      RecvTable => getRecvTable(i)
      do n = 1, RecvTable%NumComms
         redundant = .false.
         LOOP_jr: do jr = 1, NumProcWiseReceives
            if (RecvTable%RemoteProcList(n) == ProcWiseSourceProcs(jr)) then
               redundant = .true.
               exit LOOP_jr
            endif
         enddo LOOP_jr
         if (.not.redundant) then
            NumProcWiseReceives = NumProcWiseReceives + 1
            ProcWiseSourceProcs(NumProcWiseReceives) = RecvTable%RemoteProcList(n)
         endif
      enddo
   enddo
!
   Initialized = .true.
!
   end subroutine initSendRecvTmat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSendRecvTmat()
!  ===================================================================
   implicit none
!
   if (MaxNumNeighbors > 0) then
      deallocate( ProcWiseNeighbors )
      deallocate( NonLocalNeighbors )
      deallocate( ProcWiseSourceProcs )
   endif
   if (MaxNumTargets > 0) then
      deallocate( ProcWiseTargetProcs )
   endif
!
   NumProcWiseNeighbors = 0
   NumNonLocalNeighbors = 0
   MaxNumNeighbors = 0
   MaxNumTargets = 0
   NumProcWiseSends = 0
   NumProcWiseReceives = 0
!  
   Initialized = .false.
!
   end subroutine endSendRecvTmat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumProcWiseNeighbors() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumProcWiseNeighbors','Need to initialize the module first')
   endif
!
   n = NumProcWiseNeighbors
!
   end function getNumProcWiseNeighbors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumNonLocalNeighbors() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumNonLocalNeighbors','Need to initialize the module first')
   endif
!
   n = NumNonLocalNeighbors
!
   end function getNumNonLocalNeighbors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMaxNeighbors() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getMaxNeighbors','Need to initialize the module first')
   endif
!
   n = MaxNumNeighbors
!
   end function getMaxNeighbors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcWiseNeighbors() result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: p(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getProcWiseNeighbors','Need to initialize the module first')
   endif
!
   if (NumProcWiseNeighbors > 0) then
      p => ProcWiseNeighbors(1:NumProcWiseNeighbors)
   else
      nullify( p )
   endif
!
   end function getProcWiseNeighbors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNonLocalNeighbors() result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: p(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getNonLocalNeighbors','Need to initialize the module first')
   endif
!
   if (NumNonLocalNeighbors > 0) then
      p => NonLocalNeighbors(1:NumNonLocalNeighbors)
   else
      nullify( p )
   endif
!
   end function getNonLocalNeighbors
!  ===================================================================
!
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumProcWiseSends() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumProcWiseSends','Need to initialize the module first')
   endif
!
   n = NumProcWiseSends
!
   end function getNumProcWiseSends
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumProcWiseReceives() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumProcWiseReceives','Need to initialize the module first')
   endif
!
   n = NumProcWiseReceives
!
   end function getNumProcWiseReceives
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcWiseTargetProcs() result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: p(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getProcWiseTargetProcs','Need to initialize the module first')
   endif
!
   if (NumProcWiseSends > 0) then
      p => ProcWiseTargetProcs(1:NumProcWiseSends)
   else
      nullify( p )
   endif
!
   end function getProcWiseTargetProcs
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcWiseSourceProcs() result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: p(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getProcWiseSourceProcs','Need to initialize the module first')
   endif
!
   if (NumProcWiseReceives > 0) then
      p => ProcWiseSourceProcs(1:NumProcWiseReceives)
   else
      nullify( p )
   endif
!
   end function getProcWiseSourceProcs
!  ===================================================================
end module SendRecvTmatModule
