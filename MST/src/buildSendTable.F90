!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine buildSendTable(print_level)
!  ===================================================================
   use KindParamModule, only : IntKind
!
   use ErrorHandlerModule, only : WarningHandler
!
   use GroupCommModule, only : openLocalMemoryInGroup, closeLocalMemoryInGroup
   use GroupCommModule, only : writeRemoteMemoryInGroup, syncLocalMemoryInGroup
   use GroupCommModule, only : getGroupID, getMyPEinGroup
!
   use PublicTypeDefinitionsModule, only : NeighborStruct
!
   use SystemModule, only : getNumAtoms
!
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
!
   use NeighborModule, only : getNeighbor
   use NeighborModule, only : getNumReceives, getMaxReceives
   use NeighborModule, only : setSendTable, printCommunicationTable
!
   implicit none
!
   integer (kind=IntKind), intent(in) ::  print_level(*)
!
   integer (kind=IntKind) :: GroupID
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: MyPEinGroup
   integer (kind=IntKind) :: MemID
   integer (kind=IntKind) :: my_atom, pe, i, j, ig, n, numj, lid(1)
   integer (kind=IntKind) :: tlid(1), tatom, tpe
!
   integer (kind=IntKind), allocatable :: local_buf(:), SendTable(:)
   integer (kind=IntKind), allocatable :: AtomTable(:), LidTable(:)
!
   type (NeighborStruct), pointer :: Neighbor
!
   GroupID = getGroupID('Unit Cell')
   MyPEinGroup = getMyPEinGroup(GroupID)
   GlobalNumAtoms = getNumAtoms()
   LocalNumAtoms = getLocalNumAtoms()
!
   allocate( local_buf(GlobalNumAtoms) )
!
   local_buf(1:GlobalNumAtoms) = 0
!
!  -------------------------------------------------------------------
   call openLocalMemoryInGroup(GroupID,local_buf,GlobalNumAtoms,MemID)
!  -------------------------------------------------------------------
   call syncLocalMemoryInGroup(MemID)
!  -------------------------------------------------------------------
!
   do i = 1, LocalNumAtoms
      my_atom = getGlobalIndex(i)
      Neighbor => getNeighbor(i)
      do j = 1, Neighbor%NumAtoms
         pe = Neighbor%ProcIndex(j)
         lid(1) = Neighbor%LocalIndex(j)
         if (pe /= MyPEinGroup) then
!           -----------------------------------------------------------
            call writeRemoteMemoryInGroup(lid,1,pe,MemID,my_atom)
!           -----------------------------------------------------------
            call syncLocalMemoryInGroup(MemID)
!           -----------------------------------------------------------
            tlid(1) = lid(1); tpe = pe; tatom = my_atom
         endif
      enddo
!
!     ================================================================
!     The following redundant writeRemoteMemory statement is used to ensure
!     the number of syncLocalMemory processes be the same on all processes.
!     ================================================================
      n = getMaxReceives()-getNumReceives(i)
      do j = 1, n
!        -------------------------------------------------------------
         call writeRemoteMemoryInGroup(tlid,1,tpe,MemID,tatom)
         call syncLocalMemoryInGroup(MemID)
!        -------------------------------------------------------------
      enddo
   enddo
!
!  ===================================================================
!  setup the message send table
!  -------------------------------------------------------------------
   n = 0
   do ig = 1, GlobalNumAtoms
      if ( local_buf(ig) > 0 ) then
         n = n + 1
      endif
   enddo
!
   if (n > 0) then
      allocate( AtomTable(n), LidTable(n), SendTable(n) )
      j = 0
      do ig = 1, GlobalNumAtoms
         if ( local_buf(ig) > 0 ) then
            j = j + 1
            AtomTable(j) = ig
            LidTable(j) = local_buf(ig)
         endif
      enddo
      do i = 1, LocalNumAtoms
         numj = 0
         do j = 1, n
            if (LidTable(j) == i) then
               numj = numj + 1
               SendTable(numj) = AtomTable(j)
            endif
         enddo
!        -------------------------------------------------------------
         call setSendTable(i,numj,SendTable)
!        -------------------------------------------------------------
      enddo
      deallocate( AtomTable, LidTable, SendTable )
   else if (maxval(print_level(1:LocalNumAtoms)) >= 0) then
!     ----------------------------------------------------------------
      call WarningHandler('buildSendTable',                           &
           'This processor is not required to send t-matrix',MyPEinGroup)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   call closeLocalMemoryInGroup(MemID)
!  -------------------------------------------------------------------
!  call printCommunicationTable()
!  -------------------------------------------------------------------
!
!  ===================================================================
   deallocate( local_buf )
!
   end subroutine buildSendTable
