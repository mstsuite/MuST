module Atom2ProcModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initAtom2Proc,      &
          endAtom2Proc,       &
          getLocalNumAtoms,   &
          getMaxLocalNumAtoms,&
          getGlobalIndex,     &
          getAtom2ProcInGroup,&  ! this will return proc within my "Unit Cell" group
          getLocalIndex,      &
!         getNumAtomGroups,   &      ! not implemented yet
!         getNumAtomsInGroup, &      ! not implemented yet
!         getAtom2Group,      &      ! not implemented yet
!         getNumProcsInGroup, &      ! not implemented yet
!         getMyAtomGroup,     &      ! not implemented yet
!         getAtomGroupID,     &      ! not implemented yet
          printAtom2ProcTable
!
   interface getLocalNumAtoms
      module procedure getLocalNumAtoms0, getLocalNumAtoms1
   end interface ! getLocalNumAtoms
!
   interface getGlobalIndex
      module procedure getGlobalIndex0, getGlobalIndex1
   end interface ! getGlobalIndex
!
private
   integer (kind=IntKind), parameter :: MaxAtomsPerProc = 256
   integer (kind=IntKind) :: GlobalNumAtoms
   integer (kind=IntKind) :: MyProcInGroup  ! >= 0, Group means 'Unit Cell' group
   integer (kind=IntKind) :: NumProcsInGroup
   integer (kind=IntKind) :: MaxLocalNumAtoms
   integer (kind=IntKind), allocatable :: LocalNumAtoms(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:,:)
   integer (kind=IntKind), allocatable :: Atom2ProcInGroup(:)
   integer (kind=IntKind), allocatable :: LocalIndex(:)
!
   logical :: Initialized = .false.
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initAtom2Proc(num_atoms,maxp)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getMyPEinGroup,            &
                               getNumPEsInGroup, GlobalSumInGroup
   use ProcMappingModule, only : isAtomOnMyProc
   implicit none
!
   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in), optional :: maxp
   integer (kind=IntKind) :: i, j, n, ig
!
   if (num_atoms < 1) then
      call ErrorHandler('initAtom2Proc','invalid num_atoms',num_atoms)
   endif
!
   GlobalNumAtoms = num_atoms
   Initialized = .true.
!
   ig = getGroupID('Unit Cell')
   MyProcInGroup = getMyPEinGroup(ig)
   NumProcsInGroup = getNumPEsInGroup(ig)
!
   allocate( LocalNumAtoms(NumProcsInGroup) )
   allocate( Atom2ProcInGroup(GlobalNumAtoms), LocalIndex(GlobalNumAtoms) )
!
   LocalNumAtoms(1:NumProcsInGroup) = 0
   Atom2ProcInGroup(1:GlobalNumAtoms) = 0
   LocalIndex(1:GlobalNumAtoms) = 0
!
   j = 0
   do i = 1, GlobalNumAtoms
      if (isAtomOnMyProc(i)) then
         Atom2ProcInGroup(i) = MyProcInGroup
         j = j + 1
         LocalIndex(i) = j
         LocalNumAtoms(MyProcInGroup+1) = LocalNumAtoms(MyProcInGroup+1) + 1
      endif
   enddo
!
   if ( present(maxp) ) then
      if (maxp < 1) then
         n = MaxAtomsPerProc
      else
         n = maxp
      endif
   else
      n = MaxAtomsPerProc
   endif
!
   if (LocalNumAtoms(MyProcInGroup+1) > n) then
!     ----------------------------------------------------------------
      call ErrorHandler('initAtom2Proc',                              &
                        'No. local atoms exceeds upper limit',        &
                        LocalNumAtoms(MyProcInGroup+1), n)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(ig,Atom2ProcInGroup,GlobalNumAtoms)
   call GlobalSumInGroup(ig,LocalIndex,GlobalNumAtoms)
   call GlobalSumInGroup(ig,LocalNumAtoms,NumProcsInGroup)
!  -------------------------------------------------------------------
!
   MaxLocalNumAtoms = maxval(LocalNumAtoms(1:NumProcsInGroup))
!
   allocate( GlobalIndex(MaxLocalNumAtoms,NumProcsInGroup) )
!
   do i = 1, GlobalNumAtoms
      j = Atom2ProcInGroup(i) + 1
      n = LocalIndex(i)
      GlobalIndex(n,j) = i
   enddo
!
   end subroutine initAtom2Proc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endAtom2Proc()
!  ===================================================================
   implicit none
!
   GlobalNumAtoms = 0
   deallocate( Atom2ProcInGroup, LocalIndex, GlobalIndex, LocalNumAtoms )
!
   Initialized = .false.
!
   end subroutine endAtom2Proc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalNumAtoms0() result (n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getLocalNumAtoms','Need to initialize Atom2ProcModule')
   endif
   n = LocalNumAtoms(MyProcInGroup+1)
!
   end function getLocalNumAtoms0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalNumAtoms1(pe) result (n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: pe
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getLocalNumAtoms','Need to initialize Atom2ProcModule')
   else if (pe < 0 .or. pe > NumProcsInGroup-1) then
      call ErrorHandler('getLocalNumAtoms','invalid processor index',pe)
   endif
   n = LocalNumAtoms(pe+1)
!
   end function getLocalNumAtoms1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMaxLocalNumAtoms() result (n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getMaxLocalNumAtoms',                        &
                        'Need to initialize Atom2ProcModule')
   endif
   n = MaxLocalNumAtoms
!
   end function getMaxLocalNumAtoms
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtom2ProcInGroup(i) result (pe)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: pe
!
   if (.not.Initialized) then
      call ErrorHandler('getAtom2ProcInGroup','Need to initialize Atom2ProcModule')
   else if (i < 1 .or. i > GlobalNumAtoms) then
      call ErrorHandler('getAtom2ProcInGroup','invalid atom index',i)
   endif
   pe = Atom2ProcInGroup(i)
!
   end function getAtom2ProcInGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalIndex(ig) result (il)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ig
   integer (kind=IntKind) :: il
!
   if (.not.Initialized) then
      call ErrorHandler('getLocalIndex','Need to initialize Atom2ProcModule')
   else if (ig < 0 .or. ig > GlobalNumAtoms) then
      call ErrorHandler('getLocalIndex','invalid global atom index',ig)
   endif
   il = LocalIndex(ig)
!
   end function getLocalIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalIndex0(local_index) result(id)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: local_index
   integer (kind=IntKind) :: id
!
   if (.not.Initialized) then
      call ErrorHandler('getGlobalIndex','Need to initialize Atom2ProcModule')
   else if (local_index < 1 .or. local_index > LocalNumAtoms(MyProcInGroup+1)) then
      call ErrorHandler('getGlobalIndex','Invalid local index',local_index)
   endif
   id = GlobalIndex(local_index,MyProcInGroup+1)
!
   end function getGlobalIndex0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalIndex1(local_index,pe) result(id)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: local_index
   integer (kind=IntKind), intent(in) :: pe
   integer (kind=IntKind) :: id
!
   if (.not.Initialized) then
      call ErrorHandler('getGlobalIndex','Need to initialize Atom2ProcModule')
   else if (pe < 0 .or. pe > NumProcsInGroup-1) then
      call ErrorHandler('getGlobalIndex','Invalid processor index',pe)
   else if (local_index < 1 .or. local_index > LocalNumAtoms(pe+1)) then
      call ErrorHandler('getGlobalIndex','Invalid local index',local_index)
   endif
   id = GlobalIndex(local_index,pe+1)
!
   end function getGlobalIndex1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printAtom2ProcTable()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   if (.not.Initialized) then
      call WarningHandler('printAtom2ProcTable',                      &
                          'Empty table. Need to initialize Atom2ProcModule')
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,20x,a)')'***********************************'
   write(6,'( 20x,a )')'* Output from printAtom2ProcTable *'
   write(6,'(20x,a,/)')'***********************************'
   write(6,'(80(''=''))')
   write(6,'(a)')' Processor    Global Index'
   write(6,'(80(''-''))')
   do i=1,NumProcsInGroup
      write(6,'(2x,i5,5x,13i5)')i-1,(GlobalIndex(j,i),j=1,LocalNumAtoms(i))
   enddo
   write(6,'(80(''=''))')

   end subroutine printAtom2ProcTable
!  ===================================================================
end module Atom2ProcModule
