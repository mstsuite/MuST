module NeighborModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct, &
                                           CommTableStruct
!
public :: initNeighbor,            &
          endNeighbor,             &
          getNeighbor,             &
          setNeighbor,             &
          getNumNeighbors,         &
          sortNeighbors,           &
          getNeighborEIZ,          &
          setNeighborEIZ,          &
          getNumNeighborsEIZ,      &
          sortNeighborsEIZ,        &
          getNumReceives,          &
          getMaxReceives,          &
          setMaxReceives,          &
          setSendTable,            &
          getSendTable,            &
          getRecvTable,            &
          printNeighbor,           &
          printCommunicationTable
!
!  ===================================================================
!  Define NeighborStruct data type. It contains cluster/LIZ information
!  This is made public so that the neighbor data is constructed for each
!  local atom before initMST is called.
!  Moved to PublicTypeDefinitionsModule.F90
!  ===================================================================
!
private
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: MaxReceives
   integer (kind=IntKind), allocatable :: print_level(:)
!
   integer (kind=IntKind) :: MyPEinGroup, GroupID, GroupSize
!
   type (NeighborStruct), allocatable, target :: Neighbor(:)
   type (NeighborStruct), allocatable, target :: NeighborEIZ(:)
!
   type (CommTableStruct), allocatable, target :: SendTable(:)
   type (CommTableStruct), allocatable, target :: RecvTable(:)
!
   logical :: Initialized = .false.
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initNeighbor(na,iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getMyPEinGroup, getNumPEsInGroup
   implicit none
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: iprint(na)
   integer (kind=IntKind) :: id
!
   if (Initialized) then
      call ErrorHandler('initNeighbor','Module is already initialized')
   else if (na < 1) then
      call ErrorHandler('initNeighbor','invalid number of local atoms',na)
   endif
!
   GroupID = getGroupID('Unit Cell')
   MyPEinGroup = getMyPEinGroup(GroupID)
   GroupSize = getNumPEsInGroup(GroupID)
!
   allocate( Neighbor(na), NeighborEIZ(na), SendTable(na), RecvTable(na), &
             print_level(na) )
   do id=1,na
      Neighbor(id)%NumAtoms = 0
      Neighbor(id)%NumReceives = 0
      Neighbor(id)%NumShells = 0
      NeighborEIZ(id)%NumAtoms = 0
      NeighborEIZ(id)%NumReceives = 0
      NeighborEIZ(id)%NumShells = 0
      SendTable(id)%NumRequests = 0
      RecvTable(id)%NumRequests = 0
      SendTable(id)%NumComms = 0
      RecvTable(id)%NumComms = 0
      print_level(id) = iprint(id)
   enddo
   LocalNumAtoms = na
   MaxReceives = 0
!
   Initialized = .true.
!
   end subroutine initNeighbor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endNeighbor()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id
!
   if (.not.Initialized) then
      call ErrorHandler('endNeighbor','The module is not initialized')
   endif
!
   do id =1,LocalNumAtoms
      if (Neighbor(id)%NumAtoms > 0) then
         deallocate( Neighbor(id)%Z, Neighbor(id)%Lmax )
         deallocate( Neighbor(id)%ProcIndex, Neighbor(id)%LocalIndex )
         deallocate( Neighbor(id)%GlobalIndex, Neighbor(id)%Position )
         deallocate( Neighbor(id)%IndexMN, Neighbor(id)%Rmunu )
         deallocate( Neighbor(id)%ShellIndex, Neighbor(id)%ShellRad )
         Neighbor(id)%NumAtoms = 0
         Neighbor(id)%NumReceives = 0
         Neighbor(id)%NumShells = 0
      endif
   enddo
   do id =1,LocalNumAtoms
      if (NeighborEIZ(id)%NumAtoms > 0) then
         deallocate( NeighborEIZ(id)%Z, NeighborEIZ(id)%Lmax )
         deallocate( NeighborEIZ(id)%ProcIndex)
         deallocate( NeighborEIZ(id)%LocalIndex )
         deallocate( NeighborEIZ(id)%GlobalIndex )
         deallocate( NeighborEIZ(id)%Position )
         deallocate( NeighborEIZ(id)%ShellIndex,NeighborEIZ(id)%ShellRad )
         NeighborEIZ(id)%NumAtoms = 0
         NeighborEIZ(id)%NumReceives = 0
         NeighborEIZ(id)%NumShells = 0
      endif
   enddo
   deallocate(Neighbor, NeighborEIZ, print_level)
!
   do id =1,LocalNumAtoms
      if (SendTable(id)%NumRequests > 0) then
         deallocate( SendTable(id)%RemoteAtomList )
         deallocate( SendTable(id)%RemoteProcList )
      endif
      if (RecvTable(id)%NumRequests > 0) then
         deallocate( RecvTable(id)%RemoteAtomList )
         deallocate( RecvTable(id)%RemoteProcList )
      endif
   enddo
   deallocate( SendTable, RecvTable )
!
   MaxReceives = 0
!
   Initialized = .false.
!
   end subroutine endNeighbor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNeighbor(id) result(nb)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   type (NeighborStruct), pointer :: nb
!
   if (.not.Initialized) then
      call ErrorHandler('getNeighbor','Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getNeighbor','invalid local atom index',id)
   endif
!
   nb => Neighbor(id)
!
   end function getNeighbor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNeighborEIZ(id) result(nb)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   type (NeighborStruct), pointer :: nb
!
   if (.not.Initialized) then
      call ErrorHandler('getNeighbor','Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getNeighbor','invalid local atom index',id)
   endif
!
   nb => NeighborEIZ(id)
!
   end function getNeighborEIZ
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setNeighbor(id,nb)
!  ===================================================================
   use SystemModule, only : getNumAtoms
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i, n, nsa, ind_max, nshells
!
   type (NeighborStruct), intent(in) :: nb
!
   n = nb%NumAtoms
   nshells = nb%NumShells
!
   if (.not.Initialized) then
      call ErrorHandler('setNeighbor','Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setNeighbor','invalid local atom index',id)
   else if (n < 1) then
      call ErrorHandler('setNeighbor','invalid number of neighbor atoms',n)
   endif
!
   if (Neighbor(id)%NumAtoms > 0) then
      deallocate( Neighbor(id)%Z, Neighbor(id)%Lmax )
      deallocate( Neighbor(id)%ProcIndex, Neighbor(id)%LocalIndex )
      deallocate( Neighbor(id)%GlobalIndex, Neighbor(id)%Position )
      deallocate( Neighbor(id)%ShellIndex, Neighbor(id)%ShellRad )
      deallocate( Neighbor(id)%IndexMN )
      deallocate( Neighbor(id)%Rmunu )
   endif
   Neighbor(id)%NumAtoms = n
   allocate( Neighbor(id)%Z(n) )
   allocate( Neighbor(id)%Lmax(n) )
   allocate( Neighbor(id)%ProcIndex(n) )
   allocate( Neighbor(id)%LocalIndex(n) )
   allocate( Neighbor(id)%GlobalIndex(n) )
   allocate( Neighbor(id)%Position(3,n) )
   allocate( Neighbor(id)%ShellIndex(n))
   Neighbor(id)%NumShells = nshells
   allocate( Neighbor(id)%ShellRad(nshells))
!
   nsa = getNumAtoms()
!   nsa = size( nb%IndexMN, dim=1 )
   allocate( Neighbor(id)%IndexMN(nsa) )
   Neighbor(id)%IndexMN = nb%IndexMN
   ind_max = 0
   do i = 1,nsa
      ind_max = max(ind_max,nb%IndexMN(i))
   enddo
!   ind_max = size( nb%Rmunu, dim=2 )
   allocate( Neighbor(id)%Rmunu(3,ind_max ,nsa) )
   Neighbor(id)%Rmunu = nb%Rmunu
!
   do i=1,n
      Neighbor(id)%Z(i) = nb%Z(i)
      Neighbor(id)%Lmax(i) = nb%Lmax(i)
      Neighbor(id)%ProcIndex(i) = nb%ProcIndex(i)
      Neighbor(id)%LocalIndex(i) = nb%LocalIndex(i)
      Neighbor(id)%GlobalIndex(i) = nb%GlobalIndex(i)
      Neighbor(id)%Position(1:3,i) = nb%Position(1:3,i)
      Neighbor(id)%ShellIndex(i) = nb%ShellIndex(i)
   enddo
   Neighbor(id)%NumReceives = nb%NumReceives
   Neighbor(id)%NumShells = nb%NumShells
   do i=1,nshells
      Neighbor(id)%ShellRad(i)=nb%ShellRad(i)
   enddo
!
!  -------------------------------------------------------------------
   call setupRecvTable(id)
!  -------------------------------------------------------------------
!
   end subroutine setNeighbor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setNeighborEIZ(id,nb)
!  ===================================================================
   use SystemModule, only : getNumAtoms
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i, j, n, nshells
!
   type (NeighborStruct), intent(in) :: nb
!
   n = nb%NumAtoms
   nshells = nb%NumShells
!
   if (.not.Initialized) then
      call ErrorHandler( 'setNeighborEIZ', & 
                         'Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler( 'setNeighborEIZ', & 
                         'invalid local atom index',id)
   else if (n < 1) then
      call ErrorHandler( 'setNeighborEIZ', & 
                         'invalid number of neighbor atoms',n)
   endif
!
   if (NeighborEIZ(id)%NumAtoms > 0) then
      deallocate( NeighborEIZ(id)%Z, NeighborEIZ(id)%Lmax )
      deallocate( NeighborEIZ(id)%ProcIndex, NeighborEIZ(id)%LocalIndex )
      deallocate( NeighborEIZ(id)%GlobalIndex, NeighborEIZ(id)%Position )
      deallocate( NeighborEIZ(id)%ShellIndex, NeighborEIZ(id)%ShellRad )
   endif
   NeighborEIZ(id)%NumAtoms = n
   allocate( NeighborEIZ(id)%Z(n) )
   allocate( NeighborEIZ(id)%Lmax(n) )
   allocate( NeighborEIZ(id)%ProcIndex(n) )
   allocate( NeighborEIZ(id)%LocalIndex(n) )
   allocate( NeighborEIZ(id)%GlobalIndex(n) )
   allocate( NeighborEIZ(id)%Position(3,n) )
   allocate( NeighborEIZ(id)%ShellIndex(n) )
   NeighborEIZ(id)%NumShells = nshells
!
   do i=1,n
      NeighborEIZ(id)%Z(i) = nb%Z(i)
      NeighborEIZ(id)%Lmax(i) = nb%Lmax(i)
      NeighborEIZ(id)%ProcIndex(i) = nb%ProcIndex(i)
      NeighborEIZ(id)%LocalIndex(i) = nb%LocalIndex(i)
      NeighborEIZ(id)%GlobalIndex(i) = nb%GlobalIndex(i)
      do j = 1,3
         NeighborEIZ(id)%Position(j,i) = nb%Position(j,i)
      enddo
      NeighborEIZ(id)%ShellIndex(i) = nb%ShellIndex(i)
   enddo
   NeighborEIZ(id)%NumShells = nb%NumShells
   do i=1,nshells
      NeighborEIZ(id)%ShellRad(i)=nb%ShellRad(i)
   enddo
!
   end subroutine setNeighborEIZ
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumNeighbors(id) result(nn)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: nn
!
   if (.not.Initialized) then
      call ErrorHandler('getNumNeighbor','Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getNumNeighbor','invalid local atom index',id)
   endif
!
   nn = Neighbor(id)%NumAtoms
!
   end function getNumNeighbors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumNeighborsEIZ(id) result(nn)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: nn
!
   if (.not.Initialized) then
      call ErrorHandler( 'getNeighborEIZ', & 
                         'Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getNeighborEIZ','invalid local atom index',id)
   endif
!
   nn = NeighborEIZ(id)%NumAtoms
!
   end function getNumNeighborsEIZ
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printNeighbor(id)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i, j, n, nshells
   integer (kind=IntKind),allocatable:: nat(:)
!
   if (.not.Initialized) then
      call ErrorHandler('printNeighbor','Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('printNeighbor','invalid local atom index',id)
   endif
!
   n = Neighbor(id)%NumAtoms
   if (n < 1) then
      call WarningHandler('printNeighbor','Neighbor is empty')
      return
   endif
   nshells = Neighbor(id)%NumShells
!
   write(6,'(/,80(''-''))')
   write(6,'(/,25x,a)') ' *****************************'
   write(6,'(25x,a)')   ' * Output from printNeighbor *'
   write(6,'(25x,a,/)') ' *****************************'
   write(6,'(/,80(''=''))')
   write(6,'(a,i5)')' Local Atom Index = ', id
!  write(6,*) "IndexMN :: ", Neighbor(id)%IndexMN
   write(6,'(a,i5)')' Number of neighboring atoms in LIZ = ',n
   write(6,'(a)')   ' Neighboring Atom List:'
   write(6,'(a)')   &
' ============================================================================'
   write(6,'(a)')   &
      '   Z    Lmax    Proc   L-Index     G-Index shell            Position'
   write(6,'(a)')   &
' ----------------------------------------------------------------------------'
   allocate(nat(nshells))
   nat=0
   do j=1,nshells
   do i=1, n
    if(Neighbor(id)%ShellIndex(i) == j) then
      write(6,'(x,i3,5x,i2,4x,i5,5x,i5,4x,i5,3x,i4,x,3f10.5)')     &
            Neighbor(id)%Z(i), Neighbor(id)%Lmax(i),               &
            Neighbor(id)%ProcIndex(i), Neighbor(id)%LocalIndex(i), &
            Neighbor(id)%GlobalIndex(i),Neighbor(id)%ShellIndex(i),&
            Neighbor(id)%Position(1:3,i)
      nat(j)=nat(j)+1
    endif
   enddo
   enddo
   write(6,'(a)')   &
' ============================================================================'
   write(6,'(a,i5)')' Number of shells: ', nshells
   do i=1,nshells
    write(6,'(a,i3,a,3x,a,f12.5,5x,a,i5)')' shell ',i,':','radius:', &
      Neighbor(id)%ShellRad(i),'atoms: ',nat(i)
   enddo
   deallocate(nat)
!
   write(6,'(a)')   &
' ============================================================================'
   end subroutine printNeighbor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function sortNeighbors( id)                     result( pnb )
!  ===================================================================
   use MathParamModule, only : zero, ten2m8
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   type (NeighborStruct), pointer :: pnb
!
   integer (kind=IntKind) :: i, j, ztmp, ltmp, petmp, lindtmp, gindtmp
   integer (kind=IntKind) :: shelltmp
!
   real (kind=RealKind) :: pri(3), prj(3), postmp(3)
!
   if (.not.Initialized) then
      call ErrorHandler('setNeighbor','Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setNeighbor','invalid local atom index',id)
   endif
!
!  ===================================================================
!  sort neighbors such that for a specific radia the order will be
!  descending on x, y, z components of the positions
!  ===================================================================
   if (print_level(id) > 0) then
      write(6,'(a)')   &
      '   Z    Lmax    Proc    L-Index     G-Index                 Position'
      write(6,'(a)')   &
' ----------------------------------------------------------------------------'
      do i=1, Neighbor(id)%NumAtoms
         write(6,'(x,i3,5x,i2,4x,i5,5x,i5,7x,i5,5x,3f10.5)')          &
               Neighbor(id)%Z(i), Neighbor(id)%Lmax(i),               &
               Neighbor(id)%ProcIndex(i), Neighbor(id)%LocalIndex(i), &
               Neighbor(id)%GlobalIndex(i), Neighbor(id)%Position(1:3,i)
      enddo
      write(6,'(a)')   &
' ============================================================================'
      write(6,'(/,a)') &
            "the neighboring sites of Nu neighbor of my mu neighbor :: "
      write(6,'(a)')   &
' ----------------------------------------------------------------------------'
      do i=1, size( Neighbor(id)%IndexMN )
         do j=1, Neighbor(id)%IndexMN(i)
            write(6,*) " Global Atom:",i,", Contributing sites:", & 
                       Neighbor(id)%IndexMN(i)
            write(6,'(a,3f10.5)') " Rmunu:", Neighbor(id)%Rmunu(1:3,j,i)
         enddo
      enddo
      write(6,'(a,/)')   &
' ============================================================================'
   endif
!
!
   do i = 1,Neighbor(id)%NumAtoms
      pri = Neighbor(id)%Position(1:3,i)
      do j = 1,Neighbor(id)%NumAtoms
         prj = Neighbor(id)%Position(1:3,j)
         if ( prj(1) < pri(1) ) then
            postmp = Neighbor(id)%Position(1:3,i)
            ztmp = Neighbor(id)%Z(i)
            ltmp = Neighbor(id)%Lmax(i)
            petmp = Neighbor(id)%ProcIndex(i)
            lindtmp = Neighbor(id)%LocalIndex(i)
            gindtmp = Neighbor(id)%GlobalIndex(i)
            shelltmp = Neighbor(id)%ShellIndex(i)
!
            Neighbor(id)%Position(1:3,i) = Neighbor(id)%Position(1:3,j)
            Neighbor(id)%Z(i) = Neighbor(id)%Z(j)
            Neighbor(id)%Lmax(i) = Neighbor(id)%Lmax(j)
            Neighbor(id)%ProcIndex(i) = Neighbor(id)%ProcIndex(j)
            Neighbor(id)%LocalIndex(i) = Neighbor(id)%LocalIndex(j)
            Neighbor(id)%GlobalIndex(i) = Neighbor(id)%GlobalIndex(j)
            Neighbor(id)%ShellIndex(i) = Neighbor(id)%ShellIndex(j)
!
            Neighbor(id)%Position(1:3,j) = postmp
            Neighbor(id)%Z(j) = ztmp
            Neighbor(id)%Lmax(j) = ltmp
            Neighbor(id)%ProcIndex(j) = petmp
            Neighbor(id)%LocalIndex(j) = lindtmp
            Neighbor(id)%GlobalIndex(j) = gindtmp
            Neighbor(id)%ShellIndex(j) = shelltmp
         endif
         if ( prj(2) < pri(2) .and. pri(1)==prj(1) ) then
            postmp = Neighbor(id)%Position(1:3,i)
            ztmp = Neighbor(id)%Z(i)
            ltmp = Neighbor(id)%Lmax(i)
            petmp = Neighbor(id)%ProcIndex(i)
            lindtmp = Neighbor(id)%LocalIndex(i)
            gindtmp = Neighbor(id)%GlobalIndex(i)
            shelltmp = Neighbor(id)%ShellIndex(i)
!
            Neighbor(id)%Position(1:3,i) = Neighbor(id)%Position(1:3,j)
            Neighbor(id)%Z(i) = Neighbor(id)%Z(j)
            Neighbor(id)%Lmax(i) = Neighbor(id)%Lmax(j)
            Neighbor(id)%ProcIndex(i) = Neighbor(id)%ProcIndex(j)
            Neighbor(id)%LocalIndex(i) = Neighbor(id)%LocalIndex(j)
            Neighbor(id)%GlobalIndex(i) = Neighbor(id)%GlobalIndex(j)
            Neighbor(id)%ShellIndex(i) = Neighbor(id)%ShellIndex(j)
!
            Neighbor(id)%Position(1:3,j) = postmp
            Neighbor(id)%Z(j) = ztmp
            Neighbor(id)%Lmax(j) = ltmp
            Neighbor(id)%ProcIndex(j) = petmp
            Neighbor(id)%LocalIndex(j) = lindtmp
            Neighbor(id)%GlobalIndex(j) = gindtmp
            Neighbor(id)%ShellIndex(j) = shelltmp
         endif
         if ( prj(3) < pri(3) .and. pri(2) == prj(2) .and.             &
              pri(1)==prj(1)  ) then
            postmp = Neighbor(id)%Position(1:3,i)
            ztmp = Neighbor(id)%Z(i)
            ltmp = Neighbor(id)%Lmax(i)
            petmp = Neighbor(id)%ProcIndex(i)
            lindtmp = Neighbor(id)%LocalIndex(i)
            gindtmp = Neighbor(id)%GlobalIndex(i)
            shelltmp = Neighbor(id)%ShellIndex(i)
!
            Neighbor(id)%Position(1:3,i) = Neighbor(id)%Position(1:3,j)
            Neighbor(id)%Z(i) = Neighbor(id)%Z(j)
            Neighbor(id)%Lmax(i) = Neighbor(id)%Lmax(j)
            Neighbor(id)%ProcIndex(i) = Neighbor(id)%ProcIndex(j)
            Neighbor(id)%LocalIndex(i) = Neighbor(id)%LocalIndex(j)
            Neighbor(id)%GlobalIndex(i) = Neighbor(id)%GlobalIndex(j)
            Neighbor(id)%ShellIndex(i) = Neighbor(id)%ShellIndex(j)
!
            Neighbor(id)%Position(1:3,j) = postmp
            Neighbor(id)%Z(j) = ztmp
            Neighbor(id)%Lmax(j) = ltmp
            Neighbor(id)%ProcIndex(j) = petmp
            Neighbor(id)%LocalIndex(j) = lindtmp
            Neighbor(id)%GlobalIndex(j) = gindtmp
            Neighbor(id)%ShellIndex(j) = shelltmp
         endif
      enddo
   enddo
!
!   call printNeighbor(id)
!
   pnb => Neighbor(id)
!
   end function sortNeighbors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function sortNeighborsEIZ( id)                     result( pnb )
!  ===================================================================
   use MathParamModule, only : zero, ten2m8
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   type (NeighborStruct), pointer :: pnb
!
   integer (kind=IntKind) :: i, j, ztmp, ltmp, petmp, lindtmp, gindtmp
   integer (kind=IntKind) :: shelltmp
!
   real (kind=RealKind) :: pri(3), prj(3), postmp(3)
!
   if (.not.Initialized) then
      call ErrorHandler( 'setNeighborEIZ', & 
                         'Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setNeighborEIZ','invalid local atom index',id)
   endif
!  ===================================================================
!  sort neighbors such that for a specific radia the order will be
!  descending on x, y, z components of the positions
!  ===================================================================
!
   do i = 1,NeighborEIZ(id)%NumAtoms
      pri = NeighborEIZ(id)%Position(1:3,i)
      do j = 1,NeighborEIZ(id)%NumAtoms
         prj = NeighborEIZ(id)%Position(1:3,j)
         if ( prj(1) < pri(1) ) then
            postmp = NeighborEIZ(id)%Position(1:3,i)
            ztmp = NeighborEIZ(id)%Z(i)
            ltmp = NeighborEIZ(id)%Lmax(i)
            petmp = NeighborEIZ(id)%ProcIndex(i)
            lindtmp = NeighborEIZ(id)%LocalIndex(i)
            gindtmp = NeighborEIZ(id)%GlobalIndex(i)
            shelltmp = NeighborEIZ(id)%ShellIndex(i)
!
            NeighborEIZ(id)%Position(1:3,i) = NeighborEIZ(id)%Position(1:3,j)
            NeighborEIZ(id)%Z(i) = NeighborEIZ(id)%Z(j)
            NeighborEIZ(id)%Lmax(i) = NeighborEIZ(id)%Lmax(j)
            NeighborEIZ(id)%ProcIndex(i) = NeighborEIZ(id)%ProcIndex(j)
            NeighborEIZ(id)%LocalIndex(i) = NeighborEIZ(id)%LocalIndex(j)
            NeighborEIZ(id)%GlobalIndex(i) = NeighborEIZ(id)%GlobalIndex(j)
            NeighborEIZ(id)%ShellIndex(i) = NeighborEIZ(id)%ShellIndex(j)
!
            NeighborEIZ(id)%Position(1:3,j) = postmp
            NeighborEIZ(id)%Z(j) = ztmp
            NeighborEIZ(id)%Lmax(j) = ltmp
            NeighborEIZ(id)%ProcIndex(j) = petmp
            NeighborEIZ(id)%LocalIndex(j) = lindtmp
            NeighborEIZ(id)%GlobalIndex(j) = gindtmp
            NeighborEIZ(id)%ShellIndex(j) = shelltmp
         endif
         if ( prj(2) < pri(2) .and. pri(1)==prj(1) ) then
            postmp = NeighborEIZ(id)%Position(1:3,i)
            ztmp = NeighborEIZ(id)%Z(i)
            ltmp = NeighborEIZ(id)%Lmax(i)
            petmp = NeighborEIZ(id)%ProcIndex(i)
            lindtmp = NeighborEIZ(id)%LocalIndex(i)
            gindtmp = NeighborEIZ(id)%GlobalIndex(i)
            shelltmp = NeighborEIZ(id)%ShellIndex(i)
!
            NeighborEIZ(id)%Position(1:3,i) = NeighborEIZ(id)%Position(1:3,j)
            NeighborEIZ(id)%Z(i) = NeighborEIZ(id)%Z(j)
            NeighborEIZ(id)%Lmax(i) = NeighborEIZ(id)%Lmax(j)
            NeighborEIZ(id)%ProcIndex(i) = NeighborEIZ(id)%ProcIndex(j)
            NeighborEIZ(id)%LocalIndex(i) = NeighborEIZ(id)%LocalIndex(j)
            NeighborEIZ(id)%GlobalIndex(i) = NeighborEIZ(id)%GlobalIndex(j)
            NeighborEIZ(id)%ShellIndex(i) = NeighborEIZ(id)%ShellIndex(j)
!
            NeighborEIZ(id)%Position(1:3,j) = postmp
            NeighborEIZ(id)%Z(j) = ztmp
            NeighborEIZ(id)%Lmax(j) = ltmp
            NeighborEIZ(id)%ProcIndex(j) = petmp
            NeighborEIZ(id)%LocalIndex(j) = lindtmp
            NeighborEIZ(id)%GlobalIndex(j) = gindtmp
            NeighborEIZ(id)%ShellIndex(j) = shelltmp
         endif
         if ( prj(3) < pri(3) .and. pri(2) == prj(2) .and.             &
              pri(1)==prj(1)  ) then
            postmp = NeighborEIZ(id)%Position(1:3,i)
            ztmp = NeighborEIZ(id)%Z(i)
            ltmp = NeighborEIZ(id)%Lmax(i)
            petmp = NeighborEIZ(id)%ProcIndex(i)
            lindtmp = NeighborEIZ(id)%LocalIndex(i)
            gindtmp = NeighborEIZ(id)%GlobalIndex(i)
            shelltmp = NeighborEIZ(id)%ShellIndex(i)
!
            NeighborEIZ(id)%Position(1:3,i) = NeighborEIZ(id)%Position(1:3,j)
            NeighborEIZ(id)%Z(i) = NeighborEIZ(id)%Z(j)
            NeighborEIZ(id)%Lmax(i) = NeighborEIZ(id)%Lmax(j)
            NeighborEIZ(id)%ProcIndex(i) = NeighborEIZ(id)%ProcIndex(j)
            NeighborEIZ(id)%LocalIndex(i) = NeighborEIZ(id)%LocalIndex(j)
            NeighborEIZ(id)%GlobalIndex(i) = NeighborEIZ(id)%GlobalIndex(j)
            NeighborEIZ(id)%ShellIndex(i) = NeighborEIZ(id)%ShellIndex(j)
!
            NeighborEIZ(id)%Position(1:3,j) = postmp
            NeighborEIZ(id)%Z(j) = ztmp
            NeighborEIZ(id)%Lmax(j) = ltmp
            NeighborEIZ(id)%ProcIndex(j) = petmp
            NeighborEIZ(id)%LocalIndex(j) = lindtmp
            NeighborEIZ(id)%GlobalIndex(j) = gindtmp
            NeighborEIZ(id)%ShellIndex(j) = shelltmp
         endif
      enddo
   enddo
!
!   call printNeighbor(id)
!
   pnb => NeighborEIZ(id)
!
   end function sortNeighborsEIZ
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumReceives(id) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler( 'getNumReceives', & 
                         'Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getNumReceives','invalid local atom index',id)
   endif
!
   n = Neighbor(id)%NumReceives
!
   end function getNumReceives
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMaxReceives() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler( 'getMaxReceives', & 
                         'Need to initialize the module first')
   endif
!
   n = MaxReceives
!
   end function getMaxReceives
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setMaxReceives(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
!
   if (.not.Initialized) then
      call ErrorHandler( 'setMaxReceives', & 
                         'Need to initialize the module first')
   endif
!
   MaxReceives = n
!
   end subroutine setMaxReceives
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSendTable(id,n,stbl)
!  ===================================================================
   use Atom2ProcModule, only : getAtom2ProcInGroup
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: stbl(*)
!
   integer (kind=IntKind) :: i, j, p
!
   logical :: redundant
!
   if (.not.Initialized) then
      call ErrorHandler( 'setSendTable', & 
                         'Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setSendTable','invalid local atom index',id)
   else if (n < 0) then
      call ErrorHandler('setSendTable','invalid remote atom list size',n)
   endif
!
   if (SendTable(id)%NumRequests > n) then
      deallocate( SendTable(id)%RemoteAtomList )
      deallocate( SendTable(id)%RemoteProcList )
   endif
!
   SendTable(id)%NumRequests = n
   if (n > 0) then
      allocate( SendTable(id)%RemoteAtomList(n), SendTable(id)%RemoteProcList(n) )
      SendTable(id)%RemoteAtomList(1:n) = stbl(1:n)
!
      SendTable(id)%NumComms = 1
      SendTable(id)%RemoteProcList(1) = getAtom2ProcInGroup(stbl(1))
   else
      SendTable(id)%NumComms = 0
   endif
   do i = 2, n
      p = getAtom2ProcInGroup(stbl(i))
      redundant = .false.
      LOOP_j: do j = 1, SendTable(id)%NumComms
         if (SendTable(id)%RemoteProcList(j) == p) then
            redundant = .true.
            exit LOOP_j
         endif
      enddo LOOP_j
      if (.not.redundant) then
         SendTable(id)%NumComms = SendTable(id)%NumComms + 1
         j = SendTable(id)%NumComms
         SendTable(id)%RemoteProcList(j) = p
      endif
   enddo
!
   end subroutine setSendTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupRecvTable(id)
!  ===================================================================
   use Atom2ProcModule, only : getAtom2ProcInGroup
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n, i, j, p
!
   logical :: redundant
!
   n = Neighbor(id)%NumReceives
   RecvTable(id)%NumRequests = n
   if (n > 0) then
      allocate( RecvTable(id)%RemoteAtomList(n) )
      allocate( RecvTable(id)%RemoteProcList(n) )
!
      j = 0
      do i = 1, Neighbor(id)%NumAtoms
         if (Neighbor(id)%ProcIndex(i) /= MyPEinGroup) then
            j = j + 1
            if (j > RecvTable(id)%NumRequests) then
               call ErrorHandler('setupRecvTable','n out of bound',j, &
                                 RecvTable(id)%NumRequests)
            endif
            RecvTable(id)%RemoteAtomList(j) = Neighbor(id)%GlobalIndex(i)
         endif
      enddo
      n = j
!
      RecvTable(id)%NumComms = 0
      do i = 1, n
         p = getAtom2ProcInGroup(RecvTable(id)%RemoteAtomList(i))
         redundant = .false.
         LOOP_j: do j = 1, RecvTable(id)%NumComms
            if (RecvTable(id)%RemoteProcList(j) == p) then
               redundant = .true.
               exit LOOP_j
            endif
         enddo LOOP_j
         if (.not.redundant) then
            RecvTable(id)%NumComms = RecvTable(id)%NumComms + 1
            j = RecvTable(id)%NumComms
            RecvTable(id)%RemoteProcList(j) = p
         endif
      enddo
   endif
!
   end subroutine setupRecvTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSendTable(id) result (ptbl)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   type (CommTableStruct), pointer :: ptbl
!
   if (.not.Initialized) then
      call ErrorHandler( 'getSendTable', & 
                         'Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSendTable','invalid local atom index',id)
   endif
!
   ptbl => SendTable(id)
!
   end function getSendTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRecvTable(id) result (ptbl)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   type (CommTableStruct), pointer :: ptbl
!
   if (.not.Initialized) then
      call ErrorHandler( 'getRecvTable', & 
                         'Need to initialize the module first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getRecvTable','invalid local atom index',id)
   endif
!
   ptbl => RecvTable(id)
!
   end function getRecvTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printCommunicationTable()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id, i
!  
   if (GroupSize == 1) then
      return
   endif
!  
   write(6,'(/,2x,a,i8)')'Send Table of Process #',MyPEinGroup
   write(6,'(  2x,a)')'=============================================='
   write(6,'(  2x,a)')'Local Atom     Num Clients      Client Process'
   write(6,'(  2x,a)')'----------------------------------------------'
   do id = 1, LocalNumAtoms
      if (SendTable(id)%NumComms > 0) then
         write(6,'(6x,i2,12x,i4,12x,i8)')id,SendTable(id)%NumComms,  &
                                         SendTable(id)%RemoteProcList(1)
      endif
      do i = 2, SendTable(id)%NumComms
         write(6,'(36x,i8)')SendTable(id)%RemoteProcList(i)
      enddo
      write(6,'(6x,a)')   '------------------------------------------'
      write(6,'(6x,a,i4)')'Num of Requests from Remote Atoms: ',SendTable(id)%NumRequests
      write(6,'(6x,a)')   '------------------------------------------'
   enddo
   write(6,'(  2x,a)')'=============================================='
!  
   write(6,'(/,2x,a,i8)')'Receive Table of Process #',MyPEinGroup
   write(6,'(  2x,a)')'=============================================='
   write(6,'(  2x,a)')'Local Atom     Num Servers      Server Process'
   write(6,'(  2x,a)')'----------------------------------------------'
   do id =1, LocalNumAtoms
      if (RecvTable(id)%NumComms > 0) then
          write(6,'(6x,i2,12x,i4,12x,i8)')id, RecvTable(id)%NumComms,&
                                          RecvTable(id)%RemoteProcList(1)
      endif
      do i = 2, RecvTable(id)%NumComms
         write(6,'(36x,i8)')RecvTable(id)%RemoteProcList(i)
      enddo
      write(6,'(6x,a)')   '------------------------------------------'
      write(6,'(6x,a,i4)')'Num of Requests upon Remote Atoms: ',RecvTable(id)%NumRequests
      write(6,'(6x,a)')   '------------------------------------------'
   enddo
   write(6,'(  2x,a)')'=============================================='
!
   end subroutine printCommunicationTable
!  ===================================================================
end module NeighborModule
