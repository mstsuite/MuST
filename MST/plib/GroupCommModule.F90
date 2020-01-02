module GroupCommModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
!
   use MPPModule, only : MPP_BUFFER_MEM
!
   implicit none
!
public ::                   &
   initGroupComm,           &
   endGroupComm,            &
   createProcGrid,          &
   destroyProcGrid,         &
   createBoxGroupFromGrid,  &
   createGroupFromGrid,     &
   createDimGroupFromGrid,  &
   destroyGroupFromGrid,    &
   getNumGroups,            &
   getGroupID,              &
   getGroupLabel,           &
   getGroupMember,          &
   getGroupMembers,         &
   getNumClustersOfGroup,   &
   getMyClusterIndex,       & ! 1 <= index <= NumClusters with a Group ID
   getNumPEsInGroup,        &
   getMyPEinGroup,          &
   getGroupCommunicator,    &
   getNumProcGrids,         &
   getProcGridID,           &
   getProcGridDimention,    &
   getProcGridSize,         &
   getMyCoordInProcGrid,    &
   bcastMessageInGroup,     &
   GlobalSumInGroup,        &
   GlobalMaxInGroup,        &
   GlobalMinInGroup,        &
   GlobalCollectInGroup,    &
   openLocalMemoryInGroup,  &
   readRemoteMemoryInGroup, &
   writeRemoteMemoryInGroup,&
   closeLocalMemoryInGroup, &
   syncLocalMemoryInGroup,  &
   syncAllPEsInGroup,       &
   isInSameGroup,           &
   isGroupExisting
!
   interface createProcGrid
      module procedure createProcGrid0, createProcGrid1
   end interface createProcGrid
!
   interface getProcGridSize
      module procedure getProcGridSize0, getProcGridSize1
   end interface getProcGridSize
!
   interface getMyCoordInProcGrid
      module procedure getMyCoordInProcGrid0, getMyCoordInProcGrid1
   end interface getMyCoordInProcGrid
!
   interface getNumClustersOfGroup
      module procedure getNumClustersOfGroup_i, getNumClustersOfGroup_s
   end interface getNumClustersOfGroup
!
   interface getMyClusterIndex
      module procedure getMyClusterIndex_i, getMyClusterIndex_s
   end interface getMyClusterIndex
!
   interface getGroupMember
      module procedure getGroupMember_i, getGroupMember_s
   end interface getGroupMember
!
   interface getGroupMembers
      module procedure getGroupMembers_i, getGroupMembers_s
   end interface getGroupMembers
!
   interface isInSameGroup
      module procedure isInSameGroup_i, isInSameGroup_s
   end interface isInSameGroup
!
   interface bcastMessageInGroup
      module procedure bcastMessageInGroup_s0, bcastMessageInGroup_s1, &
                       bcastMessageInGroup_s2, bcastMessageInGroup_s3
      module procedure bcastMessageInGroup_i0, bcastMessageInGroup_i1, &
                       bcastMessageInGroup_i2, bcastMessageInGroup_i3
      module procedure bcastMessageInGroup_r0, bcastMessageInGroup_r1, &
                       bcastMessageInGroup_r2, bcastMessageInGroup_r3
      module procedure bcastMessageInGroup_c0, bcastMessageInGroup_c1, &
                       bcastMessageInGroup_c2, bcastMessageInGroup_c3
   end interface bcastMessageInGroup
!
   interface GlobalSumInGroup
      module procedure GlobalSumInGroup_i0, GlobalSumInGroup_i1, &
                       GlobalSumInGroup_i2, GlobalSumInGroup_i3, &
                       GlobalSumInGroup_i4
      module procedure GlobalSumInGroup_r0, GlobalSumInGroup_r1, &
                       GlobalSumInGroup_r2, GlobalSumInGroup_r3, &
                       GlobalSumInGroup_r4
      module procedure GlobalSumInGroup_c0, GlobalSumInGroup_c1, &
                       GlobalSumInGroup_c2, GlobalSumInGroup_c3, &
                       GlobalSumInGroup_c4
   end interface GlobalSumInGroup
!
   interface GlobalMaxInGroup
      module procedure GlobalMaxInGroup_i0, GlobalMaxInGroup_i1, &
                       GlobalMaxInGroup_i2
      module procedure GlobalMaxInGroup_r0, GlobalMaxInGroup_r1, &
                       GlobalMaxInGroup_r2
      module procedure GlobalMaxInGroup_c0, GlobalMaxInGroup_c1, &
                       GlobalMaxInGroup_c2
   end interface GlobalMaxInGroup
!
   interface GlobalMinInGroup
      module procedure GlobalMinInGroup_i0, GlobalMinInGroup_i1, &
                       GlobalMinInGroup_i2
      module procedure GlobalMinInGroup_r0, GlobalMinInGroup_r1, &
                       GlobalMinInGroup_r2
      module procedure GlobalMinInGroup_c0, GlobalMinInGroup_c1, &
                       GlobalMinInGroup_c2
   end interface GlobalMinInGroup
!
   interface GlobalCollectInGroup
      module procedure GlobalCollectInGroup_i0, GlobalCollectInGroup_i1, &
                       GlobalCollectInGroup_i2
      module procedure GlobalCollectInGroup_r0, GlobalCollectInGroup_r1, &
                       GlobalCollectInGroup_r2
      module procedure GlobalCollectInGroup_c0, GlobalCollectInGroup_c1, &
                       GlobalCollectInGroup_c2
      module procedure GlobalCollectInGroup_s0, GlobalCollectInGroup_s1, &
                       GlobalCollectInGroup_s2
   end interface GlobalCollectInGroup
!
   interface openLocalMemoryInGroup
      module procedure openLocalMemory_s, openLocalMemory_i,     &
                       openLocalMemory_r, openLocalMemory_c
      module procedure openLocalMemory_s1, openLocalMemory_i1,   &
                       openLocalMemory_r1, openLocalMemory_c1
   end interface
!
   interface readRemoteMemoryInGroup
      module procedure readRemoteMemory_s, readRemoteMemory_i,   &
                       readRemoteMemory_r, readRemoteMemory_c
   end interface
!
   interface writeRemoteMemoryInGroup
      module procedure writeRemoteMemory_s, writeRemoteMemory_i, &
                       writeRemoteMemory_r, writeRemoteMemory_c
   end interface
!
private
!
#ifdef MPI
   include 'mpif.h'
#else
#define MPI_STATUS_SIZE 1
#define MPI_ADDRESS_KIND IntKind
#endif
   integer status(MPI_STATUS_SIZE)
!
   logical :: Initialized = .false.
!
   type CommGroups
      character (len=48) :: Label
      integer (kind=IntKind) :: NumClusters ! number of clusters with group "Label"
      integer (kind=IntKind) :: MyCluster   ! my cluster index in the clusters of group "Label"
      integer (kind=IntKind) :: Index
      integer (kind=IntKind) :: MyRank
      integer (kind=IntKind) :: Size      ! number of processes in my group
      integer (kind=IntKind) :: Communicator
      integer (kind=IntKind), pointer :: GroupMembers(:) ! table of group members
      type (CommGroups), pointer :: next
      type (CommGroups), pointer :: prev
   end type CommGroups
!
   type ProcGrids
      character (len=48) :: Label
      integer (kind=IntKind) :: Index
      integer (kind=IntKind) :: Dim
      integer (kind=IntKind) :: Size
      integer (kind=IntKind), allocatable :: Dsize(:)
      integer (kind=IntKind), allocatable :: MyCoord(:)
      type (CommGroups), pointer :: ptr2cg ! this pointer points to the comm group formed by the grid self.
      type (ProcGrids), pointer :: next
      type (ProcGrids), pointer :: prev
   end type ProcGrids
!
   integer (kind=IntKind) :: NumProcGrids = 0
   integer (kind=IntKind) :: NumGroups = 0
   integer (kind=IntKind) :: Me
   integer (kind=IntKind) :: info
   integer (kind=IntKind) :: Globe_comm
!
   type (ProcGrids), pointer :: head2PG
   type (ProcGrids), pointer :: tail2PG
   type (ProcGrids), pointer :: ptr2PG
!
   type (CommGroups), pointer :: head2CommGr
   type (CommGroups), pointer :: tail2CommGr
   type (CommGroups), pointer :: ptr2CommGr
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initGroupComm()
!  ===================================================================
   use MPPModule, only : MyPE, NumPEs, getCommunicator
!
   implicit none
!
   if (NumPEs < 1) then
      call ErrorHandler('initGroupComm','MPPModule is not initialized')
   endif
!
   Me = MyPE
   Globe_comm = getCommunicator()
   Initialized = .true.
   nullify( head2PG, tail2PG, ptr2PG )
   nullify( head2CommGr, tail2CommGr, ptr2CommGr )
!
   end subroutine initGroupComm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endGroupComm()
!  ===================================================================
   implicit none
!
   Initialized = .false.
!
   if (NumProcGrids == 0) then
      return
   endif
!
   do while (NumGroups > 1)
#ifdef MPI
!     ----------------------------------------------------------------
      call MPI_comm_free(tail2CommGr%Communicator,info)
!     ----------------------------------------------------------------
#endif
      deallocate( tail2CommGr%GroupMembers )
      nullify( tail2CommGr%GroupMembers )
      tail2CommGr => tail2CommGr%prev
      deallocate( tail2CommGr%next )
      nullify( tail2CommGr%next )
      NumGroups = NumGroups - 1
   enddo
!
   if (NumGroups == 1) then
#ifdef MPI
!     ----------------------------------------------------------------
      call MPI_comm_free(head2CommGr%Communicator,info)
!     ----------------------------------------------------------------
#endif
      deallocate( head2CommGr )
   endif
   nullify( head2CommGr, tail2CommGr, ptr2CommGr )
   NumGroups = 0
!
   do while (NumProcGrids > 1)
      deallocate( tail2PG%Dsize, tail2PG%MyCoord )
      tail2PG => tail2PG%prev
      deallocate( tail2PG%next )
      nullify( tail2PG%next )
      NumProcGrids = NumProcGrids - 1
   enddo
!
   deallocate( head2PG%Dsize, head2PG%MyCoord )
   deallocate( head2PG )
   nullify( head2PG, tail2PG, ptr2PG )
   NumProcGrids = 0
!
   end subroutine endGroupComm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createBoxGroupFromGrid(g,box,s)
!  ===================================================================
!
!  Purpose: create subgroups in a shape defined by box from a grid
!           of processors
!
!  g : the grid of processors from which the subgroups will be created from.
!  box(1:n) : the size of the box in each dimension j = 1, 2, ..., n
!
!      In the following example, n = 2, box(1) = 4, and box(2) = 3:
!
!        x    x    x    x    o    o    o    o
!        x    x    x    x    o    o    o    o
!        x    x    x    x    o    o    o    o
!        m    m    m    m    w    w    w    w
!        m    m    m    m    w    w    w    w
!        m    m    m    m    w    w    w    w
!
!      processors "x", "o", "w", and "m" form subgroups of label "s".
!
!  s : the label of the subgroup to be created.
!
!  ===================================================================
   use MPPModule, only : MyPE
   implicit none
!
   character (len=*), intent(in) :: s
   character (len=22), parameter :: sname = 'createBoxGroupFromGrid'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: box(:)
   integer (kind=IntKind) :: i, n, color, key, ibase, jbase
   integer (kind=IntKind), allocatable :: bid(:), NumBlocks(:)
!  
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   endif
!  
!  -------------------------------------------------------------------
   call searchPG(sname,g)     ! find the right grid to work on
!  -------------------------------------------------------------------
   n = size(box(:))
   if (n /= ptr2PG%Dim) then
      call ErrorHandler(sname,'Box dimension is different from the Grid dimension',n,ptr2PG%Dim)
   endif
!
   do i = 1, n                ! check the validity of all dimensions
      if (box(i) < 0 .or. box(i) > ptr2PG%Dsize(i)) then
         call ErrorHandler(sname,'For dim j, the box size is invalid',  &
                           i,box(i))
      else if (mod(ptr2PG%Dsize(i),box(i)) /= 0) then
         call ErrorHandler(sname,'For dim i, the dim size can not be divided by the box size', &
                           i,box(i))
      endif
   enddo
!
   allocate( NumBlocks(1:n), bid(1:n) )
!
!  -------------------------------------------------------------------
   call allocateGroup(s)
!  -------------------------------------------------------------------
!
   color = 0; key = 0
   ibase = 1; jbase = 1
   ptr2CommGr%NumClusters = 1
   do i = 1, n
      NumBlocks(i) = ptr2PG%Dsize(i)/box(i)
      ptr2CommGr%NumClusters = ptr2CommGr%NumClusters * NumBlocks(i)
      bid(i) = floor(real(ptr2PG%MyCoord(i)/box(i),kind=RealKind))   ! starting from 0
      color = color + ibase*bid(i)
      key = key + jbase*mod(ptr2PG%MyCoord(i),box(i))
      ibase = ibase*NumBlocks(i)
      jbase = jbase*box(i)
   enddo
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_comm_split(ptr2PG%ptr2cg%Communicator,color,key,          &
                       ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
   call MPI_comm_size(ptr2CommGr%Communicator,ptr2CommGr%Size,info)
!  -------------------------------------------------------------------
   call MPI_comm_rank(ptr2CommGr%Communicator,ptr2CommGr%MyRank,info)
!  -------------------------------------------------------------------
#else
   ptr2CommGr%Communicator = 0
   ptr2CommGr%Size = 1
   ptr2CommGr%MyRank = 0
#endif
!
   ptr2CommGr%MyCluster = color + 1
!
   deallocate( bid, NumBlocks )
!
   n = ptr2CommGr%Size
   i = ptr2CommGr%MyRank + 1
   allocate( ptr2CommGr%GroupMembers(1:n) )
   ptr2CommGr%GroupMembers(1:n) = 0
   ptr2CommGr%GroupMembers(i) = MyPE
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,ptr2CommGr%GroupMembers,n,         &
                      MPI_INTEGER,MPI_SUM,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine createBoxGroupFromGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createGroupFromGrid(g,step,s)
!  ===================================================================
!
!  Purpose:  create subgroups, each of which is made of slides of 
!            processors cutting the d-axis with a distance "step" between
!            the slides.
!
!  g : the grid of processors from which the subgroups will be created from.
!  step: the increment in each dimension. 
!  s : the label of the subgroup to be created.
!
!  Example 1:
!
!        W    A    Y    W    A    Y    W    A    Y
!        X    O    I    X    O    I    X    O    I
!        W    A    Y    W    A    Y    W    A    Y
!        X    O    I    X    O    I    X    O    I  
!
!      processors "X", "O", "W", "I", "A", "Y", form subgroups, with
!      step(1) = 3, step(2) = 2
!
!  Example 2:
!
!        X4    O4    Y4    X4    O4    Y4    X4    O4    Y4
!        X3    O3    Y3    X3    O3    Y3    X3    O3    Y3
!        X2    O2    Y2    X2    O2    Y2    X2    O2    Y2
!        X1    O1    Y1    X1    O1    Y1    X1    O1    Y1  
!
!      processors "X1", "X2", "X3", "X4", "O1", "O2", "O3", "O4", "Y1",
!      "Y2", "Y3", and "Y4" form subgroups, with
!      step(1) = 3, step(2) = 0 (or step(2) = 4)
!
!  ===================================================================
   use MPPModule, only : MyPE
   implicit none
!
   character (len=*), intent(in) :: s
   character (len=20), parameter :: sname = 'createGroupFromGrid'
!
   integer (kind=IntKind), intent(in) :: g, step(:)
   integer (kind=IntKind) :: i, n, color, key, ibase, jbase
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   endif
!
!  -------------------------------------------------------------------
   call searchPG(sname,g)
!  -------------------------------------------------------------------
   n = size(step(:))
   if (n /= ptr2PG%Dim) then
      call ErrorHandler(sname,'Step dimension is different from Grid dimension',  &
                        n,ptr2PG%Dim)
   endif
!
   do i = 1, n
      if (step(i) < 0 .or. step(i) > ptr2PG%Dsize(i)) then
         call ErrorHandler(sname,'Invalid step size',step(i))
      else if (step(i) > 0) then
         if (mod(ptr2PG%Dsize(i),step(i)) /= 0) then
            call ErrorHandler(sname,'Dimension size is not the multiples of step',&
                              ptr2PG%Dsize(i),step(i))
         endif
      endif
   enddo
!
!  -------------------------------------------------------------------
   call allocateGroup(s)
!  -------------------------------------------------------------------
!
   ptr2CommGr%NumClusters = 1
   color = 0; key = 0
   ibase = 1; jbase = 1
   do i = 1, n
      if (step(i) == 0) then
         ptr2CommGr%NumClusters = ptr2CommGr%NumClusters*ptr2PG%Dsize(i)
         color = color + ibase*ptr2PG%MyCoord(i)
         ibase = ibase*ptr2PG%Dsize(i)
      else
         ptr2CommGr%NumClusters = ptr2CommGr%NumClusters*step(i)
         color = color + ibase*mod(ptr2PG%MyCoord(i),step(i))
         key = key + jbase*floor(real(ptr2PG%MyCoord(i)/step(i),kind=RealKind))
         ibase = ibase*step(i)
         jbase = jbase*ptr2PG%Dsize(i)/step(i)
      endif
   enddo
   ptr2CommGr%MyCluster = color + 1
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_comm_split(ptr2PG%ptr2cg%Communicator,color,key,          &
                       ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
   call MPI_comm_size(ptr2CommGr%Communicator,ptr2CommGr%Size,info)
!  -------------------------------------------------------------------
   call MPI_comm_rank(ptr2CommGr%Communicator,ptr2CommGr%MyRank,info)
!  -------------------------------------------------------------------
#else
   ptr2CommGr%Communicator = 0
   ptr2CommGr%Size = 1
   ptr2CommGr%MyRank = 0
#endif
!
   n = ptr2CommGr%Size
   i = ptr2CommGr%MyRank + 1
   allocate( ptr2CommGr%GroupMembers(1:n) )
   ptr2CommGr%GroupMembers(1:n) = 0
   ptr2CommGr%GroupMembers(i) = MyPE
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,ptr2CommGr%GroupMembers,n,         &
                      MPI_INTEGER,MPI_SUM,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine createGroupFromGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createDimGroupFromGrid(g,d,s)
!  ===================================================================
!
!  Purpose:  create subgroups, each of which is the entire set of processors
!            parallel to the d-axis.
!
!  g : the grid of processors from which the subgroups will be created from.
!  d : the dimension index along which the increment of the slide takes place
!      In the following example:
!
!        x    x    x    x    x    x    x    x    x
!        o    o    o    o    o    o    o    o    o
!        w    w    w    w    w    w    w    w    w
!        m    m    m    m    m    m    m    m    m  --->  d
!
!      processors "x", "o", "w", and "m" form subgroups
!  s : the label of the subgroup to be created.
!
!  ===================================================================
   use MPPModule, only : MyPE
   implicit none
!
   character (len=*), intent(in) :: s
   character (len=22), parameter :: sname = 'createDimGroupFromGrid'
!
   integer (kind=IntKind), intent(in) :: g, d
   integer (kind=IntKind) :: i, n
!
   logical, allocatable :: flag_dims(:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   endif
!
!  -------------------------------------------------------------------
   call searchPG(sname,g)
!  -------------------------------------------------------------------
   if (d > ptr2PG%Dim) then
      call ErrorHandler(sname,'Invalid dimension index: d > dim',d)
   endif
!
   allocate( flag_dims(1:ptr2PG%Dim) )
   flag_dims(1:ptr2PG%Dim) = .false.
   flag_dims(d) = .true.
!
!  -------------------------------------------------------------------
   call allocateGroup(s)
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_cart_sub(ptr2PG%ptr2cg%Communicator,flag_dims,            &
                     ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
   call MPI_comm_size(ptr2CommGr%Communicator,ptr2CommGr%Size,info)
   call MPI_comm_rank(ptr2CommGr%Communicator,ptr2CommGr%MyRank,info)
!  -------------------------------------------------------------------
#else
   ptr2CommGr%Communicator = 0
   ptr2CommGr%Size = 1
   ptr2CommGr%MyRank = 0
#endif
!
   ptr2CommGr%NumClusters = 1
   ptr2CommGr%MyCluster = 1
   n = 1
   do i = 1, ptr2PG%Dim
      if (i /= d) then
         ptr2CommGr%NumClusters = ptr2CommGr%NumClusters*ptr2PG%Dsize(i)
         ptr2CommGr%MyCluster = ptr2CommGr%MyCluster + n*ptr2PG%MyCoord(i)
         n = n*ptr2PG%Dsize(i)
      endif
   enddo
!
!  ===================================================================
!  check size consistency
!  ===================================================================
   if ( ptr2CommGr%Size /= ptr2PG%Dsize(d) ) then
      call ErrorHandler(sname,'Inconsistent Size',ptr2CommGr%Size, ptr2PG%Dsize(d))
   endif
   n = ptr2PG%Size/ptr2CommGr%Size
   if (n /= ptr2CommGr%NumClusters) then
      call ErrorHandler(sname,'Inconsistent number of subsets',n, ptr2CommGr%NumClusters)
   endif
!
   deallocate( flag_dims )
!
   n = ptr2CommGr%Size
   i = ptr2CommGr%MyRank + 1
   allocate( ptr2CommGr%GroupMembers(1:n) )
   ptr2CommGr%GroupMembers(1:n) = 0
   ptr2CommGr%GroupMembers(i) = MyPE
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,ptr2CommGr%GroupMembers,n,         &
                      MPI_INTEGER,MPI_SUM,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine createDimGroupFromGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine destroyProcGrid(g)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'destroyProcGrid'
!
   integer (kind=IntKind), intent(in) :: g
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   endif
!
!  -------------------------------------------------------------------
   call searchPG(sname,g)
!  -------------------------------------------------------------------
!
!
   end subroutine destroyProcGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine destroyGroupFromGrid(i,g)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'destroyGroupFromGrid'
!
   integer (kind=IntKind), intent(in) :: i,g
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (i > NumGroups .or. i < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',i)
   endif
!
!  -------------------------------------------------------------------
   call searchPG(sname,g)
   call searchGroup(sname,i)
!  -------------------------------------------------------------------
!
!
   end subroutine destroyGroupFromGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGroups() result(n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getNumGroups'
!
   integer (kind=IntKind) :: n
!
   n = NumGroups
!
   end function getNumGroups
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGroupID(s) result(n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getGroupID'
!
   character (len=*), intent(in) :: s
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   endif
!
   call searchGroup_s(sname,s)
   n = ptr2CommGr%Index
!
   end function getGroupID
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGroupLabel(i) result(s)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getGroupLabel'
   character (len=48) :: s
!
   integer (kind=IntKind), intent(in) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (i > NumGroups .or. i < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',i)
   endif
!
   call searchGroup(sname,i)
   s = ptr2CommGr%Label
!
   end function getGroupLabel
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumClustersOfGroup_s(s) result(n)
!  ===================================================================
   implicit none
!
   character (len=25), parameter :: sname = 'getNumClustersOfGroup'
!
   integer (kind=IntKind) :: n
!
   character (len=*), intent(in) :: s
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   endif
!
   call searchGroup_s(sname,s)
   n = ptr2CommGr%NumClusters
!
   end function getNumClustersOfGroup_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumClustersOfGroup_i(i) result(n)
!  ===================================================================
   implicit none
!
   character (len=25), parameter :: sname = 'getNumClustersOfGroup'
!
   integer (kind=IntKind) :: n
!
   integer (kind=IntKind), intent(in) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (i > NumGroups .or. i < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',i)
   endif
!
   call searchGroup(sname,i)
   n = ptr2CommGr%NumClusters
!
   end function getNumClustersOfGroup_i
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyClusterIndex_s(s) result(n)
!  ===================================================================
   implicit none
!
   character (len=25), parameter :: sname = 'getMyClusterIndex'
!
   integer (kind=IntKind) :: n
!
   character (len=*), intent(in) :: s
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   endif
!
   call searchGroup_s(sname,s)
   n = ptr2CommGr%MyCluster
!
   end function getMyClusterIndex_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyClusterIndex_i(i) result(n)
!  ===================================================================
   implicit none
!
   character (len=25), parameter :: sname = 'getMyClusterIndex'
!
   integer (kind=IntKind) :: n
!
   integer (kind=IntKind), intent(in) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (i > NumGroups .or. i < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',i)
   endif
!
   call searchGroup(sname,i)
   n = ptr2CommGr%MyCluster
!
   end function getMyClusterIndex_i
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumPEsInGroup(i) result(n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getNumPEsInGroup'
!
   integer (kind=IntKind) :: n
!
   integer (kind=IntKind), intent(in) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (i > NumGroups .or. i < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',i)
   endif
!
   call searchGroup(sname,i)
   n = ptr2CommGr%Size
!
   end function getNumPEsInGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyPEinGroup(i) result(n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getMyPEinGroup'
!
   integer (kind=IntKind) :: n
!
   integer (kind=IntKind), intent(in) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (i > NumGroups .or. i < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',i)
   endif
!
   call searchGroup(sname,i)
   n = ptr2CommGr%MyRank
!
   end function getMyPEinGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGroupCommunicator(g) result(comm)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getGroupCommunicator'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind) :: comm
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
   comm = ptr2CommGr%Communicator
!
   end function getGroupCommunicator
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createProcGrid0(nd,title,proc_dims)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'createProcGrid'
!
   integer (kind=IntKind), intent(in) :: nd
!
   character (len=*), intent(in) :: title
!
   integer (kind=IntKind), intent(in) :: proc_dims(1:nd)
   integer (kind=IntKind) :: i
!
   logical, allocatable :: periodic(:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (nd < 1) then
      call ErrorHandler(sname,'Dimension < 1',nd)
   endif
!
   call allocateGrid(title)
   call allocateGroup(title)
!
   ptr2PG%Dim = nd
   allocate( ptr2PG%Dsize(1:nd), ptr2PG%MyCoord(1:nd) )
   ptr2PG%Size = 1
   do i = 1, nd
      if (proc_dims(i) < 1) then
         call ErrorHandler(sname,'Dimension size < 1',proc_dims(i))
      else
         ptr2PG%Dsize(i) = proc_dims(i)
      endif
      ptr2PG%Size = proc_dims(i)*ptr2PG%Size
   enddo
!
   allocate( periodic(1:nd) )
   periodic(1:nd) = .false.
!
#ifdef MPI
!  ===================================================================
!  We set reorder to be "true" here to allow MPI library to decide the
!  best way to set up the processor grid
!  -------------------------------------------------------------------
   call MPI_cart_create(Globe_comm,nd,proc_dims,periodic,.true.,      &
                        ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
   call MPI_comm_rank(ptr2CommGr%Communicator,ptr2CommGr%MyRank,info)
   call MPI_cart_coords(ptr2CommGr%Communicator,ptr2CommGr%MyRank,nd, &
                        ptr2PG%MyCoord,info)
!  -------------------------------------------------------------------
#else
   ptr2CommGr%Communicator = 0
   ptr2CommGr%MyRank = 0
   ptr2PG%MyCoord = 0
#endif
!
   ptr2CommGr%Size = ptr2PG%Size
   ptr2PG%ptr2cg => ptr2CommGr
!
   deallocate( periodic )
!
   end subroutine createProcGrid0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createProcGrid1(comm,nd,title,proc_dims)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'createProcGrid'
!
   integer (kind=IntKind), intent(in) :: comm, nd
!
   character (len=*), intent(in) :: title
!
   integer (kind=IntKind), intent(in) :: proc_dims(1:nd)
   integer (kind=IntKind) :: i
!
   logical, allocatable :: periodic(:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (nd < 1) then
      call ErrorHandler(sname,'Dimension < 1',nd)
   endif
!
   call allocateGrid(title)
   call allocateGroup(title)
!
   ptr2PG%Dim = nd
   allocate( ptr2PG%Dsize(1:nd), ptr2PG%MyCoord(1:nd) )
   ptr2PG%Size = 1
   do i = 1, nd
      if (proc_dims(i) < 1) then
         call ErrorHandler(sname,'Dimension size < 1',proc_dims(i))
      else
         ptr2PG%Dsize(i) = proc_dims(i)
      endif
      ptr2PG%Size = proc_dims(i)*ptr2PG%Size
   enddo
!
   allocate( periodic(1:nd) )
   periodic(1:nd) = .false.
!
#ifdef MPI
!  ===================================================================
!  We set reorder to be "true" here to allow MPI library to decide the
!  best way to set up the processor grid
!  -------------------------------------------------------------------
   call MPI_cart_create(comm,nd,proc_dims,periodic,.true.,            &
                        ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
   call MPI_comm_rank(ptr2CommGr%Communicator,ptr2CommGr%MyRank,info)
   call MPI_cart_get(ptr2CommGr%Communicator,nd,proc_dims,periodic,   &
                     ptr2PG%MyCoord,info)
!  -------------------------------------------------------------------
#else
   ptr2CommGr%Communicator = 0
   ptr2CommGr%MyRank = 0
   ptr2PG%MyCoord = 0
#endif
!
   ptr2CommGr%Size = ptr2PG%Size
   ptr2PG%ptr2cg => ptr2CommGr
!
   deallocate( periodic )
!
   end subroutine createProcGrid1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumProcGrids() result (n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getNumProcGrids'
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   endif
!
   n = NumProcGrids
!
   end function getNumProcGrids
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcGridID(s) result (n)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: s
   character (len=20), parameter :: sname = 'getProcGridDimention'
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   endif
!
!  -------------------------------------------------------------------
   call searchPG_s(sname,s)
!  -------------------------------------------------------------------
!
   n = ptr2PG%Index
!
   end function getProcGridID
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcGridDimention(g) result (n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getProcGridDimention'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   endif
!
!  -------------------------------------------------------------------
   call searchPG(sname,g)
!  -------------------------------------------------------------------
!
   n = ptr2PG%Dim
!
   end function getProcGridDimention
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcGridSize0(g) result (n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getProcGridSize'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   endif
!
!  -------------------------------------------------------------------
   call searchPG(sname,g)
!  -------------------------------------------------------------------
!
   n = ptr2PG%Size
!
   end function getProcGridSize0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcGridSize1(g,d) result (n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getProcGridSize'
!
   integer (kind=IntKind), intent(in) :: g, d
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   else if (d < 1) then
      call ErrorHandler(sname,'Invalid dimension index: d < 1',d)
   endif
!
   call searchPG(sname,g)
   if (d > ptr2PG%Dim) then
      call ErrorHandler(sname,'Invalid dimension index: d > dim',d)
   endif
!
   n = ptr2PG%Dsize(d)
!
   end function getProcGridSize1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyCoordInProcGrid0(g) result (c)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getMyCoordInProcGrid'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), pointer :: c(:)
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   endif
!
   call searchPG(sname,g)
!
   n = ptr2PG%Dim
   c => ptr2PG%MyCoord(1:n)
!
   end function getMyCoordInProcGrid0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMyCoordInProcGrid1(g,d) result (c)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'getMyCoordInProcGrid'
!
   integer (kind=IntKind), intent(in) :: g, d
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumProcGrids == 0) then
      call ErrorHandler(sname,'NumProcGrids = 0')
   else if (g > NumProcGrids .or. g < 1) then
      call ErrorHandler(sname,'Invalid processor grid index',g)
   else if (d < 1) then
      call ErrorHandler(sname,'Invalid dimension index: d < 1',d)
   endif
!
   call searchPG(sname,g)
   if (d > ptr2PG%Dim) then
      call ErrorHandler(sname,'Invalid dimension index: d > dim',d)
   endif
!
   c = ptr2PG%MyCoord(d)
!
   end function getMyCoordInProcGrid1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_s0(g,s,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, root
!
   character (len=*), intent(inout) :: s
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(s,len(s),MPI_CHARACTER,root,ptr2CommGr%Communicator, &
                  info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_s0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_s1(g,s,n,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n, root
!
   character (len=*), intent(inout) :: s(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(s,len(s(1))*n,MPI_CHARACTER,root,                   &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_s1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_s2(g,s,n1,n2,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, root
!
   character (len=*), intent(inout) :: s(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(s,len(s(1,1))*n1*n2,MPI_CHARACTER,root,             &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_s2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_s3(g,s,n1,n2,n3,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, n3, root
!
   character (len=*), intent(inout) :: s(1:n1,1:n2,1:n3)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(s,len(s(1,1,1))*n1*n2*n3,MPI_CHARACTER,root,        &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_s3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_i0(g,f,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, root
!
   integer (kind=IntKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,1,MPI_INTEGER,root,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_i0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_i1(g,f,n,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n, root
!
   integer (kind=IntKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n,MPI_INTEGER,root,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_i2(g,f,n1,n2,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, root
!
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n1*n2,MPI_INTEGER,root,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_i2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_i3(g,f,n1,n2,n3,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, n3, root
!
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2,1:n3)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n1*n2*n3,MPI_INTEGER,root,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_i3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_r0(g,f,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, root
!
   real (kind=RealKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,1,MPI_DOUBLE_PRECISION,root,                      &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_r0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_r1(g,f,n,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n, root
!
   real (kind=RealKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n,MPI_DOUBLE_PRECISION,root,                      &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_r2(g,f,n1,n2,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, root
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n1*n2,MPI_DOUBLE_PRECISION,root,                  &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_r2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_r3(g,f,n1,n2,n3,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, n3, root
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2,1:n3)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n1*n2*n3,MPI_DOUBLE_PRECISION,root,               &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_r3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_c0(g,f,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, root
!
   complex (kind=CmplxKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,1,MPI_DOUBLE_COMPLEX,root,                        &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_c0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_c1(g,f,n,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n, root
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n,MPI_DOUBLE_COMPLEX,root,                        &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_c1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_c2(g,f,n1,n2,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, root
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n1*n2,MPI_DOUBLE_COMPLEX,root,                    &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_c2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine bcastMessageInGroup_c3(g,f,n1,n2,n3,root)
!  ===================================================================
   implicit none
!
   character (len=22), parameter :: sname = 'bcastMessageInGroup'
!
   integer (kind=IntKind), intent(in) :: g, n1, n2, n3, root
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2,1:n3)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   else if (root < 0 .or. root >= ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid root',root)
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_bcast(f,n1*n2*n3,MPI_DOUBLE_COMPLEX,root,                 &
                  ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine bcastMessageInGroup_c3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_i0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   integer (kind=IntKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_INTEGER,MPI_SUM,           &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalSumInGroup_i0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_i1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   integer (kind=IntKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_INTEGER,MPI_SUM,           &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalSumInGroup_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_i2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2)
   integer (kind=IntKind), allocatable :: a(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2 <= MPP_BUFFER_MEM/4) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_INTEGER,MPI_SUM,    &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2,MPI_INTEGER,MPI_SUM,               &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_i2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_i3(g,f,n1,n2,n3)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g,n1,n2,n3
!
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2,1:n3)
   integer (kind=IntKind), allocatable :: a(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2*n3 <= MPP_BUFFER_MEM/4) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2*n3,MPI_INTEGER,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2,n3))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2*n3,MPI_INTEGER,MPI_SUM,            &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_i3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_i4(g,f,n1,n2,n3,n4)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g,n1,n2,n3,n4
!
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2,1:n3,1:n4)
   integer (kind=IntKind), allocatable :: a(:,:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   else if (n4 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n4 < 1',n4)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2*n3*n4 <= MPP_BUFFER_MEM/4) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2*n3*n4,MPI_INTEGER,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2,n3,n4))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2*n3*n4,MPI_INTEGER,MPI_SUM,            &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_i4
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_r0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   real (kind=RealKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_DOUBLE_PRECISION,MPI_SUM,  &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalSumInGroup_r0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_r1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   real (kind=RealKind), intent(inout) :: f(1:n)
   real (kind=RealKind), allocatable :: a(:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n <= MPP_BUFFER_MEM/8) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_DOUBLE_PRECISION,MPI_SUM,  &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n,MPI_DOUBLE_PRECISION,MPI_SUM,          &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_r2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2)
   real (kind=RealKind), allocatable :: a(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2 <= MPP_BUFFER_MEM/8) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_r2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_r3(g,f,n1,n2,n3)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g,n1,n2,n3
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2,1:n3)
   real (kind=RealKind), allocatable :: a(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2*n3 <= MPP_BUFFER_MEM/8) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2*n3,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2,n3))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2*n3,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_r3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_r4(g,f,n1,n2,n3,n4)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g,n1,n2,n3,n4
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2,1:n3,1:n4)
   real (kind=RealKind), allocatable :: a(:,:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   else if (n4 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n4 < 1',n4)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2*n3*n4 <= MPP_BUFFER_MEM/8) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2*n3*n4,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2,n3,n4))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2*n3*n4,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_r4
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_c0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   complex (kind=CmplxKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_DOUBLE_COMPLEX,MPI_SUM,    &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalSumInGroup_c0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_c1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n)
   complex (kind=CmplxKind), allocatable :: a(:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n <= MPP_BUFFER_MEM/16) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_DOUBLE_COMPLEX,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n,MPI_DOUBLE_COMPLEX,MPI_SUM,            &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_c1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_c2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2)
   complex (kind=CmplxKind), allocatable :: a(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2 <= MPP_BUFFER_MEM/16) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_DOUBLE_COMPLEX,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2,MPI_DOUBLE_COMPLEX,MPI_SUM,        &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_c2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_c3(g,f,n1,n2,n3)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g,n1,n2,n3
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2,1:n3)
   complex (kind=CmplxKind), allocatable :: a(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2*n3 <= MPP_BUFFER_MEM/16) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2*n3,MPI_DOUBLE_COMPLEX,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2,n3))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2*n3,MPI_DOUBLE_COMPLEX,MPI_SUM,     &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_c3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalSumInGroup_c4(g,f,n1,n2,n3,n4)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalSumInGroup'
!
   integer (kind=IntKind), intent(in) :: g,n1,n2,n3,n4
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2,1:n3,1:n4)
   complex (kind=CmplxKind), allocatable :: a(:,:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   else if (n3 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n3 < 1',n3)
   else if (n4 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n4 < 1',n4)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   if (n1*n2*n3*n4 <= MPP_BUFFER_MEM/16) then
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX,MPI_SUM, &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
   else
      allocate(a(n1,n2,n3,n4))
!     ----------------------------------------------------------------
      call MPI_ALLREDUCE(f,a,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX,MPI_SUM,     &
                         ptr2CommGr%Communicator,info)
!     ----------------------------------------------------------------
      f = a
      deallocate(a)
   endif
#endif
!
   end subroutine GlobalSumInGroup_c4
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_i0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   integer (kind=IntKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_INTEGER,MPI_MAX,           &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_i0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_i1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   integer (kind=IntKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_INTEGER,MPI_MAX,           &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_i2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_INTEGER,MPI_MAX,       &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_i2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_r0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   real (kind=RealKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_DOUBLE_PRECISION,MPI_MAX,  &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_r0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_r1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   real (kind=RealKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_DOUBLE_PRECISION,MPI_MAX,  &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_r2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_DOUBLE_PRECISION,MPI_MAX, &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_r2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_c0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   complex (kind=CmplxKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_DOUBLE_COMPLEX,MPI_MAX,    &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_c0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_c1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_DOUBLE_COMPLEX,MPI_MAX,    &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_c1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMaxInGroup_c2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMaxInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_DOUBLE_COMPLEX,MPI_MAX, &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMaxInGroup_c2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_i0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   integer (kind=IntKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_INTEGER,MPI_MIN,           &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_i0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_i1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   integer (kind=IntKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_INTEGER,MPI_MIN,           &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_i2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_INTEGER,MPI_MIN,       &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_i2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_r0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   real (kind=RealKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_DOUBLE_PRECISION,MPI_MIN,  &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_r0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_r1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   real (kind=RealKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_DOUBLE_PRECISION,MPI_MIN,  &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_r2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_DOUBLE_PRECISION,MPI_MIN, &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_r2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_c0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
!
   complex (kind=CmplxKind), intent(inout) :: f
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,1,MPI_DOUBLE_COMPLEX,MPI_MIN,    &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_c0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_c1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n,MPI_DOUBLE_COMPLEX,MPI_MIN,    &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_c1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalMinInGroup_c2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalMinInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_ALLREDUCE(MPI_IN_PLACE,f,n1*n2,MPI_DOUBLE_COMPLEX,MPI_MIN, &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalMinInGroup_c2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_s0(g,a)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   character (len=*), intent(inout) :: a(*)
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(a(i),len(a(1)),                                 &
   call MPI_ALLGATHER(MPI_IN_PLACE,len(a(1)),                         &
                      MPI_CHARACTER,a,len(a(1)),                      &
                      MPI_CHARACTER,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_s0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_s1(g,a,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   character (len=*), intent(inout) :: a(n,*)
!
   integer (kind=IntKind), intent(in) :: g, n
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(a(1,i),len(a(1,1)*n,                            &
   call MPI_ALLGATHER(MPI_IN_PLACE,len(a(1,1))*n,                     &
                      MPI_CHARACTER,a,len(a(1,1))*n,                  &
                      MPI_CHARACTER,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_s1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_s2(g,a,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   character (len=*), intent(inout) :: a(n1,n2,*)
!
   integer (kind=IntKind), intent(in) :: g, n1, n2
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(a(1,1,1),len(a(1,1,1))*n1*n2                    &
   call MPI_ALLGATHER(MPI_IN_PLACE,len(a(1,1,1))*n1*n2,               &
                      MPI_CHARACTER,a,len(a(1,1,1))*n1*n2,            &
                      MPI_CHARACTER,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_s2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_i0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(inout) :: f(*)
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(i),1,MPI_INTEGER,f,1,MPI_INTEGER,             &
   call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,f,1,MPI_INTEGER,     &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_i0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_i1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(inout) :: f(1:n,*)
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(1,i),1,MPI_INTEGER,f,n,MPI_INTEGER,           &
   call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,f,n,MPI_INTEGER,     &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_i2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind), intent(inout) :: f(1:n1,1:n2,*)
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(1,1,i),n1*n2,MPI_INTEGER,f,n1*n2,MPI_INTEGER, &
   call MPI_ALLGATHER(MPI_IN_PLACE,n1*n2,MPI_INTEGER,f,n1*n2,MPI_INTEGER, &
                      ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_i2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_r0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(inout) :: f(*)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(i),1,MPI_DOUBLE_PRECISION,f,1,                &
   call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,f,1,        &
                      MPI_DOUBLE_PRECISION,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_r0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_r1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(inout) :: f(1:n,*)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(1,i),n,MPI_DOUBLE_PRECISION,f,n,              &
   call MPI_ALLGATHER(MPI_IN_PLACE,n,MPI_DOUBLE_PRECISION,f,n,        &
                      MPI_DOUBLE_PRECISION,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_r2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(inout) :: f(1:n1,1:n2,*)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(1,1,i),n1*n2,MPI_DOUBLE_PRECISION,f,n1*n2,    &
   call MPI_ALLGATHER(MPI_IN_PLACE,n1*n2,MPI_DOUBLE_PRECISION,f,n1*n2,&
                      MPI_DOUBLE_PRECISION,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_r2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_c0(g,f)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), intent(inout) :: f(*)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(i),1,MPI_DOUBLE_COMPLEX,f,1,                  &
   call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_DOUBLE_COMPLEX,f,1,          &
                      MPI_DOUBLE_COMPLEX,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_c0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_c1(g,f,n)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n,*)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n < 1',n)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(1,i),n,MPI_DOUBLE_COMPLEX,f,n,                &
   call MPI_ALLGATHER(MPI_IN_PLACE,n,MPI_DOUBLE_COMPLEX,f,n,          &
                      MPI_DOUBLE_COMPLEX,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_c1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GlobalCollectInGroup_c2(g,f,n1,n2)
!  ===================================================================
   implicit none
!
   character (len=20), parameter :: sname = 'GlobalCollectInGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), intent(inout) :: f(1:n1,1:n2,*)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (n1 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n1 < 1',n1)
   else if (n2 < 1) then
      call ErrorHandler(sname,'Invalid array dimension: n2 < 1',n2)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
   i = ptr2CommGr%MyRank + 1
#ifdef MPI
!  -------------------------------------------------------------------
!  call MPI_ALLGATHER(f(1,1,i),n1*n2,MPI_DOUBLE_COMPLEX,f,n1*n2,      &
   call MPI_ALLGATHER(MPI_IN_PLACE,n1*n2,MPI_DOUBLE_COMPLEX,f,n1*n2,  &
                      MPI_DOUBLE_COMPLEX,ptr2CommGr%Communicator,info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine GlobalCollectInGroup_c2
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_s(g,a,n,accessID)
!  ===================================================================
   use MPPModule, only : CHARACTER_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   character (len=*), intent(in) :: a(n)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = len(a(1))*n*CHARACTER_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, CHARACTER_BYTES, MPI_INFO_NULL,          &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_s
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_s1(g,a,n1,n2,accessID)
!  ===================================================================
   use MPPModule, only : CHARACTER_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   character (len=*), intent(in) :: a(n1,n2)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = len(a(1,1))*n1*n2*CHARACTER_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, CHARACTER_BYTES, MPI_INFO_NULL,          &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_s1
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_i(g,a,n,accessID)
!  ===================================================================
   use MPPModule, only : INTEGER_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   integer (kind=IntKind), intent(in) :: a(n)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = n*INTEGER_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, INTEGER_BYTES, MPI_INFO_NULL,            &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_i
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_i1(g,a,n1,n2,accessID)
!  ===================================================================
   use MPPModule, only : INTEGER_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   integer (kind=IntKind), intent(in) :: a(n1,n2)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = n1*n2*INTEGER_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, INTEGER_BYTES, MPI_INFO_NULL,            &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_i1
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_r(g,a,n,accessID)
!  ===================================================================
   use MPPModule, only : REAL_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   real (kind=RealKind), intent(in) :: a(n)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = n*REAL_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, REAL_BYTES, MPI_INFO_NULL,               &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_r
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_r1(g,a,n1,n2,accessID)
!  ===================================================================
   use MPPModule, only : REAL_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   real (kind=RealKind), intent(in) :: a(n1,n2)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = n1*n2*REAL_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, REAL_BYTES, MPI_INFO_NULL,               &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_r1
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_c(g,a,n,accessID)
!  ===================================================================
   use MPPModule, only : COMPLEX_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   complex (kind=CmplxKind), intent(in) :: a(n)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = n*COMPLEX_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, COMPLEX_BYTES, MPI_INFO_NULL,            &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_c
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine openLocalMemory_c1(g,a,n1,n2,accessID)
!  ===================================================================
   use MPPModule, only : COMPLEX_BYTES
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind), intent(out) :: accessID
   integer (kind=MPI_ADDRESS_KIND) :: m
!
   complex (kind=CmplxKind), intent(in) :: a(n1,n2)
!
   character (len=25), parameter :: sname = 'openLocalMemoryInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   m = n1*n2*COMPLEX_BYTES
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_create(a, m, COMPLEX_BYTES, MPI_INFO_NULL,            &
                       ptr2CommGr%Communicator, accessID, info)
!  -------------------------------------------------------------------
#else
   accessID = 0
#endif
!
   end subroutine openLocalMemory_c1
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine readRemoteMemory_s(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: m
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   character (len=*), intent(out) :: a(n)
!
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
   m = n*len(a(1))
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_get(a, m, MPI_CHARACTER, proc, offset,                    &
                m, MPI_CHARACTER, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!       call MPI_win_unlock(proc, accessID, info)
!! #endif
#else
   a = ' '
#endif
!
   end subroutine readRemoteMemory_s
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine readRemoteMemory_i(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   integer (kind=IntKind), intent(out) :: a(n)
!
#ifdef MPI
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
!  -------------------------------------------------------------------
   call MPI_get(a, n, MPI_INTEGER, proc, offset,                      &
                n, MPI_INTEGER, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!       call MPI_win_unlock(proc, accessID, info)
!! #endif
#else
   a = 0
#endif
!
   end subroutine readRemoteMemory_i
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine readRemoteMemory_r(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   real (kind=RealKind), intent(out) :: a(n)
!
#ifdef MPI
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
!  -------------------------------------------------------------------
   call MPI_get(a, n, MPI_DOUBLE_PRECISION, proc, offset,             &
                n, MPI_DOUBLE_PRECISION, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!      call MPI_win_unlock(proc, accessID, info)
!! #endif
#else
   a = 0.0d0
#endif
!
   end subroutine readRemoteMemory_r
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine readRemoteMemory_c(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   complex (kind=CmplxKind), intent(out) :: a(n)
!
#ifdef MPI
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
!  -------------------------------------------------------------------
   call MPI_get(a, n, MPI_DOUBLE_COMPLEX, proc, offset,               &
                n, MPI_DOUBLE_COMPLEX, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!       call MPI_win_unlock(proc, accessID, info)
!! #endif
#else
   a = (0.0d0,0.0d0)
#endif
!
   end subroutine readRemoteMemory_c
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeRemoteMemory_s(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: m
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   character (len=*), intent(in) :: a(n)
!
#ifdef MPI
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
   m = n*len(a(1))
!  -------------------------------------------------------------------
   call MPI_put(a, m, MPI_CHARACTER, proc, offset,                    &
                m, MPI_CHARACTER, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!       call MPI_win_unlock(proc, accessID, info)
!! #endif
#endif
!
   end subroutine writeRemoteMemory_s
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeRemoteMemory_i(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   integer (kind=IntKind), intent(in) :: a(n)
!
#ifdef MPI
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
!  -------------------------------------------------------------------
   call MPI_put(a, n, MPI_INTEGER, proc, offset,                      &
                n, MPI_INTEGER, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!       call MPI_win_unlock(proc, accessID, info)
!! #endif
#endif
!
   end subroutine writeRemoteMemory_i
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeRemoteMemory_r(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   real (kind=RealKind), intent(in) :: a(n)
!
#ifdef MPI
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
!  -------------------------------------------------------------------
   call MPI_put(a, n, MPI_DOUBLE_PRECISION, proc, offset,             &
                n, MPI_DOUBLE_PRECISION, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!       call MPI_win_unlock(proc, accessID, info)
!! #endif
#endif
!
   end subroutine writeRemoteMemory_r
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeRemoteMemory_c(a,n,proc,accessID,index)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: accessID
   integer (kind=IntKind), intent(in) :: proc
   integer (kind=IntKind), intent(in) :: index
   integer (kind=MPI_ADDRESS_KIND) :: offset
!
   complex (kind=CmplxKind), intent(in) :: a(n)
!
#ifdef MPI
!! #ifndef No_MPI_LOCK
!!       call MPI_win_lock(MPI_LOCK_SHARED, proc, 0, accessID, info)
!! #endif
   offset = index-1
!  -------------------------------------------------------------------
   call MPI_put(a, n, MPI_DOUBLE_COMPLEX, proc, offset,               &
                n, MPI_DOUBLE_COMPLEX, accessID, info)
!  -------------------------------------------------------------------
!! #ifndef No_MPI_LOCK
!!       call MPI_win_unlock(proc, accessID, info)
!! #endif
#endif
!
   end subroutine writeRemoteMemory_c
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine closeLocalMemoryInGroup(AccessID)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: AccessID
!
#ifdef MPI
!  -------------------------------------------------------------------
   call MPI_win_free(AccessID, info)
!  -------------------------------------------------------------------
#endif
!
   end subroutine closeLocalMemoryInGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine syncLocalMemoryInGroup(AccessID)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: AccessID
!
#ifdef MPI
!! #ifdef No_MPI_LOCK
!  -------------------------------------------------------------------
   call MPI_WIN_FENCE(0, accessID, info)
!  -------------------------------------------------------------------
!! #endif
#endif
!
   end subroutine syncLocalMemoryInGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGroupMember_s(s,i) result(p)
!  ===================================================================
!
!  Return the member processer ID in MPI_COMM_WORLD
!
!  ===================================================================
   implicit none
!
   character (len=15), parameter :: sname = 'getGroupMember'
!
   character (len=*), intent(in) :: s
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: p
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   endif
!
   call searchGroup_s(sname,s)
!
   if (i < 1 .or. i > ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid member index',i)
   endif
!
   p = ptr2CommGr%GroupMembers(i)
!
   end function getGroupMember_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGroupMember_i(g,i) result(p)
!  ===================================================================
!
!  Return the member processer ID in MPI_COMM_WORLD
!
!  ===================================================================
   implicit none
!
   character (len=15), parameter :: sname = 'getGroupMember'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: p
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if (i < 1 .or. i > ptr2CommGr%Size) then
      call ErrorHandler(sname,'Invalid member index',i)
   endif
!
   p = ptr2CommGr%GroupMembers(i)
!
   end function getGroupMember_i
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGroupMembers_s(s) result(p)
!  ===================================================================
!
!  Return the member processer ID in MPI_COMM_WORLD
!
!  ===================================================================
   implicit none
!
   character (len=15), parameter :: sname = 'getGroupMember'
!
   character (len=*), intent(in) :: s
!
   integer (kind=IntKind), pointer :: p(:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   endif
!
   call searchGroup_s(sname,s)
!
   p => ptr2CommGr%GroupMembers(1:ptr2CommGr%Size)
!
   end function getGroupMembers_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGroupMembers_i(g) result(p)
!  ===================================================================
!
!  Return the member processer ID in MPI_COMM_WORLD
!
!  ===================================================================
   implicit none
!
   character (len=15), parameter :: sname = 'getGroupMember'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), pointer :: p(:)
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   p => ptr2CommGr%GroupMembers(1:ptr2CommGr%Size)
!
   end function getGroupMembers_i
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInSameGroup_s(s,p) result(y)
!  ===================================================================
!
!  Check if processor p is in the same group of g as MyPE
!
!  ===================================================================
   use MPPModule, only : MyPE, NumPEs
   implicit none
!
   character (len=15), parameter :: sname = 'isInSameGroup'
!
   character (len=*), intent(in) :: s
!
   integer (kind=IntKind), intent(in) :: p
   integer (kind=IntKind) :: i, n
!
   logical :: y
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (p < 0 .or. p >= NumPEs) then
      call ErrorHandler(sname,'Invalid processor index',p)
   endif
!
   y = .false.
!
   if (p == MyPE) then
      y = .true.
   else
      call searchGroup_s(sname,s)
      n = ptr2CommGr%Size
      LOOP_i: do i = 1, n
         if (ptr2CommGr%GroupMembers(i) == p) then
            y = .true.
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   end function isInSameGroup_s
!  ===================================================================
!
!  *******************************************************************
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInSameGroup_i(g,p) result(y)
!  ===================================================================
!
!  Check if processor p is in the same group of g as MyPE
!
!  ===================================================================
   use MPPModule, only : MyPE, NumPEs
   implicit none
!
   character (len=15), parameter :: sname = 'isInSameGroup'
!
   integer (kind=IntKind), intent(in) :: g
   integer (kind=IntKind), intent(in) :: p
   integer (kind=IntKind) :: i, n
!
   logical :: y
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   else if (p < 0 .or. p >= NumPEs) then
      call ErrorHandler(sname,'Invalid processor index',p)
   endif
!
   y = .false.
!
   if (p == MyPE) then
      y = .true.
   else
      call searchGroup(sname,g)
      n = ptr2CommGr%Size
      LOOP_i: do i = 1, n
         if (ptr2CommGr%GroupMembers(i) == p) then
            y = .true.
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   end function isInSameGroup_i
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine syncAllPEsInGroup(g)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: g
!
   character (len=20), parameter :: sname = 'syncAllPEsInGroup'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'GroupCommModule is not initialized')
   else if (NumGroups == 0) then
      call ErrorHandler(sname,'NumGroups = 0')
   else if (g > NumGroups .or. g < 1) then
      call ErrorHandler(sname,'Invalid MPI group index',g)
   endif
!
   call searchGroup(sname,g)
!
   if ( ptr2CommGr%Size < 2) then
      return
   endif
!
#ifdef MPI
   call MPI_barrier(ptr2CommGr%Communicator,info)
#endif
!
   end subroutine syncAllPEsInGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine searchGroup_s(sname,s)
!  ===================================================================
   character (len=*), intent(in) :: sname
!
   character (len=*), intent(in) :: s
   integer (kind=IntKind) :: n
!
   logical :: found = .false.
!
   interface
      function nocaseCompare(s1,s2) result(t)
         logical :: t
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
      end function nocaseCompare
   end interface
!
   if (.not.associated(ptr2CommGr)) then
      ptr2CommGr => head2CommGr
   endif
!
   found = .false.
   if (.not.nocaseCompare(ptr2CommGr%Label,s)) then
      ptr2CommGr => head2CommGr
      LOOP_n: do n = 1, NumGroups
         if (nocaseCompare(ptr2CommGr%Label,s)) then
            found = .true.
            exit LOOP_n
         endif
         ptr2CommGr => ptr2CommGr%next
      enddo LOOP_n
      if (.not.found) then
         call ErrorHandler(sname,'MPI group does not contain the label',s)
      endif
   endif
!
   end subroutine searchGroup_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine searchGroup(sname,i)
!  ===================================================================
   character (len=*), intent(in) :: sname
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: n
!
   logical :: found = .false.
!
   if (.not.associated(ptr2CommGr)) then
      ptr2CommGr => head2CommGr
   endif
!
   found = .false.
   if (ptr2CommGr%Index /= i) then
      ptr2CommGr => head2CommGr
      LOOP_n: do n = 1, NumGroups
         if (ptr2CommGr%Index == i) then
            found = .true.
            exit LOOP_n
         endif
         ptr2CommGr => ptr2CommGr%next
      enddo LOOP_n
      if (.not.found) then
         call ErrorHandler(sname,'Corrupted MPI group data record')
      endif
   endif
!
   end subroutine searchGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine searchPG_s(sname,s)
!  ===================================================================
   character (len=*), intent(in) :: sname
!
   character (len=*), intent(in) :: s
!
   integer (kind=IntKind) :: n
!
   logical :: found = .false.
!
   interface
      function nocaseCompare(s1,s2) result(t)
         logical :: t
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
      end function nocaseCompare
   end interface
!
   if (.not.nocaseCompare(ptr2PG%Label,s)) then
      ptr2PG => head2PG
      LOOP_n: do n = 1, NumProcGrids
         if (nocaseCompare(ptr2PG%Label,s)) then
            found = .true.
            exit LOOP_n
         endif
         ptr2PG => ptr2PG%next
      enddo LOOP_n
      if (.not.found) then
         call ErrorHandler(sname,'Processor grid does not match the label',s)
      endif
   endif
!
   end subroutine searchPG_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine searchPG(sname,i)
!  ===================================================================
   character (len=*), intent(in) :: sname
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: n
!
   logical :: found = .false.
!
   if (ptr2PG%Index /= i) then
      ptr2PG => head2PG
      LOOP_n: do n = 1, NumProcGrids
         if (ptr2PG%Index == i) then
            found = .true.
            exit LOOP_n
         endif
         ptr2PG => ptr2PG%next
      enddo LOOP_n
      if (.not.found) then
         call ErrorHandler(sname,'Corrupted processor grid data record')
      endif
   endif
!
   end subroutine searchPG
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateGroup(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: s
!
   if (NumGroups == 0) then
      allocate( head2CommGr )
      nullify( head2CommGr%prev )
      ptr2CommGr => head2CommGr
   else
      allocate( tail2CommGr%next )
      ptr2CommGr => tail2CommGr%next
      ptr2CommGr%prev => tail2CommGr
   endif
!
   nullify( ptr2CommGr%next )
   tail2CommGr => ptr2CommGr
   NumGroups = NumGroups + 1
   ptr2CommGr%Index = NumGroups
   ptr2CommGr%Label = s
!
   end subroutine allocateGroup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateGrid(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: s
!
   if (NumProcGrids == 0) then
      allocate( head2PG )
      nullify( head2PG%prev )
      ptr2PG => head2PG
   else
      allocate( tail2PG%next )
      ptr2PG => tail2PG%next
      ptr2PG%prev => tail2PG
   endif
!
   nullify( ptr2PG%next )
   tail2PG => ptr2PG
   NumProcGrids = NumProcGrids + 1
   ptr2PG%Index = NumProcGrids
   ptr2PG%Label = s
!
   end subroutine allocateGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isGroupExisting(s) result(y)
!  ===================================================================
   character (len=*), intent(in) :: s
   integer (kind=IntKind) :: n
!
   logical :: y
!
   interface
      function nocaseCompare(s1,s2) result(t)
         logical :: t
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
      end function nocaseCompare
   end interface
!
   y = .false.
   ptr2CommGr => head2CommGr
   LOOP_n: do n = 1, NumGroups
      if (nocaseCompare(ptr2CommGr%Label,s)) then
         y = .true.
         exit LOOP_n
      endif
      ptr2CommGr => ptr2CommGr%next
   enddo LOOP_n
!
   end function isGroupExisting
!  ===================================================================
end module GroupCommModule
