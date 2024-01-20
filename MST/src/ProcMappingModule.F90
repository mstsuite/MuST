!  *******************************************************************
!  * Note:                                                           *
!  * After calling initProcMapping, one needs to call                *
!  * createParallelization to realize the actual processes mapping   *
!  *                                                                 *
!  *                                                                 *
!  *******************************************************************
module ProcMappingModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
!
   implicit none
!
public :: initProcMapping,  &
          endProcMapping,   &
          resetParallelParams,   &
          createParallelization, &
          isEnergyOnMyProc, &
          getNumEsOnMyProc, &
          getEnergyIndex,   &
          getNumRedundantEsOnMyProc, &
          getProcWithEnergyIndex, &     ! returns proc ID in a MPI group
          isKPointOnMyProc, &
          getNumKsOnMyProc, &
          getKPointIndex,   &
          getNumRedundantKsOnMyProc, &
          getProcWithKPointIndex, &     ! returns proc ID in a MPI group
          isAtomOnMyProc,         &
          distributeX,               &
          getNumXsOnMyProc,          &
          getXIndex,                 &
          getXCommID,                &
          getNumRedundantXsOnMyProc, &
          getProcWithXIndex,         &  ! returns proc ID in a MPI group
          isXindexOnMyProc
!
   interface initProcMapping
      module procedure initProcMapping_0, initProcMapping_1
   end interface initProcMapping
!
private
   integer (kind=IntKind) :: NumEs = 0
   integer (kind=IntKind) :: NumAtoms = 0
   integer (kind=IntKind) :: NumKs = 0
!
   integer (kind=IntKind) :: NumEsPerBox = 0
   integer (kind=IntKind) :: NumRedundantEsPerBox = 0
   integer (kind=IntKind) :: NumAtomsPerProc = 0
   integer (kind=IntKind) :: NumKsPerBox = 0
   integer (kind=IntKind) :: NumRedundantKsPerBox = 0
   integer (kind=IntKind) :: NumBoxes = 0
   integer (kind=IntKind) :: MyBoxIndex4E = 0
   integer (kind=IntKind) :: MyBoxIndex4K = 0
   integer (kind=IntKind) :: MyIndexInBox = 0
!
   integer (kind=IntKind), allocatable, target :: EsOnMyProc(:)
   integer (kind=IntKind), allocatable, target :: KsOnMyProc(:)
   integer (kind=IntKind), allocatable, target :: AsOnMyProc(:)
   integer (kind=IntKind), allocatable, target :: E2Proc(:)
   integer (kind=IntKind), allocatable, target :: K2Proc(:)
   integer (kind=IntKind), allocatable, target :: A2Proc(:)
!
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: MaxAtomsPerProc = 0
!
   real (kind=RealKind) :: Bravais(3,3)
!
   character (len=21) :: stop_routine
!
   logical :: isFullPotential = .false.
   logical :: Created = .false.
   logical :: Initialized = .false.
!
!  ===================================================================
!  MappingStruct defines a mapping information table which is a linked
!  list of each physical parameter X (atom, energy, k-points, etc) with
!  N points to processor mapping 
!  ===================================================================
   type MappingStruct
      character (len=80) :: ParamKey
      integer (kind=IntKind) :: CommGroupID
      integer (kind=IntKind) :: NumParams                  ! = N, for which X(i) has N points with i = 1, 2, ...,N
      integer (kind=IntKind) :: NumParamsOnMyProc          ! = number of i's mapped on my processor
      integer (kind=IntKind) :: NumRedundantParamsOnMyProc ! = number of redundant i's on my processor
      integer (kind=IntKind), pointer :: ParamOnMyProc(:)  ! = index i for a given local index q = 1, 2, ..., NumParamsOnMyProc
      integer (kind=IntKind), pointer :: ParamToProc(:)    ! = processor ID for a given i index
      type (MappingStruct), pointer :: next
   end type MappingStruct
!
   integer (kind=IntKind) :: NumMappingItems
   type (MappingStruct), pointer :: MappingTable
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initProcMapping_0(ifp,istop,iprint,maxp)
!  ===================================================================
   use SystemModule, only : getNumAtoms, getBravaisLattice
!
   implicit none
!
   character (len=*), intent(in) :: istop
   logical, intent(in) :: ifp
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind), intent(in), optional :: maxp
!
   isFullPotential = ifp
   Bravais(1:3,1:3) = getBravaisLattice()
   NumAtoms = getNumAtoms()
   NumEs = 0
   NumKs = 0
!
   if ( present(maxp) ) then
      if (maxp < 1) then
         if (isFullPotential) then
            MaxAtomsPerProc = 8
         else
            MaxAtomsPerProc = 12
         endif
      else
         MaxAtomsPerProc = maxp
      endif
   else
      if (isFullPotential) then
         MaxAtomsPerProc = 8
      else
         MaxAtomsPerProc = 12
      endif
   endif
!
!  -------------------------------------------------------------------
   call checkResources()
!  -------------------------------------------------------------------
!
   NumMappingItems = 0
   stop_routine = istop
   print_level = iprint
!
   Initialized = .true.
!
   end subroutine initProcMapping_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initProcMapping_1(na,ne,nk,ifp,istop,iprint,maxp)
!  ===================================================================
   use SystemModule, only : getBravaisLattice
!
   implicit none
!
   character (len=*), intent(in) :: istop
   logical, intent(in) :: ifp
   integer (kind=IntKind), intent(in) :: na, ne, nk
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind), intent(in), optional :: maxp
!
   isFullPotential = ifp
   Bravais(1:3,1:3) = getBravaisLattice()
   NumAtoms = na
   NumEs = ne
   NumKs = nk
!
   if ( present(maxp) ) then
      if (maxp < 1) then
         if (isFullPotential) then
            MaxAtomsPerProc = 8
         else
            MaxAtomsPerProc = 12
         endif
      else
         MaxAtomsPerProc = maxp
      endif
   else
      if (isFullPotential) then
         MaxAtomsPerProc = 8
      else
         MaxAtomsPerProc = 12
      endif
   endif
!
!  -------------------------------------------------------------------
   call checkResources()
!  -------------------------------------------------------------------
!
   NumMappingItems = 0
   stop_routine = istop
   print_level = iprint
!
   Initialized = .true.
!
   end subroutine initProcMapping_1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endProcMapping()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i
!
   type (MappingStruct), pointer :: p_item
!
   do i = 1, NumMappingItems
      p_item => MappingTable%next
      deallocate(MappingTable)
      MappingTable => p_item
   enddo 
   nullify(p_item)
   nullify(MappingTable)
!
   if (Created) then
      deallocate(EsOnMyProc, KsOnMyProc, AsOnMyProc)
      deallocate(E2Proc, K2Proc, A2Proc)
      Created = .false.
   endif
!
   NumEs = 0
   NumAtoms = 0
   NumKs = 0
   NumEsPerBox = 0
   NumRedundantEsPerBox = 0
   NumAtomsPerProc = 0
   NumKsPerBox = 0
   NumRedundantKsPerBox = 0
   NumMappingItems = 0
   isFullPotential = .false.
   Initialized = .false.
!
   end subroutine endProcMapping
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkResources()
!  ===================================================================
   use MPPModule, only : MyPE, NumPEs
!
   use CmdLineOptionModule, only : getCmdLineOptionValue
!
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = getCmdLineOptionValue('Number of Atoms Per Process',NumAtomsPerProc)
!
   if (NumAtomsPerProc == 0) then
      return
   else if (MyPE == 0) then
      write(6,'(a,i5)')'Number of atoms per process = ',NumAtomsPerProc
   endif
!
   if (NumAtomsPerProc > MaxAtomsPerProc) then
      call ErrorHandler('initProcMapping','NumAtomsPerProc > MaxAtomsPerProc', &
                        NumAtomsPerProc, MaxAtomsPerProc)
   else if (NumAtomsPerProc > NumAtoms) then
      call ErrorHandler('initProcMapping','NumAtomsPerProc > NumAtoms', &
                        NumAtomsPerProc, NumAtoms)
   else if (NumAtomsPerProc > 0) then
      if (mod(NumAtoms,NumAtomsPerProc) /= 0) then
         call ErrorHandler('initProcMapping',                         &
                           'NumAtoms is not divisible by NumAtomsPerProc', &
                           NumAtoms,NumAtomsPerProc)
      else if (mod(NumPEs*NumAtomsPerProc,NumAtoms) /= 0) then
         call ErrorHandler('initProcMapping',                         &
                           'NumPEs*NumAtomsPerProc is not divisible by NumAtoms', &
                           NumPEs*NumAtomsPerProc,NumAtoms)
      endif
   endif
!
   end subroutine checkResources
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetParallelParams(na,nk,ne)
!  ===================================================================
   use GroupCommModule, only : getProcGridID, destroyProcGrid
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, nk, ne
!
   character (len=20), parameter :: sname = 'resetParallelParams'
!
   integer (kind=IntKind) :: i
!
   if (na < 1) then
      call ErrorHandler(sname,'Number of Atoms < 1',na)
   else if (ne < 1) then
      call ErrorHandler(sname,'Number of Energies < 1',ne)
   else if (nk < 1) then
      call ErrorHandler(sname,'Number of KPoints < 1',nk)
   else if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   endif
!
   NumAtoms = na
   NumKs = nk
   NumEs = ne
!
   i = getProcGridID('3-D Proc Grid')
   call destroyProcGrid(i)
!
   end subroutine resetParallelParams
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isEnergyOnMyProc(ie) result(y)
!  ===================================================================
   use MPPModule, only : getMyPE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ie
   integer (kind=IntKind) :: i
!
   logical :: y
!
   character (len=15), parameter :: sname = 'isEnergyOnMyProc'
!
   if (ie < 1 .or. ie > NumEs) then
      call ErrorHandler(sname,'Energy index is out of range',ie)
   else if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   endif
!
   y = .false.
   do i = 1, NumEsPerBox
      if (ie == EsOnMyProc(i)) then
         y = .true.
         exit
      endif
   enddo
!
   end function isEnergyOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumEsOnMyProc() result(ne)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: ne
!
   ne = NumEsPerBox
!
   end function getNumEsOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEnergyIndex(i) result(ie)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: ie
!
   character (len=15), parameter :: sname = 'getEnergyIndex'
!
   if (i < 1 .or. i > NumEsPerBox) then
      call ErrorHandler(sname,'Local energy index is out of range',i)
   endif
!
   ie = EsOnMyProc(i)
!
   end function getEnergyIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcWithEnergyIndex(ie) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ie
   integer (kind=IntKind) :: p
!
   character (len=25), parameter :: sname = 'getProcWithEnergyIndex'
!
   if (ie < 1 .or. ie > NumEs) then
      call ErrorHandler(sname,'Energy index is out of range',ie)
   else if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   endif
!
   p = E2Proc(ie)
!
   end function getProcWithEnergyIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumRedundantEsOnMyProc() result(re)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: re
!
   re = NumRedundantEsPerBox
!
   end function getNumRedundantEsOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isKPointOnMyProc(ik) result(y)
!  ===================================================================
   use MPPModule, only : getMyPE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ik
   integer (kind=IntKind) :: i
!
   logical :: y
!
   character (len=15), parameter :: sname = 'isKPointOnMyProc'
!
   if (ik < 1 .or. ik > NumKs) then
      call ErrorHandler(sname,'K-Point index is out of range',ik)
   else if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   endif
!
   y = .false.
   do i = 1, NumKsPerBox
      if (ik == KsOnMyProc(i)) then
         y = .true.
         exit
      endif
   enddo
!
   end function isKPointOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumKsOnMyProc() result(nk)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: nk
!
   nk = NumKsPerBox
!
   end function getNumKsOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKPointIndex(i) result(ik)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: ik
!
   character (len=15), parameter :: sname = 'getKPointIndex'
!
   if (i < 1 .or. i > NumKsPerBox) then
      call ErrorHandler(sname,'Local k-point index is out of range',i)
   endif
!
   ik = KsOnMyProc(i)
!
   end function getKPointIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcWithKPointIndex(ik) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ik
   integer (kind=IntKind) :: p
!
   character (len=25), parameter :: sname = 'getProcWithKPointIndex'
!
   if (ik < 1 .or. ik > NumKs) then
      call ErrorHandler(sname,'K-Point index is out of range',ik)
   else if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   endif
!
   p = K2Proc(ik)
!
   end function getProcWithKPointIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumRedundantKsOnMyProc() result(rk)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: rk
!
   rk = NumRedundantKsPerBox
!
   end function getNumRedundantKsOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isAtomOnMyProc(ia) result(y)
!  ===================================================================
   use MPPModule, only : getMyPE
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind) :: i
!
   logical :: y
!
   character (len=15), parameter :: sname = 'isAtomOnMyProc'
!
   if (ia < 1 .or. ia > NumAtoms) then
      call ErrorHandler(sname,'Atom index is out of range',ia)
   else if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   endif
!
   y = .false.
   do i = 1, NumAtomsPerProc
      if (ia == AsOnMyProc(i)) then
         y = .true.
         exit
      endif
   enddo
!
   end function isAtomOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createParallelization(isAtomDistributed)
!  ===================================================================
   use MPPModule, only : getNumPEs, getMyPE, syncAllPEs, endMPP
!
   use PrimeFactorsModule, only : getSubFactors
!
   use GroupCommModule, only : createProcGrid, getProcGridID,               &
                               getProcGridSize, getProcGridDimention,       &
                               getMyCoordInProcGrid,                        &
                               createGroupFromGrid, createBoxGroupFromGrid, &
                               getGroupID, getNumPEsInGroup, getMyPEinGroup,&
                               getGroupLabel, getNumGroups,                 &
                               getNumClustersOfGroup, getMyClusterIndex,    &
                               GlobalSumInGroup
!
   implicit none
!
   logical, intent(in), optional :: isAtomDistributed
   logical :: atom_distributed = .true.
!
   integer (kind=IntKind) :: NumPEs, MyPE, NProc_E, MyPEinGroup
   integer (kind=IntKind) :: aGID, eGID, kGID
   integer (kind=IntKind) :: proc_dim(3), box(3), step(3)
   integer (kind=IntKind) :: i, j, d, n, m, grid_id, ie, ik
   integer (kind=IntKind) :: ra, rk, re
   integer (kind=IntKind) :: ra_sav, rk_sav, re_sav
   integer (kind=IntKind) :: na, nk, ne, nsum, np
   integer (kind=IntKind) :: na_sav, nk_sav, ne_sav, nsum_sav, NumAtoms_tmp
   integer (kind=IntKind), pointer :: factors(:,:)
!
   type (MappingStruct), pointer :: p_item
!
   character (len=21), parameter :: sname = 'createParallelization'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   else if (NumAtoms < 1) then
      call ErrorHandler(sname,'Number of atoms < 1',NumAtoms)
   endif
!
   NumPEs = getNumPEs()
   MyPE = getMyPE()
!
   if (present(isAtomDistributed)) then
      atom_distributed = isAtomDistributed
   else
      atom_distributed = .true.
   endif
!
!  ===================================================================
!  determine the number of atoms per processor, the number of k-points
!  per processor, and the number energy points per processor
!
!  try parallelization arrangement
!      na = the estimate of the number of processes in an imagined box,
!           in which the atoms are distributed;
!      nk = an estimate of the number of repeats of the box along k dir
!      ne = an estimate of the number of repeats of the box along e dir
!    nsum = an estimate of the number of atoms, k-points, and energy meshes
!           on each process (core or CPU)
!  ===================================================================
   if (atom_distributed) then
      NumAtoms_tmp = NumAtoms
   else
      NumAtoms_tmp = 1
   endif
   na_sav = 1; ra_sav = NumAtoms_tmp
   nk_sav = 1; rk_sav = NumKs
   ne_sav = 1; re_sav = NumEs
   nsum_sav = NumAtoms_tmp + NumEs + NumKs
   if (NumKs > 0 .and. NumEs > 0) then
      if (NumPEs > NumAtoms_tmp*NumKs*NumEs) then
!        -------------------------------------------------------------
         call WarningHandler(sname,                                      &
                             'Number of processors is more than needed', &
                             NumPEs, NumAtoms_tmp*NumKs*NumEs)
!        -------------------------------------------------------------
      endif
      factors => getSubFactors(NumPEs,3,m) ! break the total number of processes into three factors
      LOOP_i1: do i = 1, m
         if (NumAtomsPerProc > 0 .and. NumAtomsPerProc*factors(1,i) /= NumAtoms) then
            cycle LOOP_i1
         else if (NumAtoms_tmp >= factors(1,i) .and. NumKs >= factors(2,i) .and. &
                  NumEs >= factors(3,i)) then
            ra = mod(NumAtoms_tmp,factors(1,i))
            rk = mod(NumKs,factors(2,i))
            re = mod(NumEs,factors(3,i))
            if (ra == 0) then   ! we only consider the case that atoms are evenly distributed
               na = factors(1,i)
               nk = factors(2,i)
               ne = factors(3,i)
               nsum = NumAtoms_tmp/na + (NumKs-rk)/nk + rk + (NumEs-re)/ne + re
!              =======================================================
!              The following criteria for optimizing the distribution
!              needs to be further tested:
!                 The goal is to achieve the best load balancing by 
!                 requiring the smallest total residual values from
!                 distributing the k- and e- points.
!              =======================================================
               if (re+rk < re_sav+rk_sav) then
!              if (nsum < nsum_sav) then
                  na_sav = na; nk_sav = nk; ne_sav = ne
                  ra_sav = ra; rk_sav = rk; re_sav = re
                  nsum_sav = nsum
               endif
            endif
         endif
      enddo LOOP_i1
   else if (NumKs == 0 .and. NumEs > 0) then
      if (NumPEs > NumAtoms_tmp*NumEs) then
!        -------------------------------------------------------------
         call WarningHandler(sname,                                      &
                             'Number of processors is more than needed', &
                             NumPEs, NumAtoms_tmp*NumEs)
!        -------------------------------------------------------------
      endif
      factors => getSubFactors(NumPEs,2,m) ! break the total number of processes into two factors
      LOOP_i2: do i = 1, m
         if (NumAtomsPerProc > 0 .and. NumAtomsPerProc*factors(1,i) /= NumAtoms) then
            cycle LOOP_i2
         else if (NumAtoms_tmp >= factors(1,i) .and. NumEs >= factors(2,i)) then
            ra = mod(NumAtoms_tmp,factors(1,i))
            re = mod(NumEs,factors(2,i))
            if (ra == 0) then   ! we only consider the case that atoms are evenly distributed
               na = factors(1,i)
               ne = factors(2,i)
               nsum = NumAtoms_tmp/na + (NumEs-re)/ne + re
!              =======================================================
!              The following criteria for optimizing the distribution
!              needs to be further tested
!              =======================================================
               if (re < re_sav) then
!              if (nsum < nsumgetNumPEsInGroup(i)_sav) then
                  na_sav = na; ne_sav = ne
                  ra_sav = ra; re_sav = re
                  nsum_sav = nsum
               endif
            endif
         endif
      enddo LOOP_i2
   else if (NumKs > 0 .and. NumEs == 0) then
      if (NumPEs > NumAtoms_tmp*NumKs) then
!        -------------------------------------------------------------
         call WarningHandler(sname,                                      &
                             'Number of processors is more than needed', &
                             NumPEs, NumAtoms_tmp*NumKs)
!        -------------------------------------------------------------
      endif
      factors => getSubFactors(NumPEs,2,m) ! break the total number of processes into two factors
      LOOP_i3: do i = 1, m
         if (NumAtomsPerProc > 0 .and. NumAtomsPerProc*factors(1,i) /= NumAtoms) then
            cycle LOOP_i3
         else if (NumAtoms_tmp >= factors(1,i) .and. NumKs >= factors(2,i)) then
            ra = mod(NumAtoms_tmp,factors(1,i))
            rk = mod(NumKs,factors(2,i))
            if (ra == 0) then   ! we only consider the case that atoms are evenly distributed
               na = factors(1,i)
               nk = factors(2,i)
               nsum = NumAtoms_tmp/na + (NumKs-rk)/nk + rk
!              =======================================================
!              The following criteria for optimizing the distribution
!              needs to be further tested
!              =======================================================
               if (rk < rk_sav) then
!              if (nsum < nsum_sav) then
                  na_sav = na; nk_sav = nk
                  ra_sav = ra; rk_sav = rk
                  nsum_sav = nsum
               endif
            endif
         endif
      enddo LOOP_i3
   else
      ra_sav = 0
   endif
!
   if (ra_sav /= 0) then
      if (NumAtoms_tmp >= NumPEs .and. mod(NumAtoms_tmp,NumPEs) == 0) then
         na_sav = NumPEs; ra_sav = 0
      else
!        -------------------------------------------------------------
         call StopHandler(sname,'Bad load balance is found')
!        -------------------------------------------------------------
      endif
   endif
!
   if (NumAtomsPerProc == 0) then
      if (atom_distributed) then
         NumAtomsPerProc = (NumAtoms-ra_sav)/na_sav + ra_sav
      else
         NumAtomsPerProc = NumAtoms
      endif
      if (NumAtomsPerProc > MaxAtomsPerProc) then
!        -------------------------------------------------------------
         call WarningHandler(sname,                                   &
                             'Number of atoms/proc > MaxAtomsPerProc',&
                             NumAtomsPerProc, MaxAtomsPerProc)
!        -------------------------------------------------------------
      endif
   endif
!
   NumKsPerBox = (NumKs-rk_sav)/nk_sav + rk_sav
   NumEsPerBox = (NumEs-re_sav)/ne_sav + re_sav
   na = na_sav; nk = nk_sav; ne = ne_sav
   ra = ra_sav; rk = rk_sav; re = re_sav
   if (ne == 1 .and. nk == 1) then
      np = NumPEs
   else
      np = NumAtoms/NumAtomsPerProc  ! = na
   endif
!
   if (print_level >= 0 .and. MyPE == 0) then
      write(6,'(/,80(''=''))')
      write(6,'(/,12x,a)')'***************************************************'
      write(6,'(  12x,a)')'*        Output from createParallelization        *'
      write(6,'(  12x,a)')'*               in ProcMappingModule              *'
      write(6,'(  12x,a)')'*                                                 *'
      write(6,'(  12x,a)')'*  In this parellaization setup, the entire       *'
      write(6,'(  12x,a)')'*  available MPI processes are imagined to form a *'
      write(6,'(  12x,a)')'*  2-D grid, with x-dimention called e-dimension  *'
      write(6,'(  12x,a)')'*  y-dimension called k-dimention. The process    *'
      write(6,'(  12x,a)')'*  grid points are divided into equal sized 2-D   *'
      write(6,'(  12x,a)')'*  boxes. That is, the entire MPI processes are   *'
      write(6,'(  12x,a)')'*  divided into equal sized process pools, each   *'
      write(6,'(  12x,a)')'*  of which forms a MPI group. The atoms are      *'
      write(6,'(  12x,a)')'*  evenly distributed within each box, while the  *'
      write(6,'(  12x,a)')'*  energy points, and the k-points in KKR, are    *'
      write(6,'(  12x,a)')'*  distributed between these imagined boxes.      *'
      write(6,'(12x,a,/)')'***************************************************'
      write(6,'(a,i5)')'The number of processors in each box for atom parallelization:',np
      write(6,'(a,i5)')'The number of repeats of the box along the k-dimension :',nk
      write(6,'(a,i5)')'The number of repeats of the box along the e-dimension :',ne
      write(6,'(a,i5)')'The product of these three numbers   :',np*nk*ne
   endif
!
!  ===================================================================
!  At this point, we have:
!     np = the number of processes in each box
!     nk = the number of repeats of the box along the k-axis
!     ne = the number of repeats of the box along the e-axis
!     NumPEs = np*nk*ne
!  ===================================================================
   n = np*ne*nk
   if (n /= NumPEs) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'np*ne*nk /= NumPEs',n,NumPEs)
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  set up the size of the box: np = box(1)*box(2)*box(3)
!  -------------------------------------------------------------------
   call setupProcDim(np,proc_dim)
!  -------------------------------------------------------------------
   box(1:3) = proc_dim(1:3)
   if (print_level >= 0 .and. MyPE == 0) then
      write(6,'(3(a,i5))')'Dimension of the box (to be repeated to produce processor grid):', &
                          box(1),'  X',box(2),'  X',box(3)
   endif
!
   proc_dim(1) = proc_dim(1)*nk
   proc_dim(2) = proc_dim(2)*ne
!
   if (print_level >= 0 .and. MyPE == 0) then
      write(6,'(3(a,i5))')'Dimension of processors :',proc_dim(1),'  X', &
                          proc_dim(2),'  X',proc_dim(3)
   endif
!
!  -------------------------------------------------------------------
   call createProcGrid(3,'3-D Proc Grid',proc_dim)
!  -------------------------------------------------------------------
   grid_id = getProcGridID('3-D Proc Grid')
!  -------------------------------------------------------------------
   d = getProcGridDimention(grid_id)
!  -------------------------------------------------------------------
   if (MyPE == 0) then
      write(6,'(a,i5)')'Proc Grid Dimension:  ',d
      write(6,'(a,i5)')'Total Number of Procs:',getProcGridSize(grid_id)
      write(6,'(a,$)') 'Proc Grid Size in Each Dimension: ('
      do j = 1, d-1
         write(6,'(i3,a,$)')getProcGridSize(grid_id,j),','
      enddo
      write(6,'(i3,a)')getProcGridSize(grid_id,d),')'
   endif
   if (print_level >= 0 .and. MyPE < min(NumPEs,16)) then
      write(6,'(a,i5,a,$)')'MyPE = ',MyPE,' is at ('
      do j = 1, d-1
         write(6,'(i3,a,$)')getMyCoordInProcGrid(grid_id,j),','
      enddo
      write(6,'(i3,a)')getMyCoordInProcGrid(grid_id,d),')'
   endif
!  -------------------------------------------------------------------
   call createBoxGroupFromGrid(grid_id,box,'Unit Cell')
!  -------------------------------------------------------------------
   if (NumKs > 0) then
      step(1) = box(1); step(2) = 0; step(3) = 0
!     ----------------------------------------------------------------
      call createGroupFromGrid(grid_id,step,'K-Mesh')
!     ----------------------------------------------------------------
      step(1) = proc_dim(1); step(2) = box(2); step(3) = proc_dim(3)
!     ----------------------------------------------------------------
      call createBoxGroupFromGrid(grid_id,step,'A-K Plane')
!     ----------------------------------------------------------------
   endif
   if (NumEs > 0) then
      step(1) = 0; step(2) = box(2); step(3) = 0
!     ----------------------------------------------------------------
      call createGroupFromGrid(grid_id,step,'Energy Mesh')
!     ----------------------------------------------------------------
      step(1) = box(1); step(2) = proc_dim(2); step(3) = proc_dim(3)
!     ----------------------------------------------------------------
      call createBoxGroupFromGrid(grid_id,step,'A-E Plane')
!     ----------------------------------------------------------------
   endif
   if (NumKs > 0 .and. NumEs > 0) then
      step(1) = box(1); step(2) = box(2); step(3) = box(3)
!     ----------------------------------------------------------------
      call createGroupFromGrid(grid_id,step,'E-K Plane')
!     ----------------------------------------------------------------
   endif
!
   if (print_level >= 1 .and. MyPE < min(NumPEs,8)) then
      n = getNumGroups()
      if (MyPE == 0) then
         do i = 1, n
            write(6,'(3a,i5)')'For group: ',getGroupLabel(i),         &
                              ', size = ',getNumPEsInGroup(i)
         enddo
      endif
      do i = 1, n
         write(6,'(2(a,i5),2a)')'For MyPE = ',MyPE,', My rank = ',    &
                                getMyPEinGroup(i),', in group: ',     &
                                getGroupLabel(i)
      enddo
   endif
!
!  ===================================================================
!  Next, one needs to identify which box MyPE belongs to, which box each
!  energy index belongs to, and which box each k-index belongs to. And
!  finanlly, one needs to identify which atoms are mapped onto MyPE.
!  ===================================================================
   NumBoxes = nk*ne
!
   if (print_level >= 1 .and. MyPE < min(NumPEs,8)) then
      d = getProcGridDimention(grid_id)
      do j = 1, d
         write(6,'(3(a,i5))')'MyPE = ',MyPE,', for dim = ',j,         &
                             ', my coordinate in 3-D Proc Grid: ',    &
                             getMyCoordInProcGrid(grid_id,j)
      enddo
   endif
!
   aGID = getGroupID('Unit Cell')
   MyIndexInBox = getMyPEinGroup(aGID) + 1
!
   if (nk > 1) then
      kGID = getGroupID('K-Mesh')
      if (print_level >= 1 .and. MyPE < min(NumPEs,8)) then
         write(6,'(3(a,i5))')'For group K-Mesh, MyPE = ',MyPE,        &
                  ', number of members = ',getNumPEsInGroup(kGID),    &
                  ', my index in group: ',getMyPEinGroup(kGID)
      endif
      MyBoxIndex4K = getMyPEinGroup(kGID) + 1
   else
      MyBoxIndex4K = 1
      kGID = 0
   endif
!
   if (ne > 1) then
      eGID = getGroupID('Energy Mesh')
      NProc_E  = getNumPEsInGroup(eGID)
      if (print_level >= 1 .and. MyPE < min(NumPEs,8)) then
         write(6,'(3(a,i5))')'For group Energy Mesh, MyPE = ',MyPE,   &
                  ', number of members = ',getNumPEsInGroup(eGID),    &
                  ', my index in group: ',getMyPEinGroup(eGID)
      endif
      MyBoxIndex4E = getMyPEinGroup(eGID) + 1
   else
      NProc_E  = 1
      MyBoxIndex4E = 1
      eGID = 0
   endif
!
   if (Created) then
      deallocate(EsOnMyProc, KsOnMyProc, AsOnMyProc, E2Proc, K2Proc, A2Proc)
   endif
   allocate(EsOnMyProc(1:NumEsPerBox), KsOnMyProc(1:NumKsPerBox))
   allocate(AsOnMyProc(1:NumAtomsPerProc))
   allocate(E2Proc(1:NumEs), K2Proc(1:NumKs))
   allocate(A2Proc(1:NumAtoms))
   Created = .true.
!
!  ==================================================================
!  mapping atoms on processors
!  ==================================================================
   do i = 1, NumAtomsPerProc
      AsOnMyProc(i) = (MyIndexInBox-1)*NumAtomsPerProc + i
   enddo
   MyPEinGroup = getMyPEinGroup(aGID)
   A2Proc = 0
   do i = 1, NumAtomsPerProc
      j = AsOnMyProc(i)
      A2Proc(j) = MyPEinGroup
   enddo
!  ------------------------------------------------------------------
   call GlobalSumInGroup(aGID,A2Proc,NumAtoms)
!  ------------------------------------------------------------------
!
!  ==================================================================
!  mapping energy mesh on processors
!    NProc_E = the number of processes in the E-group in which NumEs
!              energy points are distributed
!  ==================================================================
   n = NumEsPerBox - re
   do i = 1, n
      EsOnMyProc(i) = NProc_E*(i-1)+MyBoxIndex4E
   enddo
   do i = 1, re
      EsOnMyProc(n+i) = NumEs - re + i
   enddo
   NumRedundantEsPerBox = re
   if (print_level >= 0 .and. NumEsPerBox > 0 .and. MyPE < min(NumPEs,8)) then
      write(6,'(a,i5,a,2i5)')'MyPE = ',MyPE,', NumEs, NumEsPerBox = ',NumEs, NumEsPerBox
      if (NumEsPerBox <= 20) then
         write(6,'(a,i5,a,20i5)')'MyPE = ',MyPE,', E index = ',        &
                  (EsOnMyProc(i),i=1,NumEsPerBox)
      endif
   endif
   m = NumEs - re
   if (NumEsPerBox > 0) then
      do ie = 1, m
         E2Proc(ie) = mod(ie-1,NProc_E)
      enddo
      if (re > 0) then
         E2Proc(m+1:NumEs) = -1 ! the redundant Es are assigned to all processors in the group
      endif
      if (MyPE == 0) then
         if (NumEs > 20) then
            write(6,'(a,20i5,a)')'Energy point index = ',(ie, ie=1, 20),', ...'
            write(6,'(a,20i5,a)')'Process in E-group = ',(E2Proc(ie), ie=1, 20),', ...'
         else
            write(6,'(a,20i5)')'Energy point index = ',(ie, ie=1, NumEs)
            write(6,'(a,20i5)')'Process in E-group = ',(E2Proc(ie), ie=1, NumEs)
         endif
      endif
   endif
!
!  ==================================================================
!  mapping k mesh on processors
!  ==================================================================
   n = NumKsPerBox - rk
   do i = 1, n
      KsOnMyProc(i) = (MyBoxIndex4K-1)*n + i
   enddo
   do i = 1, rk
      KsOnMyProc(n+i) = NumKs - rk + i
   enddo
   NumRedundantKsPerBox = rk
   if (print_level >= 0 .and. NumKsPerBox > 0 .and. MyPE < min(NumPEs,8)) then
      write(6,'(a,i5,a,2i5)')'MyPE = ',MyPE,', NumKs, NumKsPerBox = ',NumKs, NumKsPerBox
      if (NumKsPerBox > 16) then
         write(6,'(a,i5,a,16i5,a)')'MyPE = ',MyPE,', K index = ',         &
                  (KsOnMyProc(i),i=1,16),', ...,...'
      else
         write(6,'(a,i5,a,16i5)')'MyPE = ',MyPE,', K index = ',         &
                  (KsOnMyProc(i),i=1,NumKsPerBox)
      endif
   endif
   m = NumKs - rk
   if (NumKsPerBox > 0) then
      do ik = 1, m
         K2Proc(ik) = ik/n
      enddo
      if (rk > 0) then
         K2Proc(m+1:NumKs) = -1 ! the redundant Ks are assigned to all processors
      endif
      if (MyPE == 0) then
         if (NumKs > 20) then
            write(6,'(a,20i5,a)')'IBZ K point index  = ',(ik, ik=1, 20),', ...,...'
            write(6,'(a,20i5,a)')'Process in K-group = ',(K2Proc(ik), ik=1, 20),', ...,...'
         else
            write(6,'(a,20i5)')'IBZ K point index  = ',(ik, ik=1, NumKs)
            write(6,'(a,20i5)')'Process in K-group = ',(K2Proc(ik), ik=1, NumKs)
         endif
      endif
   endif
!
!  ===================================================================
!  Add each parameter-to-process mapping rule in the mapping table
!  -------------------------------------------------------------------
   call appendItem(item_name='Atoms in the Unit Cell',                &
                   comm_id = aGID,                                    &
                   num_elements=NumAtoms,                             &
                   num_local_elements=NumAtomsPerProc,                &
                   num_local_redundants=0,                            &
                   element_index=AsOnMyProc,                          &
                   index2proc=A2Proc)
!  -------------------------------------------------------------------
   call appendItem(item_name='E-Mesh along the Energy Contour',       &
                   comm_id = eGID,                                    &
                   num_elements=NumEs,                                &
                   num_local_elements=NumEsPerBox,                    &
                   num_local_redundants=NumRedundantEsPerBox,         &
                   element_index=EsOnMyProc,                          &
                   index2proc=E2Proc)
!  -------------------------------------------------------------------
   call appendItem(item_name='Normal K-Mesh in the IBZ',              &
                   comm_id = kGID,                                    &
                   num_elements=NumKs,                                &
                   num_local_elements=NumKsPerBox,                    &
                   num_local_redundants=NumRedundantKsPerBox,         &
                   element_index=KsOnMyProc,                          &
                   index2proc=K2Proc)
!  -------------------------------------------------------------------
!
   call syncAllPEs()
!
   if (print_level >= 0 .and. MyPE == 0) then
      write(6,'(/,80(''=''),/)')
   endif
!
   end subroutine createParallelization
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupProcDim(num,proc_dim)
!  ===================================================================
   use PrimeFactorsModule, only : getSubFactors
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: num
   integer (kind=IntKind), intent(out) :: proc_dim(1:3)
!
   integer (kind=IntKind), pointer :: factors(:,:)
   integer (kind=IntKind) :: m, i
   integer (kind=IntKind) :: n1, n2, n3, n1s, n2s, n3s
!
   real (kind=RealKind) :: bsize(1:3), dsqr, dsav
!
   bsize(1) = sqrt(Bravais(1,1)**2+Bravais(2,1)**2+Bravais(3,1)**2)
   bsize(2) = sqrt(Bravais(1,2)**2+Bravais(2,2)**2+Bravais(3,2)**2)
   bsize(3) = sqrt(Bravais(1,3)**2+Bravais(2,3)**2+Bravais(3,3)**2)
!
   factors => getSubFactors(num,3,m)
!
   if (m > 0) then
      dsav = 1.0d+20
      do i = 1, m
         n1 = factors(1,i)
         n2 = factors(2,i)
         n3 = factors(3,i)
         dsqr = sqrt( (n1*bsize(2)-n2*bsize(1))**2 + &
                      (n1*bsize(3)-n3*bsize(1))**2 + &
                      (n2*bsize(3)-n3*bsize(2))**2 )
         if (dsqr < dsav) then
            dsav = dsqr
            n1s = n1
            n2s = n2
            n3s = n3
         endif
      enddo
   else
      n1s = num
      n2s = 1
      n3s = 1
   endif
   proc_dim(1) = n1s; proc_dim(2) = n2s; proc_dim(3) = n3s
!
   end subroutine setupProcDim
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine distributeX(XName,CommGroupID,NumXs)
!  ===================================================================
   use MPPModule, only : MyPE, NumPEs, GlobalSUm
!
   use GroupCommModule, only : getNumPEsInGroup, getMyPEinGroup,      &
                               GlobalSumInGroup
!
   implicit none
!
   character (len=*), intent(in) :: XName
!
   integer (kind=IntKind), intent(in) :: CommGroupID
   integer (kind=IntKind), intent(in) :: NumXs
   integer (kind=IntKind) :: nproc, iproc, NumRedundants, nblock, nx_loc
   integer (kind=IntKind) :: i, id, m
   integer (kind=IntKind), allocatable :: xid(:), xid2proc(:)
!
   if (CommGroupID > 0) then
      nproc = getNumPEsInGroup(CommGroupID)
      iproc = getMyPEinGroup(CommGroupID)
   else
      nproc = NumPEs
      iproc = MyPE
   endif
   NumRedundants = mod(NumXs,nproc)
   nblock = (NumXs-NumRedundants)/nproc
   nx_loc = nblock + NumRedundants
!
   allocate(xid(nx_loc), xid2proc(NumXs))
   xid2proc = 0
!
   m = iproc*nblock
   do i = 1, nblock
      id = m + i
      xid(i) = id
      xid2proc(id) = iproc
   enddo
   m = nproc*nblock
   do i = 1, NumRedundants
      id = m + i
      xid(nblock+i) = id
      xid2proc(id) = -1
   enddo
   if (CommGroupID > 0) then
!     ---------------------------------------------------------------
      call GlobalSumInGroup(CommGroupID,xid2proc,NumXs)
!     ---------------------------------------------------------------
   else
!     ---------------------------------------------------------------
      call GlobalSum(xid2proc,NumXs)
!     ---------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   call appendItem(item_name=XName,                                   &
                   comm_id=CommGroupID,                               &
                   num_elements=NumXs,                                &
                   num_local_elements=nx_loc,                         &
                   num_local_redundants=NumRedundants,                &
                   element_index=xid,                                 &
                   index2proc=xid2proc)
!  -------------------------------------------------------------------
!
   deallocate(xid, xid2proc)
!
   end subroutine distributeX
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumXsOnMyProc(X_name) result(n)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: X_name
!
   type (MappingStruct), pointer :: p_item
!
   integer (kind=IntKind) :: n
!
   p_item => getMappingItem(X_name)
!
   n = p_item%NumParamsOnMyProc
!
   end function getNumXsOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getXIndex(X_name,n_loc) result(i)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: X_name
!
   type (MappingStruct), pointer :: p_item
!
   integer (kind=IntKind), intent(in) :: n_loc
   integer (kind=IntKind) :: i
!
   p_item => getMappingItem(X_name)
!
   if (n_loc < 1 .or. n_loc > p_item%NumParamsOnMyProc) then
      call ErrorHandler('getXIndex','The local index is out of range',n_loc)
   endif
!
   i = p_item%ParamOnMyProc(n_loc)
!
   end function getXIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getXCommID(X_name) result(i)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: X_name
!
   type (MappingStruct), pointer :: p_item
!
   integer (kind=IntKind) :: i
!
   p_item => getMappingItem(X_name)
!
   i = p_item%CommGroupID
!
   end function getXCommID
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumRedundantXsOnMyProc(X_name) result(n)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: X_name
!
   type (MappingStruct), pointer :: p_item
!
   integer (kind=IntKind) :: n
!
   p_item => getMappingItem(X_name)
!
   n = p_item%NumRedundantParamsOnMyProc
!
   end function getNumRedundantXsOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isXindexOnMyProc(X_name,i) result(y)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: X_name
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: j
!
   type (MappingStruct), pointer :: p_item
!
   logical :: y
!
   p_item => getMappingItem(X_name)
!
   if (i < 1 .or. i > p_item%NumParams) then
      call ErrorHandler('isXindexOnMyProc','The index is out of range',i)
   endif
!
   y = .false.
   LOOP_j: do j = 1, p_item%NumParamsOnMyProc
      if (i == p_item%ParamOnMyProc(j)) then
         y = .true.
         exit LOOP_j
      endif
   enddo LOOP_j
!
   end function isXindexOnMyProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcWithXIndex(X_name,i) result (n)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: X_name
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: n
!
   type (MappingStruct), pointer :: p_item
!
   p_item => getMappingItem(X_name)
!
   if (i < 1 .or. i > p_item%NumParams) then
      call ErrorHandler('getProcWithXIndex','The index is out of range',i)
   endif
!
!  ===================================================================
!  returns the proc ID in the MPI group, identified by X_name
!  ===================================================================
   n = p_item%ParamToProc(i) ! If index i is mapped on all processors, n = -1
!
   end function getProcWithXIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createMappingItem(p_item)
!  ===================================================================
   implicit none
!
   type (MappingStruct), pointer, intent(out) :: p_item
!
   integer (kind=IntKind) :: i
!
   if (NumMappingItems == 0) then
      allocate(MappingTable)
      p_item => MappingTable
   else
      p_item => MappingTable
      do i = 2, NumMappingItems
         p_item => p_item%next
      enddo
      allocate(p_item%next)
      p_item => p_item%next
   endif
   p_item%ParamKey = ' '
   p_item%CommGroupID = 0
   p_item%NumParams = 0
   p_item%NumParamsOnMyProc = 0
   p_item%NumRedundantParamsOnMyProc = 0
   nullify(p_item%ParamOnMyProc)
   nullify(p_item%ParamToProc)
   nullify(p_item%next)
   NumMappingItems = NumMappingItems + 1
!
   end subroutine createMappingItem
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMappingItem(item_name) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: item_name
!
   type (MappingStruct), pointer :: p
!
   integer (kind=IntKind) :: i
!
   logical :: found = .false.
!
   interface
      function nocaseCompare(s1,s2) result(t)
         implicit none
         logical :: t
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
      end function nocaseCompare
   end interface
!
   found = .false.
   p => MappingTable
   LOOP_i: do i = 1, NumMappingItems
      if (nocaseCompare(p%ParamKey,item_name)) then
         found = .true.
         exit LOOP_i
      endif
      p => p%next
   enddo LOOP_i
!
   if (.not.found) then
      call ErrorHandler('getMappingItem','The item is not in the mapping table',item_name)
   endif
!
   end function getMappingItem
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine appendItem(item_name,comm_id,num_elements,num_local_elements, &
                         num_local_redundants,element_index,index2proc)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: item_name
!
   integer (kind=IntKind), intent(in) :: num_elements
   integer (kind=IntKind), intent(in) :: comm_id
   integer (kind=IntKind), intent(in) :: num_local_elements
   integer (kind=IntKind), intent(in) :: num_local_redundants
   integer (kind=IntKind), intent(in) :: element_index(:)
   integer (kind=IntKind), intent(in) :: index2proc(:)
!
   type (MappingStruct), pointer :: p_item
!
!  -------------------------------------------------------------------
   call createMappingItem(p_item)
!  -------------------------------------------------------------------
!
   p_item%ParamKey = item_name
   p_item%CommGroupID = comm_id
   p_item%NumParams = num_elements
   p_item%NumParamsOnMyProc = num_local_elements
   p_item%NumRedundantParamsOnMyProc = num_local_redundants
!
   allocate(p_item%ParamOnMyProc(num_local_elements))
   p_item%ParamOnMyProc = element_index
!
   allocate(p_item%ParamToProc(num_elements))
   p_item%ParamToProc = index2proc
!
   end subroutine appendItem
!  ===================================================================
end module ProcMappingModule
