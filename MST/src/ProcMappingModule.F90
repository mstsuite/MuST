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
          getProcWithEnergyIndex, &
          isKPointOnMyProc, &
          getNumKsOnMyProc, &
          getKPointIndex,   &
          getNumRedundantKsOnMyProc, &
          getProcWithKPointIndex, &
          isAtomOnMyProc
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
   integer (kind=IntKind), allocatable :: EsOnMyProc(:)
   integer (kind=IntKind), allocatable :: KsOnMyProc(:)
   integer (kind=IntKind), allocatable :: AsOnMyProc(:)
   integer (kind=IntKind), allocatable :: E2Proc(:)
   integer (kind=IntKind), allocatable :: K2Proc(:)
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
            MaxAtomsPerProc = 4
         else
            MaxAtomsPerProc = 8
         endif
      else
         MaxAtomsPerProc = maxp
      endif
   else
      if (isFullPotential) then
         MaxAtomsPerProc = 4
      else
         MaxAtomsPerProc = 8
      endif
   endif
!
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
            MaxAtomsPerProc = 4
         else
            MaxAtomsPerProc = 8
         endif
      else
         MaxAtomsPerProc = maxp
      endif
   else
      if (isFullPotential) then
         MaxAtomsPerProc = 4
      else
         MaxAtomsPerProc = 8
      endif
   endif
!
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
   if (Created) then
      deallocate(EsOnMyProc, KsOnMyProc, AsOnMyProc)
      deallocate(E2Proc, K2Proc)
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
   isFullPotential = .false.
   Initialized = .false.
!
   end subroutine endProcMapping
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
                               getNumClustersOfGroup, getMyClusterIndex
!
   implicit none
!
   logical, intent(in), optional :: isAtomDistributed
   logical :: atom_distributed = .true.
!
   integer (kind=IntKind) :: NumPEs, MyPE, NProc_E
   integer (kind=IntKind) :: proc_dim(3), box(3), step(3)
   integer (kind=IntKind) :: i, j, d, n, m, grid_id, ie, ik
   integer (kind=IntKind) :: ra, rk, re
   integer (kind=IntKind) :: ra_sav, rk_sav, re_sav
   integer (kind=IntKind) :: na, nk, ne, nsum, np
   integer (kind=IntKind) :: na_sav, nk_sav, ne_sav, nsum_sav, NumAtoms_tmp
   integer (kind=IntKind), pointer :: factors(:,:)
!
   character (len=21), parameter :: sname = 'createParallelization'
!
   if (.not.Initialized) then
      call ErrorHandler(sname,'Module not initialized')
   else if (NumAtoms < 1) then
      call ErrorHandler(sname,'Number of atoms < 1',NumAtoms)
   endif
!
   if (present(isAtomDistributed)) then
      atom_distributed = isAtomDistributed
   else
      atom_distributed = .true.
   endif
!
   NumPEs = getNumPEs()
   MyPE = getMyPE()
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
      do i = 1, m
         if (NumAtoms_tmp >= factors(1,i) .and. NumKs >= factors(2,i) .and. &
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
!              needs to be further tested
!              =======================================================
               if (re+rk < re_sav+rk_sav) then
!              if (nsum < nsum_sav) then
                  na_sav = na; nk_sav = nk; ne_sav = ne
                  ra_sav = ra; rk_sav = rk; re_sav = re
                  nsum_sav = nsum
               endif
            endif
         endif
      enddo
   else if (NumKs == 0 .and. NumEs > 0) then
      if (NumPEs > NumAtoms_tmp*NumEs) then
!        -------------------------------------------------------------
         call WarningHandler(sname,                                      &
                             'Number of processors is more than needed', &
                             NumPEs, NumAtoms_tmp*NumEs)
!        -------------------------------------------------------------
      endif
      factors => getSubFactors(NumPEs,2,m) ! break the total number of processes into two factors
      do i = 1, m
         if (NumAtoms_tmp >= factors(1,i) .and. NumEs >= factors(2,i)) then
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
      enddo
   else if (NumKs > 0 .and. NumEs == 0) then
      if (NumPEs > NumAtoms_tmp*NumKs) then
!        -------------------------------------------------------------
         call WarningHandler(sname,                                      &
                             'Number of processors is more than needed', &
                             NumPEs, NumAtoms_tmp*NumKs)
!        -------------------------------------------------------------
      endif
      factors => getSubFactors(NumPEs,2,m) ! break the total number of processes into two factors
      do i = 1, m
         if (NumAtoms_tmp >= factors(1,i) .and. NumKs >= factors(2,i)) then
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
      enddo
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
   if (atom_distributed) then
      NumAtomsPerProc = (NumAtoms-ra_sav)/na_sav + ra_sav
   else
      NumAtomsPerProc = NumAtoms
   endif
!
   if (NumAtomsPerProc > MaxAtomsPerProc) then
!     ----------------------------------------------------------------
      call WarningHandler(sname,                                   &
                          'Number of atoms/proc > MaxAtomsPerProc',&
                          NumAtomsPerProc, MaxAtomsPerProc)
!     ----------------------------------------------------------------
   endif
!
   NumKsPerBox = (NumKs-rk_sav)/nk_sav + rk_sav
   NumEsPerBox = (NumEs-re_sav)/ne_sav + re_sav
   na = na_sav; nk = nk_sav; ne = ne_sav
   ra = ra_sav; rk = rk_sav; re = re_sav
   np = NumAtoms/NumAtomsPerProc  ! = na
!
   if (print_level >= 0 .and. MyPE == 0) then
      write(6,'(/,80(''=''))')
      write(6,'(/,12x,a)')'*********************************************'
      write(6,'(  12x,a)')'*     Output from createParallelization     *'
      write(6,'(  12x,a)')'*            in ProcMappingModule           *'
      write(6,'(12x,a,/)')'*********************************************'
      write(6,'(a,i5)')'The number of processors in each box :',np
      write(6,'(a,i5)')'The number of repeats along k-dimen. :',nk
      write(6,'(a,i5)')'The number of repeats along e-dimen. :',ne
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
      write(6,'(3(a,i5))')'Dimension of the box :',box(1),' x',box(2),' x',box(3)
   endif
!
   proc_dim(1) = proc_dim(1)*nk
   proc_dim(2) = proc_dim(2)*ne
!
!  -------------------------------------------------------------------
   call createProcGrid(3,'3-D Proc Grid',proc_dim)
!  -------------------------------------------------------------------
   grid_id = getProcGridID('3-D Proc Grid')
!  -------------------------------------------------------------------
   if (print_level >= 1 .and. MyPE < min(NumPEs,8)) then
!     ----------------------------------------------------------------
      d = getProcGridDimention(grid_id)
!     ----------------------------------------------------------------
      if (MyPE == 0) then
         write(6,'(a,i5)')'Grid Dimension: ',d
         write(6,'(a,i5)')'Grid Size: ',getProcGridSize(grid_id)
         do j = 1, d
            write(6,'(a,i5,a,i5)')'For dim = ',j,                     &
                                  ', size = ',getProcGridSize(grid_id,j)
         enddo
      endif
      do j = 1, d
         write(6,'(3(a,i5))')'MyPE = ',MyPE,', for dim = ',j,         &
                             ',my coordinate: ',getMyCoordInProcGrid(grid_id,j)
      enddo
   endif
!  -------------------------------------------------------------------
   call createBoxGroupFromGrid(grid_id,box,'Unit Cell')
!  -------------------------------------------------------------------
   if (NumKs > 0) then
      step(1) = box(1); step(2) = 0; step(3) = 0
!     ----------------------------------------------------------------
      call createGroupFromGrid(grid_id,step,'K-Mesh')
!     ----------------------------------------------------------------
   endif
   if (NumEs > 0) then
      step(1) = 0; step(2) = box(2); step(3) = 0
!     ----------------------------------------------------------------
      call createGroupFromGrid(grid_id,step,'Energy Mesh')
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
   i = getGroupID('Unit Cell')
   MyIndexInBox = getMyPEinGroup(i) + 1
!
   if (nk > 1) then
      i = getGroupID('K-Mesh')
      if (print_level >= 1 .and. MyPE < min(NumPEs,8)) then
         write(6,'(3(a,i5))')'For group K-Mesh, MyPE = ',MyPE,        &
                  ', number of members = ',getNumPEsInGroup(i),       &
                  ', my index in group: ',getMyPEinGroup(i)
      endif
      MyBoxIndex4K = getMyPEinGroup(i) + 1
   else
      MyBoxIndex4K = 1
   endif
!
   if (ne > 1) then
      i = getGroupID('Energy Mesh')
      NProc_E  = getNumPEsInGroup(i)
      if (print_level >= 1 .and. MyPE < min(NumPEs,8)) then
         write(6,'(3(a,i5))')'For group Energy Mesh, MyPE = ',MyPE,   &
                  ', number of members = ',getNumPEsInGroup(i),       &
                  ', my index in group: ',getMyPEinGroup(i)
      endif
      MyBoxIndex4E = getMyPEinGroup(i) + 1
   else
      NProc_E  = 1
      MyBoxIndex4E = 1
   endif
!
   if (Created) then
      deallocate(EsOnMyProc, KsOnMyProc, AsOnMyProc, E2Proc, K2Proc)
   endif
   allocate(EsOnMyProc(1:NumEsPerBox), KsOnMyProc(1:NumKsPerBox))
   allocate(AsOnMyProc(1:NumAtomsPerProc))
   allocate(E2Proc(1:NumEs), K2Proc(1:NumKs))
   Created = .true.
!
!  ==================================================================
!  mapping atoms on processors
!  ==================================================================
   do i = 1, NumAtomsPerProc
      AsOnMyProc(i) = (MyIndexInBox-1)*NumAtomsPerProc + i
   enddo
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
   proc_dim(1) = n1s; proc_dim(2) = n2s; proc_dim(3) = n3s
!
   end subroutine setupProcDim
!  ===================================================================
end module ProcMappingModule
