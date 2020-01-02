!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine constructDataOnGrid( grid_name, value_name, value_type, &
                                   getDataAtPoint, DenOnGrid, lmax, spin )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
   use MathParamModule, only : ZERO
!
   use GroupCommModule, only : getGroupID, GlobalSumInGroup, GlobalMaxInGroup
   use GroupCommModule, only : getGroupCommunicator
   use GroupCommModule, only : getNumPEsInGroup, getMyPEinGroup
!
   use MPPModule, only : setCommunicator, resetCommunicator
   use MPPModule, only : nbsendMessage, nbrecvMessage, waitMessage, AnyPE
!
   use PublicTypeDefinitionsModule, only : UniformGridStruct
   use PublicTypeDefinitionsModule, only : AtomOnUniformGridStruct
!
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
   use AtomModule, only : getLocalAtomPosition
!
   use ProcMappingModule, only : isAtomOnMyProc
!
   use Atom2ProcModule, only : getLocalIndex
!
   use Uniform3DGridModule, only : getUniform3DGrid, getGridIndex, getGridPosition
   use Uniform3DGridModule, only : isOnAtomicCellBoundary, isLocalGrid
   use Uniform3DGridModule, only : getSourceProc, getTargetProc
!
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: lmax, spin
!
   character(len=*), intent(in) :: grid_name, value_name, value_type
!
   real (kind=RealKind), intent(out) :: DenOnGrid(:)
!
   type (UniformGridStruct), pointer :: gp
!
   logical :: isDensityType
!
   integer (kind=IntKind) :: gCounter, ig, ig_box, jmax=0
   integer (kind=IntKind) :: gid, id, ia, i, GroupID, NumPEsInAGroup, MyPEinAGroup
   integer (kind=IntKind) :: atomsOnPoint, n_mult, n, maxn(3), ng, local_id
   integer (kind=IntKind) :: nSources, nTargets, comm
   integer (kind=IntKind), pointer :: p_source(:), p_target(:)
   integer (kind=IntKind), allocatable :: recv_msgid1(:), send_msgid1(:)
   integer (kind=IntKind), allocatable :: recv_msgid2(:), send_msgid2(:)
   integer (kind=IntKind), allocatable :: grid_local(:), grid_remote(:,:)
!
   real (kind=RealKind) :: r(3)
   real (kind=RealKind), allocatable :: den_local(:), den_remote(:,:)
!
   interface
      function getDataAtPoint( dname, id, ia, r, jmax_in, n, grad ) result(v)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         character (len=*), intent(in) :: dname
         integer (kind=IntKind), intent(in) :: id, ia
         real (kind=RealKind), intent(in) :: r(3)
         integer (kind=IntKind), intent(in), optional :: jmax_in, n
         real (kind=RealKind), intent(out), optional :: grad(3)
         real (kind=RealKind) :: v
      end function getDataAtPoint
   end interface
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   gp => getUniform3DGrid(grid_name)
!
   if (nocaseCompare(value_name,'Potential') .or.                     &
       nocaseCompare(value_name,'WaveFunct')) then
      isDensityType = .false.
   else if (nocaseCompare(value_name,'Charge') .or.                   &
            nocaseCompare(value_name,'Moment')) then
      isDensityType = .true.
   else
      call ErrorHandler("constructDataOnGrid","Invalid value name",value_name)
   endif
!
   if ( present(lmax) ) then
      if (lmax < 0) then
         call ErrorHandler("constructDataOnGrid","Invalid lmax",lmax)
      endif
      jmax = (lmax+1)*(lmax+2)/2
   endif
!
   if ( present(spin) ) then
      if (spin < 1 .or. spin > 2) then
         call ErrorHandler("constructDataOnGrid","Invalid spin index",spin)
      else if ( isDensityType ) then
         call ErrorHandler("constructDataOnGrid","No spin index is needed",spin)
      endif
   endif
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInAGroup = getNumPEsInGroup(GroupID)
   MyPEinAGroup = getMyPEinGroup(GroupID)
!
   maxn = 0
   do id = 1, gp%AtomOnGrid%NumLocalAtoms
      p_source => getSourceProc(gp,id,nSources)
      p_target => getTargetProc(gp,id,nTargets)
      maxn(1) = max(maxn(1),gp%AtomOnGrid%NumGridPointsInCell(id))
      maxn(2) = max(maxn(2),nSources)
      maxn(3) = max(maxn(3),nTargets)
   enddo
!  -------------------------------------------------------------------
   call GlobalMaxInGroup(GroupID,maxn,2)
!  -------------------------------------------------------------------
   allocate(grid_local(maxn(1)+2),den_local(maxn(1)))
   allocate(grid_remote(maxn(1)+2,maxn(2)), den_remote(maxn(1),maxn(2)))
   allocate(send_msgid1(maxn(3)),recv_msgid1(maxn(2)))
   allocate(send_msgid2(maxn(3)),recv_msgid2(maxn(2)))
!
   comm = getGroupCommunicator(GroupID)
!  -------------------------------------------------------------------
   call setCommunicator(comm,MyPEinAGroup,NumPEsInAGroup,sync=.true.)
!  -------------------------------------------------------------------
!
   DenOnGrid = ZERO
!  ===================================================================
!  loop over grid points in each atom box on my process
!  ===================================================================
   do id = 1, gp%AtomOnGrid%NumLocalAtoms
      p_source => getSourceProc(gp,id,nSources)
      p_target => getTargetProc(gp,id,nTargets)
!
!     ================================================================
!     Receive the atom box data for the grids allocated on my process. 
!     ================================================================
      den_remote = ZERO
      do i = 1, nSources
         recv_msgid1(i)=nbrecvMessage(grid_remote(:,i),maxn(1)+2,15000+id,p_source(i))
         recv_msgid2(i)=nbrecvMessage(den_remote(:,i),maxn(1),25000+id,p_source(i))
      enddo
!
      ng = gp%AtomOnGrid%NumGridPointsInCell(id)
      grid_local(1) = id
      grid_local(2) = ng
      LOOP_ig: do ig = 1, ng
!        =============================================================
!        ig_box is the grid index in AtomBox
!        =============================================================
         ig_box = gp%AtomOnGrid%InCellGridPointABIndex(ig,id)
         den_local(ig) = ZERO
         gCounter = getGridIndex(gp,id,ig_box)
         grid_local(ig+2) = gCounter
!
!        =============================================================
!        n_mult: the number of atoms that are associated with the current 
!        grid point.
!        =============================================================
         n_mult = 1
         if (isOnAtomicCellBoundary(gp,gCounter,n)) then
            n_mult = n
         endif
         r = getGridPosition(gp,id,ig_box) - gp%AtomOnGrid%AtomPosition(1:3,id)
!
         if (isDensityType) then
            if (present(lmax)) then
               do ia = 1, getLocalNumSpecies(id)
                  den_local(ig) = den_local(ig) +                          &
                                  getDataAtPoint(value_type, id, ia, r,    &
                                                 jmax_in=jmax, n=n_mult)*getLocalSpeciesContent(id,ia)
               enddo
            else
               do ia = 1, getLocalNumSpecies(id)
                  den_local(ig) = den_local(ig) +                          &
                                  getDataAtPoint(value_type, id, ia, r,    &
                                                 n=n_mult)*getLocalSpeciesContent(id,ia)
               enddo
            endif
         else
            if (present(lmax) .and. present(spin)) then
               do ia = 1, getLocalNumSpecies(id)
                  den_local(ig) = den_local(ig) +                             &
                                  getDataAtPoint(value_type, id, ia, r,       &
                                                 jmax_in=jmax, n=spin)*getLocalSpeciesContent(id,ia)
               enddo
            else if (present(lmax)) then
               do ia = 1, getLocalNumSpecies(id)
                  den_local(ig) = den_local(ig) +                             &
                                  getDataAtPoint(value_type, id, ia, r,       &
                                                 jmax_in=jmax)*getLocalSpeciesContent(id,ia)
               enddo
            else if (present(spin)) then
               do ia = 1, getLocalNumSpecies(id)
                  den_local(ig) = den_local(ig) +                             &
                                  getDataAtPoint(value_type, id, ia, r,       &
                                                 n=spin)*getLocalSpeciesContent(id,ia)
               enddo
            else
               do ia = 1, getLocalNumSpecies(id)
                  den_local(ig) = den_local(ig) +                             &
                                  getDataAtPoint(value_type, id, ia, r)*getLocalSpeciesContent(id,ia)
               enddo
            endif
         endif
!        =============================================================
!        Note: since DenOnGrid will be summed over across all atoms, den_local
!              needs to be divided by n_mult so that the density/potential 
!              value at this grid point is averaged over the assoicated atoms.
!
!              local_id = the index of the uniform grid points allocated on
!                         my process for parallel processing, e.g. FFT
!        =============================================================
         den_local(ig) = den_local(ig)/real(n_mult,kind=RealKind)
         if (isLocalGrid(gp,gCounter,local_id)) then
            DenOnGrid(local_id) = DenOnGrid(local_id) + den_local(ig)
         endif
      enddo LOOP_ig
!
!     ================================================================
!     Send the data for the atom box to the MPI processes that the atom
!     box is required for the parallel processing (e.g., FFT).
!     ================================================================
      do i = 1, nTargets
         send_msgid1(i)=nbsendMessage(grid_local,maxn(1)+2,15000+id,p_target(i))
         send_msgid2(i)=nbsendMessage(den_local,maxn(1),25000+id,p_target(i))
      enddo
!
      do i = 1, nSources
!        -------------------------------------------------------------
         call waitMessage(recv_msgid1(i))
         call waitMessage(recv_msgid2(i))
!        -------------------------------------------------------------
         if (grid_remote(1,i) /= id) then
!           ----------------------------------------------------------
            call ErrorHandler("constructDataOnGrid","Invalid atom index",grid_remote(1,i),id)
!           ----------------------------------------------------------
         endif
         do ig = 1, grid_remote(2,i)
            gCounter = grid_remote(ig+2,i)
            if (isLocalGrid(gp,gCounter,local_id)) then
               DenOnGrid(local_id) = DenOnGrid(local_id) + den_remote(ig,i)
            endif
         enddo
      enddo
!
      do i = 1, nTargets
!        -------------------------------------------------------------
         call waitMessage(send_msgid1(i))
         call waitMessage(send_msgid2(i))
!        -------------------------------------------------------------
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call resetCommunicator()
!  -------------------------------------------------------------------
!
   deallocate(den_local, den_remote)
   deallocate(grid_local, grid_remote)
   deallocate(send_msgid1,recv_msgid1)
   deallocate(send_msgid2,recv_msgid2)
!
   end subroutine constructDataOnGrid
