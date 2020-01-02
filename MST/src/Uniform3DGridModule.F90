module Uniform3DGridModule
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO, ONE
   use PublicTypeDefinitionsModule, only : UniformGridStruct
   use PublicTypeDefinitionsModule, only : AtomOnUniformGridStruct
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
!
public ::                      &
   initUniform3DGrid,          &
   endUniform3DGrid,           &
   isUniform3DGridInitialized, &
   createUniform3DGrid,        &
   createProcessorMesh,        &
   distributeUniformGrid,      &
   insertAtomsInGrid,          &
   printUniform3DGrid,         &
   getUniform3DGrid,           &
   getNumGridPoints,           &
   getGridStepSize,            &
   getGridPosition,            &
   getGridIndex,               &
   isOnAtomicCellBoundary,     &
   isLocalGrid,                &
   getSourceProc,              &
   getTargetProc
!
   interface distributeUniformGrid
      module procedure distributeUniformGrid_s, distributeUniformGrid_p
   end interface distributeUniformGrid
!
   interface isOnAtomicCellBoundary
      module procedure isOnAtomicCellBoundary_s, isOnAtomicCellBoundary_g
   end interface isOnAtomicCellBoundary
!
   interface getGridIndex
      module procedure getGridIndex_box, getGridIndex_local
   end interface getGridIndex
!
   interface getGridPosition
      module procedure getGridPosition_box, getGridPosition_global
   end interface getGridPosition
!
private
   type (UniformGridStruct), pointer :: FFTGrid
   type (UniformGridStruct), pointer :: VisualGrid
!
   logical :: Initialized = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: NumUniformGrids = 0
   integer (kind=IntKind) :: print_level
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initUniform3DGrid(istop,iprint)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: iprint
!
   if (Initialized) then
      call WarningHandler('initUniform3DGrid','module has already been initialized')
   endif
!
   stop_routine = istop
   print_level = iprint
!
   NumUniformGrids = 0
!
   Initialized = .true.
!
   end subroutine initUniform3DGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endUniform3DGrid()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
      call WarningHandler('endUniform3DGrid','module is not initialized')
      return
   endif
!
   if (associated(FFTGrid)) then
      call deleteGridSpace(FFTGrid)
      deallocate(FFTGrid)
   endif
!
   if (associated(VisualGrid)) then
      call deleteGridSpace(VisualGrid)
      deallocate(VisualGrid)
   endif
!
   NumUniformGrids = 0
!
   Initialized = .false.
!
   end subroutine endUniform3DGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isUniform3DGridInitialized() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if (Initialized) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isUniform3DGridInitialized
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteGridSpace(pUG)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   type (UniformGridStruct), intent(inout) :: pUG
!
   if (associated(pUG%AtomOnGrid)) then
      deallocate(pUG%AtomOnGrid%AtomBox)
      deallocate(pUG%AtomOnGrid%AtomPosition)
      deallocate(pUG%AtomOnGrid%NumGridPointsInAtomBox)
      deallocate(pUG%AtomOnGrid%NumGridPointsInCell)
      deallocate(pUG%AtomOnGrid%InCellGridPointABIndex)
      if (pUG%AtomOnGrid%NumGridPointsOnCellBound > 0) then
         deallocate(pUG%AtomOnGrid%CBGridPointIndex)
      endif
      if (allocated(pUG%AtomOnGrid%TargetProc)) then
         deallocate(pUG%AtomOnGrid%TargetProc)
      endif
      if (allocated(pUG%AtomOnGrid%SourceProc)) then
         deallocate(pUG%AtomOnGrid%SourceProc)
      endif
      deallocate(pUG%AtomOnGrid%NumTargetProcs, pUG%AtomOnGrid%NumSourceProcs)
      deallocate(pUG%AtomOnGrid)
   endif
!
   end subroutine deleteGridSpace
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createUniform3DGrid(gname,na,nb,nc,unit_cell,cell_origin)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind), intent(in) :: na, nb, nc
   integer (kind=IntKind) :: ic
!
   real (kind=RealKind), intent(in) :: unit_cell(3,3)
   real (kind=RealKind), optional :: cell_origin(3)
   real (kind=RealKind) :: v0(3)
!
   type (UniformGridStruct), pointer :: pUG
!
   if (.not.Initialized) then
      call ErrorHandler('createUniform3DGrid','module not initialized')
   else if (na < 1) then
      call ErrorHandler('createUniform3DGrid','na < 1',na)
   else if (nb < 1) then
      call ErrorHandler('createUniform3DGrid','nb < 1',nb)
   else if (nc < 1) then
      call ErrorHandler('createUniform3DGrid','nc < 1',nc)
   endif
!
   if (present(cell_origin)) then
      v0 = cell_origin
   else
      v0 = ZERO
   endif
!
!  ===================================================================
!  For simplicity, we only consider two uniform grids. For general purposes,
!  this should be implemented in terms of linked list.
!  ===================================================================
   if (nocaseCompare(gname,'FFT')) then
      allocate( FFTGrid )
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      allocate( VisualGrid )
      pUG => VisualGrid
   else
      call ErrorHandler('createUniform3DGrid','Unknown grid name',gname)
   endif
   NumUniformGrids = NumUniformGrids + 1
!  ===================================================================
   nullify(pUG%AtomOnGrid)
!
   pUG%nga = na
   pUG%ngb = nb
   pUG%ngc = nc
   pUG%ng = na*nb*nc
   pUG%nproc_a = 1
   pUG%nproc_b = 1
   pUG%nproc_c = 1
   pUG%vec_origin = v0
   pUG%cell(1:3,1) = unit_cell(1:3,1)-v0(1:3)
   pUG%cell(1:3,2) = unit_cell(1:3,2)-v0(1:3)
   pUG%cell(1:3,3) = unit_cell(1:3,3)-v0(1:3)
!
   pUG%grid_step_a(1:3) = pUG%cell(1:3,1)/real(na,kind=RealKind)
   pUG%grid_step_b(1:3) = pUG%cell(1:3,2)/real(nb,kind=RealKind)
   pUG%grid_step_c(1:3) = pUG%cell(1:3,3)/real(nc,kind=RealKind)
!
   end subroutine createUniform3DGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createProcessorMesh(gname,pmesh)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind), intent(in) :: pmesh(3)
!
   type (UniformGridStruct), pointer :: pUG
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('createProcessorMesh','Unknown grid name',gname)
   endif
!
   pUG%nproc_a = pmesh(1)
   pUG%nproc_b = pmesh(2)
   pUG%nproc_c = pmesh(3)
!
   end subroutine createProcessorMesh
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine distributeUniformGrid_s(gname) ! Serial implementation
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind) :: ic
!
   type (UniformGridStruct), pointer :: pUG
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('distributeUniformGrid','Unknown grid name',gname)
   endif
!
   pUG%gstart = 1
   pUG%gend(1) = pUG%nga
   pUG%gend(2) = pUG%ngb
   pUG%gend(3) = pUG%ngc
   pUG%NumLocalGridPoints = pUG%ng
!
   end subroutine distributeUniformGrid_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine distributeUniformGrid_p(gname,grid_start,grid_end)  ! parallel implementation
!  ===================================================================
!  Note: Uniform3DGridModule should NOT know which package/scheme is
!        is used to perform  parallel FFT. The information of grid
!        distribution on processors should be obtained by the caller
!        and passed in as grid_start and grid_end, via subroutine arguments.
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind), intent(in) :: grid_start(3), grid_end(3)
   integer (kind=IntKind) :: ic, ng_local
   integer (kind=IntKind) :: num_local_grid_segments, num_global_grid_segments
!
   type (UniformGridStruct), pointer :: pUG
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('distributeUniformGrid','Unknown grid name',gname)
   endif
!
   if (grid_end(1)-grid_start(1) < 0) then
      call ErrorHandler('distributeUniformGrid','grid_start(1) > grid_end(1)',grid_start(1),grid_end(1))
   else if (grid_end(2)-grid_start(2) < 0) then
      call ErrorHandler('distributeUniformGrid','grid_start(2) > grid_end(2)',grid_start(2),grid_end(2))
   else if (grid_end(3)-grid_start(3) < 0) then
      call ErrorHandler('distributeUniformGrid','grid_start(3) > grid_end(3)',grid_start(3),grid_end(3))
   endif
!
   ng_local = (grid_end(1)-grid_start(1)+1)*(grid_end(2)-grid_start(2)+1)*(grid_end(3)-grid_start(3)+1)
   pUG%gstart = grid_start
   pUG%gend = grid_end
   pUG%NumLocalGridPoints = ng_local
!
!  For real space uniform grid with 1-D indexing 1, 2, ..., ng, determine the number of segment
   num_local_grid_segments = (grid_end(2)-grid_start(2)+1)*(grid_end(3)-grid_start(3)+1)
!
   end subroutine distributeUniformGrid_p
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insertAtomsInGrid(gname, NumAtoms, AtomPosition,        &
                                getPointLocationFlag, radius)
!  ===================================================================
!  This routine identifies the uniform grid points local to the input atoms,
!  which are local to this processor. However, the identified grid points 
!  may not be local (when used with parllel FFT or parallel visualization)
!  to this processor.
!
!  Note: if the input atoms (NumAtyoms, AtomPosition) are local to the
!        this processor, the local uniform grid points are possibly not
!        associated with any local atoms, since atoms and uniform grids 
!        have independent parallel distribution.
!  ===================================================================
   use MathParamModule, only : TEN2m12, TEN2m6, TWO, PI
   use VectorModule, only : getVecLength, getDotProduct, getCrossProduct
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalCollectInGroup, GlobalMaxInGroup
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind), intent(in) :: NumAtoms
!
   real (kind=RealKind), intent(in) :: radius(:)
!
   integer (kind=IntKind) :: gCounter, ic, jc, kc, GroupID
   integer (kind=IntKind) :: ia, lp, mp, np, nshift, ig, max_incell_pts
   integer (kind=IntKind) :: icp, jcp, kcp, i, j, k, n, m, ip, nc
   integer (kind=IntKind), allocatable :: gtmp(:), ntmp(:)
   integer (kind=IntKind), allocatable :: data_collect(:,:,:)
   integer (kind=IntKind), allocatable :: size_collect(:)
   integer (kind=IntKind) :: aGID, NumPEsInAGroup, MyPEinAGroup
!
   real (kind=RealKind) :: AtomPosition(3,NumAtoms)
   real (kind=RealKind) :: r(3), rg(3), rv(3)
   real (kind=RealKind) :: a0(3), b0(3), c0(3), fac, max_rc, vol, frac
   real (kind=RealKind) :: a, b, c, step_a_len, step_b_len, step_c_len
   real (kind=RealKind) :: bxc, cxa, axb
!
   logical :: redundant
!
   type (AtomOnUniformGridStruct), pointer :: pAOG
   type (UniformGridStruct), pointer :: pUG
!
   type PointOnCellStruct
      integer (kind=IntKind) :: grid_index
      type (PointOnCellStruct), pointer :: next
   end type PointOnCellStruct
!
   type (PointOnCellStruct), pointer :: p_head, linked_list
!
   interface
      function getPointLocationFlag(i,x,y,z) result(k)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: i
         integer (kind=IntKind) :: k
         real (kind=RealKind), intent(in) :: x, y, z
      end function getPointLocationFlag
   end interface
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('getUniform3DGrid','Unknown grid name',gname)
   endif
!
   if (associated(pUG%AtomOnGrid)) then
      call ErrorHandler('insertAtomsInGrid','AtomOnGrid has already been allocated')
   endif
!
   nullify(p_head, linked_list)
!
!  ===================================================================
!  Set up AtomBox, which identifies a box of uniform grid points with
!  an atom at the center of the box. The coordinates of the box is 
!  relative to the origin of the unit cell.
!  ===================================================================
   allocate(pUG%AtomOnGrid)
   pAOG => pUG%AtomOnGrid
   allocate(pAOG%AtomBox(2,3,NumAtoms), pAOG%AtomPosition(3,NumAtoms))
   allocate(pAOG%NumGridPointsInAtomBox(NumAtoms))
   allocate(pAOG%NumGridPointsInCell(NumAtoms))
   pAOG%NumLocalAtoms = NumAtoms
!
!  ===================================================================
!  Calculate the unit vector of the cell axes and store the vectors in
!  a0, b0, and c0.
!  ===================================================================
   a0 = pUG%cell(1:3,1)/getVecLength(3,pUG%cell(:,1))
!  write(6,'(a,d15.8)')'cell a legnth = ',getVecLength(3,pUG%cell(:,1))
!  write(6,'(a,d15.8)')'a0 legnth = ',getVecLength(3,a0)
   b0 = pUG%cell(1:3,2)/getVecLength(3,pUG%cell(:,2))
!  write(6,'(a,d15.8)')'cell b legnth = ',getVecLength(3,pUG%cell(:,2))
!  write(6,'(a,d15.8)')'b0 legnth = ',getVecLength(3,b0)
   c0 = pUG%cell(1:3,3)/getVecLength(3,pUG%cell(:,3))
!  write(6,'(a,d15.8)')'cell c legnth = ',getVecLength(3,pUG%cell(:,3))
!  write(6,'(a,d15.8)')'c0 legnth = ',getVecLength(3,c0)
!
!  ===================================================================
!  Calculate fac = a0*(b0xc0).
!  ===================================================================
   rv = getCrossProduct(a0,b0)
   fac = getDotProduct(3,rv,c0,absolute=.true.)
   if (abs(fac) < TEN2m6) then
      call ErrorHandler('insertAtomsInGrid','fac = 0',fac)
   endif
!
   step_a_len = getVecLength(3,pUG%grid_step_a)
!  write(6,'(a,d15.8)')'step_a_len*nga = ',step_a_len*pUG%nga
   step_b_len = getVecLength(3,pUG%grid_step_b)
!  write(6,'(a,d15.8)')'step_b_len*ngb = ',step_b_len*pUG%ngb
   step_c_len = getVecLength(3,pUG%grid_step_c)
!  write(6,'(a,d15.8)')'step_c_len*ngc = ',step_c_len*pUG%ngc
!
   max_incell_pts = 0
   do ia = 1, NumAtoms
      pAOG%AtomPosition(1:3,ia) = AtomPosition(1:3,ia)
!     ================================================================
!     Note: Given any point r=(x,y,z), its projection onto the cell axis
!     vectors, a0, b0, c0, are given by:
!         a = r*(b0xc0)/[a0*(b0xc0)]
!         b = r*(c0xa0)/[a0*(b0xc0)]
!         c = r*(a0xb0)/[a0*(b0xc0)]
!     where a0*(b0xc0) = fac.
!
!     Determine the projection of the furtherest corner of the atomic 
!     cell onto the cell axes with the present atom as the origin. 
!     The strategy is described as follows. 
!     1. The projection of a point r onto a0, b0, and c0 is given by 
!            r*(b0xc0)/[a0*(b0xc0)]
!            r*(c0xa0)/[a0*(b0xc0)]
!            r*(a0xb0)/[a0*(b0xc0)]
!        respectively.
!     2. If r is on a sphere of radius R, the largest projection distance
!        from the atom center is
!            R*(b0xc0)**2/[a0*(b0xc0)]
!            R*(c0xa0)**2/[a0*(b0xc0)]
!            R*(a0xb0)**2/[a0*(b0xc0)]
!     3. The projection of the atom position p onto the three axis is given by
!            a = p*(b0xc0)/[a0*(b0xc0)] 
!            b = p*(c0xa0)/[a0*(b0xc0)] 
!            c = p*(a0xb0)/[a0*(b0xc0)] 
!        respectively
!     4. The uniform grid points associated with the atom is contained in 
!        AtomBox as described below:
!            AtomBox(1:2,1,atom): the index range in a0 dimension
!            AtomBox(1:2,2,atom): the index range in b0 dimension
!            AtomBox(1:2,3,atom): the index range in c0 dimension
!     5. For the index range in each dimension is the number of grid points 
!        in between
!            [a-R*(b0xc0)**2/[a0*(b0xc0)], a+R*(b0xc0)**2/[a0*(b0xc0)]]
!            [b-R*(c0xa0)**2/[a0*(b0xc0)], b+R*(c0xa0)**2/[a0*(b0xc0)]]
!            [c-R*(a0xb0)**2/[a0*(b0xc0)], c+R*(a0xb0)**2/[a0*(b0xc0)]]
!     ================================================================
      rv = getCrossProduct(b0,c0)
      a = getDotProduct(3,AtomPosition(:,ia),rv,absolute=.true.)/fac
      bxc = getVecLength(3,rv)
      max_rc = radius(ia)*bxc**2/fac
!     write(6,'(a,3d18.8)')'a,max_rc,step_a_len = ',a,max_rc,step_a_len
      pAOG%AtomBox(1,1,ia) = floor((a-max_rc)/step_a_len)
      pAOG%AtomBox(2,1,ia) = ceiling((a+max_rc)/step_a_len)
!     write(6,'(a,2f15.8,2x,2i8)')'AtomBox(1:2,1) = ',a-max_rc,a+max_rc,pAOG%AtomBox(1:2,1,ia)
!
      rv = getCrossProduct(c0,a0)
      b = getDotProduct(3,AtomPosition(:,ia),rv,absolute=.true.)/fac
      cxa = getVecLength(3,rv)
      max_rc = radius(ia)*cxa**2/fac
!     write(6,'(a,3d18.8)')'b,max_rc,step_b_len = ',b,max_rc,step_b_len
      pAOG%AtomBox(1,2,ia) = floor((b-max_rc)/step_b_len)
      pAOG%AtomBox(2,2,ia) = ceiling((b+max_rc)/step_b_len)
!     write(6,'(a,2f15.8,2x,2i8)')'AtomBox(1:2,2) = ',b-max_rc,b+max_rc,pAOG%AtomBox(1:2,2,ia)
!
      rv = getCrossProduct(a0,b0)
      c = getDotProduct(3,AtomPosition(:,ia),rv,absolute=.true.)/fac
      axb = getVecLength(3,rv)
      max_rc = radius(ia)*axb**2/fac
!     write(6,'(a,3d18.8)')'c,max_rc,step_c_len = ',c,max_rc,step_c_len
      pAOG%AtomBox(1,3,ia) = floor((c-max_rc)/step_c_len)
      pAOG%AtomBox(2,3,ia) = ceiling((c+max_rc)/step_c_len)
!     write(6,'(a,2f15.8,2x,2i8)')'AtomBox(1:2,3) = ',c-max_rc,c+max_rc,pAOG%AtomBox(1:2,3,ia)
!
      ic = pAOG%AtomBox(2,1,ia)-pAOG%AtomBox(1,1,ia)+1
      jc = pAOG%AtomBox(2,2,ia)-pAOG%AtomBox(1,2,ia)+1
      kc = pAOG%AtomBox(2,3,ia)-pAOG%AtomBox(1,3,ia)+1
      pAOG%NumGridPointsInAtomBox(ia) = ic*jc*kc
!
      rv = getCrossProduct(pUG%grid_step_a,pUG%grid_step_b)
      vol = getDotProduct(3,pUG%grid_step_c,rv,absolute=.true.)*(ic-1)*(jc-1)*(kc-1)
      frac = pAOG%NumGridPointsInAtomBox(ia)*PI*radius(ia)**3/vol
!     write(6,'(a,2i8)')'Number of grid points in AtomBox and cell = ', &
!                       pAOG%NumGridPointsInAtomBox(ia),ceiling(frac)
      max_incell_pts = max(max_incell_pts,ceiling(frac)) ! A rough estimate of the number of grid points in atomic cell
   enddo
!
   allocate(pAOG%InCellGridPointABIndex(max_incell_pts,NumAtoms))
!
   pAOG%NumGridPointsOnCellBound = 0
   do ia = 1, NumAtoms
      pAOG%NumGridPointsInCell(ia) = 0
      do ig = 1, pAOG%NumGridPointsInAtomBox(ia)
         gCounter = getGridIndex(pUG,ia,ig) ! Given the grid index defined in box,
                                            ! get the corresponding global index.
         rg = getGridPosition(pUG,ia,ig)    ! Notes on 8/22/18: this function should
                                            ! return the grid position relative to 
                                            ! the system origin.
         rv = rg - AtomPosition(1:3,ia)
         if ( radius(ia) > sqrt(rv(1)**2 + rv(2)**2 + rv(3)**2) + Ten2m6 ) then
!           ==========================================================
!           The grid point is inside the "region", but it may be
!           outside of the atomic cell, if for example the radius of
!           the region is larger than the atomic cell bounding sphere.
!
!           Note: In this implementation, a situation could arise
!                 that two equivalent uniform grid points (with 
!                 the same gCounter index) may appear in the same
!                 atomic cell.
!
!           LocationFlag = -1: grid point is outside the cell
!                        =  0: grid point is on the cell surface
!                        =  1: grid point is inside the cell
!
!           Establish a linked list of all grid points which are on atomic
!           cell boundaries.
!           ==========================================================
            n = getPointLocationFlag(ia, rv(1), rv(2), rv(3))
            if (n >= 0) then
               pAOG%NumGridPointsInCell(ia) = pAOG%NumGridPointsInCell(ia) + 1
               if (pAOG%NumGridPointsInCell(ia) > max_incell_pts) then
!                 ----------------------------------------------------
                  call ErrorHandler('insertAtomsInGrid',              &
                                    'Number of grid points > max estimated value', &
                                    pAOG%NumGridPointsInCell(ia),max_incell_pts)
!                 ----------------------------------------------------
               else
                  pAOG%InCellGridPointABIndex(pAOG%NumGridPointsInCell(ia),ia) = ig
               endif
            endif
            if (n == 0) then
               pAOG%NumGridPointsOnCellBound = pAOG%NumGridPointsOnCellBound + 1
               if (.not.associated(p_head)) then
                  allocate(p_head)
                  linked_list => p_head
               else
                  allocate(linked_list%next)
                  linked_list => linked_list%next
               endif
               nullify(linked_list%next)
               linked_list%grid_index = gCounter
            endif
         endif
      enddo
   enddo
!
!  ===================================================================
!  Copy the index data from the linked list to pAOG%CBGridPointIndex
!  ===================================================================
   n = pAOG%NumGridPointsOnCellBound
   if (n > 0) then
!     ================================================================
!     Eliminate the redundant grid points on cell boundaries that are 
!     equivalent, corresponding to the same global grid index.
!     ================================================================
      allocate(gtmp(n), ntmp(n))
      linked_list => p_head
      m = 0
      do i = 1, n
         redundant = .false.
         LOOP_j: do j = 1, m
            k = j
            if (gtmp(j) == linked_list%grid_index) then
               redundant = .true.
               exit LOOP_j
            endif
         enddo LOOP_j
         if (.not.redundant) then
            m = m + 1
            gtmp(m) = linked_list%grid_index
            ntmp(m) = 1
         else ! For redundant grid points, increment ntmp by 1
            pAOG%NumGridPointsOnCellBound = pAOG%NumGridPointsOnCellBound - 1
            ntmp(k) = ntmp(k) + 1
         endif
         linked_list => linked_list%next
         nullify(p_head%next)
         deallocate(p_head)
         p_head => linked_list
      enddo
!
!     ================================================================
!     At this point, m = the number of non-redundant grid points on the
!                        cell boundaries
!                    gtmp = the global index of each grid point on the
!                           cell boundaries
!                    ntmp = the redundancy of each grid point on the
!                           cell boundaries
!     ================================================================
      if (m /= pAOG%NumGridPointsOnCellBound) then
         call ErrorHandler('insertAtomsInGrid','m <> NumGridPointOnCellBound', &
                           m, pAOG%NumGridPointsOnCellBound)
      endif
      allocate(pAOG%CBGridPointIndex(m), pAOG%NumNeighboringAtoms(m))
      pAOG%CBGridPointIndex(1:m) = gtmp(1:m)
      pAOG%NumNeighboringAtoms(1:m) = ntmp(1:m)
      deallocate(gtmp, ntmp)
   endif
!
   nullify(p_head, linked_list)
!
!  ===================================================================
!  Collecting pAOG%CBGridPointIndex from other processors and update
!  pAOG%NumNeighboringAtoms
!  ===================================================================
   aGID = getGroupID('Unit Cell')
   NumPEsInAGroup = getNumPEsInGroup(aGID)
   if (NumPEsInAGroup > 1) then
      MyPEinAGroup = getMyPEinGroup(aGID)
      n = pAOG%NumGridPointsOnCellBound
!     ----------------------------------------------------------------
      call GlobalMaxInGroup(aGID,n)
!     ----------------------------------------------------------------
      if (n > 0) then
         allocate(size_collect(NumPEsInAGroup), data_collect(2,n,NumPEsInAGroup))
         size_collect = 0
         data_collect = 0
         m = pAOG%NumGridPointsOnCellBound
         size_collect(MyPEinAGroup+1) = m
         do i = 1, m
            data_collect(1,i,MyPEinAGroup+1)= pAOG%CBGridPointIndex(i)
            data_collect(2,i,MyPEinAGroup+1)= pAOG%NumNeighboringAtoms(i)
         enddo
!        -------------------------------------------------------------
         call GlobalCollectInGroup(aGID,size_collect)
         call GlobalCollectInGroup(aGID,data_collect,2,n)
!        -------------------------------------------------------------
         do ip = 1, NumPEsInAGroup
            if (ip /= MyPEinAGroup+1) then
               do j = 1, size_collect(ip)
                  ig = data_collect(1,j,ip)
                  k = data_collect(2,j,ip)
                  do i = 1, m
                     if (pAOG%CBGridPointIndex(i) == ig) then
                        pAOG%NumNeighboringAtoms(i) =                 &
                                       pAOG%NumNeighboringAtoms(i) + k
                     endif
                  enddo
               enddo
            endif
         enddo
         deallocate(size_collect,data_collect)
      endif
   endif
!
!  ===================================================================
!  Check if the number of uniform grid points in all atomic cells equals
!  to the number of uniform grid points in unit cell.
!  The following method for checking the consistency fails in multiple 
!  atoms case. It needs to be redesigned.
!  ===================================================================
!! n = 0
!! do ia = 1, NumAtoms
!!    n = n + pAOG%NumGridPointsInCell(ia)
!! enddo
!! do i = 1, pAOG%NumGridPointsOnCellBound
!!    n = n - pAOG%NumNeighboringAtoms(i) + 1
!! enddo
!  -------------------------------------------------------------------
!! call GlobalSumInGroup(aGID,n)
!  -------------------------------------------------------------------
!! if (n /= pUG%ng) then
!!    call ErrorHandler('insertAtomsInGrid','Failed integraty test: n <> ng',n,pUG%ng)
!! endif
!
!  ===================================================================
!  Determine the relation of each atom box with the MPI tasks, which the
!  uniform grid points in the atom box are mapped onto. Then establish
!  the send and receive list for the atom box.
!  -------------------------------------------------------------------
   call setupSRList(pUG)
!  -------------------------------------------------------------------
   nullify(pAOG)
   nullify(pUG)
!
   end subroutine insertAtomsInGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupSRList(pUG)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalCollectInGroup
   implicit none
!
   type (UniformGridStruct), intent(inout) :: pUG
   type (AtomOnUniformGridStruct), pointer :: pAOG
!
   integer (kind=IntKind) :: aGID, NumPEsInAGroup, MyPEinAGroup
   integer (kind=IntKind), allocatable :: GridParallelInfo(:,:,:)
   integer (kind=IntKind) :: ia, na, npt, nps, i, ip
!
   type Grid2ProcStruct
      integer (kind=IntKind) :: proc, atom
      type (Grid2ProcStruct), pointer :: next
   end type Grid2ProcStruct
!
   type (Grid2ProcStruct), pointer :: p_head, linked_list
!
   pAOG => pUG%AtomOnGrid
!
!  ===================================================================
!  The send and receive data activities only need to take place within 
!  each 'Unit Cell' group.
!  ===================================================================
   aGID = getGroupID('Unit Cell')
   NumPEsInAGroup = getNumPEsInGroup(aGID)
   MyPEinAGroup = getMyPEinGroup(aGID)
   na = pAOG%NumLocalAtoms
!
   allocate(GridParallelInfo(6,1+na,NumPEsInAGroup))
!
   GridParallelInfo = 0
   GridParallelInfo(1,1,MyPEinAGroup+1)= pUG%gstart(1)
   GridParallelInfo(2,1,MyPEinAGroup+1)= pUG%gend(1)
   GridParallelInfo(3,1,MyPEinAGroup+1)= pUG%gstart(2)
   GridParallelInfo(4,1,MyPEinAGroup+1)= pUG%gend(2)
   GridParallelInfo(5,1,MyPEinAGroup+1)= pUG%gstart(3)
   GridParallelInfo(6,1,MyPEinAGroup+1)= pUG%gend(3)
   do ia = 1, na
      GridParallelInfo(1,ia+1,MyPEinAGroup+1)= pAOG%AtomBox(1,1,ia)
      GridParallelInfo(2,ia+1,MyPEinAGroup+1)= pAOG%AtomBox(2,1,ia)
      GridParallelInfo(3,ia+1,MyPEinAGroup+1)= pAOG%AtomBox(1,2,ia)
      GridParallelInfo(4,ia+1,MyPEinAGroup+1)= pAOG%AtomBox(2,2,ia)
      GridParallelInfo(5,ia+1,MyPEinAGroup+1)= pAOG%AtomBox(1,3,ia)
      GridParallelInfo(6,ia+1,MyPEinAGroup+1)= pAOG%AtomBox(2,3,ia)
   enddo
!
!  -------------------------------------------------------------------
   call GlobalCollectInGroup(aGID,GridParallelInfo,6,na+1)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  determine the send list.
!  ===================================================================
   allocate(pAOG%NumTargetProcs(na))
   nullify(p_head)
   npt = 0
   do ia = 1, na
      pAOG%NumTargetProcs(ia) = 0
      do ip = 1, NumPEsInAGroup
         if (ip /= MyPEinAGroup+1) then
            if ( isOverlap(pUG%nga,pUG%ngb,pUG%ngc,GridParallelInfo(:,1,ip), &
                           GridParallelInfo(:,ia+1,MyPEinAGroup+1)) ) then
               pAOG%NumTargetProcs(ia) = pAOG%NumTargetProcs(ia) + 1
               if (.not.associated(p_head)) then
                  allocate(p_head)
                  linked_list => p_head
               else
                  allocate(linked_list%next)
                  linked_list => linked_list%next
               endif
               nullify(linked_list%next)
               linked_list%proc = ip-1
               linked_list%atom = ia
            endif
         endif
      enddo
      npt = max(npt,pAOG%NumTargetProcs(ia))
   enddo
   if (npt > 0) then
      allocate(pAOG%TargetProc(npt,na))
      pAOG%TargetProc = -3  ! An arbitrary negative number
      linked_list => p_head
      do ia = 1, na
         do ip = 1, pAOG%NumTargetProcs(ia)
            pAOG%TargetProc(ip,ia) = linked_list%proc
            if (ia /= linked_list%atom) then ! Checking data consistency
               call ErrorHandler('setupSRList','Inconsistent sending list content', &
                                 ia,linked_list%atom)
            endif
            linked_list => linked_list%next
            nullify(p_head%next)
            deallocate(p_head)
            p_head => linked_list
         enddo
      enddo
   endif
!
!  ===================================================================
!  determine the receive list.
!  ===================================================================
   allocate(pAOG%NumSourceProcs(na))
   nullify(p_head)
   nps = 0
   do ia = 1, na
      pAOG%NumSourceProcs(ia) = 0
      do ip = 1, NumPEsInAGroup
         if (ip /= MyPEinAGroup+1) then
            if ( isOverlap(pUG%nga,pUG%ngb,pUG%ngc,GridParallelInfo(:,1,MyPEinAGroup+1), &
                           GridParallelInfo(:,ia+1,ip)) ) then
               pAOG%NumSourceProcs(ia) = pAOG%NumSourceProcs(ia) + 1
               if (.not.associated(p_head)) then
                  allocate(p_head)
                  linked_list => p_head
               else
                  allocate(linked_list%next)
                  linked_list => linked_list%next
               endif
               nullify(linked_list%next)
               linked_list%proc = ip-1
               linked_list%atom = ia
            endif
         endif
      enddo
      nps = max(nps,pAOG%NumSourceProcs(ia))
   enddo
   if (nps > 0) then
      allocate(pAOG%SourceProc(nps,na))
      pAOG%SourceProc = -3 ! An arbitrary negative number
      linked_list => p_head
      do ia = 1, na
         do ip = 1, pAOG%NumSourceProcs(ia)
            pAOG%SourceProc(ip,ia) = linked_list%proc
            if (ia /= linked_list%atom) then ! Checking data consistency
               call ErrorHandler('setupSRList','Inconsistent receiving list content',ia,linked_list%atom)
            endif
            linked_list => linked_list%next
            nullify(p_head%next)
            deallocate(p_head)
            p_head => linked_list
         enddo
      enddo
   endif
!
   deallocate(GridParallelInfo)
   nullify(pAOG)
!
   end subroutine setupSRList
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isOverlap(nga,ngb,ngc,Grid,AtomBox) result(t)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nga, ngb, ngc
   integer (kind=IntKind), intent(in) :: Grid(6), AtomBox(6)
   integer (kind=IntKind) :: x0, x1, y0, y1, z0, z1
   integer (kind=IntKind) :: i, j, j0, k, n
!
   integer (kind=IntKind), parameter :: shift = 1
   integer (kind=IntKind) :: num_shifts_1d = 2*shift+1
   integer (kind=IntKind) :: num_shifts_2d = (2*shift+1)**2
   integer (kind=IntKind) :: num_shifts_3d = (2*shift+1)**3
!
   logical :: t
!
   t = .false.
!
   LOOP_n: do n = 1, num_shifts_3d
      i = mod(n,num_shifts_1d)-1
      j = mod((n-i-1)/num_shifts_1d+1,num_shifts_1d)-1
      j0 = mod((n-i-1)/num_shifts_1d,num_shifts_1d)
      k = mod((n-j0*num_shifts_1d-i-1)/num_shifts_2d+1,num_shifts_1d)-1
!
      x0 = Grid(1)+i*nga; x1 = Grid(2)+i*nga
      y0 = Grid(3)+j*ngb; y1 = Grid(4)+j*ngb
      z0 = Grid(5)+k*ngc; z1 = Grid(6)+k*ngc
!     ================================================================
!     Check if there are overlap between Grid box and AtomBox along 
!     a, b, c dimensions
!     ================================================================
      if (isOverlapIn1D(x0,x1,AtomBox(1),AtomBox(2)) .and.            &
          isOverlapIn1D(y0,y1,AtomBox(3),AtomBox(4)) .and.            &
          isOverlapIn1D(z0,z1,AtomBox(5),AtomBox(6))) then
         t = .true.
         exit LOOP_n
      endif
   enddo LOOP_n
!
   end function isOverlap
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isOverlapIn1D(a0,a1,b0,b1) result (t)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: a0, a1, b0, b1
!
   logical :: t
!
   t = .false.
!
   if ((a0 <= b0 .and. b0 <= a1) .or.                                 &
       (a0 <= b1 .and. b1 <= a1) .or.                                 &
       (b0 <= a0 .and. a0 <= b1) .or.                                 &
       (b0 <= a1 .and. a1 <= b1)) then
      t = .true.
   endif
!
   end function isOverlapIn1D
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInBox(x,y,z,Box) result (t)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: x, y, z, Box(6)
   integer (kind=IntKind) :: x0, x1, y0, y1, z0, z1
!
   logical :: t
!
   x0 = Box(1); x1 = Box(2)
   y0 = Box(3); y1 = Box(4)
   z0 = Box(5); z1 = Box(6)
!
   t = .false.
   if (x0 <= x .and. x <= x1 .and.                            &
       y0 <= y .and. y <= y1 .and.                            &
       z0 <= z .and. z <= z1) then
      t = .true.
   endif
!
   end function isInBox
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridIndex_box(pUG,ia,ig_in_box) result (grid_index)
!  ===================================================================
!
!  Given the index, ig_in_box, of a uniform grid point in AtomBox of 
!  atom ia, determine the corresponding global index of the grid point.
!
!  *******************************************************************
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in) :: ia, ig_in_box
!
   integer (kind=IntKind) :: ic, jc, kc, nk, nj, ni, icp, jcp, kcp
   integer (kind=IntKind) :: grid_index
!
   type (AtomOnUniformGridStruct), pointer :: pAOG
!
   pAOG => pUG%AtomOnGrid
!
   ni = pAOG%AtomBox(2,1,ia)-pAOG%AtomBox(1,1,ia)+1
   nj = pAOG%AtomBox(2,2,ia)-pAOG%AtomBox(1,2,ia)+1
   nk = pAOG%AtomBox(2,3,ia)-pAOG%AtomBox(1,3,ia)+1
!
   if (ig_in_box < 1 .or. ig_in_box > ni*nj*nk) then
      call ErrorHandler('getGridIndex_box','ig_in_box out of range',ig_in_box)
   endif
!
   ic = mod(ig_in_box-1,ni)+pAOG%AtomBox(1,1,ia)
   jc = mod((ig_in_box-ic+pAOG%AtomBox(1,1,ia)-1)/ni,nj)+pAOG%AtomBox(1,2,ia)
   kc = ((ig_in_box-ic+pAOG%AtomBox(1,1,ia)-1)/ni-(jc-pAOG%AtomBox(1,2,ia)))/nj+pAOG%AtomBox(1,3,ia)
!
   kcp = kc
   do while (kcp < 1)
      kcp = kcp + pUG%ngc
   enddo
   do while (kcp > pUG%ngc)
      kcp = kcp - pUG%ngc
   enddo
!
   jcp = jc
   do while (jcp < 1)
      jcp = jcp + pUG%ngb
   enddo
   do while (jcp > pUG%ngb)
      jcp = jcp - pUG%ngb
   enddo
!
   icp = ic
   do while (icp < 1)
      icp = icp + pUG%nga
   enddo
   do while (icp > pUG%nga)
      icp = icp - pUG%nga
   enddo
!
   grid_index = (kcp-1)*pUG%ngb*pUG%nga+(jcp-1)*pUG%nga+icp
!
   end function getGridIndex_box
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridIndex_local(pUG,ig_local) result (grid_index)
!  ===================================================================
!
!  Given the local index, ig_local, of a uniform grid point on local 
!  processor, determine the corresponding global index of the grid point
!
!  *******************************************************************
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in) :: ig_local
   integer (kind=IntKind) :: grid_index
   integer (kind=IntKind) :: ni, nj, nk, i, j, k
!
   ni = pUG%gend(1)-pUG%gstart(1)+1
   nj = pUG%gend(2)-pUG%gstart(2)+1
   nk = pUG%gend(3)-pUG%gstart(3)+1
!
   if (ig_local < 1 .or. ig_local > ni*nj*nk) then
      call ErrorHandler('getGridIndex_local','ig_local out of range',ig_local)
   endif
!
   i = mod(ig_local-1,ni)+pUG%gstart(1)
   j = mod((ig_local-i+pUG%gstart(1)-1)/ni,nj)+pUG%gstart(2)
   k = ((ig_local-i+pUG%gstart(1)-1)/ni-(j-pUG%gstart(2)))/nj+pUG%gstart(3)
!
   grid_index = (k-1)*pUG%ngb*pUG%nga + (j-1)*pUG%nga + i
!
   end function getGridIndex_local
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridPosition_box(pUG,ia,ig_in_box) result (pos)
!  ===================================================================
!
!  Given the index, ig_in_box, of a uniform grid point in AtomBox of 
!  atom ia, determine the position of the grid point relative to the
!  system origin, which sits on the grid point with its global index = 1.
!
!  *******************************************************************
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in) :: ia, ig_in_box
!
   integer (kind=IntKind) :: ic, jc, kc, nk, nj, ni
!
   real (kind=RealKind) :: pos(3)
!
   type (AtomOnUniformGridStruct), pointer :: pAOG
!
   pAOG => pUG%AtomOnGrid
!
   ni = pAOG%AtomBox(2,1,ia)-pAOG%AtomBox(1,1,ia)+1
   nj = pAOG%AtomBox(2,2,ia)-pAOG%AtomBox(1,2,ia)+1
   nk = pAOG%AtomBox(2,3,ia)-pAOG%AtomBox(1,3,ia)+1
   ic = mod(ig_in_box-1,ni)+pAOG%AtomBox(1,1,ia)
   jc = mod((ig_in_box-ic+pAOG%AtomBox(1,1,ia)-1)/ni,nj)+pAOG%AtomBox(1,2,ia)
   kc = ((ig_in_box-ic+pAOG%AtomBox(1,1,ia)-1)/ni-(jc-pAOG%AtomBox(1,2,ia)))/nj+pAOG%AtomBox(1,3,ia)
   pos = pUG%grid_step_a*(ic-1) + pUG%grid_step_b*(jc-1) + pUG%grid_step_c*(kc-1)
!
   end function getGridPosition_box
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridPosition_global(pUG,ig_global) result (pos)
!  ===================================================================
!
!  Given the index, ig_global, of a uniform grid point in unit cell,
!  determine the position of the grid point relative to the system 
!  origin, which sits on the grid point with its global index = 1.
!
!  *******************************************************************
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in) :: ig_global
!
   integer (kind=IntKind) :: ic, jc, kc
!
   real (kind=RealKind) :: pos(3)
!
   ic = mod(ig_global-1,pUG%nga)+1
   jc = mod((ig_global-ic)/pUG%nga,pUG%ngb)+1
   kc = ((ig_global-ic)/pUG%nga-(jc-1))/pUG%ngb+1
   pos = pUG%grid_step_a*(ic-1) + pUG%grid_step_b*(jc-1) + pUG%grid_step_c*(kc-1)
!
   end function getGridPosition_global
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getUniform3DGrid(gname) result(pUG)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   type (UniformGridStruct), pointer :: pUG
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getUniform3DGrid','need to call initUniform3DGrid first')
!     ----------------------------------------------------------------
   endif
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('getUniform3DGrid','Unknown grid name',gname)
   endif
!
   end function getUniform3DGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printUniform3DGrid(gname)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   type (UniformGridStruct), pointer :: pUG
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('printUniform3DGrid','need to call initUniform3DGrid first')
!     ----------------------------------------------------------------
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')' ************************************'
   write(6,'(24x,a)')  ' *  Output from printUniform3DGrid  *'
   write(6,'(24x,a,/)')' ************************************'
   write(6,'(a)')'============================================================'
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
      write(6,'(a)')'FFT Grid Info:'
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
      write(6,'(a)')'Visual Grid Info:'
   else
      call ErrorHandler('printUniform3DGrid','Unknown grid name',gname)
   endif
!
   write(6,'(  ''cell side a '',t40,''='',3d13.6)')pUG%cell(1:3,1)
   write(6,'(  ''cell side b '',t40,''='',3d13.6)')pUG%cell(1:3,2)
   write(6,'(  ''cell side c '',t40,''='',3d13.6)')pUG%cell(1:3,3)
   write(6,'(  ''nga        '',t40,''='',i5)')     pUG%nga
   write(6,'(  ''ngb        '',t40,''='',i5)')     pUG%ngb 
   write(6,'(  ''ngc        '',t40,''='',i5)')     pUG%ngc
   write(6,'(  ''number of grid points'',t40,''='',i10)') pUG%ng
   write(6,'(  ''step_a     '',t40,''='',3d13.6)')  pUG%grid_step_a(1:3)
   write(6,'(  ''step_b     '',t40,''='',3d13.6)')  pUG%grid_step_b(1:3)
   write(6,'(  ''step_c     '',t40,''='',3d13.6)')  pUG%grid_step_c(1:3)
   write(6,'(a)')'============================================================'
!
   end subroutine printUniform3DGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGridPoints(gname,dim,local_grid) result(n)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind), optional, intent(in) :: dim
   logical, optional, intent(in) :: local_grid
   integer (kind=IntKind) :: n
   logical :: local
!
   type (UniformGridStruct), pointer :: pUG
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumGridPoints','need to call initUniform3DGrid first')
!     ----------------------------------------------------------------
   endif
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('getUniform3DGrid','Unknown grid name',gname)
   endif
!
   if (present(local_grid)) then
      local = local_grid
   else
      local = .false.
   endif
!
   if (local) then
      if (present(dim)) then
         n = pUG%gend(dim)-pUG%gstart(dim)+1
      else
         n = pUG%NumLocalGridPoints
      endif
   else
      if (present(dim)) then
         if (dim == 1) then
            n = pUG%nga
         else if (dim == 2) then
            n = pUG%ngb
         else if (dim == 3) then
            n = pUG%ngc
         else
!           ----------------------------------------------------------
            call ErrorHandler('getNumGridPoints','Invalid dimension',dim)
!           ----------------------------------------------------------
         endif
      else
         n = pUG%ng
      endif
   endif
!
   end function getNumGridPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridStepSize(gname,dim) result(h)
!  ===================================================================
   use VectorModule, only : getVecLength
!
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind), intent(in) :: dim
!
   real (kind=RealKind) :: h
!
   type (UniformGridStruct), pointer :: pUG
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGridStepSize','need to call initUniform3DGrid first')
!     ----------------------------------------------------------------
   endif
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('getUniform3DGrid','Unknown grid name',gname)
   endif
!
   if (dim == 1) then
      h = getVecLength(3,pUG%grid_step_a)
   else if (dim == 2) then
      h = getVecLength(3,pUG%grid_step_b)
   else if (dim == 3) then
      h = getVecLength(3,pUG%grid_step_c)
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getGridStepSize','Invalid dimension',dim)
!     ----------------------------------------------------------------
   endif
!
   end function getGridStepSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isOnAtomicCellBoundary_s(gname, id, num_neighbor, local_index) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: gname
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out), optional :: num_neighbor
   logical, intent(in), optional :: local_index
!
   logical :: p
!
   integer (kind=IntKind) :: ig, n, gCounter
!
   type (UniformGridStruct), pointer :: pUG
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('isOnAtomicCellBoundary',                     &
                        'need to call initUniform3DGrid first')
!     ----------------------------------------------------------------
   endif
!
   if (nocaseCompare(gname,'FFT')) then
      pUG => FFTGrid
   else if (nocaseCompare(gname,'Visual')) then
      pUG => VisualGrid
   else
      call ErrorHandler('isOnAtomicCellBoundary','Unknown grid name',gname)
   endif
!
   if (present(local_index)) then
      if (local_index) then  ! In case id is a local grid index, it needs to be
                             ! converted to global index by calling getGridIndex
         gCounter = getGridIndex(pUG,id)
      else
         gCounter = id
      endif
   else
      gCounter = id
   endif
!
   if (gCounter < 1 .or. gCounter > pUG%ng) then
      call ErrorHandler('isOnAtomicCellBoundary','gCounter is out of range',gCounter)
   endif
!
   p = .false.
   LOOP_ig: do ig = 1, pUG%AtomOnGrid%NumGridPointsOnCellBound
      n = ig
      if (pUG%AtomOnGrid%CBGridPointIndex(ig) == gCounter) then
         p = .true.
         exit LOOP_ig
      endif
   enddo LOOP_ig
!
   if (present(num_neighbor)) then
      num_neighbor = pUG%AtomOnGrid%NumNeighboringAtoms(n)
   endif
!
   end function isOnAtomicCellBoundary_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isOnAtomicCellBoundary_g(pUG, id, num_neighbor, local_index) result(p)
!  ===================================================================
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out), optional :: num_neighbor
!
   logical, intent(in), optional :: local_index
   logical :: p
!
   integer (kind=IntKind) :: ig, n, gCounter
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('isOnAtomicCellBoundary',                     &
                        'need to call initUniform3DGrid first')
!     ----------------------------------------------------------------
   endif
!
   if (present(local_index)) then
      if (local_index) then  ! In case id is a local grid index, it needs to be
                             ! converted to global index by calling getGridIndex
         gCounter = getGridIndex(pUG,id)
      else
         gCounter = id
      endif
   else
      gCounter = id
   endif
!
   if (gCounter < 1 .or. gCounter > pUG%ng) then
      call ErrorHandler('isOnAtomicCellBoundary','gCounter is out of range',gCounter)
   endif
!
   p = .false.
   LOOP_ig: do ig = 1, pUG%AtomOnGrid%NumGridPointsOnCellBound
      n = ig
      if (pUG%AtomOnGrid%CBGridPointIndex(ig) == gCounter) then
         p = .true.
         exit LOOP_ig
      endif
   enddo LOOP_ig
!
   if (present(num_neighbor)) then
      num_neighbor = pUG%AtomOnGrid%NumNeighboringAtoms(n)
   endif
!
   end function isOnAtomicCellBoundary_g
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLocalGrid(pUG,ig,id) result(p)
!  ===================================================================
!
!  This function returns true or false value of an inquiry:
!       where the grid point of global index ig is mapped on the local 
!       process, as a result of parallization.
!  The returning value of variable "id" corresponds to the local index
!  of the grid point if it is indeed local to the current process. 
!
!  Note: id is NOT the index of the grid points in AtomBox, rather it is
!        the local index of the grid point on the local process.
!  *******************************************************************
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in)  :: ig ! global index of a grid point
   integer (kind=IntKind), intent(out) :: id ! local index of a grid point
!
   logical :: p
!
   integer (kind=IntKind) :: i, j, k
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('isLocalGrid','need to call initUniform3DGrid first')
!     ----------------------------------------------------------------
   else if (ig < 1 .or. ig > pUG%ng) then
!     ----------------------------------------------------------------
      call ErrorHandler('isLocalGrid','ig is out of range',ig)
!     ----------------------------------------------------------------
   endif
!
   i = mod(ig-1,pUG%nga)+1
   j = mod((ig-i)/pUG%nga,pUG%ngb)+1
   k = ((ig-i)/pUG%nga-(j-1))/pUG%ngb+1
!
   if (pUG%gstart(1) > i .or. i > pUG%gend(1) .or.                    &
       pUG%gstart(2) > j .or. j > pUG%gend(2) .or.                    &
       pUG%gstart(3) > k .or. k > pUG%gend(3)) then
      p = .false.
      id = 0
   else
      p = .true.
      id = (k-pUG%gstart(3))*(pUG%gend(2)-pUG%gstart(2)+1)*           &
                             (pUG%gend(1)-pUG%gstart(1)+1) +          &
           (j-pUG%gstart(2))*(pUG%gend(1)-pUG%gstart(1)+1) + i-pUG%gstart(1)+1
   endif
!
   end function isLocalGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSourceProc(pUG,id,nProc) result(p)
!  ===================================================================
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: nProc
   integer (kind=IntKind), pointer :: p(:)
!
   if (id < 1 .or. id > pUG%AtomOnGrid%NumLocalAtoms) then
!     ----------------------------------------------------------------
      call ErrorHandler('getSourceProc','atom index is out of range',id)
!     ----------------------------------------------------------------
   endif
!
   nProc = pUG%AtomOnGrid%NumSourceProcs(id)
   if (nProc > 0) then
      p => pUG%AtomOnGrid%SourceProc(:,id)
   else
      nullify(p)
   endif
!
   end function getSourceProc
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTargetProc(pUG,id,nProc) result(p)
!  ===================================================================
   implicit none
!
   type (UniformGridStruct), intent(in) :: pUG
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: nProc
   integer (kind=IntKind), pointer :: p(:)
!
   if (id < 1 .or. id > pUG%AtomOnGrid%NumLocalAtoms) then
!     ----------------------------------------------------------------
      call ErrorHandler('getTargetProc','atom index is out of range',id)
!     ----------------------------------------------------------------
   endif
!
   nProc = pUG%AtomOnGrid%NumTargetProcs(id)
   if (nProc > 0) then
      p => pUG%AtomOnGrid%TargetProc(:,id)
   else
      nullify(p)
   endif
!
   end function getTargetProc
!  ===================================================================
end module Uniform3DGridModule
