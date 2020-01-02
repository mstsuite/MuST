module PolyhedraModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
   use MathParamModule, only : ZERO, THIRD, HALF, ONE, TEN2m6, TEN2m8, TEN2m10, TEN2m12
   use Matrix3dModule, only : invm3
!
public :: initPolyhedra,        &
          genPolyhedron,        &
          endPolyhedra,         &
          getNumPolyhedra,      &
          getNeighborIndex,     &
          getSurfaceArea,       &
          getNumPlanes,         &
          getNumPlaneCorners,   &
          getPlaneCorner,       &
          getPlaneArea,         &
          getPlaneNormal,       &
          getInscrSphRadius,    &
          getOutscrSphRadius,   &
          getWignerSeitzRadius, &
          getVolume,            &
          getInscrSphVolume,    &
          getBoxVolume,         &
          getNeighborDistance,  &
          getNumCorners,        &
          getNumEdges,          &
          getPlane,             &
          getCorner,            &
          getEdge,              &
          getIndexCornerPlane,  &
          isExternalPoint,      &
          isSurfacePoint,       &
          getPointLocationFlag, &
          printPolyhedron,      &
          printPolyhedronBoundary, &
          printPolyhedraTable
!
   interface genPolyhedron
      module procedure genPoly1, genPoly2
   end interface
!
private
   type PolyhedronStruct
      integer (kind=IntKind) :: NumPlanes  ! num. of cell boundary planes
      integer (kind=IntKind) :: NumCorners ! num. of polyhedron corners
      integer (kind=IntKind) :: NumEdges   ! num. of polyhedron edges
      integer (kind=IntKind) :: GlobalIndex! the global index of the polyhedron
!
      integer (kind=IntKind), pointer :: id_planes(:) ! index of the 
                                                      ! neighboring atom 
                                                      ! who shares the 
                                                      ! same boundary plane
!
      integer (kind=IntKind), pointer :: index_corner_plane(:,:)
         ! ===========================================================
         ! index_corner_plane: gives the index of polyhedron corner for each 
         ! vertex on a boundary plane
         ! ===========================================================
!
      integer (kind=IntKind), pointer :: NumPlaneCorners(:)
         ! ===========================================================
         ! num. of vertices on a boundary plane
         ! ===========================================================
!
      real (kind=RealKind), pointer :: vplane(:,:)
         ! ===========================================================
         ! vplane: cell boundary plane vectors, which is the vector normal to 
         ! the plane. The magnitude of the vector is the distance between
         ! origin and the plane.
         ! ===========================================================
      real (kind=RealKind), pointer :: vpsq(:) ! square of cell boundary 
                                               ! plane vectors
!
      real (kind=RealKind), pointer :: corner(:,:)
         ! ===========================================================
         ! corner: polyhedron corner vectors
         ! ===========================================================
      real (kind=RealKind), pointer :: cornsq(:) ! square of polyhedron 
                                                 ! corner vectors
!
      real (kind=RealKind), pointer :: edge(:,:,:)
      real (kind=RealKind), pointer :: edgesq(:)
         ! ===========================================================
         !  edge(:,:,1): vector perpendicular to the edge
         !  edge(:,:,2): vector parallel to the edge
         !  edgesq: square of vector edge(:,:,1)
         ! ===========================================================
!
      integer (kind=IntKind) :: NumCpairs ! num. of corner pairs, each of 
                                          ! which is on the same edge
      integer (kind=IntKind), pointer :: edge_ends(:)  ! gives corner index
                                                       ! for each corner 
                                                       ! pair on the same 
                                                       ! edge
!
      real (kind=RealKind), pointer :: vp_center(:,:)
         ! ===========================================================
         ! vp_center: geometrical center of each boundary plane
         ! ===========================================================
!
      real (kind=RealKind), pointer :: PlaneArea(:) ! area of a boundary
                                                    ! plane
      real (kind=RealKind), pointer :: NeighborDist(:) ! distance from the
                                                       ! neighbor atom that
                                                       ! share the same plane
      real (kind=RealKind) :: Volume           ! volume of the polyhedron
      real (kind=RealKind) :: SurfArea         ! surface area of the polyhedron
      real (kind=RealKind) :: rmt              ! inscribed sphere radius
      real (kind=RealKind) :: rws              ! inscribed sphere radius
      real (kind=RealKind) :: rcirc            ! circumscribed sphere radius
   end type PolyhedronStruct
!
   character (len=50) :: stop_routine
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: NumLocalPolyhedra = 0
   integer (kind=IntKind) :: NumGlobalPolyhedra = 0
!
   integer (kind=IntKind), parameter :: max_planes=100
!
   logical :: Initialized = .false.
   logical :: Generated = .false.
!
   integer (kind=IntKind) :: this
!
   type (PolyhedronStruct), allocatable, target :: Polyhedron(:)
!
   real (kind=RealKind) :: UnitBox(3,3), UnitBoxM(3)
   real (kind=RealKind), parameter :: TOLERANCE = HALF*TEN2m6
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPolyhedra(n,Box,istop,iprint)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: istop
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: i
   real (kind=RealKind), intent(in) :: Box(3,3)
!
   if (n < 1) then
      call ErrorHandler('initPolyhedra','invalid num atoms on my PE',n)
   endif
!
   stop_routine = istop
   print_level = iprint
!
   NumLocalPolyhedra = n
   allocate(Polyhedron(n))
   do i=1,n
      Polyhedron(i)%NumPlanes=0
      Polyhedron(i)%NumCorners=0
      Polyhedron(i)%NumEdges=0
      Polyhedron(i)%NumCpairs=0
   enddo
!
   UnitBox(1:3,1:3)=Box(1:3,1:3)
   UnitBoxM(1) = sqrt( UnitBox(1,1)**2+UnitBox(2,1)**2+UnitBox(3,1)**2 )
   UnitBoxM(2) = sqrt( UnitBox(1,2)**2+UnitBox(2,2)**2+UnitBox(3,2)**2 )
   UnitBoxM(3) = sqrt( UnitBox(1,3)**2+UnitBox(2,3)**2+UnitBox(3,3)**2 )
   if(UnitBoxM(1).lt.TOLERANCE .or. UnitBoxM(2).lt.TOLERANCE .or.           &
      UnitBoxM(3).lt.TOLERANCE) then
      call ErrorHandler('initPolyhedra','invalid unit box vector')
   endif
   Initialized = .true.
   Generated = .false.
!
   end subroutine initPolyhedra
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endPolyhedra()
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: i
!
   do i=1,NumLocalPolyhedra
      deallocate( Polyhedron(i)%edge,                                 &
                  Polyhedron(i)%edgesq,  Polyhedron(i)%edge_ends )
      deallocate( Polyhedron(i)%vplane, Polyhedron(i)%vpsq,           &
                  Polyhedron(i)%vp_center, Polyhedron(i)%PlaneArea,   &
                  Polyhedron(i)%NumPlaneCorners, Polyhedron(i)%id_planes, &
                  Polyhedron(i)%NeighborDist )
      deallocate( Polyhedron(i)%corner, Polyhedron(i)%cornsq,         &
                  Polyhedron(i)%index_corner_plane )
   enddo
   deallocate(Polyhedron)
   Initialized = .false.
   Generated = .false.
!   
   end subroutine endPolyhedra
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genPoly1(LocalIndex,GlobalIndex,NumSeeds,SeedPosition)
!  ===================================================================
   use MathParamModule, only : THREE, FOUR, PI
   implicit none
!
   character (len=13), parameter :: sname='genPolyhedron'
!
   integer (kind=IntKind), intent(in) :: LocalIndex
   integer (kind=IntKind), intent(in) :: GlobalIndex
   integer (kind=IntKind), intent(in) :: NumSeeds
   integer (kind=IntKind) :: ip
!
   real (kind=RealKind), intent(in) :: SeedPosition(3,NumSeeds)
   real (kind=RealKind), allocatable :: SeedRad(:)
!
   if (.not.Initialized) then
      call ErrorHandler('genPolyhedron','PolyhedraModule not initialized')
   else if (LocalIndex < 1 .or. LocalIndex > NumLocalPolyhedra) then
      call ErrorHandler('genPolyhedron','Invalid local seed ID',LocalIndex)
   else if (GlobalIndex < 1 .or. GLobalIndex > NumSeeds) then
      call ErrorHandler('genPolyhedron','Invalid global seed ID',GlobalIndex)
   endif
!
   this = LocalIndex
!
   if(Polyhedron(this)%NumPlanes > 0) then
      deallocate( Polyhedron(this)%vplane, Polyhedron(this)%vpsq,       &
                  Polyhedron(this)%vp_center, Polyhedron(this)%PlaneArea,  &
                  Polyhedron(this)%NumPlaneCorners,                     &
                  Polyhedron(this)%id_planes, Polyhedron(this)%NeighborDist )
      Polyhedron(this)%NumPlanes = 0
   endif
   if(Polyhedron(this)%NumCorners > 0) then
      deallocate( Polyhedron(this)%corner, Polyhedron(this)%cornsq,     &
                  Polyhedron(this)%index_corner_plane )
      Polyhedron(this)%NumCorners = 0
   endif
   if(Polyhedron(this)%NumEdges > 0) then
      deallocate( Polyhedron(this)%edge,                                &
                  Polyhedron(this)%edgesq, Polyhedron(this)%edge_ends )
      Polyhedron(this)%NumEdges = 0
   endif
!
   NumGlobalPolyhedra = NumSeeds
!
!  ===================================================================
!  For each network mesh point get possible VP boundary planes.........
!  ===================================================================
   allocate(SeedRad(NumSeeds))
   SeedRad(1:NumSeeds) = ONE
!  -------------------------------------------------------------------
   call setupCellBound(GlobalIndex,NumSeeds,SeedPosition,SeedRad)
!  -------------------------------------------------------------------
   deallocate(SeedRad)
!
!  ===================================================================
!  For current network mesh point get the corners, edges, etc of 
!  the polyhedron.................................................
!  ------------------------------------------------------------------
   call setupPolyhedron()
!  ------------------------------------------------------------------
!
!  ==================================================================
!  For current network mesh point get the corners, edges, etc of 
!  the polyhedron.................................................
!  ------------------------------------------------------------------
   call geticp()
!  ------------------------------------------------------------------
!
!  ==================================================================
!  For current network mesh point get the area of every boundary 
!  plane of the polyhedron........................................
!  ------------------------------------------------------------------
   call getparea(GlobalIndex)
!  ------------------------------------------------------------------
!
   Polyhedron(this)%rmt   = sqrt(Polyhedron(this)%vpsq(1))
   Polyhedron(this)%rws   = (THREE*Polyhedron(this)%Volume/(FOUR*PI))**THIRD
   Polyhedron(this)%rcirc =   &
      sqrt(Polyhedron(this)%cornsq(Polyhedron(this)%NumCorners))
!
   Polyhedron(this)%SurfArea=ZERO
   do ip=1,Polyhedron(this)%NumPlanes
      Polyhedron(this)%SurfArea = &
         Polyhedron(this)%SurfArea+Polyhedron(this)%PlaneArea(ip)
   enddo
!
   Polyhedron(this)%GlobalIndex = GlobalIndex
!
   Generated = .true.
!
   if (sname==stop_routine) then
      call StopHandler(sname)
   endif
!
   end subroutine genPoly1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genPoly2(LocalIndex,GlobalIndex,NumSeeds,SeedPosition,SeedRad)
!  ===================================================================
   use MathParamModule, only : THREE, FOUR, PI
   implicit none
!
   character (len=13), parameter :: sname='genPolyhedron'
!
   integer (kind=IntKind), intent(in) :: LocalIndex
   integer (kind=IntKind), intent(in) :: GlobalIndex
   integer (kind=IntKind), intent(in) :: NumSeeds
   integer (kind=IntKind) :: ip
   real (kind=RealKind), intent(in) :: SeedPosition(3,NumSeeds)
   real (kind=RealKind), intent(in) :: SeedRad(NumSeeds)
!
   if (.not.Initialized) then
      call ErrorHandler('genPolyhedron','PolyhedraModule not initialized')
   else if (LocalIndex < 1 .or. LocalIndex > NumLocalPolyhedra) then
      call ErrorHandler('genPolyhedron','Invalid local seed ID',LocalIndex)
   else if (GlobalIndex < 1 .or. GLobalIndex > NumSeeds) then
      call ErrorHandler('genPolyhedron','Invalid global seed ID',GlobalIndex)
   endif
!
   this = LocalIndex
!
   if(Polyhedron(this)%NumPlanes > 0) then
      deallocate( Polyhedron(this)%vplane, Polyhedron(this)%vpsq,       &
                  Polyhedron(this)%vp_center, Polyhedron(this)%PlaneArea,  &
                  Polyhedron(this)%NumPlaneCorners,                     &
                  Polyhedron(this)%id_planes, Polyhedron(this)%NeighborDist )
      Polyhedron(this)%NumPlanes = 0
   endif
   if(Polyhedron(this)%NumCorners > 0) then
      deallocate( Polyhedron(this)%corner, Polyhedron(this)%cornsq,     &
                  Polyhedron(this)%index_corner_plane )
      Polyhedron(this)%NumCorners = 0
   endif
   if(Polyhedron(this)%NumEdges > 0) then
      deallocate( Polyhedron(this)%edge,                                &
                  Polyhedron(this)%edgesq, Polyhedron(this)%edge_ends )
      Polyhedron(this)%NumEdges = 0
   endif
!
   NumGlobalPolyhedra = NumSeeds
!
!  ===================================================================
!  For each network mesh point get possible VP boundary planes........
!  -------------------------------------------------------------------
   call setupCellBound(GlobalIndex,NumSeeds,SeedPosition,SeedRad)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  For current network mesh point get the corners, edges, etc of 
!  the polyhedron.................................................
!  ------------------------------------------------------------------
   call setupPolyhedron()
!  ------------------------------------------------------------------
!
!  ==================================================================
!  For current network mesh point get the corners, edges, etc of 
!  the polyhedron.................................................
!  ------------------------------------------------------------------
   call geticp()
!  ------------------------------------------------------------------
!
!  ==================================================================
!  For current network mesh point get the area of every boundary 
!  plane of the polyhedron........................................
!  ------------------------------------------------------------------
   call getparea(GlobalIndex)
!  ------------------------------------------------------------------
!
   Polyhedron(this)%rmt   = sqrt(Polyhedron(this)%vpsq(1))
   Polyhedron(this)%rws   = (THREE*Polyhedron(this)%Volume/(FOUR*PI))**THIRD
   Polyhedron(this)%rcirc =   &
      sqrt(Polyhedron(this)%cornsq(Polyhedron(this)%NumCorners))
!
   Polyhedron(this)%SurfArea=ZERO
   do ip=1,Polyhedron(this)%NumPlanes
      Polyhedron(this)%SurfArea = &
         Polyhedron(this)%SurfArea+Polyhedron(this)%PlaneArea(ip)
   enddo
!
   Polyhedron(this)%GlobalIndex = GlobalIndex
!
   Generated = .true.
!
   if (sname==stop_routine) then
      call StopHandler(sname)
   endif
!
   end subroutine genPoly2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumPolyhedra() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumPolyhedra','PolyhedraModule not initialized')
   endif
   n = NumLocalPolyhedra
   end function getNumPolyhedra
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPolyhedron(i)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: j
   type (PolyhedronStruct), pointer :: p
!
   if (.not.Generated) then
      call ErrorHandler('printPolyhedron','Polyhedron not generated')
   else if (i < 1 .or. i > NumLocalPolyhedra) then
      call ErrorHandler('printPolyhedron','invalid polyhedron index',i)
   endif
!
   p => Polyhedron(i)
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')'*********************************'
   write(6,'( 24x,a )')'*  Output from printPolyhedron  *'
   write(6,'(24x,a,/)')'*********************************'
   write(6,'(a,i5)')   'Atom index:                 ',i
   write(6,'(a,i5)')   'Number of boundary planes:  ',p%NumPlanes
   write(6,'(a,i5)')   'Number of corners:          ',p%NumCorners
   write(6,'(a,i5)')   'Number of edges:            ',p%NumEdges
   write(6,'(a,f10.5)')'Cell volume:                ',p%Volume
   write(6,'(a,f10.5)')'Cell surface area:          ',p%SurfArea
   write(6,'(a,f10.5)')'Inscribed sphere radius:    ',p%rmt
   write(6,'(a,f10.5)')'Wigner-Seitz sphere radius: ',p%rws
   write(6,'(a,f10.5)')'Bounding sphere radius:     ',p%rcirc
!
   write(6,'(/,a)')'Plane Table'
   write(6,'(80(''=''))')
   write(6,'(a)') &
'Plane  Nb Atom  No Corners          Plane Vector            Distance      Area'
   write(6,'(80(''-''))')
   do j=1,p%NumPlanes
      write(6,'(i3,3x,i5,8x,i3,4x,3f10.5,2x,f10.5,2x,f10.5)')    &
            j,p%id_planes(j),p%NumPlaneCorners(j),p%vplane(1,j), &
              p%vplane(2,j),p%vplane(3,j),sqrt(p%vpsq(j)),p%PlaneArea(j)
!     write(6,'(/,a,i5)')  'Boundary plane index:   ',j
!     write(6,'(a,i5)')    'Neighboring atom index: ',p%id_planes(j)
!     write(6,'(a,i5)')    'Number of vertices:     ',p%NumPlaneCorners(j)
!     write(6,'(a,3f10.5)')'Plane vector:           ',p%vplane(1,j),  &
!                          p%vplane(2,j),p%vplane(3,j)
!     write(6,'(a,3f10.5)')'Plane distance:         ',sqrt(p%vpsq(j))
!     write(6,'(a, f10.5)')'Plane area:             ',p%PlaneArea(j)
   enddo
   write(6,'(80(''=''))')
!
   write(6,'(/,a)')'Corner Table'
   write(6,'(80(''=''))')
   write(6,'(a)')'Corner             Corner Vector             Corner Distance'
   write(6,'(80(''-''))')
   do j=1,p%NumCorners
      write(6,'(i4,6x,3f10.5,7x,f10.5)')j,p%corner(1,j),p%corner(2,j), &
                                        p%corner(3,j),sqrt(p%cornsq(j))
!     write(6,'(/,a,i5)')  'Cell corner index:       ',j
!     write(6,'(a,3f10.5)')'Corner vector:           ',p%corner(1,j),&
!                          p%corner(2,j),p%corner(3,j)
!     write(6,'(a, f10.5)')'Corner distance:         ',sqrt(p%cornsq(j))
   enddo
   write(6,'(80(''=''),/)')
   write(6,'(80(''-''))')
!
   end subroutine printPolyhedron
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPolyhedronBoundary(i)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: j, k, n, id
   type (PolyhedronStruct), pointer :: p
!
   if (.not.Generated) then
      call ErrorHandler('printPolyhedronBoundary','Polyhedron not generated')
   else if (i < 1 .or. i > NumLocalPolyhedra) then
      call ErrorHandler('printPolyhedronBoundary','invalid polyhedron index',i)
   endif
!
   p => Polyhedron(i)
!
   write(6,'(/,80(''-''))')
   write(6,'(/,20x,a)')'*****************************************'
   write(6,'( 20x,a )')'*  Output from printPolyhedronBoundary  *'
   write(6,'(20x,a,/)')'*****************************************'
   write(6,'(a,i5)')   'Atom index:                 ',i
   write(6,'(a,i5)')   'Number of boundary planes:  ',p%NumPlanes
!
   do j=1,p%NumPlanes
      n = p%NumPlaneCorners(j)
      write(6,'(/,a,4f10.5)')'Plane: ',                               &
            p%vplane(1,j),p%vplane(2,j),p%vplane(3,j),sqrt(p%vpsq(j))
      write(6,'(a,i5)')'Number of Corners: ',n
      write(6,'(80(''=''))')
      write(6,'(a)')                                                  &
          'Corner             Corner Vector             Corner Distance'
      write(6,'(80(''-''))')
      do k = 1, n
         id = p%index_corner_plane(k,j)
         write(6,'(i4,6x,3f10.5,7x,f10.5)')id,p%corner(1,id),p%corner(2,id), &
                                           p%corner(3,id),sqrt(p%cornsq(id))
      enddo
      write(6,'(80(''=''))')
   enddo
!
   end subroutine printPolyhedronBoundary
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPolyhedraTable(gid,cid,print_pe)
!  ===================================================================
   use MathParamModule, only : PI4, ZERO
!
   use GroupCommModule, only : getMyClusterIndex, getMyPEinGroup
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: gid, cid, print_pe
!
   integer (kind=IntKind) :: i, ig
   integer (kind=IntKind), allocatable :: IntBuffer(:,:)
!
   real (kind=RealKind), allocatable :: RealBuffer(:,:)
   real (kind=RealKind) :: skewness, stretch, surfratio, volmt
!
   type (PolyhedronStruct), pointer :: p
!
   if (.not.Generated) then
      call ErrorHandler('printPolyhedraTable','Polyhedron not generated')
   endif
!
   allocate( IntBuffer(3,NumGlobalPolyhedra), RealBuffer(6,NumGlobalPolyhedra) )
!
   IntBuffer(1:3,1:NumGlobalPolyhedra) = 0
   RealBuffer(1:6,1:NumGlobalPolyhedra) = ZERO
   do i = 1, NumLocalPolyhedra
      p => Polyhedron(i)
      ig = p%GlobalIndex
      IntBuffer(1,ig) = p%NumPlanes
      IntBuffer(2,ig) = p%NumCorners
      IntBuffer(3,ig) = p%NumEdges
      RealBuffer(1,ig) = p%Volume
      RealBuffer(2,ig) = p%SurfArea
      RealBuffer(3,ig) = p%rmt
      RealBuffer(4,ig) = p%rws
      RealBuffer(5,ig) = p%rcirc
      RealBuffer(6,ig) = minval( p%NeighborDist(1:p%NumPlanes) )
   enddo
!
!  -------------------------------------------------------------------
   call GlobalSumInGroup(gid,IntBuffer,3,NumGlobalPolyhedra)
   call GlobalSumInGroup(gid,RealBuffer,6,NumGlobalPolyhedra)
!  -------------------------------------------------------------------
!
   if (getMyClusterIndex(gid) == cid .and. getMyPEinGroup(gid) == print_pe) then
      write(6,'(/,a)')'Table of Polyhedra System'
      write(6,'(124(''=''))')
      write(6,'(2x,a,a)')'Index  NumPlanes  NumCorners  NumEdges      Vol_vp        ', &
     &                   'Vol_mt     (Rc-Rmt)/Rws      Rc/Rnn      SurfRatio        Rnn'
      write(6,'(124(''-''))')
      do ig = 1, NumGlobalPolyhedra
         skewness = (RealBuffer(5,ig)-RealBuffer(3,ig))/RealBuffer(4,ig)
         stretch = RealBuffer(5,ig)/RealBuffer(6,ig)
         surfratio = RealBuffer(2,ig)/(PI4*RealBuffer(4,ig)**2)
         volmt = PI4*THIRD*RealBuffer(3,ig)**3
         write(6,'(i6,4x,i5,2(6x,i5),2x,6(2x,f12.5))')ig,IntBuffer(1:3,ig),RealBuffer(1,ig),volmt, &
                                                         skewness, stretch, surfratio, RealBuffer(6,ig)
      enddo
      write(6,'(124(''=''))')
   endif
!
   deallocate( IntBuffer, RealBuffer )
!
   end subroutine printPolyhedraTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupPolyhedron()
!  ===================================================================
   use SortModule, only : HeapSort
!
!  *******************************************************************
!  Look for those corners and edges of the polyhedron, consists of
!  boundary planes given by vplane_x, vplane_y, vplane_z.
!  *******************************************************************
!
   implicit   none
!
   character (len=15), parameter :: sname='setupPolyhedron'
!
   integer (kind=IntKind), allocatable :: indx(:)
   integer (kind=IntKind) :: NumPlanes, NumCorners, NumEdges, NumCpairs
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: ierr
   integer (kind=IntKind) :: isigma
   integer (kind=IntKind) :: iflag
!
   real (kind=RealKind), allocatable :: tmpx(:)
   real (kind=RealKind), allocatable :: tmpy(:)
   real (kind=RealKind), allocatable :: tmpz(:)
   real (kind=RealKind), allocatable :: tmp2(:)
   real (kind=RealKind) :: xm(3,3)
   real (kind=RealKind) :: xminv(3,3)
   real (kind=RealKind) :: x0
   real (kind=RealKind) :: y0
   real (kind=RealKind) :: z0
   real (kind=RealKind) :: ptest
   real (kind=RealKind) :: pxy
   real (kind=RealKind) :: pyz
   real (kind=RealKind) :: pzx
   real (kind=RealKind) :: eijz
   real (kind=RealKind) :: eijx
   real (kind=RealKind) :: eijy
   real (kind=RealKind) :: pijx
   real (kind=RealKind) :: pijy
   real (kind=RealKind) :: pijz
   real (kind=RealKind) :: e2
   real (kind=RealKind) :: ep(2)
   real (kind=RealKind) :: t
   real (kind=RealKind) :: tol
!
   logical :: corner_iterate = .true.
!
   tol = TOLERANCE
   NumPlanes=Polyhedron(this)%NumPlanes
   if(NumPlanes < 4 .or. NumPlanes > max_planes) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'Incorrect NumPlanes',NumPlanes)
!     ----------------------------------------------------------------
   else if(NumPlanes-2 > 6) then
      n=NumPlanes*(NumPlanes-1)*(NumPlanes-2)/6
   else
      n=NumPlanes*(NumPlanes-1)
   endif
!  -------------------------------------------------------------------
   allocate( tmpx(n), tmpy(n), tmpz(n), tmp2(n), indx(n) )
!  -------------------------------------------------------------------
!
!  ===================================================================
!  The equation of a plane is:
!      vplane_x * x + vplane_y * y + vplane_z * z = vpsq,
!  where:
!      vpsq = vplane_x**2+vplane_y**2+vplane_z**2
!
!  Look for corners by solving the equations:
!      vplane_x(i) * x + vplane_y(i) * y + vplane_z(i) * z = vpsq(i)
!      vplane_x(j) * x + vplane_y(j) * y + vplane_z(j) * z = vpsq(j)
!      vplane_x(k) * x + vplane_y(k) * y + vplane_z(k) * z = vpsq(k)
!  ===================================================================
   corner_iterate = .true.
   do while (corner_iterate) 
      NumCorners=0
      do k=1,NumPlanes-2
         xm(1,1)=Polyhedron(this)%vplane(1,k)/Polyhedron(this)%vpsq(k)
         xm(2,1)=Polyhedron(this)%vplane(2,k)/Polyhedron(this)%vpsq(k)
         xm(3,1)=Polyhedron(this)%vplane(3,k)/Polyhedron(this)%vpsq(k)
         do j=k+1,NumPlanes-1
            xm(1,2)=Polyhedron(this)%vplane(1,j)/Polyhedron(this)%vpsq(j)
            xm(2,2)=Polyhedron(this)%vplane(2,j)/Polyhedron(this)%vpsq(j)
            xm(3,2)=Polyhedron(this)%vplane(3,j)/Polyhedron(this)%vpsq(j)
            do i=j+1,NumPlanes
               xm(1,3)=Polyhedron(this)%vplane(1,i)/Polyhedron(this)%vpsq(i)
               xm(2,3)=Polyhedron(this)%vplane(2,i)/Polyhedron(this)%vpsq(i)
               xm(3,3)=Polyhedron(this)%vplane(3,i)/Polyhedron(this)%vpsq(i)
!              -------------------------------------------------------
               call invm3(xm,xminv,ierr)
!              -------------------------------------------------------
               if(ierr.eq.0) then
                  x0=xminv(1,1)+xminv(2,1)+xminv(3,1)
                  y0=xminv(1,2)+xminv(2,2)+xminv(3,2)
                  z0=xminv(1,3)+xminv(2,3)+xminv(3,3)
!                 ----------------------------------------------------
                  call chkpnt(x0,y0,z0,                               &
                              Polyhedron(this)%vplane,                &
                              Polyhedron(this)%vpsq,                  &
                              NumPlanes,TEN2m10,isigma)
!                 ----------------------------------------------------
                  if(isigma.eq.1) then
                     iflag=0
                     n=1
                     LOOP_n: do n = 1, NumCorners
                        if (abs(x0-tmpx(n))+abs(y0-tmpy(n))+abs(z0-tmpz(n)) &
                           .lt. TOLERANCE) then   ! Eliminate redundant corners.....
                           iflag=1
                           exit LOOP_n
                        endif
                     enddo LOOP_n
                     if (iflag.eq.0) then
                        NumCorners=NumCorners+1
                        tmpx(NumCorners)=x0
                        tmpy(NumCorners)=y0
                        tmpz(NumCorners)=z0
                        tmp2(NumCorners)=x0*x0+y0*y0+z0*z0
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
!
      corner_iterate = .false.
      LOOP_ic: do i=1,NumCorners
         do j=1,NumPlanes
            if (abs(tmpx(i)*Polyhedron(this)%vplane(1,j) +            &
                    tmpy(i)*Polyhedron(this)%vplane(2,j) +            &
                    tmpz(i)*Polyhedron(this)%vplane(3,j) - tmp2(i)) < TOLERANCE) then
               do k = j+1, NumPlanes
                  Polyhedron(this)%vplane(1,k-1) = Polyhedron(this)%vplane(1,k)
                  Polyhedron(this)%vplane(2,k-1) = Polyhedron(this)%vplane(2,k)
                  Polyhedron(this)%vplane(3,k-1) = Polyhedron(this)%vplane(3,k)
                  Polyhedron(this)%vpsq(k-1) = Polyhedron(this)%vpsq(k)
               enddo
               Polyhedron(this)%NumPlanes = Polyhedron(this)%NumPlanes - 1
               NumPlanes = NumPlanes - 1
               corner_iterate = .true.
               exit Loop_ic
            endif
         enddo
      enddo LOOP_ic
      if (NumPlanes < 4) then
!        -------------------------------------------------------------
         call ErrorHandler(sname,'NumPlanes is reduced to less than 4',NumPlanes)
!        -------------------------------------------------------------
      endif
   enddo
!
   if (NumCorners < 4 .or. &
       NumCorners > NumPlanes*(NumPlanes-1)*(NumPlanes-2)/6) then
!     ----------------------------------------------------------------
         call ErrorHandler(sname,'NumCorners out of range',NumCorners)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   allocate( Polyhedron(this)%corner(3,NumCorners),                   &
             Polyhedron(this)%cornsq(NumCorners),                     &
             Polyhedron(this)%index_corner_plane(NumCorners+1,NumPlanes) )
!  -------------------------------------------------------------------
   Polyhedron(this)%NumCorners = NumCorners
!  -------------------------------------------------------------------
   do i=1,NumCorners
      Polyhedron(this)%cornsq(i)=tmp2(i)
   enddo
!  -------------------------------------------------------------------
   call HeapSort(NumCorners,Polyhedron(this)%cornsq,indx)
!  -------------------------------------------------------------------
   do i=1,NumCorners
      j=indx(i)
      Polyhedron(this)%corner(1,i)=tmpx(j)
      Polyhedron(this)%corner(2,i)=tmpy(j)
      Polyhedron(this)%corner(3,i)=tmpz(j)
   enddo
!
!  ===================================================================
!  Look for edges. edges are represented by two vectors and one 
!  parameters as:
!       ->       ->
!       p  + t * e
!        n        n
!  where:
!       ->
!       p  = the vector starting from the origin and endind at and
!        n   perpendicular to the edge     
!       ->
!       e  = a vector parallel to the edge.
!        n
!       t  = a real parameter.
!  ===================================================================
   NumEdges=0
   do j=1,NumPlanes-1
      xm(1,1)=Polyhedron(this)%vplane(1,j)/Polyhedron(this)%vpsq(j)
      xm(2,1)=Polyhedron(this)%vplane(2,j)/Polyhedron(this)%vpsq(j)
      xm(3,1)=Polyhedron(this)%vplane(3,j)/Polyhedron(this)%vpsq(j)
      LOOP_i: do i=j+1,NumPlanes
         ptest = sqrt(Polyhedron(this)%vpsq(i)*Polyhedron(this)%vpsq(j)) &
            -abs(Polyhedron(this)%vplane(1,i)*Polyhedron(this)%vplane(1,j)+  &
                 Polyhedron(this)%vplane(2,i)*Polyhedron(this)%vplane(2,j)+  &
                 Polyhedron(this)%vplane(3,i)*Polyhedron(this)%vplane(3,j))
         if(ptest.gt.TOLERANCE) then
!           ==========================================================
!           look for vector eij first by solving equations:
!
!               ->   ->
!               e  * R = 0
!                ij   i
!
!               ->   ->
!               e  * R = 0
!                ij   j
!           ==========================================================
            pxy=Polyhedron(this)%vplane(1,i)*Polyhedron(this)%vplane(2,j) &
                  -Polyhedron(this)%vplane(1,j)*Polyhedron(this)%vplane(2,i)
            pyz=Polyhedron(this)%vplane(2,i)*Polyhedron(this)%vplane(3,j) &
                  -Polyhedron(this)%vplane(2,j)*Polyhedron(this)%vplane(3,i)
            pzx=Polyhedron(this)%vplane(3,i)*Polyhedron(this)%vplane(1,j) &
                  -Polyhedron(this)%vplane(3,j)*Polyhedron(this)%vplane(1,i)
            if(abs(pxy) .gt. TOLERANCE) then
               eijz=ONE
               eijx=(Polyhedron(this)%vplane(2,i)*Polyhedron(this)%vplane(3,j) &
                  -Polyhedron(this)%vplane(2,j)*Polyhedron(this)%vplane(3,i)) &
                  /pxy
               eijy=(Polyhedron(this)%vplane(1,j)*Polyhedron(this)%vplane(3,i) &
                  -Polyhedron(this)%vplane(1,i)*Polyhedron(this)%vplane(3,j)) &
                  /pxy
            else if(abs(pyz) .gt. TOLERANCE) then
               eijx=ONE
               eijy=(Polyhedron(this)%vplane(3,i)*Polyhedron(this)%vplane(1,j) &
                  -Polyhedron(this)%vplane(3,j)*Polyhedron(this)%vplane(1,i)) &
                  /pyz
               eijz=(Polyhedron(this)%vplane(2,j)*Polyhedron(this)%vplane(1,i) &
                  -Polyhedron(this)%vplane(2,i)*Polyhedron(this)%vplane(1,j)) &
                  /pyz
            else if(abs(pzx) .gt. TOLERANCE) then
               eijy=ONE
               eijz=(Polyhedron(this)%vplane(1,i)*Polyhedron(this)%vplane(2,j) &
                  -Polyhedron(this)%vplane(1,j)*Polyhedron(this)%vplane(2,i)) &
                  /pzx
               eijx=(Polyhedron(this)%vplane(3,i)*Polyhedron(this)%vplane(2,j) &
                  -Polyhedron(this)%vplane(3,j)*Polyhedron(this)%vplane(2,i)) &
                  /pzx
            else
!              =======================================================
!              plane i is parallel to plane j.
!              =======================================================
               cycle LOOP_i
            endif
!           ==========================================================
!           look for vector pij by solving equations:
!
!               ->   ->   2
!               p  * R = R
!                ij   j   j
!
!               ->   ->   2
!               p  * R = R
!                ij   i   i
!
!               ->   ->
!               p  * e  = 0
!                ij   ij
!           ==========================================================
            xm(1,2)=Polyhedron(this)%vplane(1,i)/Polyhedron(this)%vpsq(i)
            xm(2,2)=Polyhedron(this)%vplane(2,i)/Polyhedron(this)%vpsq(i)
            xm(3,2)=Polyhedron(this)%vplane(3,i)/Polyhedron(this)%vpsq(i)
            xm(1,3)=eijx
            xm(2,3)=eijy
            xm(3,3)=eijz
!           ----------------------------------------------------------
            call invm3(xm,xminv,ierr)
!           ----------------------------------------------------------
            if(ierr.ne.0) then
               call ErrorHandler(sname,'Cannot find pij.')
            endif
            pijx=xminv(1,1)+xminv(2,1)
            pijy=xminv(1,2)+xminv(2,2)
            pijz=xminv(1,3)+xminv(2,3)
!           ==========================================================
!           check if there is any portion of the edge inside the 
!           voronoi polyhedron.....................................
!           ----------------------------------------------------------
            call chkedge(pijx,pijy,pijz,eijx,eijy,eijz,                 &
                         Polyhedron(this)%vplane,Polyhedron(this)%vpsq, &
                         NumPlanes,TOLERANCE,isigma)
!           ----------------------------------------------------------
            if(isigma.eq.1) then
               NumEdges=NumEdges+1
               n=2*NumEdges-1
               tmp2(NumEdges)=pijx**2+pijy**2+pijz**2
               tmpx(n)=pijx
               tmpx(n+1)=eijx
               tmpy(n)=pijy
               tmpy(n+1)=eijy
               tmpz(n)=pijz
               tmpz(n+1)=eijz
            endif
         endif
      enddo LOOP_i
   enddo
   if(NumEdges < 6 .or. NumEdges > NumPlanes*(NumPlanes-1)/2) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'NumEdges out of range',NumEdges)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   allocate( Polyhedron(this)%edge(3,NumEdges,2),                     &
             Polyhedron(this)%edgesq(NumEdges),                       &
             Polyhedron(this)%edge_ends(2*NumEdges) )
!  -------------------------------------------------------------------
   Polyhedron(this)%NumEdges = NumEdges
!  -------------------------------------------------------------------
   do i=1,NumEdges
      Polyhedron(this)%edgesq(i)=tmp2(i)
   enddo
!  -------------------------------------------------------------------
   call HeapSort(NumEdges,Polyhedron(this)%edgesq,indx)
!  -------------------------------------------------------------------
   do i=1,NumEdges
      j=2*indx(i)-1
      Polyhedron(this)%edge(1,i,1)=tmpx(j)
      Polyhedron(this)%edge(1,i,2)=tmpx(j+1)
      Polyhedron(this)%edge(2,i,1)=tmpy(j)
      Polyhedron(this)%edge(2,i,2)=tmpy(j+1)
      Polyhedron(this)%edge(3,i,1)=tmpz(j)
      Polyhedron(this)%edge(3,i,2)=tmpz(j+1)
   enddo
!
!  ===================================================================
!  determine corner-pairs, each of which is on the same edge and is
!  stored in edge_ends.................................................
!  Note, edge_ends contains the index of corners.......................
!  ===================================================================
   NumCpairs=0
   do i=1,NumEdges
      e2=Polyhedron(this)%edge(1,i,2)*Polyhedron(this)%edge(1,i,2)+   &
         Polyhedron(this)%edge(2,i,2)*Polyhedron(this)%edge(2,i,2)+   &
         Polyhedron(this)%edge(3,i,2)*Polyhedron(this)%edge(3,i,2)
      iflag=0
      do j=1,NumCorners
         t=( Polyhedron(this)%corner(1,j)*Polyhedron(this)%edge(1,i,2)+      &
             Polyhedron(this)%corner(2,j)*Polyhedron(this)%edge(2,i,2)+      &
             Polyhedron(this)%corner(3,j)*Polyhedron(this)%edge(3,i,2) )/e2
         if( abs(Polyhedron(this)%corner(1,j)-t*Polyhedron(this)%edge(1,i,2)- &
                 Polyhedron(this)%edge(1,i,1)).lt.TOLERANCE .and.  &
             abs(Polyhedron(this)%corner(2,j)-t*Polyhedron(this)%edge(2,i,2)- &
                 Polyhedron(this)%edge(2,i,1)).lt.TOLERANCE .and.  &
             abs(Polyhedron(this)%corner(3,j)-t*Polyhedron(this)%edge(3,i,2)- &
                 Polyhedron(this)%edge(3,i,1)).lt.TOLERANCE ) then
            iflag=iflag+1
            if(iflag.gt.2) then
!              -------------------------------------------------------
               call ErrorHandler(sname,'more than 2 corners on this edge')
!              -------------------------------------------------------
            else
               ep(iflag)=j
            endif
         endif
      enddo
      if(iflag.eq.2) then
         NumCpairs=NumCpairs+1
         n=2*NumCpairs-1
         Polyhedron(this)%edge_ends(n)=ep(1)
         Polyhedron(this)%edge_ends(n+1)=ep(2)
      endif
   enddo
   Polyhedron(this)%NumCpairs = NumCpairs
!
!  -------------------------------------------------------------------
   deallocate( tmpx, tmpy, tmpz, tmp2, indx )
!  -------------------------------------------------------------------
!
!  ===================================================================
!  print out if needed.................................................
!  ===================================================================
   if(print_level.ge.1) then
      write(6,'(/,'' NumPlanes  ='',i5)')NumPlanes
      do i=1,NumPlanes
         write(6,'(''  plane :'',i3,4(2x,e15.8))')                     &
                       i,Polyhedron(this)%vplane(1,i),                 &
                       Polyhedron(this)%vplane(2,i),                   &
                       Polyhedron(this)%vplane(3,i),Polyhedron(this)%vpsq(i)
      enddo
      write(6,'(/,'' NumCorners ='',i5)')NumCorners
      do i=1,NumCorners
         write(6,'(''  corner:'',i3,3(2x,e15.8))')                     &
                       i,Polyhedron(this)%corner(1,i),                 &
                       Polyhedron(this)%corner(2,i),                   &
                       Polyhedron(this)%corner(3,i)
      enddo
      write(6,'(/,'' NumEdges ='',i5)')NumEdges
      do i=1,NumEdges
         write(6,'(''  edge_e: '',i3,3(2x,e15.8))')                    &
                   i,Polyhedron(this)%edge(1,i,2),                  &
                   Polyhedron(this)%edge(2,i,2),Polyhedron(this)%edge(3,i,2)
         write(6,'(''  edge_p: '',i3,3(2x,e15.8))')                    &
                   i,Polyhedron(this)%edge(1,i,1),                  &
                   Polyhedron(this)%edge(2,i,1),Polyhedron(this)%edge(3,i,1)
      enddo
      write(6,'(/,'' NumCpairs ='',i5)')NumCpairs
      do i=1,NumCpairs
         n=2*i-1
         write(6,'(''  edge_ends:  '',3i5)')                           &
                       i,Polyhedron(this)%edge_ends(n),                &
                       Polyhedron(this)%edge_ends(n+1)
      enddo
      if (print_level >= 1) then
         write(6,'(/)')
         do i=1,NumCpairs
            n=2*i-1
            write(6,'(''  set arrow'',i5,''  from '',                     &
                  &   f8.4,",",f8.4,",",f8.4,"  to ",                     &
                  &   f8.4,",",f8.4,",",f8.4,"  nohead")') i,             &
                  Polyhedron(this)%corner(1,Polyhedron(this)%edge_ends(n)), &
                  Polyhedron(this)%corner(2,Polyhedron(this)%edge_ends(n)), &
                  Polyhedron(this)%corner(3,Polyhedron(this)%edge_ends(n)), &
                  Polyhedron(this)%corner(1,Polyhedron(this)%edge_ends(n+1)), &
                  Polyhedron(this)%corner(2,Polyhedron(this)%edge_ends(n+1)), &
                  Polyhedron(this)%corner(3,Polyhedron(this)%edge_ends(n+1))
         enddo
      endif
   endif
!
   end subroutine setupPolyhedron
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupCellBound(i_seed,num_seeds,seed_pos,radical)
!  ===================================================================
   use SortModule, only : HeapSort
   implicit   none
!
   character (len=14), parameter :: sname='setupCellBound'
!
   integer (kind=IntKind), parameter :: nmax=4
   integer (kind=IntKind), parameter :: n0max=(2*nmax+1)**3
!
   integer (kind=IntKind), intent(in) :: i_seed
   integer (kind=IntKind), intent(in) :: num_seeds
!
   integer (kind=IntKind) :: tempid(max_planes)
   integer (kind=IntKind) :: indx(max_planes)
   integer (kind=IntKind) :: j_seed
   integer (kind=IntKind) :: nm1
   integer (kind=IntKind) :: isigma
   integer (kind=IntKind) :: n0
   integer (kind=IntKind) :: in(n0max)
   integer (kind=IntKind) :: jn(n0max)
   integer (kind=IntKind) :: kn(n0max)
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: n2
   integer (kind=IntKind) :: n3
   integer (kind=IntKind) :: i0
   integer (kind=IntKind) :: j0
   integer (kind=IntKind) :: k0
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: nb
   integer (kind=IntKind) :: nbt
   integer (kind=IntKind) :: checked(max_planes)
   integer (kind=IntKind) :: NumPlanes
!
   real (kind=RealKind), intent(in) :: seed_pos(3,num_seeds)
   real (kind=RealKind), intent(in) :: radical(num_seeds)
!
   real (kind=RealKind) :: shift_x
   real (kind=RealKind) :: shift_y
   real (kind=RealKind) :: shift_z
   real (kind=RealKind) :: temp(3,max_planes)
   real (kind=RealKind) :: tmpr2(max_planes)
   real (kind=RealKind) :: tmpd2(max_planes)
   real (kind=RealKind) :: x0
   real (kind=RealKind) :: y0
   real (kind=RealKind) :: z0
   real (kind=RealKind) :: x
   real (kind=RealKind) :: y
   real (kind=RealKind) :: z
   real (kind=RealKind) :: r2
   real (kind=RealKind) :: xt
   real (kind=RealKind) :: yt
   real (kind=RealKind) :: zt
   real (kind=RealKind) :: r2t
   real (kind=RealKind) :: r
   real (kind=RealKind), allocatable :: radsq(:)
   real (kind=RealKind) :: d2, radsq_dif, CutRatio, d2t
!
!  ==================================================================
!  determine i,j,k data set.......................................
!  ==================================================================
   if(num_seeds.ge.16) then
      r=max(UnitBoxM(1),UnitBoxM(2),UnitBoxM(3))
      n1=min(int(r/UnitBoxM(1)+0.8),nmax)
      n2=min(int(r/UnitBoxM(2)+0.8),nmax)
      n3=min(int(r/UnitBoxM(3)+0.8),nmax)
      if((2*n1+1)*(2*n2+1)*(2*n3+1) .gt. n0max) then
         call ErrorHandler(sname,'Needs to increase nmax')
      endif
   else
      n1=nmax
      n2=nmax
      n3=nmax
   endif
   n0=0
   do i0=-n1,n1
      do j0=-n2,n2
         do k0=-n3,n3
            n0=n0+1
            in(n0)=i0
            jn(n0)=j0
            kn(n0)=k0
         enddo
      enddo
   enddo
!
   allocate( radsq(num_seeds) )
   do n=1,num_seeds
      radsq(n) = radical(n)*radical(n)
   enddo
!
   x0=seed_pos(1,i_seed)
   y0=seed_pos(2,i_seed)
   z0=seed_pos(3,i_seed)
   j=0
   do n=1,n0
      shift_x=in(n)*UnitBox(1,1)+jn(n)*UnitBox(1,2)+kn(n)*UnitBox(1,3)-x0
      shift_y=in(n)*UnitBox(2,1)+jn(n)*UnitBox(2,2)+kn(n)*UnitBox(2,3)-y0
      shift_z=in(n)*UnitBox(3,1)+jn(n)*UnitBox(3,2)+kn(n)*UnitBox(3,3)-z0
      LOOP_j_seed: do j_seed=1,num_seeds
         x=seed_pos(1,j_seed)+shift_x
         y=seed_pos(2,j_seed)+shift_y
         z=seed_pos(3,j_seed)+shift_z
         d2=x*x+y*y+z*z
         if (d2 < TOLERANCE) then
            cycle LOOP_j_seed
         endif
         radsq_dif = radsq(i_seed)-radsq(j_seed)
         if (abs(radsq_dif) < TOLERANCE) then
            CutRatio = HALF
         else
            CutRatio = HALF*(ONE+radsq_dif/d2)
         endif
         x = CutRatio*x
         y = CutRatio*y
         z = CutRatio*z
         r2 = CutRatio*CutRatio*d2
         nb = j_seed
         id = j+1
         LOOP_k: do k=1,j
            if(abs(x-temp(1,k)).lt.TOLERANCE .and.                       &
               abs(y-temp(2,k)).lt.TOLERANCE .and.                       &
               abs(z-temp(3,k)).lt.TOLERANCE) then
               cycle LOOP_j_seed
            else if (r2.lt.tmpr2(k)) then
               id=k
               exit LOOP_k
            endif
         enddo LOOP_k
         do k=id,j
            xt=temp(1,k)
            yt=temp(2,k)
            zt=temp(3,k)
            r2t=tmpr2(k)
            nbt=tempid(k)
            d2t=tmpd2(k)
            temp(1,k)=x
            temp(2,k)=y
            temp(3,k)=z
            tmpr2(k)=r2
            tempid(k)=nb
            tmpd2(k)=d2
            x=xt
            y=yt
            z=zt
            r2=r2t
            nb=nbt
            d2=d2t
         enddo
         if(j.lt.max_planes) then
            j=j+1
            temp(1,j)=x
            temp(2,j)=y
            temp(3,j)=z
            tmpr2(j)=r2
            tempid(j)=nb
            tmpd2(j)=d2
         endif
      enddo LOOP_j_seed
   enddo
   nm1=j
   deallocate( radsq )
!
!  ===================================================================
!  reduce nm1 to speed up the process of looking for boundary planes
!  Warning: it may cause problems in some special situtaions.......
!  ===================================================================
   nm1=min(nm1,max_planes)
!
   NumPlanes=0
   k=nm1
   j=1
   id=1
   n=0
   do while(id.le.nm1)
!     ===============================================================
!     check if [temp(1,j),temp(2,j),temp(3,j)] is an possible boundary 
!     plane. It is a preliminary check before chkbnd..............
!     ===============================================================
      isigma=0
      do i=1,n
         if(id.eq.checked(i)) then
            isigma=1
         endif
      enddo
      if(isigma.eq.0) then
!        ------------------------------------------------------------
         call filter_edge(temp(1,j),temp(2,j),temp(3,j),tmpr2(j),j,  &
                          temp,tmpr2,k,isigma,i)
!        ------------------------------------------------------------
         if(isigma.eq.1 .and. i.gt.j) then
            n=n+1
            checked(n)=i+nm1-k
         endif
      endif
      if(isigma.eq.1) then
         NumPlanes=NumPlanes+1
         j=j+1
      else
         k=k-1
         do i=j,k
            temp(1,i)=temp(1,i+1)
            temp(2,i)=temp(2,i+1)
            temp(3,i)=temp(3,i+1)
         enddo
         do i=j,k
            tmpr2(i)=tmpr2(i+1)
         enddo
         do i=j,k
            tempid(i)=tempid(i+1)
         enddo
         do i=j,k
            tmpd2(i)=tmpd2(i+1)
         enddo
      endif
      id=id+1
   enddo
!
   nm1=NumPlanes
   NumPlanes=0
   k=nm1
   j=1
   id=1
   do while(id.le.nm1)
!     ===============================================================
!     check if [temp(1,j),temp(2,j),temp(3,j)] is an actual boundary 
!     plane
!     ---------------------------------------------------------------
      call chkbnd(temp(1,j),temp(2,j),temp(3,j),tmpr2(j),j,          &
                  temp,tmpr2,k,TOLERANCE,isigma)
!     ---------------------------------------------------------------
      if(isigma.eq.1) then
         NumPlanes=NumPlanes+1
         j=j+1
      else
         k=k-1
         do i=j,k
            temp(1,i)=temp(1,i+1)
            temp(2,i)=temp(2,i+1)
            temp(3,i)=temp(3,i+1)
         enddo
         do i=j,k
            tmpr2(i)=tmpr2(i+1)
         enddo
         do i=j,k
            tempid(i)=tempid(i+1)
         enddo
         do i=j,k
            tmpd2(i)=tmpd2(i+1)
         enddo
      endif
      id=id+1
   enddo
! 
   if (NumPlanes < 4 .or. NumPlanes > max_planes) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'NumPlanes out of range',NumPlanes)
!     ----------------------------------------------------------------
   endif
   Polyhedron(this)%NumPlanes = NumPlanes
!  -------------------------------------------------------------------
   allocate( Polyhedron(this)%vplane(3,NumPlanes),     &
             Polyhedron(this)%vpsq(NumPlanes),         &
             Polyhedron(this)%vp_center(3,NumPlanes),  &
             Polyhedron(this)%PlaneArea(NumPlanes),   &
             Polyhedron(this)%NumPlaneCorners(NumPlanes),      &
             Polyhedron(this)%id_planes(NumPlanes),    &
             Polyhedron(this)%NeighborDist(NumPlanes) )
!  -------------------------------------------------------------------
   call HeapSort(NumPlanes,tmpr2,indx)
!  -------------------------------------------------------------------
   do i=1,NumPlanes
      j=indx(i)
      Polyhedron(this)%vplane(1,i)=temp(1,j)
      Polyhedron(this)%vplane(2,i)=temp(2,j)
      Polyhedron(this)%vplane(3,i)=temp(3,j)
   enddo
   do i=1,NumPlanes
      j=indx(i)
      Polyhedron(this)%vpsq(i)=tmpr2(j)
   enddo
   do i=1,NumPlanes
      j=indx(i)
      Polyhedron(this)%id_planes(i)=tempid(j)
   enddo
   do i=1,NumPlanes
      j=indx(i)
      Polyhedron(this)%NeighborDist(i)=sqrt(tmpd2(j))
   enddo
!
   if (print_level.ge.1) then
      write(6,'(/,5x,60(''=''))')
      write(6,'(5x,''BOUNDARY PLANES:'')')
      write(6,'(5x,60(''=''))')
      write(6,'(5x,''Index'','' NbrTyp'',''   R_x    '',                &
            &      ''   R_y    '',''   R_z    '',''   R^2'')')
      write(6,'(5x,60(''-''))')
      write(6,'(5x,2i5,4f10.5)')(j,Polyhedron(this)%id_planes(j),   &
                                   Polyhedron(this)%vplane(1,j),    &
                                   Polyhedron(this)%vplane(2,j),    &
                                   Polyhedron(this)%vplane(3,j),    &
                                   Polyhedron(this)%vpsq(j),j=1,NumPlanes)
   endif
! 
   end subroutine setupCellBound
!  ==================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine filter_edge(x0,y0,z0,r02,i0,p,dp2,nbnd,isigma,i)
!  ===================================================================
   implicit   none
!
   character (len=11), parameter :: sname='filter_edge'
!
   integer (kind=IntKind), intent(in) :: i0
   integer (kind=IntKind), intent(in) :: nbnd
   integer (kind=IntKind), intent(inout) :: isigma
   integer (kind=IntKind), intent(out) :: i
!
   integer (kind=IntKind) :: ip
   integer (kind=IntKind) :: ierr
!
   real (kind=RealKind), intent(in) :: x0
   real (kind=RealKind), intent(in) :: y0
   real (kind=RealKind), intent(in) :: z0
   real (kind=RealKind), intent(in) :: r02
   real (kind=RealKind), intent(in) :: p(3,nbnd)
   real (kind=RealKind), intent(in) :: dp2(nbnd)
!
   real (kind=RealKind) :: xm(9)
   real (kind=RealKind) :: xminv(9)
   real (kind=RealKind) :: ptest
   real (kind=RealKind) :: pxy
   real (kind=RealKind) :: pyz
   real (kind=RealKind) :: pzx
   real (kind=RealKind) :: eijz
   real (kind=RealKind) :: eijx
   real (kind=RealKind) :: eijy
   real (kind=RealKind) :: pijx
   real (kind=RealKind) :: pijy
   real (kind=RealKind) :: pijz
!
!  *******************************************************************
!  check if the plane defined by (x0,y0,z0) is an possible boundary
!  plane...........................................................
!
!  The idea is to determine the edge formed by this plane and any
!  other one is partially inside (on the boundary of) the polyhedron...
!
!  The return values are isigma and i: if isigma=1, plane i and plane j
!  form an edge, part of which is inside (on the boundary of) the polyhedron
!  *******************************************************************
!
!  ===================================================================
!  Look for edges. edges are represented by two vectors and one 
!  parameters as:
!          ->       ->
!          p  + t * e
!           n        n
!  where:
!          ->
!          p  = the vector starting from the origin and endind at and
!           n   perpendicular to the edge     
!          ->
!          e  = a vector parallel to the edge.
!           n
!          t  = a real parameter.
!  ===================================================================
!
   xm(1)=x0/r02
   xm(2)=y0/r02
   xm(3)=z0/r02
   LOOP_ip: do ip=1,nbnd-1
      i=mod(ip+i0-1,nbnd)+1
      ptest = sqrt(dp2(i)*r02)-abs(p(1,i)*x0+p(2,i)*y0+p(3,i)*z0)
      if (ptest.gt.TOLERANCE) then
!        =============================================================
!        look for vector eij first by solving equations:
!
!              ->   ->
!              e  * R = 0
!               i    i
!
!              ->   ->
!              e  * R = 0
!               i    0
!        =============================================================
         pxy=p(1,i)*y0-x0*p(2,i)
         pyz=p(2,i)*z0-y0*p(3,i)
         pzx=p(3,i)*x0-z0*p(1,i)
         if(abs(pxy) .gt. TOLERANCE) then
            eijz=one
            eijx=(p(2,i)*z0-y0*p(3,i))/pxy
            eijy=(x0*p(3,i)-p(1,i)*z0)/pxy
         else if(abs(pyz) .gt. TOLERANCE) then
            eijx=one
            eijy=(p(3,i)*x0-z0*p(1,i))/pyz
            eijz=(y0*p(1,i)-p(2,i)*x0)/pyz
         else if(abs(pzx) .gt. TOLERANCE) then
            eijy=one
            eijz=(p(1,i)*y0-x0*p(2,i))/pzx
            eijx=(p(3,i)*y0-z0*p(2,i))/pzx
         else
!           ==========================================================
!           plane i is parallel to plane j.
!           ==========================================================
            cycle LOOP_ip
         endif
!        =============================================================
!        look for vector pij by solving equations:
!
!              ->   ->   2
!              p  * R = R
!               ij   j   j
!
!              ->   ->   2
!              p  * R = R
!               ij   i   i
!
!              ->   ->
!              p  * e  = 0
!               ij   ij
!        =============================================================
         xm(4)=p(1,i)/dp2(i)
         xm(5)=p(2,i)/dp2(i)
         xm(6)=p(3,i)/dp2(i)
         xm(7)=eijx
         xm(8)=eijy
         xm(9)=eijz
!        -------------------------------------------------------------
         call invm3(xm,xminv,ierr)
!        -------------------------------------------------------------
         if(ierr.ne.0) then
            call ErrorHandler(sname,'Cannot find pij')
         endif
         pijx=xminv(1)+xminv(2)
         pijy=xminv(4)+xminv(5)
         pijz=xminv(7)+xminv(8)
!        =============================================================
!        check if there is any portion of the edge inside the 
!        voronoi polyhedron.....................................
!        -------------------------------------------------------------
         call chkedge(pijx,pijy,pijz,eijx,eijy,eijz,p,dp2,nbnd,TOLERANCE,isigma)
!        -------------------------------------------------------------
         if(isigma.eq.1) return
      endif
   enddo LOOP_ip
!
   if (stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
!
   end subroutine filter_edge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine chkbnd(x0,y0,z0,r02,i0,p,dp2,nbnd,tol,isigma)
!  ===================================================================
   implicit   none
!
   character (len=6), parameter :: sname='chkbnd'
!
   integer (kind=IntKind), intent(in) :: i0
   integer (kind=IntKind), intent(in) :: nbnd
   integer (kind=IntKind), intent(out) :: isigma
!
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jp
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: kp
   integer (kind=IntKind) :: njoint
   integer (kind=IntKind) :: ierr
   integer (kind=IntKind) :: nbm1
   integer (kind=IntKind) :: nbm2
   integer (kind=IntKind) :: iflag
!
   real (kind=RealKind), intent(in) :: x0
   real (kind=RealKind), intent(in) :: y0
   real (kind=RealKind), intent(in) :: z0
   real (kind=RealKind), intent(in) :: r02
   real (kind=RealKind), intent(in) :: p(3,nbnd)
   real (kind=RealKind), intent(in) :: dp2(nbnd)
   real (kind=RealKind), intent(in) :: tol
!
   real (kind=RealKind) :: x1, y1, z1
   real (kind=RealKind) :: xm(9)
   real (kind=RealKind) :: xminv(9)
   real (kind=RealKind) :: xc(3),yc(3),zc(3)
!
!  tol = TOLERANCE
!  *******************************************************************
!  check if the plane defined by (x0,y0,z0) is an actual boundary
!  plane...........................................................
!  *******************************************************************
!
!  ===================================================================
!  The idea is to find the number of joint points with any other
!  two planes not being outside the polyhedron. If it is more than 
!  two, the plane is a boundary plane............................
!  ===================================================================
   xm(1)=x0/r02
   xm(2)=y0/r02
   xm(3)=z0/r02
   njoint=0
   nbm2=nbnd-2
   nbm1=nbnd-1
   kp=0
   do while(njoint.le.2 .and. kp.lt.nbm2)
      kp=kp+1
      k=mod(kp+i0-1,nbnd)+1
      xm(4)=p(1,k)/dp2(k)
      xm(5)=p(2,k)/dp2(k)
      xm(6)=p(3,k)/dp2(k)
      jp=kp
      do while(njoint.le.2 .and. jp.lt.nbm1)
         jp=jp+1
         j=mod(jp+i0-1,nbnd)+1
         xm(7)=p(1,j)/dp2(j)
         xm(8)=p(2,j)/dp2(j)
         xm(9)=p(3,j)/dp2(j)
!        -------------------------------------------------------------
         call invm3(xm,xminv,ierr)
!        -------------------------------------------------------------
         if(ierr.eq.0) then
            x1=xminv(1)+xminv(2)+xminv(3)
            y1=xminv(4)+xminv(5)+xminv(6)
            z1=xminv(7)+xminv(8)+xminv(9)
!           ----------------------------------------------------------
            call chkpnt(x1,y1,z1,p,dp2,nbnd,TEN2m10,isigma)
!           ----------------------------------------------------------
            if(isigma.eq.1) then
               iflag=0
               LOOP_i: do i=1,njoint
                  if(abs(x1-xc(i))+abs(y1-yc(i))+abs(z1-zc(i)).lt.TEN2m10) then
                     iflag=1
                     exit LOOP_i
                  end if
               end do LOOP_i
               if(iflag.eq.0) then
                  njoint=njoint+1
                  xc(njoint)=x1
                  yc(njoint)=y1
                  zc(njoint)=z1
               endif
            end if
         endif
      enddo
   enddo
!
   if(njoint.ge.3) then
      isigma=1
   else
      isigma=0
   endif
!
   if(stop_routine.eq.sname) then
      call StopHandler(sname)
   end if
!
   end subroutine chkbnd
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine geticp()
!  ===================================================================
!
!  *******************************************************************
!  determine corner index array index_corner_plane, a pointer to polyhedron 
!  corners. Its second dimension corresponds to the boundary plane.
!  Its first dimension corresponds to the corners on the same plane
!  *******************************************************************
!
   implicit none
!
   character (len=6), parameter :: sname='geticp'
!
   integer (kind=IntKind) :: ib
   integer (kind=IntKind) :: ie
   integer (kind=IntKind) :: i1
   integer (kind=IntKind) :: i2
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jp
   integer (kind=IntKind) :: nc
   integer (kind=IntKind) :: np
   integer (kind=IntKind) :: jflag, ic
   integer (kind=IntKind), allocatable :: p1(:)
   integer (kind=IntKind), allocatable :: p2(:)
!
   logical :: error_occured = .false.
!
!  -------------------------------------------------------------------
   allocate( p1(polyhedron(this)%NumPlanes), p2(polyhedron(this)%NumPlanes) )
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Loop over boundary planes and determine index_corner_plane for each 
!  one of them
!  ===================================================================
   LOOP_ib: do ib=1,Polyhedron(this)%NumPlanes
!     ================================================================
!     Looking for paired corners belonging to the plane of index ib
!     ================================================================
      np=0
      do ie=1,Polyhedron(this)%NumCpairs
         j=2*ie-1
         i1=Polyhedron(this)%edge_ends(j)
         i2=Polyhedron(this)%edge_ends(j+1)
         if (abs(Polyhedron(this)%corner(1,i1)*Polyhedron(this)%vplane(1,ib)+ &
                 Polyhedron(this)%corner(2,i1)*Polyhedron(this)%vplane(2,ib)+ &
                 Polyhedron(this)%corner(3,i1)*Polyhedron(this)%vplane(3,ib)- &
                 Polyhedron(this)%vpsq(ib)).lt.TOLERANCE .and.                   &
             abs(Polyhedron(this)%corner(1,i2)*Polyhedron(this)%vplane(1,ib)+ &
                 Polyhedron(this)%corner(2,i2)*Polyhedron(this)%vplane(2,ib)+ &
                 Polyhedron(this)%corner(3,i2)*Polyhedron(this)%vplane(3,ib)- &
                 Polyhedron(this)%vpsq(ib)).lt.TOLERANCE) then
            np=np+1
            p1(np)=i1
            p2(np)=i2
         endif
      enddo
!
      if(np.lt.3) then
         write(6,'(/,a)')'=========================='
         write(6,'(a,2i5)')'Polyhedron Index, Plane Index: ',this, ib
         write(6,'(a,i5)')'Incorrect np value: ',np
         error_occured = .true.
         exit LOOP_ib
      endif
!
      nc=1
      Polyhedron(this)%index_corner_plane(1,ib)=p1(1)
      nc=2
      Polyhedron(this)%index_corner_plane(2,ib)=p2(1)
      do while(nc.le.np)
         jflag = 0
         do jp=2,np
            if(p1(jp).eq.Polyhedron(this)%index_corner_plane(nc,ib) .and.     &
               p2(jp).ne.Polyhedron(this)%index_corner_plane(nc-1,ib)) then
               nc=nc+1
               Polyhedron(this)%index_corner_plane(nc,ib)=p2(jp)
               jflag = 1
            else if(p2(jp).eq.Polyhedron(this)%index_corner_plane(nc,ib) .and.  &
                    p1(jp).ne.Polyhedron(this)%index_corner_plane(nc-1,ib)) then
               nc=nc+1
               Polyhedron(this)%index_corner_plane(nc,ib)=p1(jp)
               jflag = 1
            endif
         enddo
         if (jflag == 0) then
            write(6,'(/,a)')'=========================='
            write(6,'(a,3i5)')'Polyhedron, Plane, Ill Corner Index: ',this, ib, &
                              Polyhedron(this)%index_corner_plane(nc,ib)
            write(6,'(a)')'Index,   Corner Pair'
            do jp=1,np
               write(6,'(i5,4x,2i5)')jp, p1(jp), p2(jp)
            enddo
            write(6,'(a)')'Ill condition: a Corner has less than three Edges attached'
            error_occured = .true.
            exit LOOP_ib
         endif
      enddo
      Polyhedron(this)%NumPlaneCorners(ib)=nc-1
      if (Polyhedron(this)%index_corner_plane(nc,ib) .ne.                 &
          Polyhedron(this)%index_corner_plane(1,ib)) then
         call ErrorHandler(sname,'Inconsistent index_corner_plane',       &
                           Polyhedron(this)%index_corner_plane(nc,ib),   &
                           Polyhedron(this)%index_corner_plane(1,ib))
      endif
   enddo LOOP_ib
!
   if (error_occured) then
      write(6,'(a)')'=========================='
      write(6,'(a)')'Index,   Corner Coordinates'
      do ic = 1, Polyhedron(this)%NumCorners
         write(6,'(i5,4x,3f14.8)') ic, Polyhedron(this)%corner(1,ic), &
                                       Polyhedron(this)%corner(2,ic), &
                                       Polyhedron(this)%corner(3,ic)
      enddo
      write(6,'(a)')'=========================='
      write(6,'(a)')'Index,   Plane  Coordinates'
      do ib = 1, Polyhedron(this)%NumPlanes
         write(6,'(i5,4x,3f14.8)') ib, Polyhedron(this)%vplane(1,ib), &
                                       Polyhedron(this)%vplane(2,ib), &
                                       Polyhedron(this)%vplane(3,ib)
      enddo
      call ErrorHandler(sname,'Error flag raised')
   endif
!
!  ==================================================================
!  find the center of each V-Plane & total number of primitive triangles
!  ==================================================================
   do ib=1,Polyhedron(this)%NumPlanes
      Polyhedron(this)%vp_center(1,ib)=ZERO
      Polyhedron(this)%vp_center(2,ib)=ZERO
      Polyhedron(this)%vp_center(3,ib)=ZERO
      do nc=1,Polyhedron(this)%NumPlaneCorners(ib)
         Polyhedron(this)%vp_center(1,ib)=Polyhedron(this)%vp_center(1,ib)+ &
         Polyhedron(this)%corner(1,Polyhedron(this)%index_corner_plane(nc,ib))
         Polyhedron(this)%vp_center(2,ib)=Polyhedron(this)%vp_center(2,ib)+ &
         Polyhedron(this)%corner(2,Polyhedron(this)%index_corner_plane(nc,ib))
         Polyhedron(this)%vp_center(3,ib)=Polyhedron(this)%vp_center(3,ib)+ &
         Polyhedron(this)%corner(3,Polyhedron(this)%index_corner_plane(nc,ib))
      enddo
      Polyhedron(this)%vp_center(1,ib)=Polyhedron(this)%vp_center(1,ib)     &
                           /dble(Polyhedron(this)%NumPlaneCorners(ib))
      Polyhedron(this)%vp_center(2,ib)=Polyhedron(this)%vp_center(2,ib)     &
                           /dble(Polyhedron(this)%NumPlaneCorners(ib))
      Polyhedron(this)%vp_center(3,ib)=Polyhedron(this)%vp_center(3,ib)     &
                           /dble(Polyhedron(this)%NumPlaneCorners(ib))
   enddo
!
!  --------------------------------------------------------------------
   deallocate( p1, p2 )
!  --------------------------------------------------------------------
!
!  ====================================================================
   if(print_level.ge.1) then
      do ib=1,Polyhedron(this)%NumPlanes
         write(6,'(/,5x,60(''=''))')
         write(6,'(5x,''Plane Index   ='',i5,5x,''Number of Corners ='',&
               &   i5)') ib,Polyhedron(this)%NumPlaneCorners(ib)
         write(6,'(5x,60(''=''))')
         write(6,'(5x,''Corner Number ='',i5,5x,''Corner Index      ='',&
               &   i5 )')(jp,Polyhedron(this)%index_corner_plane(jp,ib), &
                             jp=1,Polyhedron(this)%NumPlaneCorners(ib))
         write(6,'(5x,60(''-''))')
         write(6,'(5x,''Coordinates of Center ='',3f10.5)')             &
                   Polyhedron(this)%vp_center(1,ib),                    &
                   Polyhedron(this)%vp_center(2,ib),                    &
                   Polyhedron(this)%vp_center(3,ib)
         write(6,'(5x,60(''-''))')
      enddo
      write(6,'(/,5x,60(''=''))')
      write(6,'(5x,60(''=''))')
   endif
!
   if (stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
!
   end subroutine geticp
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getparea(i_seed)
!  ===================================================================
!
!  *******************************************************************
!  calculate the boundary plane area.
!  *******************************************************************
!
   implicit none
!
   character (len=8), parameter :: sname='getparea'
!
   integer (kind=IntKind), intent(in) :: i_seed
   integer (kind=IntKind) :: ib
   integer (kind=IntKind) :: ie
   integer (kind=IntKind) :: i1
   integer (kind=IntKind) :: i2
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: ip
   integer (kind=IntKind) :: jp
   integer (kind=IntKind) :: np
   integer (kind=IntKind), allocatable :: p1(:)
   integer (kind=IntKind), allocatable :: p2(:)
!
   real (kind=RealKind) :: psum(3)
   real (kind=RealKind) :: dp
   real (kind=RealKind) :: v1(3), v2(3), v(3)
   real (kind=RealKind) :: px
   real (kind=RealKind) :: py
   real (kind=RealKind) :: pz
!
   j=max(Polyhedron(this)%NumPlanes, Polyhedron(this)%NumEdges)
!  -------------------------------------------------------------------
   allocate( p1(j), p2(j) )
!  -------------------------------------------------------------------
!
!  ===================================================================
!  Loop over boundary plane........................................
!  ===================================================================
   Polyhedron(this)%Volume=ZERO
   psum(1:3)=ZERO
   do ib=1,Polyhedron(this)%NumPlanes
!     ================================================================
!     Looking for paired corners belonging to the plane............
!     ================================================================
      np=0
      do ie=1,Polyhedron(this)%NumCpairs
         j=2*ie-1
         i1=Polyhedron(this)%edge_ends(j)
         i2=Polyhedron(this)%edge_ends(j+1)
         if(abs(Polyhedron(this)%corner(1,i1)*Polyhedron(this)%vplane(1,ib)+ &
                Polyhedron(this)%corner(2,i1)*Polyhedron(this)%vplane(2,ib)+ &
                Polyhedron(this)%corner(3,i1)*Polyhedron(this)%vplane(3,ib)- &
                Polyhedron(this)%vpsq(ib)).lt.TOLERANCE .and.    &
            abs(Polyhedron(this)%corner(1,i2)*Polyhedron(this)%vplane(1,ib)+ &
                Polyhedron(this)%corner(2,i2)*Polyhedron(this)%vplane(2,ib)+ &
                Polyhedron(this)%corner(3,i2)*Polyhedron(this)%vplane(3,ib)- &
                Polyhedron(this)%vpsq(ib)).lt.TOLERANCE) then
            np=np+1
            p1(np)=i1
            p2(np)=i2
         endif
      enddo
!
      px=ZERO
      py=ZERO
      pz=ZERO
      i1=p1(1)
      i2=p2(1)
      do jp=0,np-1
         ip=jp+1
         do while(ip.gt.1 .and. ip.le.np) 
            if(p1(ip).eq.i2)then
               i1=p1(ip)
               i2=p2(ip)
               p1(ip)=p1(jp+1)
               p2(ip)=p2(jp+1)
               ip=np+1
            else if(p2(ip).eq.i2)then
               i1=p2(ip)
               i2=p1(ip)
               p1(ip)=p1(jp+1)
               p2(ip)=p2(jp+1)
               ip=np+1
            else
               ip=ip+1
            endif
         enddo
!        =============================================================
!        Calculate:
!
!           ->   ->
!           c  x c
!            i1   i2
!
!        and store the result in x, y, z...........................
!        =============================================================
         v1(1)=Polyhedron(this)%corner(1,i1)
         v1(2)=Polyhedron(this)%corner(2,i1)
         v1(3)=Polyhedron(this)%corner(3,i1)
         v2(1)=Polyhedron(this)%corner(1,i2)
         v2(2)=Polyhedron(this)%corner(2,i2)
         v2(3)=Polyhedron(this)%corner(3,i2)
!        -------------------------------------------------------------
         call vcross(v1,v2,v)
!        -------------------------------------------------------------
         px=px+v(1)
         py=py+v(2)
         pz=pz+v(3)
      enddo
!     ================================================================
!     Calculate the plane area using formular:
!
!                   1     ->   ->     ->   ->           ->   ->
!           Area = --- | (c1 x c2) + (c2 x c3) + ... + (cn x c1) |
!                   2 
!
!     ================================================================
      Polyhedron(this)%PlaneArea(ib)=HALF*sqrt(px*px+py*py+pz*pz)
      dp=sqrt(Polyhedron(this)%vpsq(ib))
      Polyhedron(this)%Volume=Polyhedron(this)%Volume+               &
                              THIRD*Polyhedron(this)%PlaneArea(ib)*dp
      psum(1)=psum(1)+Polyhedron(this)%PlaneArea(ib)*                &
                      Polyhedron(this)%vplane(1,ib)/dp
      psum(2)=psum(2)+Polyhedron(this)%PlaneArea(ib)*                &
                      Polyhedron(this)%vplane(2,ib)/dp
      psum(3)=psum(3)+Polyhedron(this)%PlaneArea(ib)*                &
                      Polyhedron(this)%vplane(3,ib)/dp
   enddo
!
!  ===================================================================
!  If correct, psum = 0 for each polyhedron........................
!  ===================================================================
   if( abs(psum(1)).ge.TOLERANCE .or. abs(psum(2)).ge.TOLERANCE .or.      &
       abs(psum(3)).ge.TOLERANCE ) then
      write(6,'('' SEED ID., PSUM:'',1i5,3f13.8)')                    &
            i_seed,psum(1),psum(2),psum(3)
      do ib = 1, Polyhedron(this)%NumPlanes
         write(6,'('' BOUNDARY PLANE:'',1i5,3f13.8)')                 &
                   ib,Polyhedron(this)%vplane(1,ib),                  &
                      Polyhedron(this)%vplane(2,ib),                  &
                      Polyhedron(this)%vplane(3,ib)
      enddo
      do ib = 1, Polyhedron(this)%NumCorners
         write(6,'(''         CORNER:'',1i5,3f13.8)')                 &
                   ib,Polyhedron(this)%corner(1,ib),                  &
                      Polyhedron(this)%corner(2,ib),                  &
                      Polyhedron(this)%corner(3,ib)
      enddo
      do ib = 1, Polyhedron(this)%NumEdges
         write(6,'(''          EDGES:'',1i5,3f13.8)')                 &
                   ib,Polyhedron(this)%edge(1,ib,2),                  &
                      Polyhedron(this)%edge(2,ib,2),                  &
                      Polyhedron(this)%edge(3,ib,2)
         write(6,'(''          EDGPS:'',1i5,3f13.8)')                 &
                   ib,Polyhedron(this)%edge(1,ib,1),                  &
                      Polyhedron(this)%edge(2,ib,1),                  &
                      Polyhedron(this)%edge(3,ib,1)
      enddo
      do ib = 1, Polyhedron(this)%NumCpairs
         write(6,'(''           EEND:'',3i5)')                        &
                   ib,Polyhedron(this)%edge_ends(2*ib-1),             &
                      Polyhedron(this)%edge_ends(2*ib)
      enddo
      call ErrorHandler(sname,'sum of boundary plane area vectors is not 0')
   endif
!  -------------------------------------------------------------------
   deallocate( p1, p2 )
!  -------------------------------------------------------------------
!
   if (stop_routine.eq.sname) then
      call StopHandler(sname)
   endif
   end subroutine getparea
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine chkpnt(x0,y0,z0,p,dp2,nbnd,tol,isigma)
!  ===================================================================
   implicit   none
!
   character (len=6), parameter :: sname='chkpnt'
!
   integer (kind=IntKind), intent(in) :: nbnd
   integer (kind=IntKind), intent(out) :: isigma
!
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: dp2(nbnd)
   real (kind=RealKind), intent(in) :: x0
   real (kind=RealKind), intent(in) :: y0
   real (kind=RealKind), intent(in) :: z0
   real (kind=RealKind), intent(in) :: p(3,nbnd)
   real (kind=RealKind), intent(in) :: tol
!
! ********************************************************************
!  check if (x0,y0,z0) is inside the polyhedron bounded by planes:
!  {p(1,ip),p(2,ip),p(3,ip), ip=1,nbnd}
!  return isigma = 1, inside
!                = 0, outside
! ********************************************************************
!
   isigma=1
   do i=1,nbnd
      if ( (x0*p(1,i)+y0*p(2,i)+z0*p(3,i))/dp2(i)-ONE > tol ) then
         isigma=0
         return
      end if
   end do
!
   if (stop_routine.eq.sname) then
      call StopHandler(sname)
   end if
!
   end subroutine chkpnt
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine chkedge(pijx,pijy,pijz,eijx,eijy,eijz,p,dp2,nbnd,tol,isigma)
!  ===================================================================
   implicit   none
!
   character (len=7), parameter ::  sname='chkedge'
!
   integer (kind=IntKind), intent(in) :: nbnd
   integer (kind=IntKind), intent(out) :: isigma
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: pijx
   real (kind=RealKind), intent(in) :: pijy
   real (kind=RealKind), intent(in) :: pijz
   real (kind=RealKind), intent(in) :: eijx
   real (kind=RealKind), intent(in) :: eijy
   real (kind=RealKind), intent(in) :: eijz
   real (kind=RealKind), intent(in) :: p(3,nbnd)
   real (kind=RealKind), intent(in) :: dp2(nbnd)
   real (kind=RealKind), intent(in) :: tol
!
   real (kind=RealKind) :: p2mpr
   real (kind=RealKind) :: pe
   real (kind=RealKind) :: t
   real (kind=RealKind) :: tlt
   real (kind=RealKind) :: tgt
!
!  *******************************************************************
!  check if the edge is outside the polyhedron, returns isigma:
!
!  isigma = 0, if outside
!         = 1, otherwise
!  *******************************************************************
!
   isigma = 1
   LOOP_i: do i=1,nbnd
      p2mpr=(p(1,i)*pijx+p(2,i)*pijy+p(3,i)*pijz)/dp2(i)-ONE
      pe=p(1,i)*eijx+p(2,i)*eijy+p(3,i)*eijz
      if(abs(pe).lt.tol .and. p2mpr > tol) then
         isigma=0
         exit LOOP_i
      endif
   enddo LOOP_i
   return
!
   if(sname.eq.stop_routine) then
      call StopHandler(sname)
   endif
!
!  *******************************************************************
!  The above piece of code is added to replace the following piece, 
!  which seems problematic in some cases.
!  - Yang Wang, Dec. 1st, 2009
!  *******************************************************************
   tlt= 1.0d+30
   tgt=-1.0d+30
   do i=1,nbnd
      p2mpr=dp2(i)-p(1,i)*pijx-p(2,i)*pijy-p(3,i)*pijz
      pe=p(1,i)*eijx+p(2,i)*eijy+p(3,i)*eijz
      if(abs(pe).lt.TOLERANCE .and. p2mpr.lt.-TOLERANCE) then
         isigma=0
         return
      else if(abs(pe).gt.TOLERANCE) then
         t=p2mpr/pe
         if(pe.gt.zero .and. t.lt.tlt) then
            tlt=t
         else if(pe.lt.zero .and. t.gt.tgt) then
            tgt=t
         endif
      endif
   enddo
!
   if(tlt-tgt .lt. TOLERANCE) then
      isigma=0
   else
      isigma=1
   endif
!
   end subroutine chkedge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSurfaceArea(poly) result(s)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
!
   real (kind=RealKind) :: s
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getSurfaceArea','Invalid polyhedron index',poly)
   endif
!
   s = Polyhedron(poly)%SurfArea
!
   end function getSurfaceArea
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumPlanes(poly) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: n
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getNumPlanes','Invalid polyhedron index',poly)
   endif
!
   n = Polyhedron(poly)%NumPlanes
!
   end function getNumPlanes
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumPlaneCorners(poly,plane) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: plane
   integer (kind=IntKind) :: n
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getNumPlaneCorners','Invalid polyhedron index',poly)
   else if (plane < 1 .or. plane > Polyhedron(poly)%NumPlanes) then
      call ErrorHandler('getNumPlaneCorners','Invalid plane index',plane)
   endif
!
   n = Polyhedron(poly)%NumPlaneCorners(plane)
!
   end function getNumPlaneCorners
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPlaneCorner(poly,plane,cn) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: plane
   integer (kind=IntKind), intent(in) :: cn
   integer (kind=IntKind) :: id
!
   real (kind=RealKind) :: c(3)
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getPlaneCorners','Invalid polyhedron index',poly)
   else if (plane < 1 .or. plane > Polyhedron(poly)%NumPlanes) then
      call ErrorHandler('getPlaneCorners','Invalid plane index',plane)
   else if (cn < 1 .or. cn > Polyhedron(poly)%NumPlaneCorners(plane)) then
      call ErrorHandler('getPlaneCorners','Invalid plane corner index',cn)
   endif
!
   id = Polyhedron(poly)%index_corner_plane(cn,plane)
   c(1) = Polyhedron(poly)%corner(1,id)
   c(2) = Polyhedron(poly)%corner(2,id)
   c(3) = Polyhedron(poly)%corner(3,id)
!
   end function getPlaneCorner
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIndexCornerPlane(poly,plane,cn) result(id)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: plane
   integer (kind=IntKind), intent(in) :: cn
   integer (kind=IntKind) :: id
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getPlaneCorners','Invalid polyhedron index',poly)
   else if (plane < 1 .or. plane > Polyhedron(poly)%NumPlanes) then
      call ErrorHandler('getPlaneCorners','Invalid plane index',plane)
   else if (cn < 1 .or. cn > Polyhedron(poly)%NumPlaneCorners(plane)) then
      call ErrorHandler('getPlaneCorners','Invalid plane corner index',cn)
   endif
!
   id = Polyhedron(poly)%index_corner_plane(cn,plane)
!
   end function getIndexCornerPlane
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPlaneArea(poly,plane) result(s)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: plane
!
   real (kind=RealKind) :: s
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getPlaneArea','Invalid polyhedron index',poly)
   else if (plane < 1 .or. plane > Polyhedron(poly)%NumPlanes) then
      call ErrorHandler('getPlaneArea','Invalid plane index',plane)
   endif
!
   s = Polyhedron(poly)%PlaneArea(plane)
   end function getPlaneArea
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPlaneNormal(poly,plane) result(s)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: plane
!
   real (kind=RealKind) :: s(3)
   real (kind=RealKind) :: x, y, z, p
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getPlaneNormal','Invalid polyhedron index',poly)
   else if (plane < 1 .or. plane > Polyhedron(poly)%NumPlanes) then
      call ErrorHandler('getPlaneNormal','Invalid plane index',plane)
   endif
!
   x = Polyhedron(poly)%vplane(1,plane)
   y = Polyhedron(poly)%vplane(2,plane)
   z = Polyhedron(poly)%vplane(3,plane)
   p = sqrt(Polyhedron(poly)%vpsq(plane))
!
   s(1) = x/p
   s(2) = y/p
   s(3) = z/p
!
   end function getPlaneNormal
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInscrSphRadius(poly) result(rmt)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
!
   real (kind=RealKind) :: rmt
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getInscrSphRadius','Invalid polyhedron index',poly)
   endif
!
   rmt = Polyhedron(poly)%rmt
!
   end function getInscrSphRadius
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOutscrSphRadius(poly) result(rc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
!
   real (kind=RealKind) :: rc
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getOutscrSphRadius','Invalid polyhedron index',poly)
   endif
!
   rc = Polyhedron(poly)%rcirc
!
   end function getOutscrSphRadius
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWignerSeitzRadius(poly) result(rws)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
!
   real (kind=RealKind) :: rws
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getWignerSeitzRadius','Invalid polyhedron index',poly)
   endif
!
   rws = Polyhedron(poly)%rws
!
   end function getWignerSeitzRadius
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNeighborIndex(poly,plane) result(id)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: plane
   integer (kind=IntKind) :: id
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getNeighborIndex','Invalid polyhedron index',poly)
   else if (plane < 1 .or. plane > Polyhedron(poly)%NumPlanes) then
      call ErrorHandler('getNeighborIndex','Invalid plane index',plane)
   endif
!
   id = Polyhedron(poly)%id_planes(plane)
!
   end function getNeighborIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNeighborDistance(poly,plane) result(d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: plane
!
   real (kind=RealKind) :: d
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getNeighborDistance','Invalid polyhedron index',poly)
   else if (plane < 1 .or. plane > Polyhedron(poly)%NumPlanes) then
      call ErrorHandler('getNeighborDistance','Invalid plane index',plane)
   endif
!
   d = Polyhedron(poly)%NeighborDist(plane)
!
   end function getNeighborDistance
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolume(poly) result(v)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
!
   real (kind=RealKind) :: v
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getVolume','Invalid polyhedron index',poly)
   endif
!
   v = Polyhedron(poly)%Volume
!
   end function getVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInscrSphVolume(poly) result(v)
!  ===================================================================
   use MathParamModule, only : PI4, THIRD
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
!
   real (kind=RealKind) :: v
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getInscrSphVolume','Invalid polyhedron index',poly)
   endif
!
   v = PI4*THIRD*Polyhedron(poly)%rmt**3
!
   end function getInscrSphVolume
!  ===================================================================
!
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBoxVolume() result(v)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: v
   real (kind=RealKind) :: vec(3)
!
!  -------------------------------------------------------------------
   call vcross(UnitBox(1,1),UnitBox(1,2),vec)
!  -------------------------------------------------------------------
   v = abs(vec(1)*UnitBox(1,3)+vec(2)*UnitBox(2,3)+vec(3)*UnitBox(3,3))
!
   end function getBoxVolume
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCorners(poly) result(nc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: nc
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getNumCorners','Invalid polyhedron index',poly)
   endif
!
   nc = Polyhedron(poly)%NumCorners
!
   end function getNumCorners
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumEdges(poly) result(ne)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ne
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getNumEdges','Invalid polyhedron index',poly)
   endif
!
   ne = Polyhedron(poly)%NumEdges
!
   end function getNumEdges
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPlane(poly) result(pp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: np
!
   real (kind=RealKind), pointer :: pp(:,:)
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getPlane','Invalid polyhedron index',poly)
   endif
!
   np = Polyhedron(poly)%NumPlanes
   pp => Polyhedron(poly)%vplane(1:3,1:np)
!
   end function getPlane
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCorner(poly) result(pc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: nc
!
   real (kind=RealKind), pointer :: pc(:,:)
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getCorner','Invalid polyhedron index',poly)
   endif
!
   nc = Polyhedron(poly)%NumCorners
   pc => Polyhedron(poly)%corner(1:3,1:nc)
!
   end function getCorner
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEdge(poly) result(pe)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ne
!
   real (kind=RealKind), pointer :: pe(:,:,:)
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('getEdge','Invalid polyhedron index',poly)
   endif
!
   ne = Polyhedron(poly)%NumEdges
   pe => Polyhedron(poly)%edge(1:3,1:ne,1:2)
!
   end function getEdge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isExternalPoint(poly,x,y,z) result(a)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: isigma
!
   real (kind=RealKind), intent(in) :: x, y, z
   real (kind=RealKind) :: tol
!
   logical :: a
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('isExternalPoint','Invalid polyhedron index',poly)
   endif
!
   tol = TOLERANCE
!  -------------------------------------------------------------------
   call chkpnt(x,y,z,Polyhedron(poly)%vplane,                         &
                     Polyhedron(poly)%vpsq,                           &
                     Polyhedron(poly)%NumPlanes,tol,isigma)
!  -------------------------------------------------------------------
   if (isigma == 0) then
      a = .true.
   else
      a = .false.
   endif
!
   end function isExternalPoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSurfacePoint(poly,x,y,z) result(a)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: isigma
!
   real (kind=RealKind), intent(in) :: x, y, z
   real (kind=RealKind) :: tol
!
   logical :: a
!
   if (poly < 1 .or. poly > NumLocalPolyhedra) then
      call ErrorHandler('isExternalPoint','Invalid polyhedron index',poly)
   endif
!
   tol = HALF*TOLERANCE
!  -------------------------------------------------------------------
   call chkpnt(x,y,z,Polyhedron(poly)%vplane,                         &
                     Polyhedron(poly)%vpsq,                           &
                     Polyhedron(poly)%NumPlanes,tol,isigma)
!  -------------------------------------------------------------------
   if (isigma == 0) then
      a = .false.
      return
   endif
!  -------------------------------------------------------------------
   call chkpnt(x,y,z,Polyhedron(poly)%vplane,                         &
                     Polyhedron(poly)%vpsq,                           &
                     Polyhedron(poly)%NumPlanes,-tol,isigma)
!  -------------------------------------------------------------------
   if (isigma /= 0) then
      a = .false.
      return
   endif
   a = .true.
!
   end function isSurfacePoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPointLocationFlag(poly,x,y,z) result(f)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: f
!
   real (kind=RealKind), intent(in) :: x, y, z
!
   if (isExternalPoint(poly,x,y,z)) then
      f = -1
   else if (isSurfacePoint(poly,x,y,z)) then
      f = 0
   else
      f = 1
   endif
!
   end function getPointLocationFlag
!  ===================================================================
end module PolyhedraModule
