!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDensityOnGrid(na, iprint)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, ONE
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use ChargeDensityModule, only : isSphericalChargeDensity, getChargeDensityAtPoint
!
   use ValenceDensityModule, only : getValenceElectronDensityAtPosi
!
   use PublicTypeDefinitionsModule, only : UniformGridStruct
!
   use Uniform3DGridModule, only : getUniform3DGrid, getNumGridPoints
   use Uniform3DGridModule, only : getGridIndex
!
   use MPPModule, only : MyPE, NumPEs, AnyPE
   use MPPModule, only : sendPackage, packMessage
   use MPPModule, only : recvPackage, unpackMessage
   use MPPModule, only : setCommunicator, resetCommunicator, GlobalMin
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : getGroupCommunicator
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, iprint
   integer (kind=IntKind) :: i, id, ng, n, IO_PE, comm, ng_local, ig
   integer (kind=IntKind), allocatable :: GridIndex_local(:) ! The global grid index of the uniform
                                                             ! grid points on the local process
   integer (kind=IntKind) :: GroupID, NumPEsInGroup, MyPEinGroup
   integer (kind=IntKind), allocatable :: checked(:) ! This is for sanity check
!
   real (kind=RealKind), allocatable, target :: denOnGrid(:)
   real (kind=RealKind), allocatable :: denOnGrid_local(:)   ! The density value on the uniform grid
                                                             ! points on the local process
!   
   type (UniformGridStruct), pointer :: gp
!
   interface 
      subroutine constructDataOnGrid(grid_name, value_name, value_type, getData, den, lmax, spin, tol_in)
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : UniformGridStruct
         implicit none
         character (len=*), intent(in) :: grid_name
         character (len=*), intent(in) :: value_name
         character (len=*), intent(in) :: value_type
         real (kind=RealKind), intent(out) :: den(:)
         real (kind=RealKind), intent(in), optional :: tol_in
         integer (kind=IntKind), intent(in), optional :: lmax, spin
!
         interface
            function getData( dname, id, ia, r, tol, jmax_in, n, grad, trunc ) result(v)
               use KindParamModule, only : IntKind, RealKind
               implicit none
               character (len=*), intent(in) :: dname
               integer (kind=IntKind), intent(in) :: id, ia
               real (kind=RealKind), intent(in) :: r(3), tol
               real (kind=RealKind), intent(out), optional :: grad(3)
               integer (kind=IntKind), intent(in), optional :: jmax_in, n
               real (kind=RealKind) :: v
               logical, optional, intent(in) :: trunc
            end function getData
         end interface
      end subroutine constructDataOnGrid
   end interface 
!
   interface 
      subroutine printDataOnGrid(gp, value_name, value_type, den, iprint)
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : UniformGridStruct
         implicit none
         character (len=*), intent(in) :: value_name
         character (len=*), intent(in) :: value_type
         real (kind=RealKind), intent(in), target :: den(:)
         integer (kind=IntKind), intent(in), optional :: iprint
         type (UniformGridStruct), intent(in) :: gp
      end subroutine printDataOnGrid
   end interface 
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
   comm = getGroupCommunicator(GroupID)
   if (MyPEinGroup == 0) then
      IO_PE = MyPE
   else
      IO_PE = NumPEs
   endif
!  -------------------------------------------------------------------
   call GlobalMin(IO_PE)
!  -------------------------------------------------------------------
!
   gp => getUniform3DGrid('Visual')
!
!  ng = gp%NumLocalGridPoints
   ng_local = getNumGridPoints('Visual',local_grid=.true.)
   ng = getNumGridPoints('Visual')
   allocate( denOnGrid(ng), denOnGrid_local(ng), GridIndex_local(ng) )
   allocate( checked(ng) )
!
   denOnGrid = ZERO; denOnGrid_local = ZERO; checked = 0
!  -------------------------------------------------------------------
   call constructDataOnGrid( 'Visual', 'Charge', 'Valence', getChargeDensityAtPoint, denOnGrid_local )  ! valence charge
!  -------------------------------------------------------------------
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
   if (MyPEinGroup == 0) then
      do i = 1, ng_local
         ig = getGridIndex(gp,i)     
         denOnGrid(ig) = denOnGrid_local(i)
         checked(ig) = 1
      enddo
      do n = 1, NumPEsInGroup-1
         call recvPackage(100001,AnyPE)
         call unpackMessage(ng_local)
         call unpackMessage(GridIndex_local,ng_local)
         call unpackMessage(denOnGrid_local,ng_local)
         do i = 1, ng_local
            ig = GridIndex_local(i)
            denOnGrid(ig) = denOnGrid_local(i)
            checked(ig) = 1
         enddo
      enddo
!     ================================================================
!     Performing sanity check
!     ================================================================
      do i = 1, ng
         if (checked(i) == 0) then
            call WarningHandler('printDensityOnGrid','Valence density at a grid point is missing',i)
         endif
      enddo
   else
      do i = 1, ng_local
         GridIndex_local(i) = getGridIndex(gp,i)     
      enddo
      call packMessage(ng_local)
      call packMessage(GridIndex_local,ng_local)
      call packMessage(denOnGrid_local,ng_local)
      call sendPackage(100001,0)
   endif
!  -------------------------------------------------------------------
   call resetCommunicator(sync=.true.)
!  -------------------------------------------------------------------
   if (MyPE == IO_PE) then
!     ----------------------------------------------------------------
      call printDataOnGrid(gp, 'Charge', 'Valence', denOnGrid, iprint)
!     ----------------------------------------------------------------
   endif
!
   denOnGrid = ZERO; denOnGrid_local = ZERO; checked = 0
!  -------------------------------------------------------------------
   call constructDataOnGrid( 'Visual', 'Charge', 'TotalNew', getChargeDensityAtPoint, denOnGrid_local )  ! total charge
!  -------------------------------------------------------------------
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
   if (MyPEinGroup == 0) then
      do i = 1, ng_local
         ig = getGridIndex(gp,i)     
         denOnGrid(ig) = denOnGrid_local(i)
         checked(ig) = 1
      enddo
      do n = 1, NumPEsInGroup-1
         call recvPackage(100001,AnyPE)
         call unpackMessage(ng_local)
         call unpackMessage(GridIndex_local,ng_local)
         call unpackMessage(denOnGrid_local,ng_local)
         do i = 1, ng_local
            ig = GridIndex_local(i)
            denOnGrid(ig) = denOnGrid_local(i)
            checked(ig) = 1
         enddo
      enddo
!     ================================================================
!     Performing sanity check
!     ================================================================
      do i = 1, ng
         if (checked(i) == 0) then
            call WarningHandler('printDensityOnGrid','Total density at a grid point is missing',i)
         endif
      enddo
   else
      do i = 1, ng_local
         GridIndex_local(i) = getGridIndex(gp,i)     
      enddo
      call packMessage(ng_local)
      call packMessage(GridIndex_local,ng_local)
      call packMessage(denOnGrid_local,ng_local)
      call sendPackage(100001,0)
   endif
!  -------------------------------------------------------------------
   call resetCommunicator(sync=.true.)
!  -------------------------------------------------------------------
   if (MyPE == IO_PE) then
!     ----------------------------------------------------------------
      call printDataOnGrid(gp, 'Charge', 'Total', denOnGrid, iprint)
!     ----------------------------------------------------------------
   endif
!
if (.false.) then
   do id = 1, na
      call printDataOnLine('Visual','ValenDen',id,getValenceElectronDensityAtPosi)
   enddo
!
   if ( .not.isSphericalChargeDensity() ) then
      call constructDataOnGrid( 'Visual', 'Charge', 'Pseudo', getChargeDensityAtPoint, denOnGrid )
      call printDataOnGrid(gp, 'Charge', 'Pseudo', denOnGrid, iprint)
   endif
endif
!
   deallocate( denOnGrid, denOnGrid_local, GridIndex_local, checked )
!
   call resetCommunicator(sync=.true.)
!
   end subroutine printDensityOnGrid
