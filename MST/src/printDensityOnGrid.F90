!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDensityOnGrid(na, lprint)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use ChargeDensityModule, only : isSphericalChargeDensity, getChargeDensityAtPoint
!
   use ValenceDensityModule, only : getValenceElectronDensityAtPosi
!
   use PublicTypeDefinitionsModule, only : UniformGridStruct
!
   use Uniform3DGridModule, only : getNumGridPoints
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, lprint
   integer (kind=IntKind) :: id, ng
!
   real (kind=RealKind), allocatable :: denOnGrid(:)
!
   interface 
      subroutine constructDataOnGrid(grid_name, value_name, value_type, getData, den, lmax, spin)
         use KindParamModule, only : IntKind, RealKind
         use PublicTypeDefinitionsModule, only : UniformGridStruct
         implicit none
         character (len=*), intent(in) :: grid_name
         character (len=*), intent(in) :: value_name
         character (len=*), intent(in) :: value_type
         real (kind=RealKind), intent(out) :: den(:)
         integer (kind=IntKind), intent(in), optional :: lmax, spin
!
         interface
            function getData( dname, id, ia, r, jmax_in, n, grad ) result(v)
               use KindParamModule, only : IntKind, RealKind
               implicit none
               character (len=*), intent(in) :: dname
               integer (kind=IntKind), intent(in) :: id, ia
               real (kind=RealKind), intent(in) :: r(3)
               real (kind=RealKind), intent(out), optional :: grad(3)
               integer (kind=IntKind), intent(in), optional :: jmax_in, n
               real (kind=RealKind) :: v
            end function getData
         end interface
      end subroutine constructDataOnGrid
   end interface 
!
!  ng = gp%NumLocalGridPoints
   ng = getNumGridPoints('Visual',local_grid=.true.)
   allocate( denOnGrid(ng) )
!
   call constructDataOnGrid( 'Visual', 'Charge', 'Valence', getChargeDensityAtPoint, denOnGrid )  ! total charge
   call printDataOnGrid( 'Visual', 'Charge', 'Valence', denOnGrid, lprint)
   call constructDataOnGrid( 'Visual', 'Charge', 'TotalNew', getChargeDensityAtPoint, denOnGrid )  ! total charge
   call printDataOnGrid( 'Visual', 'Charge', 'TotalNew', denOnGrid, lprint)
   do id = 1, na
      call printDataOnLine('ValenDen',id,getValenceElectronDensityAtPosi)
   enddo
   if ( .not.isSphericalChargeDensity() ) then
      call constructDataOnGrid( 'Visual', 'Charge', 'Pseudo', getChargeDensityAtPoint, denOnGrid )
      call printDataOnGrid( 'Visual', 'Charge', 'Pseudo', denOnGrid, lprint)
   endif
!
   deallocate( denOnGrid )
!
   end subroutine printDensityOnGrid
