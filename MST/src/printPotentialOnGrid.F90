!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPotentialOnGrid(na,lprint)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use PotentialGenerationModule, only : getPotentialAtPoint
!
   use PublicTypeDefinitionsModule, only : UniformGridStruct
!
   use Uniform3DGridModule, only : getUniform3DGrid, getNumGridPoints
!   
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, lprint
   integer (kind=IntKind) :: id, ng
!
   real (kind=RealKind), allocatable :: denOnGrid(:)
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
   gp => getUniform3DGrid('Visual')
!
!  ng = gp%NumLocalGridPoints
   ng = getNumGridPoints('Visual',local_grid=.true.)
   allocate( denOnGrid(ng) )
!
   call constructDataOnGrid( 'Visual', 'Potential', 'Coulomb', getPotentialAtPoint, denOnGrid )
   call printDataOnGrid( gp, 'Potential', 'Coulomb', denOnGrid, lprint)
!
   call constructDataOnGrid( 'Visual', 'Potential', 'Exchg', getPotentialAtPoint, denOnGrid )
   call printDataOnGrid( gp, 'Potential', 'Exchg', denOnGrid, lprint)
!
   call constructDataOnGrid( 'Visual', 'Potential', 'Total', getPotentialAtPoint, denOnGrid )
   call printDataOnGrid( gp, 'Potential', 'Total', denOnGrid, lprint)
!
   do id = 1, na
      call printDataOnLine('Visual','Potential',id,getPotentialAtPoint)
      call printDataOnLine('Visual','Coulomb',id,getPotentialAtPoint)
      call printDataOnLine('Visual','En_Exchg',id,getPotentialAtPoint)
   enddo
!
   deallocate( denOnGrid )
!
   end subroutine printPotentialOnGrid
