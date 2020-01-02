module StatusModule
   use KindParamModule, only : IntKind
   implicit none
!
public :: initStatus,       &
          endStatus,        &
          registerError,    &
          chechNodesStatus
!
private
   character (len=1), allocatable :: ErrorFlag(:)
!
   integer (kind=IntKind) :: PrintLevel
   integer (kind=IntKind) :: NumIndex
   integer (kind=IntKind) :: MyIndex
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initStatus(iprint)
!  ===================================================================
   use MPPModule, only : NumPEs, MyPE
   implicit none
   integer (kind=IntKind) :: i
!
   NumIndex = NumPEs
   MyIndex = MyPE+1
!
   allocate( ErrorFlag(NumPEs) )
   do i=1, NumPEs
      ErrorFlag(i) = 'n'
   enddo
!
   PrintLevel = iprint
!
   end subroutine initStatus
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endStatus()
!  ===================================================================
   implicit none
!
   deallocate( ErroFlag )
   PrintLevel = -1
!
   end subroutine endStatus
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine registerError()
!  ===================================================================
   implicit none
!
   ErrorFlag(MyIndex) = 'n'
!
   end subroutine registerError
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkNodesStatus()
!  ===================================================================
   implicit none
!
!
   end subroutine checkNodesStatus
!  ===================================================================
end module StatusModule
