!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageRFuncXProc(gid,n,f)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use GroupCommModule, only : getNumPEsInGroup
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: gid, n
   integer (kind=IntKind) :: np
!
   real (kind=RealKind), intent(inout) :: f(n)
   real (kind=RealKind) :: rfac
!
   np = getNumPEsInGroup(gid)
   rfac = real(np,kind=RealKind)
   f = f/rfac
!  -------------------------------------------------------------------
   call GlobalSumInGroup(gid,f,n)
!  -------------------------------------------------------------------
!
   end subroutine averageRFuncXProc
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageCFuncXProc(gid,n,f)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use GroupCommModule, only : getNumPEsInGroup
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: gid, n
   integer (kind=IntKind) :: np
!
   complex (kind=CmplxKind), intent(inout) :: f(n)
   complex (kind=CmplxKind) :: cfac
!
   np = getNumPEsInGroup(gid)
   cfac = real(np,kind=RealKind)
   f = f/cfac
!  -------------------------------------------------------------------
   call GlobalSumInGroup(gid,f,n)
!  -------------------------------------------------------------------
!
   end subroutine averageCFuncXProc
