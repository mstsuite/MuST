! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine vcross(v1,v2,v)
! ====================================================================
!
  use KindParamModule, only : RealKind
!
  implicit none
!
  real (kind=RealKind), intent(in) :: v1(3), v2(3)
  real (kind=RealKind), intent(out) :: v(3)
!
! ********************************************************************
! returns:
!           ->  ->   -> 
!           v = v1 x v2
!
! ********************************************************************
!
  v(1)=v1(2)*v2(3)-v1(3)*v2(2)
  v(2)=v1(3)*v2(1)-v1(1)*v2(3)
  v(3)=v1(1)*v2(2)-v1(2)*v2(1)
!
  end subroutine vcross
