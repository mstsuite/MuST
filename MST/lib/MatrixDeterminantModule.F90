module MatrixDeterminantModule
    use KindParamModule, only : IntKind, RealKind, CmplxKind
!
    use MathParamModule, only : one, zero, cone, czero
!
public :: MtxDet
    interface MtxDet
       module procedure MtxDet_r, MtxDet_c
    end interface
!
private
!
contains
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxDet_r(n,am,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n
!
    real (kind=RealKind), intent(in) :: am(:,:)
    real (kind=RealKind), intent(out) :: det
!
    integer (kind=IntKind) :: info, i, sign_P
    integer (kind=IntKind), allocatable :: ipiv(:)
!
    real (kind=RealKind), allocatable :: bm(:,:)
!
    allocate( bm(n,n), ipiv(n) )
!
    bm(1:n,1:n) = am(1:n,1:n)
    ipiv = 0
!
    call dgetrf(n,n,bm,n,ipiv,info)
!
    sign_P = 1
    det = ONE
    do i = 1,n
       if ( ipiv(i)/=i ) then
          sign_P = -sign_P
       endif
       det = det*bm(i,i)
    enddo
    det = sign_P*det
!
    deallocate( bm, ipiv )
!
    end subroutine MtxDet_r
!   ==================================================================
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxDet_c(n,am,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n
!
    complex (kind=CmplxKind), intent(in) :: am(:,:)
    complex (kind=CmplxKind), intent(out) :: det
!
    integer (kind=IntKind) :: info, i, sign_P
    integer (kind=IntKind), allocatable :: ipiv(:)
!
    complex (kind=CmplxKind), allocatable :: bm(:,:)
!
    allocate( bm(n,n), ipiv(n) )
!
    bm(1:n,1:n) = am(1:n,1:n)
    ipiv = 0
!
    call zgetrf(n,n,bm,n,ipiv,info)
!
    sign_P = 1
    det = CONE
    do i = 1,n
       if ( ipiv(i)/=i ) then
          sign_P = -sign_P
       endif
       det = det*bm(i,i)
    enddo
    det = sign_P*det
!
    deallocate( bm, ipiv )
    end subroutine MtxDet_c
!   ==================================================================


end module MatrixDeterminantModule
