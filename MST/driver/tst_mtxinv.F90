program tst_mtxinv
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ten2m8, sqrtm1, czero, one, two, three
   use MatrixInverseModule, only : MtxInv_GE
!
   implicit   none
!
   integer (kind=IntKind) :: i,j,k
   integer (kind=IntKind), parameter :: n = 81
!
   real (kind=RealKind) :: seed
   real (kind=RealKind) :: random
!
   complex (kind=CmplxKind) :: a(n,n)
   complex (kind=CmplxKind) :: b(n,n)
   complex (kind=CmplxKind) :: c(n,n)
! 
   seed=one
   do j=1,n
      do i=1,n
         a(i,j)=three*random(seed)+two*random(seed)*sqrtm1
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call MtxInv_GE(n,a,b)
!  -------------------------------------------------------------------
!
   do j=1,n
      do i=1,n
         c(i,j)=czero
         do k=1,n
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
         enddo
         if(abs(c(i,j)).ge.ten2m8) write(6,'(2i5,2d20.13)')i,j,c(i,j)
      enddo
   enddo
!
   stop 'ok'
end program tst_mtxinv
