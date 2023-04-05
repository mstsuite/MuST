subroutine wholeinv_CPU(m,a,tau00)
    use KindParamModule, only : IntKind, RealKind, CmplxKind
    integer (kind=IntKind), intent(in) :: m
    complex (kind=CmplxKind),intent(inout) :: a(m,m)
    complex (kind=CmplxKind),intent(inout) :: tau00(m,m)
    integer (kind=IntKind):: info
    integer (kind=IntKind), target :: ipvt(m)
    integer (kind=IntKind), pointer :: p_ipvt(:)
    p_ipvt => ipvt(1:m)
    call zgetrf(m,m,a,m,p_ipvt,info)
    call zgetrs("n",m,m,a,m,p_ipvt,tau00,m,info)
end subroutine wholeinv_CPU

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invertMatrixBlock_CUDA(my_atom, pBlockMatrix, kkrsz_ns, pBigMatrix, dsize)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: my_atom,kkrsz_ns,dsize
!
   complex (kind=CmplxKind), intent(inout) :: pBlockMatrix(kkrsz_ns,kkrsz_ns)
   complex (kind=CmplxKind), intent(in) :: pBigMatrix(dsize,dsize)
   complex (kind=CmplxKind) :: tmpmatrixC(kkrsz_ns,kkrsz_ns),tmpmatrixC_inv(kkrsz_ns,kkrsz_ns)
!   complex (kind=CmplxKind) :: a1(dsize,dsize),z(dsize,dsize)
   integer :: i,j
!
!  print *,'In invertMatrixBlock_CUDA'
#ifdef CUDA
   !print *,'call zblock_lu_cuda_c, nblk =',nblk
!  -------------------------------------------------------------------
!   a1 = pBigMatrix
   do i=1,kkrsz_ns
     do j=1,kkrsz_ns
       pBlockMatrix(i,j)=pBigMatrix(i,j)
     enddo
   enddo
   call cusolver_invwhole_c(dsize,pBigMatrix,kkrsz_ns,tmpmatrixC)!tempmatrixC is the first block of the inverse of BigMatrix
!   z = matmul(a1,pBigMatrix)
!   print *,"z(1,1)=",z(1,1),"z(1,10)=",z(1,10),"sum of inverted matrix =",sum(z)
   do i=1,kkrsz_ns
     do j=1,kkrsz_ns
       if (i == j) then
         tmpmatrixC_inv(i,j)=cmplx(1,0)
       else
         tmpmatrixC_inv(i,j)=cmplx(0,0)
       endif
!       tmpmatrixC(i,j)=pBigMatrix(i,j)
     enddo
   enddo
!   print *,"sum of initialized tmpmatrixC_inv =",sum(tmpmatrixC_inv)
!  inverse tmpmatrixC on CPU is efficient enough
   call wholeinv_CPU(kkrsz_ns,tmpmatrixC,tmpmatrixC_inv)
   pBlockMatrix = pBlockMatrix - tmpmatrixC_inv
!   print *,"In side invertMatrixBlock_CUDA, pBlockMatrix(:, 1) = ",pBlockMatrix(:,1)
!   call zblock_lu_cuda_c(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
!  -------------------------------------------------------------------
#else
   print *,'call zblock_lu_CPU'
!  -------------------------------------------------------------------
   call zblock_lu_CPU(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
!  -------------------------------------------------------------------
#endif
!
   end subroutine invertMatrixBlock_CUDA
