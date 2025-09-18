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
   subroutine invert_lsms_matrix(pBlockMatrix, kkrsz_ns, pBigMatrix, dsize)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kkrsz_ns,dsize
!
   complex (kind=CmplxKind), intent(inout) :: pBlockMatrix(kkrsz_ns,kkrsz_ns)
   complex (kind=CmplxKind), intent(in) :: pBigMatrix(dsize,dsize)
   integer :: i,j
!
!  print *,'In invertMatrixBlock_CUDA'
#ifdef CUDA
!  -------------------------------------------------------------------
   call cusolver_lsms_c(dsize,pBigMatrix,kkrsz_ns,pBlockMatrix)
!  -------------------------------------------------------------------
#else
   print *,'call zblock_lu_CPU'
!  -------------------------------------------------------------------
   call zblock_lu_CPU(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
!  -------------------------------------------------------------------
#endif
!
   end subroutine invert_lsms_matrix
