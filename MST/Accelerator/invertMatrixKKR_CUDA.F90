!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine invertMatrixKKR_CUDA(pBigMatrix, dsize)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: dsize
!
   complex (kind=CmplxKind), intent(inout) :: pBigMatrix(dsize,dsize)
!
!  print *,'In invertMatrixBlock_CUDA'
   call cusolver_kkr_c(dsize,pBigMatrix)
   end subroutine invertMatrixKKR_CUDA
