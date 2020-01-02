!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine zblock_lu_accel(a,lda,na,blk_sz,ld_blk,nblk,idcol,blk1,ipvt,mp,k)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lda, na, ld_blk, blk1, mp
   integer (kind=IntKind), intent(inout) :: nblk
   integer (kind=IntKind), intent(inout) :: blk_sz(ld_blk), idcol(blk1)
   integer (kind=IntKind), intent(out) :: ipvt(mp), k
!
   complex (kind=CmplxKind), intent(inout) :: a(lda,na)
!
!  print *,'In zblock_lu_accel:'
!  print *,'a(1,1) =',a(1,1)
!  print *,'a(17,16) =',a(33,32)
!  print *,'blk_sz(1) =',blk_sz(1)
!  print *,'blk_sz(2) =',blk_sz(2)
!  print *,'idcol(1) =',idcol(1)
!  print *,'idcol(16) = ',idcol(16)
!  print *,'mp = ',mp
#ifdef CUDA
   print *,'call zblock_lu_cuda_c, nblk =',nblk
!  -------------------------------------------------------------------
   call zblock_lu_cuda_c(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
!  -------------------------------------------------------------------
#elif defined(MIC)
   print *,'call zblock_lu_mic'
!  -------------------------------------------------------------------
   call zblock_lu_mic(a,lda,blk_sz,nblk_t,ipvt,mp,idcol,k)
!  -------------------------------------------------------------------
#else
   print *,'call zblock_lu_CPU'
!  -------------------------------------------------------------------
   call zblock_lu_CPU(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
!  -------------------------------------------------------------------
#endif
!
   end subroutine zblock_lu_accel
