!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initialize_matinv_accel(nthreads,num_blocks,max_block_sz)
!  ===================================================================
   use KindParamModule, only : IntKind
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nthreads
   integer (kind=IntKind), intent(in) :: num_blocks(nthreads)
   integer (kind=IntKind), intent(in) :: max_block_sz(nthreads)
!
!  print *,'Num local atoms = ',nthreads
!  print *,'LIZ size = ',num_blocks(1)
!  print *,'Max block size = ',max_block_sz(1)
#ifdef CUDA
   call initialize_matinv_cuda(nthreads,num_blocks,max_block_sz)
#elif defined(MIC)
#else
#endif
   end subroutine initialize_matinv_accel
