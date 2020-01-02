!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine finalize_matinv_accel()
!  ===================================================================
   implicit none
!
#ifdef CUDA
   call finalize_matinv_cuda()
#elif defined(MIC)
#else
#endif
   end subroutine finalize_matinv_accel
