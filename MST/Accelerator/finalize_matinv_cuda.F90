!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine finalize_matinv_cuda()
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
   integer (kind=IntKInd) :: info
!
#if defined(CULA) || defined(CUDA)
!
   integer (kind=IntKInd) :: sizeof_Z = 16
   integer (kind=IntKind) :: sizeof_I = 8
!
#if defined(CULA) || defined(CULA_FALLBACK)
   integer (kind=IntKind) :: devSZ
   common /Accelerator_CULA/ devSZ
!
   call cula_shutdown()
#endif
!
#elif defined(LIBSCI)
!  external  libsci_acc_finalize,cublas_shutdown
!  integer cublas_shutdown,ierr
!  call libsci_acc_finalize()
!  ierr = cublas_shutdown()
!  if (info.ne.0) then
!     write(*,*)'CUBLAS_SHUTDOWN',0
!  endif
#endif
!  ===================================================================
   end subroutine finalize_matinv_cuda

