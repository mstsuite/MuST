      subroutine accelerator_finalize()
      implicit none
#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_CUDA_C)
      include 'CULA_Common.h'

!      call cublas_free(devA)
!      call cublas_free(devIPVT)
#if defined(ACCELERATOR_CULA) || defined(CULA_FALLBACK)
      call cula_shutdown()
#endif
#elif defined(ACCELERATOR_LIBSCI)
!      external  libsci_acc_finalize,cublas_shutdown
!      integer cublas_shutdown,ierr
!      call libsci_acc_finalize()
!      ierr = cublas_shutdown()
!      if (ierr.ne.0) then
!        write(*,*)'CUBLAS_SHUTDOWN',0
!      endif
#else
!      write(*,*) "UNKNOWN ACCLERATOR in accelerator_finalize()"
#endif
      end subroutine

