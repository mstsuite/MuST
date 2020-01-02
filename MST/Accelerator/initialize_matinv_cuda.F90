!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initialize_matinv_cuda(nthreads,num_blocks,max_block_sz)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: nthreads
   integer (kind=IntKind), intent(in) :: num_blocks(nthreads)
   integer (kind=IntKind), intent(in) :: max_block_sz(nthreads)
!
#if defined(CULA) || defined(CUDA)
   integer (kind=IntKind) :: culaStatus, info, dev, minV, curV
   character infostr(100)

   integer*8 :: get_dev_ipvt
   integer*8 :: get_dev_a

   integer (kind=IntKind) :: cula_initialize,cublas_alloc,            &
                             cula_get_cuda_minimum_version,           &
                             cula_get_cuda_runtime_version,           &
                             cula_get_cuda_driver_version,            &
                             cula_get_cublas_minimum_version,         &
                             cula_get_cublas_runtime_version,         &
                             cula_get_device_info,                    &
                             cula_get_executing_device
!
   integer (kind=IntKInd) :: sizeof_Z = 16
   integer (kind=IntKind) :: sizeof_I = 8
!
#if defined(CULA) || defined(CULA_FALLBACK) 
   integer (kind=IntKind) :: devSZ
   common /Accelerator_CULA/ devSZ
!
   culaStatus=cula_initialize()
!  write(*,*) 'cula_initialize (',mynod,') = ',culaStatus
   if (culaStatus.ne.0) then 
      write(*,*) 'cula_initialize  = ',culaStatus
      call cula_get_status_string() 
   endif
   info = cula_get_executing_device(dev)
   info = cula_get_device_info(dev, infostr, 100)
   write(*,*) 'Accelerator found: ',infostr
   !JL Available in R12 and up
   !Check CUDA Version
   minV = cula_get_cuda_minimum_version()
   curV = cula_get_cuda_runtime_version()
   if (curV.lt.minV) then
      write(*,*)'CUDA Runtime ',minV,' required. Found ',curV
      stop
   endif
   !Check Driver Version
   curV = cula_get_cuda_driver_version()
   if (curV.lt.minV) then
      write(*,*)'CUDA Driver ',minV,' required. Found ',curV
      stop
   endif
   !Check CUBLAS Version
   minV = cula_get_cublas_minimum_version()
   curv = cula_get_cublas_runtime_version()
   if (curV.lt.minV) then
      write(*,*)'CUBLAS ',minV,' required. Found ',curV
      stop
   endif
#endif
#elif defined(LIBSCI)
   external  libsci_acc_init!,cublas_init
!  integer cublas_init,ierr
   call libsci_acc_init()
!  ierr = cublas_init()
!  if (ierr.ne.0) then
!     write(*,*)'CUBLAS_INIT',0
!  endif
#endif
#if defined(CULA) || defined(CUDA)
   call initialize_dstore(nthreads,num_blocks,max_block_sz)
#endif
!  ===================================================================
   end subroutine initialize_matinv_cuda
