! this subroutine forces the number of OpenMP threads to be the same as that
! on MPI process mpi_id (>= 0)
subroutine force_openmp(mpi_id)
  use MPPModule, only : getMyPE,getNumPEs,bcastMessage
  use KindParamModule, only: IntKind
  use ErrorHandlerModule, only : ErrorHandler
  implicit none
  integer (kind=IntKind), intent(in), optional :: mpi_id  ! MPI process id
#ifdef _OPENMP
  integer(Kind=IntKind) :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  integer(Kind=IntKind) :: n,np,ns,n_openmp,MyPE,NumPEs

  if (present(mpi_id)) then
     if (mpi_id < 0) then
        call ErrorHandler('force_openmp','Invalid MPI process ID',mpi_id)
     else
        n = mpi_id
     endif
  else
     n = 0
  endif
!
  n_openmp=0
  MyPE=getMyPE()
  NumPEs = getNumPEs()
!
  if (MyPE == n) then
     write(6,'("MyPE,NumPEs=",2i5)') MyPE,NumPEs
  endif

  if (NumPEs.le.1) return
  
  if (MyPE == n) then
!$OMP PARALLEL PRIVATE(np,ns) SHARED (n_openmp)
     np=OMP_GET_THREAD_NUM()
     ns=OMP_GET_NUM_THREADS()
     if (np.eq.0) n_openmp=ns
!$OMP BARRIER
!$OMP END PARALLEL
  endif
  
  call bcastMessage(n_openmp,n)
  
  call OMP_SET_NUM_THREADS(n_openmp)
!
  if (MyPE == n) then
     write(6,'(a,i5)')'The number openMP threads is set to ',n_openmp
  endif

#endif  
end subroutine force_openmp
