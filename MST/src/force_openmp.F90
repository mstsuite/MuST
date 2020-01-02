! this subroutine forces OpenMP threads to be the same as on the 1st processor
subroutine force_openmp()
  use MPPModule, only : getMyPE,getNumPEs,bcastMessage
   use KindParamModule, only: IntKind
  implicit none
#ifdef _OPENMP
  integer(Kind=IntKind) :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  integer(Kind=IntKind) :: np,ns,n_openmp,MyPe,NumPEs

  n_openmp=0
  MyPe=getMyPE()
  NumPEs = getNumPEs()
  write (*,'("MyPe,NumPes=",2i5)') MyPe,NumPes

  if (NumPEs.le.1) return
  
  if (MyPE == 0) then
!$OMP PARALLEL PRIVATE(np,ns) SHARED (n_openmp)
     np=OMP_GET_THREAD_NUM()
     ns=OMP_GET_NUM_THREADS()
     if (np.eq.0) n_openmp=ns
!$OMP BARRIER
!$OMP END PARALLEL
  endif
  
  call bcastMessage(n_openmp,0)
  
  call OMP_SET_NUM_THREADS(n_openmp)

#endif  
end subroutine force_openmp
