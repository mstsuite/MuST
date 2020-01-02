! this routine prints out the number of threads the cose is running with.
! i f this not not the expected number d o the following:
! - make sure you have set the system variable, e.g.
!         export OMP_NUM_THREADS=4
! - if the number of threads is different on different processors, run
!   before this call the folloging routine 
!         call force_openmp()
!
subroutine print_threads(myunit)
  use KindParamModule, only: IntKind
  implicit none
  integer(Kind=IntKind), intent(in) :: myunit
#ifdef _OPENMP
  integer(Kind=IntKind) :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  integer(Kind=IntKind) :: np,ns,n_openmp
  character(len=100)    :: ts
#endif

  write(myunit,'(/,80(''-''))')
  write(myunit,'(/,25x,a)') ' *****************************'
  write(myunit,'(25x,a)')   ' * Output from print_threads *'
  write(myunit,'(25x,a,/)') ' *****************************'
  !   write(myunit,'(/,80(''=''))')

! =============================================================  
#ifdef _OPENMP
  
!$OMP PARALLEL PRIVATE(np,ns) SHARED (n_openmp)
  np=OMP_GET_THREAD_NUM()
  ns=OMP_GET_NUM_THREADS()
  if (np.eq.0) n_openmp=ns
!$OMP BARRIER
!$OMP END PARALLEL
  
  write (ts,*) n_openmp
  ts=adjustl(ts)
  write(myunit,'(t17,3a)')  &
       & "code is compiled using threads, running with ",&
       & trim(ts)," thread(s)"
#else
  write(myunit,'(t23,a)')  "code is compiled without using threads"
#endif
! =============================================================  
  write(myunit,'(/,80(''=''))')
 
end subroutine print_threads
