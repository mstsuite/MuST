      program kkrmat_compare
      implicit none
      integer n,n2,i,j,argc
      real*8 er,ei,er2,ei2
      real*8 ar1,ai1,ar2,ai2
      integer ii,jj
      complex*16 a
      character(len=64) :: argv(2)
      real*8 tol,err

      tol=0.001

      argc=command_argument_count()
      if(argc.ne.2) then
        write(*,*) "Expected two arguments: files to compare"
        stop
      end if
      call get_command_argument(1,argv(1))
      call get_command_argument(2,argv(2))

      write(*,*) "Comapring files '",trim(argv(1)),
     &           "' and '",trim(argv(2)),"'"

      open(unit=11,file=trim(argv(1)))
      open(unit=12,file=trim(argv(2)))

      read(11,*) n,er,ei
      read(12,*) n2,er2,ei2

      if(n.ne.n2) then
        write(*,*) "n(",trim(argv(1)),") = ",n," not equal to ",
     &             "n(",trim(argv(2)),") = ",n2
        stop
      end if

      if(er.ne.er2 .or. ei.ne.ei2) then
        write(*,*) "e(",trim(argv(1)),") = (",er,", ",ei,
     &             ") not equal to ",
     &             "e(",trim(argv(2)),") = (",er2,", ",ei2,")"
!       stop
      end if
      write(*,*) "Energy = (",er,", ",ei,")"

!      allocate(a(n,n))

      do i=1,n
        do j=1,n
          read(11,*) ii,jj,ar1,ai1
          read(12,*) ii,jj,ar2,ai2
          a=cmplx(ar1,ai1)-cmplx(ar2,ai2)
          if(abs(a).ne.0.0) then
             err=abs(a)/abs(cmplx(ar1,ai1)+cmplx(ar2,ai2))
             if(err.gt.tol) then
                write(*,*) i,j,ar1,ai1,ar2,ai2
             end if
          end if
       end do
      end do

      end
