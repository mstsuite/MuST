   program tst_gaussq
!
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
!
   use MathParamModule, only : zero
   use MathParamModule, only : one
   use MathParamModule, only : two
   use MathParamModule, only : five
!
   implicit none
!
   integer (kind=IntKind) :: kind
   integer (kind=IntKind), parameter :: n=11
   integer (kind=IntKind) :: i
!
   real (kind=RealKind) :: t(n), w(n)
!
   call gauleg(-one,one,t,w,n)
   write(6,'(''Output from gauleg'')')
   do i=1,n
      write(6,'('' t, w = '',2d20.13)')t(i),w(i)
   enddo 
   write(6,'(/)')
!
   do kind=1,4
!ywg  call gaussq(kind, n, kpts, endpts, b, t, w)
      call gaussq(kind, n, 0, -one, one, t, w)
      write(6,'('' kind = '',1i5)')kind
      do i=1,n
         write(6,'('' t, w = '',2d20.13)')t(i),w(i)
      enddo 
   enddo
!
   stop 'Ok'
   end program tst_gaussq
