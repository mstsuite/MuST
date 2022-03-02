!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine decompose3Index(i,n1,n2,n3,i1,i2,i3)
!  ===================================================================
!  Purpose:
!     Given index i, determine the corresponding index i1, i2, i3,
!     where the index i in the following loop
!
!        do i = 1, n1*n2*n3
!           ...,...
!        enddo
!
!     corresponds to the index i in the following 3-level loop
!
!        i = 0
!        do i3 = 1, n3
!           do i2 = 1, n2
!              do i1 = 1, n1
!                 i = i + 1
!                 ...,...
!              enddo
!           enddo
!        enddo
!
!  *******************************************************************
!  ===================================================================
   use KindParamModule, only : IntKind
   use ErrorHandlerModule, only : ErrorHandler
   implicit none
!
   integer (kind=IntKind), intent(in)  :: i, n1, n2, n3
   integer (kind=IntKind), intent(out) :: i1, i2, i3
!
   if (i < 1) then
      call ErrorHandler('decompose3Index','index i < 1',i)
   else if (i > n1*n2*n3) then
      call ErrorHandler('decompose3Index','index i > n1*n2*n3',i,n1*n2*n3)
   endif
!
   i1 = mod(i-1,n1) + 1
   i2 = mod((i-i1)/n1,n2) + 1
   i3 = ((i-i1)/n1 - (i2-1))/n2 + 1
!
   end subroutine decompose3Index
