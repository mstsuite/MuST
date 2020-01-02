!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copyString2CharArray(a,b,n)
!     ================================================================
      use KindParamModule, only : IntKind
!
      implicit   none
!
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind) :: i,m
!
      character (len=*), intent(in) :: a
      character (len=1), intent(out) :: b(:)
!     character (len=1), intent(out) :: b(n)
!
      m = min(n,len(a))
      do i=1,m
         b(i)=a(i:i)
      enddo
      do i=m+1,n
         b(i) = ' '
      enddo
!
      end subroutine copyString2CharArray
