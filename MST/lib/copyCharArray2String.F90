!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine copyCharArray2String(a,b,n)
!     ================================================================
      use KindParamModule, only : IntKind
!
      implicit   none
!
      integer (kind=IntKind), intent(in) :: n
      integer (kind=IntKind) :: i,m
!
!     character (len=1), intent(in) :: a(n)
      character (len=1), intent(in) :: a(:)
      character (len=*), intent(out) :: b
!
      b = ' '
      m = min(n,len(b))
      do i=1,m
         b(i:i)=a(i)
      enddo
!
      end subroutine copyCharArray2String
