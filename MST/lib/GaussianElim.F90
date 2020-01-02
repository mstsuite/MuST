!  ==================================================================
   subroutine GaussianElim(A,n)
!  ==================================================================
   use KindParamModule, only : IntKind, CmplxKind
   use MathParamModule, only : ZERO
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: n
   complex(kind=CmplxKind), intent(inout) :: A(n,n)
!
   integer (kind=IntKind) :: k, i, j
   complex(kind=CmplxKind) :: fact
!
!  Convert to upper triangular form
!
   do k = 1,n-1
!
      if ( abs(A(k,k)) /= ZERO ) then
         do i = k+1,n
            fact = A(i,k)/A(k,k)
            do j = k+1,n
               A(i,j) = A(i,j) -A(k,j)*fact
            enddo
         enddo
      else
         write (6,*) 'GaussianElim:: Zero Element find in line:', k
         stop
      endif
!
   enddo
!
!  ==================================================================
   end subroutine GaussianElim
