!  *******************************************************************
!  Function getPrimeFactors(num, f) 
!  Purpose: TO FIND THE PRIME FACTORS OF A NUMBER
!  Original Author : Louisda16th a.k.a Ashwith J. Rego
!  Modified by Yang Wang
!  Description: 
!    Input: num, the number to be factorized
!    Output: f,  the number of prime factors
!    This funtion returns a pointer to the integer array which contains
!    the prime factors
!  Algorithm is quite easy:
!    Start with 2, check whether 2 is a factor by seeing if 
!    MOD(<input_number>,2) is zero. If it is zero, then 2 becomes a factor.
!    If not, check with the next number.
!    When a factor is found, divide the given number with the factor found.
!    However, do not move to the next possible factor - a number can
!    occur more than once as a factor
!
!  Function getSubFactors(num, k, ns) 
!  Purpose: TO FIND all k FACTORS OF A NUMBER
!  Original Author : Yang Wang
!  Description: 
!    Input: num, the number to be factorized into k numbers
!    Input:  k
!    Output: ns,  the number of k factor groups
!    This funtion returns a pointer to the integer array which contains
!    all possible k factors, the 1st dimension is 1:k, and the 2nd dimension
!    is 1:ns
!  *******************************************************************
module PrimeFactorsModule
   use KindParamModule, only : IntKind
!
public :: getPrimeFactors, getSubFactors, cleanPrimeFactors
!
private
   integer (kind=IntKind) :: fmax = 0
   integer (kind=IntKind) :: fsize = 0
   integer (kind=IntKind), allocatable, target :: factors(:)
!
   integer (kind=IntKind) :: tsave = 0
   integer (kind=IntKind) :: hsave = 0
   integer (kind=IntKind) :: smax = 0
   integer (kind=IntKind) :: ssize = 0
   integer (kind=IntKind) :: kmax = 0
   integer (kind=IntKind) :: ksize = 0
   integer (kind=IntKind), allocatable, target :: sub_factors(:,:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPrimeFactors(num, f) result(pf)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
!
   IMPLICIT NONE
   INTEGER (kind=IntKind), INTENT(IN)  :: num  !input number
   INTEGER (kind=IntKind), INTENT(OUT) :: f
   INTEGER (kind=IntKind) :: i, n
   INTEGER (kind=IntKind), pointer :: pf(:)
!
   i = 2  !Eligible factor
   f = 1  !Number of factors
   if (num < 0) then
      n = -num !store input number into a temporary variable
   else if (num == 0) then
      call ErrorHandler('getPrimeFactors','The number is ZERO!')
   else
      n = num
   endif
   LOOP_do0: DO
      IF (MOD(n,i) == 0) THEN !If i divides 2, it is a factor
         f = f+1
         n = n/i
      ELSE
         i = i+1     !Not a factor. Move to next number
      ENDIF
      IF (n == 1) THEN
         !Since f is incremented after a factor is found
         f = f-1        !its value will be one more than the number of factors
         !Hence the value of f is decremented
         EXIT LOOP_do0
      ENDIF
   ENDDO LOOP_do0
!
   fsize = f
   if (fsize > fmax) then
      if (fmax > 0) then
         deallocate(factors)
      endif
      allocate(factors(1:fsize))
      fmax = fsize
   endif
!
   i = 2  !Eligible factor
   f = 1  !Number of factors
   if (num < 0) then
      n = -num !store input number into a temporary variable
   else
      n = num
   endif
   LOOP_do: DO
      IF (MOD(n,i) == 0) THEN !If i divides 2, it is a factor
         factors(f) = i
         f = f+1
         n = n/i
      ELSE
         i = i+1     !Not a factor. Move to next number
      ENDIF
      IF (n == 1) THEN
         !Since f is incremented after a factor is found
         f = f-1     !its value will be one more than the number of factors
         !Hence the value of f is decremented
         EXIT LOOP_do
      ENDIF
   ENDDO LOOP_do
!
   if (num < 0) then
      factors(1) = -factors(1)
   endif
!
   pf => factors(1:f)
!
   END function getPrimeFactors
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cleanPrimeFactors()
!  ===================================================================
   implicit none
!
   if (fmax > 0) then
      deallocate(factors)
   endif
   fmax = 0; fsize = 0
!
   if (smax > 0) then
      deallocate(sub_factors)
   endif
   smax = 0; ssize = 0
   kmax = 0; ksize = 0
   hsave = 0; tsave = 0
!
   end subroutine cleanPrimeFactors
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSubFactors(num, k, ns) result(ps)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in)  :: num, k  !input number
   integer (kind=IntKind), intent(out) :: ns
   integer (kind=IntKind) :: i, j, ib, ik, n, nf, m, f(1:k), r(1:k)
   integer (kind=IntKind), pointer :: pf(:)
   integer (kind=IntKind), pointer :: ps(:,:)
!
   logical :: mtc, is_same
!
   if (num < 0) then
      n = -num !store input number into a temporary variable
   else
      n = num
   endif
!
   pf => getPrimeFactors(n,nf)
!
   m = 2**nf - 1
!
   mtc=.false.
   ns = 0
   LOOP_do0: do
      call nexcom(m,k,f(1:k),mtc)
      ib = f(1)
      do ik = 2, k
         ib =ieor(ib,f(ik))
      enddo
      if (ib == m) then
         ns = ns + 1
      endif
      if (.not.mtc) then
         exit LOOP_do0
      endif
   enddo LOOP_do0
!
   if (ns > smax .or. k > kmax) then
      if (smax > 0) then
         deallocate(sub_factors)
      endif
      allocate(sub_factors(1:k,1:ns))
      smax = ns
      kmax = k
   endif
!
   mtc=.false.
   ns = 0
   LOOP_do1: do
      call nexcom(m,k,f(1:k),mtc)
      ib = f(1)
      do ik = 2, k
         ib =ieor(ib,f(ik))
      enddo
      if (ib == m) then
         do ik = 1, k
            r(ik) = 1
            i = f(ik)
            do j = 1,nf
               if (btest(i,j-1)) then
                  r(ik) = r(ik)*pf(j)
               endif
            enddo
         enddo
         is_same = .false.
         LOOP_j: do j = 1, ns
            is_same = .true.
            LOOP_ik: do ik = 1, k
               if (sub_factors(ik,j) /= r(ik)) then
                  is_same = .false.
                  exit LOOP_ik
               endif
            enddo LOOP_ik
            if (is_same) then
               exit LOOP_j
            endif
         enddo LOOP_j
         if (.not.is_same) then
            ns = ns + 1
            sub_factors(1:k,ns) = r(1:k)
         endif
      endif
      if (.not.mtc) then
         exit LOOP_do1
      endif
   enddo LOOP_do1
   ssize = ns
   ksize = k
!
   if (num < 0) then
      sub_factors(1,1:ns) = -sub_factors(1,1:ns)
   endif
!
   ps => sub_factors(1:k,1:ns)
!
   end function getSubFactors
!  ===================================================================
!
!  *******************************************************************
!  subroutine nexcom is copied from:
!     http://www.cs.sunysb.edu/~algorith/implement/wilf/distrib/processed/
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine nexcom(n,k,r,mtc)
!  ===================================================================
   implicit none
!
   logical, intent(inout) :: mtc
!
   integer (kind=IntKind), intent(in) :: n, k
   integer (kind=IntKind), intent(inout) :: r(k)
   integer (kind=IntKind) :: i
!
   if (.not.mtc) then
      r(1) = n
      tsave = n
      hsave = 0
      if (k /= 1) then
         do i = 2, k
            r(i) = 0
         enddo
      endif
   else
      if (tsave > 1) hsave = 0
      hsave = hsave + 1
      tsave = r(hsave)
      r(hsave) = 0
      r(1) = tsave - 1
      r(hsave+1) = r(hsave+1) + 1
   endif
   mtc = (r(k) /= n)
!
   end subroutine nexcom
!  ===================================================================
end module PrimeFactorsModule
