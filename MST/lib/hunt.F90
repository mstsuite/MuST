subroutine hunt(n,xx,x,jlo)
   use KindParamModule, only : IntKind, RealKind
!
   implicit none
!
   integer (kind=IntKind), intent(inout) :: jlo
   integer (kind=IntKind),intent(in) :: n
   integer (kind=IntKind) :: inc, jhi, jm
!
   real (kind=RealKind), intent(in) :: xx(n)
   real (kind=RealKind), intent(in) :: x
!
   logical :: ascnd
! 
!  *******************************************************************
!  Taken from numerical recipes. Modified by Yang Wang
!  Modified according to the Fortran 90 version
!
!  Given an array xx(1:n), and given a value x, returns a value jlo
!  such that x is between xx(jlo) and xx(jlo+1). xx must be monotonic,
!  either increasing or decreasing. jlo=0 or jlo=n is returned to
!  indicate that x is out of range. jlo on input is taken as the initial
!  guess for jlo on output.
!  *******************************************************************
!  n=size(xx,1)
   ascnd=(xx(n) >= xx(1))
   if (jlo <= 0 .or. jlo > n) then
      jlo=0
      jhi=n+1
   else
      inc=1
      if (x >= xx(jlo) .eqv. ascnd) then
         do
            jhi=jlo+inc
            if(jhi > n)then
               jhi=n+1
               exit
            else if (x < xx(jhi) .eqv. ascnd) then
               exit
            else 
               jlo=jhi
               inc=inc+inc
            endif
         enddo
      else
         jhi=jlo
         do
            jlo=jhi-inc
            if (jlo < 1) then
               jlo=0
               exit
            else if (x >= xx(jlo) .eqv. ascnd) then
               exit
            else
               jhi=jlo
               inc=inc+inc
            endif
         enddo
      endif
   endif
   do
      if (jhi-jlo <= 1) then
         if (x == xx(1)) then
            jlo=1
         else if (x == xx(n)) then
            jlo=n-1
         endif
         exit
      else
         jm=(jhi+jlo)/2
         if (x >= xx(jm) .eqv. ascnd) then
            jlo=jm
         else
            jhi=jm
         endif
      endif
   enddo
end subroutine hunt
