module IntegerFactorsModule
   use KindParamModule, only : IntKind
!
public :: initIntegerFactors,          &
          endIntegerFactors,           &
          pushIntegerFactorsToAccel,   &
          deleteIntegerFactorsOnAccel, &
          isIntegerFactorsInitialized
!
public :: lofk, mofk, jofk, lofj, mofj, kofj, m1m, bofk
!
   integer (kind=IntKind), allocatable :: m1m(:),                     &
                                          lofk(:), mofk(:), jofk(:),  &
                                          lofj(:), mofj(:), kofj(:), bofk(:)
private
   logical :: Initialized = .false.
   integer (kind=IntKind) :: lmax
   integer (kind=IntKind) :: NumInits
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initIntegerFactors(lmax_in)
!  ===================================================================
   integer(kind=IntKind), intent(in) :: lmax_in
!
   logical :: isResetOn
!
   if ( Initialized ) then
      NumInits = NumInits + 1             ! count initialization calls
      if ( lmax_in<=lmax ) then
         return
      else
         isResetOn=.true.
!        -------------------------------------------------------------
         call endIntegerFactors(isResetOn)
!        -------------------------------------------------------------
      endif
   else
      NumInits = 1
   endif
!
   lmax = lmax_in
!  -------------------------------------------------------------------
   call genFactors(lmax)
!  -------------------------------------------------------------------
   Initialized = .true.
!
   end subroutine initIntegerFactors
!  =================================================================== 
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endIntegerFactors(isResetOn)
!  ===================================================================
!
   logical, optional :: isResetOn
!
   if ( present(isResetOn) ) then
!     ----------------------------------------------------------------
      deallocate( lofk, mofk, jofk, lofj, mofj, kofj, m1m, bofk )
!     ----------------------------------------------------------------
      return
   endif
!
   NumInits = NumInits - 1
   if (NumInits > 0) then
      return
   else if ( .not.Initialized ) then
      return
   endif
!
!  -------------------------------------------------------------------
   deallocate( lofk, mofk, jofk, lofj, mofj, kofj, m1m, bofk )
!  -------------------------------------------------------------------
   Initialized = .false.
!
   end subroutine endIntegerFactors
!  =================================================================== 
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genFactors(l0)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: l0
!
   integer (kind=IntKind) :: kmax, jmax, l, m, kl, jl, n
!
   kmax=(l0+1)*(l0+1)
   jmax=(l0+1)*(l0+2)/2
!
   allocate( lofk(kmax), mofk(kmax), jofk(kmax), bofk(kmax) )
   allocate( lofj(jmax), mofj(jmax), kofj(jmax) )
   allocate( m1m(-l0:l0) )
!  ===================================================================
!  calculate the factors: lofk, mofk, jofk, lofj, mofj, kofj, m1m, and bofk.
!  ===================================================================
   jl = 0
   kl = 0
   do l = 0,l0
      n = (l+1)*(l+2)/2 - l
!
      do m = -l,l,1
         kl = kl+1
         lofk(kl) = l
         mofk(kl) = m
         bofk(kl) = kl-2*m
         if ( m>=0 ) then
            jofk(kl) = n + m
         else
            jofk(kl) = n - m
         endif
      enddo
!
      do m = 0,l
         jl=jl+1
         lofj(jl)=l
         mofj(jl)=m
         kofj(jl)=(l+1)*(l+1)-l+m
      enddo
   enddo
!
!  ===================================================================
!  calculate the factor (-1)**m and store in m1m(-lmax:lmax)..........
!  ===================================================================
   m1m(0)=1
   do m = 1, l0
      m1m(m)  = -m1m(m-1)
      m1m(-m) =  m1m(m)
   enddo
!
   end subroutine genFactors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pushIntegerFactorsToAccel()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: jmax, kmax
!
   jmax = (lmax+1)*(lmax+2)/2
   kmax = (lmax+1)**2
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call initialize_integerfactors(lmax,jmax,kmax)
   call push_integerfactors(lofk,mofk,lofj,mofj,kofj,m1m,bofk)
!  -------------------------------------------------------------------
#endif
!
   end subroutine pushIntegerFactorsToAccel
!  =================================================================== 
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteIntegerFactorsOnAccel()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: jmax, kmax
!
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call delete_integerfactors()
!  -------------------------------------------------------------------
#endif
!
   end subroutine deleteIntegerFactorsOnAccel
!  =================================================================== 
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isIntegerFactorsInitialized() result(t)
!  =================================================================== 
   implicit none
!
   logical :: t
!
   t = Initialized
!
   end function isIntegerFactorsInitialized
!  =================================================================== 
end module IntegerFactorsModule
