module SystemSymmetryModule
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, ONE, TWO, TEN2m10, TEN2m8, SQRTm1, HALF
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public:: initSystemSymmetry,  &
         getSymmetryFlags,    &
         setSymmetryFlags,    &
         resetSymmetryFlags,  &
         calSymmetryFlags,    &
         endSystemSymmetry
!
   interface getSymmetryFlags
      module procedure getSymmetryFlags0, getSymmetryFlags1
   end interface
!
   interface setSymmetryFlags
      module procedure setSymmetryFlags0, setSymmetryFlags1
   end interface
!
private
!
   logical :: Initialized = .false.
!
   type SymmetryFlags
      logical :: isFlagSet
      integer (kind=IntKind) :: lmax_mad
      integer (kind=IntKind) :: lmax_step
      integer (kind=IntKind) :: lmax
      integer (kind=IntKind) :: jmax
      integer (kind=IntKind), pointer :: flags(:)
   end type SymmetryFlags
!
   type (SymmetryFlags), allocatable, target :: SymmFlags(:)
!
   integer (kind=IntKind) :: NumLocalSites
   integer (kind=IntKind) :: NumGlobalSites
   integer (kind=IntKind) :: NumSiteTypes
!
   integer (kind=IntKind) :: jmax_max
   integer (kind=IntKind) :: lmax_max
!
   integer (kind=IntKind), allocatable :: print_level(:)
!
   real (kind=RealKind), parameter :: PrimeValue = 11.0d0
   real (kind=RealKind), parameter :: tol_abs = TEN2m8
   real (kind=RealKind), parameter :: tol_fraction = TEN2m8
   real (kind=RealKind), allocatable :: tol_l(:)
!
contains
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSystemSymmetry( nga, nla, lmax_mad, lmax_step, iprint )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nga, nla
   integer (kind=IntKind), intent(in) :: lmax_mad(nla), lmax_step(nla), iprint(nla)
!
   integer (kind=IntKind) :: id, lmax_id, jmax_id, ll, powerl
!
   NumGlobalSites = nga
   NumLocalSites  = nla
!
   allocate( SymmFlags(1:NumLocalSites) )
   allocate( print_level(1:NumLocalSites) )
!
   lmax_max = 0
   jmax_max = 0
   do id = 1,NumLocalSites
      SymmFlags(id)%isFlagSet = .false.
      lmax_id = lmax_step(id) + lmax_mad(id)
      jmax_id = (lmax_id+1)*(lmax_id+2)/2
      SymmFlags(id)%lmax_mad = lmax_mad(id)
      SymmFlags(id)%lmax_step = lmax_step(id)
      SymmFlags(id)%lmax = lmax_id
      SymmFlags(id)%jmax = jmax_id
      allocate( SymmFlags(id)%flags(1:jmax_id) )
      SymmFlags(id)%flags(1:jmax_id) = 0
      jmax_max = max(jmax_max,jmax_id)
      lmax_max = max(lmax_max,lmax_id)
      print_level(id) = iprint(id)
   enddo
!
   allocate( tol_l(0:lmax_max) )
   tol_l(0) = tol_abs
   do ll = 1, lmax_max
      powerl = ll/2
      tol_l(ll) = tol_abs/(10.0d0**powerl)
   enddo
!
   Initialized = .true.
!
   end subroutine initSystemSymmetry
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSystemSymmetry()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id
!
   do id = 1, NumLocalSites
      deallocate(SymmFlags(id)%flags)
   enddo
   deallocate( SymmFlags )
   deallocate( print_level )
!
   deallocate( tol_l )
!
   Initialized = .false.
!
   end subroutine endSystemSymmetry
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSymmetryFlags(des)
!  ===================================================================
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
   implicit none
!
   character (len=*), optional :: des
!
   integer (kind=IntKind) :: ni, p, jmax_step, jmax_mad, jmax
   integer (kind=IntKind) :: jls, kls, jlm, klm, i3, kl3, l3, jl3, nj3_i3
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
   integer (kind=IntKind), pointer :: kj3_i3(:)
   integer (kind=IntKind), pointer :: pflag(:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
   real (kind=RealKind), pointer :: cgnt_i3(:)
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (.not.present(des)) then
      call calTruncSymmetryFlags()
      call calMadSymmetryFlags()
!
      nj3 => getNumK3()
      kj3 => getK3()
      cgnt => getGauntFactor()
      do ni = 1, NumLocalSites
         pflag => SymmFlags(ni)%flags(:)
         jmax_step = (SymmFlags(ni)%lmax_step+1)*(SymmFlags(ni)%lmax_step+2)/2
         jmax_mad = (SymmFlags(ni)%lmax_mad+1)*(SymmFlags(ni)%lmax_mad+2)/2
         do jls = 1, jmax_step
            kls = kofj(jls)
            do jlm = 1, jmax_mad
               klm = kofj(jlm)
               p = min(pflag(jlm),pflag(jls))
               nj3_i3 = nj3(kls,klm)
               kj3_i3 => kj3(1:nj3_i3,kls,klm)
               cgnt_i3 => cgnt(1:nj3_i3,kls,klm)
               do i3 = 1, nj3_i3
                  kl3 = kj3_i3(i3)
                  l3 = lofk(kl3)
                  if (mofk(kl3) >= 0 .and. l3 > lofj(jls)+lofj(jlm)) then
                     jl3 = (l3+1)*(l3+2)/2
                     if (abs(p*cgnt_i3(i3)) > tol_l(l3)) then
                        pflag(jl3) = max(pflag(jl3),pflag(jlm),pflag(jls))
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
!
   else if ( nocaseCompare(des,'Madelung') ) then
      call calMadSymmetryFlags()
   else
      call calTruncSymmetryFlags()  ! Using truncation function to determine the symmetry
   endif
!
   if ( print_level(1) >= 1 ) then
      write(6,'(/,a,/)')"calSymmetryFlags::   Flags"
      do ni = 1, NumLocalSites
         jmax = size(SymmFlags(ni)%flags)
         write(6,*)"  l    m   flag"
         do jls = 1, jmax
            write(6,'(3i5)') lofj(jls), mofj(jls), SymmFlags(ni)%flags(jls)
         enddo
         write(6,'(60(''-''),/)')
      enddo
   endif
!
   end subroutine calSymmetryFlags
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calTruncSymmetryFlags()
!  ===================================================================
   use IntegerFactorsModule, only : lofj, mofj
   use StepFunctionModule, only : getStepFunction
   use SystemModule, only : getAtomType
!
   implicit none
!
   integer (kind=IntKind) :: ni, nr, ir, jmax, jl, l, sf, izamax
   integer (kind=IntKind), pointer :: pSymmFlags(:)
!
   real (kind=RealKind) :: sr, si
   complex (kind=CmplxKind), pointer :: stepf(:,:)
!
   do ni = 1, NumLocalSites
      jmax = (SymmFlags(ni)%lmax_step+1)*(SymmFlags(ni)%lmax_step+2)/2
      pSymmFlags => SymmFlags(ni)%flags(1:jmax)
      pSymmFlags(1) = 1
      stepf => getStepFunction(ni)
      if ( jmax > size(stepf,2) ) then
         call ErrorHandler('calTruncSymmetryFlags','In initSystemSymmetry, jmax > jmax_step', &
                             jmax,size(stepf,2))
      endif
!
      nr = size(stepf,1)
      do jl = 1, jmax
!        ============================================================
!        set the flags for the non-zero component of
!        spherical harmonic expansion:
!           flags == 0 the component is zero - no contribution
!           flags == 1 the component has the real part non-zero
!           flags == 2 the component has the imaginary part non-zero
!           flags == 3 the component contributes with
!                      both real and imaginary parts
!        ============================================================
         l = lofj(jl)
         sf = 0
         LOOP_ir: do ir = 1, nr
!           if (abs(stepf(ir,jl)) > tol_l(l)) then
            if (abs(stepf(ir,jl)) > tol_abs*10.0d0) then
               sr = real(stepf(ir,jl),kind=RealKind); si = aimag(stepf(ir,jl))
               if ( abs(sr) > tol_l(l) .and. abs(si) <= TEN2m10 ) then
!              if ( abs(sr) > tol_abs  .and. abs(si) <= TEN2m10 ) then
                  sf = 1
               else if ( abs(si) <= TEn2m10 .and. abs(si) > tol_l(l) ) then
!              else if ( abs(si) <= TEn2m10 .and. abs(si) > tol_abs ) then
                  sf = 2
               else
                  sf = 3
               endif
               exit LOOP_ir
            endif
         enddo LOOP_ir
         pSymmFlags(jl) = max(pSymmFlags(jl),sf)
      enddo
      SymmFlags(ni)%isFlagSet = .true.
!
      if ( print_level(ni) >= 0 ) then
         write(6,'(/,a,/)')"calSymmetryFlags::   Flags, StepFunction"
         write(6,*)"  l   m flag          Re[stepf]                  Im[stepf]        "
         do jl = 1, jmax
            ir = izamax(nr,stepf(1:nr,jl),1)
            write(6,'(3i4,2(1x,d23.14))') lofj(jl), mofj(jl), pSymmFlags(jl), stepf(ir,jl)
         enddo
         write(6,'(60(''-''),/)')
      endif
   enddo
!
   end subroutine calTruncSymmetryFlags
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calMadSymmetryFlags()
!  ===================================================================
   use MadelungModule, only : getDLMatrix, getDLFactor, getMadelungMatrix
!
   use SystemModule, only : getAtomType, getAtomicNumber
!
   implicit none
!
   integer (kind=IntKind) :: ni, id, ll, ml, jl, kl, m1m, lmax, jmax, sf
   integer (kind=IntKind), pointer :: pAtomType(:)
   integer (kind=IntKind), pointer :: pSymmFlags(:)
!
   real (kind=RealKind) :: q0, ql_r, ql_i, a0
   real (kind=RealKind), pointer :: dlf(:), mad(:)
   real (kind=RealKind), target  :: dlf_0(1:1)
!
   complex (kind=CmplxKind) :: ql, sumat
   complex (kind=CmplxKind), allocatable :: ql_ni(:)
   complex (kind=CmplxKind), pointer :: dlm(:,:)
   complex (kind=CmplxKind), allocatable, target :: dlm_0(:,:)
!
   allocate(ql_ni(jmax_max))
!
   pAtomType => getAtomType()
!
   do ni  = 1, NumLocalSites
      if (SymmFlags(ni)%lmax_mad == 0) then
         allocate( dlm_0(NumGlobalSites,1) )
         exit
      endif
   enddo
!
   do ni  = 1, NumLocalSites
      lmax = SymmFlags(ni)%lmax_mad
      jmax = (lmax+1)*(lmax+2)/2
      pSymmFlags => SymmFlags(ni)%flags(1:jmax)
      pSymmFlags(1) = 1
      if ( lmax==0 ) then
         mad => getMadelungMatrix(ni)
         dlm_0 (1:NumGlobalSites,1) = cmplx(mad,ZERO,kind=RealKind)
         a0 = ONE
         dlm => dlm_0
         dlf_0(1) = ONE
      else
         dlm => getDLMatrix(ni,a0)
      endif
!
      m1m = 1
      ql_ni = CZERO
      do ll = 0, lmax
         do ml = 0,ll
            kl = (ll+1)*(ll+1) - ll + ml
            jl = (ll+1)*(ll+2)/2 -ll + ml
            sf = 0
!
            if ( lmax==0 ) then
               dlf => dlf_0
            else
               dlf => getDLFactor(jl)
            endif
            sumat = CZERO
            do id = 1,NumGlobalSites
               q0 = getAtomicNumber(id)
               if (abs(q0) < 1.0) then
!                 q0 = q0 + id/100.0d0
                  q0 = q0 + ONE/100.0d0
               endif
!              q0 = max(q0,HALF)
!               q0 = real(PrimeValue+100*pAtomType(id),kind=RealKind)
!               q0 = real(PrimeValue+100*id,kind=RealKind)
               sumat = sumat + dlm(id,kl)*q0/a0**ll
            enddo
            ql = dlf(1)*sumat
            ql_ni(jl) = ql
!
!           ---------------------------------------------------------
!           set the flags for the non-zero component of
!           spherical harmonic expansion:
!           flags == 0 the component is zero - no contribution
!           flags == 1 the component has the real part non-zero
!           flags == 2 the component has the imaginary part non-zero
!           flags == 3 the component contributes with
!                                    both real and imaginary parts
!           ---------------------------------------------------------
!
            if ( abs(ql)/tol_l(ll) > 1 ) then
               sf = 3
            endif
!
            if ( pSymmFlags(jl) == 3 ) then
               ql_r = real(ql,kind=RealKind)
               ql_i = real(-sqrtm1*ql,kind=RealKind)
!              if ( abs(ql_r) /= ZERO  .and. abs(ql_i)/abs(ql_r) < tol_fraction ) then
               if ( abs(ql_r) > TEN2m10  .and. abs(ql_i)/abs(ql_r) < tol_fraction ) then
                  sf = 1
               else if ( abs(ql_i) > TEN2m10  .and. abs(ql_r)/abs(ql_i) < tol_fraction ) then
                  sf = 2
               endif
            endif
            pSymmFlags(jl) = max(pSymmFlags(jl),sf)
         enddo
         m1m = -m1m
      enddo
      SymmFlags(ni)%isFlagSet = .true.
!
      if ( print_level(ni) >= 0 ) then
         jl = 0
         write(6,'(/,a,/)')"calSymmetryFlags::   Flags, q_Value"
         write(6,*)"  l   m flag          Re[ql]                     Im[ql]           "
         do ll = 0,lmax
            do ml = 0,ll
               jl = jl+1
               write(6,'(3i4,2(1x,d23.14))') ll, ml, pSymmFlags(jl), ql_ni(jl)
            enddo
         enddo
         write(6,'(60(''-''),/)')
      endif
   enddo
!
   nullify(dlm)
   if ( allocated(dlm_0) ) then
      deallocate( dlm_0 )
   endif
   deallocate( ql_ni )
!
   end subroutine calMadSymmetryFlags
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSymmetryFlags0(id)                         result(flag)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind), pointer :: flag(:)
!
   if (id<1 .or. id>NumLocalSites) then
      call ErrorHandler('getSymmFlags0','Invalid atom index',id)
   endif
   if ( .not.SymmFlags(id)%isFlagSet ) then
      call WarningHandler('getSymmFlags0',"Flags have not been initialized")
   endif
   jmax = SymmFlags(id)%jmax
   flag => SymmFlags(id)%flags(1:jmax)
!
   end function getSymmetryFlags0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSymmetryFlags1(id,jl)                         result(flag)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, jl
!
   integer (kind=IntKind) :: flag
!
   if (id<1 .or. id>NumLocalSites) then
      call ErrorHandler('getSymmFlags1','Invalid atom index',id)
   endif
   if ( jl<0 .or. jl>SymmFlags(id)%jmax ) then
      call ErrorHandler('getSymmFlags1','Invalid jl',jl)
   endif
   if ( .not.SymmFlags(id)%isFlagSet ) then
      call WarningHandler('getSymmFlags1',"Flags have not been initialized")
   endif
!
   flag = SymmFlags(id)%flags(jl)
!
   end function getSymmetryFlags1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSymmetryFlags0(id,jmax,flags)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, jmax
   integer (kind=IntKind), intent(in) :: flags(:)
!
   integer (kind=IntKind) :: jl
!
   if (id<1 .or. id>NumLocalSites) then
      call ErrorHandler('setSymmFlags','Invalid atom index',id)
   endif
   if ( jmax < 0 .or. jmax > SymmFlags(id)%jmax ) then
      call WarningHandler('setSymmFlags',                             &
                       'Input jmax size greater than local one',jmax)
   endif
!
   SymmFlags(id)%flags = 0
   do jl = 1, min(jmax,SymmFlags(id)%jmax)
      SymmFlags(id)%flags(jl) = flags(jl)
   enddo
   SymmFlags(id)%isFlagSet = .true.
!
   end subroutine setSymmetryFlags0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSymmetryFlags1(id,jl,flag)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, jl
   integer (kind=IntKind), intent(in) :: flag
!
   if (id<1 .or. id>NumLocalSites) then
      call ErrorHandler('setSymmFlags','Invalid atom index',id)
   endif
   if ( jl < 0 .or. jl > SymmFlags(id)%jmax ) then
      call WarningHandler('setSymmFlags',                             &
                       'Input jmax size greater than local one',jl)
   endif
!
   SymmFlags(id)%flags(jl) = flag
!
   end subroutine setSymmetryFlags1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetSymmetryFlags()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id
!
   do id = 1,NumLocalSites
      SymmFlags(id)%flags = 3
      SymmFlags(id)%isFlagSet = .false.
   enddo
!
   end subroutine resetSymmetryFlags
!  ===================================================================
end module SystemSymmetryModule
