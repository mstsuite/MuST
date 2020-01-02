module GauntFactorsModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
   use IntegerFactorsModule, only : lofk, mofk, m1m
   use MathParamModule, only : ZERO, PI4, ten2m10, ten2m14
!
public :: initGauntFactors, &
          getNumK3,         &
          getK3,            &
          getGauntFactor,   &
          endGauntFactors,  &
          pushGauntFactorsToAccel,   &
          deleteGauntFactorsOnAccel, &
          isInitialized
!
   interface getNumK3
      module procedure getNumK3_0, getNumK3_2
   end interface
!
   interface getK3
      module procedure getK3_0, getK3_2
   end interface
!
   interface getGauntFactor
      module procedure getGauntFactor0, getGauntFactor2, getGauntFactor3
   end interface
!
private
!
   integer (kind=IntKind) :: lmax_sav
   integer (kind=IntKind) :: kmax_sav
   integer (kind=IntKind) :: MaxJ3
   integer (kind=IntKind), allocatable, target :: nj3(:,:)
   integer (kind=IntKind), allocatable, target :: kj3(:,:,:)
!
   real (kind=RealKind) :: tol = ten2m14
   real (kind=RealKind), pointer :: clm(:)
   real (kind=RealKind), allocatable, target :: cgnt(:,:,:)
!
   logical :: Initialized = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: NumInits
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initGauntFactors(lmax, istop, iprint)
!  ===================================================================
   use SphericalHarmonicsModule, only : getClm
   use IntegerFactorsModule, only : initIntegerFactors
   implicit   none
!
   character (len=*) :: istop
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind) :: jmax,kmax,lmax2,jmax2,kmax2
!
   if (Initialized) then
      NumInits = NumInits + 1      ! count for number of times initGauntFactors
                                   ! have been called.
      if (lmax <= lmax_sav) then
         return
      else
!        -------------------------------------------------------------
         call endGauntFactors()
!        -------------------------------------------------------------
      endif
   else
      NumInits = 1
   endif
!
   Initialized=.true.
   print_level = iprint
   stop_routine = istop
!
   kmax=(lmax+1)**2
   jmax=(lmax+1)*(lmax+2)/2
!
   lmax_sav=lmax
   kmax_sav=kmax
   MaxJ3 = lmax+1
!
   lmax2=2*lmax
   kmax2=(lmax2+1)**2
   jmax2=(lmax2+1)*(lmax+1)
!
!  -------------------------------------------------------------------
   allocate( nj3(kmax,kmax)        )
   allocate( kj3(MaxJ3,kmax,kmax)  )
   allocate( cgnt(MaxJ3,kmax,kmax) )
   nj3(1:kmax,1:kmax)          = 0
   kj3(1:MaxJ3,1:kmax,1:kmax)  = 0
   cgnt(1:MaxJ3,1:kmax,1:kmax) = ZERO
!  -------------------------------------------------------------------
!
!  ===================================================================
!  set up factors lofk, mofk, etc..., for l <= lmax2
!  -------------------------------------------------------------------
   call initIntegerFactors(lmax2)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  get prefators of complex spherical harmonics for l <= lmax2
!  -------------------------------------------------------------------
   clm => getClm(lmax2)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  set up gaunt numbers.[N.B. this routine needs the clm above],
!
!                   L       4pi  ->   ->   * ->     ->
!  cgnt(j,L',L") = C      = int do*Y (o )*Y (o )*Y (o )
!                   L',L"    0      L'     L"     L
!
!  L', L" <= lmax
!
!  L = kj3(j,L',L"), j=1,2,...,nj3(L',L") <= lmax+1
!  -------------------------------------------------------------------
   call GauntFactor(lmax,lmax)
!  -------------------------------------------------------------------
!
   nullify( clm )
!
   end subroutine initGauntFactors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumK3_0() result(pnj3)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: pnj3(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getNumK3','GauntFactorsModule is not initialized')
   endif
!
   pnj3 => nj3(:,:)
!
   end function getNumK3_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumK3_2(k1,k2) result(pnj3)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: k1, k2
   integer (kind=IntKind) :: pnj3
!
   if (.not.Initialized) then
      call ErrorHandler('getNumK3','GauntFactorsModule is not initialized')
   else if (k1 > kmax_sav .or. k1 < 1) then
      call ErrorHandler('getNumK3','K1 index out of bound',k1)
   else if (k2 > kmax_sav .or. k2 < 1) then
      call ErrorHandler('getNumK3','K2 index out of bound',k2)
   endif
!
   pnj3 = nj3(k1,k2)
!
   end function getNumK3_2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getK3_0() result(pkj3)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: pkj3(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getK3','GauntFactorsModule is not initialized')
   endif
!
   pkj3 => kj3(:,:,:)
!
   end function getK3_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getK3_2(k1,k2) result(pkj3)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: k1, k2
   integer (kind=IntKind), pointer :: pkj3(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getK3','GauntFactorsModule is not initialized')
   else if (k1 > kmax_sav .or. k1 < 1) then
      call ErrorHandler('getK3','K1 index out of bound',k1)
   else if (k2 > kmax_sav .or. k2 < 1) then
      call ErrorHandler('getK3','K2 index out of bound',k2)
   endif
!
   pkj3 => kj3(1:MaxJ3,k1,k2)
!
   end function getK3_2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGauntFactor0() result(pcgnt)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: pcgnt(:,:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getGauntFactor',                             &
                        'GauntFactorsModule is not initialized')
   endif
!
   pcgnt => cgnt(:,:,:)
!
   end function getGauntFactor0
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGauntFactor2(k1,k2) result(pcgnt)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: k1, k2
!
   real (kind=RealKind), pointer :: pcgnt(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getGauntFactor',                             &
                        'GauntFactorsModule is not initialized')
   else if (k1 > kmax_sav .or. k1 < 1) then
      call ErrorHandler('getGauntFactor','K1 index out of bound',k1)
   else if (k2 > kmax_sav .or. k2 < 1) then
      call ErrorHandler('getGauntFactor','K2 index out of bound',k2)
   endif
!
   pcgnt => cgnt(1:MaxJ3,k1,k2)
!
   end function getGauntFactor2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGauntFactor3(k1,k2,k3) result(cg)
!  ===================================================================
!
!  given k1, k2, and k3, returns:
!
!    cg = int_{4pi} dr { Y_{k1}(r) * Y^{*}_{k2}(r) * Y_{k3}(r) }
!
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: k1, k2, k3
   integer (kind=IntKind) :: j3
!
   real (kind=RealKind) :: cg
!
   if (.not.Initialized) then
      call ErrorHandler('getGauntFactor',                             &
                        'GauntFactorsModule is not initialized')
   else if (k1 > kmax_sav .or. k1 < 1) then
      call ErrorHandler('getGauntFactor','K1 index out of bound',k1)
   else if (k2 > kmax_sav .or. k2 < 1) then
      call ErrorHandler('getGauntFactor','K2 index out of bound',k2)
   endif
!
   cg = ZERO
   do j3 = 1, nj3(k1,k2)
      if (kj3(j3,k1,k2) == k3) then
         cg = cgnt(j3,k1,k2)
         return
      endif
   enddo
!
   end function getGauntFactor3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endGauntFactors()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
   implicit   none
!
   NumInits = NumInits - 1
   if (NumInits > 0) then
      return
   else if ( .not.Initialized ) then
      return
   endif
!
!  -------------------------------------------------------------------
   deallocate( nj3, kj3, cgnt )
!  -------------------------------------------------------------------
   call endIntegerFactors()
!  -------------------------------------------------------------------
!
   Initialized = .false.
!
   end subroutine endGauntFactors
!  ===================================================================
!
!  *******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genFactors(lmax)
!  ====================================================================
   implicit   none 
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: kl
!
!  ===================================================================
!  set up factors lofk and mofk for l up to lmax...................
!  ===================================================================
   kl=0
   do l=0,lmax
      do m=-l,l
         kl=kl+1
         lofk(kl)=l
         mofk(kl)=m
      enddo
   enddo
!
!  ===================================================================
!  load the factors (-1)**m into m1m(-lmax:lmax)...................
!  ===================================================================
   m1m(0) = 1
   do m= 1,lmax
      m1m(m)= - m1m(m-1)
      m1m(-m)= m1m(m)
   enddo
!
   end subroutine genFactors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine GauntFactor(lmax_1,lmax_2)
!  ===================================================================
   use LegendreModule, only : legendre
!
   implicit   none
!
   character (len=11), parameter :: sname = 'GauntFactor'
!
   integer (kind=IntKind), intent(in) :: lmax_1
   integer (kind=IntKind), intent(in) :: lmax_2
   integer (kind=IntKind) :: lmax_3
   integer (kind=IntKind) :: kmax_1
   integer (kind=IntKind) :: kmax_2
   integer (kind=IntKind) :: jmax_3
   integer (kind=IntKind) :: kmax_3
   integer (kind=IntKind) :: lm3p1
   integer (kind=IntKind) :: ng
   integer (kind=IntKind) :: h3
   integer (kind=IntKind) :: kl1, l1, m1
   integer (kind=IntKind) :: kl2, l2, m2
   integer (kind=IntKind) :: l3, m3
   integer (kind=IntKind) :: j3count
! 
   real (kind=RealKind) :: endpts(2)
   real (kind=RealKind) :: ctemp
!
   real (kind=RealKind), allocatable :: plmg(:,:)
   real (kind=RealKind), allocatable :: tg(:)
   real (kind=RealKind), allocatable :: wg(:)
!
   data       endpts/-1.d0,+1.d0/
!
   lmax_3 = lmax_1+lmax_2
   kmax_1 = (lmax_1+1)**2
   kmax_2 = (lmax_2+1)**2
   kmax_3 = (lmax_3+1)**2
   jmax_3 = (lmax_3+1)*(lmax_3+2)/2
!
   if (min(lmax_1,lmax_2)+1 > MaxJ3) then
      call ErrorHandler(sname,'1st dim of cgnt is out of bound')
   endif
!
   lm3p1 = lmax_3+1
   allocate( plmg(jmax_3,lm3p1) )
   allocate( tg(2*lm3p1) )
   allocate( wg(2*lm3p1) )
!
!  ===================================================================
!  generate the gaussian integration pts and weights...............
!  -------------------------------------------------------------------
   call gaussq(1,2*lm3p1,0,endpts(1),endpts(2),tg,wg)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  generate the plm's at the gaussian nodes........................
!  ===================================================================
   do ng=1,lm3p1
!    -----------------------------------------------------------------
      call legendre(lmax_3,tg(ng),plmg(1:jmax_3,ng))
!    -----------------------------------------------------------------
   enddo
!
!  ===================================================================
!
!         4pi  _  m1 _     m2* _     m3 _
!  clll = int do Y  (o) * Y   (o) * Y  (o)
!          0      l1       l2        l3
!
!          L3
!       = C
!          L1,L2
!
!      l1 = 0, 1, ..., lmax_1
!      l2 = 0, 1, ..., lmax_2
!      l3 = 0, 1, ..., lmax_1+lmax_2
!
!  store clll in cgnt(j,L1,L2), j=1,2,...,nj3(L1,L2), L3=kj3(j,L1,L2)
!  ===================================================================
   do kl2=1,kmax_2
      l2=lofk(kl2)
      m2=mofk(kl2)
      do kl1=1,kmax_1
         l1=lofk(kl1)
         m1=mofk(kl1)
         m3=m2-m1
         h3=max(abs(m3),abs(l1-l2))
         j3count=0
         do l3=l1+l2,h3,-2
!           ----------------------------------------------------------
            ctemp=clll(l1,m1,l2,m2,l3,m3,wg,lm3p1,plmg,lmax_3)
!           ----------------------------------------------------------
            if(abs(ctemp).gt.tol) then
               if(j3count < MaxJ3) then
                  j3count=j3count+1
                  cgnt(j3count,kl1,kl2)=ctemp
                  kj3(j3count,kl1,kl2)=(l3+1)*(l3+1)-l3+m3
               else
                  call ErrorHandler(sname,'j3count > MaxJ3',j3count+1,MaxJ3)
               endif
            endif
         enddo
         nj3(kl1,kl2)=j3count
!        do j3count=1,nj3(kl1,kl2)
!           kl3=kj3(j3count,kl1,kl2)
!           write(6,'('' l1,m1,l2,m2,j3,cgnt ='',5i5,1d15.8)') &
!                        lofk(kl1),mofk(kl1),lofk(kl2),mofk(kl2), &
!                        j3count,cgnt(j3count,kl1,kl2)
!        enddo
      enddo
   enddo
!
   deallocate(plmg)
   deallocate(tg  )
   deallocate(wg  )
!
   end subroutine GauntFactor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function clll(l1,m1,l2,m2,l3,m3,wg,ngauss,plmg,lmax)
!  ===================================================================
!
   implicit  none
!
   integer (kind=IntKind), intent(in) :: l1
   integer (kind=IntKind), intent(in) :: l2
   integer (kind=IntKind), intent(in) :: l3
   integer (kind=IntKind), intent(in) :: m1
   integer (kind=IntKind), intent(in) :: m2
   integer (kind=IntKind), intent(in) :: m3
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: ngauss
!
   integer (kind=IntKind) :: ng
   integer (kind=IntKind) :: ifac1
   integer (kind=IntKind) :: jl1
   integer (kind=IntKind) :: ifac2
   integer (kind=IntKind) :: jl2
   integer (kind=IntKind) :: ifac3
   integer (kind=IntKind) :: jl3
!
   real (kind=RealKind), intent(in) :: wg(ngauss)
   real (kind=RealKind), intent(in) :: plmg((lmax+1)*(lmax+2)/2,ngauss)
   real (kind=RealKind) :: clll
!
!  *******************************************************************
!
!               4pi  _  m1 _     m2* _     m3 _
!        clll = int do Y  (o) * Y   (o) * Y  (o)
!                0      l1       l2        l3   
!
!                L3
!             = C
!                L1,L2
!
!  *******************************************************************
!
   if(l1.gt.lmax .or. l2.gt.lmax .or. l3.gt.lmax ) then
      call ErrorHandler('CLLL','bad parameters: l1,l2,l3,lmax',l1,l2,l3,lmax)
   endif
   clll=zero
   if(mod(l1+l2+l3,2).ne.0) then
      return
   else if(m1-m2+m3 .ne. 0) then
      return
   else if(l1+l2.lt.l3 .or. l2+l3.lt.l1 .or. l3+l1.lt.l2) then
      return
   else
!     ----------------------------------------------------------------
      call defac(l1, m1,ifac1,jl1)
      call defac(l2,-m2,ifac2,jl2)
      call defac(l3, m3,ifac3,jl3)
!     ----------------------------------------------------------------
      do ng=1,ngauss
         clll=clll+wg(ng)*plmg(jl1,ng)*plmg(jl2,ng)*plmg(jl3,ng)
      enddo
      if (abs(clll).lt.tol) then
         clll=zero
         return
      else
         clll=pi4*ifac1*clm(jl1)*m1m(m2)*ifac2*clm(jl2)*ifac3*clm(jl3)*clll
      endif
   endif
!  ===================================================================
   if (print_level.ge.2) then
      write(6,'(''l1,m1,l2,m2,l3,m3,clll = '',6i3,2x,d20.13)')        &
                  l1,m1,l2,m2,l3,m3,clll
   endif
!
   end function clll
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine defac(l,m,ifac,jl)
!  ===================================================================
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: l
   integer (kind=IntKind), intent(in) :: m
   integer (kind=IntKind), intent(out) :: ifac
   integer (kind=IntKind), intent(out) :: jl
!
!  ===================================================================
   if (m.ge.0) then
      ifac= 1
      jl = (l+1)*(l+2)/2-l+m
   else
      ifac= m1m(m)
      jl = (l+1)*(l+2)/2-l-m
   end if
!
   end subroutine defac
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isInitialized() result(y)
!  ===================================================================
!
   implicit   none
!
   logical :: y
!
   y = Initialized
!
   end function isInitialized
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pushGauntFactorsToAccel()
!  ===================================================================
   implicit none
!
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call initialize_gauntfactors(MaxJ3,kmax_sav)
   call push_gauntfactors(nj3,kj3,cgnt)
!  -------------------------------------------------------------------
#endif
!
   end subroutine pushGauntFactorsToAccel
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteGauntFactorsOnAccel()
!  ===================================================================
   implicit none
!
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call delete_gauntfactors()
!  -------------------------------------------------------------------
#endif
!
   end subroutine deleteGauntFactorsOnAccel
!  ===================================================================
end module GauntFactorsModule
