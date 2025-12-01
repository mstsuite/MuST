!  *******************************************************************
!  * Module for calculating KKR structure constant matrix            *
!  * Public functions:                                               *
!  *                                                                 *
!  *    call initRSpaceStrConst(lmax,istop,iprint)                   *
!  *    Purpose: initialize the module for calculating real space    *
!  *             KKR structure constant matrix                       *
!  *    Note:    GauntFactorsModule needs to be initialized first    *
!  *    Input:   lmax = largest possible l-cut off for gij matrix    *
!  *                    (integer)                                    *
!  *             istop = routine name to stop (character string)     *
!  *             iprint= print level (integer)                       *
!  *                                                                 *
!  *    getStrConstMatrix(kappa,rij,lmaxi,lmaxj)                     *
!  *    Purpose: calculate real space KKR structure constant matrix  *
!  *    Note:    need to call initRSpaceStrConst first               *
!  *    Input:   kappa  = sqrt of energy (complex)                   *
!  *             rij    = vector from atom i to j (real array(3))    *
!  *             lmaxi  = l-cut off for row index of gij (integer)   *
!  *             lmaxj  = l-cut off for column index of gij (integer)*
!  *             Both lmaxi and lmaxj can not be larger than lmax,   *
!  !             which is previously defined when call initStrConst. *
!  *    Result:  gij    = real space KKR structure constant matrix   *
!  *                      (complex array(:,:)), where                *
!  *                      1st dim = kmaxi=(lmaxi+1)**2               *
!  *                      2nd dim = kmaxj=(lmaxj+1)**2               *
!  *                                                                 *
!  *    call delStrConstMatrix()                                     *
!  *    Purpose: deallocate the structure constanct matrix gij       *
!  *                                                                 *
!  *    call endRspaceStrConst()                                     *
!  *    Purpose: deallocate the internal allocated arrays and clear  *
!  *             the storage                                         *
!  *    Usage:   should be used whenever one or more atoms change    *
!  *             their position or before the end of the program     *
!  *******************************************************************
module RSpaceStrConstModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use KindParamModule, only : CmplxKind
!
   use MathParamModule, only : zero
   use MathParamModule, only : half
   use MathParamModule, only : one
   use MathParamModule, only : two
   use MathParamModule, only : ten2m6
   use MathParamModule, only : ten2m8
   use MathParamModule, only : pi
   use MathParamModule, only : pi2
   use MathParamModule, only : pi4
   use MathParamModule, only : czero
   use MathParamModule, only : cone
   use MathParamModule, only : sqrtm1
!
   use ErrorHandlerModule, only : StopHandler
   use ErrorHandlerModule, only : ErrorHandler
   use ErrorHandlerModule, only : WarningHandler
!
   use IntegerFactorsModule, only : lofk, mofk
!
   use SphericalHarmonicsModule, only : calYlmConjg
!
public :: initRSpaceStrConst, &
          endRSpaceStrConst,  &
!         delStrConstMatrix,  &
          getStrConstMatrix
!
private
!
   logical :: Initialized = .false.
!
   character (len=40) :: stop_routine
!
   integer (kind=IntKind) :: lmax_max
   integer (kind=IntKind) :: kmax_max
   integer (kind=IntKind) :: lmax_dlm
   integer (kind=IntKind) :: kmax_dlm
   integer (kind=IntKind) :: print_level
!
   complex (kind=CmplxKind), allocatable :: illp(:,:)
   complex (kind=CmplxKind), allocatable :: ilp1(:)
!
   complex (kind=CmplxKind), allocatable :: ylmcc(:)
   complex (kind=CmplxKind), allocatable :: dlm(:)
   complex (kind=CmplxKind), allocatable :: hfn(:)
   complex (kind=CmplxKind), allocatable, target :: gijp(:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genFactors(lmax)
!  ===================================================================
   implicit   none
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l, m, kl, i, lp, klp, mp
!
   ilp1(0)=sqrtm1
   i=lmax
   do l=1,i
      ilp1(l)=ilp1(l-1)*sqrtm1
   enddo
   do l=0,lmax
      do m = -l,l
         kl = l*(l+1)+1+m
         do lp=0,lmax
            do mp = -lp,lp
               klp = lp*(lp+1)+1+mp
               illp(kl,klp)=ilp1(l)/ilp1(lp)
            enddo
         enddo
      enddo
   enddo
!
!  do l=0,lmax
!     do lp=0,lmax
!        illp(kl,klp)=ilp1(l)/ilp1(lp)
!     enddo
! enddo
!
   end subroutine genFactors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initRSpaceStrConst(lmax,istop,iprint)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
!
   implicit   none
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: lmax_max2, kmax_max2
!
   Initialized=.true.
   stop_routine=istop
   print_level=iprint
!
   lmax_max=lmax
   kmax_max=(lmax+1)*(lmax+1) 
   lmax_max2=2*lmax
   kmax_max2=(lmax_max2+1)*(lmax_max2+1)
   call initIntegerFactors(lmax_max2)
   allocate( illp(1:kmax_max2,1:kmax_max2), ilp1(0:lmax_max2) )
   allocate( gijp(kmax_max*kmax_max) )
   allocate( ylmcc(kmax_max2) )
   allocate( dlm(kmax_max2) )
   allocate( hfn(0:lmax_max2) )
!  -------------------------------------------------------------------
   call genFactors(lmax_max2)
!  -------------------------------------------------------------------
!
   lmax_dlm=0
   kmax_dlm=1
!
   end subroutine initRSpaceStrConst
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endRSpaceStrConst()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
   implicit   none
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endRSpaceStrConst',                     &
                        'need to initialize RSpaceStrConstModule first')
!     ----------------------------------------------------------------
   endif
!
   deallocate( gijp, ylmcc, dlm, hfn )
!
   deallocate( illp, ilp1 )
!  -------------------------------------------------------------------
   call endIntegerFactors()
!  -------------------------------------------------------------------
!
   Initialized=.false.
!
   end subroutine endRSpaceStrConst
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStrConstMatrix(kappa,rij,lmaxi,lmaxj,iprint) result(gij)
!  ===================================================================
!
!  *******************************************************************
!  *                                                                 *
!  *  calculate KKR real space structure constant matrix: gij        *
!  *                                                                 *
!  *  Note:                                                          *
!  *        * the returned value strcon_matrix is in the same units  *
!  *          as the input quantities.                               *
!  *                                                                 *
!  *        * L = 1,2,...,kmaxi, L1=1,2,...,kmaxj                    *
!  *                                                                 *
!  *                                                                 *
!  *   ij         l+1                                m ->  *         *
!  *  D  [E]  = -i    * sqrt(E) * h [sqrt(E)*R  ] * Y [R  ]          *
!  *   L                           l          ij     l  ij           *
!  *                                                                 *
!  *                                                                 *
!  *              -l+1                               m  ->  *        *
!  *          = -i    * sqrt(E) * h [sqrt(E)*R  ] * Y [-R  ]         *
!  *                               l          ij     l   ij          *
!  *                                                                 *
!  *                                                                 *
!  *   m       m     [ (2*l+1)*(l-|m|)!]  m                          *
!  *  Y  = (-1) *sqrt[-----------------]*P [cos(theta)]*exp(i*m*phi) *
!  *   l             [   4*pi*(l+|m|)! ]  l                          *
!  *                                                                 *
!  *  for m>=0                                                       *
!  *                                                                 *
!  *   m                  m  -m           *                          *
!  *  Y [theta,phi] = (-1)  Y  [theta,phi]      for m<0              *
!  *   l                     l                                       *
!  *                                                                 *
!  *  ->    ->   ->                                                  *
!  *  R   = R  - R  ==> [theta,phi]                                  *
!  *   ij    j    i                                                  *
!  *                                                                 *
!  *  cos(theta)=Rij(3)/sqrt(Rij(1)**2+Rij(2)**2+Rij(3)**2)          *
!  *                                                                 *
!  *  cos(phi)  =Rij(1)/sqrt(Rij(1)**2+Rij(2)**2)                    *
!  *                                                                 *
!  *             m     m                                             *
!  *  Indexing: P  :  P(l*(l+1)/2+m+1)  only +m are calculated       *
!  *             l     l                                             *
!  *                                                                 *
!  *             m     m                                             *
!  *  Indexing: C  :  C(l*(l+1)/2+m+1)  only +m are calculated       *
!  *             l     l                                             *
!  *                                                                 *
!  *             m     m                                             *
!  *  Indexing: D  :  D(l*(l+1)+m+1)    all   m are calculated       *
!  *             l     l                                             *
!  *                                                                 *
!  *  Now for the real space structure contant we have :             *
!  *                                                                 *
!  *   ij             l-l1         L2     ij                         *
!  *  G   (E) = 4*pi*i    * SUM { C    * D  (E) }                    *
!  *   L,L1                  L2    L1,L   L2                         *
!  *                                                                 *
!  *******************************************************************
!
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
use LegendreModule, only : legendre
use SphericalHarmonicsModule, only : getClm
!
   implicit   none
!
   character (len=23), parameter :: sname='getRSpaceStrConstMatrix'
!
   integer (kind=IntKind), intent(in) :: lmaxi
   integer (kind=IntKind), intent(in) :: lmaxj
   integer (kind=IntKind), intent(in), optional :: iprint
!
   integer (kind=IntKind) :: kmaxi
   integer (kind=IntKind) :: kmaxj
   integer (kind=IntKind) :: kl
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: klp
   integer (kind=IntKind) :: lp
   integer (kind=IntKind) :: nnj3
#ifdef DEBUG
   integer (kind=IntKind) :: mp
#endif
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
   integer (kind=IntKind), pointer :: pkj3(:)
!
   real (kind=RealKind), intent(in) :: rij(3)
!
   real (kind=RealKind) :: rmag
   real (kind=RealKind), pointer :: cgnt(:,:,:)
   real (kind=RealKind), pointer :: pcgnt(:)
real (kind=RealKind), pointer :: clm(:)
real (kind=RealKind) :: plm(138), cos_theta
integer (kind=IntKind) :: jmax_dlm
!
   complex (kind=CmplxKind), intent(in) :: kappa
!
   complex (kind=CmplxKind) :: z
   complex (kind=CmplxKind) :: fac, gij_llp
   complex (kind=CmplxKind), pointer :: gij(:,:)
!
   rmag=sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
!
#ifdef DEBUG
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'need to call initRSpaceStrConst first')
!     ----------------------------------------------------------------
   else if(abs(kappa) < ten2m8) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'energy = 0.0')
!     ----------------------------------------------------------------
   else if(rmag < ten2m6) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'rmag = 0',rmag)
!     ----------------------------------------------------------------
   else if(lmaxi > lmax_max) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'lmaxi > lmax_max',lmaxi,lmax_max)
!     ----------------------------------------------------------------
   else if(lmaxj > lmax_max) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'lmaxj > lmax_max',lmaxj,lmax_max)
!     ----------------------------------------------------------------
   endif
#endif
!
   kmaxi=(lmaxi+1)*(lmaxi+1)
   kmaxj=(lmaxj+1)*(lmaxj+1)
!
   gijp = CZERO
   gij => aliasArray2_c(gijp,kmaxi,kmaxj)
!
   lmax_dlm=lmaxi+lmaxj
   kmax_dlm=(lmax_dlm+1)*(lmax_dlm+1)
!  ====================================================================
!  calculate the hankel function. [dangerous code if z is close to 0]
!  hankel function is hfn*exp(i*z)/z
!  ====================================================================
   z=kappa*rmag
   hfn(0)=-sqrtm1
   if ( lmax_dlm > 0 ) then
      hfn(1)=-(cone+sqrtm1/z)
      do l=2,lmax_dlm
         hfn(l)=(2*l-1)*hfn(l-1)/z - hfn(l-2)
      enddo
   endif
!
!  ====================================================================
!  calculate Y(l,m)...................................................
!  --------------------------------------------------------------------
   call calYlmConjg(rij,lmax_dlm,ylmcc)
!  --------------------------------------------------------------------
   if (present(iprint)) then
      if (iprint > 0) then
         write(6,'(a,3f12.5)')'rij = ',rij(1:3)
clm => getClm(lmax_dlm)
cos_theta=rij(3)/rmag
plm = ZERO
call legendre(lmax_dlm,cos_theta,plm)
jmax_dlm=(lmax_dlm+1)*(lmax_dlm+2)/2
         write(6,'(a,1i5,2d15.8)')'jl, Plm, Clm ',jmax_dlm-lmax_dlm+1, &
                                   plm(jmax_dlm-lmax_dlm+1),clm(jmax_dlm-lmax_dlm+1)
         write(6,'(a,1i5,2d15.8)')'kl, Ylmcc = ',kmax_dlm-lmax_dlm+1,  &
                                  ylmcc(kmax_dlm-lmax_dlm+1)
         write(6,'(a,1i5,2d15.8)')'jl, Plm, Clm ',jmax_dlm-1, &
                                   plm(jmax_dlm-1),clm(jmax_dlm-1)
         write(6,'(a,1i5,2d15.8)')'kl, Ylmcc = ',kmax_dlm-1,  &
                                  ylmcc(kmax_dlm-1)
         write(6,'(a,1i5,2d15.8)')'jl, Plm, Clm ',jmax_dlm, &
                                   plm(jmax_dlm),clm(jmax_dlm)
         write(6,'(a,1i5,2d15.8)')'kl, Ylmcc = ',kmax_dlm,  &
                                  ylmcc(kmax_dlm)
      endif
   endif
!
!  ===================================================================
!  generate the KKR real space lattice structure matrix for the energy
!  and store result in gij
!
!             l+1
!     fac = -i   *h (k*R  )*sqrt(E)
!                  l    ij
!  ===================================================================
   z = exp(sqrtm1*z)/rmag
   do kl = 1,kmax_dlm
      l =lofk(kl)
!     fac = -hfn(l)*z*ilp1(l)   ! 08/28/14
      fac =  hfn(l)*z/ilp1(l)
      dlm(kl) = fac*ylmcc(kl)
   enddo
   if (present(iprint)) then
      if (iprint > 0) then
         write(6,'(a,3f12.5)')'rij(1:3) = ',rij(1:3)
         write(6,'(a,2d15.8)')'kappa = ',kappa
         write(6,'(a,2d15.8)')'hfn   = ',hfn(4)
         write(6,'(a,2d15.8)')'hfn   = ',hfn(7)
         write(6,'(a,2d15.8)')'hfn   = ',hfn(8)
         write(6,'(a,2d15.8)')'ylmcc(37) = ',ylmcc(37)
         write(6,'(a,2d15.8)')'ylmcc(66) = ',ylmcc(66)
         write(6,'(a,2d15.8)')'ylmcc(81) = ',ylmcc(81)
         write(6,'(a,2d15.8)')'dlm(37)   = ',dlm(37)
         write(6,'(a,2d15.8)')'dlm(66)   = ',dlm(66)
         write(6,'(a,2d15.8)')'dlm(81)   = ',dlm(81)
      endif
   endif
!
!  ===================================================================
!  calculate g(R_ij)...................................................
!  ===================================================================
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
!  ===================================================================
!  loop over klp......................................................
!  ===================================================================
   do klp=1,kmaxj
!      lp=lofk(klp)
!     ================================================================
!     loop over kl....................................................
!     ================================================================
      do kl=1,kmaxi
!         l=lofk(kl)
!        =============================================================
!                      l-lp
!        illp(l,lp) = i
!
!        perform sum over j with gaunt # ..............................
!        =============================================================
         nnj3 = nj3(klp,kl)
         pkj3 => kj3(1:nnj3,klp,kl)
         pcgnt=> cgnt(1:nnj3,klp,kl)
         gij_llp = CZERO
         do j = 1,nnj3
            gij_llp = gij_llp+pcgnt(j)*dlm(pkj3(j))
   if (present(iprint)) then
      if (iprint > 0 .and. kl == 25 .and. klp == 17) then
         write(6,'(a,i5,2x,2d15.8)')'pkj3, dlm = ',pkj3(j),dlm(pkj3(j))
         write(6,'(a,d15.8,2x,2d15.8)')'pcgnt, gllp = ',pcgnt(1), gij_llp
      endif
   endif
         enddo
         gij(kl,klp)=pi4*illp(kl,klp)*gij_llp
      enddo
   enddo
   if (present(iprint)) then
      if (iprint > 0) then
         write(6,'(a,2i5)')'nj3(25,17),nj3(17,25) = ',nj3(25,17),nj3(17,25)
         write(6,'(a,2d15.8)')'gij(25,17) = ', gij(25,17)
         write(6,'(a,2d15.8)')'gij(17,25) = ', gij(17,25)
      endif
   endif
!
#ifdef DEBUG
   if(print_level >= 1) then
!     ================================================================
!     loop over lp,mp.................................................
!     ================================================================
      do klp=1,kmaxj
         lp=lofk(klp)
         mp=mofk(klp)
!        =============================================================
!        loop over l,m................................................
!        =============================================================
         do kl=1,kmaxi
            l=lofk(kl)
            m=mofk(kl)
            write(6,'('' klp,kl,lp,mp,l,m, gij:'',6i3,2d12.5)')       &
                         klp,kl,lp,mp,l,m,gij(kl,klp)
         enddo
      enddo
   endif
#endif
!  -------------------------------------------------------------------
   nullify( nj3, kj3, cgnt )
!  -------------------------------------------------------------------
   end function getStrConstMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray_c(size1,size2,array)              result(parray)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: size1,size2
   complex(kind=CmplxKind), target :: array(size1,size2)
!
   complex(kind=CmplxKind), pointer :: parray(:,:)
!
   parray => array
!
   end function aliasArray_c
!  ===================================================================
end module RSpaceStrConstModule
