!  *******************************************************************
!  *  Title: Step Function Module                                    *
!  *                                                                 *
!  *  Version History:                                               *
!  *         v2.0     October 21, 2005        Y. Wang                *
!  *         v1.0     October  7, 2003        Y. Wang                *
!  *                  First release.                                 *
!  *                                                                 *
!  *  Interfaces and Usage:                                          *
!  *       # initStepFunction(na,lmax,ngr,ngt,istop,iprint)          *
!  *             Purpose: initialize the StepFunction module.        *
!  *             Input:: na: no. of polyhedra to be considered.      *
!  *                     lmax(1:na): lmax cutoff of step function    *
!  *                                 for each polyhedron.            *
!  *                     ngr(1:na): no. of radial Gaussian points    *
!  *                                within each critical interval    *
!  *                                for each polyhedron. This is     *
!  *                                only an estimate.                *
!  *                     ngt(1:na): no. of cos(Theta) Gaussian       *
!  *                                points within each cos(Theta)    *
!  *                                interval for each polyhedron.    *
!  *                     istop: the stop routine name.               *
!  *                     iprint: the print out instruction.          *
!  *       # initStepFunction(na,lmax,istop,iprint)                  *
!  *             Purpose:: initialize the StepFunction module.       *
!  *             Input:: na: no. of polyhedra to be considered.      *
!  *                     lmax(1:na): lmax cutoff of step function    *
!  *                                 for each polyhedron.            *
!  *                     istop: the stop routine name.               *
!  *                     iprint: the print out instruction.          *
!  *             Note:: if this routine is used for initialization,  *
!  *                    the number of Gaussian points are taking     *
!  *                    from the internal default parameters.        *
!  *       # endStepFunction()                                       *
!  *             Purpose:: clean memory spaces.                      *
!  *       # getNumGaussRs(n)                                        *
!  *             Purpose:: given a polyhedron index, returns the     *
!  *                       number of radial Gaussian points.         *
!  *             Input:: n: the polyhedron index.                    *
!  *             Return:: the number of radial Gaussian points.      *
!  *       # getGaussR(n)                                            *
!  *             Purpose:: given a polyhedron index, returns the     *
!  *                       radial Gaussian points.                   *
!  *             Input:: n: the polyhedron index.                    *
!  *             Return:: a pointer to the radial Gaussian points.   *
!  *       # getGaussRWeight(n)                                      *
!  *             Purpose:: given a polyhedron index, returns the     *
!  *                       weight of radial Gaussian points.         *
!  *             Input:: n: the polyhedron index.                    *
!  *             Return:: a pointer to the weight of radial Gausian  *
!  *                      points.                                    *
!  *       # getNumCriticalRs(n)                                     *
!  *             Purpose:: given a polyhedron index, returns the     *
!  *                       number of radial critical points.         *
!  *             Input:: n: the polyhedron index.                    *
!  *             Return:: the number of radial critical points.      *
!  *       # getCriticalR(n)                                         *
!  *             Purpose:: given a polyhedron index, returns the     *
!  *                       radial critical points.                   *
!  *             Input:: n: the polyhedron index.                    *
!  *             Return:: a pointer to the radial critical points.   *
!  *       # getStepFunction(n)                                      *
!  *             Purpose:: given a polyhedron index, returns the     *
!  *                       step function of the radial grid points.  *
!  *             Input:: n: the polyhedron index.                    *
!  *             Return:: a pointer to the step function of the      *
!  *                      radial Gaussian points for l upto lmax.    *
!  *       # getStepFunction(n,jl)                                   *
!  *             Purpose:: given a polyhedron index and the {l,m}    *
!  *                       (m >= 0) index, returns the step function *
!  *                       of the radial grid points.                *
!  *             Input:: n: the polyhedron index.                    *
!  *                     jl: the {l, m, m>=0} index.                 *
!  *             Return:: a pointer to the step function of the      *
!  *                      radial Gaussian points of jl.              *
!  *       # getStepFunction(n,jl,r)                                 *
!  *             Purpose:: given a polyhedron index, the {l,m, m>=0} *
!  *                       index, and r, returns the step function.  *
!  *             Input:: n: the polyhedron index.                    *
!  *                     jl: the {l, m, m>=0} index.                 *
!  *                     r: a radial value.                          *
!  *             Return:: the step function of r and jl.             *
!  *       # interpolateStepFunction(n,nr,r,kmax1,kmax2,s)           *
!  *             Purpose:: given a polyhedron index, kmax1, kmax2,   *
!  *                       nr, and r-mesh r(nr), returns the inter-  *
!  *                       polated step function s(nr,kmax1,kmax2).  *
!  *             Input:: n: the polyhedron index.                    *
!  *                     nr: the number if r-mesh.                   *
!  *                     r: a radial mesh.                           *
!  *                     kmax1, kmax2: {l,m} index dimension         *
!  *             Return:: the interpolated step function.            *
!  *       # getVolumeIntegration(n,nr,r,k,f,truncate)               *
!  *             Purpose:: given a polyhedron index and a real       *
!  *                       function data, calculates the volume      *
!  *                       integration as follows:                   *
!  *                          v = int_{cell n} d^3r r^k*f(r)         *
!  *             Input:: n: the polyhedron index.                    *
!  *                     nr: number of radial points.                *
!  *                     r: r-mesh array, on which the function is   *
!  *                        defined.                                 *
!  *                     k: the input function contains a factor     *
!  *                        r^{k}.                                   *
!  *                     f: the function multiplied by r^k for the   *
!  *                          volume integration. Note: the function *
!  *                          is defined on r-mesh.                  *
!  *                     truncate: optional, if it is present and is *
!  *                               false, the function will not be   *
!  *                               truncated.                        *
!  *             Return:: the volume integration.                    *
!  *       # getVolumeIntegration(n,nr,r,kmax,k,frk,truncate)        *
!  *             Purpose:: given a polyhedron index and a function   *
!  *                       data, calculates the volume integration.  *
!  *             Input:: n: the polyhedron index.                    *
!  *                     nr: number of radial points.                *
!  *                     r: r-mesh array, on which the function is   *
!  *                        defined.                                 *
!  *                     kmax: the upper limit of {l,m} index of the *
!  *                           the input function.                   *
!  *                     k: the input function contains a factor     *
!  *                        r^{k}.                                   *
!  *                     frk: the function multiplied by r^k for the *
!  *                          volume integration. Note: the function *
!  *                          is defined on r-mesh and is expanded   *
!  *                          on complex spherical harmonics.        *
!  *                     truncate: optional, if it is present and is *
!  *                               false, the function will not be   *
!  *                               truncated.                        *
!  *             Return:: a complex value of the volume integration. *
!  *       # getVolumeIntegration(n,nr,r,kmaxf,kf,frk,kmaxg,kg,grk,t)*
!  *             Purpose:: given a polyhedron index and two function *
!  *                       data, calculates the volume integration   *
!  *                       of the two function product.              *
!  *             Input:: n: the polyhedron index.                    *
!  *                     nr: number of radial points.                *
!  *                     r: r-mesh array, on which the function is   *
!  *                        defined.                                 *
!  *                     kmaxf: the upper limit of {l,m} index of    *
!  *                            the input function frk.              *
!  *                     kf: the input function frk contains a       *
!  *                         factor r^{kf}.                          *
!  *                     frk: f function multiplied by r^kf. It is   *
!  *                          defined on r-mesh and is expanded on   *
!  *                          complex spherical harmonics.           *
!  *                     kmaxg: the upper limit of {l,m} index of    *
!  *                            the input function grk.              *
!  *                     kg: the input function grk contains a       *
!  *                         factor r^{kg}.                          *
!  *                     grk: g function multiplied by r^kg. It is   *
!  *                          defined on r-mesh and is expanded on   *
!  *                          complex spherical harmonics.           *
!  *                     t:  optional, if it is present and is false,*
!  *                         the function will not be truncated      *
!  *       # testStepFunction(n)                                     *
!  *             Purpose:: given a polyhedron index, test the        *
!  *                       step function.                            *
!  *             Input:: n: the polyhedron index.                    *
!  *             Note:: What has been tested here is the integration *
!  *                    of the step function over the polyhedron.    *
!  *                    The result should be the polyhedron volume.  *
!  *                    Other test cases can be implemented later.   *
!  *       # printStepFunction(n)                                    *
!  *             Purpose:: given a polyhedron index, print out the   *
!  *                       step function information.                *
!  *             Input:: n: the polyhedron index.                    *
!  *                                                                 *
!  *******************************************************************
module StepFunctionModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   use MathParamModule, only : ZERO, THIRD, HALF, ONE, TWO, THREE
   use MathParamModule, only : PI, PI2, PI4, Y0inv
   use MathParamModule, only : TEN2m3, TEN2m4, TEN2m6, TEN2m8, TEN2m10, TEN2m12, TEN
   use MathParamModule, only : CZERO, CONE, SQRTm1
   use IntegerFactorsModule, only : lofk, mofk, jofk, lofj, mofj, m1m
   use InterpolationModule, only : FitInterp
!
public :: initStepFunction,     &
          endStepFunction,      &
          getNumGaussRs,        &
          getGaussR,            &
          getGaussRWeight,      &
          getNumCriticalRs,     &
          getCriticalR,         &
          getStepFunction,      &
          getVolumeIntegration, &
          testStepFunction,     &
          interpolateStepFunction, &
          getSigmaLL,           &
          truncate,             &
          getRadialStepFunction,&
          printStepFunction
!
   interface initStepFunction
      module procedure initStepFunction1, initStepFunction2
   end interface
!
   interface getStepFunction
      module procedure getStepFunction1, getStepFunction2, getStepFunction3
   end interface
!
   interface getVolumeIntegration
      module procedure getVolumeIntegration0
      module procedure getVolumeIntegration1, getVolumeIntegration1_jl, &
                       getVolumeIntegration1r_jl
      module procedure getVolumeIntegration2
   end interface
!
   interface truncate
      module procedure truncate_all, truncate_jl, truncate_symm
   end interface
!
   interface testStepFunction
      module procedure testStepFunction_i, testStepFunction_t
   end interface
!
private
   character (len=50) :: stop_routine
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: NumPolyhedra
!
!  *******************************************************************
!  * Note:                                                           *
!  *                                                                 *
!  *    sigma_L  = sigma(r,L)                                        *
!  *             = int_{4pi} ds * Y^{*}_{L}(vec{r}) * sigma(vec{r})  *
!  *                                                                 *
!  *    sigma_LL = sigma(r,L,L1)                                     *
!  *             = int_{4pi} ds * Y_{L}(vec{r}) * sigma(vec{r})      *
!  *                            * Y^{*}_{L1}(vec{r})                 *
!  *******************************************************************
!
   type StepFunctionStruct
      integer (kind=IntKind) :: lmax_step
      integer (kind=IntKind) :: jmax_step
      integer (kind=IntKind) :: kmax_step
      integer (kind=IntKind) :: lmax_half
      integer (kind=IntKind) :: kmax_half
      integer (kind=IntKind) :: NumGaussRs
      integer (kind=IntKind) :: NumGaussCosThetas
      integer (kind=IntKind) :: NumCriticalRs
      integer (kind=IntKind) :: np_index(2)
      integer (kind=IntKind), pointer :: SigmaCompFlag(:)
!
      real (kind=RealKind) :: rmin
      real (kind=RealKind) :: rmax
      real (kind=RealKind) :: rp_index(2)
      real (kind=RealKind), pointer :: CriticalR(:)
      real (kind=RealKind), pointer :: GaussCosTheta(:)
      real (kind=RealKind), pointer :: GaussCosThetaWeight(:)
      real (kind=RealKind), pointer :: GaussR(:)
      real (kind=RealKind), pointer :: GaussRWeight(:)
!
      complex (kind=CmplxKind), pointer :: sigma_L(:,:)
      complex (kind=CmplxKind), pointer :: sigma_LL(:,:,:)
   end type StepFunctionStruct
!
   type (StepFunctionStruct), allocatable, target :: StepFunction(:) 
!
   integer (kind=IntKind), allocatable :: fflag(:), gflag(:)
   integer (kind=IntKind) :: kmax_ff, kmax_gf, max_nr
   integer (kind=IntKind) :: current_polyhedron
!
   integer (kind=IntKind), parameter :: mode = 0
   integer (kind=IntKind), parameter :: NumGaussRs = 80
   integer (kind=IntKind), parameter :: NumGaussCosThetas = 60
!
   real (kind=RealKind) :: current_r
   real (kind=RealKind), pointer :: clm(:)
!
   real (kind=RealKind), allocatable, target :: sqrt_r(:), rr0(:)
!
   real (kind=RealKind), parameter :: sigma_tol = TEN2m8
!
   real (kind=RealKind), allocatable :: rfg(:), res(:)
   complex (kind=CmplxKind), allocatable :: cfg(:), ces(:)
   complex (kind=CmplxKind), allocatable :: wylm(:), ga(:)
   complex (kind=CmplxKind), allocatable :: fr2(:)
!
   logical :: Initialized = .false.
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initStepFunction1(na,lmax_in,lmax,ngr,ngt,istop,iprint)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
!
   use SphericalHarmonicsModule, only : getClm
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: lmax_in
   integer (kind=IntKind), intent(in) :: lmax(na)
   integer (kind=IntKind), intent(in) :: ngr(na)
   integer (kind=IntKind), intent(in) :: ngt(na)
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: lmax_max, jmax_max
   integer (kind=IntKind) :: i
!
   stop_routine = istop
   print_level = iprint
!
   if (na < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initStepFunction','invalid number of polyhedra',na)
!     ----------------------------------------------------------------
   endif
!
   NumPolyhedra = na
   allocate( StepFunction(na) )
!
   lmax_max = 0
   do i = 1, na
      if (lmax(i) < 0) then
!        -------------------------------------------------------------
         call ErrorHandler('initStepFunction','invalid lmax',lmax(i))
!        -------------------------------------------------------------
      else
!        lmax_max = max(lmax_max, 2*lmax(i))
         lmax_max = max(lmax_max, lmax(i))
      endif
   enddo
!
!  -------------------------------------------------------------------
!  call initIntegerFactors(2*max(lmax_max,lmax_in))
   call initIntegerFactors(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  allocate the working space for step function calculation
!  ===================================================================
   jmax_max = (lmax_max+1)*(lmax_max+2)/2
   allocate( wylm(jmax_max) )
!  -------------------------------------------------------------------
   clm => getClm(lmax_max)
!  -------------------------------------------------------------------
!
   Initialized = .true.
!
   do i = 1, na
      if (ngt(i) < 1) then
!        -------------------------------------------------------------
         call ErrorHandler('initStepFunction',                        &
                           'invalid number of theta gaussian points',ngt(i))
!        -------------------------------------------------------------
      else if (ngr(i) < 1) then
!        -------------------------------------------------------------
         call ErrorHandler('initStepFunction',                        &
                           'invalid number of radial gaussian points',ngr(i))
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
      call setupStepFunction(i,lmax(i),ngr(i),ngt(i))
!     ----------------------------------------------------------------
   enddo
!
   kmax_ff = 0; kmax_gf = 0; max_nr = 0
!
   end subroutine initStepFunction1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initStepFunction2(na,lmax_in,lmax,istop,iprint)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
!
   use SphericalHarmonicsModule, only : getClm
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: lmax_in
   integer (kind=IntKind), intent(in) :: lmax(na)
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: lmax_max, jmax_max
   integer (kind=IntKind) :: i
!
   stop_routine = istop
   print_level = iprint
!
   if (na < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initStepFunction','invalid number of polyhedra',na)
!     ----------------------------------------------------------------
   endif
!
   NumPolyhedra = na
   allocate( StepFunction(na) )
!
   lmax_max = 0
   do i = 1, na
      if (lmax(i) < 0) then
!        -------------------------------------------------------------
         call ErrorHandler('initStepFunction','invalid lmax',lmax(i))
!        -------------------------------------------------------------
      else
         lmax_max = max(lmax_max, lmax(i))
      endif
   enddo
!
!  -------------------------------------------------------------------
!  call initIntegerFactors(2*max(lmax_max,lmax_in))
   call initIntegerFactors(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  allocate the working space for step function calculation
!  ===================================================================
   jmax_max = (lmax_max+1)*(lmax_max+2)/2
   allocate( wylm(jmax_max) )
!  -------------------------------------------------------------------
   clm => getClm(lmax_max)
!  -------------------------------------------------------------------
!
   Initialized = .true.
!
   do i = 1, na
!     ----------------------------------------------------------------
      call setupStepFunction(i,lmax(i),NumGaussRs,NumGaussCosThetas)
!     ----------------------------------------------------------------
   enddo
!
   kmax_ff = 0; kmax_gf = 0; max_nr = 0
!
   end subroutine initStepFunction2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupStepFunction(poly,lmax,ngr,ngt)
!  ===================================================================
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
   use PolyhedraModule, only : getNumPlanes, getPlane
   use PolyhedraModule, only : getNumCorners, getCorner
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: ngr
   integer (kind=IntKind), intent(in) :: ngt
   integer (kind=IntKind), allocatable :: ng(:)
   integer (kind=IntKind) :: ic, ig, n, j, np, ncorn, nbnd, nmax
   integer (kind=IntKind) :: lmax_step, jmax_step, kmax_step
   integer (kind=IntKind) :: lmax_half, kmax_half
   integer (kind=IntKind) :: m, jl, kl, kl1, kl2
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
!
   real (kind=RealKind) :: r1, r2, fac
   real (kind=RealKind) , allocatable:: rg(:,:), wg(:,:)
   real (kind=RealKind), pointer :: xpt(:,:), corn(:,:)
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!
   complex (kind=CmplxKind) :: ssum
!
   current_polyhedron = 0
!
!  lmax_half = lmax/2
   lmax_half = lmax
   kmax_half = (lmax_half+1)*(lmax_half+1)
   lmax_step = lmax
   jmax_step = (lmax_step+1)*(lmax_step+2)/2
   kmax_step = (lmax_step+1)*(lmax_step+1)
!
   StepFunction(poly)%lmax_step = lmax_step
   StepFunction(poly)%jmax_step = jmax_step
   StepFunction(poly)%kmax_step = kmax_step
   StepFunction(poly)%lmax_half = lmax_half
   StepFunction(poly)%kmax_half = kmax_half
   StepFunction(poly)%rmin = getInscrSphRadius(poly)
   StepFunction(poly)%rmax = getOutscrSphRadius(poly)
   if ( StepFunction(poly)%rmin < TEN2m6) then
!     ----------------------------------------------------------------
      call ErrorHandler('setupStepFunction',                          &
                        'invalid inscribed sphere radius',            &
                        StepFunction(poly)%rmin)
!     ----------------------------------------------------------------
   else if ( StepFunction(poly)%rmax/StepFunction(poly)%rmin < 1.01) then
!     ----------------------------------------------------------------
      call ErrorHandler('setupStepFunction',                          &
                        'invalid outscribed sphere radius',           &
                        StepFunction(poly)%rmax)
!     ----------------------------------------------------------------
   endif
!
   StepFunction(poly)%NumGaussCosThetas = ngt
   allocate( StepFunction(poly)%GaussCosTheta(1:ngt) )
   allocate( StepFunction(poly)%GaussCosThetaWeight(1:ngt) )
!  -------------------------------------------------------------------
   call genGaussianPoints( ngt,-ONE,ONE,                              &
                           StepFunction(poly)%GaussCosTheta,          &
                           StepFunction(poly)%GaussCosThetaWeight )
!  -------------------------------------------------------------------
!
!  ===================================================================
!  determine the critical points of the atomic cell.
!  -------------------------------------------------------------------
   call computeCriticalPoints(poly)
!  -------------------------------------------------------------------
!
   fac = StepFunction(poly)%rmin*(sqrt(TWO)-ONE)
!
!  ===================================================================
!  The number Gaussian points in each interval is taken proportionally.
!  It is assumed that for a square with side length = 2, the number of
!  Gaussian points between the inscribed circle and the circumscribed
!  circle is ngr.
!  ===================================================================
   n = StepFunction(poly)%NumCriticalRs
   allocate( ng(n-1) )
   nmax = 0
   do ic = 2,n
      ng(ic-1) = ceiling(ngr*(StepFunction(poly)%CriticalR(ic)-       &
                              StepFunction(poly)%CriticalR(ic-1))/fac)
!       ng(ic-1)=20
      nmax = max(nmax, ng(ic-1))
   enddo
   allocate( rg(nmax,n-1), wg(nmax,n-1) )
   StepFunction(poly)%NumGaussRs = 0
   do ic = 2,n
      r1 = StepFunction(poly)%CriticalR(ic-1)
      r2 = StepFunction(poly)%CriticalR(ic)
      StepFunction(poly)%NumGaussRs = StepFunction(poly)%NumGaussRs + ng(ic-1)
!     ----------------------------------------------------------------
      call genGaussianPoints(ng(ic-1),r1,r2,rg(1,ic-1),wg(1,ic-1))
!     ----------------------------------------------------------------
   enddo
   n = StepFunction(poly)%NumGaussRs
   allocate( StepFunction(poly)%GaussR(1:n) )
   allocate( StepFunction(poly)%GaussRWeight(1:n) )
   allocate( StepFunction(poly)%SigmaCompFlag(1:jmax_step) )
   allocate( StepFunction(poly)%sigma_L(1:n,1:jmax_step) )
   allocate( StepFunction(poly)%sigma_LL(1:kmax_half,1:kmax_half,1:n) )
   StepFunction(poly)%sigma_L = CZERO
   StepFunction(poly)%sigma_LL = CZERO
!
!  ===================================================================
!  determine np_index and rp_index.
!  ===================================================================
   StepFunction(poly)%np_index(1:2)=0
   StepFunction(poly)%rp_index(1:2)=ZERO
!  -------------------------------------------------------------------
   nbnd = getNumPlanes(poly)
   ncorn = getNumCorners(poly)
   xpt => getPlane(poly)
   corn => getCorner(poly)
!  -------------------------------------------------------------------
   n=0
   do np=1,nbnd
      if(abs(xpt(1,np)) <= TEN2m10 .and. abs(xpt(2,np)) <= TEN2m10) then
         n=n+1
         StepFunction(poly)%np_index(n)=np
         if (sigma(ZERO,ZERO,xpt(3,np),xpt,nbnd,1) == 1) then
            do j = 1, ncorn
               if(abs(corn(3,j)-xpt(3,np)) < TEN2m10) then
                  StepFunction(poly)%rp_index(n)=sqrt(corn(1,j)*corn(1,j)+ &
                                                      corn(2,j)*corn(2,j)+ &
                                                      corn(3,j)*corn(3,j))
               endif
            enddo
         endif
      endif
   enddo
   nullify( xpt, corn )
!
!  ===================================================================
!  compute and store the step function for each radial Gaussian 
!  points and jl component.
!  ===================================================================
   StepFunction(poly)%SigmaCompFlag = 0
   current_r = -ONE
   n = 0
   do ic = 2,StepFunction(poly)%NumCriticalRs
      do ig = 1, ng(ic-1)
         n = n + 1
         StepFunction(poly)%GaussR(n) = rg(ig,ic-1)
         StepFunction(poly)%GaussRWeight(n) = wg(ig,ic-1)
         do jl = 1, jmax_step
            StepFunction(poly)%sigma_L(n,jl) =                        &
                getStepFunction(poly,jl,rg(ig,ic-1))
            if ( abs(StepFunction(poly)%sigma_L(n,jl)) > sigma_tol) then
               StepFunction(poly)%SigmaCompFlag(jl) = 1
            endif
         enddo
      enddo
   enddo
   if ( n /= StepFunction(poly)%NumGaussRs ) then
      call ErrorHandler('setupStepFunction','n <> StepFunction(poly)%NumGaussRs', &
                        n,StepFunction(poly)%NumGaussRs)
   endif
!
   do jl = 1, jmax_step
      if (StepFunction(poly)%SigmaCompFlag(jl) == 0) then
         StepFunction(poly)%sigma_L(:,jl) = CZERO
      endif
   enddo
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
   do n = 1, StepFunction(poly)%NumGaussRs
      do kl2 = 1, kmax_half
         do kl1 = 1, kmax_half
            ssum = CZERO
            LoopKl : do j=1,nj3(kl1,kl2)
               kl=kj3(j,kl1,kl2)
               if ( kl>StepFunction(poly)%kmax_step ) cycle LoopKl
               jl = jofk(kl)
               m = mofk(kl)
               if (m >= 0) then
                  ssum = ssum + cgnt(j,kl1,kl2)*StepFunction(poly)%sigma_L(n,jl)
               else
                  ssum = ssum + m1m(m)*cgnt(j,kl1,kl2)*conjg(StepFunction(poly)%sigma_L(n,jl))
               endif
            enddo LoopKl
            StepFunction(poly)%sigma_LL(kl1,kl2,n) = ssum
         enddo
      enddo
   enddo
!
   nullify( nj3, kj3, cgnt )
   deallocate( rg, wg, ng )
!
   end subroutine setupStepFunction
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endStepFunction()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
   implicit none
!
   integer (kind=IntKind) :: i
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endStepFunction',                            &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   endif
!
   do i = 1, NumPolyhedra
      deallocate( StepFunction(i)%SigmaCompFlag )
      deallocate( StepFunction(i)%sigma_L )
      deallocate( StepFunction(i)%sigma_LL )
      deallocate( StepFunction(i)%CriticalR )
      deallocate( StepFunction(i)%GaussR )
      deallocate( StepFunction(i)%GaussRWeight )
      deallocate( StepFunction(i)%GaussCosTheta )
      deallocate( StepFunction(i)%GaussCosThetaWeight )
!     deallocate( StepFunction(i)%sigma_L, StepFunction(i)%sigma_LL,  &
!                 StepFunction(i)%CriticalR,                          &
!                 StepFunction(i)%GaussR,                             &
!                 StepFunction(i)%GaussRWeight,                       &
!                 StepFunction(i)%GaussCosTheta,                      &
!                 StepFunction(i)%GaussCosThetaWeight )
   enddo
   deallocate( StepFunction, wylm )
   call endIntegerFactors()
   nullify( clm )
!
   if (allocated(fflag)) then
      deallocate(fflag)
   endif
   if (allocated(gflag)) then
      deallocate(gflag, ga)
   endif
   kmax_ff = 0; kmax_gf = 0
!
   if ( max_nr/=0 ) then
      deallocate( sqrt_r, fr2 )
      deallocate( res, ces, rfg, cfg, rr0 )
   endif
   max_nr = 0
!
   Initialized = .false.
!
   end subroutine endStepFunction
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genFactors(lmax)
!  ===================================================================
   implicit   none
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: kmax, jmax, l, m, kl, jl, n
!
   kmax=(lmax+1)*(lmax+1)
   jmax=(lmax+1)*(lmax+2)/2
!
   allocate( lofk(kmax), mofk(kmax), jofk(kmax), lofj(jmax), mofj(jmax), &
             m1m(-lmax:lmax) )
!
!  ===================================================================
!  calculate the factors: lofk, mofk, jofk, lofj, and mofj............
!  ===================================================================
   jl=0
   kl=0
   do l=0,lmax
      n=((l+1)*(l+2))/2-l
      do m=-l,l
         kl=kl+1
         lofk(kl)=l
         mofk(kl)=m
         jofk(kl)=n+abs(m)
      enddo
      do m=0,l
         jl=jl+1
         lofj(jl)=l
         mofj(jl)=m
      enddo
   enddo
!
!  ===================================================================
!  calculate the factor (-1)**m and store in m1m(-lmax:lmax)..........
!  ===================================================================
   m1m(-lmax)=1-2*mod(lmax,2)
   do m=-lmax+1,lmax
      m1m(m)=-m1m(m-1)
   enddo
!
   end subroutine genFactors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStepFunction1(n,jl,r) result(sf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, jl
!
   real (kind=RealKind), intent(in) :: r
!
   complex (kind=CmplxKind) :: sf
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction',                            &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (jl < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction','invalid jl index',jl)
!     ----------------------------------------------------------------
   else if (r < ZERO) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction','invalid r',r)
!     ----------------------------------------------------------------
   endif 
!
   if (n /= current_polyhedron .or. abs(r-current_r) > TEN2m12) then
!     ================================================================
!     given n, jl, and r, calculates the step function, sf.
!     ----------------------------------------------------------------
      call stepyll(n,StepFunction(n)%jmax_step,r)
!     ----------------------------------------------------------------
      current_polyhedron = n
      current_r = r
   endif
   sf = wylm(jl)
!
   end function getStepFunction1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStepFunction2(n,jl) result(sf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, jl
   integer (kind=IntKind) :: nr
!
   complex (kind=CmplxKind), pointer :: sf(:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction',                            &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (jl < 1 .or. jl > StepFunction(n)%jmax_step) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction','invalid jl index',jl)
!     ----------------------------------------------------------------
   endif 
!
   nr = StepFunction(n)%NumGaussRs
   sf => StepFunction(n)%sigma_L(1:nr,jl)
!
   end function getStepFunction2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStepFunction3(n) result(sf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: nr, jmax
!
   complex (kind=CmplxKind), pointer :: sf(:,:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction',                            &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getStepFunction','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   endif 
!
   nr = StepFunction(n)%NumGaussRs
   jmax = StepFunction(n)%jmax_step
   sf => StepFunction(n)%sigma_L(1:nr,1:jmax)
!
   end function getStepFunction3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSigmaLL(n) result(sig)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: nr, kmax
!
   complex (kind=CmplxKind), pointer :: sig(:,:,:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getSigmaLL',                                 &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getSigmaLL','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   endif 
!
   nr = StepFunction(n)%NumGaussRs
   kmax = StepFunction(n)%kmax_half
   sig => StepFunction(n)%sigma_LL(1:kmax,1:kmax,1:nr)
!
   end function getSigmaLL
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine interpolateStepFunction(n,nr,r,kmax1,kmax2,sig)
!  ===================================================================
!  this is a fast way to get the step function, that is to interpolate
!  the step function calculated on Gaussian points
!  *******************************************************************
!
   use InterpolationModule, only : getInterpolation
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, nr, kmax1, kmax2
   integer (kind=IntKind) :: kl, kl1, kl2, ir, ig, ic, mr, kr, ng, nc, kmin
   integer (kind=IntKind) :: rinc, jc, jg, jr
   integer (kind=IntKind), allocatable :: kg1(:), kg2(:), ginc(:), nginc(:)
!
   real (kind=RealKind), intent(in) :: r(nr)
   real (kind=RealKind), pointer :: rg(:), rc(:)
   real (kind=RealKind) :: rmin, rmax, err
!
   complex (kind=CmplxKind), intent(out) :: sig(nr,kmax1,kmax2)
   complex (kind=CmplxKind), allocatable :: fac(:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('interpolateStepFunction',                    &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('interpolateStepFunction','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (nr < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('interpolateStepFunction','invalid number of r mesh',nr)
!     ----------------------------------------------------------------
   else if (kmax1 < 1 .or. kmax1 > StepFunction(n)%kmax_half) then
!     ----------------------------------------------------------------
      call ErrorHandler('interpolateStepFunction','invalid kmax1',kmax1)
!     ----------------------------------------------------------------
   else if (kmax2 < 1 .or. kmax2 > StepFunction(n)%kmax_half) then
!     ----------------------------------------------------------------
      call ErrorHandler('interpolateStepFunction','invalid kmax2',kmax2)
!     ----------------------------------------------------------------
   endif
!
   sig(1:nr,1:kmax1,1:kmax2) = CZERO
!
   rmin = StepFunction(n)%rmin
   rmax = StepFunction(n)%rmax
!
   if (r(1) >= rmax) return
!
   mr = nr
   LOOP_ir1: do ir = 1, nr
      if (r(ir) > rmin+TEN2m12) then
         mr = ir - 1
         exit LOOP_ir1
      endif
   enddo LOOP_ir1
!
   kr = nr
   LOOP_ir2: do ir = nr, 1, -1
      if (r(ir) <= rmax) then
         kr = ir
         exit LOOP_ir2
      endif
   enddo LOOP_ir2
!
   if (mr > 0) then
      kmin = min(kmax1, kmax2)
      do kl = 1,kmin
         sig(1:mr,kl,kl) = CONE
      enddo
   endif
!
   if (mr == nr) then
      return
   endif
!
!  mesh points on r(mr + 1), ..., r(kr) will be interpolated
!
   ng = StepFunction(n)%NumGaussRs
   rg => StepFunction(n)%GaussR(1:ng)
!
   nc = StepFunction(n)%NumCriticalRs
   rc => StepFunction(n)%CriticalR(1:nc)
!
   allocate(fac(ng), kg1(kr-mr), kg2(kr-mr), ginc(ng), nginc(nc-1))
!
   jc = 1
   nginc(1:nc-1) = 0
   do ig = 1, ng
      LOOP_ic1: do ic = jc, nc-1
         if (rg(ig) >= rc(ic) .and. rg(ig) <= rc(ic+1)) then
            ginc(ig) = ic
            nginc(ic) = nginc(ic) + 1
            jc = ic
            exit LOOP_ic1
         endif
      enddo LOOP_ic1
   enddo
!
   jc = 1
   jg = 1
   do ir = mr+1,kr
      rinc = 0
      LOOP_ic2: do ic = jc, nc-1
         if (r(ir) >= rc(ic) .and. r(ir) <= rc(ic+1)) then
            rinc = ic
            jc = ic
            exit LOOP_ic2
         endif
      enddo LOOP_ic2
      if (rinc == 0) then
         call ErrorHandler('interpolateStepFunction','bad r-mesh',r(ir))
      endif
      LOOP_ig: do ig = jg, ng
         if (ginc(ig) == rinc) then
            jg = ig
            exit LOOP_ig
         endif
      enddo LOOP_ig
      kg1(ir-mr) = jg
      kg2(ir-mr) = nginc(rinc) + jg - 1
   enddo
!
   do kl2 =1, kmax2
      do kl1 =1, kmax1
         do ig = 1, ng
            fac(ig) = StepFunction(n)%sigma_LL(kl1,kl2,ig)
         enddo
         do ir = mr+1, kr
            jr = ir - mr
!           ----------------------------------------------------------
            sig(ir,kl1,kl2) =                                           &
                getInterpolation(kg2(jr)-kg1(jr)+1,rg(kg1(jr):kg2(jr)), &
                                 fac(kg1(jr):kg2(jr)),r(ir),err)
!           ----------------------------------------------------------
         enddo
      enddo
   enddo
!
   deallocate(kg1, kg2, fac, ginc, nginc)
   nullify(rg, rc)
!
   end subroutine interpolateStepFunction
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegration0(n,nr,r,k,f,v_mr,truncated) result(v)
!  ===================================================================
!
!  returns the volume integration of atomic cell n:
!
!     v = int_{cell n} d^3r r^k*f(r)
!
!  where f(r) is a real function
!  *******************************************************************
   use InterpolationModule, only : getInterpolation
!
   use IntegrationModule, only : calIntegration
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, nr, k
!
   real (kind=RealKind), intent(in) :: r(nr)
   real (kind=RealKind), intent(in) :: f(nr)
!
   real (kind=RealKind), optional :: v_mr
!
   logical, intent(in), optional :: truncated
!
   logical :: isTruncated
!
   integer (kind=IntKind) :: ig, ng, ir, mr, mt
!
   real (kind=RealKind) :: err, rt, rmax, rmin
   real (kind=RealKind) :: fg, v, g0, vng
   real (kind=RealKind), pointer :: rg(:), wg(:)
   real (kind=RealKind), allocatable :: g(:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration',                       &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration0','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (nr < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration0','invalid number of r mesh',nr)
!     ----------------------------------------------------------------
!   else if (k < 0) then
!     ----------------------------------------------------------------
!      call WarningHandler('getVolumeIntegration0','k < 0',k)
!     ----------------------------------------------------------------
   endif
!
   call setupIntWkSpace(nr,r)
!
   rmin = StepFunction(n)%rmin
   rmax = StepFunction(n)%rmin          ! By default, will truncate
   isTruncated =.true.
!
   if ( present(truncated) ) then
      if (.not.truncated) then
         rmax = r(nr)                      ! will not truncate
         if ( r(nr)<rmin ) then
            rmin = r(nr)
         endif
         isTruncated = .false.
      endif
   endif
!
   mr = nr
   if ( isTruncated ) then
      LOOP_ir: do ir = nr, 1, -1
         if ( r(ir) < rmax + TEN2m12) then
            mr = ir
            mt = mr
            exit LOOP_ir
         endif
      enddo LOOP_ir
      if ( mr == nr ) then
         isTruncated = .false.
      endif
   else
      LOOP_irt: do ir = nr, 1, -1
         if ( r(ir) < rmin + TEN2m12) then
            mt = ir
            exit LOOP_irt
         endif
      enddo LOOP_irt
   endif
!
   allocate( g(0:mr+1) )
!
   g(0) = ZERO
   v = ZERO
   if (mr > 0) then
      if ( abs( rmax-r(mr) ) < TEN2m12 ) then
         do ir = 1,mr
            g(ir) = f(ir)*r(ir)
         enddo
!        -------------------------------------------------------------
         call calIntegration(mr+1,sqrt_r(0:mr),g(0:mr),g0,2*k+3)
!        -------------------------------------------------------------
         v = g0*PI4*TWO
!
         if ( present(v_mr) ) then
!           ----------------------------------------------------------
            call calIntegration(mt+1,sqrt_r(0:mt),g(0:mt),g0,2*k+3)
!           ----------------------------------------------------------
            v_mr = g0*PI4*TWO
        endif
      else
!        -------------------------------------------------------------
         fg = getInterpolation(nr,r(1:nr),f(1:nr),rmax,err)
!        -------------------------------------------------------------
         rt = sqrt_r(mr+1)
         sqrt_r(mr+1) = sqrt(rmax)
         do ir = 1,mr
            g(ir) = f(ir)*r(ir)
         enddo
         g(mr+1) = fg
!        -------------------------------------------------------------
         call calIntegration(mr+2,sqrt_r(0:mr+1),g(0:mr+1),g0,2*k+3)
!        -------------------------------------------------------------
         v = g0*PI4*TWO
         sqrt_r(mr+1) = rt
         if ( present(v_mr) ) then
            if ( abs(rmin-r(mt))<TEN2m12 ) then
               v_mr = v
            else
!              -------------------------------------------------------
               fg = getInterpolation(nr,r(1:nr),f(1:nr),rmin,err)
!              -------------------------------------------------------
               rt = sqrt_r(mt+1)
               sqrt_r(mt+1) = sqrt(rmin)
               do ir = 1,mt
                  g(ir) = f(ir)*r(ir)
               enddo
               g(mt+1) = fg
!              -------------------------------------------------------
               call calIntegration(mt+2,sqrt_r(0:mt+1),g(0:mt+1),g0,2*k+3)
!              -------------------------------------------------------
               v_mr = g0*PI4*TWO
               sqrt_r(mt+1) = rt
            endif
         endif
      endif
   endif
!
   deallocate( g )
!
   if ( .not.isTruncated ) then
      return
   endif
!
   ng = StepFunction(n)%NumGaussRs
   rg => StepFunction(n)%GaussR(1:ng)
   wg => StepFunction(n)%GaussRWeight(1:ng)
!
   if ( r(nr)<StepFunction(n)%rmax-Ten2m12 ) then
      call hunt(ng,rg,r(nr),ir)
   else
      ir = ng
   endif
!
   vng = ZERO
   if ( ir<ng ) then
      vng = vng + TWO*sqrt(PI)*wg(ir+1)*f(nr)*(r(nr)**(2+k))*             &
                  real(StepFunction(n)%sigma_L(ir+1,1),kind=RealKind)
   endif
   do ig = ir,1,-1
!     ----------------------------------------------------------------
      fg = getInterpolation(nr-mr,r(mr+1:nr),f(mr+1:nr),rg(ig),err)
!     ----------------------------------------------------------------
      vng = vng + TWO*sqrt(PI)*wg(ig)*fg*(rg(ig)**(2+k))*             &
                  real(StepFunction(n)%sigma_L(ig,1),kind=RealKind)
   enddo
   v = v +vng
   nullify(rg, wg)
!
   end function getVolumeIntegration0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegration1(n,nr,r,kmax,k,frk,v_mr,truncated) result(v)
!  ===================================================================
!
!  returns the volume integration of atomic cell n:
!
!     v = int_{cell n} d^3r f(vec{r})
!
!  where f(vec{r}) = sum_{L} f(r,L) * Y_{L}
!
!  the input frk = f(r,L)*r^{k}
!  *******************************************************************
   use InterpolationModule, only : getInterpolation, FitInterp
!
   use IntegrationModule, only : calIntegration
   implicit none
!
   logical, intent(in), optional :: truncated
!
   integer (kind=IntKind), intent(in) :: n, nr, kmax, k
!
   real (kind=RealKind), intent(in) :: r(1:nr)
!
   complex (kind=CmplxKind), optional :: v_mr
   complex (kind=CmplxKind), dimension(:,:) :: frk
!
   logical :: isTruncated
!
   integer (kind=IntKind) :: jl, ml, kl, ig, ng, ir, mr, mt, lmax, l
   integer (kind=IntKind) :: izamax, i
!
   real (kind=RealKind) :: fac, err, rmax, rmin, fr, fi
   real (kind=RealKind), pointer :: rg(:), wg(:)
!
   complex (kind=CmplxKind) :: fg, v, vig, vigm, vng, vtmp, fgmm, v_m0
   complex (kind=CmplxKind), pointer :: sigma_L(:,:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration',                       &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (nr < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid number of r mesh',nr)
!     ----------------------------------------------------------------
   else if (kmax < 1 ) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid kmax',kmax)
!     ----------------------------------------------------------------
   endif
!
   if (.not.allocated(fflag)) then
      allocate(fflag(kmax))
      kmax_ff = kmax
   else if(kmax_ff < kmax) then
      deallocate(fflag)
      allocate(fflag(kmax))
      kmax_ff = kmax
   endif
!
   do kl = 1, kmax
      fflag(kl) = 0
      ir = izamax(nr,frk(1:nr,kl),1)
      if (abs(frk(ir,kl)) > TEN2m8*(TEN**(-2*lofk(kl)))) then
         fflag(kl) = 1
      endif
   enddo
!
   if (maxval(fflag(1:kmax)) == 0) then
      v = CZERO
      if (present(v_mr)) then
         v_mr = CZERO
      endif
      return
   endif
!
   rmin = StepFunction(n)%rmin
   rmax = StepFunction(n)%rmin          ! By default, will truncate
   isTruncated =.true.
!
   if ( present(truncated) ) then
      if (.not.truncated) then
         rmax = r(nr)                      ! will not truncate
         if ( r(nr)<rmin ) then
            rmin = r(nr)
         endif
         isTruncated = .false.
      endif
   endif
!
   mr = nr
   if ( isTruncated ) then
      LOOP_ir: do ir = nr, 1, -1
         if ( r(ir) < rmax + TEN2m12) then
            mr = ir
            mt = mr
            exit LOOP_ir
         endif
      enddo LOOP_ir
      if ( mr == nr ) then
         isTruncated = .false.
      endif
   else if ( present(v_mr) ) then
      LOOP_irt: do ir = nr, 1, -1
         if ( r(ir) < rmin + TEN2m12) then
            mt = ir
            exit LOOP_irt
         endif
      enddo LOOP_irt
   endif
!
   v = CZERO
!
   if (mr > 0) then
!     ----------------------------------------------------------------
      call setupIntWkSpace(nr,r)
!     ----------------------------------------------------------------
      cfg(1:nr) = frk(1:nr,1)
!     ----------------------------------------------------------------
      call calIntegrationToR( mode, nr, r(1:nr), k, cfg(0:nr),        &
                              rmax, mr, ces(0:nr) )
!     ----------------------------------------------------------------
      if ( abs(rmax - r(mr)) < TEN2m12 ) then
         v = ces(mr)*TWO*sqrt(PI)
      else
         v = ces(mr+1)*TWO*sqrt(PI)
      endif
      if ( present(v_mr) ) then
         if ( abs( r(mt) - rmin ) < TEN2m12 .and. rmax >= rmin ) then
            v_mr = ces(mt)*TWO*sqrt(PI)
         else
!           ----------------------------------------------------------
            call calIntegrationToR( mode, nr, r(1:nr), k, cfg(0:nr),  &
                                    rmin, mt, ces(0:nr) )
!           ----------------------------------------------------------
            v_mr = ces(mt+1)*TWO*sqrt(PI)
         endif
      endif
   endif
!
   if ( .not.isTruncated ) then
      return
   endif
!
!   write(6,*) "getVolumeIntegration1:: Step"
   ng = StepFunction(n)%NumGaussRs
   rg => StepFunction(n)%GaussR(1:ng)
   wg => StepFunction(n)%GaussRWeight(1:ng)
   sigma_L => StepFunction(n)%sigma_L
!
   if ( r(nr)<StepFunction(n)%rmax-Ten2m10 ) then
      call hunt(ng,rg,r(nr),ir)
   else
      ir = ng
   endif
!
   cfg = CZERO
   kl = min(kmax,StepFunction(n)%kmax_step)
   lmax = lofk(kl)
   vng  = CZERO
   v_m0 = CZERO
   if (ir<ng) then
      vig  = CZERO
      do l = lmax,0,-1
         vigm = CZERO
         do ml = l, 0, -1
            kl = (l+1)*(l+1) - l + ml
            jl = ((l+1)*(l+2))/2 -l + ml
            fg = CZERO
            if ( fflag(kl)==1 ) then 
               fg = frk(nr,kl)
            endif
            if (ml == 0) then
               vig  = vig  + conjg(sigma_L(ir+1,jl))*fg
            else
               fac = m1m(ml)
               fgmm = CZERO
               if (fflag(kl-2*ml)==1) then
                  fgmm = frk(nr,kl-2*ml)
               endif
               vtmp = conjg(sigma_L(ir+1,jl))*fg+     &
                      fac*sigma_L(ir+1,jl)*fgmm
               vigm = vigm + vtmp
            endif
         enddo
         vig = vig + vigm
      enddo
      if ( k==2 ) then
         fac = wg(ir+1)
      else
         fac = (r(nr)**(2-k))*wg(ir+1)
      endif
      vng = vng + fac*vig
   endif
   do ig = ir,1,-1
      vig  = CZERO
      LOOP_l: do l = lmax,0,-1
         vigm = CZERO
         do ml = l, 0, -1
            kl = (l+1)*(l+1) - l + ml
            jl = ((l+1)*(l+2))/2 -l + ml
            fg = CZERO
            if ( fflag(kl)==1 ) then 
!              -------------------------------------------------------
               fg = getInterpolation( nr-mr,r(mr+1:nr),               &
                                      frk(mr+1:nr,kl),rg(ig),err)
!              -------------------------------------------------------
            endif
            if (ml == 0) then
               vigm = vigm + conjg(sigma_L(ig,jl))*fg
            else
               fac = real(m1m(ml),kind=RealKind)
               fgmm = CZERO
               if ( fflag(kl-2*ml)==1 ) then 
!                 ----------------------------------------------------
                  fgmm = getInterpolation( nr-mr,r(mr+1:nr),          &
                                      frk(mr+1:nr,kl-2*ml),rg(ig),err)
!                 ----------------------------------------------------
               endif
               vtmp = conjg(sigma_L(ig,jl))*fg +                    &
                      fac*sigma_L(ig,jl)*fgmm
               vigm = vigm + vtmp
            endif
         enddo
         vig = vig + vigm
      enddo LOOP_l
      if ( k==2 ) then
         fac = wg(ig)
      else
         fac = (rg(ig)**(2-k))*wg(ig)
      endif
      vng = vng + fac*vig
   enddo
!
   v = v + vng
!
   end function getVolumeIntegration1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegration1_jl(n,nr,r,kmax,jmax,k,frk,v_mr,truncated,tol_in) result(v)
!  ===================================================================
!
!  returns the volume integration of atomic cell n:
!
!     v = int_{cell n} d^3r f(vec{r})
!
!  where f(vec{r}) = sum_{L} f(r,L) * Y_{L}
!
!  the input frk = f(r,L)*r^{k}
!  *******************************************************************
   use InterpolationModule, only : getInterpolation, FitInterp
!
   use IntegrationModule, only : calIntegration
   implicit none
!
   logical, intent(in), optional :: truncated
!
   integer (kind=IntKind), intent(in) :: n, nr, kmax, k, jmax
!
   real (kind=RealKind), intent(in) :: r(1:nr)
   real (kind=RealKind), intent(in), optional :: tol_in
!
   real (kind=RealKind), intent(out) :: v_mr
!
   complex (kind=CmplxKind), dimension(:,:) :: frk
!
   integer (kind=IntKind) :: jl, ml, ig, ng, ir, mr, mt, l
   integer (kind=IntKind) :: izamax, jmin, lmax
!
   logical :: isTruncated
!
   real (kind=RealKind) :: fac, err, rmax, rmin, tol, v
   real (kind=RealKind), pointer :: rg(:), wg(:)
!
   complex (kind=CmplxKind) :: fg, vig, vigm, vtmp, vng
   complex (kind=CmplxKind), pointer :: sigma_L(:,:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration',                       &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (nr < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid number of r mesh',nr)
!     ----------------------------------------------------------------
   else if (kmax < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid kmax',kmax)
!     ----------------------------------------------------------------
   endif
!
   if (.not.allocated(fflag)) then
      allocate(fflag(kmax))
      kmax_ff = kmax
   else if(kmax_ff < kmax) then
      deallocate(fflag)
      allocate(fflag(kmax))
      kmax_ff = kmax
   endif
!
   if ( present(tol_in) ) then
      tol = tol_in
   else
      tol = sigma_tol
   endif
!
   do jl = 1, jmax
      fflag(jl) = 0
      ir = izamax(nr,frk(1:nr,jl),1)
      if (abs(frk(ir,jl)) > tol*(TEN**(-2*lofj(jl)))) then
         fflag(jl) = 1
      endif
   enddo
!
   if (maxval(fflag(1:jmax)) == 0) then
      v = CZERO
      v_mr = CZERO
      return
   endif
!
   rmin = StepFunction(n)%rmin
   rmax = StepFunction(n)%rmin          ! By default, will truncate
   isTruncated =.true.
!
   if ( present(truncated) ) then
      if (.not.truncated) then
         rmax = r(nr)                      ! will not truncate
         if ( r(nr)<rmin ) then
            rmin = r(nr)
         endif
         isTruncated = .false.
      endif
   endif
!
   mr = nr
   if ( isTruncated ) then
      LOOP_ir: do ir = nr, 1, -1
         if ( r(ir) < rmax + TEN2m12) then
            mr = ir
            mt = mr
            exit LOOP_ir
         endif
      enddo LOOP_ir
      if ( mr == nr ) then
         isTruncated = .false.
      endif
   else
      LOOP_irt: do ir = nr, 1, -1
         if ( r(ir) < rmin + TEN2m12) then
            mt = ir
            exit LOOP_irt
         endif
      enddo LOOP_irt
   endif
!
   v = CZERO
!
   if (mr > 0) then
!
!     ----------------------------------------------------------------
      call setupIntWkSpace(nr,r)
!     ----------------------------------------------------------------
      cfg(1:nr) = frk(1:nr,1)
!     ----------------------------------------------------------------
      call calIntegrationToR( mode, nr, r(1:nr), k, cfg(0:nr),      &
                              rmax, mr, ces(0:nr) )
!     ----------------------------------------------------------------
      if ( abs(rmax - r(mr)) < TEN2m12 ) then
         v = ces(mr)*TWO*sqrt(PI)
      else
         v = ces(mr+1)*TWO*sqrt(PI)
      endif
      if ( abs( r(mt) - rmin ) < TEN2m12 .and. rmax >= rmin ) then
         v_mr = ces(mt)*TWO*sqrt(PI)
      else
!        -------------------------------------------------------------
         call calIntegrationToR( mode, nr, r(1:nr), k, cfg(0:nr),&
                                 rmin, mt, ces(0:nr) )
!        -------------------------------------------------------------
         v_mr = ces(mt+1)*TWO*sqrt(PI)
      endif
   endif
!
   if ( .not.isTruncated ) then
      return
   endif
!
   ng = StepFunction(n)%NumGaussRs
   rg => StepFunction(n)%GaussR(1:ng)
   wg => StepFunction(n)%GaussRWeight(1:ng)
   sigma_L => StepFunction(n)%sigma_L
!
   if ( r(nr)<StepFunction(n)%rmax-Ten2m10 ) then
      call hunt(ng,rg,r(nr),ir)
   else
      ir = ng
   endif
!
   jmin = (StepFunction(n)%lmax_step+1)*(StepFunction(n)%lmax_step+2)/2
   jmin = min(jmax,jmin)
   lmax = lofj(jmin)
   vng = CZERO
   if (ir<ng) then
      vig  = CZERO
      jl = jmin
      do l = lmax,0,-1
         vigm = CZERO
         do ml = l,0,-1
            fg = frk(nr,jl)
            if ( ml==0 ) then
               vigm = vigm + conjg(sigma_L(ir+1,jl))*fg
            else
               vtmp = conjg(sigma_L(ir+1,jl))*fg +       &
                      sigma_L(ir+1,jl)*conjg(fg)
               vigm = vigm + vtmp
            endif
            jl = jl - 1
         enddo
         vig = vig + vigm
      enddo
!
      if ( k==2 ) then
         vng = vng + wg(ir+1)*vig
      else
         vng = vng + ((r(nr)**(2-k))*wg(ir+1))*vig
      endif
   endif
   do ig = ir,1,-1
      vig  = CZERO
      jl = jmin
      do l = lmax,0,-1
         vigm = CZERO
         do ml = l,0,-1
            if ( fflag(jl)==1 ) then
!              -------------------------------------------------------
               fg = getInterpolation(nr-mr,r(mr+1:nr),frk(mr+1:nr,jl),&
                                     rg(ig),err)
!              -------------------------------------------------------
               if ( ml==0 ) then
                  vigm = vigm + conjg(sigma_L(ig,jl))*fg
               else
                  vtmp = conjg(sigma_L(ig,jl))*fg +         &
                         sigma_L(ig,jl)*conjg(fg)
                  vigm = vigm + vtmp
              endif
            endif
            jl = jl - 1
         enddo
         vig = vig + vigm
      enddo
!
      if ( k==2 ) then
         fac = wg(ig)
      else
         fac = (rg(ig)**(2-k))*wg(ig)
      endif
      vng = vng + fac*vig
!
   enddo
!
   v = v + vng
!
   end function getVolumeIntegration1_jl
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegration1r_jl(n,nr,r,jmax,k,frk,v_mr,v_jl,    &
                                      truncated,tol_in) result(v)
!  ===================================================================
!
!  returns the volume integration of atomic cell n:
!
!     v = int_{cell n} d^3r f(vec{r})
!
!  where f(vec{r}) = sum_{L} f(r,L) * Y_{L}
!
!  the input frk = f(r,L)*r^{k}
!
!  returns v_mr, the volume integration over the spherical volume whose
!          radius is at the mesh point closest to the Muffin-tin radius
!  Optionally, it returns v_jl, the volume integration of each jl component
!              where,     v_jl = int_{cell n} d^3r f(r,l,m) * Y_{l,m}
!                               + int_{cell n} d^3r f(r,l,-m) * Y_{l,-m}
!                              = 2*Re[int_{cell n} d^3r f(r,l,m) * Y_{l,m}], if m <> 0
!                          or  = Re[int_{cell n} d^3r f(r,l,0) * Y_{l,0}], if m = 0
!  *******************************************************************
   use InterpolationModule, only : getInterpolation, FitInterp
!
   use IntegrationModule, only : calIntegration
   implicit none
!
   logical, intent(in), optional :: truncated
!
   integer (kind=IntKind), intent(in) :: n, nr, k, jmax
!
   real (kind=RealKind), intent(in) :: r(1:nr)
   real (kind=RealKind), intent(in), optional :: tol_in
   real (kind=RealKind), intent(out) :: v_mr
   real (kind=RealKind), intent(out), optional :: v_jl(jmax)
   real (kind=RealKind) :: v
!
   complex (kind=CmplxKind), intent(in), dimension(:,:) :: frk
!
   integer (kind=IntKind) :: jl, ml, ig, ng, ir, mr, mt, l
   integer (kind=IntKind) :: izamax, jmin, lmax
!
   logical :: isTruncated
!
   real (kind=RealKind) :: fac, err, rmax, rmin, tol
   real (kind=RealKind), pointer :: rg(:), wg(:)
   real (kind=RealKind) :: vig, vigm, vtmp, vng
!
   complex (kind=CmplxKind) :: fg
   complex (kind=CmplxKind), pointer :: sigma_L(:,:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration',                       &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (nr < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid number of r mesh',nr)
!     ----------------------------------------------------------------
   else if ( r(nr)<StepFunction(n)%rmax-Ten2m10 ) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','r(nr) <  bounding sphere radius', &
                        r(nr),StepFunction(n)%rmax)
!     ----------------------------------------------------------------
   else if (jmax < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid jmax',jmax)
!     ----------------------------------------------------------------
   endif
!
   if (.not.allocated(fflag)) then
      allocate(fflag(jmax))
      kmax_ff = jmax
   else if(kmax_ff < jmax) then
      deallocate(fflag)
      allocate(fflag(jmax))
      kmax_ff = jmax
   endif
!
   if ( present(tol_in) ) then
      tol = tol_in
   else
      tol = sigma_tol
   endif
!
   do jl = 1, jmax
      fflag(jl) = 0
      ir = izamax(nr,frk(1:nr,jl),1)
      if (abs(frk(ir,jl)) > tol*(TEN**(-2*lofj(jl)))) then
         fflag(jl) = 1
      endif
   enddo
!
   if (present(v_jl)) then
      v_jl = ZERO
   endif
!
   if (maxval(fflag(1:jmax)) == 0) then
      v = ZERO
      v_mr = ZERO
      return
   endif
!
   rmin = StepFunction(n)%rmin
   rmax = StepFunction(n)%rmin          ! By default, will truncate
   isTruncated =.true.
!
   if ( present(truncated) ) then
      if (.not.truncated) then
         rmax = r(nr)                      ! will not truncate
         if ( r(nr)<rmin ) then
            rmin = r(nr)
         endif
         isTruncated = .false.
      endif
   endif
!
   mr = nr
   if ( isTruncated ) then
      LOOP_ir: do ir = nr, 1, -1
         if ( r(ir) < rmax + TEN2m12) then
            mr = ir
            mt = mr
            exit LOOP_ir
         endif
      enddo LOOP_ir
      if ( mr == nr ) then
         isTruncated = .false.
      endif
   else
      LOOP_irt: do ir = nr, 1, -1
         if ( r(ir) < rmin + TEN2m12) then
            mt = ir
            exit LOOP_irt
         endif
      enddo LOOP_irt
   endif
!
   v = CZERO
!
!  ===================================================================
!  Volume integration up to mesh point mr. No truncation is needed up 
!  to this point.
!  ===================================================================
   if (mr > 0) then
!     ----------------------------------------------------------------
      call setupIntWkSpace(nr,r)
!     ----------------------------------------------------------------
      cfg(1:nr) = frk(1:nr,1)
!     ----------------------------------------------------------------
      call calIntegrationToR( mode, nr, r(1:nr), k, cfg(0:nr),      &
                              rmax, mr, ces(0:nr) )
!     ----------------------------------------------------------------
      if ( abs(rmax - r(mr)) < TEN2m12 ) then
         v = real(ces(mr),kind=RealKind)*TWO*sqrt(PI)
      else
         v = real(ces(mr+1),kind=RealKind)*TWO*sqrt(PI)
      endif
      if ( abs( r(mt) - rmin ) < TEN2m12 .and. rmax >= rmin ) then
         v_mr = real(ces(mt),kind=RealKind)*TWO*sqrt(PI)
      else
!        -------------------------------------------------------------
         call calIntegrationToR( mode, nr, r(1:nr), k, cfg(0:nr),     &
                                 rmin, mt, ces(0:nr) )
!        -------------------------------------------------------------
         v_mr = real(ces(mt+1),kind=RealKind)*TWO*sqrt(PI)
      endif
      if (present(v_jl)) then
         v_jl(1) = v
      endif
   endif
!
   if ( .not.isTruncated ) then
      return
   endif
!
   ng = StepFunction(n)%NumGaussRs
   rg => StepFunction(n)%GaussR(1:ng)
   wg => StepFunction(n)%GaussRWeight(1:ng)
   sigma_L => StepFunction(n)%sigma_L
!
   jmin = (StepFunction(n)%lmax_step+1)*(StepFunction(n)%lmax_step+2)/2
   jmin = min(jmax,jmin)
   lmax = lofj(jmin)
   vng = ZERO
!
!  ===================================================================
!  Volume integration from mesh point mr to the bounding sphere. 
!  ===================================================================
   do ig = ng,1,-1
      if ( k==2 ) then
         fac = wg(ig)
      else
         fac = (rg(ig)**(2-k))*wg(ig)
      endif
!
      vig  = ZERO
      jl = jmin
      do l = lmax,0,-1
         vigm = ZERO
         do ml = l,0,-1
            if ( fflag(jl)==1 ) then
!              -------------------------------------------------------
               fg = getInterpolation(nr-mr,r(mr+1:nr),frk(mr+1:nr,jl),&
                                     rg(ig),err)
!              -------------------------------------------------------
               if ( ml==0 ) then
                  vtmp = conjg(sigma_L(ig,jl))*fg
               else
                  vtmp = TWO*real(conjg(sigma_L(ig,jl))*fg,kind=RealKind)
               endif
               vigm = vigm + vtmp
               if (present(v_jl)) then
                  v_jl(jl) = v_jl(jl) + fac*vtmp
               endif
            endif
            jl = jl - 1
         enddo
         vig = vig + vigm
      enddo
      vng = vng + fac*vig
   enddo
!
   v = v + vng
!
   end function getVolumeIntegration1r_jl
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegration2(n,nr,r,kmaxf,kf,frk,kmaxg,kg,grk,   &
                                  v_mr,truncated,tol_in) result(v)
!  ===================================================================
!
!  returns the volume integration of atomic cell n:
!
!     v = int_{cell n} d^3r f^{*}(vec{r}) * g(vec{r})
!
!  where f(vec{r}) = sum_{L} f(r,L) * Y_{L},  L <= kmaxf
!        g(vec{r}) = sum_{L} g(r,L) * Y_{L},  L <= kmaxg
!
!  the input frk = f(r,L)*r^{kf}, grk = g(r,L)*r^{kg}
!  *******************************************************************
   use InterpolationModule, only : getInterpolation
!
   use IntegrationModule, only : calIntegration
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, nr
   integer (kind=IntKind), intent(in) :: kmaxf, kf
   integer (kind=IntKind), intent(in) :: kmaxg, kg
!
   complex (kind=CmplxKind), intent(out), optional :: v_mr
!
   logical, intent(in), optional :: truncated
!
   real (kind=RealKind), intent(in), optional :: tol_in
!
   logical :: isTruncated
!
   integer (kind=IntKind) :: kmax1, kmax0, k, kls, klfs, m, mf
   integer (kind=IntKind) :: kl, ig, ng, ir, mr, mt, klf, klg
   integer (kind=IntKind) :: min_kf1, min_kg1, kminf, kming
   integer (kind=IntKind) :: izamax
!
   real (kind=RealKind), intent(in) :: r(1:nr)
   real (kind=RealKind) :: fac, err, rmax, rmin, tol
   real (kind=RealKind), pointer :: rg(:), wg(:)
!
   complex (kind=CmplxKind), intent(in) :: frk(:,:)
   complex (kind=CmplxKind), intent(in) :: grk(:,:)
   complex (kind=CmplxKind) :: fi, v, vig, v0, vng, dps
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   kmax1 = max(kmaxf, kmaxg)
   kmax0 = min(kmaxf, kmaxg)
   k = kf + kg
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration',                       &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   else if (nr < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid number of r mesh',nr)
!     ----------------------------------------------------------------
   else if (kmax0 < 1 .or. kmax1 <1 ) then
!     ----------------------------------------------------------------
      call ErrorHandler('getVolumeIntegration','invalid kmax',kmaxf,kmaxg)
!     ----------------------------------------------------------------
   endif
!
   rmin = StepFunction(n)%rmin
   rmax = StepFunction(n)%rmin          ! By default, will truncate
   isTruncated =.true.
!
   if ( present(truncated) ) then
      if (.not.truncated) then
         rmax = r(nr)                      ! will not truncate
         if ( r(nr)<rmin ) then
            rmin = r(nr)
         endif
         isTruncated = .false.
      endif
   endif
!
   mr = nr
   if ( isTruncated ) then
      LOOP_ir: do ir = nr, 1, -1
         if ( r(ir) < rmax + TEN2m12) then
            mr = ir
            mt = mr
            exit LOOP_ir
         endif
      enddo LOOP_ir
      if ( mr == nr ) then
         isTruncated = .false.
      endif
   else
      LOOP_irt: do ir = nr, 1, -1
         if ( r(ir) < rmin + TEN2m12) then
            mt = ir
            exit LOOP_irt
         endif
      enddo LOOP_irt
   endif
!
   if ( present(tol_in) ) then
      tol = tol_in
   else
      tol = sigma_tol
   endif
!
   if (.not.allocated(fflag)) then
      allocate( fflag(kmaxf) )
      kmax_ff = kmaxf
   else if (kmax_ff < kmaxf) then
      deallocate( fflag )
      allocate( fflag(kmaxf) )
      kmax_ff = kmaxf
   endif
   if (.not.allocated(gflag)) then
      allocate( gflag(kmaxg), ga(kmaxg) )
      kmax_gf = kmaxg
   else if (kmax_gf < kmaxg) then
      deallocate( gflag, ga )
      allocate( gflag(kmaxg), ga(kmaxg) )
      kmax_gf = kmaxg
   endif
!
   min_kf1 = min(kmaxf,kmax1)
   min_kg1 = min(kmaxg,kmax1)
!
   v = CZERO
   v_mr = CZERO
!
   if (mr > 0) then
      fflag(1:min_kf1) = 0
      do kl = 1, min_kf1
         ir = izamax(nr,frk(1:nr,kl),1)
         if (abs(frk(ir,kl)) > TEN2m8*(TEN**(-2*lofk(kl)))) then
            fflag(kl) = 1
         endif
      enddo
      gflag(1:min_kg1) = 0
      do kl = 1,min_kg1
         ir = izamax(nr,grk(1:nr,kl),1)
         if (abs(grk(ir,kl)) > TEN2m8*(TEN**(-2*lofk(kl)))) then
            gflag(kl) = 1
         endif
      enddo
!
!     ----------------------------------------------------------------
      call setupIntWkSpace(nr,r)
!     ----------------------------------------------------------------
!
      do kl = kmax0, 1, -1
         m = mofk(kl)
         kls = kl-2*m
         if (fflag(kls) == 1 .and. gflag(kl) == 1) then
            do ir = 1, mr
               cfg(ir) = frk(ir,kls)*grk(ir,kl)*m1m(m)
            enddo
!           ----------------------------------------------------------
            call calIntegrationToR(mode,nr,r(1:nr),k,cfg(0:nr),rmax,mr,ces(0:nr))
!           ----------------------------------------------------------
            if ( abs(rmax - r(mr)) < TEN2m12 ) then
               v = v + ces(mr)
            else
               v = v + ces(mr+1)
            endif
            if ( abs( r(mt) - rmin ) < TEN2m12 .and. rmax >= rmin ) then
               v_mr = v_mr + ces(mt)
            else
!              -------------------------------------------------------
               call calIntegrationToR(mode,nr,r(1:nr),k,cfg(0:nr),rmin,mt,ces(0:nr))
!              -------------------------------------------------------
               v_mr = v_mr + ces(mt+1)
            endif
         endif
      enddo
   endif
!
   if ( .not.isTruncated ) then
!      write(6,*)"Truncation"
      return
   endif
!
   kminf = min(kmaxf,StepFunction(n)%kmax_half)
   kming = min(kmaxg,StepFunction(n)%kmax_half)
!
   ng = StepFunction(n)%NumGaussRs
   rg => StepFunction(n)%GaussR(1:ng)
   wg => StepFunction(n)%GaussRWeight(1:ng)
!
   if ( r(nr)<StepFunction(n)%rmax-Ten2m12 ) then
      call hunt(ng,rg,r(nr),ir)
   else
      ir = ng
   endif
!
   vng = CZERO
   if (k == 2) then
      fac = ONE
      if (ir<ng) then
         vig = CZERO
         ga(1:kmaxg) = CZERO
         do klg = kming,1,-1
            if ( gflag(klg) == 1 ) then
!              -------------------------------------------------------
               ga(klg)=grk(nr,klg)
!              -------------------------------------------------------
            endif
         enddo
         do klf = kminf, 1, -1
            mf = mofk(klf)
            klfs = klf - 2*mf
            if (fflag(klfs) == 1) then
!              -------------------------------------------------------
               fi=frk(nr,klfs)*m1m(mf)
!              -------------------------------------------------------
               v0 = CZERO
               do klg = kming,1,-1
                  if (gflag(klg) == 1) then
                     v0 = v0 + StepFunction(n)%sigma_LL(klg,klf,ir+1) &
                              *ga(klg)
                  endif
               enddo
               vig = vig + v0*fi
            endif
         enddo
         vng = vng + wg(ir+1)*vig
      endif
      do ig = ir,1,-1
         vig = CZERO
         ga(1:kmaxg) = CZERO
         do klg = kming,1,-1
            if ( gflag(klg) == 1 ) then
!              -------------------------------------------------------
               ga(klg)=getInterpolation(nr-mr,r(mr+1:nr),grk(mr+1:nr,klg),rg(ig),err)
!              -------------------------------------------------------
            endif
         enddo
         do klf = kminf, 1, -1
            mf = mofk(klf)
            klfs = klf - 2*mf
            if (fflag(klfs) == 1) then
!              -------------------------------------------------------
               fi=getInterpolation(nr-mr,r(mr+1:nr),frk(mr+1:nr,klfs),rg(ig),err)*m1m(mf)
!              -------------------------------------------------------
               v0 = CZERO
               do klg = kming,1,-1
                  if (gflag(klg) == 1) then
                     v0 = v0 + StepFunction(n)%sigma_LL(klg,klf,ig)   &
                              *ga(klg)
                  endif
               enddo
               vig = vig + v0*fi
            endif
         enddo
         vng = vng + wg(ig)*vig
      enddo
   else
      if (ir<ng) then
         vig = CZERO
         fac = wg(ir+1)*(r(nr)**(2-k))
         ga(1:kmaxg) = CZERO
         do klg = kming,1,-1
            if ( gflag(klg) == 1 ) then
!              -------------------------------------------------------
               ga(klg)=grk(nr,klg)
!              -------------------------------------------------------
            endif
         enddo
         do klf = kminf, 1, -1
            mf = mofk(klf)
            klfs = klf - 2*mf
            if (fflag(klfs) == 1) then
!              -------------------------------------------------------
               fi=frk(nr,klfs)*m1m(mf)
!              -------------------------------------------------------
               v0 = CZERO
               do klg = kming,1,-1
                  if (gflag(klg) == 1) then
                     v0 = v0 + StepFunction(n)%sigma_LL(klg,klf,ir+1) &
                              *ga(klg)
                  endif
               enddo
               vig = vig + v0*fi
            endif
         enddo
         vng = vng + fac*vig
      endif
      do ig = ir,1,-1
         vig = CZERO
         fac = wg(ig)*(rg(ig)**(2-k))
         ga(1:kmaxg) = CZERO
         do klg = kming,1,-1
            if ( gflag(klg) == 1 ) then
!              -------------------------------------------------------
               ga(klg)=getInterpolation(nr-mr,r(mr+1:nr),grk(mr+1:nr,klg),rg(ig),err)
!              -------------------------------------------------------
            endif
         enddo
         do klf = kminf, 1, -1
            mf = mofk(klf)
            klfs = klf - 2*mf
            if (fflag(klfs) == 1) then
!              -------------------------------------------------------
               fi=getInterpolation(nr-mr,r(mr+1:nr),frk(mr+1:nr,klfs),rg(ig),err)*m1m(mf)
!              -------------------------------------------------------
               v0 = CZERO
               do klg = kming,1,-1
                  if (gflag(klg) == 1) then
                     v0 = v0 + StepFunction(n)%sigma_LL(klg,klf,ig)   &
                              *ga(klg)
                  endif
               enddo
               vig = vig + v0*fi
            endif
         enddo
         vng = vng + fac*vig
      enddo
   endif
   nullify(rg, wg)
!
   v = v + vng
!
   end function getVolumeIntegration2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCriticalRs(n) result(nc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: nc
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumCriticalRs',                           &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumCriticalRs','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   endif
!
   nc = StepFunction(n)%NumCriticalRs
!
   end function getNumCriticalRs
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCriticalR(n) result(rc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: nc
!
   real (kind=RealKind), pointer :: rc(:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCriticalR',                               &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCriticalR','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   endif
!
   nc = StepFunction(n)%NumCriticalRs
   rc => StepFunction(n)%CriticalR(1:nc)
!
   end function getCriticalR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGaussRs(n) result(ng)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: ng
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumGaussRs',                              &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumGaussRs','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   endif
!
   ng = StepFunction(n)%NumGaussRs
!
   end function getNumGaussRs
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussR(n) result(rg)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: ng
!
   real (kind=RealKind), pointer :: rg(:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGaussR',                               &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGaussR','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   endif
!
   ng = StepFunction(n)%NumGaussRs
   rg => StepFunction(n)%GaussR(1:ng)
!
   end function getGaussR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussRWeight(n) result(wg)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: ng
!
   real (kind=RealKind), pointer :: wg(:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGaussRWeight',                         &
                        'StepFunctionModule is not initialized')
!     ----------------------------------------------------------------
   else if (n < 1 .or. n > NumPolyhedra) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGaussRWeight','invalid polyhedron index',n)
!     ----------------------------------------------------------------
   endif
!
   ng = StepFunction(n)%NumGaussRs
   wg => StepFunction(n)%GaussRWeight(1:ng)
!
   end function getGaussRWeight
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeCriticalPoints(n)
!  ===================================================================
   use SortModule, only : QuickSort
!
   use PolyhedraModule, only : getNumPlanes, getNumCorners, getNumEdges
   use PolyhedraModule, only : getPlane, getCorner, getEdge
   use PolyhedraModule, only : isExternalPoint
!
   use MathParamModule, only : TEN2m9, TEN2m8
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: nump, numc, nume, i, j, nc, ncrit
!
   real (kind=RealKind), allocatable :: rc2(:), rcrit(:)
   real (kind=RealKind), pointer :: pp(:,:)
   real (kind=RealKind), pointer :: pc(:,:)
   real (kind=RealKind), pointer :: pe(:,:,:)
!
   nump = getNumPlanes(n)
   numc = getNumCorners(n)
   nume = getNumEdges(n)
!
   allocate ( rc2(nump+numc+nume) )
!
   pp => getPlane(n)
   do i = 1, nump
      rc2(i) = pp(1,i)*pp(1,i)+pp(2,i)*pp(2,i)+pp(3,i)*pp(3,i)
   enddo
   nullify( pp )
!
   pc => getCorner(n)
   do i = 1, numc
      rc2(i+nump) = pc(1,i)*pc(1,i)+pc(2,i)*pc(2,i)+pc(3,i)*pc(3,i)
   enddo
   nullify( pc )
!
   nc = nump + numc
   pe => getEdge(n)
   j = 0
   do i = 1, nume
      if (.not.isExternalPoint(n,pe(1,i,1),pe(2,i,1),pe(3,i,1))) then
         j = j + 1
         rc2(j+nc) = pe(1,i,1)*pe(1,i,1)+pe(2,i,1)*pe(2,i,1)+         &
                     pe(3,i,1)*pe(3,i,1)
      endif
   enddo
   nullify( pe )
!
   nc = nump + numc + j
!
!  ===================================================================
!  sort the table of possible critical points
!  -------------------------------------------------------------------
   call QuickSort(nc,rc2)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  obtain list of unique critical points
!  ===================================================================
   allocate( rcrit(nc) )
   ncrit=1
   rcrit(1)=rc2(1)
   do i=2,nc
      if(abs(rcrit(ncrit)-rc2(i)) > TEN2m8) then
         ncrit=ncrit+1
         rcrit(ncrit)=rc2(i)
      endif
   enddo
   StepFunction(n)%NumCriticalRs = ncrit
   allocate( StepFunction(n)%CriticalR(ncrit) )
   do i=1,ncrit
      StepFunction(n)%CriticalR(i)=sqrt(rcrit(i))
   enddo
!
   deallocate( rc2, rcrit )
!
!  ===================================================================
!  check for consistency
!  ===================================================================
   if (abs(StepFunction(n)%CriticalR(1)-StepFunction(n)%rmin) > TEN2m6) then
!     ----------------------------------------------------------------
      call ErrorHandler('computeCriticalPoints',                      &
                        'the first critical point is not rmin',       &
                        StepFunction(n)%CriticalR(1),StepFunction(n)%rmin)
!     ----------------------------------------------------------------
   else if (abs(StepFunction(n)%CriticalR(ncrit)-                     &
                StepFunction(n)%rmax) > TEN2m6) then ! change TEn2m10 to TEN2m9
                                                     ! this is a temporary fix
!     ----------------------------------------------------------------
      call ErrorHandler('computeCriticalPoints',                      &
                        'the last critical point is not rmax',        &
                        StepFunction(n)%CriticalR(ncrit),StepFunction(n)%rmax)
!     ----------------------------------------------------------------
   endif
!
   end subroutine computeCriticalPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine stepyll(poly,jmax,r)
!  ===================================================================
   use PolyhedraModule, only : getNumPlanes, getPlane
   use PolyhedraModule, only : getNumEdges
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius
!
!  *******************************************************************
!  *  This subroutine calculates the expansion coefficients,         *
!  *                                    _                            *
!  *  sigma (r),  of stepfunction sigma(r) on complex spherical      *
!  *       L                                                         *
!  *             m _                                                 *
!  *  harmonics Y (r) for given jmax and r                           *
!  *             l                                                   *
!  *              ___                                                *
!  *         _    \                   m _                            *
!  *  sigma( r ) = >   sigma ( r ) * Y (r)                           *
!  *              /__       L         l                              *
!  *              l,m                                                *
!  *                                                                 *
!  *  It calls calsig for calculating the expansion coefficients     *
!  *  of sigma ( r ) through complex spherical harmonics.            *
!  *          L                                                      *
!  *                                                                 *
!  *  We only calculate and store sigma (r) for m >= 0 and {l,m} up  *
!  *                                   L                             *
!  *  to jmax.                                                       *
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: jmax
   integer (kind=IntKind) :: nbnd, lmax
   integer (kind=IntKind) :: ip, np, l, jl, it
   integer (kind=IntKind) :: node
!
   real (kind=RealKind), intent(in) :: r
!
   real (kind=RealKind) :: begth
   real (kind=RealKind) :: endth
   real (kind=RealKind), pointer :: xpt(:,:)
   real (kind=RealKind), allocatable :: xp(:,:)
   real (kind=RealKind), allocatable :: tnode(:)
   real (kind=RealKind), allocatable :: plm(:)
!
   real (kind=RealKind), parameter :: tol = TEN2m10
!
!  ===================================================================
!  special case when r is inside the muffin-tin radius or outside 
!  the bounding sphere radius.
!  ===================================================================
   if (jmax < 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('stepyll','invalid jmax',jmax)
!     ----------------------------------------------------------------
   else if (r < ZERO) then
!     ----------------------------------------------------------------
      call ErrorHandler('stepyll','invalid r',r)
!     ----------------------------------------------------------------
!  else if ( r <= getInscrSphRadius(poly)+TEN2m12 ) then
   else if ( r <= getInscrSphRadius(poly)+TEN2m8 ) then
      wylm(1) = TWO*sqrt(PI)
      if (jmax > 1) then
         wylm(2:jmax) = CZERO
      endif
      return
!  else if ( r >= getOutscrSphRadius(poly)-TEN2m12 ) then
   else if ( r >= getOutscrSphRadius(poly)-TEN2m8 ) then
      wylm(1:jmax) = CZERO
      return
   endif
!
   wylm(1:jmax) = CZERO
!
   nbnd = getNumPlanes(poly)
   xpt => getPlane(poly)
!
   if (sigma(ZERO,ZERO,r,xpt,nbnd,1) == 1) then
      endth=ONE
   else
      endth=TWO
   endif
!
   if(sigma(ZERO,ZERO,-r,xpt,nbnd,1) == 1) then
      begth=-ONE
   else
      begth=-TWO
   endif
!
   allocate( xp(3,nbnd) )
   ip=0
   do np=1,nbnd
      if( (StepFunction(poly)%np_index(1) == np .and.                 &
           r < StepFunction(poly)%rp_index(1))  .or.                  &
          (StepFunction(poly)%np_index(2) == np .and.                 &
           r < StepFunction(poly)%rp_index(2)) ) then
         if (xpt(3,np) > ZERO .and. xpt(3,np) <= r) then
            endth=xpt(3,np)/r
         else if(xpt(3,np) < ZERO .and. -xpt(3,np) <= r) then
            begth=xpt(3,np)/r
         endif
      else if((StepFunction(poly)%np_index(1) /= np) .and.            &
              (StepFunction(poly)%np_index(2) /= np)) then
!        =============================================================
!        only store those boundaries planes not perpendicular to the z-axis.
!        =============================================================
         ip=ip+1
         xp(1,ip)=xpt(1,np)
         xp(2,ip)=xpt(2,np)
         xp(3,ip)=xpt(3,np)
      endif
   enddo
!
   allocate( tnode(nbnd+2*getNumEdges(poly)+2) )
!  -------------------------------------------------------------------
   call caltnode(poly,r,begth,endth,tnode,node)
!  -------------------------------------------------------------------
!
   lmax = lofj(jmax)
!
   if(abs(tnode(1)+ONE) <= tol) then
      allocate( plm(0:lmax) )
!     ----------------------------------------------------------------
      call intpl0(-ONE,tnode(2),lmax,plm)
!     ----------------------------------------------------------------
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         wylm(jl)=wylm(jl)+plm(l)*PI2
      enddo
      deallocate( plm )
      do it=1,node-1
         tnode(it)=tnode(it+1)
      enddo
      node=node-1
   endif
!
   if(abs(tnode(node)-ONE) <= tol .and. node > 1) then
      allocate( plm(0:lmax) )
!     ----------------------------------------------------------------
      call intpl0(tnode(node-1),ONE,lmax,plm)
!     ----------------------------------------------------------------
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         wylm(jl)=wylm(jl)+plm(l)*PI2
      enddo
      deallocate( plm )
      node=node-1
   endif
!
!  ===================================================================
!  it calls calsig to calculate  sigma   ( r ) ,
!                                     l,m
!  which are the expansion coefficients of step function on the 
!  complex spherical harmonics. the data are stored in wylm(jl).
!  -------------------------------------------------------------------
   call calsig(poly,lmax,r,xp,ip,tnode,node)
!  -------------------------------------------------------------------
   do jl=1,jmax
      wylm(jl)=wylm(jl)*clm(jl)
   enddo
   deallocate( xp, tnode )
!
   end subroutine stepyll
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine caltnode(poly,r,begth,endth,tnode,node)
!  ===================================================================
   use PolyhedraModule, only : getNumPlanes, getPlane
   use PolyhedraModule, only : getNumEdges, getEdge
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(out) :: node
   integer (kind=IntKind) :: nbnd, nedge
   integer (kind=IntKind) :: i, n, i1
!
   real (kind=RealKind), intent(in) :: r
   real (kind=RealKind), intent(in) :: begth
   real (kind=RealKind), intent(in) :: endth
   real (kind=RealKind), intent(out) :: tnode(:)
   real (kind=RealKind), pointer :: edget(:,:,:)
   real (kind=RealKind), pointer :: edge(:,:)
   real (kind=RealKind), pointer :: edgp(:,:)
   real (kind=RealKind), pointer :: xp(:,:)
   real (kind=RealKind) :: dn2, dp2, rhop, t, cost2, e2, x0, y0, z0, a, r2
   real (kind=RealKind), parameter :: tol = TEN2m12
!
   if(begth+ONE > -tol) then
      n=1
      tnode(1)=begth
   else
      n=0
   endif
!
   nbnd = getNumPlanes(poly)
   xp => getPlane(poly)
   nedge = getNumEdges(poly)
   edget => getEdge(poly)
   edgp => edget(1:3,1:nedge,1)
   edge => edget(1:3,1:nedge,2)
!
   r2=r*r
   do i=1,nbnd
      rhop=xp(1,i)*xp(1,i)+xp(2,i)*xp(2,i)
      dp2=rhop+xp(3,i)*xp(3,i)
      if (r2 > dp2) then
         if(StepFunction(poly)%np_index(1) /= i .and.                    &
            StepFunction(poly)%np_index(2) /= i) then
            cost2 = sqrt((r2-dp2)*rhop/dp2)
            do i1=1,2
               z0=xp(3,i)+cost2
               if (z0 > r) then
                   z0=r
               else if (z0 < -r) then
                   z0=-r
               endif
               a=sqrt((r2-z0*z0)/rhop)
               x0=xp(1,i)*a
               y0=xp(2,i)*a
               if (sigma(x0,y0,z0,xp,nbnd,1) == 1 .or.                 &
                   sigma(-x0,-y0,z0,xp,nbnd,1) == 1) then
                  n=n+1
                  tnode(n)=z0/r
               endif
               cost2=-cost2
            enddo
         end if
      end if
   enddo
   if (n > 1) then
!     ----------------------------------------------------------------
      call cleanup(tnode,n,tol,endth)
!     ----------------------------------------------------------------
   end if
!
!  ===================================================================
!  look for tnode by solving:
!
!            ->       ->  2    2
!          ( p  + t * e  )  = r
!             i        i
!
!                        ->       ->
!          r * tnode = ( p  + t * e  )
!                         i        i  z
!  ===================================================================
   do i=1,nedge
      dn2=edgp(1,i)*edgp(1,i)+edgp(2,i)*edgp(2,i)+edgp(3,i)*edgp(3,i)
      if (r2 > dn2) then
         e2=edge(1,i)*edge(1,i)+edge(2,i)*edge(2,i)+edge(3,i)*edge(3,i)
         t=sqrt((r2-dn2)/e2)
         do i1=1,2
            x0=edgp(1,i)+t*edge(1,i)
            y0=edgp(2,i)+t*edge(2,i)
            z0=edgp(3,i)+t*edge(3,i)
            if(sigma(x0,y0,z0,xp,nbnd,1) == 1) then
               n=n+1
               tnode(n)=z0/r
            endif
            t=-t
         enddo
      end if
   enddo
!
   if (n > 1) then
!     ----------------------------------------------------------------
      call cleanup(tnode,n,tol,endth)
!     ----------------------------------------------------------------
   end if
   if(endth-one > tol) then
      node=n
   else
      n=n+1
      tnode(n)=endth
      node=n
   endif
!
   nullify( xp, edget, edge, edgp )
!
   end subroutine caltnode
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cleanup(tnode,node,tol,endth)
!  ===================================================================
   use SortModule, only : QuickSort
   implicit none
   integer (kind=IntKind), intent(inout) :: node 
   integer (kind=IntKind) :: n , i
!
   real (kind=RealKind), intent(in) :: tol, endth
   real (kind=RealKind), intent(inout) :: tnode(node)
!
   n=1
!  -------------------------------------------------------------------
   call QuickSort(node,tnode)
!  -------------------------------------------------------------------
   do i=2,node
      if ( abs(tnode(i)-tnode(n)) > tol .and. tnode(i) < endth) then
         n=n+1
         tnode(n)=tnode(i)
      end if
   end do
   node=n
!
   end subroutine cleanup
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine intpl0(x1,x2,lmax,yl)
!  ===================================================================
!
!  *******************************************************************
!                               x2       0
!  calculate the integration:  inte {dx P (x)}
!                               x1       l
!
!  and store in intpl0.
!  *******************************************************************
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: x1
   real (kind=RealKind), intent(in) :: x2
   real (kind=RealKind), intent(out) :: yl(0:lmax)
   real (kind=RealKind) :: p0x1
   real (kind=RealKind) :: p0x2
   real (kind=RealKind) :: p1x1
   real (kind=RealKind) :: p1x2
   real (kind=RealKind) :: p2x1
   real (kind=RealKind) :: p2x2
!
   p0x1=ONE
   p0x2=ONE
   p1x1=x1
   p1x2=x2
   yl(0)=x2-x1
!
   do l=2,lmax+1
      p2x1=((2*l-1)*x1*p1x1-(l-1)*p0x1)/dble(l)
      p2x2=((2*l-1)*x2*p1x2-(l-1)*p0x2)/dble(l)
      yl(l-1)=(p2x2-p2x1-x2*p1x2+x1*p1x1)/dble(l-1)
      p0x1=p1x1
      p0x2=p1x2
      p1x1=p2x1
      p1x2=p2x2
   enddo
!
   end subroutine intpl0
!  ===================================================================
!
!  *******************************************************************  
! 
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine intphi(lmax,cosrth,r,xp,nbnd,sumfi)
!  ===================================================================
   use SortModule, only : QuickSort
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: nbnd
   integer (kind=IntKind) :: mp, ip, m
!
   real (kind=RealKind), intent(in) :: cosrth
   real (kind=RealKind), intent(in) :: r
   real (kind=RealKind), intent(in) :: xp(3,nbnd)
   real (kind=RealKind), allocatable :: ford(:)
   real (kind=RealKind) :: dp2, rhop, sinrth, z, rho
   real (kind=RealKind) :: fh, x, y, a, fi1, fi2, f, phi
!
   real (kind=RealKind), parameter :: HalfPI=HALF*PI
   real (kind=RealKind), parameter :: macherr=TEN2m12
!
   complex (kind=CmplxKind), intent(out) :: sumfi(0:lmax)
   complex (kind=CmplxKind) :: dummy1, dummy2, exp1, exp2
!
   sinrth=sqrt(ONE-cosrth*cosrth)
   z=r*cosrth
   rho=r*sinrth
!
!  ===================================================================
!  find those circle portions within the w-s boundary, and only
!  integrate over phi through them. the value of phi for both
!  ending points of any portion on the boundary are stored in
!  ford.
!  ===================================================================
!
   allocate( ford(2*nbnd+2) )
!
   mp=1
   ford(mp)=ZERO
   do ip=1,nbnd
!     ================================================================
!     rp = the distance of the plane from the origin.
!     if r < rp the curve does not touch the plane.
!     if r = rp the curve touches the plane.
!     if r > rp the curve intersects with the plane.
!     mp = the number of intersection points.
!     ================================================================
      rhop=xp(1,ip)*xp(1,ip)+xp(2,ip)*xp(2,ip)
      dp2=rhop+xp(3,ip)*xp(3,ip)
      if(r*r > dp2 .and. rhop > macherr*macherr) then
         rhop=sqrt(rhop)
         a=(dp2-xp(3,ip)*z)/(rho*rhop)
         if(abs(a) <= ONE+macherr) then
            if(abs(a) > ONE) then
               a=a/abs(a)
            endif
            if(abs(xp(2,ip)) <=  macherr) then
               fi1=halfpi*(one-sign(one,xp(1,ip)))
            else if(xp(2,ip) .ge. ZERO) then
               fi1=acos(xp(1,ip)/rhop)
            else
               fi1=PI2-acos(xp(1,ip)/rhop)
            endif
            fi2=acos(a)
            mp=mp+1
            f=fi1-fi2
            if (f >= ZERO) then
               ford(mp)=f
            else
               ford(mp)=PI2+f
            end if
            mp=mp+1
            f=fi1+fi2
            if (f < PI2) then
               ford(mp)=f
            else
               ford(mp)=f-PI2
            endif
         endif
      endif
   enddo
   mp=mp+1
   ford(mp)=PI2
!
!  ===================================================================
!  sort ford so that :   ford(1) < ford(2) < ... < ford(mp).
!  -------------------------------------------------------------------
   call QuickSort(mp,ford)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  start the phi integration.
!  ===================================================================
   do m=0,lmax
      sumfi(m)=CZERO
   enddo
   exp1=exp(-SQRTm1*ford(1))
   exp2=CONE
   do ip=1,mp-1
      phi=ford(ip+1)-ford(ip)
      exp1=exp1*conjg(exp2)
      exp2=exp(SQRTm1*HALF*phi)
      exp1=exp1*conjg(exp2)
!
!     ================================================================
!     sigma:    = 1       if curve between ford(i) and
!                         ford(i+1) is within the cell.
!               = 0       otherwise.
!
!     i = 1,2,...,mp-1
!     ================================================================
      fh=HALF*(ford(ip)+ford(ip+1))
      x=rho*cos(fh)
      y=rho*sin(fh)
!
!     ================================================================
!     check if the point (x,y,z) is inside or outside the polyhedron.
!     ================================================================
!
      if(sigma(x,y,z,xp,nbnd,0) == 1) then
         sumfi(0)=sumfi(0)+phi
         dummy1=CONE
         dummy2=CONE
         do m=1,lmax
            dummy1=dummy1*exp1
            dummy2=dummy2*exp2
            sumfi(m) = sumfi(m) + dummy1*aimag(dummy2)
         enddo
      endif  ! sigma=1
   enddo
   do m=1,lmax
      sumfi(m)=TWO*sumfi(m)/m
   enddo
!
   deallocate( ford )
!
   end subroutine intphi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function sigma(x0,y0,z0,xp,nbnd,m) result(s)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: nbnd
   integer (kind=IntKind), intent(in) :: m
   integer (kind=IntKind) :: s
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: x0
   real (kind=RealKind), intent(in) :: y0
   real (kind=RealKind), intent(in) :: z0
   real (kind=RealKind), intent(in) :: xp(3,nbnd)
   real (kind=RealKind) :: tol
!
   tol = TEN2m10*m
!
   s=1
   do i=1,nbnd
      if((x0-xp(1,i))*xp(1,i)+(y0-xp(2,i))*xp(2,i)+(z0-xp(3,i))*xp(3,i) &
         > tol) then
         s=0
         return
      end if
   end do
!
   end function sigma
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calsig(poly,lmax,r,xp,nbnd,tnode,node)
!  ===================================================================
!
!  *******************************************************************
!  *                                                                 *
!  *  This is the new program to expand sigma(r) on complex spherical*
!  *  harmonics. It only requires the position of boundary planes.   *
!  *                                                                 *
!  *  the integration :                                              *
!  *                     /                                           *
!  *                     |  _          _      *   _                  *
!  *                     | do * sigma( r ) * y  ( r ) = wylm         *
!  *                     |                    l,m                    *
!  *                     /                                           *
!  *                                                                 *
!  *  for certain r and m => 0, is carried through this subroutine.  *
!  *                                                                 *
!  *                           m        *                            *
!  *  For m < 0, w(l,-m) = (-1) * w(l,m)                             *
!  *                                                                 *
!  *                             *                                   *
!  *                     = w(l,m)         if m = even                *
!  *                                                                 *
!  *  The array is stored in the following way:                      *
!  *                                                                 *
!  *      jl       l       m              jl      l      m           *
!  *                                                                 *
!  *    w( 1 )     0       0            w( 7 )    3      0           *
!  *    w( 2 )     1       0            w( 8 )    3      1           *
!  *    w( 3 )     1       1            w( 9 )    3      2           *
!  *    w( 4 )     2       0            w(10 )    3      3           *
!  *    w( 5 )     2       1               .      .      .           *
!  *    w( 6 )     2       2               .      .      .           *
!  *                                       .      .      .           *
!  *                                                                 *
!  *  input:                                                         *
!  *                                                                 *
!  *        lmax     = the maximum value of l.                       *
!  *        xgq      = the Gaussian quadrature points.               *
!  *        wgq      = the Gaussian quadrature weight associating    *
!  *                   each value of xgq.                            *
!  *        npts     = the number of Gaussian qurdrature points.     *
!  *        r        = the ratial value.                             *
!  *        xp       = the vector normal to a boundary plane and     *
!  *                   ending on the plane (from origin).            *
!  *        nbnd     = the number of boundary planes.                *
!  *        tnode    = the nodes along theta integration which is the*
!  *                   non-regular points of function: f(theta).     *
!  *                   The integration over theta is broken into     *
!  *                   several pieces, each of which has two nodes as*
!  *                   the terminal points.                          *
!  *        node     = the number of tnode.                          *
!  *                                                                 *
!  *                                                                 *
!  *  output:                                                        *
!  *                                                                 *
!  *        wylm     = step function expansion value for current     *
!  *                   value of r.                                   *
!  *                                                                 * 
!  *******************************************************************
!
   use LegendreModule, only : legendre
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(in) :: nbnd
   integer (kind=IntKind), intent(in) :: node
   integer (kind=IntKind) :: jmax, i, l, m, jl, n, npts
!
   real (kind=RealKind), intent(in) :: tnode(node)
   real (kind=RealKind), intent(in) :: xp(3,nbnd)
   real (kind=RealKind), allocatable :: plm(:)
   real (kind=RealKind), pointer :: xgq(:), wgq(:)
   real (kind=RealKind) :: wt, cosrth, offset, d, halfd, h0, sh0, r, u, t(2)
   real (kind=RealKind), parameter :: wtol = TEN2m12
!
   complex (kind=CmplxKind), allocatable :: sumfi(:,:)
   complex (kind=CmplxKind) :: af1
!
!  ===================================================================
!  start to integrate over xn=cos(theta), using gaussian method. 
!
!             2pi              ->   -i*m*phi
!     sumfi = int d(phi) sigma(r )*e        ,   for certain cos(theta)
!              0
!
!  ===================================================================
!
   jmax = ((lmax+1)*(lmax+2))/2
   allocate( plm(jmax), sumfi(0:lmax,3) )
!
   npts = StepFunction(poly)%NumGaussCosThetas
   xgq => StepFunction(poly)%GaussCosTheta(1:npts)
   wgq => StepFunction(poly)%GaussCosThetaWeight(1:npts)
!
   do i=1,node-1
      d=tnode(i+1)-tnode(i)
      h0=d*TEN2m6
      if(h0 < TEN2m8) then
         h0=min(TEN2m8,d*TEN2m3)
      endif
      sh0=sqrt(h0)
!     ----------------------------------------------------------------
      call intphi(lmax,tnode(i),r,xp,nbnd,sumfi(0,3))
      call intphi(lmax,tnode(i)+h0,r,xp,nbnd,sumfi(0,1))
!     ----------------------------------------------------------------
      do m=0,lmax
         sumfi(m,1)=(sumfi(m,1)-sumfi(m,3))/sh0
      enddo
!
!     ----------------------------------------------------------------
      call intphi(lmax,tnode(i+1),r,xp,nbnd,sumfi(0,3))
      call intphi(lmax,tnode(i+1)-h0,r,xp,nbnd,sumfi(0,2))
!     ----------------------------------------------------------------
      do m=0,lmax
         sumfi(m,2)=(sumfi(m,2)-sumfi(m,3))/sh0
      enddo
!
      offset=HALF*(tnode(i)+tnode(i+1))
      halfd=HALF*d
      do n=1,npts
         cosrth=offset+halfd*xgq(n)
!
!        =============================================================
!        call legendre to generate the Legendre functions up to
!        l = lmax for each cos(theta) value.
!        -------------------------------------------------------------
         call legendre(lmax,cosrth,plm)
!        -------------------------------------------------------------
         call intphi(lmax,cosrth,r,xp,nbnd,sumfi(0,3))
!        -------------------------------------------------------------
!
         t(1)=sqrt(halfd+halfd*xgq(n))
         t(2)=sqrt(halfd-halfd*xgq(n))
         wt=wgq(n)*halfd
         do m=0,lmax
            af1= wt*(sumfi(m,3)-(sumfi(m,1)*t(1)+sumfi(m,2)*t(2)))
            jl=(m*(m+1))/2+1
            do l=m,lmax
               jl=jl+l
               wylm(jl) = wylm(jl)+af1*plm(jl)
            enddo
         enddo
!
         if(xgq(n) >= ZERO) then
            u=sqrt(d)*xgq(n)
            wt=two*sqrt(d)*wgq(n)*u*u
!           ----------------------------------------------------------
            call legendre(lmax,u*u+tnode(i),plm)
!           ----------------------------------------------------------
            do m=0,lmax
               af1=wt*sumfi(m,1)
               jl=(m*(m+1))/2+1
               do l=m,lmax
                  jl=jl+l
                  wylm(jl)=wylm(jl)+af1*plm(jl)
               enddo
            enddo
!           ----------------------------------------------------------
            call legendre(lmax,tnode(i+1)-u*u,plm)
!           ----------------------------------------------------------
            do m=0,lmax
               af1=wt*sumfi(m,2)
               jl=(m*(m+1))/2+1
               do l=m,lmax
                  jl=jl+l
                  wylm(jl)=wylm(jl)+af1*plm(jl)
               enddo
            enddo
         endif
      enddo
   enddo
!
   do jl=1,jmax
      if (abs(wylm(jl)) < wtol) then
         wylm(jl)=CZERO
      endif
   enddo
!
   deallocate( plm, sumfi )
   nullify( xgq, wgq )
!
   end subroutine calsig
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getRadialStepFunction(id, mr, nr, r, lmax, sigma)
!  ===================================================================
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, mr, nr, lmax
!
   real (kind=RealKind), target :: r(nr)
!
   complex (kind=CmplxKind), intent(inout) :: sigma(:,:)
!
   integer (kind=IntKind) :: ir, mrt, jl
   integer (kind=IntKind) :: lmax_sigma, jmax_sigma
   integer (kind=IntKind) :: lmax_step, jmax_step
!
!  mr = nr
!  LOOP_ir: do ir = nr-1, 1, -1
!     if ( r(ir) < StepFunction(id)%rmin  - TEN2m8 ) then
!        mr = ir+1
!        exit LOOP_ir
!     endif
!  enddo LOOP_ir
!
   lmax_sigma = lmax
   lmax_step = min(lmax,StepFunction(id)%lmax_step)
   jmax_step = ((lmax_step+1)*(lmax_step+2))/2
   jmax_sigma = ((lmax_sigma+1)*(lmax_sigma+2))/2
!
!  if ( abs(StepFunction(id)%rmin-r(mr)) > TEN2m8 ) then
!     call WarningHandler("StepFunction:: getRadialStepFunction:",                 &
!                         "rmin doesn't fall on the radial mesh ")
!  endif
!
   mrt = nr-mr+1
!
   if ( mr == nr ) then
      sigma = CZERO ! nothing to truncate
      return
   else if (mr > nr) then
      call ErrorHandler('StepFunction:: getRadialStepFunction:','mr > nr',mr,nr)
   else if (mrt > size(sigma,1)) then
      call ErrorHandler('StepFunction:: getRadialStepFunction:','mrt > 1st dim size of sigma',mrt,size(sigma,1))
   endif
!
   sigma = CZERO
   do ir = mr,nr
      sigma(ir-mr+1,1) = getStepFunction(id,1,r(ir))
      if ( jmax_step>1 .and. jmax_sigma>=2 ) then
         do jl = 2,min(jmax_sigma,jmax_step)
            sigma(ir-mr+1,jl) = wylm(jl)
         enddo
      endif
   enddo
!
   end subroutine getRadialStepFunction
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine truncate_symm(id, mr, nr, r, f, fsymm, jl_f, nrt, ftrnc, jl_ft)
!  ===================================================================
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: mr, nr, nrt
   integer (kind=IntKind), intent(in) :: jl_f, jl_ft
   integer (kind=IntKind), intent(in) :: fsymm(jl_f)
!
   real (kind=RealKind), target :: r(nr)
!
   complex (kind=CmplxKind), intent(in)  :: f(nr,jl_f)
   complex (kind=CmplxKind), intent(out) :: ftrnc(1:nrt,1:jl_ft)
!
   integer (kind=IntKind) :: ir, mrt, izamax
   integer (kind=IntKind) :: lmax_sigma, jmax_sigma, kmax_sigma
   integer (kind=IntKind) :: lmax_step, jmax_step, lmax_f, lmax_ft
   integer (kind=IntKind) :: l, ml, jl, kl, lp, mlp, jlp, klp
   integer (kind=IntKind) :: mlpp, jlpp, klpp, lpp, j, kmax_f, jmax_f
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!
   complex (kind=CmplxKind) :: sigma_lpp
   complex (kind=CmplxKind) :: cfac, cfacp, f_ir
!
   complex (kind=CmplxKind), allocatable :: sigma(:,:)
   complex (kind=CmplxKind), allocatable :: sum_i3(:), sum_kl(:)
!
!  mr = nr
!  LOOP_ir: do ir = nr-1, 1, -1
!     if ( r(ir) < StepFunction(id)%rmin  - TEN2m8 ) then
!        mr = ir+1
!        exit LOOP_ir
!     endif
!  enddo LOOP_ir
!
   mrt = nr-mr+1
!
   if ( mr == nr ) then
      ftrnc = CZERO ! nothing to truncate
      return
   else if (mr > nr) then
      call ErrorHandler('StepFunction:: truncate:','mr > nr',mr,nr)
   else if (mrt > nrt) then
      call WarningHandler('StepFunction:: truncate:','mrt > nrt',mrt,nrt)
      mrt = nrt
   endif
!
   lmax_f  = lofj(jl_f)
   lmax_ft = lofj(jl_ft)
   jmax_f = (lmax_f+1)*(lmax_f+2)/2
   kmax_f  = (lmax_f+1)*(lmax_f+1)
   lmax_step = lmax_f+lmax_ft
   lmax_sigma = min(lmax_step,StepFunction(id)%lmax_step)
   jmax_step = ((lmax_step+1)*(lmax_step+2))/2
   jmax_sigma = ((lmax_sigma+1)*(lmax_sigma+2))/2
   kmax_sigma = (lmax_sigma+1)*(lmax_sigma+1)
!
!  if ( abs(StepFunction(id)%rmin-r(mr)) > TEN2m8 ) then
!     call WarningHandler("StepFunction:: truncate:",                 &
!                         "rmin doesn't fall on the radial mesh ")
!  endif
!
   allocate( sigma(mrt,jmax_step) )
   allocate( sum_i3(mrt), sum_kl(mrt))
!
   sigma = CZERO
   do ir = mr,nr
      sigma(ir-mr+1,1) = getStepFunction(id,1,r(ir))
      if ( jmax_step>1 .and. jmax_sigma>=2 ) then
         do jl = 2,jmax_sigma
            sigma(ir-mr+1,jl) = wylm(jl)
         enddo
      endif
   enddo
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   ftrnc = CZERO
!
   LOOP_jl: do jl = 1,jl_ft
      if ( jl <= min(jmax_f,jmax_sigma) ) then
         if (fsymm(jl) == 0 .and. StepFunction(id)%SigmaCompFlag(jl) == 0) cycle LOOP_jl
      else if (jl <= StepFunction(id)%jmax_step) then
         if (StepFunction(id)%SigmaCompFlag(jl) == 0) cycle LOOP_jl
      endif
      l  = lofj(jl)
      ml = mofj(jl)
      kl = (l+1)*(l+1) - l + ml
      sum_kl(1:mrt) = CZERO
      do klp = 1,kmax_f
         lp = lofk(klp)
         mlp = mofk(klp)
         jlp = jofk(klp)
         sum_i3(1:mrt) = CZERO
         Loop_i3: do j=nj3(klp,kl),1,-1
            klpp = kj3(j,klp,kl)
            if ( klpp>kmax_sigma ) cycle Loop_i3
            lpp = lofk(klpp)
            cfac=cgnt(j,klp,kl)
            jlpp = jofk(klpp)
            if (StepFunction(id)%SigmaCompFlag(jlpp) == 0) cycle Loop_i3
            mlpp = mofk(klpp)
            cfacp=cfac*m1m(mlpp)
            do ir = mr,mr+mrt-1
               if ( mlp>=0 ) then
                  f_ir = f(ir,jlp)
               else
                  f_ir = m1m(mlp)*conjg(f(ir,jlp))
               endif
               if ( mlpp>=0 ) then
                  sigma_lpp = cfac*sigma(ir-mr+1,jlpp)
               else
                  sigma_lpp = cfacp*conjg(sigma(ir-mr+1,jlpp))
               endif
               sum_i3(ir-mr+1) = sum_i3(ir-mr+1) + f_ir*sigma_lpp
            enddo
         enddo Loop_i3
         do ir = 1,mrt
            sum_kl(ir) = sum_kl(ir) + sum_i3(ir)
         enddo
      enddo
      do ir = 1,mrt
         ftrnc(ir,jl) = sum_kl(ir)
      enddo
   enddo LOOP_jl
!
   deallocate( sigma, sum_i3, sum_kl )
!
   end subroutine truncate_symm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine truncate_all(id, mr, nr, r, f, jl_f, nrt, ftrnc, jl_ft)
!  ===================================================================
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: mr, nr, nrt
   integer (kind=IntKind), intent(in) :: jl_f, jl_ft
!
   real (kind=RealKind), target :: r(nr)
!
   complex (kind=CmplxKind), intent(in)  :: f(nr,jl_f)
   complex (kind=CmplxKind), intent(out) :: ftrnc(1:nrt,1:jl_ft)
!
   integer (kind=IntKind) :: ir, mrt, izamax
   integer (kind=IntKind) :: lmax_sigma, jmax_sigma, kmax_sigma
   integer (kind=IntKind) :: lmax_step, jmax_step, lmax_f, lmax_ft
   integer (kind=IntKind) :: l, ml, jl, kl, lp, mlp, jlp, klp
   integer (kind=IntKind) :: mlpp, jlpp, klpp, lpp, j, kmax_f
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!
   complex (kind=CmplxKind) :: sigma_lpp
   complex (kind=CmplxKind) :: cfac, cfacp, f_ir
!
   complex (kind=CmplxKind), allocatable :: sigma(:,:)
   complex (kind=CmplxKind), allocatable :: sum_i3(:), sum_kl(:)
!
!  mr = nr
!  LOOP_ir: do ir = nr-1, 1, -1
!     if ( r(ir) < StepFunction(id)%rmin  - TEN2m8 ) then
!        mr = ir+1
!        exit LOOP_ir
!     endif
!  enddo LOOP_ir
!
   mrt = nr-mr+1
!
   if ( mr == nr ) then
      ftrnc = CZERO ! nothing to truncate
      return
   else if (mr > nr) then
      call ErrorHandler('StepFunction:: truncate:','mr > nr',mr,nr)
   else if (mrt > nrt) then
      call WarningHandler('StepFunction:: truncate:','mrt > nrt',mrt,nrt)
      mrt = nrt
   endif
!
   lmax_f  = lofj(jl_f)
   lmax_ft = lofj(jl_ft)
   kmax_f  = (lmax_f+1)*(lmax_f+1)
   lmax_step = lmax_f+lmax_ft
   lmax_sigma = min(lmax_step,StepFunction(id)%lmax_step)
   jmax_step = ((lmax_step+1)*(lmax_step+2))/2
   jmax_sigma = ((lmax_sigma+1)*(lmax_sigma+2))/2
   kmax_sigma = (lmax_sigma+1)*(lmax_sigma+1)
!
!  if ( abs(StepFunction(id)%rmin-r(mr)) > TEN2m8 ) then
!     call WarningHandler("StepFunction:: truncate:",                 &
!                         "rmin doesn't fall on the radial mesh ")
!  endif
!
   allocate( sigma(mrt,jmax_step) )
   allocate( sum_i3(mrt), sum_kl(mrt))
!
   sigma = CZERO
   do ir = mr,nr
      sigma(ir-mr+1,1) = getStepFunction(id,1,r(ir))
      if ( jmax_step>1 .and. jmax_sigma>=2 ) then
         do jl = 2,jmax_sigma
            sigma(ir-mr+1,jl) = wylm(jl)
         enddo
      endif
   enddo
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   ftrnc = CZERO
!
   do jl = 1,jl_ft
      l  = lofj(jl)
      ml = mofj(jl)
      kl = (l+1)*(l+1) - l + ml
      sum_kl(1:mrt) = CZERO
      do klp = 1,kmax_f
         lp = lofk(klp)
         mlp = mofk(klp)
         jlp = jofk(klp)
         sum_i3(1:mrt) = CZERO
         Loop_i3: do j=nj3(klp,kl),1,-1
            klpp = kj3(j,klp,kl)
            if ( klpp>kmax_sigma ) cycle Loop_i3
            lpp = lofk(klpp)
            cfac=cgnt(j,klp,kl)
            jlpp = jofk(klpp)
            if (StepFunction(id)%SigmaCompFlag(jlpp) == 0) cycle Loop_i3
            mlpp = mofk(klpp)
            cfacp=cfac*m1m(mlpp)
            do ir = mr,mr+mrt-1
               if ( mlp>=0 ) then
                  f_ir = f(ir,jlp)
               else
                  f_ir = m1m(mlp)*conjg(f(ir,jlp))
               endif
               if ( mlpp>=0 ) then
                  sigma_lpp = cfac*sigma(ir-mr+1,jlpp)
               else
                  sigma_lpp = cfacp*conjg(sigma(ir-mr+1,jlpp))
               endif
               sum_i3(ir-mr+1) = sum_i3(ir-mr+1) + f_ir*sigma_lpp
            enddo
         enddo Loop_i3
         do ir = 1,mrt
            sum_kl(ir) = sum_kl(ir) + sum_i3(ir)
         enddo
      enddo
      do ir = 1,mrt
         ftrnc(ir,jl) = sum_kl(ir)
      enddo
   enddo
!
   deallocate( sigma, sum_i3, sum_kl )
!
   end subroutine truncate_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine truncate_jl(id, mr, nr, r, f, jl_f, nrt, ftrnc, jl_ft)
!  ===================================================================
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: mr, nr, nrt
   integer (kind=IntKind), intent(in) :: jl_f, jl_ft
!
   real (kind=RealKind), target :: r(nr)
!
   complex (kind=CmplxKind), intent(in)  :: f(nr,jl_f)
   complex (kind=CmplxKind), intent(out) :: ftrnc(1:nrt)
!
   integer (kind=IntKind) :: ir, mrt
   integer (kind=IntKind) :: lmax_sigma, jmax_sigma, kmax_sigma
   integer (kind=IntKind) :: lmax_step, jmax_step, lmax_f, lmax_ft
   integer (kind=IntKind) :: l, ml, jl, kl, lp, mlp, jlp, klp
   integer (kind=IntKind) :: mlpp, jlpp, klpp, lpp, j, kmax_f
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!
   complex (kind=CmplxKind) :: sigma_lpp
   complex (kind=CmplxKind) :: cfac, cfacp, f_ir
!
   complex (kind=CmplxKind), allocatable :: sigma(:,:)
   complex (kind=CmplxKind), allocatable :: sum_i3(:)
!
!  mr = nr
!  LOOP_ir: do ir = nr-1, 1, -1
!     if ( r(ir) < StepFunction(id)%rmin  - TEN2m8 ) then
!        mr = ir+1
!        exit LOOP_ir
!     endif
!  enddo LOOP_ir
!
   mrt = nr-mr+1
!
   if ( mr == nr ) then
      ftrnc = CZERO             ! nothing to truncate
      return
   else if (mr > nr) then
      call ErrorHandler('StepFunction:: truncate:','mr > nr',mr,nr)
   else if (mrt > nrt) then
      call WarningHandler('StepFunction:: truncate:','mrt > nrt',mrt,nrt)
      mrt = nrt
   endif
!
   lmax_f  = lofj(jl_f)
   lmax_ft = lofj(jl_ft)
   kmax_f  = (lmax_f+1)*(lmax_f+1)
   lmax_step = lmax_f+lmax_ft
   lmax_sigma = min(lmax_step,StepFunction(id)%lmax_step)
   jmax_step = ((lmax_step+1)*(lmax_step+2))/2
   jmax_sigma = ((lmax_sigma+1)*(lmax_sigma+2))/2
   kmax_sigma = (lmax_sigma+1)**2
!
!  if ( abs(StepFunction(id)%rmin-r(mr)) > TEN2m8 ) then
!     call WarningHandler("StepFunction:: truncate:",                 &
!                         "rmin doesn't fall on the radial mesh ")
!  endif
!
   allocate( sigma(mrt,jmax_step) )
   allocate( sum_i3(mrt) )
!
   sigma = CZERO
   do ir = mr,nr
      sigma(ir-mr+1,1) = getStepFunction(id,1,r(ir))
      if ( jmax_step>1 .and. jmax_sigma>=2 ) then
         do jl = 2,jmax_sigma
            sigma(ir-mr+1,jl) = wylm(jl)
         enddo
      endif
   enddo
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   l  = lofj(jl_ft)
   ml = mofj(jl_ft)
   kl = (l+1)*(l+1) - l + ml
!
   ftrnc = CZERO
!
   do klp = 1,kmax_f
      lp = lofk(klp)
      mlp = mofk(klp)
      jlp = jofk(klp)
      sum_i3(1:mrt) = CZERO
      Loop_i3: do j=nj3(klp,kl),1,-1
         klpp = kj3(j,klp,kl)
         if ( klpp>kmax_sigma ) cycle Loop_i3
         lpp = lofk(klpp)
         cfac=cgnt(j,klp,kl)
         jlpp = jofk(klpp)
         mlpp = mofk(klpp)
         cfacp=cfac*m1m(mlpp)
         do ir = mr,mr+mrt-1
            if ( mlp>=0 ) then
               f_ir = f(ir,jlp)
            else
               f_ir = m1m(mlp)*conjg(f(ir,jlp))
            endif
            if ( mlpp>=0 ) then
               sigma_lpp = cfac*sigma(ir-mr+1,jlpp)
            else
               sigma_lpp = cfacp*conjg(sigma(ir-mr+1,jlpp))
            endif
            sum_i3(ir-mr+1) = sum_i3(ir-mr+1) + f_ir*sigma_lpp
         enddo
      enddo Loop_i3
      do ir = 1,mrt
         ftrnc(ir) = ftrnc(ir) + sum_i3(ir)
      enddo
   enddo
!
   deallocate( sigma, sum_i3 )
!
   end subroutine truncate_jl
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegrationToR(m,nr,r,k,fr,rmax,mr,g)
!  ===================================================================
   use InterpolationModule, only : getInterpolation
!
   use IntegrationModule, only : calIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: m, nr, k, mr
!
   real(kind=RealKind), intent(in) :: rmax
   real(kind=RealKind), intent(in) :: r(1:nr)
!
   complex (kind=CmplxKind), intent(in) :: fr(0:nr)
   complex (kind=CmplxKind), intent(out) :: g(0:nr)
!
   integer (kind=IntKind) :: ir, k_int
!
   real (kind=RealKind) :: err, rt, fac
   real(kind=RealKind), pointer :: rr(:)
!
   complex (kind=CmplxKind) :: fg, dps
!
   if ( m == 0 ) then
      k_int = 1
      rr => sqrt_r(0:nr+1)
      fac = TWO
   else
      k_int = 0
      rr => rr0(0:nr)
      fac = ONE
   endif
!
   if ( abs(rmax-r(mr) ) < TEN2m12 ) then
      if ( k==2 ) then
!        -------------------------------------------------------------
         call FitInterp( 4, rr(2:5), fr2(1:4), ZERO, fr2(0), dps )
!        ------------------------------------------------------------
         call calIntegration(mr+1,rr(1:mr+1),fr(0:mr),g(0:mr),k_int)
!        ------------------------------------------------------------
      else
         do ir = 1, mr
            fr2(ir) = fr(ir)*r(ir)**(2-k)
         enddo
!        -------------------------------------------------------------
         call FitInterp( 4, rr(2:5), fr2(1:4), ZERO, fr2(0), dps )
!        -------------------------------------------------------------
         call calIntegration(mr+1,rr(1:mr+1),fr2(0:mr),g(0:mr),k_int)
!        -------------------------------------------------------------
      endif
      g(0:mr)= fac*g(0:mr)
   else
      rt = rr(mr+2)
      if ( m==0 ) then
         rr(mr+2) = sqrt(rmax)
      else
         rr(mr+2) = rmax
      endif
!     ----------------------------------------------------------------
      fg = getInterpolation(nr,rr(1:nr+1),fr(1:nr),rr(mr+2),err)
!     ----------------------------------------------------------------
      if ( k==2 ) then
         do ir = 1, mr
            fr2(ir) = fr(ir)
         enddo
         fr2(mr+1) = fg
!        -------------------------------------------------------------
         call FitInterp( 4, rr(2:5), fr2(1:4), ZERO, fr2(0), dps )
!        ------------------------------------------------------------
         call calIntegration(mr+2,rr(1:mr+2),fr2(0:mr+1),g(0:mr+1),k_int)
!        ------------------------------------------------------------
      else
         do ir = 1, mr
            fr2(ir) = fr(ir)*r(ir)**(2-k)
         enddo
         fr2(mr+1) = fg*rmax**(2-k)
!        -------------------------------------------------------------
         call FitInterp( 4, rr(2:5), fr2(1:4), ZERO, fr2(0), dps )
!        ------------------------------------------------------------
         call calIntegration(mr+2,rr(1:mr+2),fr2(0:mr+1),g(0:mr+1),k_int)
!        ------------------------------------------------------------
      endif
      g(0:mr+1)= fac*g(0:mr+1)
      rr(mr+2) = rt
   endif
!
   end subroutine calIntegrationToR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupIntWkSpace(nr,r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr
   real (kind=RealKind), intent(in) :: r(nr)
!
   integer (kind=IntKind) :: ir
!
   if (max_nr>=nr+1) then
      sqrt_r(0) = ZERO
      rr0(0)    = ZERO
      do ir = 1,nr
         sqrt_r(ir) = sqrt(r(ir))
         rr0(ir)    = r(ir)
      enddo
   else
      if ( max_nr/=0 ) then
         deallocate( sqrt_r, fr2 )
         deallocate( res, rfg, ces, cfg, rr0 )
      endif
      allocate( sqrt_r(0:nr+1), fr2(0:nr+1), rr0(0:nr+1) )
      allocate( res(0:nr+1), ces(0:nr+1),  rfg(0:nr+1), cfg(0:nr+1) )
      sqrt_r(0) = ZERO
      rr0(0)    = ZERO
      do ir = 1,nr
         sqrt_r(ir)=sqrt(r(ir))
         rr0(ir)    = r(ir)
      enddo
      max_nr = nr+1
   endif
!
   fr2(0:nr+1) = CZERO
   cfg(0:nr+1) = CZERO
   ces(0:nr+1) = CZERO
   res(0:nr+1) = ZERO
   rfg(0:nr+1) = ZERO
!
   end subroutine setupIntWkSpace
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printStepFunction(poly)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: j, jl
!
   type (StepFunctionStruct), pointer :: p
!
   if (.not.Initialized) then
      call ErrorHandler('printStepFunction',                          &
                        'StepFunction module not initialized')
   else if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('printStepFunction','invalid polyhedron index',poly)
   endif
!
   p => StepFunction(poly)
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')'*********************************'
   write(6,'( 24x,a )')'* Output from printStepFunction *'
   write(6,'(24x,a,/)')'*********************************'
   write(6,'(a,i5)')   'Polyhedron Index:                 ',poly
   write(6,'(a,i5)')   'Lmax for the Step Function:       ',p%lmax_step
   write(6,'(a,i5)')   'Number of Critical Points:        ',p%NumCriticalRs
   write(6,'(a,i5)')   'Number of Radial Gaussian Points: ',p%NumGaussRs
   write(6,'(a,i5)')   'Number of Theta  Gaussian Points: ',p%NumGaussCosThetas
   write(6,'(a,f10.5)')   'Inner Bounding Sphere Radius:  ',p%rmin
   write(6,'(a,f10.5)')   'Outer Bounding Sphere Radius:  ',p%rmax
!
   write(6,'(/,a)')'Non-zero (l,m) components of the step function'
   j = 0
   do jl = 1, StepFunction(poly)%jmax_step
      if ( StepFunction(poly)%SigmaCompFlag(jl) /= 0 ) then
         if (j == 0) then 
            write(6,'(a,i2,a,i2,a,$)')'(',lofj(jl),',',mofj(jl),')'
            j = 7
         else if (j < 88) then
            write(6,'(a,i2,a,i2,a,$)')', (',lofj(jl),',',mofj(jl),')'
            j = j + 9
         else
            write(6,'(a)')','
            write(6,'(a,i2,a,i2,a,$)')'(',lofj(jl),',',mofj(jl),')'
            j = 7
         endif
      endif
   enddo
   write(6,'(a)')' '
!
   write(6,'(/,40(''=''))')
   write(6,'(a)')'Critical Radial Point Table'
   write(6,'(40(''=''))')
   write(6,'(a)')'Critical R Point         Radius'
   write(6,'(40(''-''))')
   do j=1,p%NumCriticalRs
      write(6,'(4x,i5,13x,f10.5)')j,p%CriticalR(j)
   enddo
   write(6,'(40(''=''))')
!
   nullify( p )
!
   end subroutine printStepFunction
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine testStepFunction_i(poly)
!  ===================================================================
   use PolyhedraModule, only : getVolume
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ig, ng
!
   real (kind=RealKind) :: exactV, calcV
   real (kind=RealKind), pointer :: rg(:)
   real (kind=RealKind), pointer :: wg(:)
!
   complex (kind=CmplxKind) :: vint
!
   if (.not.Initialized) then
      call ErrorHandler('printStepFunction',                          &
                        'StepFunction module not initialized')
   else if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('printStepFunction','invalid polyhedron index',poly)
   endif
!
   ng = StepFunction(poly)%NumGaussRs
   rg => StepFunction(poly)%GaussR(1:ng)
   wg => StepFunction(poly)%GaussRWeight(1:ng)
!
   vint = CZERO
   do ig = 1, ng
      vint = vint + wg(ig)*StepFunction(poly)%sigma_L(ig,1)*rg(ig)*rg(ig)
   enddo
   nullify(rg, wg)
!
   calcV = PI4*(StepFunction(poly)%rmin)**3/THREE + real(vint,RealKind)*Y0inv
   exactV = getVolume(poly)
!
   if (print_level >= 0) then
      write(6,'(/,80(''-''))')
      write(6,'(/,23x,a)')'**********************************'
      write(6,'( 23x,a )')'*  Output from testStepFunction  *'
      write(6,'(23x,a,/)')'**********************************'
      write(6,'(a,i5)')   'Polyhedron Index:       ',poly
      write(6,'(a,d22.12)')'Exact Cell Volume:      ',exactV
      write(6,'(a,d22.12)')'Calculated Cell Volume: ',calcV
      write(6,'(80(''=''))')
   endif
!
   end subroutine testStepFunction_i
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine testStepFunction_t()
!  ===================================================================
   use PolyhedraModule, only : getVolume
   implicit none
!
   integer (kind=IntKind) :: poly
   integer (kind=IntKind) :: ig, ng
!
   real (kind=RealKind) :: exactV, calcV
   real (kind=RealKind), pointer :: rg(:)
   real (kind=RealKind), pointer :: wg(:)
!
   complex (kind=CmplxKind) :: vint
!
   if (.not.Initialized) then
      call ErrorHandler('printStepFunction',                          &
                        'StepFunction module not initialized')
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,23x,a)')'**********************************'
   write(6,'( 23x,a )')'*  Output from testStepFunction  *'
   write(6,'(23x,a,/)')'**********************************'
!
   write(6,'(a)')'============================================================='
   write(6,'(a)')'Polyhedron Index       Exact Volume         Calculated Volume'
   write(6,'(a)')'-------------------------------------------------------------'
   do poly = 1, NumPolyhedra
      ng = StepFunction(poly)%NumGaussRs
      rg => StepFunction(poly)%GaussR(1:ng)
      wg => StepFunction(poly)%GaussRWeight(1:ng)
!
      vint = CZERO
      do ig = 1, ng
         vint = vint + wg(ig)*StepFunction(poly)%sigma_L(ig,1)*rg(ig)*rg(ig)
      enddo
      nullify(rg, wg)
!
      calcV = PI4*(StepFunction(poly)%rmin)**3/THREE + real(vint,RealKind)*Y0inv
      exactV = getVolume(poly)
      write(6,'(5x,i5,8x,d20.12,3x,d20.12)')poly,exactV,calcV
   enddo
   write(6,'(a)')'============================================================='
!
   end subroutine testStepFunction_t
!  ===================================================================
end module StepFunctionModule
