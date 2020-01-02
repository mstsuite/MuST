module RelGreenFunctionModule

!********************************************************************************
! Name: RelGreenFunctionModule                                                  *
!                                                                               *
! Version 1.0: June 2011 by wkb                                                 *
!                                                                               *
! Description                                                                   *
!    Calculates the relativistic Green's function and integrates it along       *
!    a radial mesh (all in kappa,mu form).                                      *
!                                                                               *
! ----------------------------------------------------------------------------- *
! subroutine initRelGreenFunction(na, atomic_index, lmaxkkr, lmaxphi, &         *
!                                 lmaxgreen, posi, pola, cant,          &       *
!                                 i_stop, i_print)                              *
!    Initizalizes variables and allocates the memory used by the routines       *
!    in this module. Must be called first.                                      *
!  Input: na - Integer, number of local atoms.                                  *
!         atomic_index - 1D integer array of size na, containing the global     *
!                        index of each local atom.                              *
!         lmaxkkr - 1D integer array of size na, containing the single-site     *
!                   l-index cutoff for each local atom.                         *
!         lmaxphi - 1D integer array of size na, containing the single-site     *
!                   wavefunction expansion cutoff for each local atom.          *
!         lmaxgreen - 1D integer array of size na, containing the single-site   *
!                   green function expansion cutoff for each local atom.        *
!         posi - 2D real array of size (3,na) containing the position           *
!                coordinates of each local atom.                                *
!         pola - Integer, spin polarization index. Must equal 2 for             *
!                relativistic calculations.                                     *
!         cant - Integer, spin canting index. Must equal 2 for relativistic     *
!                calculations.                                                  *
!         i_stop - Character string, name of the routine in which to end        *
!                  processing.                                                  *
!         i_print - Integer, level of printed output                            *
! ----------------------------------------------------------------------------- *
! subroutine endRelGreenFunction()                                              *
!    Deallocates the memory used by routines in this module. This must be the   *
!    last routine called, unless initRelGreenFunction is used again.            *
! ----------------------------------------------------------------------------- *
! function getRelGreenFunction(atom, e) result(green)                           *
!    Calculates and returns the relativistic Green's function for a given       *
!    atom.                                                                      *
!  Input: atom - Integer, local atom index                                      *
!         e - Complex energy                                                    *
! Output: green - Pointer to a 3D complex array of size (jws,kmax,4) containing *
!                 the four-component relativistic Green's function, where jws   *
!                 is the number of radial points in the W-S sphere. Currently,  *
!                 kmax is assumed to be 1 (lmax_green=0).                       *
! ----------------------------------------------------------------------------- *
! function getRelGreenIntegral(ia, e) result(integral)                          *
!    Calculates and returns the integral of the Green's function along the      *
!    radial mesh from the center to the W-S sphere.                             *
!  Input: ia - Integer, local atom index                                        *
!         e - Complex energy                                                    *
!  Output: integral - Pointer to a 1D complex array of size 4, containing the   *
!                     integral of the 4-component Green's function from the     *
!                     center of the W-S sphere to its edge.                     *
! ----------------------------------------------------------------------------- *
! function getRelGreenIntegralMT(ia, e) result(integral)                        *
!    Calculates and returns the integra; of the Green's function along the      *
!    radial mesh from the center to the MT sphere.                              *
!  Input: Same as above.                                                        *
!  Output: integral - Pointer to a 1D complex array of size 4, containing the   *
!                     integral of the 4-component Green's function from the     *
!                     center of the MT sphere to its edge.                      *
!-------------------------------------------------------------------------------*
! The following functions return a vector or scalar containing zero. They can   *
! implemented at a later date if needed.                                        *
!        getSpinTorqueMomentR                                                   *
!        getStonerParamMomentR                                                  *
!        getOneSiteStonerParamMomentR                                           *
!        getOneSiteSusceptibilityMomentR                                        *
!        getOneSiteExchangeMomentR                                              *
!        getPairExchangeMomentR                                                 *
!********************************************************************************

use KindParamModule, only : IntKind, RealKind, CmplxKind
use MathParamModule, only : ZERO, HALF, ONE, TWO, FOUR, PI, PI4, SQRT_PI, &
                            SQRTM1, TEN2m10, CZERO, CONE, Y0
use ErrorHandlerModule, only : ErrorHandler, StopHandler

public :: initRelGreenFunction,            &
          endRelGreenFunction,             &
          getRelGreenFunction,             &
          getRelGreenIntegral,             &
          getRelGreenIntegralMT,           &
          getSpinTorqueMomentR,            &
          getStonerParamMomentR,           &
          getOneSiteStonerParamMomentR,    &
          getOneSiteSusceptibilityMomentR, &
          getOneSiteExchangeParamMomentR,  &
          getPairExchangeParamMomentR

! Relativistic Gaunt coefficients, for gfill and gafill
!-----------------------------------------------------------------
integer (kind=IntKind) :: kmymaxp

complex (kind=CmplxKind), allocatable :: sxcoeff(:,:), sxbcoeff(:,:)
complex (kind=CmplxKind), allocatable :: sycoeff(:,:), sybcoeff(:,:)
complex (kind=CmplxKind), allocatable :: szcoeff(:,:), szbcoeff(:,:)

complex (kind=CmplxKind), allocatable :: lxcoeff(:,:), lxbcoeff(:,:)
complex (kind=CmplxKind), allocatable :: lycoeff(:,:), lybcoeff(:,:)
complex (kind=CmplxKind), allocatable :: lzcoeff(:,:), lzbcoeff(:,:)

complex (kind=CmplxKind), allocatable :: gacoeff(:,:,:)
complex (kind=CmplxKind), allocatable :: rgacoeff(:,:,:)
!-------------------------------------------------------------------

logical :: Initialized = .false.
logical, allocatable :: done(:)

character (len=32) :: istop

real (kind=RealKind), allocatable :: xstep(:)

integer (kind=IntKind), allocatable :: j_mt(:), j_ws(:)
integer (kind=IntKind), allocatable :: l_max(:), kkr_sz(:)
integer (kind=IntKind) :: jmt, jws     ! The current values, used for
integer (kind=IntKind) :: lmax, kkrsz  ! easier access
integer (kind=IntKind) :: lmax_max
integer (kind=IntKind) :: num_rpts
integer (kind=IntKind) :: LocalNumAtoms
integer (kind=IntKind), pointer :: iprint(:)
integer (kind=IntKind), parameter :: nuzp = 2

real (kind=RealKind), parameter :: small=1.d-15

complex (kind=CmplxKind) :: Energy
complex (kind=CmplxKind), allocatable, target :: greenfunc(:,:,:,:)
complex (kind=CmplxKind), allocatable, target :: greenint(:,:)
complex (kind=CmplxKind), allocatable, target :: greenintmt(:,:)

contains

!=====================================================================
   subroutine initRelGreenFunction(na, atomic_index, lmaxkkr,    &
                                   lmaxphi, lmaxgreen, posi, pola, &
                                   cant, i_stop, i_print)
!=====================================================================

   use PublicTypeDefinitionsModule, only : GridStruct
   use RadialGridModule, only : getGrid, getMaxNumRmesh
   use ScfDataModule, only : isKKR, isLSMS
   use SystemModule, only : getBravaisLattice
!v2.0   use MultipleScatteringModule, only : initMultipleScattering
   use MSSolverModule, only : initMSSolver

   implicit none

   character (len=32) :: i_stop

   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: NumRmesh
   integer (kind=IntKind), intent(in) :: na, lmaxkkr(na)
   integer (kind=IntKind), intent(in) :: lmaxphi(na), lmaxgreen(na)
   integer (kind=IntKind), intent(in) :: pola, cant
   integer (kind=IntKind), intent(in) :: atomic_index(na)
   integer (kind=IntKind), intent(in), target :: i_print(na)

   real (kind=RealKind), intent(in) :: posi(3,na)
   real (kind=RealKind) :: bravais(3,3)

   type (GridStruct), pointer :: Grid

   integer (kind=IntKind), parameter :: lamp=1, lammp=2*(lamp+1)*(lamp+1)

   if(pola.ne.2 .or. cant.ne.2) then
      call ErrorHandler('initRelGreenFunction', &
                        'pola and cant must equal 2 for relativistic calculations')
   endif

   allocate( done(na) )
   allocate( j_mt(na) )
   allocate( j_ws(na) )
   allocate( l_max(na) )
   allocate( kkr_sz(na) )
   allocate( greenint(4, na) )
   allocate( greenintmt(4, na) )
   allocate( xstep(na) )

   done(1:na) = .false.
   LocalNumAtoms = na
   NumRmesh = getMaxNumRmesh()
   lmax_max = 1

   do i = 1,na
      Grid => getGrid(i)
      j_mt(i) = Grid%jmt
      j_ws(i) = Grid%jend
      xstep(i) = Grid%hin
      l_max(i) = lmaxkkr(i)
      kkr_sz(i) = (lmaxkkr(i)+1)*(lmaxkkr(i)+1)
      lmax_max = max(lmaxkkr(i), lmax_max)
   enddo

   allocate( greenfunc(NumRmesh,1,4,na) )

   kmymaxp = 2*(lmax_max+1)*(lmax_max+1)

   allocate( sxcoeff(kmymaxp,kmymaxp) )
   allocate( sxbcoeff(kmymaxp,kmymaxp) )
   allocate( sycoeff(kmymaxp,kmymaxp) )
   allocate( sybcoeff(kmymaxp,kmymaxp) )
   allocate( szcoeff(kmymaxp,kmymaxp) )
   allocate( szbcoeff(kmymaxp,kmymaxp) )
   allocate( lxcoeff(kmymaxp,kmymaxp) )
   allocate( lxbcoeff(kmymaxp,kmymaxp) )
   allocate( lycoeff(kmymaxp,kmymaxp) )
   allocate( lybcoeff(kmymaxp,kmymaxp) )
   allocate( lzcoeff(kmymaxp,kmymaxp) )
   allocate( lzbcoeff(kmymaxp,kmymaxp) )
   allocate( gacoeff((lmax_max+1)*(lmax_max+1), &
                     (lmax_max+1)*(lmax_max+1),lammp) )
   allocate( rgacoeff(kmymaxp,kmymaxp,lammp) )

   call gfill
   call gafill

   iprint => i_print
   istop = i_stop

   nullify( Grid )

   if (isLSMS()) then
!v2.0      call initMultipleScattering(na, atomic_index, lmaxkkr, &
!v2.0                                  lmaxphi, lmaxpot,  posi,   &
!v2.0                                  pola, cant, 2, i_stop, i_print)
      call initMSSolver(na, atomic_index, lmaxkkr, &
                        lmaxphi, lmaxgreen,  posi,   &
                        pola, cant, 2, i_stop, i_print)
   else if (isKKR()) then
      bravais(1:3,1:3) = getBravaisLattice()
!v2.0      call initMultipleScattering(bravais, na, lmaxkkr, lmaxphi, &
!v2.0                                  lmaxpot, posi, pola, cant,     &
!v2.0                                  2, i_stop, i_print)
      call initMSSolver(na, atomic_index, lmaxkkr, &
                        lmaxphi, lmaxgreen,  posi,   &
                        pola, cant, 2, i_stop, i_print)
   else
      call ErrorHandler('initRelGreenFunction', &
            'Relativistic calculations only support LSMS and KKR methods')
   endif

   Energy = CZERO
   Initialized = .true.

   end subroutine initRelGreenFunction
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endRelGreenFunction()
!  ===================================================================

!v2.0   use MultipleScatteringModule, only : endMultipleScattering
   use MSSolverModule, only : endMSSolver

   implicit none

   if (.not.Initialized) then
      call ErrorHandler('endRelGreen', 'module has not been initialized')
   endif

   Initialized = .false.

   deallocate( sxcoeff, sxbcoeff )
   deallocate( sycoeff, sybcoeff )
   deallocate( szcoeff, szbcoeff )
   deallocate( lxcoeff, lxbcoeff )
   deallocate( lycoeff, lybcoeff )
   deallocate( lzcoeff, lzbcoeff )
   deallocate( gacoeff, rgacoeff )

   deallocate( j_mt, j_ws )
   deallocate( l_max, kkr_sz )
   deallocate( greenint, greenintmt )
   deallocate( greenfunc )
   deallocate( xstep )
   deallocate( done )

!v2.0   call endMultipleScattering()
   call endMSSolver()
!
   Energy = CZERO
   Initialized = .false.

   end subroutine endRelGreenFunction
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRelGreenFunction(e)
!  ===================================================================
!v2.0   use MultipleScatteringModule, only : solveMSTeqns
   use MSSolverModule, only : computeMSGreenFunction
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e

   character (len=20), parameter :: sname = 'calRelGreenFunction'
!

   if ( abs(e-Energy) > TEN2m10 ) then
!     ----------------------------------------------------------------
!v2.0      call solveMSTeqns(e, 1)
      call computeMSGreenFunction(1, e)
!     ----------------------------------------------------------------
   endif
!
   Energy = e
!
   end subroutine calRelGreenFunction
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRelGreenFunction(atom, e) result(green)
!  ===================================================================

! Original arguments included dos, dosck, dos_orb, dosck_orb, dens_orb

! original output:
!                dos     (wigner-seitz cell density of states)
!                dosck   (muffin-tin density of states)
!                green   (Green's function) (actually complex
!                              charge & spin - magnetization densities)
!                dos_orb    )
!                dosck_orb  ) orbital dos & magn. densities
!                dens_orb   )

   use MathParamModule, only : TEN2m10
   use AtomModule, only : getLocalEvecOld
   use DiracSolverModule, only : getSizeOfRMesh, getDmat, getDmatP
   use RelScattererModule, only : getRelTmat, getNumRadialPts, &
                                  getGZ, getFZ, getGJ, getFJ,  &
                                  get_Nuz, get_Indz
   use RadialGridModule, only : getRmesh, getNumRmesh, getMaxNumRmesh
   use RelativityToolsModule, only : repl, doubmt, submat, rsimp, &
                                     rotr, tripmt
!   use SpinRotationModule, only : TransformDensityMatrix
   use ScfDataModule, only : isKKR, isLSMS
!v2.0   use MultipleScatteringModule, only : getTauijRInGlobalSpinFrame, &
!v2.0                                        getTauijKInGlobalSpinFrame

   implicit none

   character (len=32), parameter :: sname = 'getRelGreenFunction'

   integer (kind=IntKind), intent(in) :: atom
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: ir, irmax
   integer (kind=IntKind) :: NumRmesh
   integer (kind=IntKind) :: i, j
   integer (kind=IntKind) :: kmy, kmyp
   integer (kind=IntKind) :: nu, nup
   integer (kind=IntKind) :: kmy1, kmyp1
   integer (kind=IntKind), pointer :: nuz(:)
   integer (kind=IntKind), pointer :: indz(:,:)
   integer (kind=IntKind) :: lambda
   integer (kind=IntKind) :: kmymax
!
!   real (kind=RealKind), allocatable :: rrr(:), rir(:)
!   real (kind=RealKind), allocatable :: rri(:), rii(:)
!   real (kind=RealKind) :: r1rr, r1ir
!   real (kind=RealKind) :: r1ri, r1ii
   real (kind=RealKind) :: rlam
   real (kind=RealKind) :: fac
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind) :: h, phi, a
   real (kind=RealKind) :: axis(3), rot(3,3), evec_r(3)

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: tau00(:,:)
   complex (kind=CmplxKind), pointer :: gz(:,:,:)
   complex (kind=CmplxKind), pointer :: fz(:,:,:)
   complex (kind=CmplxKind), pointer :: gj(:,:,:)
   complex (kind=CmplxKind), pointer :: fj(:,:,:)
   complex (kind=CmplxKind), pointer :: dmat(:,:), dmatp(:,:)
   complex (kind=CmplxKind), allocatable, target :: tau00L(:,:)
!   complex (kind=CmplxKind), allocatable :: sum(:,:,:)
!   complex (kind=CmplxKind), allocatable :: gf(:,:)
   complex (kind=CmplxKind) :: cxr,cxi
   complex (kind=CmplxKind) :: cox
   complex (kind=CmplxKind) :: t1, t2, t3
!   complex (kind=CmplxKind) :: dos(4) ! Output, these are already calculated in ValenceStatesModule
!   complex (kind=CmplxKind) :: dosck(4) ! Output
   complex (kind=CmplxKind), pointer :: green(:,:,:)
!   complex (kind=CmplxKind) :: dos_orb(3) ! Output
!   complex (kind=CmplxKind) :: dosck_orb(3) ! Output
!  complex (kind=CmplxKind), allocatable :: dens_orb(:,:) ! Output

   if (.not.Initialized) then
      call ErrorHandler('getRelGreenFunction', 'Module has not been initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelGreenFunction', 'Invalid local atom index', atom)
   else if ( abs(e-Energy) > TEN2m10 ) then
      call ErrorHandler('getRelGreenFunction','Need to call to call calRelGreenFunction first')
   endif

   jmt = j_mt(atom)
   jws = j_ws(atom)
   num_rpts = getNumRmesh(atom)
   irmax = getNumRadialPts(atom)
   lmax = l_max(atom)
   kkrsz = kkr_sz(atom)
   h = xstep(atom)
   lambda = 0
   fac = -ONE/pi
   kmymax=2*(lmax+1)*(lmax+1)

   allocate( tau00L(kmymax,kmymax) )

   if ( isLSMS() ) then
!v2.0      tau00 => getTauijRInGlobalSpinFrame(atom,1,0,0)
stop 'Under construction...'
   else if ( isKKR() ) then
!v2.0      tau00 => getTauijKInGlobalSpinFrame(atom,1)
stop 'Under construction...'
   else
      call ErrorHandler('getRelGreenFunction', 'Must use LSMS or KKR')
   endif

!  Transform tau00 into the local frame:
   dmat => getDmat(atom)
   dmatp => getDmatP(atom)
   call repl(tau00L, tau00, kmymax, kmymax)
   call tripmt(dmatp,tau00L,dmat,kmymax,kmymax,kmymax)
   tau00 => tau00L

   NumRmesh = getMaxNumRmesh()
   green => greenfunc(1:NumRmesh, 1:1, 1:4, atom)
   r_mesh => getRmesh(atom)
   gz => getGZ(atom)
   fz => getFZ(atom)
   gj => getGJ(atom)
   fj => getFJ(atom)
   nuz => get_Nuz(atom)
   indz => get_Indz(atom)

!   allocate( rrr(num_rpts) )
!   allocate( rri(num_rpts) )
!   allocate( rir(num_rpts) )
!   allocate( rii(num_rpts) )
!   allocate( sum(kkrsz*2,kkrsz*2,2) )
!   allocate( gf(kkrsz*2,kkrsz*2) )
!
!   do kmy=1,kmymax
!      do kmyp=1,kmymax
!         call zeroout(rrr,jmt)
!         call zeroout(rir,jmt)
!         call zeroout(rri,jmt)
!         call zeroout(rii,jmt)
!         do nu=1,nuz(kmy)
!            kmy1=indz(nu,kmy)
!            do nup=1,nuz(kmyp)
!               kmyp1=indz(nup,kmyp)
!               cox=fac*rgacoeff(kmy1,kmyp1,1)
!               if(cdabs(cox).gt.small) then
!                  do i=1,jmt
!!                    rlam=r_mesh(i)**lambda
!                     rlam=ONE
!                     cxr=(gz(i,nu,kmy)*gz(i,nup,kmyp)+ &
!                         fz(i,nu,kmy)*fz(i,nup,kmyp))*cox*rlam
!                     cxi=(gz(i,nu,kmy)*gj(i,nup,kmyp)+ &
!                         fz(i,nu,kmy)*fj(i,nup,kmyp))*cox*rlam
!                     rrr(i)=rrr(i)+dreal(cxr)
!                     rir(i)=rir(i)+dimag(cxr)
!                     rri(i)=rri(i)+dreal(cxi)
!                     rii(i)=rii(i)+dimag(cxi)
!                  end do
!               end if
!            end do ! nup
!         end do   ! nu
!         r1rr = rsimp(rrr,r_mesh,jmt,h)
!         r1ir = rsimp(rir,r_mesh,jmt,h)
!         r1ri = rsimp(rri,r_mesh,jmt,h)
!         r1ii = rsimp(rii,r_mesh,jmt,h)
!         sum(kmy,kmyp,1) = dcmplx(r1rr,r1ir)
!         sum(kmy,kmyp,2) = dcmplx(r1ri,r1ii)
!      end do ! kmyp
!   end do   ! kmy

!   call repl(gf,sum(:,:,1),kmymax,kmymax)
!   call doubmt(gf,tau00,kmymax,kmymax)
!   call submat(gf,sum(:,:,2),kmymax,kmymax)

!   dos(1) = ZERO
!   do i=1,kmymax
!      dos(1) = dos(1)+gf(i,i)
!   end do

!   call magnet(tau00, r_mesh, h,     &
!               gz,fz,gj,fj,nuz,indz, &
!               sxcoeff,sxbcoeff,     &
!               dos(2))
!
!   call magnet(tau00, r_mesh, h,     &
!               gz,fz,gj,fj,nuz,indz, &
!               sycoeff,sybcoeff,     &
!               dos(3))
!
!   call magnet(tau00, r_mesh, h,     &
!               gz,fz,gj,fj,nuz,indz, &
!               szcoeff,szbcoeff,     &
!               dos(4))
!
!   call magnet(tau00, r_mesh, h,     &
!               gz,fz,gj,fj,nuz,indz, &
!               lxcoeff,lxbcoeff,     &
!               dos_orb(1))
!
!
!   call magnet(tau00, r_mesh, h,     &
!               gz,fz,gj,fj,nuz,indz, &
!               lycoeff,lybcoeff,     &
!               dos_orb(2))
!
!   call magnet(tau00, r_mesh, h,     &
!               gz,fz,gj,fj,nuz,indz, &
!               lzcoeff,lzbcoeff,     &
!               dos_orb(3))

!!!
!!! The following changes have been made by Yang Wang, following William Balunas' work
!!!
   call new_dens(irmax,gz,fz,gj,fj,nuz,indz, &       ! changes: jmt -> irmax
                 tau00, green(1:irmax,1,1),r_mesh)   ! changes: jmt -> irmax

   call magnetic_dens(irmax,gz,fz,gj,fj,   &    ! add irmax as the 1st argument
                      nuz,indz,tau00,      &
                      sxcoeff,sxbcoeff,    &
                      green(1:irmax,1,2))       ! changes: jmt -> irmax
 
   call magnetic_dens(irmax,gz,fz,gj,fj,   &    ! add irmax as the 1st argument
                      nuz,indz,tau00,      &
                      sycoeff,sybcoeff,    &
                      green(1:irmax,1,3))       ! changes: jmt -> irmax

   call magnetic_dens(irmax,gz,fz,gj,fj,   &
                      nuz,indz,tau00,      &
                      szcoeff,szbcoeff,    &
                      green(1:irmax,1,4))       ! changes: jmt -> irmax

!     call magnetic_dens(lmax,kkrsz,
!    >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
!    >                lxcoeff,lxbcoeff,
!    >                dens_orb(1,1),
!    >                iprint,istop)
!
!     call magnetic_dens(lmax,kkrsz,
!    >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
!    >                lycoeff,lybcoeff,
!    >                dens_orb(1,2),
!    >                iprint,istop)
!
!     call magnetic_dens(lmax,kkrsz,
!    >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
!    >                lzcoeff,lzbcoeff,
!    >                dens_orb(1,3),
!    >                iprint,istop)


! currently we are have only implemented ASA for the relativistic case,
! therefore dos == dosck

!   dosck(1) = dos(1)
!   dosck(2) = dos(2)  
!   dosck(3) = dos(3)
!   dosck(4) = dos(4)
!   dosck_orb(1) = dos_orb(1)
!   dosck_orb(2) = dos_orb(2)
!   dosck_orb(3) = dos_orb(3)

!   if(iprint(atom).ge.2) then
!      write(6,*) 'dos(1)=',dos(1)
!      write(6,*) 'dos(2)=',dos(2)
!      write(6,*) 'dos(3)=',dos(3)
!      write(6,*) 'dos(4)=',dos(4)
!      write(6,*) 'dos_orb(1)=',dos_orb(1)
!      write(6,*) 'dos_orb(2)=',dos_orb(2)
!      write(6,*) 'dos_orb(3)=',dos_orb(3)
!   end if

   done(atom) = .true.

   nullify( tau00 )
!   deallocate( rrr, rri )
!   deallocate( rir, rii )
!   deallocate( sum )
!   deallocate( gf )
   deallocate( tau00L )

   if(istop.eq.sname) then
      call StopHandler(sname)
   end if

! we now have to rotate into the global frame of reference:
!   evec_r(1:3) = getLocalEvecOld(atom)
!   a = sqrt(evec_r(1)**2+evec_r(2)**2)
!   phi = acos(evec_r(3))
!   if(a.eq.0.d0) then
!      axis(1) =  0.d0
!      axis(2) =  0.d0
!      axis(3) =  1.d0
!      axis(3) = sign(1.d0,evec_r(3))
!      a = 1.d0
!   else
!      axis(1) = -evec_r(2)
!      axis(2) =  evec_r(1)
!      axis(3) =  0.d0
!      a = 1.d0/sqrt(axis(1)**2+axis(2)**2+axis(3)**2)
!      a = 1.d0/a
!   end if
!   axis(1) = axis(1)*a
!   axis(2) = axis(2)*a
!   axis(3) = axis(3)*a
!   call rotr(rot,axis,phi)
!   t1 = rot(1,1)*dos(2)+rot(1,2)*dos(3)+rot(1,3)*dos(4)
!   t2 = rot(2,1)*dos(2)+rot(2,2)*dos(3)+rot(2,3)*dos(4)
!   t3 = rot(3,1)*dos(2)+rot(3,2)*dos(3)+rot(3,3)*dos(4)
!   dos(2) = t1
!   dos(3) = t2
!   dos(4) = t3
!   t1 = rot(1,1)*dosck(2)+rot(1,2)*dosck(3)+rot(1,3)*dosck(4)
!   t2 = rot(2,1)*dosck(2)+rot(2,2)*dosck(3)+rot(2,3)*dosck(4)
!   t3 = rot(3,1)*dosck(2)+rot(3,2)*dosck(3)+rot(3,3)*dosck(4)
!   dosck(2) = t1
!   dosck(3) = t2
!   dosck(4) = t3
!   t1 = rot(1,1)*dos_orb(1)+rot(1,2)*dos_orb(2)+ &
!        rot(1,3)*dos_orb(3)
!   t2 = rot(2,1)*dos_orb(1)+rot(2,2)*dos_orb(2)+ &
!        rot(2,3)*dos_orb(3)
!   t3 = rot(3,1)*dos_orb(1)+rot(3,2)*dos_orb(2)+ &
!        rot(3,3)*dos_orb(3)
!   dos_orb(1) = t1
!   dos_orb(2) = t2
!   dos_orb(3) = t3
!   t1 = rot(1,1)*dosck_orb(1)+rot(1,2)*dosck_orb(2)+ &
!        rot(1,3)*dosck_orb(3)
!   t2 = rot(2,1)*dosck_orb(1)+rot(2,2)*dosck_orb(2)+ &
!        rot(2,3)*dosck_orb(3)
!   t3 = rot(3,1)*dosck_orb(1)+rot(3,2)*dosck_orb(2)+ &
!        rot(3,3)*dosck_orb(3)
!   dosck_orb(1) = t1
!   dosck_orb(2) = t2
!   dosck_orb(3) = t3
!   do j=1,jws
!      t1 = rot(1,1)*green(j,1,2)+rot(1,2)*green(j,1,3)+ &
!           rot(1,3)*green(j,1,4)
!      t2 = rot(2,1)*green(j,1,2)+rot(2,2)*green(j,1,3)+ &
!           rot(2,3)*green(j,1,4)
!      t3 = rot(3,1)*green(j,1,2)+rot(3,2)*green(j,1,3)+ &
!           rot(3,3)*green(j,1,4)
!      green(j,1,2) = t1
!      green(j,1,3) = t2
!      green(j,1,4) = t3
!!      t1 = rot(1,1)*dens_orb(j,1)+rot(1,2)*dens_orb(j,2)+ &
!!           rot(1,3)*dens_orb(j,3)
!!      t2 = rot(2,1)*dens_orb(j,1)+rot(2,2)*dens_orb(j,2)+ &
!!           rot(2,3)*dens_orb(j,3)
!!      t3 = rot(3,1)*dens_orb(j,1)+rot(3,2)*dens_orb(j,2)+ &
!!           rot(3,3)*dens_orb(j,3)
!!      dens_orb(j,1) = t1
!!      dens_orb(j,2) = t2
!!      dens_orb(j,3) = t3
!   enddo

   green(1:irmax,1,1:4) = green(1:irmax,1,1:4)/Y0/PI4

   end function getRelGreenFunction

!======================================================================
   function getRelGreenIntegral(ia, e) result(integral)
!======================================================================

   use RadialGridModule, only : getRmesh
   use RelScattererModule, only : getNumRadialPts
   use IntegrationModule, only : calIntegration

   implicit none

   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind) :: iend
   real (kind=RealKind), pointer :: rmesh(:)
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: temp(1:jws)
   complex (kind=CmplxKind), pointer :: green(:)
   complex (kind=CmplxKind), pointer :: integral(:)

   if (.not.Initialized) then
      call ErrorHandler('getRelGreenIntegral', 'Module has not been initialized')
   else if (ia < 1 .or. ia > LocalNumAtoms) then
      call ErrorHandler('getRelGreenIntegral', 'Invalid atom index', ia)
   else if (.not.done(ia)) then
      call ErrorHandler('getRelGreenIntegral', 'Must call RelGreenFunction for this atom first')
   else if (abs(e-Energy) > TEN2m10) then
      call ErrorHandler('getRelGreenIntegral','Need to call to call calRelGreenFunction first')
   endif

   rmesh => getRmesh(ia)
   iend = getNumRadialPts(ia)
   green => greenfunc(1:iend, 1, 1, ia)
   integral => greenint(1:4, ia)

   call calIntegration(0, jws, rmesh, green, temp)
   integral(1) = TWO*SQRT_PI*temp(jws)
   green => greenfunc(1:iend, 1, 2, ia)
   call calIntegration(0, jws, rmesh, green, temp)
   integral(2) = TWO*SQRT_PI*temp(jws)
   green => greenfunc(1:iend, 1, 3, ia)
   call calIntegration(0, jws, rmesh, green, temp)
   integral(3) = TWO*SQRT_PI*temp(jws)
   green => greenfunc(1:iend, 1, 4, ia)
   call calIntegration(0, jws, rmesh, green, temp)
   integral(4) = TWO*SQRT_PI*temp(jws)

   end function getRelGreenIntegral

!======================================================================
   function getRelGreenIntegralMT(ia, e) result(integral)
!======================================================================

   use RadialGridModule, only : getRmesh
   use RelScattererModule, only : getNumRadialPts
   use IntegrationModule, only : calIntegration

   implicit none

   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind) :: iend
   real (kind=RealKind), pointer :: rmesh(:)
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: temp(1:jmt)
   complex (kind=CmplxKind), pointer :: green(:)
   complex (kind=CmplxKind), pointer :: integral(:)

   if (.not.Initialized) then
      call ErrorHandler('getRelGreenIntegralMT', 'Module has not been initialized')
   else if (ia < 1 .or. ia > LocalNumAtoms) then
      call ErrorHandler('getRelGreenIntegralMT', 'Invalid atom index', ia)
   else if (.not.done(ia)) then
      call ErrorHandler('getRelGreenIntegralMT', 'Must call RelGreenFunction for this atom first')
   else if (abs(e-Energy) > TEN2m10) then
      call ErrorHandler('getRelGreenIntegralMT','Need to call calRelGreenFunction first')
   endif

   rmesh => getRmesh(ia)
   iend = getNumRadialPts(ia)
   green => greenfunc(1:iend, 1, 1, ia)
   integral => greenintmt(1:4, ia)

   call calIntegration(0, jmt, rmesh, green, temp)
   integral(1) = TWO*SQRT_PI*temp(jmt)
   green => greenfunc(1:iend, 1, 2, ia)
   call calIntegration(0, jmt, rmesh, green, temp)
   integral(2) = TWO*SQRT_PI*temp(jmt)
   green => greenfunc(1:iend, 1, 3, ia)
   call calIntegration(0, jmt, rmesh, green, temp)
   integral(3) = TWO*SQRT_PI*temp(jmt)
   green => greenfunc(1:iend, 1, 4, ia)
   call calIntegration(0, jmt, rmesh, green, temp)
   integral(4) = TWO*SQRT_PI*temp(jmt)

   end function getRelGreenIntegralMT

!======================================================================
   function getSpinTorqueMomentR(id, de) result(td)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: id
   real (kind=RealKind) :: td(3)
   complex (kind=CmplxKind), intent(in) :: de

   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSpinTorqueMomentR', 'Invalid local atom index', id)
   endif

   td(1) = ZERO
   td(2) = ZERO
   td(3) = ZERO

   end function getSpinTorqueMomentR

!======================================================================
   function getStonerParamMomentR(id, de) result(st)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: id
   real (kind=RealKind) :: st(3)
   complex (kind=CmplxKind), intent(in) :: de

   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getStonerParamMoment','Invalid local atom index', id)
   endif

   st(1) = ZERO
   st(2) = ZERO
   st(3) = ZERO

   end function getStonerParamMomentR

!======================================================================
   function getOneSiteStonerParamMomentR(id,de) result(os)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: id
   real (kind=RealKind) :: os
   complex (kind=CmplxKind), intent(in) :: de

   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getOneSiteStonerParamMoment', 'Invalid local atom index', id)
   endif

   os = ZERO

   end function getOneSiteStonerParamMomentR

!=====================================================================
   function getOneSiteSusceptibilityMomentR(id,de) result(os)
!=====================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i
   real (kind=RealKind) :: os(0:9)
   complex (kind=CmplxKind), intent(in) :: de

   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getOneSiteSusceptibilityMoment', 'Invalid local atom index',id)
   endif

   os(0) = ZERO
   do i=1,9
      os(i) = ZERO
   enddo

   end function getOneSiteSusceptibilityMomentR

!======================================================================
   function getOneSiteExchangeParamMomentR(id,de) result(osex)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: id
   real (kind=RealKind) :: osex
   complex (kind=CmplxKind), intent(in) :: de

   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getOneSiteExchangeParamMoment', 'Invalid local atom index', id)
   endif

   osex = ZERO

   end function getOneSiteExchangeParamMomentR

!=======================================================================
   function getPairExchangeParamMomentR(id,j,de,jid,r0j) result(exch)
!=======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: id, j
   integer (kind=IntKind), intent(out) :: jid
   integer (kind=IntKind) :: i

   real (kind=RealKind) :: exch(9)
   real (kind=RealKind), intent(out) :: r0j(3)

   complex (kind=CmplxKind), intent(in) :: de

   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPairExchangeParamMoment', 'Invalid local atom index', id)
   endif

   jid = 0
   r0j(1:3) = ZERO
   do i=1,9
      exch(i) = ZERO
   enddo

   end function getPairExchangeParamMomentR

!!======================================================================
!   subroutine magnet(tau00, r_mesh, h,      &
!                     gz,fz,gj,fj,nuz,indz,  &
!                     gcoeff,gbcoeff,        &
!                     out)
!!======================================================================
!
!   use RelativityToolsModule, only : repl, doubmt, submat, rsimp
!
!   implicit none
!
!   character (len=32), parameter :: sname = 'magnet'
!
!   integer (kind=IntKind) :: kmymax
!   integer (kind=IntKind) :: m
!   integer (kind=IntKind) :: ir
!   integer (kind=IntKind) :: i
!   integer (kind=IntKind) :: kmy, kmyp
!   integer (kind=IntKind) :: nu, nup
!   integer (kind=IntKind) :: kmy1, kmyp1
!
!   integer (kind=IntKind), intent(in) :: nuz(:)
!   integer (kind=IntKind), intent(in) :: indz(:,:)
!
!   real (kind=RealKind) :: fac
!   real (kind=RealKind), intent(in) :: r_mesh(:)
!   real (kind=RealKind), intent(in) :: h
!   real (kind=RealKind) :: rrr(num_rpts), rir(num_rpts)
!   real (kind=RealKind) :: rri(num_rpts), rii(num_rpts)
!   real (kind=RealKind) :: r1rr, r1ri
!   real (kind=RealKind) :: r1ir, r1ii
!
!   complex (kind=CmplxKind), intent(in) :: tau00(:,:)
!   complex (kind=CmplxKind), intent(in) :: gz(:,:,:)
!   complex (kind=CmplxKind), intent(in) :: fz(:,:,:)
!   complex (kind=CmplxKind), intent(in) :: gj(:,:,:)
!   complex (kind=CmplxKind), intent(in) :: fj(:,:,:)
!   complex (kind=CmplxKind) :: sum(kkrsz*2,kkrsz*2,2)
!   complex (kind=CmplxKind) :: gf(kkrsz*2,kkrsz*2)
!   complex (kind=CmplxKind), intent(out) :: out
!   complex (kind=CmplxKind) :: cxr,cxi
!   complex (kind=CmplxKind), intent(in) :: gcoeff(:,:)
!   complex (kind=CmplxKind), intent(in) :: gbcoeff(:,:)
!   complex (kind=CmplxKind) :: gcox, gbcox
!
!   fac = -ONE/pi
!   kmymax=2*(lmax+1)*(lmax+1)
!
!   do kmy=1,kmymax
!      do kmyp=1,kmymax
!         call zeroout(rrr,jmt)
!         call zeroout(rir,jmt)
!         call zeroout(rri,jmt)
!         call zeroout(rii,jmt)
!         do nu=1,nuz(kmy)
!            kmy1=indz(nu,kmy)
!            do nup=1,nuz(kmyp)
!               kmyp1=indz(nup,kmyp)
!               gcox=fac*gcoeff(kmy1,kmyp1)
!               gbcox=fac*gbcoeff(kmy1,kmyp1)
!               if(abs(gcox).gt.small) then
!                  do i=1,jmt
!                     cxr=gz(i,nu,kmy)*gz(i,nup,kmyp)*gcox
!                     cxi=gz(i,nu,kmy)*gj(i,nup,kmyp)*gcox
!                     rrr(i)=rrr(i)+dreal(cxr)
!                     rir(i)=rir(i)+dimag(cxr)
!                     rri(i)=rri(i)+dreal(cxi)
!                     rii(i)=rii(i)+dimag(cxi)
!                  end do
!               end if
!               if(abs(gbcox).gt.small) then
!                  do i=1,jmt
!                     cxr=fz(i,nu,kmy)*fz(i,nup,kmyp)*gbcox
!                     cxi=fz(i,nu,kmy)*fj(i,nup,kmyp)*gbcox
!                     rrr(i)=rrr(i)-dreal(cxr)
!                     rir(i)=rir(i)-dimag(cxr)
!                     rri(i)=rri(i)-dreal(cxi)
!                     rii(i)=rii(i)-dimag(cxi)
!                  end do
!               end if
!            end do ! nup
!         end do   ! nu
!         r1rr = rsimp(rrr,r_mesh,jmt,h)
!         r1ir = rsimp(rir,r_mesh,jmt,h)
!         r1ri = rsimp(rri,r_mesh,jmt,h)
!         r1ii = rsimp(rii,r_mesh,jmt,h)
!         sum(kmy,kmyp,1) = dcmplx(r1rr,r1ir)
!         sum(kmy,kmyp,2) = dcmplx(r1ri,r1ii)
!      end do ! kmyp
!   end do   ! kmy
!
!   call repl(gf,sum(:,:,1),kmymax,kmymax)
!   call doubmt(gf,tau00,kmymax,kmymax)
!   call submat(gf,sum(:,:,2),kmymax,kmymax)
!
!   out = CZERO
!   do i=1,kmymax
!      out = out + gf(i,i)
!   end do
!
!   if(istop.eq.sname) then
!      call StopHandler(sname)
!   end if
!
!   end subroutine magnet

!   add irmax to the 1st argument, by Yang Wang
!=====================================================================
   subroutine magnetic_dens(irmax, gz, fz, gj, fj, nuz, indz, tau00, gcoeff, gbcoeff, out)
!======================================================================

   use RelativityToolsModule, only : repl

   implicit none

   character (len=32), parameter :: sname = 'magnetic_dens'

   integer (kind=IntKind) :: kmymax
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: kmy, kmyp
   integer (kind=IntKind) :: nu, nup
   integer (kind=IntKind) :: kmy1, kmyp1
   integer (kind=IntKind), intent(in) :: irmax
   integer (kind=IntKind), intent(in) :: nuz(:), indz(:,:)

   complex (kind=CmplxKind) :: cxr, cxi
   complex (kind=CmplxKind), intent(in) :: gcoeff(:,:)
   complex (kind=CmplxKind), intent(in) :: gbcoeff(:,:)
   complex (kind=CmplxKind), intent(in) :: gz(:,:,:), fz(:,:,:)
   complex (kind=CmplxKind), intent(in) :: gj(:,:,:), fj(:,:,:)
   complex (kind=CmplxKind), intent(in) :: tau00(kkrsz*2,kkrsz*2)
   complex (kind=CmplxKind) :: gf(kkrsz*2,kkrsz*2)
   complex (kind=CmplxKind) :: sum1(kkrsz*2,kkrsz*2)
   complex (kind=CmplxKind) :: sum2(kkrsz*2,kkrsz*2)
   complex (kind=CmplxKind), intent(out) :: out(:)
   complex (kind=CmplxKind) :: gcox, gbcox

   kmymax=2*(lmax+1)*(lmax+1)

!  do i=1,jmt       ! look up to irmax istead of jmt. By Yang Wang
   do i=1,irmax
      do kmyp=1,kmymax
         do kmy=1,kmymax
            cxr = CZERO
            cxi = CZERO
            do nu=1,nuz(kmy)
               kmy1=indz(nu,kmy)
               do nup=1,nuz(kmyp)
                  kmyp1=indz(nup,kmyp)
                  gcox=gcoeff(kmy1,kmyp1)
                  gbcox=gbcoeff(kmy1,kmyp1)
                  if(abs(gcox).gt.small) then
                     cxr=cxr+(gz(i,nu,kmy)*gz(i,nup,kmyp))*gcox
                     cxi=cxi+(gz(i,nu,kmy)*gj(i,nup,kmyp))*gcox
                  end if
                  if(abs(gbcox).gt.small) then
                     cxr=cxr-(fz(i,nu,kmy)*fz(i,nup,kmyp))*gbcox
                     cxi=cxi-(fz(i,nu,kmy)*fj(i,nup,kmyp))*gbcox
                  end if
               end do ! nup
            end do   ! nu
            sum1(kmy,kmyp) = cxr
            sum2(kmy,kmyp) = cxi
         end do ! kmyp
      end do   ! kmy

      call repl(gf,sum2,kmymax,kmymax)

      call zgemm('N','N',kmymax,kmymax,kmymax,CONE,sum1,kmymax, &
                 tau00,kmymax,-CONE,gf,kmymax)

      out(i) = CZERO
      do kmy=1,kmymax
         out(i) = out(i) + gf(kmy,kmy)
      end do

   end do ! i

   if(istop.eq.sname) then
      do i=1,irmax
         write(6,*) i,out(i)
      end do
      call StopHandler(sname)
   end if

   end subroutine magnetic_dens

!=========================================================================
   subroutine new_dens(ns,gz,fz,gj,fj,nuz,indz, &
                       tau,zrho,r_mesh)
!=========================================================================

! output: zrho - radial density distribution

   use RelativityToolsModule, only : repl, doubmt, submat

   implicit none

   character (len=32), parameter :: sname = 'dens'

   real (kind=RealKind) :: fac
   real (kind=RealKind) :: rlam
   real (kind=RealKind) :: r_mesh(:)

   integer (kind=IntKind) :: lmmax,kmymax
   integer (kind=IntKind), intent(in) :: nuz(:), indz(:,:)
   integer (kind=IntKind), intent(in) :: ns
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: kmy, kmyp
   integer (kind=IntKind) :: nu, nup
   integer (kind=IntKind) :: kmy1, kmyp1

   complex (kind=CmplxKind), intent(in) :: gz(:,:,:), fz(:,:,:)
   complex (kind=CmplxKind), intent(in) :: gj(:,:,:), fj(:,:,:)
   complex (kind=CmplxKind), intent(out) :: zrho(:)
   complex (kind=CmplxKind), intent(in) :: tau(:,:)
   complex (kind=CmplxKind) :: sum(2*kkrsz,2*kkrsz,2)
   complex (kind=CmplxKind) :: gf(2*kkrsz,2*kkrsz)
   complex (kind=CmplxKind) :: cox

   fac=-ONE/PI
   lmmax=(lmax+1)*(lmax+1)
   kmymax=2*lmmax
!  call zeroout(zrho,2*ns)

   do i=1,ns
      call zeroout(sum,2*2*kmymax*kmymax)
      do kmyp=1,kmymax
         do kmy=1,kmymax
            do nup=1,nuz(kmyp)
               kmyp1=indz(nup,kmyp)
               do nu=1,nuz(kmy)
                  kmy1=indz(nu,kmy)
                  cox=rgacoeff(kmy1,kmyp1,1)
                  if(cdabs(cox).gt.small) then
!                    rlam=r_mesh(i)**lambda
                     rlam=ONE
                     sum(kmy,kmyp,1)=sum(kmy,kmyp,1)+ &
                        (gz(i,nu,kmy)*gz(i,nup,kmyp)+ &
                         fz(i,nu,kmy)*fz(i,nup,kmyp))*cox*rlam
                     sum(kmy,kmyp,2)=sum(kmy,kmyp,2)+ &
                        (gz(i,nu,kmy)*gj(i,nup,kmyp)+ &
                         fz(i,nu,kmy)*fj(i,nup,kmyp))*cox*rlam
                  end if
               end do  ! nup
            end do    ! nu
         end do      ! kmyp
      end do      ! kmy

      call repl(gf,sum(1:2*kkrsz,1:2*kkrsz,1),kmymax,kmymax)
      call doubmt(gf,tau,kmymax,kmymax)
      call submat(gf,sum(1:2*kkrsz,1:2*kkrsz,2),kmymax,kmymax)

      zrho(i)=CZERO
      do kmy=1,kmymax
         zrho(i)=zrho(i)+gf(kmy,kmy)
      end do

   end do
!
   end subroutine new_dens

!============================================================
   subroutine gfill
!============================================================

! matrices of spin and orbital momentum operators
! on the basis of the spinor spherical harmonics

   use DiracSolverModule, only : getCGCu1, getCGCu2

   implicit none

! Clebsch-Gordan coeffs
   real (kind=RealKind), pointer :: u1(:)
   real (kind=RealKind), pointer :: u2(:)

   complex (kind=CmplxKind) :: sx(98,98), sxb(98,98)
   complex (kind=CmplxKind) :: sy(98,98), syb(98,98)
   complex (kind=CmplxKind) :: sz(98,98), szb(98,98)
   complex (kind=CmplxKind) :: jx(98,98), jxb(98,98)
   complex (kind=CmplxKind) :: jy(98,98), jyb(98,98)
   complex (kind=CmplxKind) :: jz(98,98), jzb(98,98)

   integer (kind=IntKind) :: l1, l2
   integer (kind=IntKind) :: j1, j2
   integer (kind=IntKind) :: my1, my2
   integer (kind=IntKind) :: kap1, kap2
   integer (kind=IntKind) :: kmy1, kmy2
   integer (kind=IntKind) :: kapp1, kapp2
   integer (kind=IntKind) :: kmyp1, kmyp2

   real (kind=RealKind) :: gamm, gamp
   real (kind=RealKind) :: x, y
   real (Kind=RealKind) :: c1d, c1u, c2d, c2u

   integer (kind=IntKind) :: kapdex(98), ldex(98)
   integer (kind=IntKind) :: jdex(98), mydex(98)

   data kapdex/                          &
    -1,-1,                               &
     1, 1,                               &
    -2,-2,-2,-2,                         &
     2, 2, 2, 2,                         &
    -3,-3,-3,-3,-3,-3,                   &
     3, 3, 3, 3, 3, 3,                   &
    -4,-4,-4,-4,-4,-4,-4,-4,             &
     4, 4, 4, 4, 4, 4, 4, 4,             &
    -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,       &
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       &
    -6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6,-6, &
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
    -7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7/
   data ldex/                            &
     0, 0,                               &
     1, 1,                               &
     1, 1, 1, 1,                         &
     2, 2, 2, 2,                         &
     2, 2, 2, 2, 2, 2,                   &
     3, 3, 3, 3, 3, 3,                   &
     3, 3, 3, 3, 3, 3, 3, 3,             &
     4, 4, 4, 4, 4, 4, 4, 4,             &
     4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       &
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       &
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6/
   data jdex/                            &
     1, 1,                               &
     1, 1,                               &
     3, 3, 3, 3,                         &
     3, 3, 3, 3,                         &
     5, 5, 5, 5, 5, 5,                   &
     5, 5, 5, 5, 5, 5,                   &
     7, 7, 7, 7, 7, 7, 7, 7,             &
     7, 7, 7, 7, 7, 7, 7, 7,             &
     9, 9, 9, 9, 9, 9, 9, 9, 9, 9,       &
     9, 9, 9, 9, 9, 9, 9, 9, 9, 9,       &
    11,11,11,11,11,11,11,11,11,11,11,11, &
    11,11,11,11,11,11,11,11,11,11,11,11, &
    13,13,13,13,13,13,13,13,13,13,13,13,13,13/
   data mydex/                                        &
    -1, 1,                                            &
    -1, 1,                                            &
    -3,-1, 1, 3,                                      &
    -3,-1, 1, 3,                                      &
    -5,-3,-1, 1, 3, 5,                                &
    -5,-3,-1, 1, 3, 5,                                &
    -7,-5,-3,-1, 1, 3, 5, 7,                          &
    -7,-5,-3,-1, 1, 3, 5, 7,                          &
    -9,-7,-5,-3,-1, 1, 3, 5, 7, 9,                    &
    -9,-7,-5,-3,-1, 1, 3, 5, 7, 9,                    &
    -11, -9, -7, -5, -3, -1,  1,  3,  5,  7,  9,  11, &
    -11, -9, -7, -5, -3, -1,  1,  3,  5,  7,  9,  11, &
    -13,-11, -9, -7, -5, -3, -1,  1,  3,  5,  7,  9,  11, 13/

   call zeroout(sx,98*98*2)
   call zeroout(sy,98*98*2)
   call zeroout(sz,98*98*2)
   call zeroout(sxb,98*98*2)
   call zeroout(syb,98*98*2)
   call zeroout(szb,98*98*2)
   call zeroout(jx,98*98*2)
   call zeroout(jy,98*98*2)
   call zeroout(jz,98*98*2)
   call zeroout(jxb,98*98*2)
   call zeroout(jyb,98*98*2)
   call zeroout(jzb,98*98*2)

   u1 => getCGCu1()
   u2 => getCGCu2()

   do kmy2=1,98
      c2d=u1(kmy2)
      c2u=u2(kmy2)
      kap2=kapdex(kmy2)
      l2=ldex(kmy2)
      j2=jdex(kmy2)
      my2=mydex(kmy2)
      x=j2*(j2+2.0d0)-my2*(my2+2.0d0)
      y=j2*(j2+2.0d0)-my2*(my2-2.0d0)
      if(x.lt.small) then
         gamp=ZERO
      else
         gamp=0.25d0*sqrt(x)
      end if
      if(y.lt.small) then
         gamm=ZERO
      else
         gamm=0.25*sqrt(y)
      end if
      do 10 kmy1=1,98
         c1d=u1(kmy1)
         c1u=u2(kmy1)
         kap1=kapdex(kmy1)
         l1=ldex(kmy1)
         j1=jdex(kmy1)
         my1=mydex(kmy1)

         if(l1.ne.l2) goto 10

! 'x' and 'y' matrices

         if(my1.eq.my2+2) then
            sx(kmy1,kmy2)=c1u*c2d
            sy(kmy1,kmy2)=-sqrtm1*c1u*c2d
            if(j1.eq.j2) then
               jx(kmy1,kmy2)=gamp
               jy(kmy1,kmy2)=-sqrtm1*gamp
            end if
         else if(my1.eq.my2-2) then
            sx(kmy1,kmy2)=c1d*c2u
            sy(kmy1,kmy2)=sqrtm1*c1d*c2u
            if(j1.eq.j2) then
               jx(kmy1,kmy2)=gamm
               jy(kmy1,kmy2)=sqrtm1*gamm
            end if
         end if

! 'z' matrices

         if(my1.eq.my2) then
            sz(kmy1,kmy2)=c1u*c2u-c1d*c2d
            if(j1.eq.j2) jz(kmy1,kmy2)=0.5d0*my1
         end if

  10  continue 
      enddo
 
   do kmy2=1,72
      kap2=kapdex(kmy2)
      j2=jdex(kmy2)
      my2=mydex(kmy2)
      kapp2=-kap2
      kmyp2=2*kapp2*kapp2+kapp2+(my2+1)/2
      do kmy1=1,72
         kap1=kapdex(kmy1)
         j1=jdex(kmy1)
         my1=mydex(kmy1)
         kapp1=-kap1
         kmyp1=2*kapp1*kapp1+kapp1+(my1+1)/2
         sxb(kmy1,kmy2)=sx(kmyp1,kmyp2)
         syb(kmy1,kmy2)=sy(kmyp1,kmyp2)
         szb(kmy1,kmy2)=sz(kmyp1,kmyp2)
         jxb(kmy1,kmy2)=jx(kmyp1,kmyp2)
         jyb(kmy1,kmy2)=jy(kmyp1,kmyp2)
         jzb(kmy1,kmy2)=jz(kmyp1,kmyp2)
      end do
   end do

   do kmy2=1,kmymaxp
      do kmy1=1,kmymaxp
        sxcoeff(kmy1,kmy2)=sx(kmy1,kmy2)
        sycoeff(kmy1,kmy2)=sy(kmy1,kmy2)
        szcoeff(kmy1,kmy2)=sz(kmy1,kmy2)
        sxbcoeff(kmy1,kmy2)=sxb(kmy1,kmy2)
        sybcoeff(kmy1,kmy2)=syb(kmy1,kmy2)
        szbcoeff(kmy1,kmy2)=szb(kmy1,kmy2)
        lxcoeff(kmy1,kmy2)=jx(kmy1,kmy2)-0.5d0*sx(kmy1,kmy2)
        lycoeff(kmy1,kmy2)=jy(kmy1,kmy2)-0.5d0*sy(kmy1,kmy2)
        lzcoeff(kmy1,kmy2)=jz(kmy1,kmy2)-0.5d0*sz(kmy1,kmy2)
        lxbcoeff(kmy1,kmy2)=jxb(kmy1,kmy2)-0.5d0*sxb(kmy1,kmy2)
        lybcoeff(kmy1,kmy2)=jyb(kmy1,kmy2)-0.5d0*syb(kmy1,kmy2)
        lzbcoeff(kmy1,kmy2)=jzb(kmy1,kmy2)-0.5d0*szb(kmy1,kmy2)
      end do
   end do

   end subroutine gfill

!=====================================================================
   subroutine gafill
!=====================================================================

   use RelativityToolsModule, only : getclm, plglmax, gaunt

   implicit none

   character (len=32), parameter :: sname='gafill'

   integer (kind=IntKind) :: i, l, m
   integer (kind=IntKind) :: lam, lp, lmp
   integer (kind=IntKind) :: lpmax, lpmin
   integer (kind=IntKind) :: lm, mp
   integer (kind=IntKind) :: nu
   integer (kind=IntKind) :: lmmaxp
   integer (kind=IntKind) :: kkr1

   real (kind=RealKind) :: plmg((2*lmax_max+1)*(2*lmax_max+2)/2,2*lmax_max+1)
   real (kind=RealKind) :: clm((2*lmax_max+1)*(2*lmax_max+2)/2)
   real (kind=RealKind) :: tg(2*(2*lmax_max+1))
   real (kind=RealKind) :: wg(2*(2*lmax_max+1))
   real (kind=RealKind) :: fac

!     set up things for the gaunt coeff:s
   call getclm(2*lmax_max,clm)
   call gauleg(-1.d0,1.d0,tg,wg,2*(2*lmax_max+1))
   do i=1,2*lmax_max+1
      call plglmax(2*lmax_max,tg(i),plmg(:,i))
   end do

   lmmaxp =(lmax_max+1)*(lmax_max+1)

   i=0
   do lam=0,1  ! originally "do lam=0,lamp"
      fac=TWO*sqrt_pi/(2*lam+1)
      do nu=-lam,lam
         i=i+1
         call zeroout(gacoeff(1,1,i),2*lmmaxp*lmmaxp)
         call zeroout(rgacoeff(1,1,i),2*kmymaxp*kmymaxp)
         do l=0,lmax_max
            do m=-l,l
               lm=l*l+l+m+1
               mp=m-nu
               lpmin=iabs(l-lam)
               lpmax=l+lam
               do 1 lp=lpmin,lpmax,2
                  if(lp.gt.lmax_max) goto 1
                  if(iabs(mp).gt.lp) goto 1
                  lmp=lp*lp+lp+mp+1
                  gacoeff(lm,lmp,i)=fac*gaunt(lam,nu,l,m,lp,mp, &
                                    2*lmax_max,2*lmax_max+1,wg,plmg,clm)
    1          continue
            end do
         end do
         kkr1=(lmax_max+1)*(lmax_max+1)
         call relmtrx(gacoeff(:,:,i),rgacoeff(:,:,i),kkr1,kkr1)

      end do
   end do

   end subroutine gafill

!=====================================================================
   subroutine relmtrx(a,b,kkr1,kkr2)
!=====================================================================

! *****************************************************
! * transformation from a non-relativistic matrix 'a' *
! *                to a relativistic matrix 'b'       *
! *****************************************************

   use DiracSolverModule, only : getCGCu1, getCGCu2, &
                                 getCGCind1, getCGCind2

   implicit none

   integer (kind=IntKind), intent(in) :: kkr1, kkr2
   integer (kind=IntKind) :: i, i1, i2, j, j1, j2
   integer (kind=IntKind), pointer :: ind1(:), ind2(:)
   real (kind=RealKind), pointer :: u1(:), u2(:)
   complex (kind=CmplxKind), intent(in) :: a(:,:)
   complex (kind=CmplxKind), intent(out) :: b(:,:)

   ind1 => getCGCind1()
   ind2 => getCGCind2()
   u1 => getCGCu1()
   u2 => getCGCu2()

   do j=1,2*kkr2
      j1=ind1(j)
      j2=ind2(j)
      do i=1,2*kkr1
         i1=ind1(i)
         i2=ind2(i)
         b(i,j)=u1(i)*a(i1,j1)*u1(j)+ &
                u2(i)*a(i2,j2)*u2(j)
      end do
   end do

   end subroutine relmtrx

end module RelGreenFunctionModule
