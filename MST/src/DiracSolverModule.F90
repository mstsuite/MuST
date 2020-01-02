module DiracSolverModule

!**************************************************************************
! Name: DiracSolverModule                                                 *
!                                                                         *
! Version 1.0: June 2011 by wkb                                           *
!                                                                         *
! Description:                                                            *
!    Solves the single-site Dirac equation and generates relativistic     *
!    t-matrices. Adapted from LSMS 1.9 routines.                          *
!                                                                         *
! ----------------------------------------------------------------------- *
! subroutine initDiracSolver(na, lmax, z, i_print, stop_routine)          *
!    Initializes variables and allocates the memory needed for the        *
!    routines in this module. Must be called first.                       *
!  Input: na - integer number of atoms on this processor                  *
!         lmax - 1D integer array, spin index cutoff for each local atom  *
!         i_print - integer, print output level                           *
!         stop_routine - character string, name of routine in which to    *
!                        end the program                                  *
! ----------------------------------------------------------------------- *
! subroutine endDiracSolver()                                             *
!    Deallocates the memory used in this module. This is the last routine *
!    that can be called unless initDiracSolver is called again.           *
! ----------------------------------------------------------------------- *
! subroutine solveSingleDirac(ia, energy)                                 *
!    Calculates the t-matrix and solutions to the Dirac equation for a    *
!    given atom and stores them in memory.                                *
!  Input: ia - integer, local index of the atom to be solved              *
!         energy - complex energy                                         *
! ----------------------------------------------------------------------- *
! function getRelTMatrix(atom, e) result(t_mat)                           *
!    Returns the relativistic t-matrix (kappa,mu form) for a given atom.  *
!    solveSingleDirac must first be called on that atom.                  *
!  Input: atom - integer, local atom index                                *
!         e - complex energy                                              *
!  Output: t_mat - pointer to 2D complex array containing relativistic    *
!                  t-matrix                                               *
! ----------------------------------------------------------------------- *
! function getGReg(atom, e) result(gz)                                    *
!    Returns the large component of the regular solution to the Dirac     *
!    equation for a given atom.                                           *
!  Input: same as above                                                   *
!  Output: gz - pointer to 3D complex array containing the large          *
!               component of the regular solution                         *            
! ----------------------------------------------------------------------- *
! function getFReg(atom, e) result(fz)                                    *
!    Returns the small component of the regular solution to the Dirac     *
!    equation for a given atom.                                           *
!  Input: same as above                                                   *
!  Output: fz - pointer to 3D complex array containing the small          *
!               component of the regular solution                         *
! ----------------------------------------------------------------------- *
! function getGIrr(atom, e) result(gj)                                    *
!    Returns the large component of the irregular solution to the Dirac   *
!    equation for a given atom.                                           *
!  Input: same as above                                                   *
!  Output: gj - pointer to 3D complex array containing the large          *
!               component of the irregular solution                       *
! ----------------------------------------------------------------------- *
! function getFIrr(atom, e) result(fj)                                    *
!    Returns the small component of the irregular solution to the Dirac   *
!    equation for a given atom.                                           *
!  Input: same as above                                                   *
!  Output: fj - pointer to 3D complex array containing the small          *
!               componentof the irregular solution                        *
! ----------------------------------------------------------------------- *
! function getnuz(atom, e) result(nuz)                                    *
!    Returns the number of (kappa',mu') components for (kappa,mu).        *
!  Input: same as above                                                   *
!  Output: nuz - pointer to 1D integer array containing the number of     *
!                (kappa',mu') components for a given atom                 *
! ----------------------------------------------------------------------- *
! function getindz(atom, e) result(indz)                                  *
!    Returns the (kappa',mu') indexes for (kappa,mu).                     *
!  Input: same as above                                                   *
!  Output: nuz - pointer to 2D integer array containing the               *
!                (kappa',mu') indexes for a given atom                    *
! ----------------------------------------------------------------------- *
! function getCGCu1() result(cgc1)                                        *
!    Returns the Clebsch-Gordan coefficient array u1.                     *
!  Output: cgc1 - pointer to 1D real array containing the Clebsch-Gordan  *
!                 coefficients for a given atom (spin down)               *
! ----------------------------------------------------------------------- *
! function getCGCu2() result(cgc2)                                        *
!    Returns the Clebsch-Gordan coefficient array u2.                     *
!  Input: same as above                                                   *
!  Output: cgc2 - pointer to 2D real array containing the Clebsch-Gordan  *
!                 coefficients for a given atom (spin up)                 *
! ----------------------------------------------------------------------- *
! function getCGCind1() result(pind1)                                     *
!    Returns the Clebsch-Gordan index array ind1.                         *
!  Input: same as above                                                   *
!  Output: pind1 - pointer to 1D integer array containing the             *
!                  Clebsch-Gordan indexes for a given atom (spin down)    *
! ----------------------------------------------------------------------- *
! function getCGCind2() result(pind2)                                     *
!    Returns the Clebsch-Gordan index array ind2.                         *
!  Input: same as above                                                   *
!  Output: pind2 - pointer to 1D integer array containing the             *
!                  Clebsch-Gordan indexes for a given atom (spin up)      *
! ----------------------------------------------------------------------- *
! function getSizeOfRMesh(atom) result(nr)                                *
!    Returns the number of points in a given atom's radial mesh.          *
!  Input: atom - integer, local atom index                                *
!  Output: nr - integer, number of radial points on the atom's mesh       *
! ----------------------------------------------------------------------- *
! function getDmat(atom) result(dmat)                                     *
!    Returns the matrix to transform between the global frame and an      *
!    atom's local frame.                                                  *
!  Input: atom - integer, local atom index                                *
!  Output: dmat - pointer to 2D complex array containing the              *
!                 transformation matrix for the given atom                *
! ----------------------------------------------------------------------- *
! function getDmatP(atom) result(dmatp)                                   *
!    Returns the transposed conjugate of the matrix to transform between  *
!    the global frame and an atom's local frame.                          *
!  Input: atom - integer, local atom index                                *
!  Output: dmat - pointer to 2D complex array containing the transposed   *
!                 conjugate of the transformation matrix for the          *
!                 given atom                                              *
!**************************************************************************

use ErrorHandlerModule, only : ErrorHandler, StopHandler
use KindParamModule, only : IntKind, CmplxKind, RealKind
use MathParamModule, only : ZERO, HALF, ONE, TWO, PI, SIX, CZERO, TEN2m12, SQRTm1
use PhysParamModule, only : LightSpeed
use PublicTypeDefinitionsModule, only : GridStruct

public :: initDiracSolver,  &
          endDiracSolver,   &
          solveSingleDirac, &
          getRelTMatrix,    &
          getGReg,          &
          getFReg,          &
          getGIrr,          &
          getFIrr,          &
          getnuz,           &
          getindz,          &
          getCGCu1,         & ! These 4 could be calculated elsewhere
          getCGCu2,         &
          getCGCind1,       &
          getCGCind2,       &
          getSizeOfRMesh,   &
          getDmat,          &
          getDmatP

private

integer (kind=IntKind) :: LocalNumAtoms
integer (kind=IntKind), allocatable :: AtomicNumber(:)
integer (kind=IntKind), parameter :: n_spin_pola = 2, n_spin_cant = 2
integer (kind=IntKind), parameter :: nuzp = 2
integer (kind=IntKind) :: iplmax
integer (kind=IntKind) :: iprpts
integer (kind=IntKind) :: currentAtom

type RelScatterStruct
   logical :: done
   integer (kind=IntKind) :: l_max
   integer (kind=IntKind) :: NumRs
   integer (kind=IntKind) :: NumRs_cs
   integer (kind=IntKind), pointer :: nuz(:)
   integer (kind=IntKind), pointer :: indz(:,:)
   real (kind=RealKind) :: r_end ! in general, should = r_ws
   real (kind=RealKind) :: h
   complex (kind=CmplxKind), pointer :: gz(:,:,:)
   complex (kind=CmplxKind), pointer :: fz(:,:,:)
   complex (kind=CmplxKind), pointer :: gj(:,:,:)
   complex (kind=CmplxKind), pointer :: fj(:,:,:)
   complex (kind=CmplxKind), pointer :: t_mat(:,:) 
   complex (kind=CmplxKind), pointer :: d_mat(:,:)
   complex (kind=CmplxKind), pointer :: d_matp(:,:)
end type RelScatterStruct

type(RelScatterStruct), allocatable :: Scatter(:)
type(GridStruct), pointer :: Grid

logical :: Initialized = .false.

real (kind=RealKind) :: rs

character (len=32) :: istop

complex (kind=CmplxKind) :: last_energy
complex (kind=CmplxKind), allocatable, target :: wks_dmat(:)
complex (kind=CmplxKind), allocatable, target :: wks_dmatp(:)
complex (kind=CmplxKind), allocatable, target :: wks_gz(:)
complex (kind=CmplxKind), allocatable, target :: wks_fz(:)
complex (kind=CmplxKind), allocatable, target :: wks_gj(:)
complex (kind=CmplxKind), allocatable, target :: wks_fj(:)
complex (kind=CmplxKind), allocatable, target :: wks_tmat(:)
integer (kind=IntKind), allocatable, target :: wks_nuz(:)
integer (kind=IntKind), allocatable, target :: wks_indz(:)

! for Clebsch-Gordan coeffs:
real (kind=RealKind), allocatable, target :: u1(:), u2(:)
integer (kind=IntKind), allocatable, target :: ind1(:), ind2(:)

contains

include '../lib/arrayTools.F90'

!==========================================================================
   subroutine initDiracSolver(na, lmax, z, i_print, stop_routine)
!==========================================================================

! Relativistic solver currently supports only spherically symmetric potentials.
! lmax_phi and lmax_pot need to be added as arguments to implement full solver.

   use AtomModule, only : getLocalEvecOld
   use RadialGridModule, only : getGrid

   implicit none

   character (len=*), intent(in) :: stop_routine
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: lmax(:)
   integer (kind=IntKind), intent(in) :: z(:)
   integer (kind=IntKind), intent(in) :: i_print(:)

   integer (kind=IntKind) :: ia
   integer (kind=IntKind) :: iend
   integer (kind=IntKind) :: kkrsz
   integer (kind=IntKind) :: sz_ind_solutions=0, sz_ind_tmat=0
   integer (kind=IntKind) :: sz_ind_nuz=0, sz_ind_indz=0
   integer (kind=IntKind) :: ind_solutions(na), ind_tmat(na)
   integer (kind=IntKind) :: ind_nuz(na), ind_indz(na)
   integer (kind=IntKind) :: inda, indb
   integer (kind=IntKind), pointer :: pi1(:)

   real (kind=RealKind) :: evec_l(3), r_global(3)
!
   complex (kind=CmplxKind), pointer :: p1(:)

   if (Initialized) then
      call ErrorHandler('initDiracSolver', 'module is already initialized')
   else if (na < 1) then
      call ErrorHandler('initDiracSolver', 'LocalNumAtoms < 1', na)
   endif

   LocalNumAtoms = na
   last_energy = CZERO
   istop = stop_routine
   r_global(1) = ZERO
   r_global(2) = ZERO
   r_global(3) = ONE

   allocate( AtomicNumber(LocalNumAtoms) )
   allocate( Scatter(LocalNumAtoms) )
   allocate( ind1(500) )
   allocate( ind2(500) )
   allocate( u1(500) )
   allocate( u2(500) )

   iprpts = 0
   sz_ind_solutions=0; sz_ind_tmat=0
   sz_ind_nuz=0; sz_ind_indz=0
   do ia=1, LocalNumAtoms
      kkrsz = (lmax(ia)+1)*(lmax(ia)+1)
      AtomicNumber(ia) = z(ia)
      Grid => getGrid(ia)
      if (Grid%nmult .ne. 1) then
         call ErrorHandler('initDiracSolver', &
              'nmult must equal 1 for relativistic calculation')
      endif
      iend = max(Grid%jmt+11, Grid%jend)
      iend = min(iend, Grid%jend_plus_n)
      Scatter(ia)%done = .false.
      Scatter(ia)%l_max = lmax(ia)
      Scatter(ia)%NumRs = Grid%jend ! iend
      Scatter(ia)%NumRs_cs = Grid%jend
      Scatter(ia)%r_end = Grid%rend
      Scatter(ia)%h = Grid%hin
      ind_solutions(ia) = sz_ind_solutions + 1
      ind_tmat(ia) = sz_ind_tmat + 1
      ind_nuz(ia) = sz_ind_nuz + 1
      ind_indz(ia) = sz_ind_indz + 1
      sz_ind_solutions = sz_ind_solutions + (2*Scatter(ia)%NumRs*nuzp*kkrsz)
      sz_ind_tmat = sz_ind_tmat + (kkrsz*kkrsz*n_spin_cant*n_spin_cant)
      sz_ind_nuz = sz_ind_nuz + (2*kkrsz)
      sz_ind_indz = sz_ind_indz + (2*nuzp*kkrsz)
      iprpts = max(iprpts,Scatter(ia)%NumRs + 30)
   enddo

   allocate( wks_gz(sz_ind_solutions) ); wks_gz = CZERO
   allocate( wks_fz(sz_ind_solutions) ); wks_fz = CZERO
   allocate( wks_gj(sz_ind_solutions) ); wks_gj = CZERO
   allocate( wks_fj(sz_ind_solutions) ); wks_fj = CZERO
   allocate( wks_dmat(sz_ind_tmat) ); wks_dmat = CZERO
   allocate( wks_dmatp(sz_ind_tmat) ); wks_dmatp = CZERO
   allocate( wks_tmat(sz_ind_tmat) ); wks_tmat = CZERO
   allocate( wks_nuz(sz_ind_nuz) ); wks_nuz = 0
   allocate( wks_indz(sz_ind_indz) ); wks_indz = 0

   iplmax = 0
   do ia=1, LocalNumAtoms
      iplmax = max(lmax(ia),iplmax)
      kkrsz = (lmax(ia)+1)*(lmax(ia)+1)
      inda = ind_solutions(ia)
      indb = inda + (2*Scatter(ia)%NumRs*nuzp*kkrsz) - 1
      p1 => wks_gz(inda:indb)
      Scatter(ia)%gz => aliasArray3_c(p1, Scatter(ia)%NumRs, nuzp, 2*kkrsz)
      p1 => wks_fz(inda:indb)
      Scatter(ia)%fz => aliasArray3_c(p1, Scatter(ia)%NumRs, nuzp, 2*kkrsz)
      p1 => wks_gj(inda:indb)
      Scatter(ia)%gj => aliasArray3_c(p1, Scatter(ia)%NumRs, nuzp, 2*kkrsz)
      p1 => wks_fj(inda:indb)
      Scatter(ia)%fj => aliasArray3_c(p1, Scatter(ia)%NumRs, nuzp, 2*kkrsz)
      inda = ind_tmat(ia)
      indb = inda + (kkrsz*kkrsz*n_spin_cant*n_spin_cant) - 1
      p1 => wks_tmat(inda:indb)
      Scatter(ia)%t_mat => aliasArray2_c(p1, kkrsz*n_spin_cant, kkrsz*n_spin_cant)
      p1 => wks_dmat(inda:indb)
      Scatter(ia)%d_mat => aliasArray2_c(p1, kkrsz*n_spin_cant, kkrsz*n_spin_cant)
      p1 => wks_dmatp(inda:indb)
      Scatter(ia)%d_matp => aliasArray2_c(p1, kkrsz*n_spin_cant, kkrsz*n_spin_cant)
      evec_l(1:3) = getLocalEvecOld(ia)
      call matrot1(r_global, evec_l, lmax(ia), Scatter(ia)%d_mat, Scatter(ia)%d_matp)
      inda = ind_nuz(ia)
      indb = inda + (2*kkrsz) - 1
      pi1 => wks_nuz(inda:indb)
      Scatter(ia)%nuz => aliasArray1_i(pi1, 2*kkrsz)
      inda = ind_indz(ia)
      indb = inda + (2*nuzp*kkrsz) - 1
      pi1 => wks_indz(inda:indb)
      Scatter(ia)%indz => aliasArray2_i(pi1, nuzp, 2*kkrsz)
   enddo

   call clebsch

   nullify( Grid )

   Initialized = .true.

   end subroutine initDiracSolver

!==========================================================================
   subroutine endDiracSolver()
!==========================================================================

   implicit none

   integer (kind=IntKind) :: ia

   if (.not. Initialized) then
      call ErrorHandler('endDiracSolver', 'module not initialized')
   endif

   Initialized = .false.

   deallocate( AtomicNumber )
   deallocate( ind1 )
   deallocate( ind2 )
   deallocate( u1 )
   deallocate( u2)

   do ia=1, LocalNumAtoms
      Scatter(ia)%done = .false.
      nullify( Scatter(ia)%gz )
      nullify( Scatter(ia)%fz )
      nullify( Scatter(ia)%gj )
      nullify( Scatter(ia)%fj )
      nullify( Scatter(ia)%t_mat )
      nullify( Scatter(ia)%nuz )
      nullify( Scatter(ia)%indz )
   enddo
   deallocate( Scatter )
   deallocate( wks_gz )
   deallocate( wks_fz )
   deallocate( wks_gj )
   deallocate( wks_fj )
   deallocate( wks_tmat )
   deallocate( wks_dmat )
   deallocate( wks_dmatp )
   deallocate( wks_nuz )
   deallocate( wks_indz )

   end subroutine endDiracSolver

!==========================================================================
   subroutine solveSingleDirac(ia, energy)
!==========================================================================

! tmat_g is the t-matrix in the global frame of reference.

   use AtomModule, only : getLocalEvecOld
   use PhysParamModule, only : LightSpeed
   use PotentialModule, only : getSphPotr
   use RelativityToolsModule, only : gjinv, tripmt, outmat1

   implicit none

   character (len=32), parameter :: sname='solveSingleDirac'
   character (len=10) :: idpot

   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind) :: jws
   integer (kind=IntKind) :: iend
   integer (kind=IntKind) :: lmax
   integer (kind=IntKind) :: kkrsz
   integer (kind=IntKind), pointer :: nuz(:)
   integer (kind=IntKind), pointer :: indz(:,:)
   integer (kind=IntKind) :: kmymax
   integer (kind=IntKind) ::  j
   integer (kind=IntKind) :: info
   integer (kind=IntKind) :: isq, isp

   real (kind=RealKind), pointer :: vrup(:), vrdown(:)
   real (kind=RealKind) :: v0
   real (kind=RealKind) :: dx
   real (kind=RealKind), allocatable :: vrr(:)
   real (kind=RealKind), allocatable :: brr(:)
   real (kind=RealKind), allocatable :: boprr(:,:)
   real (kind=RealKind), parameter :: soscal=1.d0 !spin-orbit scaling

   complex (kind=CmplxKind), pointer ::  tmat_g(:,:)
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind) :: psq
   complex (kind=CmplxKind), pointer :: gz(:,:,:)
   complex (kind=CmplxKind), pointer :: fz(:,:,:)
   complex (kind=CmplxKind), pointer :: gj(:,:,:)
   complex (kind=CmplxKind), pointer :: fj(:,:,:)
   complex (kind=CmplxKind), allocatable :: wbig(:)
   complex (kind=CmplxKind), pointer :: dmat(:,:)
   complex (kind=CmplxKind), pointer :: dmatp(:,:)
   complex (kind=CmplxKind) :: detl

   if(.not.Initialized) then
      call ErrorHandler('solveSingleDirac', 'module not initialized')
   endif

   jws = Scatter(ia)%NumRs_cs
   iend = Scatter(ia)%NumRs
   lmax = Scatter(ia)%l_max
   kkrsz = (lmax+1)*(lmax+1)
   rs = Scatter(ia)%r_end
   dx = Scatter(ia)%h
   currentAtom = ia

   vrup => getSphPotr(ia, 1, 1)
   vrdown => getSphPotr(ia, 1, 2)

   idpot='x123456789' !needs to be modified, have it passed in as argument
   v0=0.d0
   last_energy = energy

   allocate( vrr(iend) )
   allocate( brr(iend) )
   allocate( boprr(iend,2) )
   allocate( wbig(4*kkrsz*kkrsz) )

   gz => Scatter(ia)%gz(1:iend,1:nuzp,1:2*kkrsz)
   fz => Scatter(ia)%fz(1:iend,1:nuzp,1:2*kkrsz)
   gj => Scatter(ia)%gj(1:iend,1:nuzp,1:2*kkrsz)
   fj => Scatter(ia)%fj(1:iend,1:nuzp,1:2*kkrsz)
   nuz => Scatter(ia)%nuz(1:2*kkrsz)
   indz => Scatter(ia)%indz(1:nuzp,1:2*kkrsz)
   tmat_g => Scatter(ia)%t_mat(1:kkrsz*n_spin_cant,1:kkrsz*n_spin_cant)
   dmat => Scatter(ia)%d_mat(1:kkrsz*n_spin_cant,1:kkrsz*n_spin_cant)
   dmatp => Scatter(ia)%d_matp(1:kkrsz*n_spin_cant,1:kkrsz*n_spin_cant)
   kmymax = kkrsz*n_spin_cant

!  currently the relativistic part works only for ASA

   do j=1,iend
      vrr(j) = 0.5d0*(vrup(j)+vrdown(j))
      if(vrup(1).gt.-1.d-4) vrr(j)=vrr(j)-1.d-4
      brr(j) = 0.5d0*(vrup(j)-vrdown(j))
      boprr(j,1) = 0.d0
      boprr(j,2) = 0.d0
   end do
   200   format(4d20.12)
   psq = energy  +energy*energy/(LightSpeed*LightSpeed)

   call single_scatterer_rel(energy,psq,lmax,kmymax, &
             'x123456789',v0,                        &
             vrr,brr,boprr,dx,jws,                   &
             tmat_g,gz,fz,gj,fj,nuz,indz,1,soscal)

!  now we have tinv but we need t:

   call gjinv(tmat_g,kmymax,kmymax,detl)

!  t is in the local frame of reference, now we have to
!  transform it into the global frame:

   call tripmt(dmat,tmat_g,dmatp,kmymax,kmymax,kmymax)

   Scatter(ia)%done = .true.

   deallocate( vrr, brr )
   deallocate( boprr )
   deallocate( wbig )

   if(istop.eq.sname) then
      write(6,*) 'tmat_g:'
      call outmat1(tmat_g,kmymax,kmymax,kmymax,1.d-15,6)
      call gjinv(tmat_g,kmymax,kmymax,detl)
      call replms(wbig,tmat_g,lmax,kmymax)
      write(6,*) 'diag(tmat_lms):'
      do j=1,2*(lmax+1)**2
         write(6,*) j,j,wbig(j)
      end do
      write(6,*) 't^(-1) :'
      call outmat1(tmat_g,kmymax,kmymax,kmymax,1.d-15,6)
      call StopHandler(sname)
   end if

   end subroutine solveSingleDirac

!==================================================================
   subroutine single_scatterer_rel(ce,psq,lmax,kmymax, &
                   idpot,v0,vr,br,bopr,dx,ns,          &
                   tminv,gz,fz,gj,fj,nuz,indz,iflag,socsc)        
!==================================================================

! tminv is the inverse of the single-site scattering matrix

   use RelativityToolsModule, only : csbf

   implicit none

   character (len=10) :: idpot
   character (len=32), parameter :: sname='single_scatterer_rel'

   integer (kind=IntKind), intent(in) :: lmax, kmymax, ns, iflag
   integer (kind=IntKind), intent(out) :: nuz(:), indz(:,:)
   integer (kind=IntKind) :: k, l, kmax, kap, kap1, kap2, j, lb
   integer (kind=IntKind) :: my, kapmy, kapmy1, kapmy2, i

   real (kind=RealKind), intent(in) :: vr(:), br(:), v0, dx, socsc
   real (kind=RealKind), intent(in) :: bopr(:,:)

   complex (kind=CmplxKind) :: tminv2(2,2)
   complex (kind=CmplxKind), pointer :: tminv(:,:)
   complex (kind=CmplxKind) :: gz2(iprpts,2,2), fz2(iprpts,2,2)
   complex (kind=CmplxKind) :: gj2(iprpts,2,2), fj2(iprpts,2,2)
   complex (kind=CmplxKind), intent(out) :: gz(:,:,:) ,fz(:,:,:) ! dimension (iprpts,nuzp,kmymax)
   complex (kind=CmplxKind), intent(out) :: gj(:,:,:), fj(:,:,:)
   complex (kind=CmplxKind) :: fb(0:lmax+1), fn(0:lmax+1), fh(0:lmax+1)
   complex (kind=CmplxKind) :: fb1(0:lmax+1), fn1(0:lmax+1), fh1(0:lmax+1)
   complex (kind=CmplxKind), intent(in) :: ce, psq
   complex (kind=CmplxKind) :: p, sk, dk, xk, react, cevac, psqvac, pvac

   kmax=2*lmax+1
   call zeroout(tminv,2*kmymax*kmymax)

   if(idpot.eq.'Va        ') then
      p=cdsqrt(psq)
      cevac=ce-v0
      pvac=cdsqrt(cevac)
      call csbf(lmax+1,p,rs,fb,fn,fh)
      call csbf(lmax+1,pvac,rs,fb1,fn1,fh1)
      do k=1,kmax
         l=k/2
         if(k.eq.2*l) then
            kap=l
            lb=l-1
            j=2*l-1
         else
            kap=-l-1
            lb=l+1
            j=2*l+1
         end if
         sk=dcmplx(dfloat(l-lb),0.d0)
         xk=sk*cevac*fb1(lb)/pvac
         xk=xk/fb1(l)
         sk=sk*ce/p
         dk=(xk*fb(l)-sk*fb(lb))/(xk*fn(l)-sk*fn(lb))
         react=-dk/p
         do my=-j,j,2
            kapmy=2*kap*kap+kap+(my+1)/2
            tminv(kapmy,kapmy)=1.d0/react+sqrtm1*p
         end do
      end do
   else  
      do l=0,lmax
         kap1=l
         kap2=-l-1
         do my=-2*l-1,2*l+1,2
!---------------------------------------------------------------
            call spzwafu(socsc,ce,psq,l,my,vr,br,bopr,dx,ns, &
                         tminv2,gz2,fz2,gj2,fj2,iflag)
!---------------------------------------------------------------
            if(abs(my).eq.2*l+1) then
               kapmy=2*kap2*kap2+kap2+(my+1)/2
               tminv(kapmy,kapmy)=tminv2(2,2)
               if(iflag.eq.1) then
                  nuz(kapmy)=1
                  indz(1,kapmy)=kapmy
		  indz(2,kapmy)=0
                  do i=1,ns
                     gz(i,1,kapmy)=gz2(i,2,2)
                     fz(i,1,kapmy)=fz2(i,2,2)
                     gj(i,1,kapmy)=gj2(i,2,2)
                     fj(i,1,kapmy)=fj2(i,2,2)
                  end do
               end if
            else
               kapmy1=2*kap1*kap1+kap1+(my+1)/2
               kapmy2=2*kap2*kap2+kap2+(my+1)/2
               tminv(kapmy1,kapmy1)=tminv2(1,1)
               tminv(kapmy1,kapmy2)=tminv2(1,2)
               tminv(kapmy2,kapmy1)=tminv2(2,1)
               tminv(kapmy2,kapmy2)=tminv2(2,2)
               if(iflag.eq.1) then
                  nuz(kapmy1)=2
                  nuz(kapmy2)=2
                  indz(1,kapmy1)=kapmy1
                  indz(2,kapmy1)=kapmy2
                  indz(1,kapmy2)=kapmy2
                  indz(2,kapmy2)=kapmy1
                  do i=1,ns
                     gz(i,1,kapmy1)=gz2(i,1,1)
                     fz(i,1,kapmy1)=fz2(i,1,1)
                     gz(i,2,kapmy1)=gz2(i,2,1)
                     fz(i,2,kapmy1)=fz2(i,2,1)
                     gz(i,1,kapmy2)=gz2(i,2,2)
                     fz(i,1,kapmy2)=fz2(i,2,2)
                     gz(i,2,kapmy2)=gz2(i,1,2)
                     fz(i,2,kapmy2)=fz2(i,1,2)
                     gj(i,1,kapmy1)=gj2(i,1,1)
                     fj(i,1,kapmy1)=fj2(i,1,1)
                     gj(i,2,kapmy1)=gj2(i,2,1)
                     fj(i,2,kapmy1)=fj2(i,2,1)
                     gj(i,1,kapmy2)=gj2(i,2,2)
                     fj(i,1,kapmy2)=fj2(i,2,2)
                     gj(i,2,kapmy2)=gj2(i,1,2)
                     fj(i,2,kapmy2)=fj2(i,1,2)
                  end do
               end if
            end if
         end do
      end do
   end if

   if(istop.eq.sname) then
      call StopHandler(sname)
   end if

   end subroutine single_scatterer_rel

!===============================================================
   subroutine spzwafu(socsc,ce,psq,l,my,vr,br,bopr,dx,ns, &
                      tminv,gz,fz,gj,fj,iflag)
!===============================================================

! tminv is (l,my)-like 2x2 block of inverse t-matrix

   use RelativityToolsModule, only : gjinv, doubmt, submat, &
                                     repl, addmat, csbf

   implicit none

   logical scale

   character (len=32), parameter :: sname='spzwafu'

   integer (kind=IntKind) :: i, ii, jj
   integer (kind=IntKind), intent(in) :: l, my, ns, iflag

   real (kind=RealKind), intent(in) :: vr(:), br(:), socsc, dx
   real (kind=RealKind), intent(in) :: bopr(:,:)
   real (kind=RealKind) :: c, xl, xnot

   complex (kind=CmplxKind) :: tminv(:,:)
   complex (kind=CmplxKind) :: fz(:,:,:), gz(:,:,:)
   complex (kind=CmplxKind) :: fj(:,:,:), gj(:,:,:)
   complex (kind=CmplxKind) :: f1(2,2,iprpts),g1(2,2,iprpts),gp1(2,2)
   complex (kind=CmplxKind) :: f2(2,2,iprpts),g2(2,2,iprpts),gp2(2,2)
   complex (kind=CmplxKind) :: f11(iprpts),g11(iprpts),f12(iprpts),g12(iprpts)
   complex (kind=CmplxKind) :: gp11,gp12
   complex (kind=CmplxKind) :: fb(0:iplmax+1),fn(0:iplmax+1),fh(0:iplmax+1)
   complex (kind=CmplxKind) :: ce,psq,p,pl,sl,smlm1,detl
   complex (kind=CmplxKind) :: aa1,aa2,bb1,bb2,gg1,ff1,ggam1,gg2,ff2,ggam2,ggamb   
   complex (kind=CmplxKind) :: a1(2,2),a2(2,2),b1(2,2),b2(2,2)
   complex (kind=CmplxKind) :: g1mat(2,2),f1mat(2,2),gam1(2,2)
   complex (kind=CmplxKind) :: g2mat(2,2),f2mat(2,2),gam2(2,2)
   complex (kind=CmplxKind) :: jmat(2,2),jbmat(2,2),nmat(2,2),nbmat(2,2),gamb(2,2)
   complex (kind=CmplxKind) :: x1(2,2),x2(2,2),x3(2,2),x4(2,2)

   scale=abs(1.0d0-socsc).gt.1.0d-3

   call zeroout(gz,8*iprpts)
   call zeroout(fz,8*iprpts)
   call zeroout(gj,8*iprpts)
   call zeroout(fj,8*iprpts)
   call zeroout(tminv,8)
   call zeroout(jmat,8)
   call zeroout(jbmat,8)
   call zeroout(nmat,8)
   call zeroout(nbmat,8)
   call zeroout(gamb,8)

   xnot=dlog(rs)-(ns-1)*dx

!  c in rydberg units:
   c=274.072d0

   p=cdsqrt(psq)
   xl=dfloat(l)
   sl=1.d0
   smlm1=-1.d0
   sl=sl*ce/psq
   smlm1=smlm1*ce/psq

   call csbf(l+1,p,rs,fb,fn,fh)

!  if(abs(my)-(2*l+1)) 1,2,3
   if(abs(my)-(2*l+1) < 0) then
      goto 1
   else if (abs(my)-(2*l+1) == 0) then
      goto 2
   else
      goto 3
   endif
   3 stop ' spzwafu: my is out of range!'
   2 continue
!-----------------------------------------------------------------
   call dirmago1op(socsc,ce,l,my,vr,br,bopr,dx,xnot,ns, &
                   g11,f11,gp11)
!-----------------------------------------------------------------
   if(scale) then
      gg1=g11(ns)/rs
      ff1=gp11/rs
      ggam1=ff1/gg1
      aa2=p*((ggam1-xl/rs)*fn(l)+p*fn(l+1))/ &
            ((ggam1-xl/rs)*fb(l)+p*fb(l+1))
   else
      gg1=g11(ns)/rs
      ff1=f11(ns)/rs
      ggam1=ff1/gg1
      aa2=p*(ggam1*fn(l)-smlm1*p*fn(l+1)/c)/ &
            (ggam1*fb(l)-smlm1*p*fb(l+1)/c)
   end if

!  This is now the inverse reactance!
   tminv(2,2)=-aa2

!  This is now the inverse t-matrix!
   tminv(2,2)=tminv(2,2)+sqrtm1*p

   if(iflag.eq.0) return
!----------------------------------------------------------------
   call dirmagi1op(socsc,ce,l,my,vr,br,bopr,dx,xnot,ns, &
                   g12,f12,gp12)
!----------------------------------------------------------------
   if(scale) then
      gg2=g12(ns)/rs
      ff2=gp12/rs
      ggam2=ff2/gg2
      ggamb=(xl*fb(l)/rs-p*fb(l+1))/fb(l)
      aa1=p*((ggamb-xl/rs)*fn(l)+p*fn(l+1))/(ggamb*gg1-ff1)
      bb1=  ((ggam2-xl/rs)*fb(l)+p*fb(l+1))/(ggam2*gg1-ff1)
      bb2=  ((ggam1-xl/rs)*fb(l)+p*fb(l+1))/(ggam1*gg2-ff2)
      do i=1,ns
         gz(i,2,2)=g11(i)*aa1
         gj(i,2,2)=g11(i)*bb1+g12(i)*bb2
      end do
   else
      gg2=g12(ns)/rs
      ff2=f12(ns)/rs
      ggam2=ff2/gg2
      ggamb=smlm1*p*fb(l+1)/(c*fb(l))
      aa1=p*(ggamb*fn(l)-smlm1*p*fn(l+1)/c)/(ggamb*gg1-ff1)
      bb1=  (ggam2*fb(l)-smlm1*p*fb(l+1)/c)/(ggam2*gg1-ff1)
      bb2=  (ggam1*fb(l)-smlm1*p*fb(l+1)/c)/(ggam1*gg2-ff2)
      do i=1,ns
         gz(i,2,2)=g11(i)*aa1
         fz(i,2,2)=f11(i)*aa1
         gj(i,2,2)=g11(i)*bb1+g12(i)*bb2
         fj(i,2,2)=f11(i)*bb1+f12(i)*bb2
      end do
   end if
   
   return 

   1 continue
!---------------------------------------------------------------
   call dirmago2op(socsc,ce,l,my,vr,br,bopr,dx,xnot,ns, &
                   g1,f1,gp1)
!---------------------------------------------------------------
   jmat(1,1)=fb(l)
   jmat(2,2)=fb(l)
   nmat(1,1)=p*fn(l)
   nmat(2,2)=p*fn(l)

   if(scale) then
      jbmat(1,1)=xl*fb(l)/rs-p*fb(l+1)
      jbmat(2,2)=jbmat(1,1)
      nbmat(1,1)=xl*p*fn(l)/rs-p*p*fn(l+1)
      nbmat(2,2)=nbmat(1,1)
      do ii=1,2
         do jj=1,2
            g1mat(ii,jj)=g1(ii,jj,ns)/rs 
            f1mat(ii,jj)=gp1(ii,jj)/rs 
         end do
      end do
   else
      jbmat(1,1)=sl*p*fb(l-1)/c
      jbmat(2,2)=smlm1*p*fb(l+1)/c
      nbmat(1,1)=sl*p*p*fn(l-1)/c
      nbmat(2,2)=smlm1*p*p*fn(l+1)/c
      do ii=1,2
         do jj=1,2
            g1mat(ii,jj)=g1(ii,jj,ns)/rs 
            f1mat(ii,jj)=f1(ii,jj,ns)/rs 
         end do
      end do
   end if

   call repl(x1,g1mat,2,2)
   call gjinv(x1,2,2,detl)
   call repl(gam1,f1mat,2,2)
   call doubmt(gam1,x1,2,2)

   call repl(x1,gam1,2,2) 
   call doubmt(x1,jmat,2,2) 
   call submat(x1,jbmat,2,2)
   call repl(a2,x1,2,2) 
   call gjinv(a2,2,2,detl)
   call repl(x2,gam1,2,2) 
   call doubmt(x2,nmat,2,2) 
   call submat(x2,nbmat,2,2) 
   call doubmt(a2,x2,2,2)

! This is now the inverse reactance matrix!
   do ii=1,2                 
      do jj=1,2                 
         tminv(ii,jj)=-a2(ii,jj)
      end do
   end do

! This is now the inverse t-matrix!
   tminv(1,1)=tminv(1,1)+sqrtm1*p
   tminv(2,2)=tminv(2,2)+sqrtm1*p

   if(iflag.eq.0) return
!--------------------------------------------------------------
   call dirmagi2op(socsc,ce,l,my,vr,br,bopr,dx,xnot,ns, &
                   g2,f2,gp2)
!--------------------------------------------------------------
   if(scale) then
      do ii=1,2
         do jj=1,2
            g2mat(ii,jj)=g2(ii,jj,ns)/rs 
            f2mat(ii,jj)=gp2(ii,jj)/rs 
         end do
      end do
   else
      do ii=1,2
         do jj=1,2
            g2mat(ii,jj)=g2(ii,jj,ns)/rs 
            f2mat(ii,jj)=f2(ii,jj,ns)/rs 
         end do
      end do
   end if

   call repl(x2,g2mat,2,2)
   call gjinv(x2,2,2,detl)
   call repl(gam2,f2mat,2,2)
   call doubmt(gam2,x2,2,2)
   gamb(1,1)=jbmat(1,1)/jmat(1,1)
   gamb(2,2)=jbmat(2,2)/jmat(2,2)

   call repl(b2,gam1,2,2) 
   call doubmt(b2,g2mat,2,2) 
   call submat(b2,f2mat,2,2) 
   call gjinv(b2,2,2,detl)
   call doubmt(b2,x1,2,2) 

   call repl(a1,gamb,2,2)
   call doubmt(a1,g1mat,2,2)
   call submat(a1,f1mat,2,2)
   call gjinv(a1,2,2,detl)
   call repl(x1,gamb,2,2)
   call doubmt(x1,nmat,2,2)
   call submat(x1,nbmat,2,2) 
   call doubmt(a1,x1,2,2) 

   call repl(b1,gam2,2,2)
   call doubmt(b1,g1mat,2,2)
   call submat(b1,f1mat,2,2) 
   call gjinv(b1,2,2,detl)
   call repl(x1,gam2,2,2)
   call doubmt(x1,jmat,2,2)
   call submat(x1,jbmat,2,2) 
   call doubmt(b1,x1,2,2)

   do i=1,ns
      call repl(x1,g1(1:2,1:2,i),2,2)
      call repl(x2,f1(1:2,1:2,i),2,2)
      call doubmt(x1,a1,2,2)
      call doubmt(x2,a1,2,2)
      do ii=1,2
         do jj=1,2
            gz(i,ii,jj)=x1(ii,jj)
            fz(i,ii,jj)=x2(ii,jj)
         end do
      end do
      call repl(x1,g1(1:2,1:2,i),2,2)
      call repl(x2,f1(1:2,1:2,i),2,2)
      call repl(x3,g2(1:2,1:2,i),2,2)
      call repl(x4,f2(1:2,1:2,i),2,2)
      call doubmt(x1,b1,2,2)
      call doubmt(x2,b1,2,2)
      call doubmt(x3,b2,2,2)
      call doubmt(x4,b2,2,2)
      call addmat(x1,x3,2,2)
      call addmat(x2,x4,2,2)
      do ii=1,2
         do jj=1,2
            gj(i,ii,jj)=x1(ii,jj)
            fj(i,ii,jj)=x2(ii,jj)
         end do
      end do
   end do

   if(scale) then
      call zeroout(fz,8*iprpts)
      call zeroout(fj,8*iprpts)
   endif

   end subroutine spzwafu

!================================================================
   subroutine dirmago1op(socsc,e,l,my,vr,bspr,bopr, &
                         dx,xnot,ns,g,f,gp)
!================================================================
   implicit none

!  ***********************************************************
!  * integration of relativistic radial dirac equation       *
!  * in the presence of an internal field by adams method    *
!  * integrate outward!                                      *
!  * strange et al., j.phys.c.solid state phys.17,3355(1984) *
!  * scaling of SOC and OP-term ala Hubert Ebert included    *
!  *                                                         *
!  * my=+/-(l+1/2), kap=-l-1  case!!!                        *
!  ***********************************************************

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: p(iprpts), q(iprpts)
   complex (kind=CmplxKind) :: pp(iprpts), qp(iprpts)
   complex (kind=CmplxKind) :: pn,qn,k1,k2,k3,k4,m1,m2,m3,m4
   complex (kind=CmplxKind) :: psum,qsum,p0,q0,p1,q1,ppnp1,qpnp1
   complex (kind=CmplxKind) :: g(:), gp, f(:)

   integer (kind=IntKind), intent(in) :: l, my, ns
   integer (kind=IntKind), parameter :: imode=0, nitmax=100
   integer (kind=IntKind) :: n, nit

   real (kind=RealKind), intent(in) :: socsc, dx, xnot
   real (kind=RealKind) :: vr(:)
   real (kind=RealKind), intent(in) :: bspr(:), bopr(:,:)
   real (kind=RealKind) :: ba11(iprpts), ba12(iprpts)
   real (kind=RealKind) :: ba22(iprpts), bb11(iprpts), bb22(iprpts)
   real (kind=RealKind), parameter :: test=1.0d10, c=LightSpeed
   real (kind=RealKind), parameter :: dkoef1=475.d0/502.0
   real (kind=RealKind), parameter :: dkoef2=27.d0/502.0
   real (kind=RealKind) :: xl, xkap, xk, facl, dx2, cin, hoc, u
   real (kind=RealKind) :: r, rn, x, vrn, ba22n, bb22n, rnp1
   real (kind=RealKind) :: vrnp1, ba22np1, bb22np1, rmid
   real (kind=RealKind) :: vrmid, ba22mid, bb22mid
   real (kind=RealKind) :: vrc, ba22c, bb22c

   if(abs(my).lt.2*l+1) stop 'dirmag1: -l-1/2 < my < l+1/2'

   xl=dfloat(l)
   xkap=-xl-ONE
   xk=socsc*(ONE+xkap)-ONE
   facl=xl*(xl+ONE)-xk*(xk+ONE)

   dx2=HALF*dx
   cin = ONE/(c*c)

!  calculate boundary conditions

   if(vr(1).lt.-1.d-3)then
      hoc=-vr(1)/c
   else
      hoc=-vr(2)/c
      vr(1)=vr(2)
   endif
   u=(xk+sqrt(xk*xk-hoc*hoc+facl))/hoc
   p(1) = 1.0d-20
   q(1) = c*u*1.0d-20

!  get combined effective fields

   call brmat(l,my,ns,bspr,bopr,ba11,ba12,ba22,bb11,bb22)

   r=dexp(xnot)
!------------------------------------------------------------
   call dmag1op(pp(1),qp(1),p(1),q(1), &
                xk,facl,cin,e,r,vr(1),ba22(1),bb22(1))
!------------------------------------------------------------

!  start runge-kutta procedure (points 2, ... , 6)

   do n=1,5
      x=xnot+(n-1)*dx
      rn=dexp(x)
      vrn=vr(n)
      ba22n=ba22(n)
      bb22n=bb22(n)

      rnp1=dexp(x+dx)
      vrnp1=vr(n+1)
      ba22np1=ba22(n+1)
      bb22np1=bb22(n+1)

      rmid=dexp(x+dx2)
      vrmid=HALF*(vrn+vrnp1)
      ba22mid=HALF*(ba22n+ba22np1)
      bb22mid=HALF*(bb22n+ba22np1)

      pn=p(n)
      qn=q(n)
!-------------------------------------------------------------------
      call dmag1op(k1,m1,pn,qn, &
                   xk,facl,cin,e,rn,vrn,ba22n,bb22n)
!-------------------------------------------------------------------
      call dmag1op(k2,m2,pn+dx2*k1,qn+dx2*m1, &
                   xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)
!-------------------------------------------------------------------
      call dmag1op(k3,m3,pn+dx2*k2,qn+dx2*m2, &
                   xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)
!-------------------------------------------------------------------
      call dmag1op(k4,m4,pn+dx*k3,qn+dx*m3, &
                   xk,facl,cin,e,rnp1,vrnp1,ba22np1,bb22np1)
!-------------------------------------------------------------------
      p(n+1)=pn+dx*(k1+TWO*k2+TWO*k3+k4)/SIX
      q(n+1)=qn+dx*(m1+TWO*m2+TWO*m3+m4)/SIX
!-------------------------------------------------------------------
      call dmag1op(pp(n+1),qp(n+1),p(n+1),q(n+1), &
                   xk,facl,cin,e,rnp1,vrnp1,ba22np1,bb22np1)
!-------------------------------------------------------------------
   end do                 

!  begin adams procedure (points 7, 8, ... , ns)

   do n=6,ns-1
      x=xnot+(n-1)*dx
      r=dexp(x+dx)
      vrc=vr(n+1)
      ba22c=ba22(n+1)
      bb22c=bb22(n+1)

      psum=646.d0*pp(n)-264.d0*pp(n-1)+106.d0*pp(n-2)-19.d0*pp(n-3)
      qsum=646.d0*qp(n)-264.d0*qp(n-1)+106.d0*qp(n-2)-19.d0*qp(n-3)

!     predict for point n+1

      p0=p(n)+dx*(251.d0*pp(n-4)-1274.d0*pp(n-3)+2616.d0*pp(n-2)- &
                  2774.d0*pp(n-1)+1901.d0*pp(n))/720.d0
      q0=q(n)+dx*(251.d0*qp(n-4)-1274.d0*qp(n-3)+2616.d0*qp(n-2)- &
                  2774.d0*qp(n-1)+1901.d0*qp(n))/720.d0

!     correct

      if(imode.eq.1) then

         do nit=1,nitmax
!--------------------------------------------------------------------
            call dmag1op(ppnp1,qpnp1,p0,q0, &
                         xk,facl,cin,e,r,vrc,ba22c,bb22c)
!--------------------------------------------------------------------
            p1=p(n)+dx*(251.d0*ppnp1+psum)/720.d0
            q1=q(n)+dx*(251.d0*qpnp1+qsum)/720.d0

!           compare predictor with corrector

            if(test*abs(p1-p0).gt.abs(p0)) goto 15
            if(test*abs(q1-q0).gt.abs(q0)) goto 15 
            goto 10
     15     p0=p1   
            q0=q1   
         end do  
         write(6,*) n+1,r,nit,' not converged'
      else
!--------------------------------------------------------------------
         call dmag1op(ppnp1,qpnp1,p0,q0, &
                      xk,facl,cin,e,r,vrc,ba22c,bb22c)
!--------------------------------------------------------------------
         p1=p(n)+dx*(251.d0*ppnp1+psum)/720.d0
         q1=q(n)+dx*(251.d0*qpnp1+qsum)/720.d0
         q1=dkoef1*q1+dkoef2*q0
         p1=dkoef1*p1+dkoef2*p0
      end if

 10   q(n+1)=q1  
      p(n+1)=p1
!--------------------------------------------------------------------
      call dmag1op(pp(n+1),qp(n+1),p(n+1),q(n+1), &
                   xk,facl,cin,e,r,vrc,ba22c,bb22c)
!--------------------------------------------------------------------
   end do

!  store radial amplitudes times radius

   do n=1,ns
      g(n)=p(n)
      f(n)=q(n)/c
   end do
   gp=(pp(ns)-p(ns))/rs

   end subroutine dirmago1op

!=====================================================================
   subroutine dirmagi1op(socsc,e,l,my,vr,bspr,bopr, &
                         dx,xnot,ns,g,f,gp)
!=====================================================================

   use RelativityToolsModule, only : csbf

   implicit none

!  ***********************************************************
!  * integration of relativistic radial dirac equation       *
!  * in the presence of an internal field by adams method    *
!  * integrate inward!                                       *
!  * strange et al., j.phys.c.solid state phys.17,3355(1984) *
!  * scaling of SOC and OP-term ala Hubert Ebert included    *
!  *                                                         *
!  * my=+/-(l+1/2), kap=-l-1  case!!!                        *
!  ***********************************************************

   integer (kind=IntKind), intent(in) :: ns, l, my
   integer (kind=IntKind), parameter :: nitmax=100, imode=0, nout=5
   integer (kind=IntKind) :: lb, sk, ilag, i, nlag, n, nit, iex

   real (kind=RealKind) :: rs
   real (kind=RealKind), intent(in) :: socsc, dx, xnot
   real (kind=RealKind), intent(in) :: vr(:), bspr(:), bopr(:,:)
   real (kind=RealKind) :: xlag(iprpts), vrlag(iprpts)
   real (kind=RealKind) :: bsprlag(iprpts), boprlag(iprpts,2)
   real (kind=RealKind) :: ba11(iprpts), ba12(iprpts)
   real (kind=RealKind) :: ba22(iprpts), bb11(iprpts), bb22(iprpts)
   real (kind=RealKind), parameter :: test=1.0d10, c=LightSpeed
   real (kind=RealKind), parameter :: dkoef1=475.d0/502.d0
   real (kind=RealKind), parameter :: dkoef2= 27.d0/502.d0
   real (kind=RealKind) :: xl, xkap, xk, facl, dx2, cin, xs, x, r
   real (kind=RealKind) :: rn, vrn, ba22n, bb22n, rmid, vrmid
   real (kind=RealKind) :: ba22mid, bb22mid, rnm1, vrnm1
   real (kind=RealKind) :: ba22nm1, bb22nm1, rc, vrc
   real (kind=RealKind) :: ba22c, bb22c

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: psq, pe
   complex (kind=CmplxKind) :: p(iprpts), q(iprpts)
   complex (kind=CmplxKind) :: pp(iprpts), qp(iprpts)
   complex (kind=CmplxKind) :: pn,qn,k1,k2,k3,k4,m1,m2,m3,m4
   complex (kind=CmplxKind) :: psum,qsum,p0,q0,p1,q1,ppnm1,qpnm1
   complex (kind=CmplxKind), intent(out) :: g(:), gp, f(:)
   complex (kind=CmplxKind) :: fb(0:5), fn(0:5), fh(0:5)

   rs = Scatter(currentAtom)%r_end

   if(abs(my).lt.2*l+1) stop 'dirmag1: -l-1/2 < my < l+1/2'

   xl=dfloat(l)
   xkap=-xl-ONE
   xk=socsc*(ONE+xkap)-ONE
   facl=xl*(xl+ONE)-xk*(xk+ONE)
   lb=l+1
   sk=l-lb

   dx2=HALF*dx
   cin = ONE/(c*c)
   ilag=3
   xs=dlog(rs)
   do i=1,ns
      xlag(i)=xs+(i-ns)*dx
      vrlag(i)=vr(i)
      bsprlag(i)=bspr(i)
      boprlag(i,1)=bopr(i,1)
      boprlag(i,2)=bopr(i,2)
   end do
   nlag=ns+1
   vrlag(nlag)=ZERO
   bsprlag(nlag)=ZERO
   boprlag(nlag,1)=ZERO
   boprlag(nlag,2)=ZERO
   xlag(nlag)=xs+nout*dx

!  get combined effective fields

   call brmat(l,my,nlag,bsprlag,boprlag,ba11,ba12,ba22,bb11,bb22)

   psq=e+e*e*cin
   pe=cdsqrt(psq)

!  set up the starting values outside the muffin tin
!  corresponding to the boundary condition

   x=xs+nout*dx
   r=dexp(x)
   call csbf(l+1,pe,r,fb,fn,fh)
   n=ns+nout
   p(n)= r*fb(l)
   q(n)= ( (xk-xkap)*fb(l) + sk*pe*r*fb(lb) ) * e/psq
   call dmag1op(pp(n),qp(n),p(n),q(n), &
             xk,facl,cin,e,r,vrlag(nlag),ba22(nlag),bb22(nlag))

!  start runge-kutta procedure (points ns+nout-1, ... , ns+nout-5)

   do n=ns+nout,ns+nout-4,-1
      x=xs+(n-ns)*dx
      rn=dexp(x)
      vrn=ylag(x,xlag,vrlag,0,ilag,nlag,iex)
      ba22n=ylag(x,xlag,ba22,0,ilag,nlag,iex)
      bb22n=ylag(x,xlag,bb22,0,ilag,nlag,iex)

      rmid=dexp(x-dx2)
      vrmid=ylag(x-dx2,xlag,vrlag,0,ilag,nlag,iex)
      ba22mid=ylag(x-dx2,xlag,ba22,0,ilag,nlag,iex)
      bb22mid=ylag(x-dx2,xlag,bb22,0,ilag,nlag,iex)

      rnm1=dexp(x-dx)
      vrnm1=ylag(x-dx,xlag,vrlag,0,ilag,nlag,iex)
      ba22nm1=ylag(x-dx,xlag,ba22,0,ilag,nlag,iex)
      bb22nm1=ylag(x-dx,xlag,bb22,0,ilag,nlag,iex)

      pn=p(n)
      qn=q(n)
      call dmag1op(k1,m1,pn,qn, &
                xk,facl,cin,e,rn,vrn,ba22n,bb22n)

      call dmag1op(k2,m2,pn-dx2*k1,qn-dx2*m1, &
                xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)

      call dmag1op(k3,m3,pn-dx2*k2,qn-dx2*m2, &
                xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)

      call dmag1op(k4,m4,pn-dx*k3,qn-dx*m3, &
                xk,facl,cin,e,rnm1,vrnm1,ba22nm1,bb22nm1)

      p(n-1)=pn-dx*(k1+TWO*k2+TWO*k3+k4)/SIX
      q(n-1)=qn-dx*(m1+TWO*m2+TWO*m3+m4)/SIX
      call dmag1op(pp(n-1),qp(n-1),p(n-1),q(n-1), &
                 xk,facl,cin,e,rnm1,vrnm1,ba22nm1,bb22nm1)
   end do

!  begin adams procedure (points ns+nout-6,  ... , 1)

   do n=ns+nout-5,2,-1
      x=xs+(n-ns)*dx
      r=dexp(x-dx)
      vrc=vr(n-1)
      ba22c=ba22(n-1)
      bb22c=bb22(n-1)

      psum=646.d0*pp(n)-264.d0*pp(n+1)+106.d0*pp(n+2)-19.d0*pp(n+3)
      qsum=646.d0*qp(n)-264.d0*qp(n+1)+106.d0*qp(n+2)-19.d0*qp(n+3)

!     predict for point n-1

      p0=p(n)-dx*(251.d0*pp(n+4)-1274.d0*pp(n+3)+2616.d0*pp(n+2)- &
                  2774.d0*pp(n+1)+1901.d0*pp(n))/720.d0
      q0=q(n)-dx*(251.d0*qp(n+4)-1274.d0*qp(n+3)+2616.d0*qp(n+2)- &
                  2774.d0*qp(n+1)+1901.d0*qp(n))/720.d0

!     correct

      if(imode.eq.1) then
         do nit=1,nitmax
            call dmag1op(ppnm1,qpnm1,p0,q0, &
                       xk,facl,cin,e,r,vrc,ba22c,bb22c)
            p1=p(n)-dx*(251.d0*ppnm1+psum)/720.d0
            q1=q(n)-dx*(251.d0*qpnm1+qsum)/720.d0

!           compare predictor with corrector

            if(test*abs(p1-p0).gt.abs(p0)) goto 25
            if(test*abs(q1-q0).gt.abs(q0)) goto 25 
            goto 20
   25       p0=p1   
            q0=q1   
         end do  
         write(6,*) n-1,r,nit,' not converged'
      else
         call dmag1op(ppnm1,qpnm1,p0,q0, &
                    xk,facl,cin,e,r,vrc,ba22c,bb22c)
         p1=p(n)-dx*(251.d0*ppnm1+psum)/720.d0
         q1=q(n)-dx*(251.d0*qpnm1+qsum)/720.d0
         q1=dkoef1*q1+dkoef2*q0
         p1=dkoef1*p1+dkoef2*p0
      end if

 20   q(n-1)=q1  
      p(n-1)=p1
      call dmag1op(pp(n-1),qp(n-1),p(n-1),q(n-1), &
                 xk,facl,cin,e,r,vrc,ba22c,bb22c)
   end do

!  store radial amplitudes times radius

   do n=1,ns
      g(n)=p(n)
      f(n)=q(n)/c
   end do   
   gp=(pp(ns)-p(ns))/rs

   end subroutine dirmagi1op

!=================================================================
   subroutine dmag1op(pp,qp,p,q,xk,facl,cin,e,r,vr,ba22,bb22)
!=================================================================

   implicit none

   real (kind=RealKind), intent(in) :: xk, facl, cin, r
   real (kind=RealKind), intent(in) :: vr, ba22, bb22

   complex (kind=CmplxKind), intent(out) :: pp, qp
   complex (kind=CmplxKind), intent(in) :: p, q, e
   complex (kind=CmplxKind) :: tr, sr

   tr=r*e-vr
   sr=cin*tr+r+cin*bb22

   qp= xk*q - tr*p + ba22*p + facl*p/sr
   pp=-xk*p + sr*q

   end subroutine dmag1op

!================================================================
   subroutine dirmago2op(socsc,e,l,my,vr,bspr,bopr, &
                         dx,xnot,ns,g,f,gp)
!================================================================
   implicit none

!  ***********************************************************
!  * integration of relativistic radial dirac equation       *
!  * in the presence of an internal field by adams method    *
!  * integrate outward!                                      *
!  *                                                         *
!  * strange et al., j.phys.c.solid state phys.17,3355(1984) *
!  * scaling of SOC and OP-term ala Hubert Ebert included    *
!  ***********************************************************

   integer (kind=IntKind), parameter :: nitmax=100, imode=0
   integer (kind=IntKind), intent(in) :: l, my, ns
   integer (kind=IntKind) :: ii, n, nit

   real (kind=RealKind) :: rs, xl, dx2, cin, hoc, xkap, xk1, facl1
   real (kind=RealKind) :: u1a, xk2, facl2, u2a, r, rn, vrn, x
   real (kind=RealKind) :: ba11n, ba12n, ba22n, bb11n, bb22n
   real (kind=RealKind) :: rnp1, vrnp1, ba11np1, ba12np1
   real (kind=RealKind) :: ba22np1, bb11np1, bb22np1, rmid, vrmid
   real (kind=RealKind) :: ba11mid, ba12mid, ba22mid, bb11mid, bb22mid
   real (kind=RealKind) :: vrc, ba11c, ba12c, ba22c, bb11c, bb22c
   real (kind=RealKind), intent(in) :: socsc, dx, xnot
   real (kind=RealKind), intent(in) :: bspr(:), bopr(:,:)
   real (kind=RealKind) :: vr(:)
   real (kind=RealKind) :: ba11(iprpts), ba12(iprpts), ba22(iprpts)
   real (kind=RealKind) :: bb11(iprpts), bb22(iprpts)
   real (kind=RealKind) :: fk1(2), fk2(2), gk1(2), gk2(2)
   real (kind=RealKind), parameter :: test = 1.0d10, c=LightSpeed
   real (kind=RealKind), parameter :: dkoef1=475.d0/502.d0
   real (kind=RealKind), parameter :: dkoef2= 27.d0/502.d0

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: p1(iprpts), q1(iprpts)
   complex (kind=CmplxKind) :: p2(iprpts), q2(iprpts)
   complex (kind=CmplxKind) :: pp1(iprpts), qp1(iprpts)
   complex (kind=CmplxKind) :: pp2(iprpts), qp2(iprpts)
   complex (kind=CmplxKind) :: p1n,q1n,k11,k12,k13,k14,m11,m12,m13,m14       
   complex (kind=CmplxKind) :: p2n,q2n,k21,k22,k23,k24,m21,m22,m23,m24       
   complex (kind=CmplxKind) :: p1sum,q1sum,p10,q10,p11,q11,pp1np1,qp1np1          
   complex (kind=CmplxKind) :: p2sum,q2sum,p20,q20,p21,q21,pp2np1,qp2np1          
   complex (kind=CmplxKind), intent(out) :: g(:,:,:), gp(:,:), f(:,:,:)

   rs = Scatter(currentAtom)%r_end

   if(abs(my).eq.2*l+1) stop 'dirmag2: my=+/-(l+1/2)'
   xl = dfloat(l)
   dx2 = HALF*dx
   cin = ONE/(c*c)

!  calculate boundary conditions

   if(vr(1).lt.-1.d-3)then
      hoc=-vr(1)/c
   else
      hoc=-vr(2)/c
      vr(1)=vr(2)
   endif

   xl=dfloat(l)
   xkap=xl
   xk1=socsc*(ONE+xkap)-ONE
   facl1=xl*(xl+ONE)-xk1*(xk1+ONE)
   u1a=(xk1+sqrt(xk1*xk1-hoc*hoc+facl1))/hoc
   xkap=-xl-ONE
   xk2=socsc*(ONE+xkap)-ONE
   facl2=xl*(xl+ONE)-xk2*(xk2+ONE)
   u2a=(xk2+sqrt(xk2*xk2-hoc*hoc+facl2))/hoc

!  get combined effective fields

   call brmat(l,my,ns,bspr,bopr,ba11,ba12,ba22,bb11,bb22)

   do ii=1,2
      if(ii.eq.2) go to 40

      p1(1) = 1.0d-20
      q1(1) = c*u1a*1.0d-20
      p2(1) = CZERO
      q2(1) = CZERO
      go to 50

   40 continue
      p1(1) = CZERO
      q1(1) = CZERO
      p2(1) = 1.0d-20
      q2(1) = c*u2a*1.0d-20
   50 continue

   r=dexp(xnot)
   call dmag2op(pp1(1),qp1(1),pp2(1),qp2(1), &
              p1(1),q1(1),p2(1),q2(1),       &
              xk1,xk2,facl1,facl2,cin,e,r,   &
              vr(1),ba11(1),ba12(1),ba22(1),bb11(1),bb22(1))

! start runge-kutta procedure (points 2, ... , 6)

   do n=1,5
      x=xnot+(n-1)*dx
      rn=dexp(x)
      vrn=vr(n)
      ba11n=ba11(n)
      ba12n=ba12(n)
      ba22n=ba22(n)
      bb11n=bb11(n)
      bb22n=bb22(n)

      rnp1=dexp(x+dx)
      vrnp1=vr(n+1)
      ba11np1=ba11(n+1)
      ba12np1=ba12(n+1)
      ba22np1=ba22(n+1)
      bb11np1=bb11(n+1)
      bb22np1=bb22(n+1)

      rmid=dexp(x+dx2)
      vrmid=0.5d0*(vrn+vrnp1)
      ba11mid=0.5d0*(ba11n+ba11np1)
      ba12mid=0.5d0*(ba12n+ba12np1)
      ba22mid=0.5d0*(ba22n+ba22np1)
      bb11mid=0.5d0*(bb11n+bb11np1)
      bb22mid=0.5d0*(bb22n+bb22np1)

      p1n=p1(n)
      q1n=q1(n)
      p2n=p2(n)
      q2n=q2(n)
      call dmag2op(k11,m11,k21,m21,p1n,q1n,p2n,q2n, &
                 xk1,xk2,facl1,facl2,cin,e,rn,      &
                 vrn,ba11n,ba12n,ba22n,bb11n,bb22n) 

      call dmag2op(k12,m12,k22,m22,                               &
                 p1n+dx2*k11,q1n+dx2*m11,p2n+dx2*k21,q2n+dx2*m21, &
                 xk1,xk2,facl1,facl2,cin,e,rmid,                  &
                 vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)

      call dmag2op(k13,m13,k23,m23,                               &
                 p1n+dx2*k12,q1n+dx2*m12,p2n+dx2*k22,q2n+dx2*m22, &
                 xk1,xk2,facl1,facl2,cin,e,rmid,                  &
                 vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)

      call dmag2op(k14,m14,k24,m24,                           &
                 p1n+dx*k13,q1n+dx*m13,p2n+dx*k23,q2n+dx*m23, &
                 xk1,xk2,facl1,facl2,cin,e,rnp1,              &
                 vrnp1,ba11np1,ba12np1,ba22np1,bb11np1,bb22np1)

      p1(n+1)=p1n+dx*(k11+2.d0*k12+2.d0*k13+k14)/6.d0
      q1(n+1)=q1n+dx*(m11+2.d0*m12+2.d0*m13+m14)/6.d0
      p2(n+1)=p2n+dx*(k21+2.d0*k22+2.d0*k23+k24)/6.d0
      q2(n+1)=q2n+dx*(m21+2.d0*m22+2.d0*m23+m24)/6.d0
      call dmag2op(pp1(n+1),qp1(n+1),pp2(n+1),qp2(n+1), &
                 p1(n+1),q1(n+1),p2(n+1),q2(n+1),       &
                 xk1,xk2,facl1,facl2,cin,e,rnp1,        &
                 vrnp1,ba11np1,ba12np1,ba22np1,bb11np1,bb22np1)
   end do

!  begin adams procedure (points 7, 8, ... , ns)

   do n=6,ns-1
      x=xnot+(n-1)*dx
      r=dexp(x+dx)
      vrc=vr(n+1)
      ba11c=ba11(n+1)
      ba12c=ba12(n+1)
      ba22c=ba22(n+1)
      bb11c=bb11(n+1)
      bb22c=bb22(n+1)

      p1sum=646.d0*pp1(n)-264.d0*pp1(n-1)+106.d0*pp1(n-2)- &
            19.d0*pp1(n-3)
      q1sum=646.d0*qp1(n)-264.d0*qp1(n-1)+106.d0*qp1(n-2)- &
            19.d0*qp1(n-3)
      p2sum=646.d0*pp2(n)-264.d0*pp2(n-1)+106.d0*pp2(n-2)- &
            19.d0*pp2(n-3)
      q2sum=646.d0*qp2(n)-264.d0*qp2(n-1)+106.d0*qp2(n-2)- &
            19.d0*qp2(n-3)

!     predict for point n+1

      p10=p1(n)+dx*(251.d0*pp1(n-4)-1274.d0*pp1(n-3)+ &
          2616.d0*pp1(n-2)-2774.d0*pp1(n-1)+1901.d0*pp1(n))/720.d0
      q10=q1(n)+dx*(251.d0*qp1(n-4)-1274.d0*qp1(n-3)+ &
          2616.d0*qp1(n-2)-2774.d0*qp1(n-1)+1901.d0*qp1(n))/720.d0
      p20=p2(n)+dx*(251.d0*pp2(n-4)-1274.d0*pp2(n-3)+ &
          2616.d0*pp2(n-2)-2774.d0*pp2(n-1)+1901.d0*pp2(n))/720.d0
      q20=q2(n)+dx*(251.d0*qp2(n-4)-1274.d0*qp2(n-3)+ &
          2616.d0*qp2(n-2)-2774.d0*qp2(n-1)+1901.d0*qp2(n))/720.d0

!     correct

      if(imode.eq.1) then
         do nit=1,nitmax
            call dmag2op(pp1np1,qp1np1,pp2np1,qp2np1, &
                       p10,q10,p20,q20, &
                       xk1,xk2,facl1,facl2,cin,e,r, &
                       vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
            p11=p1(n)+dx*(251.d0*pp1np1+p1sum)/720.d0
            q11=q1(n)+dx*(251.d0*qp1np1+q1sum)/720.d0
            p21=p2(n)+dx*(251.d0*pp2np1+p2sum)/720.d0
            q21=q2(n)+dx*(251.d0*qp2np1+q2sum)/720.d0

!           compare predictor with corrector

            if(test*abs(p11-p10).gt.abs(p10)) goto 35
            if(test*abs(q11-q10).gt.abs(q10)) goto 35
            if(test*abs(p21-p20).gt.abs(p20)) goto 35
            if(test*abs(q21-q20).gt.abs(q20)) goto 35
            goto 30

   35       p10=p11
            q10=q11
            p20=p21
            q20=q21
         end do
         write(6,*) n+1,r,nit,' not converged'
      else
         call dmag2op(pp1np1,qp1np1,pp2np1,qp2np1, &
                    p10,q10,p20,q20, &
                    xk1,xk2,facl1,facl2,cin,e,r, &
                    vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
         p11=p1(n)+dx*(251.d0*pp1np1+p1sum)/720.d0
         q11=q1(n)+dx*(251.d0*qp1np1+q1sum)/720.d0
         p21=p2(n)+dx*(251.d0*pp2np1+p2sum)/720.d0
         q21=q2(n)+dx*(251.d0*qp2np1+q2sum)/720.d0

         q11=dkoef1*q11+dkoef2*q10
         p11=dkoef1*p11+dkoef2*p10
         q21=dkoef1*q21+dkoef2*q20
         p21=dkoef1*p21+dkoef2*p20

      end if

 30      q1(n+1)=q11
         p1(n+1)=p11
         q2(n+1)=q21
         p2(n+1)=p21
         call dmag2op(pp1(n+1),qp1(n+1),pp2(n+1),qp2(n+1), &
                    p1(n+1),q1(n+1),p2(n+1),q2(n+1),       &
                    xk1,xk2,facl1,facl2,cin,e,r,           &
                    vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
         end do

!     store radial amplitudes times radius

      do n=1,ns
         g(1,ii,n)=p1(n)
         g(2,ii,n)=p2(n)
         f(1,ii,n)=q1(n)/c
         f(2,ii,n)=q2(n)/c
      end do
      gp(1,ii)=(pp1(ns)-p1(ns))/rs
      gp(2,ii)=(pp2(ns)-p2(ns))/rs

   end do  

   end subroutine dirmago2op

!=================================================================
   subroutine dirmagi2op(socsc,e,l,my,vr,bspr,bopr, &
                         dx,xnot,ns,g,f,gp)
!=================================================================

   use RelativityToolsModule, only : csbf

   implicit none

!  ***********************************************************
!  * integration of relativistic radial dirac equation       *
!  * in the presence of an internal field by adams method    *
!  * integrate inward!                                       *
!  *                                                         *
!  * strange et al., j.phys.c.solid state phys.17,3355(1984) *
!  * scaling of SOC and OP-term ala Hubert Ebert included    *
!  ***********************************************************

   integer (kind=IntKind), intent(in) :: l, my, ns
   integer (kind=IntKind), parameter :: nitmax=100, imode=0, nout=5
   integer (kind=IntKind) :: lb1, lb2, ilag, i, ii, nlag, n, nit, iex

   real (kind=RealKind), intent(in) :: socsc, dx, xnot
   real (kind=RealKind), intent(in) :: bspr(:), bopr(:,:)
   real (kind=RealKind), parameter :: test=1.0d10, c=LightSpeed
   real (kind=RealKind), parameter :: dkoef1=475.d0/502.d0
   real (kind=RealKind), parameter :: dkoef2= 27.d0/502.d0
   real (kind=RealKind) :: vr(:), xlag(iprpts), vrlag(iprpts)
   real (kind=RealKind) :: bsprlag(iprpts), boprlag(iprpts,2)
   real (kind=RealKind) :: ba11(iprpts), ba12(iprpts)
   real (kind=RealKind) :: ba22(iprpts), bb11(iprpts), bb22(iprpts)
   real (kind=RealKind) :: fk1(2), fk2(2), gk1(2), gk2(2)
   real (kind=RealKind) :: rs, sk1, sk2, xl, xkap, xk1, facl1, xk2
   real (kind=RealKind) :: facl2, dx2, cin, xs, x, r, rn, vrn
   real (kind=RealKind) :: ba11n, ba12n, ba22n, bb11n, bb22n
   real (kind=RealKind) :: rmid, vrmid, ba11mid, ba12mid, ba22mid
   real (kind=RealKind) :: bb11mid, bb22mid, rnm1, vrnm1
   real (kind=RealKind) :: ba11nm1, ba12nm1, ba22nm1, bb11nm1, bb22nm1
   real (kind=RealKind) :: vrc, ba11c, ba12c, ba22c, bb11c, bb22c

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: psq, pe
   complex (kind=CmplxKind) :: p1(iprpts), q1(iprpts)
   complex (kind=CmplxKind) :: p2(iprpts), q2(iprpts)
   complex (kind=CmplxKind) :: pp1(iprpts), qp1(iprpts)
   complex (kind=CmplxKind) :: pp2(iprpts), qp2(iprpts)
   complex (kind=CmplxKind) :: p1n,q1n,k11,k12,k13,k14,m11,m12,m13,m14       
   complex (kind=CmplxKind) :: p2n,q2n,k21,k22,k23,k24,m21,m22,m23,m24       
   complex (kind=CmplxKind) :: p1sum,q1sum,p10,q10,p11,q11,pp1nm1,qp1nm1          
   complex (kind=CmplxKind) :: p2sum,q2sum,p20,q20,p21,q21,pp2nm1,qp2nm1          
   complex (kind=CmplxKind), intent(out) :: g(:,:,:), gp(:,:), f(:,:,:)
   complex (kind=CmplxKind) :: fb(0:5),fn(0:5),fh(0:5)

   rs = Scatter(currentAtom)%r_end

   if(abs(my).eq.2*l+1) stop 'dirmag2: my=+/-(l+1/2)'

   lb1=l-1
   lb2=l+1
   sk1=l-lb1
   sk2=l-lb2

   xl = dfloat(l)
   xkap=xl
   xk1=socsc*(1.d0+xkap)-1.d0
   facl1=xl*(xl+1.d0)-xk1*(xk1+1.d0)
   xkap=-xl-1.d0
   xk2=socsc*(1.d0+xkap)-1.d0
   facl2=xl*(xl+1.d0)-xk2*(xk2+1.d0)

   dx2=0.5d0*dx
   cin = 1.0d00/(c*c)
   ilag=3
   xs=dlog(rs)
   do i=1,ns
      xlag(i)=xs+(i-ns)*dx
      vrlag(i)=vr(i)
      bsprlag(i)=bspr(i)
      boprlag(i,1)=bopr(i,1)
      boprlag(i,2)=bopr(i,2)
   end do
   nlag=ns+1
   vrlag(nlag)=0.d0
   bsprlag(nlag)=0.d0
   boprlag(nlag,1)=0.d0
   boprlag(nlag,2)=0.d0
   xlag(nlag)=xs+nout*dx

!  get combined effective fields

   call brmat(l,my,nlag,bsprlag,boprlag,ba11,ba12,ba22,bb11,bb22)

   psq=e+e*e*cin
   pe=cdsqrt(psq)

   do ii=1,2
!     set up the starting values outside the muffin tin
!     corresponding to the boundary condition

      x=xs+nout*dx
      r=dexp(x)
      call csbf(l+1,pe,r,fb,fn,fh)
      n=ns+nout

      if(ii.eq.2) go to 60
      p1(n) = r*fb(l)
      q1(n) = ( (xk1-xl)*fb(l) + sk1*pe*r*fb(lb1) ) * e/psq
      p2(n) = 0.d0             
      q2(n) = 0.d0                     
      call dmag2op(pp1(n),qp1(n),pp2(n),qp2(n),                &
                 p1(n),q1(n),p2(n),q2(n),                      &
                 xk1,xk2,facl1,facl2,cin,e,r,                  &
                 vrlag(nlag),ba11(nlag),ba12(nlag),ba22(nlag), &
                             bb11(nlag),bb22(nlag))
      go to 70
60 continue
      p1(n) = 0.d0    
      q1(n) = 0.d0                    
      p2(n) = r*fb(l)          
      q2(n) = ( (xk2+xl+1.d0)*fb(l) + sk2*pe*r*fb(lb2) ) * e/psq 
      call dmag2op(pp1(n),qp1(n),pp2(n),qp2(n),                &
                 p1(n),q1(n),p2(n),q2(n),                      &
                 xk1,xk2,facl1,facl2,cin,e,r,                  &
                 vrlag(nlag),ba11(nlag),ba12(nlag),ba22(nlag), &
                             bb11(nlag),bb22(nlag))
70 continue

! start runge-kutta procedure (points ns+nout-1, ... , ns+nout-5)

   do n=ns+nout,ns+nout-4,-1
      x=xs+(n-ns)*dx
      rn=dexp(x)
      vrn=ylag(x,xlag,vrlag,0,ilag,nlag,iex)
      ba11n=ylag(x,xlag,ba11,0,ilag,nlag,iex)
      ba12n=ylag(x,xlag,ba12,0,ilag,nlag,iex)
      ba22n=ylag(x,xlag,ba22,0,ilag,nlag,iex)
      bb11n=ylag(x,xlag,bb11,0,ilag,nlag,iex)
      bb22n=ylag(x,xlag,bb22,0,ilag,nlag,iex)

      rmid=dexp(x-dx2)
      vrmid=ylag(x-dx2,xlag,vrlag,0,ilag,nlag,iex)
      ba11mid=ylag(x-dx2,xlag,ba11,0,ilag,nlag,iex)
      ba12mid=ylag(x-dx2,xlag,ba12,0,ilag,nlag,iex)
      ba22mid=ylag(x-dx2,xlag,ba22,0,ilag,nlag,iex)
      bb11mid=ylag(x-dx2,xlag,bb11,0,ilag,nlag,iex)
      bb22mid=ylag(x-dx2,xlag,bb22,0,ilag,nlag,iex)

      rnm1=dexp(x-dx)
      vrnm1=ylag(x-dx,xlag,vrlag,0,ilag,nlag,iex)
      ba11nm1=ylag(x-dx,xlag,ba11,0,ilag,nlag,iex)
      ba12nm1=ylag(x-dx,xlag,ba12,0,ilag,nlag,iex)
      ba22nm1=ylag(x-dx,xlag,ba22,0,ilag,nlag,iex)
      bb11nm1=ylag(x-dx,xlag,bb11,0,ilag,nlag,iex)
      bb22nm1=ylag(x-dx,xlag,bb22,0,ilag,nlag,iex)

      p1n=p1(n)
      q1n=q1(n)
      p2n=p2(n)
      q2n=q2(n)
      call dmag2op(k11,m11,k21,m21,p1n,q1n,p2n,q2n, &
                 xk1,xk2,facl1,facl2,cin,e,rn,      &
                 vrn,ba11n,ba12n,ba22n,bb11n,bb22n)

      call dmag2op(k12,m12,k22,m22,                               &
                 p1n-dx2*k11,q1n-dx2*m11,p2n-dx2*k21,q2n-dx2*m21, &
                 xk1,xk2,facl1,facl2,cin,e,rmid,                  &
                 vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)

      call dmag2op(k13,m13,k23,m23,                               &
                 p1n-dx2*k12,q1n-dx2*m12,p2n-dx2*k22,q2n-dx2*m22, &
                 xk1,xk2,facl1,facl2,cin,e,rmid,                  &
                 vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)

      call dmag2op(k14,m14,k24,m24,                           &
                 p1n-dx*k13,q1n-dx*m13,p2n-dx*k23,q2n-dx*m23, &
                 xk1,xk2,facl1,facl2,cin,e,rnm1,              &
                 vrnm1,ba11nm1,ba12nm1,ba22nm1,bb11nm1,bb22nm1)

      p1(n-1)=p1n-dx*(k11+2.d0*k12+2.d0*k13+k14)/6.d0
      q1(n-1)=q1n-dx*(m11+2.d0*m12+2.d0*m13+m14)/6.d0
      p2(n-1)=p2n-dx*(k21+2.d0*k22+2.d0*k23+k24)/6.d0
      q2(n-1)=q2n-dx*(m21+2.d0*m22+2.d0*m23+m24)/6.d0
      call dmag2op(pp1(n-1),qp1(n-1),pp2(n-1),qp2(n-1), &
                 p1(n-1),q1(n-1),p2(n-1),q2(n-1),       &
                 xk1,xk2,facl1,facl2,cin,e,rnm1,        &
                 vrnm1,ba11nm1,ba12nm1,ba22nm1,bb11nm1,bb22nm1)

   end do


!  begin adams procedure (points ns+nout-6,  ... , 1)

   do n=ns+nout-5,2,-1
      x=xs+(n-ns)*dx
      r=dexp(x-dx)
      vrc=vr(n-1)
      ba11c=ba11(n-1)
      ba12c=ba12(n-1)
      ba22c=ba22(n-1)
      bb11c=bb11(n-1)
      bb22c=bb22(n-1)

      p1sum=646.d0*pp1(n)-264.d0*pp1(n+1)+106.d0*pp1(n+2)- &
            19.d0*pp1(n+3)
      q1sum=646.d0*qp1(n)-264.d0*qp1(n+1)+106.d0*qp1(n+2)- &
            19.d0*qp1(n+3)
      p2sum=646.d0*pp2(n)-264.d0*pp2(n+1)+106.d0*pp2(n+2)- &
            19.d0*pp2(n+3)
      q2sum=646.d0*qp2(n)-264.d0*qp2(n+1)+106.d0*qp2(n+2)- &
            19.d0*qp2(n+3)

!     predict for point n-1

      p10=p1(n)-dx*(251.d0*pp1(n+4)-1274.d0*pp1(n+3)+ &
          2616.d0*pp1(n+2)-2774.d0*pp1(n+1)+1901.d0*pp1(n))/720.d0
      q10=q1(n)-dx*(251.d0*qp1(n+4)-1274.d0*qp1(n+3)+ &
          2616.d0*qp1(n+2)-2774.d0*qp1(n+1)+1901.d0*qp1(n))/720.d0
      p20=p2(n)-dx*(251.d0*pp2(n+4)-1274.d0*pp2(n+3)+ &
          2616.d0*pp2(n+2)-2774.d0*pp2(n+1)+1901.d0*pp2(n))/720.d0
      q20=q2(n)-dx*(251.d0*qp2(n+4)-1274.d0*qp2(n+3)+ &
          2616.d0*qp2(n+2)-2774.d0*qp2(n+1)+1901.d0*qp2(n))/720.d0

!     correct

      if(imode.eq.1) then

      do nit=1,nitmax
         call dmag2op(pp1nm1,qp1nm1,pp2nm1,qp2nm1, &
                    p10,q10,p20,q20,               &
                    xk1,xk2,facl1,facl2,cin,e,r,   &
                    vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
         p11=p1(n)-dx*(251.d0*pp1nm1+p1sum)/720.d0
         q11=q1(n)-dx*(251.d0*qp1nm1+q1sum)/720.d0
         p21=p2(n)-dx*(251.d0*pp2nm1+p2sum)/720.d0
         q21=q2(n)-dx*(251.d0*qp2nm1+q2sum)/720.d0

!        compare predictor with corrector

         if(test*abs(p11-p10).gt.abs(p10)) goto 85
         if(test*abs(q11-q10).gt.abs(q10)) goto 85
         if(test*abs(p21-p20).gt.abs(p20)) goto 85
         if(test*abs(q21-q20).gt.abs(q20)) goto 85
         goto 80

  85     p10=p11
         q10=q11
         p20=p21
         q20=q21

      end do
      write(6,*) n-1,r,nit,' not converged'

      else
         call dmag2op(pp1nm1,qp1nm1,pp2nm1,qp2nm1, &
                    p10,q10,p20,q20,               &
                    xk1,xk2,facl1,facl2,cin,e,r,   &
                    vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
         p11=p1(n)-dx*(251.d0*pp1nm1+p1sum)/720.d0
         q11=q1(n)-dx*(251.d0*qp1nm1+q1sum)/720.d0
         p21=p2(n)-dx*(251.d0*pp2nm1+p2sum)/720.d0
         q21=q2(n)-dx*(251.d0*qp2nm1+q2sum)/720.d0

         q11=dkoef1*q11+dkoef2*q10
         p11=dkoef1*p11+dkoef2*p10
         q21=dkoef1*q21+dkoef2*q20
         p21=dkoef1*p21+dkoef2*p20

      end if

    80   q1(n-1)=q11
         p1(n-1)=p11
         q2(n-1)=q21
         p2(n-1)=p21
         call dmag2op(pp1(n-1),qp1(n-1),pp2(n-1),qp2(n-1), &
                    p1(n-1),q1(n-1),p2(n-1),q2(n-1),       &
                    xk1,xk2,facl1,facl2,cin,e,r,           &
                    vrc,ba11c,ba12c,ba22c,bb11c,bb22c)

      end do

!  store radial amplitudes times radius

      do n=1,ns+nout
         g(1,ii,n)=p1(n)
         g(2,ii,n)=p2(n)
         f(1,ii,n)=q1(n)/c
         f(2,ii,n)=q2(n)/c
      end do
      gp(1,ii)=(pp1(ns)-p1(ns))/rs
      gp(2,ii)=(pp2(ns)-p2(ns))/rs

   end do  

   end subroutine dirmagi2op

!==================================================================
   subroutine dmag2op(pp1,qp1,pp2,qp2,p1,q1,p2,q2, &
                    xk1,xk2,facl1,facl2,cin,e,r,   &
                    vr,ba11,ba12,ba22,bb11,bb22)
!==================================================================

   implicit none

   real (kind=RealKind), intent(in) :: xk1, xk2, facl1, facl2, r
   real (kind=RealKind), intent(in) :: vr, ba11, ba12, ba22
   real (kind=RealKind), intent(in) :: bb11, bb22, cin

   complex (kind=CmplxKind), intent(out) :: pp1, qp1, pp2, qp2
   complex (kind=CmplxKind), intent(in) :: p1, q1, p2, q2, e
   complex (kind=CmplxKind) :: tr, sr1, sr2

   tr=r*e-vr
   sr1=cin*tr+r+cin*bb11     
   sr2=cin*tr+r+cin*bb22     

   qp1= xk1*q1 - tr*p1 + ba11*p1 + ba12*p2 + facl1*p1/sr1 
   pp1=-xk1*p1 + sr1*q1
   qp2= xk2*q2 - tr*p2 + ba22*p2 + ba12*p1 + facl2*p2/sr2 
   pp2=-xk2*p2 + sr2*q2

   end subroutine dmag2op

!=================================================================
   subroutine clebsch
!=================================================================
   implicit none

   integer (kind=IntKind) :: i, l, m, inr, kap, ir

   real (kind=RealKind), parameter :: tiny = 1.0d-15
   real (kind=RealKind) :: twolp1

!  clebsch-gordan rectangular matrices to transform from (lm) to
! (kappa,my) basis

   do 101 i=1,500
      u1(i)=ZERO
      u2(i)=ZERO
      ind1(i)=1
      ind2(i)=1
   101 continue

   inr=0
   do 103 l=0,12
      twolp1=dfloat(2*l+1)
      do m=-l,l
         inr=inr+1

! j=l-1/2
         kap=l
         if(kap.eq.0) goto 102

! ms=-1/2
         ir=2*kap*kap+kap+m
         u1(ir)=sqrt((l+m)/twolp1)
         ind1(ir)=inr

! ms=+1/2
         ir=2*kap*kap+kap+m+1
         u2(ir)=-sqrt((l-m)/twolp1)
         ind2(ir)=inr
    102  continue

! j=l+1/2
         kap=-l-1

! ms=-1/2
         ir=2*kap*kap+kap+m
         u1(ir)=sqrt((l-m+1)/twolp1)
         ind1(ir)=inr

! ms=+1/2
         ir=2*kap*kap+kap+m+1
         u2(ir)=sqrt((l+m+1)/twolp1)
         ind2(ir)=inr
      enddo
  103 continue

  do ir=1,inr
     if(abs(u1(ir)).lt.tiny) ind1(ir)=1
  end do
  do ir=1,inr
     if(abs(u2(ir)).lt.tiny) ind2(ir)=1
  end do

   end subroutine clebsch

!========================================================================
   function getCGCu1() result(cgc1)
!========================================================================

   implicit none

   real (kind=RealKind), pointer :: cgc1(:)

   if (.not.Initialized) then
      call ErrorHandler('getCGCu1', 'module is not initialized')
   endif

   cgc1 => u1(1:500)

   end function getCGCu1

!=========================================================================
   function getCGCu2() result(cgc2)
!=========================================================================

   implicit none

   real (kind=RealKind), pointer :: cgc2(:)

   if(.not.Initialized) then
      call ErrorHandler('getCGCu2', 'module is not initialized')
   endif

   cgc2 => u2(1:500)

   end function getCGCu2

!=========================================================================
   function getCGCind1() result(pind1)
!=========================================================================

   implicit none

   integer (kind=IntKind), pointer :: pind1(:)

   if(.not.Initialized) then
      call ErrorHandler('getCGCind1', 'module is not intialized')
   endif

   pind1 => ind1(1:500)

   end function getCGCind1

!=========================================================================
   function getCGCind2() result(pind2)
!=========================================================================

   implicit none

   integer (kind=IntKind), pointer :: pind2(:)

   if(.not.Initialized) then
      call ErrorHandler('getCGCind2', 'module is not intialized')
   endif

   pind2 => ind2(1:500)

   end function getCGCind2


!====================================================================
   subroutine brmat(l,my,n,bspr,bopr, &
                    ba11,ba12,ba22,bb11,bb22)
!====================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: l, my, n
   integer (kind=IntKind) :: kap1, kapb1, kap1my, kapb1my, i
   integer (kind=IntKind) :: kap2, kapb2, kap2my, kapb2my

   real (kind=RealKind), intent(in) :: bspr(:), bopr(:,:)
   real (kind=RealKind), intent(out) :: ba11(:), ba12(:), ba22(:)
   real (kind=RealKind), intent(out) :: bb11(:), bb22(:)
   real (kind=RealKind) :: br(iprpts,2)
   real (kind=RealKind) :: ca11_u, ca11_d, ca12_u, ca12_d, cb11_u
   real (kind=RealKind) :: ca22_u, ca22_d, cb22_u, cb22_d, cb11_d

   call zeroout(ba11,iprpts)
   call zeroout(ba12,iprpts)
   call zeroout(ba22,iprpts)
   call zeroout(bb11,iprpts)
   call zeroout(bb22,iprpts)

   do i=1,n
      br(i,1)=-bspr(i)
      br(i,2)=bspr(i)
   end do
   if(l.eq.2) then
      do i=1,n
         br(i,1)=br(i,1)+bopr(i,1)*(my+1.0d0)/2.0d0
         br(i,2)=br(i,2)+bopr(i,2)*(my-1.0d0)/2.0d0
      end do
   end if

   if(abs(my).eq.2*l+1) then
      kap2=-l-1
      kapb2=l+1
      kap2my=2*kap2*kap2+kap2+(my+1)/2
      kapb2my=2*kapb2*kapb2+kapb2+(my+1)/2
      ca22_d=u1(kap2my)*u1(kap2my)
      ca22_u=u2(kap2my)*u2(kap2my)
      cb22_d=u1(kapb2my)*u1(kapb2my)
      cb22_u=u2(kapb2my)*u2(kapb2my)
      do i=1,n
         ba22(i)=ca22_d*br(i,1)+ca22_u*br(i,2)
         bb22(i)=cb22_d*br(i,1)+cb22_u*br(i,2)
      end do
   else 
      kap1=l
      kap2=-l-1
      kapb1=-l 
      kapb2=l+1
      kap1my=2*kap1*kap1+kap1+(my+1)/2
      kapb1my=2*kapb1*kapb1+kapb1+(my+1)/2
      kap2my=2*kap2*kap2+kap2+(my+1)/2
      kapb2my=2*kapb2*kapb2+kapb2+(my+1)/2
      ca11_d=u1(kap1my)*u1(kap1my)
      ca11_u=u2(kap1my)*u2(kap1my)
      ca12_d=u1(kap1my)*u1(kap2my)
      ca12_u=u2(kap1my)*u2(kap2my)
      ca22_d=u1(kap2my)*u1(kap2my)
      ca22_u=u2(kap2my)*u2(kap2my)
      cb11_d=u1(kapb1my)*u1(kapb1my)
      cb11_u=u2(kapb1my)*u2(kapb1my)
      cb22_d=u1(kapb2my)*u1(kapb2my)
      cb22_u=u2(kapb2my)*u2(kapb2my)
      do i=1,n
         ba11(i)=ca11_d*br(i,1)+ca11_u*br(i,2)
         ba12(i)=ca12_d*br(i,1)+ca12_u*br(i,2)
         ba22(i)=ca22_d*br(i,1)+ca22_u*br(i,2)
         bb11(i)=cb11_d*br(i,1)+cb11_u*br(i,2)
         bb22(i)=cb22_d*br(i,1)+cb22_u*br(i,2)
      end do
   endif

   end subroutine brmat

!==============================================================
   subroutine replms(alms, akmy, lmax, kmymax)
!==============================================================

! find diagonal elements of a matrix in (l,m,s) representation
! to be given in (kappa,my) representation

   implicit none

   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(inout) :: kmymax
   integer (kind=IntKind) :: l, lm, lmmax, m, lmm, lmp
   integer (kind=IntKind) :: kmy, kmyp

   complex (kind=CmplxKind), intent(in) :: akmy(:,:)
   complex (kind=CmplxKind), intent(out) :: alms(:)

   lmmax=(lmax+1)*(lmax+1)
   kmymax=2*lmmax

   lm=0
   do l=0,lmax
      do m=-l,l
         lm=lm+1
         lmm=lm
         lmp=lm+lmmax
         alms(lmm)=(0.d0,0.d0)
         alms(lmp)=(0.d0,0.d0)
         do kmy=1,kmymax
            do kmyp=1,kmymax
               if(ind1(kmy).eq.lm.and.ind1(kmyp).eq.lm) &
                  alms(lmm)=alms(lmm)+u1(kmy)*akmy(kmy,kmyp)*u1(kmyp)
               if(ind2(kmy).eq.lm.and.ind2(kmyp).eq.lm) &
                  alms(lmp)=alms(lmp)+u2(kmy)*akmy(kmy,kmyp)*u2(kmyp)
            end do
         end do
      end do
   end do

   end subroutine replms

!==============================================================
   subroutine matrot1(r0,r1,lmax,dmat,dmatp)
!==============================================================

! if r1 = R * r0, dmat = D(R) and dmatp = D(R)+ 

   use RelativityToolsModule, only : matr, rotmat
   implicit none

   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: i

   real (kind=RealKind), intent(inout) :: r0(3), r1(3)
   real (kind=RealKind) :: rr(3), drot(3,3), tvec(3), tlm(4)
   real (kind=RealKind) :: sth, r0mod, r1mod, cosphi, phi
   real (kind=RealKind) :: r0xy, tmod, cosp2, sinp2, trd

   complex (kind=CmplxKind), intent(out) :: dmat(:,:)
   complex (kind=CmplxKind), intent(out) :: dmatp(:,:)
   complex (kind=CmplxKind) :: cmat(2*(lmax+1)*(lmax+1),2*(lmax+1)*(lmax+1))
   complex (kind=CmplxKind) :: d(2*lmax+2,2*lmax+2), dz(3)

   sth=sqrt(3.d0)

   r0mod=0.d0
   r1mod=0.d0
   do i=1,3
      r0mod=r0mod+r0(i)*r0(i)
      r1mod=r1mod+r1(i)*r1(i)
   end do
   r0mod=sqrt(r0mod)
   r1mod=sqrt(r1mod)
   do i=1,3
      r0(i)=r0(i)/r0mod
      r1(i)=r1(i)/r1mod
   end do

!  Find normal vector and angle of rotation

   cosphi=r0(1)*r1(1)+r0(2)*r1(2)+r0(3)*r1(3)
   if(abs(cosphi-1.d0).lt.1.d-10) then
      phi=0.d0
      tvec(1)=r0(1) 
      tvec(2)=r0(2)
      tvec(3)=r0(3) 
      goto 10
   else if(abs(cosphi+1.d0).lt.1.d-10) then
      phi=pi
      r0xy=sqrt(r0(1)**2+r0(2)**2)
      if(r0xy.gt.1.d-10) then
         tvec(1)=-r0(2)/r0xy
         tvec(2)=r0(1)/r0xy
         tvec(3)=0.d0
      else
         tvec(1)=1.d0
         tvec(2)=0.d0
         tvec(3)=0.d0
      end if
      goto 10
   end if
   phi=dacos(cosphi)

!  T=(R0xR1)/|R0xR1|

   tvec(1)=r0(2)*r1(3)-r0(3)*r1(2)
   tvec(2)=r0(3)*r1(1)-r0(1)*r1(3)
   tvec(3)=r0(1)*r1(2)-r0(2)*r1(1)
   tmod=sqrt(tvec(1)**2+tvec(2)**2+tvec(3)**2)
   tvec(1)=tvec(1)/tmod
   tvec(2)=tvec(2)/tmod
   tvec(3)=tvec(3)/tmod

   10 continue

   if(0.gt.3) then
      write(6,'(''   r0= '',3f10.4)') r0
      write(6,'(''   r1= '',3f10.4)') r1
      write(6,'('' Normal vector: '',3f10.4)') tvec
      write(6,'('' Angle in degree:'',f10.4)') phi*180.d0/pi
   end if

   cosp2=dcos(phi/2.d0)
   sinp2=dsin(phi/2.d0)
   tlm(1)=cosp2
   tlm(2)=sinp2*tvec(1)
   tlm(3)=sinp2*tvec(2)
   tlm(4)=sinp2*tvec(3)

   call matr(tlm,lmax,dmat,dmatp)

!  Set up matrix of rotation for spherical harmonics with l=1 and m=0

   call rotmat(d,trd,3,tlm,1,lmax)
   dz(1)=dconjg(d(2,3))*sth
   dz(2)=dconjg(d(2,2))*sth
   dz(3)=dconjg(d(2,1))*sth

   if(0.gt.3) then
   write(6,*) ' Rotation matrix for (l,m)=(1,0)'
   write(6,'(6d13.5)') dz(1),dz(2),dz(3)
   end if

   end subroutine matrot1

!======================================================================
   function getRelTMatrix(atom, e) result(t_mat)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: t_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelTMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelTMat', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if(abs(e-last_energy) > TEN2m12) then
      call resetRelScatter()
   endif

   if (.not.Scatter(atom)%done) then
      call solveSingleDirac(atom, e)
   endif

   t_mat => Scatter(atom)%t_mat

   end function getRelTMatrix

!======================================================================
   function getGReg(atom, e) result(gz)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: gz(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getGReg', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getGReg', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if(abs(e-last_energy) > TEN2m12) then
      call resetRelScatter()
   endif

   if (.not.Scatter(atom)%done) then
      call solveSingleDirac(atom, e)
   endif

   gz => Scatter(atom)%gz

   end function getGReg

!======================================================================
   function getFReg(atom, e) result(fz)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: fz(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getFReg', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getFReg', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if(abs(e-last_energy) > TEN2m12) then
      call resetRelScatter()
   endif

   if (.not.Scatter(atom)%done) then
      call solveSingleDirac(atom, e)
   endif

   fz => Scatter(atom)%fz

   end function getFReg

!======================================================================
   function getGIrr(atom, e) result(gj)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: gj(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getGIrr', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getGIrr', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if(abs(e-last_energy) > TEN2m12) then
      call resetRelScatter()
   endif

   if (.not.Scatter(atom)%done) then
      call solveSingleDirac(atom, e)
   endif

   gj => Scatter(atom)%gj

   end function getGIrr

!======================================================================
   function getFIrr(atom, e) result(fj)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: fj(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getFIrr', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getFIrr', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if(abs(e-last_energy) > TEN2m12) then
      call resetRelScatter()
   endif

   if (.not.Scatter(atom)%done) then
      call solveSingleDirac(atom, e)
   endif

   fj => Scatter(atom)%fj

   end function getFIrr

!======================================================================
   function getnuz(atom, e) result(nuz)
!======================================================================

  implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), intent(in) :: e
   integer (kind=IntKind), pointer :: nuz(:)

   if (.not.Initialized) then
      call ErrorHandler('getnuz', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getnuz', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if(abs(e-last_energy) > TEN2m12) then
      call resetRelScatter()
   endif

   if (.not.Scatter(atom)%done) then
      call solveSingleDirac(atom, e)
   endif

   nuz => Scatter(atom)%nuz

   end function getnuz

!======================================================================
   function getindz(atom, e) result(indz)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), intent(in) :: e
   integer (kind=IntKind), pointer :: indz(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getindz', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getindz', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if(abs(e-last_energy) > TEN2m12) then
      call resetRelScatter()
   endif

   if (.not.Scatter(atom)%done) then
      call solveSingleDirac(atom, e)
   endif

   indz => Scatter(atom)%indz

   end function getindz

!======================================================================
   function getSizeOfRMesh(atom) result(nr)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom
   integer (kind=IntKind) :: nr

   if (.not.Initialized) then
      call ErrorHandler('getSizeOfRmesh', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getSizeOfRmesh', 'invalid number of local atoms', LocalNumAtoms)
   endif

   nr = Scatter(atom)%NumRs

   end function getSizeOfRMesh

!======================================================================
   function getDmat(atom) result(dmat)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: dmat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getDmat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getDMat', 'invalid number of local atoms', LocalNumAtoms)
   endif

   dmat => Scatter(atom)%d_mat

   end function getDMat

!======================================================================
   function getDmatP(atom) result(dmatp)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: dmatp(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getDmatP', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getDMatP', 'invalid number of local atoms', LocalNumAtoms)
   endif

   dmatp => Scatter(atom)%d_matp

   end function getDMatP

!======================================================================
   subroutine resetRelScatter()
!======================================================================

   implicit none

   integer (kind=IntKind) :: ia

   do ia=1, LocalNumAtoms
      Scatter(ia)%done= .false.
   enddo

   end subroutine resetRelScatter

!===========================================================================
   function ylag(xi,x,y,index,n1,imax,iex) result(lag)
!===========================================================================
!     lagrangian interpolation
!   modified from subroutine polint from Numerical Recipe, by Press et al.
!     xi is intepolated entry into x-array
!     n is the order of lagrangran interpolation
!     y is array from which ylag is obtained by interpolation
!     ind is the min-i for x(i).gt.xi
!     if ind=0,x-array will be searched
!     imax is max index of x-and y-arrays
!     extrapolation can occur,iex=-1 or +1

   implicit none

   integer (kind=IntKind), parameter :: nmax=10
   integer (kind=IntKind), intent(in) :: index, n1, imax
   integer (kind=IntKind), intent(out) :: iex
   integer (kind=IntKind) :: ind, n, j, j1, inl, inu, ns, i, ij

   real (kind=RealKind), intent(in) :: xi, x(:), y(:)
   real (kind=RealKind) :: c(nmax), d(nmax), dx(nmax)
   real (kind=RealKind) :: dift, dif, den, lag

   ind=index
   n=n1
   iex=0
   if (n.le.imax) go to 10
   n=imax
   iex=n
   10 if (ind.gt.0) go to 40
   do 20 j = 1,imax
!     if (xi-x(j)) 30, 130, 20
      if (xi-x(j) < 0.0d0) then
        goto 30
      else if (xi-x(j) > 0.0d0) then
        goto 20
      else
        goto 130
      endif
   20 continue
      iex=1
      go to 70
   30 ind=j
   40 if (ind.le.1) iex=-1
      inl=max(0,ind-ishft(n+1,-1))
      inu=inl+n
      if (inu.le.imax) go to 80
   70 inl=imax-n
      inu=imax
   80 dif=abs(xi-x(inl+1))
      ns=1
      do 105 j=1,inu-inl
	 j1=j+inl
	 dx(j)=xi-x(j1)
         dift=abs(dx(j))
         if(dift.lt.dif) then
           ns=j
           dif=dift
         endif
         c(j)=y(j1)
         d(j)=y(j1)
  105 continue
      lag=y(inl+ns)
      ns=ns-1
      do 110 j=inl+1,inu-1
         do 100 i=1,inu-j
	    ij=i+j
            den=x(i+inl)-x(ij)
            if(den.eq.0.d0) stop 'Two xs are the same in ylag.'
            den=(d(i)-c(i+1))/den
            d(i)=dx(ij-inl)*den
            c(i)=dx(i)*den
  100    continue
         if(ishft(ns,1).lt.inu-j) then
           lag=lag+c(2*ns+1)
         else
           lag=lag+d(ns)
           ns=ns-1
         endif
  110 continue
      return
  130 lag=y(j)
   end function ylag

end module DiracSolverModule
