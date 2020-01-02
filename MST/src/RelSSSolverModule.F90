!**************************************************************************
! Name: RelSSSolverModule                                                 *
!                                                                         *
! Version 1.0: June 2016 by Xianglin Liu                                  *
!                                                                         *
! Description:                                                            *
!    Solves the full-potential single-site Dirac equation using sine and  *
!    cosine matrices formalism. Adapoted for LSMS library.                *
!    Relativistic version of SSSolverModule                               *
!                                                                         *
!-------------------------------------------------------------------------*
! subroutine initRelSSSolver(na, lmax, lmax_pot, z, isBxyz, iprint,&      *
!            stop_routine)                                                *
!    Initializes variables and allocates the memory needed for the        *
!    routines in this module. Must be called first.                       *
!    input: na - integer number of atoms on this processor                *
!           lmax - 1D integer array, lamx for square sine matrices        *
!           lmax_pot - 1D integer array, lmax of potential expansion      *
!           iprint - integer, print output level                          *
!           stop_routine - character string, name of routine in which to  * 
!                          end the program                                *
! ----------------------------------------------------------------------- *
! subroutine endRelSSSolver()                                             *
!    Deallocates the memory used in this module. This is the last routine *
!    that can be called unless initDiracSolver is called again.           *
! ----------------------------------------------------------------------- *  
! subroutine SingleDiracScattering(ia,Ek)                                 *
!    Solve the Dirac equation correspond to atom ia, with relativistic    *
!    kinetic energy Ek. First it allocates and assign values to gloabl    *
!    variabls, such as r(nr), vr_new(kl,nr), kyk(num_vl,2,lamax,lamax);   *
!    Then it calls SolveDirac() to solve the Dirac equation.              *
!    Note: kyk correspond to <chi_\Lambda|Y_Lp|chi_\Lambda'>, it finds    *
!    potential component coupled to given indices \Lambda and \Lambda'.   *
!    num_vl is number of non-zero potential component, with index in      *
!    list_vl the "2" corresponds to upper and lower part.                 *
! ----------------------------------------------------------------------- *
! subroutine SolveDirac(ia,Ek,lmax,lmax_p2,lamax,num_vl)                  *
!    Solve the Dirac equation from r=0                                    *
!    boundary condtion: sine(ir=1)=0, cosine(ir=1)=-1                     *
!    Logrithm grid; RK4 for the first few points, then 4th order predictor*
!    corrector for the rest points.                                       *
!    getffgg are called, with d/dx sine(r)=ff; d/dx cosine(r)=gg          *
!    output: sine matrix, cosine matrix, tau matrix, etc.                 *
! ----------------------------------------------------------------------- *  
! subroutine cal_tmat(ia,kappa,lamax,s_matrix,c_matrix,s_matrix_d,&       *
!            c_matrix_d,t_matrix,jinv_matrix,OmegaHat,OmegaHatInv)        *
!    Calculate all kinds of matrices needed. (not limited to t_matrix !)  *
!    t_mat: -1/k*S*(C-i*S)^{-1}                                           *
!    S_mat: 1-2*i*k*t_mat: (1 here is unit matrix)                        *
!    jost_mat: i*S-C                                                      *
!    jinv_mat: (i*S-C)^{-1}                                               *
!    Omega_mat: (S^t*S+C^t*C) to be defined when anti-pathology used      *
!    OmegaHat_mat: {S^T(I*S-C)}^{-1}                                      *
!    OmegaHatinv_mat: OmegaHatInv=S^T(I*S-C)                              *
! ----------------------------------------------------------------------- *                     
! subroutine computeRelSingleSiteDOS(id)                                  *
!     return the single site DOS, which is calculated in                   *
!     SingleDiracScattering by calling:                                    *
!     subroutine computeRelGF(ia): calculate the Green's function,         *
!     expanded on spherical harmonics and stored as green(ir,kl1)          *
!     subroutine computeRelDOS(atom,add_highl_fec): calculate the          *
!     r-dependent dos and stored as dos(ir,jl)                             *
! computeRelSingleSiteDOS use getVolumeintegration to integrate on space  *
!**************************************************************************
module RelSSSolverModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use PhysParamModule, only : LightSpeed
   use BesselModule,only: SphericalBessel, SphericalNeumann
   use IntegerFactorsModule, only : lofk, mofk, &
       jofk, lofj, mofj, kofj, m1m
   use IntegrationModule, only : calIntegration

   use PublicTypeDefinitionsModule, only : GridStruct


!  use TimerModule, only : getTime

!==============================================================================
!charge density use
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
   use GauntFactorsModule, only : getGauntFactor
!==============================================================================
!
   use MathParamModule, only : TEN2m12, TEN2m10, TEN2m8, TEN2m6, TEN2m4, TEN2m11

   use MathParamModule, only : PI, PI2, SQRT_PI, HALF, TEN2m5, TEN2m16
   use MathParamModule, only : CZERO, CONE, SQRTm1, Y0
!
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
!
!   use IntegrationModule, only : calIntegration
!
public :: initRelSSSolver, &
          endRelSSSolver, &
          SolveDirac, &
          SingleDiracScattering,    &
!          RelSSScattering, &
          getRelSineMatrix,         &
          getRelCosineMatrix,       &
          getRelJostMatrix,         &
          getRelJostInvMatrix,      &
!          getRelOmegaMatrix,        &
          getRelOmegaHatMatrix,     &
          getRelOmegaHatInvMatrix,  &
          getRelTMatrix,            &
          getRelSMatrix,            &
          getRelScatteringMatrix,   &
          getRelSolutionRmeshSize,  &
          getRelRegSolutionSine,       &
          getRelRegSolutionCosine,  &
          computeRelSingleSiteDOS,  &
          computeRelSingleSiteMMDOS,&
          calNum_kyk,                  &
          getRelkyk,                  &
          getRelkyk_Bz,                  &
          getRelkyk_Bx,                  &
          getRelkyk_By,               &
          getldex,                     &
          getkapdex,                  &
          getmydex,                  &
          getldex_bar,               &
          getRelSSDOS,             &
                        getRelSSGreen,            &
          getRelSSGreenOut
         

private
   logical :: isRealSpace = .false. !MST variable to choose betweeen LSMS and KKR
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind), allocatable :: AtomicNumber(:)
   integer (kind=IntKind) :: lmax_max, lmax_pot_max

   real (kind=RealKind), parameter :: Me = 0.5d0

   type RelScatterStruct
      logical :: done
      integer (kind=IntKind) :: l_max
      integer (kind=IntKind) :: lmax_pot
      integer (kind=IntKind) :: lamax !lamax=2*(lmax+1)**2
      integer (kind=IntKind) :: lmax_green
      integer (kind=IntKind) :: kmax_green
      integer (kind=IntKind) :: num_vl
      integer (kind=IntKind) :: numrs ! number of r points, not necessarily the same as numrs_cs
      integer (kind=IntKind) :: numrs_cs !number of circumscribed sphere points
      integer (kind=IntKind) :: numrs_mt
      integer (kind=IntKind) :: numrs_trunc
      real (kind=RealKind) :: dx
      real (kind=RealKind) :: evec(3)
      complex (kind=CmplxKind), pointer :: sin_mat(:,:)
      complex (kind=CmplxKind), pointer :: cos_mat(:,:)
      complex (kind=CmplxKind), pointer :: t_mat(:,:)
      complex (kind=CmplxKind), pointer :: S_mat(:,:)
      complex (kind=CmplxKind), pointer :: jost_mat(:,:)
      complex (kind=CmplxKind), pointer :: jinv_mat(:,:)
!      complex (kind=CmplxKind), pointer :: Omega_mat(:,:)
      complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)
      complex (kind=CmplxKind), pointer :: OmegaHatInv_mat(:,:)
      complex (kind=CmplxKind), pointer :: ss_mat_r(:,:,:)
      complex (kind=CmplxKind), pointer :: cc_mat_r(:,:,:)
      complex (kind=CmplxKind), pointer :: green(:,:)
      complex (kind=CmplxKind), pointer :: green_Bxyz(:,:,:)
      complex (kind=CmplxKind), pointer :: dos(:,:)
      complex (kind=CmplxKind), pointer :: dos_Bxyz(:,:,:)
      real (kind=RealKind), pointer :: kyk(:,:,:,:)
      real (kind=RealKind), pointer :: kyk_Bz(:,:,:,:)
      real (kind=RealKind), pointer :: kyk_Bx(:,:,:,:)
      real (kind=RealKind), pointer :: kyk_By(:,:,:,:)
!      complex (kind=CmplxKind), pointer :: pdos(:,:,:)
   end type RelScatterStruct

   type (RelScatterStruct), allocatable :: Scatter(:)
   type(GridStruct), pointer :: Grid
   !R Grid is obtained from data import, need to be revised after test finished

   logical :: Initialized = .false.
   logical, allocatable :: is_Bxyz(:)
   logical :: is_Bxyz_max

   character (len=32) :: istop
   
   real (kind=RealKind), allocatable :: u1(:),u2(:)  !Clebsch-Gordan coeffs:
   integer (kind=IntKind), allocatable, target:: kapdex(:), ldex(:), jdex(:), mydex(:) !set indices
   integer (kind=IntKind), allocatable, target:: ldex_bar(:)

   integer (kind=IntKind), allocatable :: print_instruction(:)

   complex (kind=CmplxKind), allocatable, target :: wks_sinmat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_cosmat(:)
!
   complex (kind=CmplxKind), allocatable, target :: wks_tmat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_Smat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_jostmat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_jinvmat(:)
!   complex (kind=CmplxKind), allocatable, target :: wks_Omega(:)
   complex (kind=CmplxKind), allocatable, target :: wks_OmegaHat(:)
   complex (kind=CmplxKind), allocatable, target :: wks_OmegaHatInv(:)
   complex (kind=CmplxKind), allocatable, target :: wks_ssmat_r(:)
   complex (kind=CmplxKind), allocatable, target :: wks_ccmat_r(:)
   complex (kind=CmplxKind), allocatable, target :: wks_green(:)
   complex (kind=CmplxKind), allocatable, target :: wks_green_Bxyz(:)
   complex (kind=CmplxKind), allocatable, target :: wks_dos(:)
   complex (kind=CmplxKind), allocatable, target :: wks_dos_Bxyz(:)
!   complex (kind=CmplxKind), allocatable, target :: wks_pdos(:)
   complex (kind=CmplxKind), allocatable, target :: TmpSpace(:)
   real (kind=RealKind), allocatable, target :: wks_kyk(:)
   real (kind=RealKind), allocatable, target :: wks_kyk_Bx(:)
   real (kind=RealKind), allocatable, target :: wks_kyk_By(:)
   real (kind=RealKind), allocatable, target :: wks_kyk_Bz(:)
!====================================================================
!  Gaunt Factors
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
   real (kind=RealKind), pointer :: cgnt(:,:,:)
!====================================================================
! input Grid variables !be careful about the global variables
   character*80 title
   integer (kind=IntKind) :: nr
   real (kind=RealKind) :: dx, zz
   real (kind=RealKind),allocatable:: r(:)

   complex(kind=CmplxKind), allocatable:: vr_s1(:,:),vr_s2(:,:) 
   complex(kind=CmplxKind), allocatable:: vr_new(:,:),br_new(:,:,:) 
!====================================================================
!  kyk variables
   real (kind=RealKind),allocatable, target :: kyk_wks(:,:,:,:)
   real (kind=RealKind), allocatable, target :: kyk_short_wks(:,:,:,:)
   real (kind=RealKind),allocatable, target :: kyk_Bxyz_wks(:,:,:,:,:)
   real (kind=RealKind), allocatable, target :: kyk_short_Bxyz_wks(:,:,:,:,:)

   real (kind=RealKind), pointer :: kyk(:,:,:,:)
   real (kind=RealKind), pointer :: kyk_short(:,:,:,:)
   real (kind=RealKind), pointer :: kyk_Bxyz(:,:,:,:,:)
   real (kind=RealKind), pointer :: kyk_short_Bxyz(:,:,:,:,:)
!====================================================================
!  Bessel functions
   complex (kind=CmplxKind), allocatable :: bjl(:,:), bnl(:,:)
!====================================================================
!  global total energy and momentum
   complex (kind=CmplxKind) :: kappa_g, energy_g
!====================================================================
   integer (kind=IntKind) :: iend, iend_max
!====================================================================
   !single-site solver accelerated by finding out non-zero matrix elements
   logical :: Speedup = .true.
   integer (kind=IntKind) :: n_sscc
   integer (kind=IntKind), allocatable :: ind_sscc(:)
!real (kind=RealKind) :: time_ie
!====================================================================
   ! finding the nonvanishing potential components
   integer (kind=IntKind), allocatable, target   :: flags_jl(:)
   integer (kind=IntKind), allocatable, target   :: flags_trunc_jl(:)
!====================================================================
contains

!include './arrayTools.F90'
include '../lib/arrayTools.F90'
!need to be revised after test finished

   !=================================================================
   subroutine initRelSSSolver(na, lmax, lmax_pot, lmax_green,        &
                              getNumSpecies, getAtomicNumber,        &
                              isBxyz, iprint, stop_routine, evec_in)
   !=================================================================
   use RadialGridModule, only : getGrid

   implicit none

!   integer (kind=IntKind), parameter :: num_vl=12 
  !True only for fcc and bcc, lmax<=4, change with list_vl
   integer (kind=IntKind) :: list_vl(12)


   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: lmax(na)
   integer (kind=IntKind), intent(in) :: lmax_pot(na)
   integer (kind=IntKind), intent(in) :: lmax_green(na)
   integer (kind=IntKind), intent(in) :: iprint(na)
   character (len=*), intent(in) :: stop_routine
   logical, intent(in) :: isBxyz(na)
   real (kind=RealKind), intent(in), optional :: evec_in(3,na)

   integer (kind=IntKind) :: ia, jl,ir,inr,kl
   integer (kind=IntKind) :: kkrsz
   integer (kind=IntKind) :: sz_ind_solutions=0, sz_ind_tmat=0
   integer (kind=IntKind) :: ind_solutions(na), ind_tmat(na)
   integer (kind=IntKind) :: inda, indb
   real (kind=RealKind) :: evec(3)
   integer (kind=IntKind) :: id,lamax,sz,lpot,loc1,loc2,lgreen
!
   real (kind=RealKind), pointer :: p1r(:)
   complex (kind=CmplxKind), pointer :: p1c(:)
!
   interface
      function getNumSpecies(id) result(n)
         use KindParamModule, only : IntKind
         implicit none
         integer (kind=IntKind), intent(in) :: id
         integer (kind=IntKind) :: n
      end function getNumSpecies
!
      function getAtomicNumber(id,ic) result(n)
         use KindParamModule, only : IntKind
         implicit none
         integer (kind=IntKind), intent(in) :: id
         integer (kind=IntKind), intent(in), optional :: ic
         integer (kind=IntKind) :: n
      end function getAtomicNumber
   end interface

   evec=(/1.d0,0.d0,0.d0/) !default direction vector, (z,x,y)
!   list_vl=(/1,17,21,25,39,43,47,65,69,73,77,81/)

   if (Initialized) then
      call ErrorHandler('initRelSSSolver','intend to initialize twice')
   else if (na < 1) then
      call ErrorHandler('initRelSSSolver','LocalNumAtoms < 1',na)
   endif

   LocalNumAtoms = na
   istop=stop_routine

   allocate( AtomicNumber(LocalNumAtoms) )
   allocate( Scatter(LocalNumAtoms) )
   allocate( print_instruction(LocalNumAtoms) )
   allocate( u1(500) )
   allocate( u2(500) )
   allocate( kapdex(98), ldex(98), jdex(98), mydex(98), ldex_bar(98))
   allocate( is_Bxyz(LocalNumAtoms))

   call set_indices(kapdex,ldex,jdex,mydex,ldex_bar,98)
   call clebsch(u1,u2,500)

   do ia=1,na
!      if ( lmax(ia) < 2 ) then
!         Scatter(ia)%num_vl = 1       
!      else if ( lmax(ia) == 2) then
!         Scatter(ia)%num_vl = 4
!      else if ( lmax(ia) == 3 ) then
!         Scatter(ia)%num_vl = 7
!      else if ( lmax(ia) == 4 ) then
!         Scatter(ia)%num_vl = 12
!      else
!         call ErrorHandler ('initRelSSSolver','lmax>4 has not been implemented in Rel: lmax=',lmax(ia))
!      endif
      if ( present(evec_in) ) then !note evec changed from x,y,z to z,x,y in scatter
         Scatter(ia)%evec(1) = evec_in(3,ia)
         Scatter(ia)%evec(2) = evec_in(1,ia)
         Scatter(ia)%evec(3) = evec_in(2,ia)
      else
         Scatter(ia)%evec = evec
      endif
   enddo
   sz_ind_tmat=0
   sz_ind_solutions=0
   do ia=1,LocalNumAtoms
      print_instruction(ia) = iprint(ia)
      AtomicNumber(ia) = getAtomicNumber(ia) ! Needs to be fixed for the CPA case
!     if (AtomicNumber(ia) < 0 .or. AtomicNumber(ia) > 108) then
!         call ErrorHandler('initSSSolver','invalid atomic number',AtomicNumber(ia))
!     endif
      Grid => getGrid(ia)
      if (Grid%nmult .ne. 1) then  !numlt=hin/hout here only single log grid is used
         call ErrorHandler('initRelSSSolver', &
              'nmult must equal 1 for relativistic calculation')
      endif

! determine iend, iend is global, note that it's different from scatter%jend, which depends on atom 
!      iend = max(Grid%jmt+11, Grid%jend)
      iend = Grid%jend !numrs and numrs_cs are the same
!      iend = min(iend, Grid%jend_plus_n)
      iend_max = max(iend_max,iend)
      Scatter(ia)%numrs    = iend
      Scatter(ia)%numrs_cs = Grid%jend
      Scatter(ia)%numrs_mt = Grid%jmt
      Scatter(ia)%numrs_trunc = Grid%jend-Grid%jmt+1
      Scatter(ia)%lmax_pot = lmax_pot(ia)
      Scatter(ia)%lmax_green = lmax_green(ia)
      Scatter(ia)%dx = Grid%hin
      ind_solutions(ia) = sz_ind_solutions + 1

      is_Bxyz(ia) = isBxyz(ia)
      if (is_Bxyz(ia)) then
         is_Bxyz_max=.true.
      endif 
      Scatter(ia)%done = .false.
      Scatter(ia)%l_max=lmax(ia)
      Scatter(ia)%lamax=2*(lmax(ia)+1)**2
      Scatter(ia)%kmax_green=(lmax_green(ia)+1)**2
!print*,"lmax_green, kmax_green=", lmax_green(ia), Scatter(ia)%kmax_green
      kkrsz = 2*(lmax(ia)+1)**2
      ind_tmat(ia) = sz_ind_tmat + 1
      sz_ind_solutions = sz_ind_solutions + iend*kkrsz*kkrsz
      sz_ind_tmat = sz_ind_tmat + kkrsz*kkrsz
   enddo

   allocate( TmpSpace(1:5*iend_max) )

! note here we didn't distinguish kmax_kkr and kmax_int as in SSSolver. We consider only square matrix
   allocate( wks_tmat(sz_ind_tmat) ); wks_tmat = CZERO
   allocate( wks_sinmat(sz_ind_tmat) ); wks_sinmat = CZERO
   allocate( wks_cosmat(sz_ind_tmat) ); wks_cosmat = CZERO
   allocate( wks_Smat(sz_ind_tmat) ); wks_Smat = CZERO
   allocate( wks_jostmat(sz_ind_tmat) ); wks_jostmat = CZERO
   allocate( wks_jinvmat(sz_ind_tmat) ); wks_jinvmat = CZERO
!   allocate( wks_Omega(sz_ind_tmat) ); wks_Omega = CZERO
   allocate( wks_OmegaHat(sz_ind_tmat) ); wks_OmegaHat = CZERO
   allocate( wks_OmegaHatInv(sz_ind_tmat) ); wks_OmegaHatInv = CZERO
   allocate( wks_ssmat_r(sz_ind_solutions) ); wks_ssmat_r = CZERO
   allocate( wks_ccmat_r(sz_ind_solutions) ); wks_ccmat_r = CZERO

   lmax_max = 0
   lmax_pot_max = 0
   do ia=1,LocalNumAtoms
      kkrsz = 2*(lmax(ia)+1)**2
      lmax_max = max(lmax(ia),lmax_max)
      lmax_pot_max = max(lmax_pot(ia),lmax_pot_max)
      inda = ind_solutions(ia)
      indb = inda + (Scatter(ia)%numrs*kkrsz*kkrsz) - 1
      p1c => wks_ssmat_r(inda:indb)
      Scatter(ia)%ss_mat_r => aliasArray3_c(p1c, kkrsz, kkrsz, Scatter(ia)%numrs)
      p1c => wks_ccmat_r(inda:indb)
      Scatter(ia)%cc_mat_r => aliasArray3_c(p1c, kkrsz, kkrsz, Scatter(ia)%numrs)
      inda = ind_tmat(ia) 
!print*,"ind_tmat=",ind_tmat(ia)
      indb = inda + kkrsz*kkrsz - 1
!print*,"indb=",indb
      p1c => wks_tmat(inda:indb)
      Scatter(ia)%t_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
      p1c => wks_Smat(inda:indb)
      Scatter(ia)%S_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
      p1c => wks_sinmat(inda:indb)
      Scatter(ia)%sin_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
      p1c => wks_cosmat(inda:indb)
      Scatter(ia)%cos_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
      p1c => wks_jostmat(inda:indb)
      Scatter(ia)%jost_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
      p1c => wks_jinvmat(inda:indb)
      Scatter(ia)%jinv_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
!     p1c => wks_Omega(inda:indb)
!     Scatter(ia)%Omega_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
      p1c => wks_OmegaHat(inda:indb)
      Scatter(ia)%OmegaHat_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
      p1c => wks_OmegaHatInv(inda:indb)
      Scatter(ia)%OmegaHatInv_mat => aliasArray2_c(p1c,kkrsz,kkrsz)
   enddo
!=============================================================================
! init kyk, note that wks_kyk is different from kyk_wks: wks_kyk is for output 
! to get_kyk
   sz = 0
   do id = 1, LocalNumAtoms !kyk( (lgreen+1)**2,2,lamax,lamax))
      lamax = Scatter(id)%lamax
      lgreen = Scatter(id)%lmax_pot
      sz = sz + 2*(lgreen+1)**2*lamax*lamax
   enddo
   allocate(wks_kyk(sz)); wks_kyk=0.d0
   allocate(wks_kyk_Bz(sz)); wks_kyk_Bz=0.d0
   allocate(wks_kyk_Bx(sz)); wks_kyk_Bx=0.d0
   allocate(wks_kyk_By(sz)); wks_kyk_By=0.d0
   loc1 = 0
   do id = 1, LocalNumAtoms
      lamax = Scatter(id)%lamax
      lgreen = Scatter(id)%lmax_green
      loc2 = loc1 + 2*(lgreen+1)**2*lamax*lamax
      p1r => wks_kyk(loc1+1:loc2)
      Scatter(id)%kyk => aliasArray4_r(p1r,(lgreen+1)**2,2,lamax,lamax)
      p1r => wks_kyk_Bz(loc1+1:loc2)
      Scatter(id)%kyk_Bz => aliasArray4_r(p1r,(lgreen+1)**2,2,lamax,lamax)
      p1r => wks_kyk_Bx(loc1+1:loc2)
      Scatter(id)%kyk_Bx => aliasArray4_r(p1r,(lgreen+1)**2,2,lamax,lamax)
      p1r => wks_kyk_By(loc1+1:loc2)
      Scatter(id)%kyk_By => aliasArray4_r(p1r,(lgreen+1)**2,2,lamax,lamax)
      loc1 = loc2
   enddo

   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()

   allocate( ind_sscc(2*(lmax_max+1)**2) )
!

   Initialized = .true.
   end subroutine initRelSSSolver


   !=================================================================
   subroutine endRelSSSolver()
   !=================================================================

   implicit none

   integer (kind=IntKind) :: ia

!   if (.not. Initialized) then
!      call ErrorHandler('endDiracSolver', 'module not initialized')
!   endif

   Initialized = .false.

   deallocate( AtomicNumber, print_instruction )
   deallocate( kapdex,ldex,jdex,mydex,ldex_bar )
   deallocate( is_Bxyz)
   deallocate( u1 )
   deallocate( u2)

   do ia=1, LocalNumAtoms
      Scatter(ia)%done = .false.
      nullify( Scatter(ia)%t_mat)
      nullify( Scatter(ia)%S_mat )
      nullify( Scatter(ia)%sin_mat )
      nullify( Scatter(ia)%cos_mat )
      nullify( Scatter(ia)%jost_mat )
      nullify( Scatter(ia)%jinv_mat )
!      nullify( Scatter(ia)%Omega_mat )
      nullify( Scatter(ia)%OmegaHat_mat )
      nullify( Scatter(ia)%OmegaHatInv_mat )
      nullify( Scatter(ia)%ss_mat_r)
      nullify( Scatter(ia)%cc_mat_r)
      nullify( Scatter(ia)%green)
      nullify( Scatter(ia)%dos)
      nullify( Scatter(ia)%dos_Bxyz)
      nullify( Scatter(ia)%kyk)
      nullify( Scatter(ia)%kyk_Bx)
      nullify( Scatter(ia)%kyk_By)
      nullify( Scatter(ia)%kyk_Bz)
   enddo
   nullify( nj3, kj3, cgnt )
   deallocate( Scatter )
   deallocate( wks_tmat )
   deallocate( wks_Smat )
   deallocate( wks_sinmat )
   deallocate( wks_cosmat )
   deallocate( wks_jostmat )
   deallocate( wks_jinvmat )
!   deallocate( wks_Omega )
   deallocate( wks_OmegaHat )
   deallocate( wks_OmegaHatInv )
   deallocate( TmpSpace)
   deallocate( wks_green)
   deallocate( wks_dos)
   if (allocated(wks_dos_Bxyz)) then
      deallocate(wks_dos_Bxyz)
   endif
   if (allocated(wks_green_Bxyz)) then
      deallocate(wks_green_Bxyz)
   endif
   deallocate(wks_ccmat_r)
   deallocate(wks_ssmat_r)
!   deallocate( r,vr_s1,vr_s2)
!   deallocate( vr_new,br_new)
!   deallocate(kyk_wks,kyk_short_wks,kyk_Bxyz_wks,kyk_short_Bxyz_wks)
   deallocate( wks_kyk, wks_kyk_Bx, wks_kyk_By, wks_kyk_Bz)
   deallocate( ind_sscc)
   end subroutine endRelSSSolver

!  =================================================================
   subroutine SingleDiracScattering(ia,Ek)
!  =================================================================
   use RadialGridModule, only : getGrid
   use PotentialModule, only : getPotComponentFlag
   use PotentialModule, only : getTruncatedPotComponentFlag
   use PotentialModule, only : getPotential, getTruncatedPotential
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ia
   complex (kind=CmplxKind), intent(in) :: Ek
   integer (kind=IntKind) :: list_vl(12)
   integer (kind=IntKind) :: lmax_kkr
!  integer (kind=IntKind) :: num_vl
   integer (kind=IntKind) :: lpot, jmax_pot, ir, jl, kl
   integer (kind=IntKind) :: lamax
   real (kind=RealKind) :: evec(3)
   complex (kind=CmplxKind), pointer :: pot_jl(:,:)
   complex (kind=CmplxKind), pointer :: pot_trunc_jl(:,:)
   integer (kind=IntKind), pointer :: pflag(:)
   integer (kind=IntKind) :: num_vl_new
   integer (kind=IntKind), allocatable :: list_vl_new(:)
!
   list_vl=(/1,17,21,25,39,43,47,65,69,73,77,81/)

   lmax_kkr=Scatter(ia)%l_max
!   num_vl=Scatter(ia)%num_vl

   lpot = Scatter(ia)%lmax_pot
   jmax_pot = ((lpot+1)*(lpot+2))/2
   Grid => getGrid(ia)
   evec =Scatter(ia)%evec

   pflag => getPotComponentFlag(ia)
   num_vl_new = 0
   if ( SIZE(pflag) .ne. (lpot+1)*(lpot+2)/2 ) then
      call ErrorHandler ("SingleDiracScattering", "length of pflag .ne. (lpot+1)*(lpot+2)/2",&
                         SIZE(pflag),(lpot+1)*(lpot+2)/2)
   endif
!
   do kl = 1, (lpot+1)**2
      if ( pflag(jofk(kl)).ne. 0  ) then
         num_vl_new = num_vl_new + 1
      endif
   enddo
   allocate ( list_vl_new(num_vl_new) )
   jl = 0
   do kl = 1, (lpot+1)**2
      if ( pflag(jofk(kl)) .ne. 0  ) then
         jl= jl+1
         list_vl_new (jl) = kl
      endif
   enddo
   Scatter(ia)%num_vl = num_vl_new
!  =================================================================
!  set up global variables
   nr = Scatter(ia)%numrs_cs
   dx = Scatter(ia)%dx

   allocate( r(nr),vr_s1(jmax_pot,nr),vr_s2(jmax_pot,nr) )
   allocate( vr_new(num_vl_new,nr), br_new(num_vl_new,nr,3) )
   r = Grid%r_mesh
!  allocate bessel functions
      allocate( bjl(0:lmax_kkr+2,1:nr), bnl(0:lmax_kkr+2,1:nr) )

!  Spin up
   pot_jl => getPotential(ia,1,1)
   pot_trunc_jl => getTruncatedPotential(ia,1,1)
!print*,"pot_jl ", pot_jl(1,:)
!print*,"pot_trunc_jl ", pot_trunc_jl(1,:)
!print*, jmax_pot, nr, Scatter(ia)%numrs_trunc
!print*, "Size of r_mesh", SIZE(r)
!   open (500,file="FullPotential_S1",action="write")
   do jl=1,jmax_pot
         do ir=1,nr
         if (ir<= (nr-Scatter(ia)%numrs_trunc)) then
            vr_s1(jl,ir)=r(ir)*pot_jl(ir,jl)
         else
            vr_s1(jl,ir)=r(ir)*pot_trunc_jl(ir-(Scatter(ia)%numrs_mt-1),jl)
         endif
!         write(500,*) ir, r(ir), vr_s1(jl,ir)
         enddo
   enddo
!   close (500)

!  Spin down
   if (is_Bxyz(ia)) then  !if no B field, vr_s2=vr_s1
      pot_jl => getPotential(ia,1,2) !ia is atom number, 1 is species, for cpa
      pot_trunc_jl => getTruncatedPotential(ia,1,2)
   else
      pot_jl => getPotential(ia,1,1)
      pot_trunc_jl => getTruncatedPotential(ia,1,1)
   endif

   do jl=1,jmax_pot
         do ir=1,nr
         if (ir<= (nr-Scatter(ia)%numrs_trunc)) then
            vr_s2(jl,ir)=r(ir)*pot_jl(ir,jl)
         else
            vr_s2(jl,ir)=r(ir)*pot_trunc_jl(ir-(Scatter(ia)%numrs_mt-1),jl)
         endif
         enddo
   enddo

!print*,"nr,jmax_pot, num_vl,lpot=",nr,jmax_pot,num_vl,lpot
!print*,"list_vl",list_vl
!print*,"evec", evec
!print*,"size of vr_new",SIZE(vr_new,1), SIZE(vr_new,2)
!print*,"start getvrbr"
!time_ie = getTime()
   call get_vrbr(vr_s1,vr_s2,vr_new,br_new,evec,jmax_pot,num_vl_new,&
        list_vl_new,nr,lpot)

!print*,"Time usef for get_vrbr",getTime()-time_ie
!print*,"after get_vrbr"
   !======================================
   lamax=2*(lmax_kkr+1)**2
   allocate(kyk_wks( (lpot+1)**2,2,lamax,lamax))
   allocate(kyk_short_wks(num_vl_new,2,lamax,lamax))
   allocate(kyk_Bxyz_wks( (lpot+1)**2,2,lamax,lamax,3)) !the last argument '3' is for 3 magnetic direction
   allocate(kyk_short_Bxyz_wks(num_vl_new,2,lamax,lamax,3))
!print*,"after kyk_wks"
   kyk => kyk_wks
   kyk_short => kyk_short_wks
   kyk_Bxyz => kyk_Bxyz_wks
   kyk_short_Bxyz => kyk_short_Bxyz_wks

   kyk=0.0d0 !never forget to initialize kyk
   call get_kyk(lmax_kkr,lpot,lamax,kyk)
   do kl=1,num_vl_new
      kyk_short(kl,:,:,:)=kyk(list_vl_new(kl),:,:,:)
   enddo
!print*,"kyk_short"
!   open (100,file="kyk",action="write")
!      write(100,*) kyk(:,1,:,3)
!   close(100)

   kyk => kyk_short_wks

   if (is_Bxyz_max) then
      kyk_Bxyz=0.0d0 !never forget to initialize kyk
      call get_kyk_Bxyz(lmax_kkr,lpot,lamax,kyk_Bxyz)
      do kl=1,num_vl_new
         kyk_short_Bxyz(kl,:,:,:,:)=kyk_Bxyz(list_vl_new(kl),:,:,:,:)
      enddo
      kyk_Bxyz => kyk_short_Bxyz
   endif

!  store kyk to Scatter(ia)%kyk
!print*,"kmax_green(ia)=", Scatter(ia)%kmax_green
!note lmax_green (size of wks_kyk) different from lmax_pot (size of kyk_wks), in general
   do kl=1,Scatter(ia)%kmax_green
      Scatter(ia)%kyk(kl,:,:,:) = kyk_wks(kl,:,:,:)
   enddo
   if (is_Bxyz_max) then
      do kl=1, Scatter(ia)%kmax_green
         Scatter(ia)%kyk_Bz(kl,:,:,:) = kyk_Bxyz_wks(kl,:,:,:,1)
         Scatter(ia)%kyk_Bx(kl,:,:,:) = kyk_Bxyz_wks(kl,:,:,:,2)
         Scatter(ia)%kyk_By(kl,:,:,:) = kyk_Bxyz_wks(kl,:,:,:,3)
      enddo
   endif
   if (.not.Scatter(ia)%done ) then 
   ! output when the first time SingleDiracScattering is called for atom ia
   ! Scatter(ia)%done = .true. when SolveDirac is called
      write(6,'(/,80(''-''))')
      write(6,'(/,25x,a)')'*************************************************************'
      write(6,'( 25x,a )')'* Output from SingleDiracScattering (Symmetry of Potential) *'
      write(6,'( 25x,a )')'*************************************************************'
      write(6,'( 4x,a)') 'getPotComponentFlag   :' 
      write(6,'( 4x,t20,a,t40,a)') 'kl','pflag'
      do kl = 1, (lpot+1)**2     
         write(6,'( 4x,t20,i5,t40,i5)') kl, pflag(jofk(kl))
      enddo
      write(6,'(4x,a,t40,a,i3)')'num_vl_new','=', num_vl_new
      write(6,'(4x,t20,a,t40,a)')'i','list_vl_new'
      do kl = 1, num_vl_new
         write(6,'(4x,t20,i5,t40,i5)') kl, list_vl_new(kl)
      enddo
      write(6,'(80(''-''))')
   endif
   !======================================
!print*,"before solveDirac"
   call SolveDirac(ia,Ek,lmax_kkr,lmax_kkr+2,2*(lmax_kkr+1)**2,num_vl_new)
!print*,"before computeRelGF"
   call computeRelGF(ia)
!print*,"after computeRelGF"
   !open (100,file="dos",action="write")
!      do kl=1,(Scatter(ia)%lmax_green+1)*(Scatter(ia)%lmax_green+2)/2
 !           do ir=1,nr
  !             write(100,*) ir,r(ir),Scatter(ia)%dos(ir,kl)
!         enddo
!      enddo
!   close(100)
!
!   open (100,file="green",action="write")
!      do kl=1,(Scatter(ia)%lmax_green+1)*(Scatter(ia)%lmax_green+1)
!            do ir=1,nr
!               write(100,*) ir,r(ir),Scatter(ia)%green(ir,kl)
!         enddo
!      enddo
!   close(100)

   !if (is_Bxyz_max) then
      !open (100,file="green_mm",action="write")
      !   do kl=1,1!(Scatter(ia)%lmax_green+1)*(Scatter(ia)%lmax_green+1)
      !         do ir=1,nr
      !            write(100,*) ir,r(ir),Scatter(ia)%green_Bxyz(ir,kl,3)
   !         enddo
   !      enddo
   !   close(100)
   !endif
!
   deallocate(r,vr_s1,vr_s2,vr_new,br_new)
   deallocate(kyk_wks,kyk_short_wks,kyk_Bxyz_wks,kyk_short_Bxyz_wks)
   deallocate(bjl, bnl)
   deallocate(list_vl_new)
   end subroutine SingleDiracScattering

   !=================================================================
   subroutine get_vrbr(vr_s1,vr_s2,vr_new,br_new,evec,jmax_pot,num_vl,list_vl,nr,lmax_pot)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none

   integer (kind=IntKind) :: jmax_pot,num_vl,nr,list_vl(num_vl),lmax_pot
   real (kind=RealKind) :: evec(3)
   complex (kind=CmplxKind) :: vr_s1(jmax_pot,nr), vr_s2(jmax_pot,nr)
   complex (kind=CmplxKind) :: vr_new(num_vl,nr), br_new(num_vl,nr,3)

   complex (kind=CmplxKind) :: vr_tmp((lmax_pot+1)**2,nr), br_tmp((lmax_pot+1)**2,nr,3)
   complex (kind=CmplxKind) :: vr(jmax_pot,nr), br(jmax_pot,nr,3)
   integer (kind=IntKind) :: ir,kl,jl
!print*,"getvrbr 1"
   do jl=1,jmax_pot
      do ir=1,nr
         vr(jl,ir)=.5d0*(vr_s1(jl,ir) + vr_s2(jl,ir))
         br(jl,ir,1)=.5d0*(vr_s1(jl,ir) - vr_s2(jl,ir))*evec(1)
         br(jl,ir,2)=.5d0*(vr_s1(jl,ir) - vr_s2(jl,ir))*evec(2)
         br(jl,ir,3)=.5d0*(vr_s1(jl,ir) - vr_s2(jl,ir))*evec(3)
      enddo
   enddo
!print*,"getvrbr 2"

   do ir=1,nr
      do kl=1,(lmax_pot+1)**2
         if (mofk(kl)>=0) then
            vr_tmp(kl,ir)=vr(jofk(kl),ir)
            br_tmp(kl,ir,:)=br(jofk(kl),ir,:)

         else
            vr_tmp(kl,ir)=m1m(mofk(kl))*conjg(vr(jofk(kl), ir))
            br_tmp(kl,ir,:)=m1m(mofk(kl))*conjg(br(jofk(kl), ir,:))

         endif
      enddo
   enddo
!print*,"getvrbr 3"
   do ir=1,nr
      do kl=1,num_vl
         vr_new(kl,ir)=vr_tmp(list_vl(kl),ir)
         br_new(kl,ir,:)=br_tmp(list_vl(kl),ir,:)
      enddo
   enddo
!print*,"getvrbr 4"
   end subroutine get_vrbr


   !=================================================================
   subroutine get_kyk(lmax,lmax_pot,lamax,kyk)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none
   character (len=50) :: istop="kyk_FP"
   integer (kind=IntKind) :: lmax,lmax_pot,i,j,m,n,lamax,nr,lp,lpp,kl1
   integer (kind=IntKind) :: lmax_gaunt
   real (kind=RealKind) :: kyk( (lmax_pot+1)**2,2,lamax,lamax )
   integer (kind=IntKind):: lam,lam_new,lambar,lambar_new
   integer (kind=IntKind):: kap_lp,kap_new_lp,my_lp,l_lp,lbar_lp,kapbar_lp,jbar_lp,kapbar_new_lp
   integer (kind=IntKind):: kap_lpp,kap_new_lpp,my_lpp,l_lpp,lbar_lpp,kapbar_lpp,jbar_lpp,kapbar_new_lpp
   integer (kind=IntKind):: k1,k2,lambar_lp,lambar_lpp
   real (kind=RealKind), parameter :: tiny=TEN2m6

   integer (kind=IntKind) :: nj3_i3
   integer (kind=IntKind), pointer :: kj3_i3(:)
   real (kind=RealKind), pointer :: cgnt_i3(:)

   lmax_gaunt=2*lmax !important

   do lp=1,lamax
       kap_lp=kapdex(lp)
       my_lp=mydex(lp)
       l_lp=ldex(lp)

       kapbar_lp=-kap_lp
       lambar_lp=2*kapbar_lp*kapbar_lp+kapbar_lp+(my_lp+1)/2
       lbar_lp=ldex(lambar_lp)
      do lpp=1, lamax
         kap_lpp=kapdex(lpp)
         my_lpp=mydex(lpp)
         l_lpp=ldex(lpp)
   !need to figure out how to get cgnt(:,k1,k2), where max(k1) is half of lamax since no spin
           if (abs(u2(lp))>tiny .and. abs(u2(lpp))>tiny) then !avoid 0 term. l=0 and s=1/2 only gives j=1/2
              k1=l_lp*l_lp+l_lp+(my_lp-1)/2+1
              k2=l_lpp*l_lpp+l_lpp+(my_lpp-1)/2+1

              nj3_i3 = nj3(k2,k1)
              kj3_i3 => kj3(1:nj3_i3,k2,k1)
              cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

              do i=1,nj3_i3
                 kl1=kj3_i3(i)
                 if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
           !k=l^2+(l_m+l+1), when spin up, l_m=(my-1)/2
                    kyk(kl1,1,lpp,lp)=u2(lpp)*u2(lp)*cgnt(i,k2,k1)
                 endif
              enddo
   !           print*,"Spin up: l_lpp=",l_lpp,"my_lpp=",my_lpp,"k2=",k2
   !           print*,"l_lp=",l_lp,"my_lp=",my_lp,"k1=",k1
           endif
           ! when spin down, l_m=(my+1)/2
           if (abs(u1(lp))>tiny .and. abs(u1(lpp))>tiny) then
              k1=l_lp*l_lp+l_lp+(my_lp+1)/2+1
              k2=l_lpp*l_lpp+l_lpp+(my_lpp+1)/2+1

              nj3_i3 = nj3(k2,k1)
              kj3_i3 => kj3(1:nj3_i3,k2,k1)
              cgnt_i3 => cgnt(1:nj3_i3,k2,k1)
              do i=1,nj3_i3
                 kl1=kj3_i3(i)
                 if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
                    kyk(kl1,1,lpp,lp)=kyk(kl1,1,lpp,lp)+u1(lpp)*u1(lp)*cgnt(i,k2,k1)
                 endif
              enddo
   !           print*,"Spin down: l_lpp=",l_lpp,"my_lpp=",my_lpp,"k2=",k2
              !print*,"l_lp=",l_lp,"my_lp=",my_lp,"k1=",k1
           endif

          kapbar_lpp=-kap_lpp
          lambar_lpp=2*kapbar_lpp*kapbar_lpp+kapbar_lpp+(my_lpp+1)/2
   !kappamy = 2*kappa*kappa + kappa + (my+1)/2
          lbar_lpp=ldex(lambar_lpp)

   !first, make sure lbar<=lmax. Then, make sure abs(m_l)<lbar. Finally, remove null CG coefficient.
          if (lbar_lpp<=lmax .and. lbar_lp<=lmax ) then ! avoid l>lmax
             if( abs((my_lp-1)/2)<=lbar_lp .and. abs((my_lpp-1)/2)<=lbar_lpp ) then
                if( abs(u2(lambar_lp))>tiny .and. abs(u2(lambar_lpp))>tiny  ) then
                   k1=lbar_lp*lbar_lp+lbar_lp+(my_lp-1)/2+1
                   k2=lbar_lpp*lbar_lpp+lbar_lpp+(my_lpp-1)/2+1

                   nj3_i3 = nj3(k2,k1)
                   kj3_i3 => kj3(1:nj3_i3,k2,k1)
                   cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

                   do i=1,nj3_i3
                      kl1=kj3_i3(i)
                      if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
   !print*, "k1,k2=",k1,k2
                      kyk(kl1,2,lpp,lp)=u2(lambar_lpp)*u2(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)
                      endif
                   enddo
                endif
             endif

             if( abs(my_lp+1)/2<=lbar_lp .and. abs(my_lpp+1)/2<=lbar_lpp ) then
                if( abs(u1(lambar_lp))>tiny .and. abs(u1(lambar_lpp))>tiny ) then
                   k1=lbar_lp*lbar_lp+lbar_lp+(my_lp+1)/2+1
                   k2=lbar_lpp*lbar_lpp+lbar_lpp+(my_lpp+1)/2+1

                   nj3_i3 = nj3(k2,k1)
                   kj3_i3 => kj3(1:nj3_i3,k2,k1)
                   cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

                   do i=1,nj3_i3
                      kl1=kj3_i3(i)
                      if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
   !print*, "k1,k2=",k1,k2
                      kyk(kl1,2,lpp,lp)=kyk(kl1,2,lpp,lp)+u1(lambar_lpp)*u1(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)
                      endif
                   enddo
                endif
             endif
          endif

      enddo
   enddo

   end subroutine get_kyk


   !=================================================================
   subroutine get_kyk_Bxyz(lmax,lmax_pot,lamax,kyk)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none
   character (len=50) :: istop="get_kyk_Bxyz"
   integer (kind=IntKind) :: lmax,lmax_pot,i,j,m,n,lamax,nr,lp,lpp,kl1
   integer (kind=IntKind) :: lmax_gaunt
   real (kind=RealKind) :: kyk( (lmax_pot+1)**2,2,lamax,lamax,3 )
   integer (kind=IntKind):: lam,lam_new,lambar,lambar_new
   integer (kind=IntKind):: kap_lp,kap_new_lp,my_lp,l_lp,lbar_lp,kapbar_lp,jbar_lp,kapbar_new_lp
   integer (kind=IntKind):: kap_lpp,kap_new_lpp,my_lpp,l_lpp,lbar_lpp,kapbar_lpp,jbar_lpp,kapbar_new_lpp
   integer (kind=IntKind):: k1,k2,lambar_lp,lambar_lpp
   real (kind=RealKind), parameter :: tiny=TEN2m6

   integer (kind=IntKind) :: nj3_i3
   integer (kind=IntKind), pointer :: kj3_i3(:)
   real (kind=RealKind), pointer :: cgnt_i3(:)

   lmax_gaunt=2*lmax !important

   kyk=0. !in case forgot kyk=0 in main function
   do lp=1,lamax
       kap_lp=kapdex(lp)
       my_lp=mydex(lp)
       l_lp=ldex(lp)

       kapbar_lp=-kap_lp
       lambar_lp=2*kapbar_lp*kapbar_lp+kapbar_lp+(my_lp+1)/2
       lbar_lp=ldex(lambar_lp)
      do lpp=1, lamax
         kap_lpp=kapdex(lpp)
         my_lpp=mydex(lpp)
         l_lpp=ldex(lpp)
   !need to figure out how to get cgnt(:,k1,k2), where max(k1) is half of lamax since no spin
           if (abs(u2(lp))>tiny .and. abs(u2(lpp))>tiny) then !avoid 0 term. l=0 and s=1/2 only gives j=1/2
              k1=l_lp*l_lp+l_lp+(my_lp-1)/2+1
              k2=l_lpp*l_lpp+l_lpp+(my_lpp-1)/2+1

              nj3_i3 = nj3(k2,k1)
              kj3_i3 => kj3(1:nj3_i3,k2,k1)
              cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

              do i=1,nj3_i3
                 kl1=kj3_i3(i)
                 if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
           !k=l^2+(l_m+l+1), when spin up, l_m=(my-1)/2
                    kyk(kl1,1,lpp,lp,1)=u2(lpp)*u2(lp)*cgnt(i,k2,k1)
                 endif
              enddo
   !           print*,"Spin up: l_lpp=",l_lpp,"my_lpp=",my_lpp,"k2=",k2
   !           print*,"l_lp=",l_lp,"my_lp=",my_lp,"k1=",k1
           endif
           ! when spin down, l_m=(my+1)/2
           if (abs(u1(lp))>tiny .and. abs(u1(lpp))>tiny) then
              k1=l_lp*l_lp+l_lp+(my_lp+1)/2+1
              k2=l_lpp*l_lpp+l_lpp+(my_lpp+1)/2+1

              nj3_i3 = nj3(k2,k1)
              kj3_i3 => kj3(1:nj3_i3,k2,k1)
              cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

              do i=1,nj3_i3
                 kl1=kj3_i3(i)
                 if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
                    kyk(kl1,1,lpp,lp,1)=kyk(kl1,1,lpp,lp,1) - u1(lpp)*u1(lp)*cgnt(i,k2,k1) !change due to Bz
                 endif
              enddo
   !           print*,"Spin down: l_lpp=",l_lpp,"my_lpp=",my_lpp,"k2=",k2
              !print*,"l_lp=",l_lp,"my_lp=",my_lp,"k1=",k1
           endif
   !==============================================================================
   !Added component for Bx and By. U2(lp) and U1(lpp)
           if (abs(u2(lp))>tiny .and. abs(u1(lpp))>tiny) then !avoid 0 term. l=0 and s=1/2 only gives j=1/2
              k1=l_lp*l_lp+l_lp+(my_lp - 1)/2+1
              k2=l_lpp*l_lpp+l_lpp+(my_lpp + 1)/2+1

              nj3_i3 = nj3(k2,k1)
              kj3_i3 => kj3(1:nj3_i3,k2,k1)
              cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

              do i=1,nj3_i3
                 kl1=kj3_i3(i)
                 if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
           !k=l^2+(l_m+l+1), when spin up, l_m=(my-1)/2
                    kyk(kl1,1,lpp,lp,2)=u2(lp)*u1(lpp)*cgnt(i,k2,k1)
                    kyk(kl1,1,lpp,lp,3)=u2(lp)*u1(lpp)*cgnt(i,k2,k1)*(-1.)
                 endif
              enddo
   !           print*,"Spin up: l_lpp=",l_lpp,"my_lpp=",my_lpp,"k2=",k2
   !           print*,"l_lp=",l_lp,"my_lp=",my_lp,"k1=",k1
           endif
           if (abs(u1(lp))>tiny .and. abs(u2(lpp))>tiny) then !avoid 0 term. l=0 and s=1/2 only gives j=1/2
              k1=l_lp*l_lp+l_lp+(my_lp + 1)/2+1
              k2=l_lpp*l_lpp+l_lpp+(my_lpp - 1)/2+1

              nj3_i3 = nj3(k2,k1)
              kj3_i3 => kj3(1:nj3_i3,k2,k1)
              cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

              do i=1,nj3_i3
                 kl1=kj3_i3(i)
                 if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
           !k=l^2+(l_m+l+1), when spin up, l_m=(my-1)/2
                    kyk(kl1,1,lpp,lp,2)=kyk(kl1,1,lpp,lp,2)+u1(lp)*u2(lpp)*cgnt(i,k2,k1)
                    kyk(kl1,1,lpp,lp,3)=kyk(kl1,1,lpp,lp,3)+u1(lp)*u2(lpp)*cgnt(i,k2,k1)
                 endif
              enddo
   !           print*,"Spin up: l_lpp=",l_lpp,"my_lpp=",my_lpp,"k2=",k2
   !           print*,"l_lp=",l_lp,"my_lp=",my_lp,"k1=",k1
           endif
   !==============================================================================

          kapbar_lpp=-kap_lpp
          lambar_lpp=2*kapbar_lpp*kapbar_lpp+kapbar_lpp+(my_lpp+1)/2
   !kappamy = 2*kappa*kappa + kappa + (my+1)/2
          lbar_lpp=ldex(lambar_lpp)


   !first, make sure lbar<=lmax. Then, make sure abs(m_l)<lbar. Finally, remove null CG coefficient.
          if (lbar_lpp<=lmax .and. lbar_lp<=lmax ) then ! avoid l>lmax
             if( abs((my_lp-1)/2)<=lbar_lp .and. abs((my_lpp-1)/2)<=lbar_lpp ) then
                if( abs(u2(lambar_lp))>tiny .and. abs(u2(lambar_lpp))>tiny  ) then
                   k1=lbar_lp*lbar_lp+lbar_lp+(my_lp-1)/2+1
                   k2=lbar_lpp*lbar_lpp+lbar_lpp+(my_lpp-1)/2+1

                   nj3_i3 = nj3(k2,k1)
                   kj3_i3 => kj3(1:nj3_i3,k2,k1)
                   cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

                   do i=1,nj3_i3
                      kl1=kj3_i3(i)
                      if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
   !print*, "k1,k2=",k1,k2
                      kyk(kl1,2,lpp,lp,1)=u2(lambar_lpp)*u2(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)
                      endif
                   enddo
                endif
             endif

             if( abs(my_lp+1)/2<=lbar_lp .and. abs(my_lpp+1)/2<=lbar_lpp ) then
                if( abs(u1(lambar_lp))>tiny .and. abs(u1(lambar_lpp))>tiny ) then
                   k1=lbar_lp*lbar_lp+lbar_lp+(my_lp+1)/2+1
                   k2=lbar_lpp*lbar_lpp+lbar_lpp+(my_lpp+1)/2+1

                   nj3_i3 = nj3(k2,k1)
                   kj3_i3 => kj3(1:nj3_i3,k2,k1)
                   cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

                   do i=1,nj3_i3
                      kl1=kj3_i3(i)
                      if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
   !print*, "k1,k2=",k1,k2
   !change due to Bz, because VIII=-sigma*B,will be done in get_ffgg
                      kyk(kl1,2,lpp,lp,1)=kyk(kl1,2,lpp,lp,1) - u1(lambar_lpp)*u1(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)
                      endif
                   enddo
                endif
             endif
   !===================================================================================
   !Added component for Bx and By. U2(lp) and U1(lpp)
             if( abs((my_lp - 1)/2)<=lbar_lp .and. abs((my_lpp + 1)/2)<=lbar_lpp ) then
                if( abs(u2(lambar_lp))>tiny .and. abs(u1(lambar_lpp))>tiny  ) then
                   k1=lbar_lp*lbar_lp+lbar_lp+(my_lp - 1)/2+1
                   k2=lbar_lpp*lbar_lpp+lbar_lpp+(my_lpp + 1)/2+1

                   nj3_i3 = nj3(k2,k1)
                   kj3_i3 => kj3(1:nj3_i3,k2,k1)
                   cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

                   do i=1,nj3_i3
                      kl1=kj3_i3(i)
                      if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
   !print*, "k1,k2=",k1,k2
                      kyk(kl1,2,lpp,lp,2)=u1(lambar_lpp)*u2(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)
                      kyk(kl1,2,lpp,lp,3)=u1(lambar_lpp)*u2(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)*(-1.)
                      endif
                   enddo
                endif
             endif

             if( abs((my_lp + 1)/2)<=lbar_lp .and. abs((my_lpp - 1)/2)<=lbar_lpp ) then
                if( abs(u1(lambar_lp))>tiny .and. abs(u2(lambar_lpp))>tiny  ) then
                   k1=lbar_lp*lbar_lp+lbar_lp+(my_lp + 1)/2+1
                   k2=lbar_lpp*lbar_lpp+lbar_lpp+(my_lpp - 1)/2+1

                   nj3_i3 = nj3(k2,k1)
                   kj3_i3 => kj3(1:nj3_i3,k2,k1)
                   cgnt_i3 => cgnt(1:nj3_i3,k2,k1)

                   do i=1,nj3_i3
                      kl1=kj3_i3(i)
                      if ( kj3_i3(i)<=(lmax_pot+1)**2 ) then
   !print*, "k1,k2=",k1,k2
                      kyk(kl1,2,lpp,lp,2)=kyk(kl1,2,lpp,lp,2)+u2(lambar_lpp)*u1(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)
                      kyk(kl1,2,lpp,lp,3)=kyk(kl1,2,lpp,lp,3)+u2(lambar_lpp)*u1(lambar_lp)*cgnt(i,k2,k1)&
                      *sign(1,kap_lp)*sign(1,kap_lpp)
                      endif
                   enddo
                endif
             endif
   !===================================================================================
          endif

      enddo
   enddo
   end subroutine get_kyk_Bxyz

   !=================================================================
   subroutine calNum_kyk(lmax_pot,lamax,kyk_in,n_kyk,ind_kyk)
   !=================================================================
   ! return n_kyk and ind_kyk. n_kyk is the number of third component of 
   ! the gaunt factors, ind_kyk are the indices
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none

   integer (kind=IntKind), intent(in) :: lmax_pot,lamax
   real (kind=RealKind), intent(in) :: kyk_in((lmax_pot+1)**2,lamax,lamax)
   integer (kind=IntKind), intent(out) :: n_kyk(lamax,lamax)
   integer (kind=IntKind), intent(out) :: ind_kyk((lmax_pot+1)**2,lamax,lamax)
   integer (kind=IntKind) :: lp,lpp,i
   real (kind=RealKind), parameter :: tiny=TEN2m6

   n_kyk=0
   ind_kyk=0
   do lp=1,lamax
      do lpp=1,lamax
         n_kyk(lpp,lp)=0
         do i=1,(lmax_pot+1)**2
            if (Abs(kyk_in(i,lpp,lp))>tiny) then
               n_kyk(lpp,lp)=n_kyk(lpp,lp)+1
               ind_kyk(n_kyk(lpp,lp),lpp,lp)=i
            endif
         enddo
      enddo
   enddo


   end subroutine calNum_kyk

   !=================================================================
   subroutine get_ffgg(bj_tmp,bn_tmp,kappa,lam,lamax,r,ss_vec,cc_vec,vr,lmax_p2,&
   ff,gg,kyk_vec,num_vl,coeff)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none

   integer (kind=IntKind) :: lam,l,lbar,lambar,lam_new,lambar_new,lamax,lmax_p2
   integer (kind=IntKind) :: lpp,l_lpp,kj3,l_lppbar
   integer (kind=IntKind) :: num_vl

   complex (kind=CmplxKind) :: ss_vec(lamax), cc_vec(lamax), ff,gg, kappa
   complex (kind=CmplxKind) :: vr(num_vl)
   complex (kind=CmplxKind) :: coeff(2)
   complex (kind=CmplxKind) :: bj_tmp(0:lmax_p2),bn_tmp(0:lmax_p2)
   complex (kind=CmplxKind) :: nvkykSa,nvkykSb,nvkykSc,nvkykSd, vkyk

   real (kind=RealKind) :: r
   integer (kind=IntKind):: kap,my,kapbar
   real (kind=RealKind) :: kyk_vec(num_vl,2,lamax)


   l=ldex(lam)
   kap=kapdex(lam)
   my=mydex(lam)

   if (kap>0) then
      lbar=l-1
   else if (kap<0) then
      lbar=l+1
   else
      print*,"error:kapdex(lam)==0 in get_ff"
      stop
   endif

   ff=(0.d0,0.d0)
   gg=(0.d0,0.d0)

   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
   !upper part
   do lpp=1,lamax
      if (cc_vec(lpp) .ne. CZERO) then
         l_lpp=ldex(lpp)
         vkyk=(0.d0,0.d0)
         vkyk=vkyk+dot_product(vr,kyk_vec(:,1,lpp)) 
         nvkykSa=nvkykSa+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
         nvkykSb=nvkykSb+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
         nvkykSc=nvkykSc+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
         nvkykSd=nvkykSd+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
      endif
   enddo
   ff=ff+coeff(1)*r*r*bj_tmp(l)*(nvkykSa-nvkykSb)
   gg=gg+coeff(1)*r*r*bn_tmp(l)*(nvkykSc-nvkykSd)

!print*,"check upper: bn_tmp=",bn_tmp(0),"cc_vec(1)=",cc_vec(1),"kyk_vec=",kyk_vec(1,1,1),"vr=",vr(1)


   !lower part, lbar
   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
   do lpp=1,lamax
!      l_lpp=ldex(lpp)
!      if (kapdex(lpp)>0) then
!         l_lppbar=l_lpp-1
!      else
!         l_lppbar=l_lpp+1
!      endif
      if (cc_vec(lpp) .ne. CZERO) then
         l_lppbar=ldex_bar(lpp) !a faster way to find l_lppbar
         vkyk=(0.d0,0.d0)
         vkyk=vkyk+dot_product(vr,kyk_vec(:,2,lpp)) 
         nvkykSa=nvkykSa+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
         nvkykSb=nvkykSb+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
         nvkykSc=nvkykSc+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
         nvkykSd=nvkykSd+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
      endif
   enddo
   ff=ff+coeff(2)*r*r*bj_tmp(lbar)*(nvkykSa-nvkykSb)
   gg=gg+coeff(2)*r*r*bn_tmp(lbar)*(nvkykSc-nvkykSd)

!   print*,"check lower: gg=",gg,"r=",r,"nvkykSc=",nvkykSc,"nvkykSd=",nvkykSd

   end subroutine get_ffgg
   !=================================================================

   !=================================================================
   subroutine get_ffgg_speed(bj_tmp,bn_tmp,kappa,lam,lamax,r,ss_vec,cc_vec,vr,lmax_p2,&
   ff,gg,kyk_vec,num_vl,coeff)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none

   integer (kind=IntKind) :: lam,l,lbar,lambar,lam_new,lambar_new,lamax,lmax_p2
   integer (kind=IntKind) :: lpp,l_lpp,kj3,l_lppbar
   integer (kind=IntKind) :: num_vl

   complex (kind=CmplxKind) :: ss_vec(lamax), cc_vec(lamax), ff,gg, kappa
   complex (kind=CmplxKind) :: vr(num_vl)
   complex (kind=CmplxKind) :: coeff(2)
   complex (kind=CmplxKind) :: bj_tmp(0:lmax_p2),bn_tmp(0:lmax_p2)
   complex (kind=CmplxKind) :: nvkykSa,nvkykSb,nvkykSc,nvkykSd, vkyk

   real (kind=RealKind) :: r
   integer (kind=IntKind):: kap,my,kapbar, lam_loc
   real (kind=RealKind) :: kyk_vec(num_vl,2,lamax)


   l=ldex(lam)
   kap=kapdex(lam)
   my=mydex(lam)

   if (kap>0) then
      lbar=l-1
   else if (kap<0) then
      lbar=l+1
   else
      print*,"error:kapdex(lam)==0 in get_ff"
      stop
   endif

   ff=(0.d0,0.d0)
   gg=(0.d0,0.d0)

   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
   !upper part
   do lam_loc=1,n_sscc
      lpp=ind_sscc(lam_loc)
      if (cc_vec(lpp) .ne. CZERO) then
         l_lpp=ldex(lpp)
         vkyk=(0.d0,0.d0)
         vkyk=vkyk+dot_product(vr,kyk_vec(:,1,lpp)) 
         nvkykSa=nvkykSa+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
         nvkykSb=nvkykSb+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
         nvkykSc=nvkykSc+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
         nvkykSd=nvkykSd+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
      endif
   enddo
   ff=ff+coeff(1)*r*r*bj_tmp(l)*(nvkykSa-nvkykSb)
   gg=gg+coeff(1)*r*r*bn_tmp(l)*(nvkykSc-nvkykSd)

!print*,"check upper: bn_tmp=",bn_tmp(0),"cc_vec(1)=",cc_vec(1),"kyk_vec=",kyk_vec(1,1,1),"vr=",vr(1)


   !lower part, lbar
   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
   do lam_loc=1,n_sscc
      lpp=ind_sscc(lam_loc)
!      l_lpp=ldex(lpp)
!      if (kapdex(lpp)>0) then
!         l_lppbar=l_lpp-1
!      else
!         l_lppbar=l_lpp+1
!      endif
      if (cc_vec(lpp) .ne. CZERO) then
         l_lppbar=ldex_bar(lpp) !a faster way to find l_lppbar
         vkyk=(0.d0,0.d0)
         vkyk=vkyk+dot_product(vr,kyk_vec(:,2,lpp)) 
         nvkykSa=nvkykSa+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
         nvkykSb=nvkykSb+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
         nvkykSc=nvkykSc+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
         nvkykSd=nvkykSd+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
      endif
   enddo
   ff=ff+coeff(2)*r*r*bj_tmp(lbar)*(nvkykSa-nvkykSb)
   gg=gg+coeff(2)*r*r*bn_tmp(lbar)*(nvkykSc-nvkykSd)

!   print*,"check lower: gg=",gg,"r=",r,"nvkykSc=",nvkykSc,"nvkykSd=",nvkykSd

   end subroutine get_ffgg_speed
   !=================================================================


   !=================================================================
   subroutine get_ffgg_B(bj_tmp,bn_tmp,kappa,lam,lamax,r,ss_vec,cc_vec,vr,br,lmax_p2,&
   ff,gg,kyk_vec, kyk_vec_Bxyz, num_vl,coeff)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none

   integer (kind=IntKind) :: lam,l,lbar,lambar,lam_new,lambar_new,lamax,lmax_p2
   integer (kind=IntKind) :: lpp,l_lpp,kj3,l_lppbar
   integer (kind=IntKind) :: num_vl

   complex (kind=CmplxKind) :: ss_vec(lamax), cc_vec(lamax), ff,gg, kappa
   complex (kind=CmplxKind) :: vr(num_vl),br(num_vl,3)
   complex (kind=CmplxKind) :: coeff(2)
   complex (kind=CmplxKind) :: bj_tmp(0:lmax_p2),bn_tmp(0:lmax_p2)
   complex (kind=CmplxKind) :: nvkykSa,nvkykSb,nvkykSc,nvkykSd, vkyk

   real (kind=RealKind) :: r
   integer (kind=IntKind):: kap,my,kapbar
   real (kind=RealKind) :: kyk_vec(num_vl,2,lamax),kyk_vec_Bxyz(num_vl,2,lamax,3)

   integer (kind=IntKind) :: lam_loc !note lam is input, so here use lam_loc

   l=ldex(lam)
   kap=kapdex(lam)
   my=mydex(lam)

   if (kap>0) then
      lbar=l-1
   else if (kap<0) then
      lbar=l+1
   else
      print*,"error:kapdex(lam)==0 in get_ff"
      stop
   endif

   ff=(0.d0,0.d0)
   gg=(0.d0,0.d0)

   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
   !upper part
   if ( .false.) then
      do lam_loc=1,n_sscc
         lpp=ind_sscc(lam_loc)
            l_lpp=ldex(lpp)
            vkyk=(0.d0,0.d0)
            vkyk=vkyk+dot_product(vr,kyk_vec(:,1,lpp)) + dot_product(br(:,1),kyk_vec_Bxyz(:,1,lpp,1))&
                 +dot_product(br(:,2),kyk_vec_Bxyz(:,1,lpp,2))+dot_product(br(:,3),kyk_vec_Bxyz(:,1,lpp,3))*sqrtm1
            nvkykSa=nvkykSa+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
            nvkykSb=nvkykSb+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
            nvkykSc=nvkykSc+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
            nvkykSd=nvkykSd+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
      enddo
   else
      do lpp=1,lamax
         if (cc_vec(lpp) .ne. CZERO) then
            l_lpp=ldex(lpp)
            vkyk=(0.d0,0.d0)
            vkyk=vkyk+dot_product(vr,kyk_vec(:,1,lpp)) + dot_product(br(:,1),kyk_vec_Bxyz(:,1,lpp,1))&
                 +dot_product(br(:,2),kyk_vec_Bxyz(:,1,lpp,2))+dot_product(br(:,3),kyk_vec_Bxyz(:,1,lpp,3))*sqrtm1
            nvkykSa=nvkykSa+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
            nvkykSb=nvkykSb+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
            nvkykSc=nvkykSc+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
            nvkykSd=nvkykSd+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
         endif
      enddo
   endif
   ff=ff+coeff(1)*r*r*bj_tmp(l)*(nvkykSa-nvkykSb)
   gg=gg+coeff(1)*r*r*bn_tmp(l)*(nvkykSc-nvkykSd)


   !lower part, lbar
   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
   if ( .false.) then
      do lam_loc=1,n_sscc
         lpp=ind_sscc(lam_loc)
   !      l_lpp=ldex(lpp)
   !      if (kapdex(lpp)>0) then
   !         l_lppbar=l_lpp-1
   !      else
   !         l_lppbar=l_lpp+1
   !      endif
   !      if (cc_vec(lpp) .ne. CZERO) then
            l_lppbar=ldex_bar(lpp) !a faster way to find l_lppbar
            vkyk=(0.d0,0.d0)
            vkyk=vkyk+dot_product(vr,kyk_vec(:,2,lpp)) - dot_product(br(:,1),kyk_vec_Bxyz(:,2,lpp,1))&
                 - dot_product(br(:,2),kyk_vec_Bxyz(:,2,lpp,2))- dot_product(br(:,3),kyk_vec_Bxyz(:,2,lpp,3))*sqrtm1
            nvkykSa=nvkykSa+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
            nvkykSb=nvkykSb+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
            nvkykSc=nvkykSc+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
            nvkykSd=nvkykSd+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
   !      endif
      enddo

   else
      do lpp=1,lamax
   !      l_lpp=ldex(lpp)
   !      if (kapdex(lpp)>0) then
   !         l_lppbar=l_lpp-1
   !      else
   !         l_lppbar=l_lpp+1
   !      endif
         if (cc_vec(lpp) .ne. CZERO) then
            l_lppbar=ldex_bar(lpp) !a faster way to find l_lppbar
            vkyk=(0.d0,0.d0)
            vkyk=vkyk+dot_product(vr,kyk_vec(:,2,lpp)) - dot_product(br(:,1),kyk_vec_Bxyz(:,2,lpp,1))&
                 - dot_product(br(:,2),kyk_vec_Bxyz(:,2,lpp,2))- dot_product(br(:,3),kyk_vec_Bxyz(:,2,lpp,3))*sqrtm1
            nvkykSa=nvkykSa+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
            nvkykSb=nvkykSb+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
            nvkykSc=nvkykSc+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
            nvkykSd=nvkykSd+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
         endif
      enddo
   endif
   ff=ff+coeff(2)*r*r*bj_tmp(lbar)*(nvkykSa-nvkykSb)
   gg=gg+coeff(2)*r*r*bn_tmp(lbar)*(nvkykSc-nvkykSd)

   !print*,"check lower: gg=",gg,"r=",r,"nvkykSc=",nvkykSc,"nvkykSd=",nvkykSd

   end subroutine get_ffgg_B
   !=================================================================


   !=================================================================
   subroutine get_ffgg_B_speed(bj_tmp,bn_tmp,kappa,lam,lamax,r,ss_vec,cc_vec,vr,br,lmax_p2,&
   ff,gg,kyk_vec, kyk_vec_Bxyz, num_vl,coeff)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none

   integer (kind=IntKind) :: lam,l,lbar,lambar,lam_new,lambar_new,lamax,lmax_p2
   integer (kind=IntKind) :: lpp,l_lpp,kj3,l_lppbar
   integer (kind=IntKind) :: num_vl

   complex (kind=CmplxKind) :: ss_vec(lamax), cc_vec(lamax), ff,gg, kappa
   complex (kind=CmplxKind) :: vr(num_vl),br(num_vl,3)
   complex (kind=CmplxKind) :: coeff(2)
   complex (kind=CmplxKind) :: bj_tmp(0:lmax_p2),bn_tmp(0:lmax_p2)
   complex (kind=CmplxKind) :: nvkykSa,nvkykSb,nvkykSc,nvkykSd, vkyk

   real (kind=RealKind) :: r
   integer (kind=IntKind):: kap,my,kapbar
   real (kind=RealKind) :: kyk_vec(num_vl,2,lamax),kyk_vec_Bxyz(num_vl,2,lamax,3)

   integer (kind=IntKind) :: lam_loc !note lam is input, so here use lam_loc

   l=ldex(lam)
   kap=kapdex(lam)
   my=mydex(lam)

   if (kap>0) then
      lbar=l-1
   else if (kap<0) then
      lbar=l+1
   else
      print*,"error:kapdex(lam)==0 in get_ff"
      stop
   endif

   ff=(0.d0,0.d0)
   gg=(0.d0,0.d0)

   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
   !upper part
      do lam_loc=1,n_sscc
         lpp=ind_sscc(lam_loc)
            l_lpp=ldex(lpp)
            vkyk=(0.d0,0.d0)
            vkyk=vkyk+dot_product(vr,kyk_vec(:,1,lpp)) + dot_product(br(:,1),kyk_vec_Bxyz(:,1,lpp,1))&
                 +dot_product(br(:,2),kyk_vec_Bxyz(:,1,lpp,2))+dot_product(br(:,3),kyk_vec_Bxyz(:,1,lpp,3))*sqrtm1
            nvkykSa=nvkykSa+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
            nvkykSb=nvkykSb+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
            nvkykSc=nvkykSc+bn_tmp(l_lpp)*vkyk*ss_vec(lpp)
            nvkykSd=nvkykSd+bj_tmp(l_lpp)*vkyk*cc_vec(lpp)
      enddo
   ff=ff+coeff(1)*r*r*bj_tmp(l)*(nvkykSa-nvkykSb)
   gg=gg+coeff(1)*r*r*bn_tmp(l)*(nvkykSc-nvkykSd)


   !lower part, lbar
   nvkykSa=(0.d0,0.d0)
   nvkykSb=(0.d0,0.d0)
   nvkykSc=(0.d0,0.d0)
   nvkykSd=(0.d0,0.d0)
      do lam_loc=1,n_sscc
         lpp=ind_sscc(lam_loc)
   !      l_lpp=ldex(lpp)
   !      if (kapdex(lpp)>0) then
   !         l_lppbar=l_lpp-1
   !      else
   !         l_lppbar=l_lpp+1
   !      endif
   !      if (cc_vec(lpp) .ne. CZERO) then
            l_lppbar=ldex_bar(lpp) !a faster way to find l_lppbar
            vkyk=(0.d0,0.d0)
            vkyk=vkyk+dot_product(vr,kyk_vec(:,2,lpp)) - dot_product(br(:,1),kyk_vec_Bxyz(:,2,lpp,1))&
                 - dot_product(br(:,2),kyk_vec_Bxyz(:,2,lpp,2))- dot_product(br(:,3),kyk_vec_Bxyz(:,2,lpp,3))*sqrtm1
            nvkykSa=nvkykSa+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
            nvkykSb=nvkykSb+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
            nvkykSc=nvkykSc+bn_tmp(l_lppbar)*vkyk*ss_vec(lpp)
            nvkykSd=nvkykSd+bj_tmp(l_lppbar)*vkyk*cc_vec(lpp)
   !      endif
      enddo
   ff=ff+coeff(2)*r*r*bj_tmp(lbar)*(nvkykSa-nvkykSb)
   gg=gg+coeff(2)*r*r*bn_tmp(lbar)*(nvkykSc-nvkykSd)

   !print*,"check lower: gg=",gg,"r=",r,"nvkykSc=",nvkykSc,"nvkykSd=",nvkykSd

   end subroutine get_ffgg_B_speed
   !=================================================================


   !=================================================================
   subroutine SolveDirac(ia,Ek,lmax,lmax_p2,lamax,num_vl)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind

   implicit none

   integer (kind=IntKind), intent(in) :: ia, lmax, lmax_p2, lamax, num_vl
   complex (kind=CmplxKind), intent(in) :: Ek
   integer (kind=IntKind) :: ir, lam, lam_0, lp
   real (kind=RealKind) :: r_tmp
   complex (kind=CmplxKind) :: kappa, energy, coeff(2) 
   complex (kind=CmplxKind) :: ss_mat(lamax,nr), cc_mat(lamax,nr) !note nr is a global variable
!   complex (kind=CmplxKind) :: ss_mat_r(lamax,nr,lamax), cc_mat_r(lamax,nr,lamax)
   complex (kind=CmplxKind) :: bjtmp(0:lmax_p2), bntmp(0:lmax_p2) !,bjl(0:lmax_p2,1:nr), bnl(0:lmax_p2,1:nr)
   complex (kind=CmplxKind) :: vr_tmp(num_vl), br_tmp(num_vl,3)
   complex (kind=CmplxKind) :: K1(lamax),L1(lamax),K2(lamax),L2(lamax),K3(lamax),L3(lamax),K4(lamax),L4(lamax)
   complex (kind=CmplxKind) :: ff,gg 
   complex (kind=CmplxKind) :: ff_list(0:4),gg_list(0:4)
   complex (kind=CmplxKind) :: t_matrix(lamax,lamax),cone_diag(lamax,lamax),jinv_matrix(lamax,lamax)
   complex (kind=CmplxKind) :: s_matrix(lamax,lamax),c_matrix(lamax,lamax)
   complex (kind=CmplxKind) :: s_matrix_d(lamax,lamax),c_matrix_d(lamax,lamax)
   complex (kind=CmplxKind) :: OmegaHat(lamax,lamax),OmegaHatInv(lamax,lamax)
!  speedup variables

   !k=sqrt( E_k*(2m)+E_k**2/c**2)
   kappa=sqrt(2.d0*Me*Ek + Ek**2/LightSpeed**2)
   do ir=1,nr
      call SphericalBessel(lmax_p2,kappa*r(ir),bjl(:,ir))
      call SphericalNeumann(lmax_p2,kappa*r(ir),bnl(:,ir))
   enddo
   energy=sqrt(kappa**2*LightSpeed**2+(Me*LightSpeed*LightSpeed)**2)
   coeff(1)=kappa*(energy + Me*LightSpeed**2)/LightSpeed**2
   coeff(2)=kappa*(energy - Me*LightSpeed**2)/LightSpeed**2
!  assign globle values kappa_g, energy_g
   kappa_g = kappa
   energy_g = energy !energy_g=Ek+Me*C^2


!   ss_mat_r=CZERO
!   cc_mat_r=CZERO
   s_matrix=CZERO
   c_matrix=CZERO
   do lam_0=1,lamax
   !   print*, "lam_0=",lam_0
      ss_mat=( 0.0d0,0.0d0)
      cc_mat=( 0.0d0,0.0d0)
      cc_mat(lam_0,1)=(-1.0d0,0.0d0) !boundry condition at origin
      if (is_Bxyz(ia)) then
         do ir=2,5
            bjtmp=.5d0*(bjl(:,ir)+bjl(:,ir-1))
            bntmp=.5d0*(bnl(:,ir)+bnl(:,ir-1))
            do lp=1,lamax
               call get_ffgg_B (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),cc_mat(:,ir-1),&
           vr_new(:,ir-1),br_new(:,ir-1,:),lmax_p2,ff,gg,kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
               K1(lp)=ff*dx
               L1(lp)=gg*dx
            enddo

            r_tmp=r(ir-1)*exp(0.5d0*dx)
            vr_tmp(:)=vr_new(:,ir-1)+(r_tmp-r(ir-1))/(r(ir)-r(ir-1))*(vr_new(:,ir)-vr_new(:,ir-1)) 
            br_tmp(:,:)=br_new(:,ir-1,:)+(r_tmp-r(ir-1))/(r(ir)-r(ir-1))*(br_new(:,ir,:)-br_new(:,ir-1,:))
            !not exact, but fine for the first five points
            do lp=1,lamax
               call get_ffgg_B (bjtmp,bntmp,kappa,lp,lamax,r_tmp,(ss_mat(:,ir-1)+0.5*K1(:)),&
           (cc_mat(:,ir-1)+0.5d0*L1(:)),vr_tmp(:),br_tmp(:,:),lmax_p2,ff,gg,kyk(:,:,:,lp),&
           kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
               K2(lp)=ff*dx
               L2(lp)=gg*dx
            enddo
         
            do lp=1,lamax
               call get_ffgg_B (bjtmp,bntmp,kappa,lp,lamax,r_tmp,(ss_mat(:,ir-1)+0.5*K2(:)),&
           (cc_mat(:,ir-1)+0.5d0*L2(:)),vr_tmp(:),br_tmp(:,:),lmax_p2,ff,gg,kyk(:,:,:,lp),&
           kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
               K3(lp)=ff*dx
               L3(lp)=gg*dx
            enddo

            do lp=1,lamax
               call get_ffgg_B (bjl(:,ir),bnl(:,ir),kappa,lp,lamax,r(ir),(ss_mat(:,ir-1)+K3(:)), &
            (cc_mat(:,ir-1)+L3(:)),vr_new(:,ir),br_new(:,ir,:),lmax_p2,ff,gg,kyk(:,:,:,lp),&
            kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
               K4(lp)=ff*dx
               L4(lp)=gg*dx
            enddo

            do lam=1,lamax
               ss_mat(lam,ir)=ss_mat(lam,ir-1)+1.d0/6.d0*(K1(lam)+2.d0*K2(lam)+2.d0*K3(lam)+k4(lam))
               cc_mat(lam,ir)=cc_mat(lam,ir-1)+1.d0/6.d0*(L1(lam)+2.d0*L2(lam)+2.d0*L3(lam)+L4(lam))
            enddo

         enddo
! Find the Speedup flags 
         if ( Speedup ) then
            n_sscc = 0
            ind_sscc = 0
            do lam=1,lamax
               if ( (cc_mat(lam, ir-1) .ne. CZERO ) .or. (ss_mat(lam,ir-1) .ne. CZERO) ) then
                  n_sscc = n_sscc + 1
                  ind_sscc(n_sscc) = lam
               endif
            enddo
         !   if(lam_0 == 1 ) then
         !   open (102,file="flag_sscc",action="write")
         !   write (102,*) n_sscc 
         !   write (102,*) ind_sscc
         !   close (102)
         !   endif
         endif

         do ir=6,nr
!           bjtmp=bjl(:,ir)
!           bntmp=bnl(:,ir)
            if ( Speedup ) then
                do lam=1,n_sscc !change made for speedup
                  lp=ind_sscc(lam) !change made for speedup
                  call get_ffgg_B_speed (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),br_new(:,ir-1,:),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B_speed (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),br_new(:,ir-2,:),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B_speed (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),br_new(:,ir-3,:),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B_speed (bjl(:,ir-4),bnl(:,ir-4),kappa,lp,lamax,r(ir-4),ss_mat(:,ir-4),&
                   cc_mat(:,ir-4),vr_new(:,ir-4),br_new(:,ir-4,:),lmax_p2,&
                   ff_list(4),gg_list(4),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(55.d0*ff_list(1)-59.d0*ff_list(2)+&
                            37.d0*ff_list(3)-9.d0*ff_list(4))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(55.d0*gg_list(1)-59.d0*gg_list(2)+&
                            37.d0*gg_list(3)-9.d0*gg_list(4))
                  enddo
                do lam=1,n_sscc
                  lp=ind_sscc(lam)
                  call get_ffgg_B_speed (bjl(:,ir),bnl(:,ir),kappa,lp,lamax,r(ir),ss_mat(:,ir),&
                   cc_mat(:,ir),vr_new(:,ir),br_new(:,ir,:),lmax_p2,&
                   ff_list(0),gg_list(0),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B_speed (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),br_new(:,ir-1,:),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B_speed (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),br_new(:,ir-2,:),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B_speed (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),br_new(:,ir-3,:),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(9.d0*ff_list(0)+19.d0*ff_list(1)&
                            -5.d0*ff_list(2)+1.d0*ff_list(3))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(9.d0*gg_list(0)+19.d0*gg_list(1)&
                            -5.d0*gg_list(2)+1.d0*gg_list(3))
                  enddo
            else !whether speedup
                  do lp=1,lamax
                  call get_ffgg_B (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),br_new(:,ir-1,:),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),br_new(:,ir-2,:),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),br_new(:,ir-3,:),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B (bjl(:,ir-4),bnl(:,ir-4),kappa,lp,lamax,r(ir-4),ss_mat(:,ir-4),&
                   cc_mat(:,ir-4),vr_new(:,ir-4),br_new(:,ir-4,:),lmax_p2,&
                   ff_list(4),gg_list(4),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(55.d0*ff_list(1)-59.d0*ff_list(2)+&
                            37.d0*ff_list(3)-9.d0*ff_list(4))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(55.d0*gg_list(1)-59.d0*gg_list(2)+&
                            37.d0*gg_list(3)-9.d0*gg_list(4))
               enddo

               do lp=1,lamax
                  call get_ffgg_B (bjl(:,ir),bnl(:,ir),kappa,lp,lamax,r(ir),ss_mat(:,ir),&
                   cc_mat(:,ir),vr_new(:,ir),br_new(:,ir,:),lmax_p2,&
                   ff_list(0),gg_list(0),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),br_new(:,ir-1,:),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),br_new(:,ir-2,:),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
                  call get_ffgg_B (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),br_new(:,ir-3,:),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),kyk_Bxyz(:,:,:,lp,:),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(9.d0*ff_list(0)+19.d0*ff_list(1)&
                            -5.d0*ff_list(2)+1.d0*ff_list(3))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(9.d0*gg_list(0)+19.d0*gg_list(1)&
                            -5.d0*gg_list(2)+1.d0*gg_list(3))
               enddo
            endif !speedup
         enddo

      else !isBxyz or not

         do ir=2,5
            bjtmp=.5d0*(bjl(:,ir)+bjl(:,ir-1))
            bntmp=.5d0*(bnl(:,ir)+bnl(:,ir-1))
            do lp=1,lamax
               call get_ffgg (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),cc_mat(:,ir-1),&
           vr_new(:,ir-1),lmax_p2,ff,gg,kyk(:,:,:,lp),num_vl,coeff)
               K1(lp)=ff*dx
               L1(lp)=gg*dx
            enddo

            r_tmp=r(ir-1)*exp(0.5d0*dx)
            vr_tmp(:)=vr_new(:,ir-1)+(r_tmp-r(ir-1))/(r(ir)-r(ir-1))*(vr_new(:,ir)-vr_new(:,ir-1)) 
            !not exact, but fine for the first five points
            do lp=1,lamax
               call get_ffgg (bjtmp,bntmp,kappa,lp,lamax,r_tmp,(ss_mat(:,ir-1)+0.5*K1(:)),&
           (cc_mat(:,ir-1)+0.5d0*L1(:)),vr_tmp(:),lmax_p2,ff,gg,kyk(:,:,:,lp),num_vl,coeff)
               K2(lp)=ff*dx
               L2(lp)=gg*dx
            enddo
         
            do lp=1,lamax
               call get_ffgg (bjtmp,bntmp,kappa,lp,lamax,r_tmp,(ss_mat(:,ir-1)+0.5*K2(:)),&
           (cc_mat(:,ir-1)+0.5d0*L2(:)),vr_tmp(:),lmax_p2,ff,gg,kyk(:,:,:,lp),num_vl,coeff)
               K3(lp)=ff*dx
               L3(lp)=gg*dx
            enddo

            do lp=1,lamax
               call get_ffgg (bjl(:,ir),bnl(:,ir),kappa,lp,lamax,r(ir),(ss_mat(:,ir-1)+K3(:)), &
            (cc_mat(:,ir-1)+L3(:)),vr_new(:,ir),lmax_p2,ff,gg,kyk(:,:,:,lp),num_vl,coeff)
               K4(lp)=ff*dx
               L4(lp)=gg*dx
            enddo

            do lam=1,lamax
               ss_mat(lam,ir)=ss_mat(lam,ir-1)+1.d0/6.d0*(K1(lam)+2.d0*K2(lam)+2.d0*K3(lam)+k4(lam))
               cc_mat(lam,ir)=cc_mat(lam,ir-1)+1.d0/6.d0*(L1(lam)+2.d0*L2(lam)+2.d0*L3(lam)+L4(lam))
            enddo

         enddo

         if ( Speedup ) then
            n_sscc = 0
            ind_sscc = 0
            do lam=1,lamax
               if ( (cc_mat(lam, ir-1) .ne. CZERO ) .or. (ss_mat(lam,ir-1) .ne. CZERO) ) then
                  n_sscc = n_sscc + 1
                  ind_sscc(n_sscc) = lam
               endif
            enddo
         endif

         do ir=6,nr
!           bjtmp=bjl(:,ir)
!           bntmp=bnl(:,ir)

            if ( Speedup ) then
                do lam=1,n_sscc
                  lp=ind_sscc(lam)
                  call get_ffgg_speed (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg_speed (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg_speed (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg_speed (bjl(:,ir-4),bnl(:,ir-4),kappa,lp,lamax,r(ir-4),ss_mat(:,ir-4),&
                   cc_mat(:,ir-4),vr_new(:,ir-4),lmax_p2,&
                   ff_list(4),gg_list(4),kyk(:,:,:,lp),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(55.d0*ff_list(1)-59.d0*ff_list(2)+&
                            37.d0*ff_list(3)-9.d0*ff_list(4))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(55.d0*gg_list(1)-59.d0*gg_list(2)+&
                            37.d0*gg_list(3)-9.d0*gg_list(4))
               enddo

                do lam=1,n_sscc
                  lp=ind_sscc(lam)
                  call get_ffgg_speed (bjl(:,ir),bnl(:,ir),kappa,lp,lamax,r(ir),ss_mat(:,ir),&
                   cc_mat(:,ir),vr_new(:,ir),lmax_p2,&
                   ff_list(0),gg_list(0),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg_speed (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg_speed (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg_speed (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(9.d0*ff_list(0)+19.d0*ff_list(1)&
                            -5.d0*ff_list(2)+1.d0*ff_list(3))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(9.d0*gg_list(0)+19.d0*gg_list(1)&
                            -5.d0*gg_list(2)+1.d0*gg_list(3))
              enddo
           else !speedup or not
               do lp=1,lamax
                  call get_ffgg (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg (bjl(:,ir-4),bnl(:,ir-4),kappa,lp,lamax,r(ir-4),ss_mat(:,ir-4),&
                   cc_mat(:,ir-4),vr_new(:,ir-4),lmax_p2,&
                   ff_list(4),gg_list(4),kyk(:,:,:,lp),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(55.d0*ff_list(1)-59.d0*ff_list(2)+&
                            37.d0*ff_list(3)-9.d0*ff_list(4))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(55.d0*gg_list(1)-59.d0*gg_list(2)+&
                            37.d0*gg_list(3)-9.d0*gg_list(4))
               enddo

               do lp=1,lamax
                  call get_ffgg (bjl(:,ir),bnl(:,ir),kappa,lp,lamax,r(ir),ss_mat(:,ir),&
                   cc_mat(:,ir),vr_new(:,ir),lmax_p2,&
                   ff_list(0),gg_list(0),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg (bjl(:,ir-1),bnl(:,ir-1),kappa,lp,lamax,r(ir-1),ss_mat(:,ir-1),&
                   cc_mat(:,ir-1),vr_new(:,ir-1),lmax_p2,&
                   ff_list(1),gg_list(1),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg (bjl(:,ir-2),bnl(:,ir-2),kappa,lp,lamax,r(ir-2),ss_mat(:,ir-2),&
                   cc_mat(:,ir-2),vr_new(:,ir-2),lmax_p2,&
                   ff_list(2),gg_list(2),kyk(:,:,:,lp),num_vl,coeff)
                  call get_ffgg (bjl(:,ir-3),bnl(:,ir-3),kappa,lp,lamax,r(ir-3),ss_mat(:,ir-3),&
                   cc_mat(:,ir-3),vr_new(:,ir-3),lmax_p2,&
                   ff_list(3),gg_list(3),kyk(:,:,:,lp),num_vl,coeff)
              ss_mat(lp,ir)=ss_mat(lp,ir-1)+dx/24.d0*(9.d0*ff_list(0)+19.d0*ff_list(1)&
                            -5.d0*ff_list(2)+1.d0*ff_list(3))
              cc_mat(lp,ir)=cc_mat(lp,ir-1)+dx/24.d0*(9.d0*gg_list(0)+19.d0*gg_list(1)&
                            -5.d0*gg_list(2)+1.d0*gg_list(3))
              enddo
           endif !speedup
        enddo

      endif
      do ir=1,nr
         Scatter(ia)%ss_mat_r(:,lam_0,ir)=ss_mat(:,ir)
         Scatter(ia)%cc_mat_r(:,lam_0,ir)=cc_mat(:,ir)
      enddo
      s_matrix(:,lam_0)=ss_mat(:,nr)
      c_matrix(:,lam_0)=cc_mat(:,nr)
   enddo

! dot notation. This is not true when there is By component
   s_matrix_d=s_matrix
   c_matrix_d=c_matrix
   
!   do ir=1,5
!      print*, "ss and cc", ss_mat(50,ir), cc_mat(50,ir)
!      print*, "ff,gg", ff, gg
!   enddo


!   open(111,file="reg_ss_mat",action="write")
!   write(111,*),Scatter(ia)%ss_mat_r(1,1,:)
!   close(111)
!   open(111,file="reg_cc_mat",action="write")
!   write(111,*),Scatter(ia)%cc_mat_r(1,1,:)
!   close(111)

   Scatter(ia)%sin_mat=s_matrix!ss_mat_r(:,nr,:)
   Scatter(ia)%cos_mat=c_matrix!cc_mat_r(:,nr,:)
!   open(111,file="ss_mat_nr",action="write")
!   write(111,*), Scatter(ia)%sin_mat
!   close(111)
!   open(111,file="cc_mat_nr",action="write")
!   write(111,*), Scatter(ia)%cos_mat
!   close(111)

   call cal_tmat(ia,kappa,lamax,s_matrix,c_matrix,s_matrix_d,c_matrix_d,&
   t_matrix,jinv_matrix,OmegaHat,OmegaHatInv)
   Scatter(ia)%t_mat=t_matrix
   cone_diag=CZERO
   do lam=1,lamax
      cone_diag(lam,lam)=CONE
   enddo
   Scatter(ia)%S_mat=cone_diag-2*sqrtm1*kappa*t_matrix
   Scatter(ia)%jost_mat=sqrtm1*s_matrix-c_matrix
   Scatter(ia)%jinv_mat=jinv_matrix
   Scatter(ia)%OmegaHat_mat=OmegaHat
   Scatter(ia)%OmegaHatInv_mat=OmegaHatInv
   Scatter(ia)%done=.true.

 !  open(131,file="t_mat",action="write")
 !  write(131,*) Scatter(ia)%t_mat
 !  close(131)

 !  open(131,file="jinv",action="write")
 !  write(131,*) Scatter(ia)%jinv_mat
 !  close(131)

 !  open(131,file="OmegaHat",action="write")
 !  write(131,*) Scatter(ia)%OmegaHat_mat
 !  close(131)

 !  open(131,file="OmegaHatInv",action="write")
 !  write(131,*) Scatter(ia)%OmegaHatInv_mat
 !  close(131)
!   call cal_tmat_0(ia,kappa,lamax,s_matrix,c_matrix,t_matrix)

   end subroutine SolveDirac


   !=================================================================
   subroutine cal_tmat(ia,kappa,lamax,s_matrix,c_matrix,s_matrix_d,c_matrix_d,&
   t_matrix,jinv_matrix,OmegaHat,OmegaHatInv)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none
   !calculate t_mat, t_mat=1/k*S*(I*S-C)^-1

   integer (kind=IntKind),intent(in) :: lamax,ia
   complex (kind=CmplxKind),intent(in) :: kappa
   complex (kind=CmplxKind),intent(in) :: s_matrix(lamax,lamax),c_matrix(lamax,lamax)
   complex (kind=CmplxKind),intent(in) :: s_matrix_d(lamax,lamax),c_matrix_d(lamax,lamax)
   complex (kind=CmplxKind),intent(out) :: t_matrix(lamax,lamax),jinv_matrix(lamax,lamax)
   complex (kind=CmplxKind),intent(out) :: OmegaHat(lamax,lamax),OmegahatInv(lamax,lamax)
   complex (kind=CmplxKind) :: mat_tmp(lamax,lamax),kappa_dig(lamax,lamax),kinv

   integer (kind=IntKind) :: INFO
   integer (kind=IntKind), allocatable :: IPVT(:)
   integer (kind=IntKind) :: LWORK
   complex (kind=CmplxKInd), allocatable :: WORK(:)

   mat_tmp=(sqrtm1*s_matrix-c_matrix)

   allocate( IPVT(1:2*lamax) ) !2 times larger incase we need to use scalapack
   call ZGETRF(lamax,lamax,mat_tmp,lamax,IPVT,INFO)

   LWORK=lamax*lamax
   allocate( WORK(1:LWORK))
   WORK=(0.d0,0.d0)

   call ZGETRI(lamax,mat_tmp,lamax,IPVT,WORK,LWORK,INFO)
   jinv_matrix=mat_tmp

!   open(131,file="mat_tmp",action="write")
!   write(131,*) mat_tmp
!   close(131)


   kinv=1/kappa
   call zgemm('n','n',lamax,lamax,lamax,kinv,&
   s_matrix, lamax, mat_tmp, lamax, czero, t_matrix, lamax)

!  calculate OmegaHat and OmegaHatInv
   OmegaHatInv=CZERO
!  OmegaHatInv=S^T(I*S-C)
   call zgemm( 't', 'n', lamax, lamax, lamax, sqrtm1,&
   s_matrix_d, lamax, s_matrix, lamax, CZERO, OmegaHatInv, lamax)

!   open(131,file="STS",action="write")
!   write(131,*) OmegaHatInv
!   close(131)

   call zgemm( 't', 'n', lamax, lamax, lamax, -CONE,&
   s_matrix_d, lamax, c_matrix, lamax, CONE, OmegaHatInv, lamax) !be careful with zgemm. CONR or CZERo

!   open(131,file="OmegaHatInv",action="write")
!   write(131,*) OmegaHatInv
!   close(131)

   mat_tmp=OmegaHatInv
      call ZGETRF(lamax,lamax,mat_tmp,lamax,IPVT,INFO)
   call ZGETRI(lamax,mat_tmp,lamax,IPVT,WORK,LWORK,INFO)
   OmegaHat=mat_tmp

!   open(131,file="OmegaHat02",action="write")
!   write(131,*) OmegaHat
!   close(131)
   
   deallocate(IPVT,WORK)
   end subroutine cal_tmat


   !=================================================================
   subroutine cal_tmat_0(ia,kappa,lamax,s_matrix,c_matrix,t_matrix)
   !=================================================================
   !a different way to calculate tmatrix, just for test
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none
   !calculate t_mat, first calculate t_mat^-1=-k*C_mat*S_mat^-1+i*k

   integer (kind=IntKind),intent(in) :: lamax,ia
   complex (kind=CmplxKind),intent(in) :: kappa
   integer (kind=IntKind) :: i
   complex (kind=CmplxKind) :: s_matrix(lamax,lamax),c_matrix(lamax,lamax),t_matrix(lamax,lamax)
   complex (kind=CmplxKind) :: mat_tmp(lamax,lamax),kappa_dig(lamax,lamax),kappa_tmp

   integer (kind=IntKind) :: INFO
   integer (kind=IntKind), allocatable :: IPVT(:)
   integer (kind=IntKind) :: LWORK
   complex (kind=CmplxKInd), allocatable :: WORK(:)
   complex (kind=CmplxKind), parameter :: sqrtm1 = ( 0.0d0,1.0d0)

   kappa_dig=(0.d0,0.d0)
   do i=1,lamax
   kappa_dig(i,i)=kappa
   enddo

   !print*,kappa_dig

   mat_tmp=s_matrix

   allocate( IPVT(1:2*lamax) ) !2 times larger incase we need to use scalapack
   call ZGETRF(lamax,lamax,mat_tmp,lamax,IPVT,INFO)



   LWORK=lamax*lamax
   allocate( WORK(1:LWORK))
   WORK=(0.d0,0.d0)

   call ZGETRI(lamax,mat_tmp,lamax,IPVT,WORK,LWORK,INFO)


   !print*, mat_tmp
   !mat_tmp=c_mat*s_mat^-1

   !print*,"mat-tmp(3,3)=",mat_tmp(3,3)
   !print*,"mat-tmp(6,3)=",mat_tmp(6,3)


   kappa_tmp=-kappa

   call zgemm( 'n', 'n', lamax, lamax, lamax, kappa_tmp,       &
   c_matrix, lamax, mat_tmp, lamax, sqrtm1, kappa_dig, lamax)
   t_matrix=kappa_dig


   WORK=(0.d0,0.d0)
   call ZGETRF(lamax,lamax,t_matrix,lamax,IPVT,INFO)
   call ZGETRI(lamax,t_matrix,lamax,IPVT,WORK,LWORK,INFO)

   !output t_matrix
!   open(131,file="t_mat_01",action="write")
!   write(131,*) t_matrix
!   close(131)

   deallocate(IPVT,WORK)
   end subroutine cal_tmat_0


   !=================================================================
   subroutine clebsch(u1,u2,n)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none

   integer (kind=IntKind):: i, l, m, inr, kap, ir
   integer (kind=IntKind), intent(in):: n
   real (kind=RealKind):: twolp1
   real (kind=RealKind):: u1(n),u2(n)

   !  clebsch-gordan rectangular matrices to transform from (lm) to
   ! (kappa,my) basis

   do i=1,n
      u1(i)=0
      u2(i)=0
   enddo

   do l=0,12
      twolp1=dfloat(2*l+1)
      ! j=l-1/2
      kap=l

      do m=-l,l
         if(kap.eq.0) goto 102
         ! ms=-1/2
         ir=2*kap*kap+kap+m
         u1(ir)=sqrt((l+m)/twolp1)
         ! ms=+1/2
         ir=2*kap*kap+kap+m+1
         u2(ir)=-sqrt((l-m)/twolp1)
      enddo

      102  continue
      ! j=l+1/2
      kap=-l-1
      do m=-l,l
         ! ms=-1/2
         ir=2*kap*kap+kap+m
         u1(ir)=sqrt((l-m+1)/twolp1)
         ! ms=+1/2
         ir=2*kap*kap+kap+m+1
         u2(ir)=sqrt((l+m+1)/twolp1)
      enddo

   enddo
   end subroutine clebsch

   !=================================================================
   subroutine set_indices(kapdex,ldex,jdex,mydex,ldex_bar,n)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none
   integer (kind=IntKind),intent(in)::n
   integer (kind=IntKind):: kapdex(n), ldex(n), ldex_bar(n)
   integer (kind=IntKind):: jdex(n), mydex(n)
   kapdex=(/                          &
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
   -7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7/)
   ldex=(/                            &
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
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6/)
   jdex=(/                            &
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
   13,13,13,13,13,13,13,13,13,13,13,13,13,13/)
   mydex=(/                                        &
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
   -13,-11, -9, -7, -5, -3, -1,  1,  3,  5,  7,  9,  11, 13/)
   ldex_bar=(/                            &
   1, 1,                               &
   0, 0,                               &
   2, 2, 2, 2,                         &
   1, 1, 1, 1,                         &
   3, 3, 3, 3, 3, 3,                   &
   2, 2, 2, 2, 2, 2,                   &
   4, 4, 4, 4, 4, 4, 4, 4,             &
   3, 3, 3, 3, 3, 3, 3, 3,             &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       &
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       &
   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
   5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
   7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7/)
   end subroutine set_indices

!  ===================================================================
   function computeRelSingleSiteDOS(id) result(dos)
!  ===================================================================
!  
!  This function returns single site DOS for a given energy on the real
!  energy axis.
!  
!  ===================================================================
   use IntegrationModule, only : calIntegration
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration

   implicit none
!  
   integer (kind=IntKind), intent(in) :: id
!  
   real (kind=RealKind) :: dos, dos_mt
   integer (kind=IntKind) :: ir, nr, jmax_dos
   real (kind=RealKind) :: dos_r(Scatter(id)%numrs_cs), dos_int(Scatter(id)%numrs_cs)
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind) :: energy !Total energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   logical :: SphericalInt=.false., truncated !to test getVolumeIntegration, .false. by default

   kappa = kappa_g
   energy = energy_g
   nr=Scatter(id)%numrs_cs


   Grid => getGrid(id)
   dos_r_jl => Scatter(id)%dos


   if (SphericalInt) then
      do ir=1,nr
         dos_r(ir)=Grid%r_mesh(ir)*Grid%r_mesh(ir)*Real(dos_r_jl(ir,1))*Sqrt(4.d0*PI)
         !r^2*dos*Y_0/Y_0^2, note scatter%dos are basically GF expanded on spherical harmonics
      enddo
      dos_int=0.d0
      call calIntegration(0, nr, Grid%r_mesh, dos_r(:), dos_int(:))
      dos=dos_int(nr)
   else
      jmax_dos = size(dos_r_jl,2)
      dos = getVolumeIntegration( id, nr, Grid%r_mesh,   &
                           jmax_dos, 2, dos_r_jl, dos_mt, truncated=.true. )
   endif

   end function computeRelSingleSiteDOS
!  ===================================================================

!  ===================================================================
   function computeRelSingleSiteMMDOS(id) result(dos)
!  ===================================================================
!  
!  This function returns single site magnetic moment DOS for a given energy on the real
!  energy axis.
!  
!  ===================================================================
   use IntegrationModule, only : calIntegration
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration

   implicit none
!  
   integer (kind=IntKind), intent(in) :: id
!  
   real (kind=RealKind) :: dos(3), dos_mt(3)
   integer (kind=IntKind) :: ir, nr, jmax_dos, i
   real (kind=RealKind) :: dos_r(Scatter(id)%numrs_cs), dos_int(Scatter(id)%numrs_cs)
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind) :: energy !Total energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:,:)
   logical :: SphericalInt=.false., truncated !to test getVolumeIntegration, .false. by default

   kappa = kappa_g
   energy = energy_g
   nr=Scatter(id)%numrs_cs


   Grid => getGrid(id)
   dos_r_jl => Scatter(id)%dos_Bxyz


   if (is_Bxyz_max) then
      do i=1,3
         jmax_dos = size(dos_r_jl,2)
         dos(i) = getVolumeIntegration( id, nr, Grid%r_mesh,   &
                              jmax_dos, 2, dos_r_jl(:,:,i), dos_mt(i), truncated=.true. )
      enddo
   else
      call ErrorHandler ('computeRelSingleSiteMMDOS',&
      'non-magnetic calculation, but magnetic moment density DOS is called')
   endif

   end function computeRelSingleSiteMMDOS
!  ===================================================================

   !=================================================================
   subroutine computeRelDOS(atom,add_highl_fec)
   !=================================================================
   use GroupCommModule, only : GlobalSumInGroup
   use RadialGridModule, only : getGrid
!
   implicit none

   logical, intent(in), optional :: add_highl_fec   ! Add high L contribution by free electron
   integer (kind=IntKind), intent(in) :: atom
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: ir, jl, kl, klc, l, m, kl1, kl1p, kl2, kl2c, kl2p, &
         kl2pc, m1, m2, m2p, ma, i, l2
   integer (kind=IntKind) :: js, ia, nr, lm, jm, sz, loc1, loc2
   integer (kind=IntKind) :: lmax_dos, jmax_dos, kmax_kkr_loc
   complex (kind=CmplxKind), pointer :: dos(:,:),dos_Bxyz(:,:,:)
   complex (kind=CmplxKind) :: cfac, kappa, energy
   complex (kind=CmplxKind), pointer :: green(:,:),green_Bxyz(:,:,:), p1c(:)

   id = atom
   kappa=kappa_g
   energy=energy_g

   if (.not. allocated(wks_dos)) then
      sz = 0
      do ia = 1, LocalNumAtoms
         nr = Scatter(ia)%numrs_cs
         lm = Scatter(ia)%lmax_green
         jm = (lm+1)*(lm+2)/2
         sz = sz + nr*jm
      enddo
      allocate( wks_dos(sz) )
      if ( is_Bxyz_max ) then
            allocate( wks_dos_Bxyz(sz*3) )
      endif
      loc1 = 0
      do ia = 1, LocalNumAtoms
         nr = Scatter(ia)%numrs_cs
         lm = Scatter(ia)%lmax_green
         jm = (lm+1)*(lm+2)/2
         loc2 = loc1 + nr*jm
         p1c => wks_dos(loc1+1:loc2)
         Scatter(ia)%dos => aliasArray2_c(p1c,nr,jm)
         Scatter(ia)%dos = CZERO
         loc1 = loc2
      enddo

      if ( is_Bxyz_max ) then
         loc1 = 0
            do ia = 1, LocalNumAtoms
            nr = Scatter(ia)%numrs_cs
            lm = Scatter(ia)%lmax_green
            jm = (lm+1)*(lm+2)/2
               loc2 = loc1 + nr*jm*3
            p1c => wks_dos_Bxyz(loc1+1:loc2)
            Scatter(ia)%dos_Bxyz => aliasArray3_c(p1c,nr,jm,3)
            Scatter(ia)%dos_Bxyz = CZERO
            loc1 = loc2
         enddo
      endif
   endif

   lmax_dos = Scatter(id)%lmax_green
   jmax_dos = (lmax_dos+1)*(lmax_dos+2)/2
   green => Scatter(id)%green
   dos => Scatter(id)%dos
   nr = Scatter(id)%numrs_cs
   Grid => getGrid(id)
   do jl = 1, jmax_dos
      kl = kofj(jl)
      m = mofj(jl)
      klc = kl-2*m
      do ir = 1, nr
         dos(ir,jl) = green(ir,kl) - m1m(m)*conjg(green(ir,klc))
      enddo
   enddo
   cfac = sqrtm1/PI2 !divid by 2*Pi, the 2 is from 1/2*(green-conjg(green))
!
   dos = cfac*dos

   if (is_Bxyz_max) then
      green_Bxyz => Scatter(id)%green_Bxyz
      dos_Bxyz => Scatter(id)%dos_Bxyz
      do jl = 1, jmax_dos
         kl = kofj(jl)
         m = mofj(jl)
         klc = kl-2*m
         do i = 1, 3
            do ir = 1, nr   
               dos_Bxyz(ir,jl,i) = green_Bxyz(ir,kl,i) - m1m(m)*conjg(green_Bxyz(ir,klc,i))
            enddo
         enddo
      enddo
      cfac = sqrtm1/PI2 !divid by 2*Pi, the 2 is from 1/2*(green-conjg(green))
   !
      dos_Bxyz = cfac*dos_Bxyz
   endif
!
   if ( present(add_highl_fec) ) then
      if (add_highl_fec) then
         TmpSpace = CZERO
         do l = lofk(kmax_kkr_loc), 0, -1
            do ir = 1, nr
               TmpSpace(ir) = TmpSpace(ir) + (2*l+1)*bjl(ir,l)**2/kappa**2
            enddo
         enddo
         cfac = kappa*sqrt(PI)/(2.0d0*PI**2)*(4.d0*energy/LightSpeed**2)
         ! =(4*W/c^2)*kappa/(4.0d0*PI**2)/Y00
         do ir = 1, nr
            dos(ir,1) = dos(ir,1) + cfac*(Grid%r_mesh(ir)**2-TmpSpace(ir))
         enddo
      endif
   endif

   end subroutine computeRelDOS

   !=================================================================
   subroutine computeRelGF(ia)
   !=================================================================
   use RadialGridModule, only : getGrid
   implicit none

   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: sz, loc1, loc2
   integer (kind=IntKind) :: nofr, km
   integer (kind=IntKind) :: kap_lp, my_lp, l_lp, kapbar_lp, lambar_lp, lbar_lp
   integer (kind=IntKind) :: l_lpp, lbar_lpp
   integer (kind=IntKind) :: i,ir, lp, lpp, lamax, k_1, k_2, kl1, lmax, lmax_pot
!   integer (kind=IntKind) :: nj3_i3
!   integer (kind=IntKind), pointer :: kj3_i3(:)
!   real (kind=RealKind), pointer :: cgnt_i3(:)
   complex (kind=CmplxKind), pointer :: ss_mat_r(:,:,:)
   complex (kind=CmplxKind), pointer :: cc_mat_r(:,:,:)
   complex (kind=CmplxKind), pointer :: green(:,:), green_Bxyz(:,:,:)
   complex (kind=CmplxKind) :: kappa
   complex (kind=CmplxKind) :: energy !Total energy
   complex (kind=CmplxKind) :: bjtmp( 0:Scatter(ia)%l_max+2 )
   complex (kind=CmplxKind) :: bntmp( 0:Scatter(ia)%l_max+2 )
   complex (kind=CmplxKind) :: mid_matrix( Scatter(ia)%lamax,Scatter(ia)%lamax,4)
   complex (kind=CmplxKind) :: prefactor,prefactor_bar
   !complex (kind=CmplxKind), pointer :: t_matrix(:,:)
   !complex (kind=CmplxKind), pointer :: s_matrix(:,:)
   !complex (kind=CmplxKind), pointer :: c_matrix(:,:)
   complex (kind=CmplxKind), pointer :: OmegaHat(:,:)
!   complex (kind=CmplxKind), pointer :: d_s_matrix(:,:)
!   complex (kind=CmplxKind), pointer :: d_c_matrix(:,:)
   complex (kind=CmplxKind), pointer :: d_ss_mat_r(:,:,:)
   complex (kind=CmplxKind), pointer :: d_cc_mat_r(:,:,:)
   complex (kind=CmplxKind), pointer :: p1c(:)

   integer (kind=IntKind) :: n_kyk(Scatter(ia)%lamax,Scatter(ia)%lamax)

   integer (kind=IntKind) :: ind_kyk((Scatter(ia)%lmax_pot+1)**2,Scatter(ia)%lamax,Scatter(ia)%lamax)

   integer (kind=IntKind) :: n_kyk_bar(Scatter(ia)%lamax,Scatter(ia)%lamax)

   integer (kind=IntKind) :: ind_kyk_bar((Scatter(ia)%lmax_pot+1)**2,&
                             Scatter(ia)%lamax,Scatter(ia)%lamax)
   integer (kind=IntKind), allocatable :: n_kyk_Bxyz(:,:,:)
   integer (kind=IntKind), allocatable :: ind_kyk_Bxyz(:,:,:,:)
   integer (kind=IntKind), allocatable :: n_kyk_bar_Bxyz(:,:,:)
   integer (kind=IntKind), allocatable :: ind_kyk_bar_Bxyz(:,:,:,:)

   complex (kind=CmplxKind) :: tmpI

   if (is_Bxyz_max) then
      allocate( n_kyk_Bxyz(Scatter(ia)%lamax,Scatter(ia)%lamax,3) )
      allocate( ind_kyk_Bxyz((Scatter(ia)%lmax_pot+1)**2,Scatter(ia)%lamax,Scatter(ia)%lamax,3) )
      allocate( n_kyk_bar_Bxyz(Scatter(ia)%lamax,Scatter(ia)%lamax,3) )
      allocate( ind_kyk_bar_Bxyz((Scatter(ia)%lmax_pot+1)**2,Scatter(ia)%lamax,Scatter(ia)%lamax,3) )
   endif

   OmegaHat => Scatter(ia)%OmegaHat_mat
   !s_matrix => Scatter(ia)%sin_mat
   !c_matrix => Scatter(ia)%cos_mat
   ss_mat_r => Scatter(ia)%ss_mat_r
   cc_mat_r => Scatter(ia)%cc_mat_r

   !Right now no distinction between dot and undot sine matrix, to be generalized
   !d_s_matrix => Scatter(ia)%sin_mat
   !d_c_matrix => Scatter(ia)%cos_mat
   d_ss_mat_r => Scatter(ia)%ss_mat_r
   d_cc_mat_r => Scatter(ia)%cc_mat_r

   kappa = kappa_g
   energy = energy_g
   lmax = Scatter(ia)%l_max
   lamax = Scatter(ia)%lamax
   lmax_pot = Scatter(ia)%lmax_pot

!   open (101,file="kyk_ss",action="write")
!   write (101,*) kyk_wks(:,1,:,2)
!   close(101)

   call calNum_kyk(lmax_pot,lamax,kyk_wks(:,1,:,:),n_kyk,ind_kyk)
   call calNum_kyk(lmax_pot,lamax,kyk_wks(:,2,:,:),n_kyk_bar,ind_kyk_bar)

   if (is_Bxyz_max) then
      do i=1,3
         call calNum_kyk(lmax_pot,lamax,kyk_Bxyz_wks(:,1,:,:,i),n_kyk_Bxyz(:,:,i),ind_kyk_Bxyz(:,:,:,i))
         call calNum_kyk(lmax_pot,lamax,kyk_Bxyz_wks(:,2,:,:,i),n_kyk_bar_Bxyz(:,:,i),ind_kyk_bar_Bxyz(:,:,:,i))
      enddo
   endif
!      print*,"SS n_kyk="
!      print*, n_kyk(2,2)
!      print*,"ind_kyk="
!      print*, ind_kyk(:,2,2)

   if (.not. allocated(wks_green)) then
      sz = 0
      do id = 1, LocalNumAtoms
         sz = sz + Scatter(id)%numrs_cs*Scatter(id)%kmax_green
      enddo
      allocate( wks_green(sz))
      loc1 = 0
      do id = 1, LocalNumAtoms
         nofr = Scatter(id)%numrs_cs
         km = Scatter(id)%kmax_green !kmax_green=(lgreen+1)**2
         loc2 = loc1 + nofr*km
         p1c => wks_green(loc1+1:loc2)
         Scatter(id)%green => aliasArray2_c(p1c,nofr,km)
         loc1 = loc2
      enddo
   endif
   if (.not. allocated(wks_green_Bxyz) .and. is_Bxyz_max) then
      sz = 0
      do id = 1, LocalNumAtoms
         sz = sz + Scatter(id)%numrs_cs*Scatter(id)%kmax_green*3
      enddo
      allocate( wks_green_Bxyz(sz))
      loc1 = 0
      do id = 1, LocalNumAtoms
         nofr = Scatter(id)%numrs_cs
         km = Scatter(id)%kmax_green !kmax_green=(lgreen+1)**2
         loc2 = loc1 + nofr*km*3
         p1c => wks_green_Bxyz(loc1+1:loc2)
         Scatter(id)%green_Bxyz => aliasArray3_c(p1c,nofr,km,3)
         loc1 = loc2
      enddo
   endif

   nr = Scatter(ia)%numrs_cs !nr is a global variable
   green => Scatter(ia)%green
   if (is_Bxyz_max) then
      green_Bxyz => Scatter(ia)%green_Bxyz
      green_Bxyz = CZERO
   endif
   green=CZERO
   Grid => getGrid(ia)
   tmpI=CZERO
   do ir=1, nr
   !note divided by pi in calculating DOS or density is not included in prefactor
      prefactor=kappa*(energy+Me*LightSpeed**2)/LightSpeed**2
      prefactor_bar=kappa*(energy+Me*LightSpeed**2)/LightSpeed**2*&
                    (kappa**2*LightSpeed**2)/(energy+Me*LightSpeed**2)**2
      bjtmp=bjl(:,ir)
      bntmp=bnl(:,ir)
      call get_midM(lamax,ss_mat_r(:,:,ir),&
           cc_mat_r(:,:,ir),OmegaHat,mid_matrix,d_ss_mat_r(:,:,ir),&
           d_cc_mat_r(:,:,ir)) !mid_matrix has been initialized in get_midM
      do lp=1,lamax
   !      kap_lp=kapdex(lp)
   !      my_lp=mydex(lp)
         l_lp=ldex(lp)
         lbar_lp=ldex_bar(lp)
         do lpp=1,lamax
   !         kap_lpp=kapdex(lpp)
   !         my_lpp=mydex(lpp)
            l_lpp=ldex(lpp)
            lbar_lpp=ldex_bar(lpp)
!
            tmpI=bntmp(l_lpp)*mid_matrix(lpp,lp,1)*bntmp(l_lp)
!if(ir==nr) then
!print*,"SS tmpI 1", tmpI, kyk_wks(kl1,1,lpp,lp), n_kyk(lpp,lp)
!endif
            tmpI=tmpI - bntmp(l_lpp)*mid_matrix(lpp,lp,2)*bjtmp(l_lp)
!if(ir==nr) then
!print*,"SS tmpI 2", tmpI
!endif
            tmpI=tmpI - bjtmp(l_lpp)*mid_matrix(lpp,lp,3)*bntmp(l_lp)
!if(ir==nr) then
!print*,"SS tmpI 3", tmpI
!endif
            tmpI=tmpI + bjtmp(l_lpp)*mid_matrix(lpp,lp,4)*bjtmp(l_lp)
!if(ir==nr) then
!print*,"SS tmpI 4", tmpI
!endif
            tmpI=tmpI*prefactor
            do i=1,n_kyk(lpp,lp)
               kl1=ind_kyk(i,lpp,lp)
               green(ir,kl1)=green(ir,kl1) + tmpI*kyk_wks(kl1,1,lpp,lp)
            enddo
            if (is_Bxyz_max) then
               do i=1,n_kyk_Bxyz(lpp,lp,2) !note change form zxy of n_kyk_Bxyz to xyz of green_Bxyz
                  kl1=ind_kyk_Bxyz(i,lpp,lp,2)
                  green_Bxyz(ir,kl1,1)=green_Bxyz(ir,kl1,1) + tmpI*kyk_Bxyz_wks(kl1,1,lpp,lp,2)
               enddo
               do i=1,n_kyk_Bxyz(lpp,lp,3) 
                  kl1=ind_kyk_Bxyz(i,lpp,lp,3)
                  green_Bxyz(ir,kl1,2)=green_Bxyz(ir,kl1,2) + tmpI*kyk_Bxyz_wks(kl1,1,lpp,lp,3)
               enddo
               do i=1,n_kyk_Bxyz(lpp,lp,1) 
                  kl1=ind_kyk_Bxyz(i,lpp,lp,1)
                  green_Bxyz(ir,kl1,3)=green_Bxyz(ir,kl1,3) + tmpI*kyk_Bxyz_wks(kl1,1,lpp,lp,1)
               enddo
            endif
!if(ir==nr) then
!print*,"green(1,1)", green(1,1)
!endif
   ! the lower "bar" part
            tmpI=bntmp(lbar_lpp)*mid_matrix(lpp,lp,1)*bntmp(lbar_lp)
            tmpI=tmpI - bntmp(lbar_lpp)*mid_matrix(lpp,lp,2)*bjtmp(lbar_lp)
            tmpI=tmpI - bjtmp(lbar_lpp)*mid_matrix(lpp,lp,3)*bntmp(lbar_lp)
            tmpI=tmpI + bjtmp(lbar_lpp)*mid_matrix(lpp,lp,4)*bjtmp(lbar_lp)
            tmpI=tmpI*prefactor_bar
            do i=1,n_kyk_bar(lpp,lp)
               kl1=ind_kyk_bar(i,lpp,lp)
               green(ir,kl1)=green(ir,kl1) + tmpI*kyk_wks(kl1,2,lpp,lp)
            enddo
            if (is_Bxyz_max) then
               do i=1,n_kyk_bar_Bxyz(lpp,lp,2) !note change form zxy of nkyk_Bxyz to xyz of green_Bxyz
                  kl1=ind_kyk_bar_Bxyz(i,lpp,lp,2)
                  green_Bxyz(ir,kl1,1)=green_Bxyz(ir,kl1,1) + tmpI*kyk_Bxyz_wks(kl1,2,lpp,lp,2)
               enddo
               do i=1,n_kyk_bar_Bxyz(lpp,lp,3) 
                  kl1=ind_kyk_bar_Bxyz(i,lpp,lp,3)
                  green_Bxyz(ir,kl1,2)=green_Bxyz(ir,kl1,2) + tmpI*kyk_Bxyz_wks(kl1,2,lpp,lp,3)
               enddo
               do i=1,n_kyk_bar_Bxyz(lpp,lp,1) 
                  kl1=ind_kyk_bar_Bxyz(i,lpp,lp,1)
                  green_Bxyz(ir,kl1,3)=green_Bxyz(ir,kl1,3) + tmpI*kyk_Bxyz_wks(kl1,2,lpp,lp,1)
               enddo
            endif
         enddo
      enddo
      green(ir,:) = green(ir,:)*Grid%r_mesh(ir)**2
      if (is_Bxyz_max) then
         green_Bxyz(ir,:,:) = green_Bxyz(ir,:,:)*Grid%r_mesh(ir)**2
      endif
   enddo

!   open (101,file="mid_M",action="write")
!   write (101,*) green(:,1)!mid_matrix(:,:,1)
!   close(101)
!print*,"before computeRelDOS"
   call computeRelDOS(ia)

   if (is_Bxyz_max) then
      deallocate( n_kyk_Bxyz )
      deallocate( ind_kyk_Bxyz )
      deallocate( n_kyk_bar_Bxyz )
      deallocate( ind_kyk_bar_Bxyz )
   endif
   end subroutine computeRelGF

   !=================================================================
   subroutine get_midM(lamax,ss_mat_r,cc_mat_r,OmegaHat,mid_matrix,&
              d_ss_mat_r,d_cc_mat_r)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : CZERO, CONE, SQRTm1

   implicit none

   integer (kind=IntKind),intent(in) :: lamax
   integer (kind=IntKind) :: i
   complex (kind=CmplxKind) :: OmegaHat(lamax,lamax)
   complex (kind=CmplxKind) :: ss_mat_r(lamax,lamax),cc_mat_r(lamax,lamax),Omat(lamax,lamax)
   complex (kind=CmplxKind) :: d_ss_mat_r(lamax,lamax),d_cc_mat_r(lamax,lamax)
   complex (kind=CmplxKind) :: mid_matrix(lamax,lamax,4)
   complex (kind=CmplxKind) :: mat_tmp(lamax,lamax)

   mat_tmp=OmegaHat
   !mat_tmp=(S^t*(i*S-C))^-1

   Omat=(0.d0,0.d0)
   call zgemm('n','n',lamax,lamax,lamax,cone,&
   ss_mat_r, lamax, mat_tmp, lamax, czero, Omat, lamax)
   !Omat=ss_r*OmegaHat

   mid_matrix(:,:,1)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_ss_mat_r, lamax, czero, mid_matrix(:,:,1), lamax)
   !mid_matrix(1)=ss_r*OmegaHat*d_ss_r

   mid_matrix(:,:,2)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_cc_mat_r, lamax, czero, mid_matrix(:,:,2), lamax)
   !mid_matrix(2)=ss_r*OmegaHat*d_cc_r

   Omat=(0.d0,0.d0)
   call zgemm('n','n',lamax,lamax,lamax,cone,&
   cc_mat_r, lamax, mat_tmp, lamax, czero, Omat, lamax)
   !Omat=cc_r*OmegaHat

   mid_matrix(:,:,3)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_ss_mat_r, lamax, czero, mid_matrix(:,:,3), lamax)
   !mid_matrix(3)=cc_r*OmegaHat*d_ss_r

   mid_matrix(:,:,4)=(0.d0,0.d0)
   call zgemm('n','t',lamax,lamax,lamax,cone,&
   Omat, lamax, d_cc_mat_r, lamax, czero, mid_matrix(:,:,4), lamax)
   !mid_matrix(4)=cc_r*OmegaHat*d_cc_r

   end subroutine get_midM


   ! This function has the same interface as the non-relativistic case
   !=================================================================
   function getRelScatteringMatrix(sm_type,spin,site,atom,dsize) result(sm)
   !=================================================================

   implicit none
 
   character (len=*), intent(in) :: sm_type

   integer (kind=IntKind), intent(in), optional :: spin, site, atom
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: id

   complex (kind=CmplxKind), pointer :: sm(:,:)

   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface

   if (.not.Initialized) then
      call ErrorHandler('getRelScatteringMatrix', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelScatteringMatrix',                    &
                        'invalid number of local atoms', LocalNumAtoms)
   endif

!  This needs to be fixed for the CPA case
   if (present(site) .and. present(atom)) then
      id = max(site,atom) ! this is just a temporary fix
   else if (present(site)) then
      id = site
   else if (present(atom)) then
      id = atom
   else
      id = 1
   endif
!
   if (nocaseCompare(sm_type,'T-Matrix')) then
      sm => Scatter(id)%t_mat
   else if (nocaseCompare(sm_type,'S-Matrix')) then
      sm => Scatter(id)%S_mat
   else if (nocaseCompare(sm_type,'Sine-Matrix')) then
      sm => Scatter(id)%sin_mat
   else if (nocaseCompare(sm_type,'Cosine-Matrix')) then
      sm => Scatter(id)%cos_mat
   else if (nocaseCompare(sm_type,'Jost-Matrix')) then
      sm => Scatter(id)%jost_mat
   else if (nocaseCompare(sm_type,'JostInv-Matrix')) then
      sm => Scatter(id)%jinv_mat
!  else if (nocaseCompare(sm_type,'Omega-Matrix')) then
!     sm => Scatter(id)%Omega_mat
   else if (nocaseCompare(sm_type,'OmegaHat-Matrix')) then
      sm => Scatter(id)%OmegaHat_mat
   else if (nocaseCompare(sm_type,'OmegaHatInv-Matrix')) then
      sm => Scatter(id)%OmegaHatInv_mat
   else
      call ErrorHandler('getRelScatteringMatrix',                    &
                        'Scattering matrix type is invalid',sm_type)
   endif
!
   if (present(dsize)) then
      dsize = size(sm,1)
   endif

   end function getRelScatteringMatrix

   !=================================================================
   function getRelTMatrix(atom) result(t_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: t_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelTMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelTMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   t_mat => Scatter(atom)%t_mat

   end function getRelTMatrix

   !=================================================================
   function getRelSineMatrix(atom) result(sin_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: sin_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelSinMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSinMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   sin_mat => Scatter(atom)%sin_mat

   end function getRelSineMatrix

   !=================================================================
   function getRelCosineMatrix(atom) result(cos_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: cos_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelSinMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSinMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   cos_mat => Scatter(atom)%cos_mat

   end function getRelCosineMatrix

   !=================================================================
   function getRelSMatrix(atom) result(S_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: S_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelSinMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSinMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   S_mat => Scatter(atom)%S_mat

   end function getRelSMatrix

   !=================================================================
   function getRelJostMatrix(atom) result(jost_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: jost_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelSinMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSinMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   jost_mat => Scatter(atom)%jost_mat

   end function getRelJostMatrix

   !=================================================================
   function getRelJostInvMatrix(atom) result(jinv_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: jinv_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelSinMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSinMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   jinv_mat => Scatter(atom)%jinv_mat

   end function getRelJostInvMatrix

   !=================================================================
   function getRelOmegaHatMatrix(atom) result(OmegaHat_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: OmegaHat_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelSinMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSinMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   OmegaHat_mat => Scatter(atom)%OmegaHat_mat

   end function getRelOmegaHatMatrix

   !=================================================================
   function getRelOmegaHatInvMatrix(atom) result(OmegaHatInv_mat)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom

   complex (kind=CmplxKind), pointer :: OmegaHatInv_mat(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelSinMat', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSinMat', 'invalid number of local atoms', LocalNumAtoms)
   endif


   OmegaHatInv_mat => Scatter(atom)%OmegaHatInv_mat

   end function getRelOmegaHatInvMatrix

!  ===================================================================
   function getRelSolutionRmeshSize(atom) result(nr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: atom
   integer (kind=IntKind) :: nr
!
   if (.not.Initialized) then
      call ErrorHandler('getRelSolutionRmeshSize','module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelSolutionRmeshSize','invalid number of local atoms', &
                        LocalNumAtoms)
   endif
!
   nr = Scatter(atom)%numrs
!
   end function getRelSolutionRmeshSize

!======================================================================
   function getRelRegSolutionSine(atom) result(ss_mat_r)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom
   complex (kind=CmplxKind), pointer :: ss_mat_r(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelRegSolutionSine', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelRegSolutionSine', 'invalid number of local atoms', &
                        LocalNumAtoms)
   endif

   ss_mat_r => Scatter(atom)%ss_mat_r

   end function getRelRegSolutionSine

!======================================================================
   function getRelRegSolutionCosine(atom) result(cc_mat_r)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: atom
   complex (kind=CmplxKind), pointer :: cc_mat_r(:,:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelRegSolutionSine', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelRegSolutionSine', 'invalid number of local atoms', LocalNumAtoms)
   endif

   cc_mat_r => Scatter(atom)%cc_mat_r

   end function getRelRegSolutionCosine

!======================================================================
   function getRelkyk(atom) result(kyk_loc)
!======================================================================

   implicit none
   integer (kind=IntKind), intent(in) :: atom
   real (kind=RealKind), pointer :: kyk_loc(:,:,:,:)

   if (.not. Initialized) then
      call ErrorHandler('getRelkyk_Bz','module not initialized')
   else
      kyk_loc => Scatter(atom)%kyk
   endif
   end function getRelkyk

!======================================================================
   function getRelkyk_Bz(atom) result(kyk_loc)
!======================================================================

   implicit none
   integer (kind=IntKind), intent(in) :: atom
   real (kind=RealKind), pointer :: kyk_loc(:,:,:,:)

   if (.not. Initialized) then
      call ErrorHandler('getRelkyk_Bz','module not initialized')
   else
      if (is_Bxyz_max) then
         kyk_loc => Scatter(atom)%kyk_Bz
      else
         call ErrorHandler('getRelkyk_Bz','is_Bxyz is NOT .true. and kyk_Bxyz not calculated')
      endif
   endif
   end function getRelkyk_Bz

!======================================================================
   function getRelkyk_Bx(atom) result(kyk_loc)
!======================================================================

   implicit none
   integer (kind=IntKind), intent(in) :: atom
   real (kind=RealKind), pointer :: kyk_loc(:,:,:,:)

   if (.not. Initialized) then
      call ErrorHandler('getRelkyk_Bx','module not initialized')
   else
      if (is_Bxyz_max) then
         kyk_loc => Scatter(atom)%kyk_Bx
      else
         call ErrorHandler('getRelkyk_Bx','is_Bxyz is NOT .true. and kyk_Bxyz not calculated')
      endif
   endif
   end function getRelkyk_Bx

!======================================================================
   function getRelkyk_By(atom) result(kyk_loc)
!======================================================================

   implicit none
   integer (kind=IntKind), intent(in) :: atom
   real (kind=RealKind), pointer :: kyk_loc(:,:,:,:)

   if (.not. Initialized) then
      call ErrorHandler('getRelkyk_By','module not initialized')
   else
      if (is_Bxyz_max) then
         kyk_loc => Scatter(atom)%kyk_By
      else
         call ErrorHandler('getRelkyk_By','is_Bxyz is NOT .true. and kyk_Bxyz not calculated')
      endif
   endif
   end function getRelkyk_By

!======================================================================
   function getkapdex() result(kapdex_loc)
!======================================================================
   implicit none
   integer (kind=IntKind),pointer :: kapdex_loc(:)
   if (.not. Initialized) then
      call ErrorHandler('getkapdex','module not initialized')
   else
      kapdex_loc => kapdex   
   endif
   end function getkapdex
   

!======================================================================
   function getldex() result(ldex_loc)
!======================================================================
   implicit none
   integer (kind=IntKind),pointer :: ldex_loc(:)
   if (.not. Initialized) then
      call ErrorHandler('getldex','module not initialized')
   else
      ldex_loc => ldex   
   endif
   end function getldex
   
!======================================================================
   function getjdex() result(jdex_loc)
!======================================================================
   implicit none
   integer (kind=IntKind),pointer :: jdex_loc(:)
   if (.not. Initialized) then
      call ErrorHandler('getjdex','module not initialized')
   else
      jdex_loc => jdex   
   endif
   end function getjdex

!======================================================================
   function getmydex() result(mydex_loc)
!======================================================================
   implicit none
   integer (kind=IntKind),pointer :: mydex_loc(:)
   if (.not. Initialized) then
      call ErrorHandler('getldex','module not initialized')
   else
      mydex_loc => mydex   
   endif
   end function getmydex

!======================================================================
   function getldex_bar() result(ldex_bar_loc)
!======================================================================
   implicit none
   integer (kind=IntKind),pointer :: ldex_bar_loc(:)
   if (.not. Initialized) then
      call ErrorHandler('getldex','module not initialized')
   else
      ldex_bar_loc => ldex_bar   
   endif
   end function getldex_bar
   
   !=================================================================
   function getRelSSDOS(is, atom) result(dos_loc)
   !=================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: is, atom
   !is from 1:4, four parts of green function or DOS

   complex (kind=CmplxKind), pointer :: dos_loc(:,:)

   if (.not.Initialized) then
      call ErrorHandler('getRelDOS', 'module not initialized')
   else if (atom < 1 .or. atom > LocalNumAtoms) then
      call ErrorHandler('getRelDOS', 'invalid number of local atoms', LocalNumAtoms)
   endif

   if ( is >1 .and. (.not.is_Bxyz_max) ) then
      call ErrorHandler('getRelSSDOS','is > 1 in non-magnetic calculation, invalid case')
   endif

   if (is == 1) then
         dos_loc => Scatter(atom)%dos
   else if (is == 2) then
         dos_loc => Scatter(atom)%dos_Bxyz(:,:,1)
   else if (is == 3) then
         dos_loc => Scatter(atom)%dos_Bxyz(:,:,2)
   else if (is == 4) then
         dos_loc => Scatter(atom)%dos_Bxyz(:,:,3)
   else
      call ErrorHandler('getRelSSDOS','is is not between 1 and 4')
   endif

   end function getRelSSDOS

        !=================================================================
        function getRelSSGreen(is, atom) result(dos_loc)
        !=================================================================
   !note that Scatter(atom)%dos has been changed
        implicit none

        integer (kind=IntKind), intent(in) :: is, atom
        !is from 1:4, four parts of green function or DOS
        integer (kind=IntKind) :: lmax_dos, jmax_dos, kmax_kkr_loc
        integer (kind=IntKind) :: ir, jl, kl, klc, l, m, kl1, kl1p, kl2, kl2c, kl2p,&
          kl2pc, m1, m2, m2p, ma, i, l2

        complex (kind=CmplxKind), pointer :: dos_loc(:,:)

        if (.not.Initialized) then
                call ErrorHandler('getRelGreen', 'module not initialized')
        else if (atom < 1 .or. atom > LocalNumAtoms) then
                call ErrorHandler('getRelGreen', 'invalid number of local atoms', LocalNumAtoms)
        endif

        if ( is >1 .and. (.not.is_Bxyz_max) ) then
                call ErrorHandler('getRelSSGreen','is > 1 in non-magnetic calculation, invalid case')
        endif

        lmax_dos = Scatter(atom)%lmax_green
        jmax_dos = (lmax_dos+1)*(lmax_dos+2)/2
        nr = Scatter(atom)%numrs_cs
!        do jl = 1, jmax_dos
!           kl = kofj(jl)
!                m = mofj(jl)
!                klc = kl-2*m
!                do ir = 1, nr
!                   Scatter(atom)%dos(ir,jl)=Scatter(atom)%green(ir,kl)/PI
!                enddo
!        enddo
!        if (is_Bxyz_max) then
!           kl = kofj(jl)
!                m = mofj(jl)
!                klc = kl-2*m
!                do i = 1, 3
!                   do ir = 1, nr
!                      Scatter(atom)%dos_Bxyz(ir,jl,i)=Scatter(atom)%green_Bxyz(ir,kl,i)/PI
!                   enddo
!                enddo
!        endif

        if (is == 1) then
                dos_loc => Scatter(atom)%green
        else if (is == 2) then
                dos_loc => Scatter(atom)%green_Bxyz(:,:,1)
        else if (is == 3) then
                dos_loc => Scatter(atom)%green_Bxyz(:,:,2)
        else if (is == 4) then
                dos_loc => Scatter(atom)%green_Bxyz(:,:,3)
        else
                call ErrorHandler('getRelSSDOS','is is not between 1 and 4')
        endif

        end function getRelSSGreen

   !=================================================================
   function getRelSSGreenOut(ia) result(dos)
   !=================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use RadialGridModule, only : getGrid
   implicit none

   integer (kind=IntKind), intent(in) :: ia
   complex (kind=CmplxKind) :: dos
   integer (kind=IntKind) :: lmax,lamax,lam,l
   complex (kind=CmplxKind) :: fint(0:Scatter(ia)%l_max)

   if (.not. Initialized) then
      call ErrorHandler('getRelSSGreenOut','module not initialized')
   else if (ia < 1 .or. ia > LocalNumAtoms) then
      call ErrorHandler('getRelGreenOut', 'invalid number of local atoms', LocalNumAtoms)
   endif

   nr=Scatter(ia)%numrs_cs
   lmax=Scatter(ia)%l_max
   lamax=Scatter(ia)%lamax
   Grid => getGrid(ia)

   call calIntSphHankelSq0(lmax,Grid%r_mesh(nr),kappa_g**2,fint)
   dos=CZERO
   do lam=1,lamax
      l = ldex(lam)
      dos=dos-2.d0*energy_g/(LightSpeed**2)*kappa_g*kappa_g*&
          Scatter(ia)%t_mat(lam,lam)*fint(l)/PI
   enddo
   end function getRelSSGreenOut

!  ===================================================================
   subroutine calIntSphHankelSq0(lmax,rc,energy,fint)
!  ===================================================================
      use KindParamModule, only : IntKind, RealKind, CmplxKind
      use MathParamModule, only : Half, SQRTm1
      use BesselModule, only : SphericalBessel, SphericalNeumann
!
      implicit none
!
      integer (kind=IntKind), intent(in) :: lmax
      integer (kind=IntKind) :: l
!
      real (kind=RealKind), intent(in) :: rc
!
      complex (kind=CmplxKind), intent(in) :: energy
      complex (kind=CmplxKind), intent(out) :: fint(0:lmax)
      complex (kind=CmplxKind) :: bjl(0:lmax+1), bnl(0:lmax+1), x, bhl, bhlp1, bhlm1
      complex (kind=CmplxKind) :: Ul, Ulp1
!
      x = rc*sqrt(energy) !+SQRTm1*0.00001d0 !Changed by Xianglin
!     ----------------------------------------------------------------
      call SphericalBessel(lmax+1,x,bjl)
      call SphericalNeumann(lmax+1,x,bnl)
!     ----------------------------------------------------------------
      bhl = bjl(0)+SQRTm1*bnl(0)
      bhlp1 = bjl(1)+SQRTm1*bnl(1)
      fint(0) = (bhl*bhlp1/x-bhl**2-bhlp1**2)*HALF
      do l = 1, lmax
         bhlm1 = bhl
         bhl = bhlp1
         bhlp1 = bjl(l+1)+SQRTm1*bnl(l+1)
!        fint(l) = ((2*l+1)*bhl*bhlp1/x-bhl**2-bhlp1**2)*HALF
         fint(l) = (bhlm1*bhlp1-bhl**2)*HALF
      enddo
!
!     ================================================================
!     Another way of calculating fint is to use recursive relation.
!     Both ways have shown to give the same results.
!     ================================================================
!!    Ul = -SQRTm1*exp(SQRTm1*TWO*x)/(TWO*x**3)
!!    write(6,'(a,2d16.8,2x,2d16.8)')'fint(0), Ul = ',fint(0), Ul
!!    fint(0) = Ul
!!    do l = 1, lmax
!!       bhlm1 = bjl(l-1)+SQRTm1*bnl(l-1)
!!       bhl = bjl(l)+SQRTm1*bnl(l)
!!       Ulp1 = (bhl**2+bhlm1**2+(2*l+1)*fint(l-1))/(2*l-1.0d0)
!!       write(6,'(a,2d16.8,2x,2d16.8)')'fint(l), Ulp1 = ',fint(l), Ulp1
!!       fint(l) = Ulp1
!!    enddo
!     ================================================================
!
      do l = 0, lmax
         fint(l) = fint(l)*rc**3
      enddo
!
   end subroutine calIntSphHankelSq0

end module RelSSSolverModule

