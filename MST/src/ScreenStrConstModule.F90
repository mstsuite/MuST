!  *******************************************************************
!  * Module for calculating Screened KKR structure constant matrix   *
!  * External Modules used:                                          *
!  *   KindParamModule                                               *
!  *   MathParamModule                                               *
!  *   WriteMatrixModule                                             *
!  *   ErrorHandlerModule                                            *
!  *   RSpaceStrConstModule                                          *
!  *   GauntFactorsModule                                            *
!  *   PhysParamModule                                               *
!  *   BesselModule                                                  *
!  *   SingleScatteringModule                                        *
!  *                                                                 *
!  * Public functions:                                               *
!  *                                                                 *
!  *    initScreenStrConst(Bravais,NumAtoms,AtomPosition,lmax_swiss, *
!  *                       rmt_swiss,v_swiss,rcut_swiss,iprint,istop)*
!  *    Purpose: initialize the module for calculating screened      *
!  *             KKR structure constant matrix for n atom per cell   *
!  *             Define work space size, clusters sizes and setup    *
!  *             the neighbors around the atoms in the unit cell     *
!  *    Note:    GauntFactorsModule needs to be initialized first    *
!  *    Input:   Bravais  = Bravais lattice vector (real array(3,3)) *
!  *             NumAtoms = number of atoms in the unit cell         *
!  *             AtomPosition = the position of the atoms in the unit*
!  *                    cell, a 2-d array (1:3,1:NumAtoms)           *
!  *             lmax_swiss = largest possible l-cut off for gij     *
!  *                    matrix (integer), 1-d array of NumAtoms size *
!  *             rmt_swiss = muffintin radia of the square well      *
!  *                    potentials v_swiss 1-d array of NumAtoms size*
!  *             v_swiss = the hight of the repulsive potential      *
!  *                    an 1-d array of NumAtoms size                *
!  *             rcut_swiss = radia of the screening cluster         *
!  *                    an 1-d array of NumAtoms size                *
!  *             istop = routine name to stop (character string)     *
!  *             iprint= print level (integer)                       *
!  *                                                                 *
!  *    endScreenStrConst()                                          *
!  *    Purpose: deallocate the internal allocated arrays and clear  *
!  *             the storage                                         *
!  *    Usage:   should be used whenever one or more atoms change    *
!  *             their position or before the end of the program     *
!  *                                                                 *
!  *    calScreenTau0J( myatom, epsilon, isRelativistic )            *
!  *    Purpose: Calculates the Tau0J matrices, arow of block        *
!  *             matrices of length equals to the number of neighbors*
!  *             in the screening cluster                            *
!  *    Input:   myatom = the central site, integre which cannot be  *
!  *                    larger than the NumAtoms                     *
!  *             epsilon = the energy at which the calculations are  *
!  *                    performed                                    *
!  *             isRelativistic = logical value, sets the type of    *
!  *                    calculations, false(true)=non-(relativistic) *
!  *                                                                 *
!  *    getTauBlockJ( mu, nrs )                                      *
!  *    Purpose: retrieve the nrs block of TAu0J for site mu         *
!  *    Input:   mu = Site index, cannot be largerthan NumAtoms      *
!  *             nrs = index of the neighbor in the screening cluster*
!  *                 of site mu, cannot be larger than NrSiteClu_Max,*
!  *                 wich is setup in initScreenStrConst             *
!  *    Result:  ptau0j = the nrs block of screen structure constant *
!  *                    {i, j} matrix (complex array(:,:))  block    *
!  *                    where,  1st dim = kmaxi=(lmaxi+1)**2         *
!  *                            2nd dim = kmaxj=(lmaxj+1)**2         *
!  *                                                                 *
!  *    SpinSpace( InputMatrix )                                     *
!  *    Purpose: Brings the InputMatrix to spin space                *
!  *    Input:   InputMatrix - square matrix of dim kmax_max         *
!  *    Output:  OutMatrix - square matrix of dim nSpinCant*kmax_max *
!  *                   where the Input matrix is copyed in each block*
!  *                                                                 *
!  *    getDaltaM( ns, myAtom )                                      *
!  *    Note     The SingleScatteringModule needs to be initalized   *
!  *    Purpose: Retrive the matrix {T_swissSpin**(-1)-Tmatrix**(-1)]*
!  *    Input:   ns - spin dimension                                 *
!  *             myAtom - site index(cannot be greater than NumAtoms)*
!  *    Output:  pDelta - the Delta matrix [M_swiss-M]               *
!  *                                                                 *
!  *    getScreenLmax()                                              *
!  *    Purpose: returns the maximum value of cutt-off lmax          *
!  *                                                                 *
!  *    printScreenTau0J(ns,nb)                                      *
!  *    Purpose: Prints out screen cluster info and                  *
!  *             Screen Matrices if they are allocated               *
!  *    Input: ns - spin dim                                         *
!  *           nb - number of Tau0J blocks to be printed out         *
!  *                                                                 *
!  *******************************************************************
module ScreenStrConstModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : zero, one, czero, cone
   use MathParamModule, only : two, half, third
   use MathParamModule, only : ten2m8, ten2m6, ten2m14
   use MathParamModule, only : sqrtm1, pi2, pi
   use WriteMatrixModule, only  : writeMatrix
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler,        &
                                  StopHandler
!
public :: initScreenStrConst, &
          endScreenStrConst,  &
          calScreenTauK,      &
          getTau00K,          &
          getTauK,            &
          calScreenTau0J,     &
          SpinSpace,          &
          getDeltaM,          &
          getTauBlockJ,       &
          getScreenLmax,      &
          printScreenTau0J
!
   interface initScreenStrConst
      module procedure initScreenStrConst0, initScreenStrConst1
   end interface
!
private
!
   type ScreenMatrices
      complex (kind=CmplxKind), pointer :: tauS0jp(:,:,:)    ! tau0j swiss
      complex (kind=CmplxKind), pointer :: tauS00p(:,:)      ! 00 block of tau0j
      complex (kind=CmplxKind), pointer :: deltaMp(:,:)      ! [T_swiss^(-1)-T_syst^(-1)]
      complex (kind=CmplxKind), pointer :: Tswiss(:,:)
      complex (kind=CmplxKind), pointer :: TswissSpin(:,:)   ! spin space
   end type ScreenMatrices
!
   type (ScreenMatrices), allocatable, target :: ScreenM(:)
!
   logical :: isSpinSpace = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: NumAtoms, NumUniqPots, nSpinCant
   integer (kind=IntKind) :: lmax_max, kmax_max
   integer (kind=IntKind) :: NrClu_Max, NrSiteClu_Max, NrSiteMuNu_Max
   integer (kind=IntKind), allocatable :: UniqPots(:)
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: UniqV(:)
   real (kind=RealKind), allocatable :: UniqRmt(:)
!
   complex (kind=CmplxKind) :: energy
!
   type NeighborScreenStruct
      integer (kind=IntKind) :: NumLocalSites
      integer (kind=IntKind), pointer :: AtomType(:)
      integer (kind=IntKind), pointer :: LmaxClu(:)
      integer (kind=IntKind), pointer :: IndexIJK(:,:)
      real (kind=RealKind), pointer :: Position(:,:)
   end type NeighborScreenStruct
!
   type ScreenStrStruct
      integer (kind=IntKind) :: AtomType
      integer (kind=IntKind) :: NrMat
      integer (kind=IntKind) :: Lmax
      real (kind=RealKind) :: Position(3)
      real (kind=RealKind) :: Rmt
      real (kind=RealKind) :: Rcut
      real (kind=RealKind) :: V
      integer (kind=IntKind), pointer :: NrNu(:)
      real (kind=RealKind), pointer   :: R_MuNu(:,:,:)
      type (NeighborScreenStruct)     :: Neighbors
   end type ScreenStrStruct
!
   type (ScreenStrStruct), allocatable :: ScreenSite(:)
!
   complex (kind=CmplxKind), allocatable, target :: TauK_Big(:,:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initScreenStrConst0(br,na,ns,ap,lmax_swiss,rmt_swiss,       &
                                  v_swiss,rcut_swiss,iprint,istop)
!  ===================================================================
   use RSpaceStrConstModule, only : initRSpaceStrConst
   use GauntFactorsModule, only : initGauntFactors
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: na, ns
   integer (kind=IntKind), intent(in) :: lmax_swiss
   integer (kind=IntKind), intent(in) :: iprint
!
   real (kind=RealKind), intent(in) :: br(3,3)
   real (kind=RealKind), intent(in) :: ap(3,na)
   real (kind=RealKind), intent(in) :: rmt_swiss
   real (kind=RealKind), intent(in) :: v_swiss
   real (kind=RealKind), intent(in) :: rcut_swiss
!
   integer (kind=IntKind) :: i
!
   print_level = iprint
   stop_routine = istop
!
   NumAtoms = na
   bravais = br
   nSpinCant = ns
   if ( nSpinCant == 2 ) then
      isSpinSpace = .true.
   endif
!
   allocate ( ScreenSite(NumAtoms) )
   do i=1,NumAtoms
      ScreenSite(i)%Lmax = lmax_swiss
      ScreenSite(i)%Position(1:3) = ap(1:3,i)
      ScreenSite(i)%Rmt = rmt_swiss
      ScreenSite(i)%V = v_swiss
      ScreenSite(i)%Rcut = rcut_swiss
      ScreenSite(i)%AtomType = 1
   enddo
!  ===================================================================
!  find the maximun dimension of real space cluster
!  ===================================================================
   call getClustDim()
!  ===================================================================
!  locates the neighbors in the cluster associated to each atom in the
!  system
!  ===================================================================
   do i = 1,NumAtoms
      call neighbors(rcut_swiss,i)
   enddo
!
   energy = -(100000.d0,0.d0)
!
   lmax_max = 0
   do i = 1,NumAtoms
      lmax_max = max( lmax_max,ScreenSite(i)%Lmax )
   enddo
   kmax_max = (lmax_max+1)**2
!
   NumUniqPots = 1
   allocate(UniqPots(NumAtoms), UniqV(1), UniqRmt(1))
   UniqPots(1:NumAtoms) = 1
   UniqV(1) = v_swiss
   UniqRmt(1) = rmt_swiss
!
   allocate( ScreenM(NumUniqPots) )
!
   call initRSpaceStrConst(2*lmax_max,stop_routine,print_level)
   call initGauntFactors(lmax_max,stop_routine,print_level)
!
   end subroutine initScreenStrConst0
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initScreenStrConst1(br,na,ns,ap,lmax_swiss,rmt_swiss,       &
                                  v_swiss,rcut_swiss,iprint,istop)
!  ===================================================================
   use RSpaceStrConstModule, only : initRSpaceStrConst
   use GauntFactorsModule, only : initGauntFactors
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: na, ns
   integer (kind=IntKind), intent(in) :: lmax_swiss(na)
   integer (kind=IntKind), intent(in) :: iprint
!
   real (kind=RealKind), intent(in) :: br(3,3)
   real (kind=RealKind), intent(in) :: ap(3,na)
   real (kind=RealKind), intent(in) :: rmt_swiss(na)
   real (kind=RealKind), intent(in) :: v_swiss(na)
   real (kind=RealKind), intent(in) :: rcut_swiss(na)
!
   integer (kind=IntKind) :: i,j
!
   print_level = iprint
   stop_routine = istop
!
   NumAtoms = na
   bravais = br
   nSpinCant = ns
   if ( nSpinCant == 2 ) then
      isSpinSpace = .true.
   endif
!
   allocate ( ScreenSite(NumAtoms) )
   do i=1,NumAtoms
      ScreenSite(i)%Lmax = lmax_swiss(i)
      ScreenSite(i)%Position(1:3) = ap(1:3,i)
      ScreenSite(i)%Rmt = rmt_swiss(i)
      ScreenSite(i)%V = v_swiss(i)
      ScreenSite(i)%Rcut = rcut_swiss(i)
      ScreenSite(i)%AtomType = i
      allocate( ScreenSite(i)%NrNu(NumAtoms))
   enddo
!
   allocate(UniqPots(NumAtoms), UniqV(NumAtoms), UniqRmt(NumAtoms))
   NumUniqPots = 1
   ScreenSite(1)%AtomType = 1
   UniqPots(1) = NumUniqPots
   UniqV(1) = ScreenSite(1)%V
   UniqRmt(1) = ScreenSite(1)%Rmt
   LOOP_i: do i = 2,NumAtoms
      do j = 1,NumUniqPots
         if ( ScreenSite(i)%V==UniqV(j) .and.                           &
              ScreenSite(i)%Rmt==UniqRmt(j) ) then
            UniqPots(i)=j
            ScreenSite(i)%AtomType = j
            cycle LOOP_i
         endif
      enddo
      NumUniqPots = NumUniqPots+1
      UniqPots(i) = NumUniqPots
      UniqV(NumUniqPots) = ScreenSite(i)%V
      UniqRmt(NumUniqPots) = ScreenSite(i)%Rmt
   enddo LOOP_i
!
!  ===================================================================
!  finds the maximun dimension of the real space clusters
!  ===================================================================
   call getClustDim()
!  ===================================================================
!  locates the neighbors in the cluster associated to each atom
!
   do i = 1,NumAtoms
      call neighbors(rcut_swiss(i),i)
   enddo
!
   energy = -100000.d0
!
   lmax_max = 0
   do i = 1,NumAtoms
      lmax_max = max( lmax_max,ScreenSite(i)%Lmax )
   enddo
   kmax_max = (lmax_max+1)**2
!
   allocate( ScreenM(NumUniqPots) )
!
   call initRSpaceStrConst(2*lmax_max,stop_routine,print_level)
   call initGauntFactors(lmax_max,stop_routine,print_level)
!
   end subroutine initScreenStrConst1
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endScreenStrConst()
!  ===================================================================
   use RSpaceStrConstModule, only : endRSpaceStrConst
!
   implicit none
!
   integer (kind=IntKind) :: i
!
   do i=1,NumAtoms
      if (ScreenSite(i)%Neighbors%NumLocalSites > 0) then
         deallocate( ScreenSite(i)%Neighbors%LmaxClu)
         deallocate( ScreenSite(i)%Neighbors%AtomType )
         deallocate( ScreenSite(i)%Neighbors%Position )
         deallocate( ScreenSite(i)%Neighbors%IndexIJK )
      endif
!     deallocate( ScreenSite(i)%Neighbors )
   enddo
   deallocate( ScreenSite )
!
   deallocate( UniqPots, UniqV, UniqRmt )
!
   do i = 1,NumUniqPots
      if ( associated( ScreenM(i)%tauS0jp ) ) then
           call WarningHandler('endScreenStrConst',                    &
                        'Call delScreenMatrices() to clean up memory')
           call delScreenMatrices(i)
      endif
      if ( associated( ScreenM(i)%deltaMp ) ) then
         call WarningHandler('endScreenStrConst',                      &
                          'Call delDeltaM() to clean up memory')
         call delDeltaM(i)
      endif
      if ( associated( ScreenM(i)%Tswiss ) ) then
         deallocate( ScreenM(i)%Tswiss )
         nullify( ScreenM(i)%Tswiss )
      endif
   enddo
   deallocate( ScreenM, TauK_big )
!
   call endRSpaceStrConst()
!
   end subroutine endScreenStrConst
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delScreenMatrices(i)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!
   if ( associated( ScreenM(i)%tauS0jp ) ) then
      deallocate( ScreenM(i)%tauS0jp )
      nullify( ScreenM(i)%tauS0jp )
   endif
!
   if ( associated( ScreenM(i)%tauS00p ) ) then
      deallocate( ScreenM(i)%tauS00p )
      nullify( ScreenM(i)%tauS00p )
   endif
!
   end subroutine delScreenMatrices
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delDeltaM(i)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!

   if ( associated( ScreenM(i)%deltaMp ) ) then
      deallocate( ScreenM(i)%deltaMp )
      nullify( ScreenM(i)%deltaMp )
   endif
!
   end subroutine delDeltaM
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calScreenTauK( kvec )
!  ===================================================================
!
   implicit none
!
   real (kind=RealKind), intent(in) :: kvec(3)
!   complex (kind=CmplxKind), intent(in) :: epsilon
!
   integer (kind=IntKind) :: i, j, i_ind, j_ind
   integer (kind=IntKind) :: mu, nu, mu_ind, nu_ind
   integer (kind=IntKind) :: mdim, ndim, kmax
!
   complex (kind=CmplxKind), pointer :: Tau00K(:,:)
   complex (kind=CmplxKind), pointer :: TauK_g(:,:)
   complex (kind=CmplxKind), allocatable :: TauSW_Big(:,:)
   complex (kind=CmplxKind), allocatable :: M_Big(:,:)
   complex (kind=CmplxKind), pointer :: DeltaM_mu(:,:)
   complex (kind=CmplxKind), allocatable :: Work_g(:,:)
!
   kmax = kmax_max
   if ( isSpinSpace ) then
      mdim = 2*kmax_max
   else
      mdim = kmax_max
   endif
   ndim = mdim*NumAtoms
!
   allocate( TauK_g(mdim,mdim),TauSW_Big(ndim,ndim),Work_g(mdim,mdim), &
             M_Big(ndim,ndim), DeltaM_mu(mdim,mdim) )
!
   TauK_g = czero
   M_Big  = czero
   do i = 1,ndim
      M_Big(i,i) = cone
   enddo
   TauSW_Big = czero
!
   mu_ind = 0
!
   MuLoop: do mu = 1,NumAtoms
!
!     ================================================================
!     obtain delta_m for current sublattice
!     ================================================================
      DeltaM_mu = czero
      if ( .not.associated(ScreenM(UniqPots(mu))%deltaMp) ) then
         call calDeltaM(mu)
      endif
      DeltaM_mu = ScreenM(UniqPots(mu))%deltaMp
!
      nu_ind = 0
      NuLoop: do nu = 1,NumAtoms
!
         TauK_g => getTauK_SpinMuNu( mu, nu, kvec )
!
!        =============================================================
!        calculate TauK*DeltaM
!        =============================================================
!        -------------------------------------------------------------
         call zgemm( 'n','n',mdim,mdim,mdim,cone,TauK_g,mdim,          &
                     DeltaM_mu,mdim,czero,Work_g,mdim )
!        -------------------------------------------------------------
!        =============================================================
!        add [tau_k*Delta_m] in correct place in m_big
!        Atom on my node is located in the 1,1 block
!        =============================================================
!
         do i = 1,mdim
            i_ind = nu_ind + i
            do j = 1,mdim
               j_ind = mu_ind + j
               M_Big(i_ind,j_ind) = M_Big(i_ind,j_ind) + Work_g(i,j)
               TauSW_Big(i_ind,j_ind) = TauK_g(i,j)
            enddo
         enddo
         nu_ind = nu_ind + mdim
!
      enddo NuLoop
      mu_ind = mu_ind + mdim
!
   enddo MuLoop
!
   deallocate( DeltaM_mu, TauK_g, Work_g  )
!
!  ===================================================================
!  check that nu_ind and mu_ind are correct
!  ===================================================================
   if ( nu_ind /= ndim .or. mu_ind /= ndim ) then
      call ErrorHandler('Tauk','error: nu_ind,mu_ind',nu_ind,mu_ind)
   endif
!
!  ===================================================================
!  invert [1+delta_m*tau(k)]_{mu,nu}
!  ===================================================================
!  -------------------------------------------------------------------
   call lu_inverse( M_Big,ndim )
!  -------------------------------------------------------------------
!
!  ===================================================================
!  compute tau(k)_{mu,nu} = [1+delta_m*tau(k)]^{-1}[tau(k)]_{mu,nu}
!  ===================================================================
!  -------------------------------------------------------------------
   if ( .not.allocated(TauK_Big) ) then
      allocate( TauK_Big(ndim,ndim) )
   endif
   TauK_Big = czero
   call zgemm( 'n','n',ndim,ndim,ndim,cone,M_Big,ndim,TauSW_Big,       &
               ndim,czero,TauK_Big,ndim )
!  -------------------------------------------------------------------
!  ===================================================================
!  deallocate temporary storage for big matrices
!  ===================================================================
   deallocate( M_Big, TauSW_Big )
!
   end subroutine calScreenTauK
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTauK_SpinMuNu( mu, nu, kvec )      result(TauK_Spin)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: mu, nu
   real (kind=RealKind), intent(in)   :: kvec(3)
!
   integer (kind=IntKind) :: nrs, NrsClu
   integer (kind=IntKind) :: indexLFT, indexLFT_max
   integer (kind=IntKind) :: kmax, mdim
   integer (kind=IntKind) :: i
!
   real (kind=RealKind) :: KdotRij, rdif
   real (kind=RealKind), allocatable :: RsCluSwiss(:,:)
   real (kind=RealKind), allocatable :: R_munuLFT(:,:)
!
   complex (kind=CmplxKind) :: e2iKRij
   complex (kind=CmplxKind), pointer :: TauK(:,:)
   complex (kind=CmplxKind), pointer :: TauK_Spin(:,:)
!
   kmax = kmax_max
   if ( isSpinSpace ) then
      mdim = 2*kmax_max
   else
      mdim = kmax_max
   endif
   allocate( TauK(kmax,kmax), TauK_Spin(mdim,mdim))
!
   NrsClu = ScreenSite(mu)%Neighbors%NumLocalSites
   indexLFT_max = ScreenSite(mu)%NrNu(nu)
   allocate( RsCluSwiss(3,NrsClu), R_munuLFT(3,indexLFT_max) )
   R_munuLFT = ScreenSite(mu)%R_MuNu(:,:,nu)
   RsCluSwiss = ScreenSite(mu)%Neighbors%Position(:,:)
!
   allocate( TauK(kmax,kmax) )
!
   TauK = czero
   rdif = zero
!
!  ===================================================================
!  Calculates the [mu,nu] block of the lattice Fourier Transform
!  of tau_0j_swiss. Puts the result in TauK_swiss
!  ===================================================================
!
   NuClu_Loop: do indexLFT = 1,indexLFT_max
      SwissClu_Loop: do nrs = 1,NrsClu
         do i = 1,3
            rdif = rdif+( R_munuLFT(i,indexLFT) -                      &
                          (RsCluSwiss(i,1) - RsCluSwiss(i,nrs)) )**2
         enddo
         if (rdif < ten2m6) then
!           ==========================================================
!           calculate exp(i*k.R_ij)
!           ==========================================================
!
            KdotRij = kvec(1)*R_munuLFT(1,indexLFT) +                  &
                      kvec(2)*R_munuLFT(2,indexLFT) +                  &
                      kvec(3)*R_munuLFT(3,indexLFT)
            e2iKRij = exp( sqrtm1*KdotRij )
!           ==========================================================
!           calculates exp(i*K*R_0j)*tau0j
!           ==========================================================
            TauK(:,:) = TauK(:,:) +                                    &
                        e2iKRij*ScreenM(UniqPots(mu))%tauS0jp(:,:,nrs)
         endif
         rdif = zero
      enddo SwissClu_Loop
   enddo NuClu_Loop
!
   deallocate( R_munuLFT, RsCluSwiss )
!
!  ===================================================================
!  put TauK into spin space: NB it is diagonal
!  ===================================================================
!
   TauK_Spin => SpinSpace( TauK )
   deallocate( TauK )
!
   end function getTauK_SpinMuNu
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calScreenTau0J(myatom,epsilon,isRelativistic)
!  ===================================================================
   use RSpaceStrConstModule, only : getGij => getStrConstMatrix
   use PhysParamModule, only : LightSpeed
!
   implicit none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
!
   logical, intent(in) :: isRelativistic
!
   integer (kind=IntKind), intent(in) :: myatom
!
   complex (kind=CmplxKind), intent(in) :: epsilon
!
!  ===================================================================
!  Local variables
!  ===================================================================
!  *******************************************************************
!  Tau0jSwiss - row of block matrices with length equals to
!                 the number of neighbors in the cluster
!  Tau00Swiss  - saves the 00-block of Tau0jSwiss
!                 kmax by kmax matrix
!  WsmallSwiss - kmax by kmax matrix used to save t_swiss*gij
!  WBig         - the local cluster matrix 1-Tswiss*G
!                 NrMat by NrMat matrix
!
!  *******************************************************************
   character (len=14), parameter :: sname = 'getScreenTauIJ'
!
   integer (kind=IntKind) :: n, i, j
   integer (kind=IntKind) :: nrst, ncst
   integer (kind=IntKind) :: lmaxi, lmaxj
   integer (kind=IntKind) :: ir1
   integer (kind=IntKind) :: ir2
   integer (kind=IntKind) :: nrmat_swiss
   integer (kind=IntKind) :: nrs_swiss
   integer (kind=IntKind) :: atom_type
   integer (kind=IntKind) :: iprint
!
   real (kind=RealKind) :: rij(3)
   real (kind=RealKind) :: rs_swissi(3), rs_swissj(3)
!
   complex (kind=CmplxKind), pointer :: gij(:,:)
   complex (kind=CmplxKind) :: prel
!
   complex (kind=CmplxKind), pointer :: Tau0jSwiss(:,:,:)
   complex (kind=CmplxKind), pointer :: Tau00Swiss(:,:)
   complex (kind=CmplxKind), allocatable :: WsmallSwiss(:,:)
!
   complex (kind=CmplxKind), allocatable :: WBig(:,:)
!
!  *******************************************************************
!  This subroutine calculates Tau0jSwiss
!  *******************************************************************
!
   nrs_swiss = ScreenSite(myatom)%Neighbors%NumLocalSites
   nrmat_swiss = ScreenSite(myatom)%NrMat
!
   if ( abs(energy-epsilon) > ten2m8 ) then
      energy = epsilon
!     ================================================================
!     calculates the Tswiss of each unique site
!     ================================================================
      do n = 1,NumUniqPots
         call calTswiss(n)
      enddo
   endif
!
   if ( associated( ScreenM(UniqPots(myatom))%tauS0jp ) ) then
      call WarningHandler(sname,'call delScreenMatrices() to clear space')
      call delScreenMatrices(UniqPots(myatom))
   endif
!
!  ===================================================================
!  setup the dimensions of the cluster matrices of myatom
!     nrs_swiss   - number of neighbors in the cluster
!     nrmat_swiss - dimension of the cluster matrix
!  ===================================================================
!
   if ( UniqPots(myatom) == 1 ) then
      allocate( Tau0jSwiss(kmax_max,kmax_max,nrs_swiss),               &
                WsmallSwiss(kmax_max,kmax_max),                        &
                WBig(nrmat_swiss,nrmat_swiss) )
   else
      allocate( WsmallSwiss(kmax_max,kmax_max),                        &
                WBig(nrmat_swiss,nrmat_swiss) )
   endif
   Tau0jSwiss = czero
   WsmallSwiss = czero
!
!  ===================================================================
!  setup unit matrix
!  ===================================================================
!  -------------------------------------------------------------------
   WBig = czero
   do n = 1, nrmat_swiss
      WBig(n,n) = cone
   enddo
!
!  ===================================================================
!  for each Rij set up Tau_swiss=[1-t_swiss*G]** -1 matrix
!  ===================================================================
!
   if ( isRelativistic ) then
      prel=sqrt(energy*(cone+energy/(LightSpeed**2)))
   else
      prel=sqrt(energy)
   endif
!
   nrst=0
!   write(6,*)"getTau0J :: before ij loop "
   do ir1 = 1,nrs_swiss
      ncst = 0
      lmaxi = ScreenSite(myatom)%Neighbors%LmaxClu(ir1)
      rs_swissi = ScreenSite(myatom)%Neighbors%Position(:,ir1)
      atom_type = ScreenSite(myatom)%Neighbors%AtomType(ir1)
      do ir2 = 1,nrs_swiss
         if (ir1 /= ir2) then
            rs_swissj = ScreenSite(myatom)%Neighbors%Position(:,ir2)
            lmaxj = ScreenSite(myatom)%Neighbors%LmaxClu(ir2)
            rij(1:3) = rs_swissi(1:3)-rs_swissj(1:3)
            WsmallSwiss = czero
!            write(6,*)"getTau0J :: rs_swissi rs_swissj"
!            write(6,*) rs_swissi, rs_swissj
!
!           ==========================================================
!           g(Rij) calculation
!           ==========================================================
            if (ir1 == 1 .and. ir2 == 2) then
               iprint = 2
            else
               iprint = 0
            endif
!           ----------------------------------------------------------
            allocate( gij((lmaxi+1)**2,(lmaxj+1)**2) )
            gij => getGij( prel, rij, lmaxi, lmaxj)
!           ----------------------------------------------------------
!            write(6,'(''getTau0J :: gij'',2i3)')ir1,ir2
!            call writeMatrix('Gij ::',gij,(lmaxi+1)**2,(lmaxj+1)**2)
!           ==========================================================
!           form t_swiss * g_ij for current Rij: => WsmallSwiss
!           ==========================================================
!           ----------------------------------------------------------
            call zgemm( 'n', 'n', kmax_max, kmax_max, kmax_max, cone,  &
                        ScreenM(UniqPots(atom_type))%Tswiss(:,:),      &
                        kmax_max, gij, kmax_max, czero,                &
                        WsmallSwiss, kmax_max )
!           ----------------------------------------------------------
!           ==========================================================
!           load the current block into WBig 1-t*G
!           ==========================================================
!           ----------------------------------------------------------
            call cmtrins( kmax_max, nrmat_swiss, kmax_max, kmax_max,   &
                          nrst, ncst, -CONE, WsmallSwiss, kmax_max,    &
                          WBig, nrmat_swiss)
!           ----------------------------------------------------------
            deallocate( gij )
            nullify(gij)
         endif
         ncst = ncst+kmax_max
      enddo
      nrst = nrst+kmax_max
   enddo
!
!  ===================================================================
   deallocate( WsmallSwiss )
   if ( nrst /= ncst .or. nrst /= nrmat_swiss ) then
      call ErrorHandler( 'getScreenTau0J',                             &
                         'nrst ne ncst or nrst ne nrmat_swiss',        &
                         ncst, nrmat_swiss)
      call StopHandler(sname,'nrst differ on ncst or nrmat',nrst)
   endif
!
!  ===================================================================
!  create Tau0jSwiss => {[1-t*G]**(-1)}*t ( block matrices row)
!  ===================================================================
!  -------------------------------------------------------------------
   call getTAU0j_RowBlk( myatom, kmax_max, nrmat_swiss, nrs_swiss,     &
                         WBig, Tau0jSwiss )
!   write(6,*)"getTau0J:: end of calculations"
!  -------------------------------------------------------------------
!  ===================================================================
!  load 00-block of Tau0jSwiss to Tau00Swiss: used in n(e) calc.
!  ===================================================================
!
   if ( .not.associated(ScreenM(UniqPots(myatom))%tauS0jp) ) then
      allocate(ScreenM(UniqPots(myatom))%                              &
                                  tauS0jp(kmax_max,kmax_max,nrs_swiss))
      allocate(ScreenM(UniqPots(myatom))%tauS00p(kmax_max,kmax_max))
   endif
!
   ScreenM(UniqPots(myatom))%tauS0jp = Tau0jSwiss
   ScreenM(UniqPots(myatom))%tauS00p = Tau0jSwiss(:,:,1)
!
   deallocate( WBig, Tau0jSwiss )
!
   end subroutine calScreenTau0J
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calTswiss(n)
!  ===================================================================
!
!  *******************************************************************
!  set up matrix t_swiss as the non-relativistic t-matrix of a
!  square well
!  *******************************************************************
   use BesselModule, only : SphericalBessel, SphericalNeumann
!
   implicit none
!
   character (len=9), parameter :: sname = 'calTswiss'
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: lm
   integer (kind=IntKind) :: l1max
!
   complex (kind=CmplxKind) :: eprime
   complex (kind=CmplxKind) :: rte
   complex (kind=CmplxKind) :: rteprime
   complex (kind=CmplxKind) :: tand
   complex (kind=CmplxKind) :: react
   complex (kind=CmplxKind) :: fmb(0:lmax_max+1)
   complex (kind=CmplxKind) :: fmn(0:lmax_max+1)
   complex (kind=CmplxKind) :: fmh(0:lmax_max+1)
   complex (kind=CmplxKind) :: fmodb(0:lmax_max+1)
   complex (kind=CmplxKind) :: fmodn(0:lmax_max+1)
   complex (kind=CmplxKind) :: fmodh(0:lmax_max+1)
!
   rte = sqrt(energy)
   eprime = energy-UniqV(n)
   rteprime = sqrt(eprime)
   if(print_level > 0) then
      write(6,*) 'energy    in calTswiss:',energy
      write(6,*) 'rmt_swiss in calTswiss:',UniqRmt(n)
      write(6,*) 'pot_swiss in calTswiss:',UniqV(n)
   endif
!
!  -------------------------------------------------------------------
   call SphericalBessel(lmax_max+1,rte*UniqRmt(n),fmb)
   call SphericalNeumann(lmax_max+1,rte*UniqRmt(n),fmn)
!  -------------------------------------------------------------------
   call SphericalBessel(lmax_max+1,rteprime*UniqRmt(n),fmodb)
   call SphericalNeumann(lmax_max+1,rteprime*UniqRmt(n),fmodn)
!  -------------------------------------------------------------------
   allocate( ScreenM(n)%Tswiss(kmax_max,kmax_max) )
   ScreenM(n)%Tswiss(1:kmax_max,1:kmax_max) = czero
!  -------------------------------------------------------------------
   lm=0
   do l=0,lmax_max
      tand = (rte*fmb(l+1)*fmodb(l)-rteprime*fmb(l)*fmodb(l+1))/       &
             (rte*fmn(l+1)*fmodb(l)-rteprime*fmn(l)*fmodb(l+1))
      react = -tand/rte
      do m = -l,l
         lm = lm+1
         ScreenM(n)%Tswiss(lm,lm) = react/(one+sqrtm1*rte*react)
      end do
   end do
!  ===================================================================
   if ( print_level > 0 ) then
      write(6,*)'  calTswiss ::'
!     ----------------------------------------------------------------
      call writeMatrix( 'Tswiss:: block 1', ScreenM(n)%Tswiss(:,:),    &
                        kmax_max, kmax_max)
!     ----------------------------------------------------------------
   endif
!
!
   if ( sname.eq.stop_routine ) then
      call StopHandler(sname)
   endif
   end subroutine calTswiss
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function SpinSpace(InputMatrix) result(OutMatrix)
!  ===================================================================
!
   implicit   none
!
!  ===================================================================
!   I/O variables...................................................
!  ===================================================================
!
   complex (kind=CmplxKind), intent(in) :: InputMatrix(:,:)
!
!  ===================================================================
!  Local variables.................................................
!  ===================================================================
   integer (kind=IntKind) :: i, j, n
   complex (kind=CmplxKind), pointer :: OutMatrix(:,:)
!
!  *******************************************************************
!  copies a matrix into spin space : ..............................
!  ******************************************************************
!
!  ===================================================================
!  printout if needed: ............................................
!  ===================================================================
   if ( isSpinSpace ) then
      allocate( OutMatrix(2*kmax_max,2*kmax_max) )
      OutMatrix = czero
      OutMatrix(1:kmax_max,1:kmax_max)=InputMatrix(:,:)
      do i = 1,kmax_max
         do j = 1,kmax_max
            OutMatrix(kmax_max+j,kmax_max+i) = InputMatrix(j,i)
         enddo
      enddo
   else
      allocate( OutMatrix(kmax_max,kmax_max) )
      OutMatrix(:,:) = InputMatrix(:,:)
   endif
!
   return
   end function SpinSpace
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDeltaM( myAtom )
!  ===================================================================
   use SingleScatteringModule, only : getTMatrix
!
   implicit   none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
!
   integer (kind=IntKind), intent(in) :: myAtom
!
!  ====================================================================
!  Local variables
!  ====================================================================
!
   integer (kind=IntKind) :: nsize
   integer (kind=IntKind) :: i,n
   integer (kind=IntKind),allocatable :: lu_piv(:)
!
   complex (kind=CmplxKind), pointer :: DeltaM(:,:)
   complex (kind=CmplxKind), allocatable :: Mswiss(:,:)
   complex (kind=CmplxKind), allocatable :: Msystm(:,:)
   complex (kind=CmplxKind), allocatable :: Tsystm(:,:)
   complex (kind=CmplxKind), allocatable :: lu_wrk(:)
!
!  ===================================================================
!  calculates : DeltaM = [Msystm -  Mswiss ]
!  ===================================================================
!
   nsize = kmax_max*nSpinCant
!  ===================================================================
!  allocate work arrays:
!  ===================================================================
   allocate( DeltaM(nsize,nsize), Tsystm(nsize,nsize) )
   allocate( Mswiss(nsize,nsize), Msystm(nsize,nsize) )
!
   Tsystm = czero
   do i = 1,nSpinCant
      Tsystm((i-1)*kmax_max+1:i*kmax_max,(i-1)*kmax_max+1:i*kmax_max)= &
                           getTMatrix(UniqPots(myAtom),energy,i,czero)
   enddo
!
   Msystm = Tsystm
!
   do n = 1,NumUniqPots
      DeltaM = czero
      allocate( ScreenM(n)%TswissSpin(nsize,nsize) )
      ScreenM(n)%TswissSpin => SpinSpace(ScreenM(n)%Tswiss(:,:))
!
!     ================================================================
!     printout if needed: TswissSpin , Tsystm
!     ================================================================
!
      if (print_level > 1) then
         write(6,'('' GET_DeltaM:: n_spin_cant ='',i5)') nSpinCant
         write(6,'('' GET_DeltaM::       kmax_max ='',i5)') kmax_max
         write(6,'('' GET_DeltaM::       nsize ='',i5)') nsize
         write(6,'('' GET_DeltaM:: [TswissSpin]'')')
!        -------------------------------------------------------------
         call writeMatrix( 'TswissSpin',ScreenM(n)%TswissSpin,         &
                           nsize,nsize )
!        -------------------------------------------------------------
         write(6,'('' GET_DeltaM:: [Tsystm]'')')
!        -------------------------------------------------------------
         call writeMatrix('calDeltaM : Tsystm', Tsystm,nsize,nsize )
!        -------------------------------------------------------------
      endif
!
!     ================================================================
!     calculate Mswiss [TswissSpin^(-1)] and Msystm [Tsystm^(-1)]
!     ================================================================
!
      Mswiss = czero
      Mswiss = ScreenM(n)%TswissSpin
!
!      call writeMatrix('TswissSpin',ScreenM(n)%TswissSpin,nsize,nsize)
!     ----------------------------------------------------------------
      call lu_inverse( Mswiss, nsize )
!     ----------------------------------------------------------------
      call lu_inverse( Msystm, nsize )
!     ----------------------------------------------------------------
!     ================================================================
!     calculate DeltaM(i,j)
!     ================================================================
!
      DeltaM =  Msystm - Mswiss
!
!     ================================================================
!     printout if needed: Mswiss , m_matrix and delta-m
!     ================================================================
      if (print_level > 1) then
         write(6,'(//,'' GET_DeltaM:: [Mswiss]'')')
!        ----------------------------------------------------------
         call writeMatrix( 'Mswiss',Mswiss, nsize, nsize)
!        ----------------------------------------------------------
         write(6,'('' GET_DeltaM:: [Msystm]'')')
!        ----------------------------------------------------------
         call writeMatrix( 'Msystm',Msystm, nsize, nsize )
!        ----------------------------------------------------------
         write(6,'('' GET_DeltaM:: [DeltaM]'')')
!        ----------------------------------------------------------
         call writeMatrix( 'DeltaM', DeltaM, nsize, nsize)
!        ----------------------------------------------------------
      endif
!
      allocate ( ScreenM(n)%deltaMp(nsize,nsize) )
      ScreenM(n)%deltaMp = DeltaM
   enddo
!  ===================================================================
!  deallocate work arrays:
!  ===================================================================
!
   deallocate( Mswiss, Msystm, Tsystm )
!
!  ===================================================================
   return
!
   end subroutine calDeltaM
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getClustDim()
!  ===================================================================
!  needs NrClu_Max and NrSiteClu_Max to be defined as module variables
!  and updates their values
!
   implicit none
!
   integer (kind=IntKind) :: m, n, n1
   integer (kind=IntKind) :: i0, j0, k0, n0
   integer (kind=IntKind) :: nn(3)
   integer (kind=IntKind) :: nn_size
   integer (kind=IntKind) :: n_max, nrsclu, nrsclu_mn
   integer (kind=IntKind), parameter :: n_max0 = 5
   integer (kind=IntKind), allocatable :: i(:), j(:), k(:)
!
   real (kind=RealKind) :: r(3), posv(3), posv_mn(3), shift(3)
   real (kind=RealKind) :: atdist, atdist_mn, Rcirclu
!
   NrClu_Max = 0
   NrSiteClu_Max = 0
   NrSiteMuNu_Max = 0
!
   r(1:3) = bravais(1,1:3)**2 + bravais(2,1:3)**2 + bravais(3,1:3)**2
   r(1:3) = sqrt(r(1:3))
!
   MainLoop: do n = 1,NumAtoms
!
      Rcirclu = ScreenSite(n)%Rcut
      nn(1:3) = int(Rcirclu/r(1:3)+0.9)
      nn_size = (2*nn(1)+1)*(2*nn(2)+1)*(2*nn(3)+1)
      n_max=int( half*((nn_size)**third-1)+0.9)
      if ( n_max > n_max0 ) then
         call ErrorHandler('getClustDim','n_max exceeds the upper limit', &
                           n_max,(2*n_max0+1)**3)
      endif
      allocate( i(nn_size), j(nn_size), k(nn_size))
      n0 = 0
      do i0 = -nn(1),nn(1)
         do j0 = -nn(2),nn(2)
            do k0 = -nn(3),nn(3)
               n0 = n0+1
               i(n0) = i0
               j(n0) = j0
               k(n0) = k0
            enddo
         enddo
      enddo

!     ================================================================
!     calculate neighbors that are closer than Rcirclu
!
      if ( abs(Rcirclu) < ten2m8 ) then
!        =============================================================
!        special case Rcirclu=0.0 then default to single site
!
         nrsclu = 1
      else
         nrsclu = 0
         nrsclu_mn = 0
         do m = 1,n0
            shift(1:3) = i(m)*bravais(1:3,1) + j(m)*bravais(1:3,2)     &
                   + k(m)*bravais(1:3,3)- ScreenSite(n)%Position(1:3)
            Num_Atoms: do n1 = 1,NumAtoms
               atdist = zero
               posv(1:3) = shift(1:3) + ScreenSite(n1)%Position(1:3)   &
                                      - ScreenSite(n)%Position(1:3)
               posv_mn(:) =shift(:) + ScreenSite(n)%Position(1:3)    &
                                       -ScreenSite(n1)%Position(1:3)
               atdist = posv(1)**2 + posv(2)**2 + posv(3)**2
               atdist = sqrt(atdist)
               atdist_mn = posv_mn(1)**2 + posv_mn(2)**2 + posv_mn(3)**2
               atdist_mn = sqrt(atdist_mn)
               if ( atdist <= Rcirclu ) then
                  nrsclu = nrsclu + 1
               endif
               if (atdist_mn.lt.Rcirclu) then
                  nrsclu_mn = nrsclu_mn + 1
               endif
            enddo Num_Atoms
         enddo
      endif
!
      deallocate( i,j,k )
!     ================================================================
!     updates the cluster dimensions
!
      NrClu_Max = max(NrClu_Max,n_max)
      NrSiteClu_Max = max(NrSiteClu_Max,nrsclu)
      NrSiteMuNu_Max = max(NrSiteMuNu_Max,nrsclu_mn)
!
   enddo MainLoop
!
   end subroutine getClustDim
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine neighbors( Rcirclu, myatom )
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), parameter :: n_max = 5
   integer (kind=IntKind), intent(in) ::  myatom
   integer (kind=IntKind) ::  nrsclu, nrsclu_mn(NumAtoms)
   integer (kind=IntKind) ::  NrMat, max_NrNu
   integer (kind=IntKind) ::  n, m
   integer (kind=IntKind) ::  i0, j0, k0, n0
   integer (kind=IntKind) ::  nn(3)
   integer (kind=IntKind) ::  nn_size
   integer (kind=IntKind) ::  i(3,(2*n_max+1)**3)
   integer (kind=IntKind) ::  index(3)
   integer (kind=IntKind) ::  AtomType(NrSiteClu_MAX)
   integer (kind=IntKind) ::  AtomMode(NrSiteClu_MAX)
   integer (kind=IntKind) ::  AtomIndexIJK(3,NrSiteClu_Max)
   integer (kind=IntKind), allocatable ::  LmaxClu(:)
!
   real (kind=RealKind), intent(in) ::  Rcirclu
   real (kind=RealKind) ::  R_munu(3,NrSiteMuNu_Max,NumAtoms)
   real (kind=RealKind) ::  atdist, atdist_mn
   real (kind=RealKind) ::  atdsqr, atdsqr_mn
   real (kind=RealKind) ::  r(3)
   real (kind=RealKind) ::  r0
   real (kind=RealKind) ::  shift(3)
   real (kind=RealKind) ::  posv(3), posv_mn(3)
   real (kind=RealKind) ::  rsclu(3,NrSiteClu_Max)
   real (kind=RealKind) ::  rsclusq(NrSiteClu_Max)
!
!  *******************************************************************
!  Calculate neighbors to be used in the cluster calculation
!  *******************************************************************
!
!  ===================================================================
!  determine i,j,k data set
!
   r(1:3) = sqrt( bravais(1,1:3)**2 + bravais(2,1:3)**2                &
                + bravais(3,1:3)**2 )
   nn(1:3) = max(n_max,int(Rcirclu/r(1:3)+0.9))
   nn_size = (2*nn(1)+1)*(2*nn(2)+1)*(2*nn(3)+1)
   if ( nn_size .gt. (2*n_max+1)**3 ) then
      call ErrorHandler( 'neighbors','n0 exceeds the upper limit.',    &
                         nn_size,(2*n_max+1)**3)
   endif
   n0 = 0
   do i0 = -nn(1),nn(1)
      do j0 = -nn(2),nn(2)
         do k0 = -nn(3),nn(3)
            n0 = n0+1
            i(1,n0) = i0
            i(2,n0) = j0
            i(3,n0) = k0
         enddo
      enddo
   enddo
!  ===================================================================
!  zero out strorage arrays
!  -------------------------------------------------------------------
   rsclu = zero
   rsclusq = zero
   R_munu = zero
!  -------------------------------------------------------------------
   AtomMode = 1
!
!  ===================================================================
!  calculate neighbors that are closer than Rcirclu
!
   if ( abs(Rcirclu) < ten2m8 ) then
!     ================================================================
!     special case Rcirclu=0.0 then default to single site
!
      AtomMode(1) = myatom
      nrsclu = 1
   else
      nrsclu = 0
      nrsclu_mn = 0
      do m = 1,n0
         shift(1:3) = i(1,m)*bravais(1:3,1)+i(2,m)*bravais(1:3,2)      &
                    + i(3,m)*bravais(1:3,3)
!         index(1) = i(m)
!         index(2) = j(m)
!         index(3) = k(m)
         do n = 1,NumAtoms
            posv(1:3) = shift(1:3) + ScreenSite(n)%Position(1:3)       &
                        -ScreenSite(myatom)%Position(1:3)
            atdsqr = posv(1)**2 + posv(2)**2 + posv(3)**2
            atdist = sqrt(atdsqr)
            posv_mn(1:3) = shift(1:3) - ScreenSite(n)%Position(1:3)       &
                           + ScreenSite(myatom)%Position(1:3)
            atdsqr_mn =posv_mn(1)**2 + posv_mn(2)**2 + posv_mn(3)**2
            atdist_mn = sqrt(atdsqr_mn)
            if ( atdist <= Rcirclu ) then
               if ( nrsclu+1 > NrSiteClu_Max ) then
                  call ErrorHandler('neighbors','',myatom,nrsclu)
               endif
!              =======================================================
!              insert posv(3) in rsclu(3*nrsclu) list such that they
!              are in order of increasing length
!              -------------------------------------------------------
               call sort_rsclu( posv, atdsqr, n, AtomType, AtomMode,   &
                        AtomIndexIJK, i(:,m), rsclu, rsclusq, nrsclu )
!              -------------------------------------------------------
            endif
            if ( atdist_mn < Rcirclu ) then
               if ( nrsclu_mn(n)+1 > NrSiteMuNu_Max ) then
                  call ErrorHandler('neighbors','',myatom,nrsclu_mn(n))
               endif
               nrsclu_mn(n) = nrsclu_mn(n) + 1
               R_munu(1:3,nrsclu_mn(n),n) = posv_mn(1:3)
            endif
         enddo
      enddo
   endif
   allocate ( LmaxClu(nrsclu) )
   call size_clu_mat( Rcirclu, rsclusq, nrsclu,                        &
                      NrMat, LmaxClu, ScreenSite(myatom)%Lmax)
!
!  ===================================================================
!  update the neighbors of myatom
!
   ScreenSite(myatom)%NrMat = NrMat
   ScreenSite(myatom)%NrNu = nrsclu_mn
   max_NrNu = 0
   do i0 = 1,NumAtoms
      max_NrNu = max(max_NrNu,ScreenSite(myatom)%NrNu(i0))
   enddo
   ScreenSite(myatom)%Neighbors%NumLocalSites=nrsclu
   allocate ( ScreenSite(myatom)%Neighbors%LmaxClu(nrsclu),             &
              ScreenSite(myatom)%Neighbors%AtomType(nrsclu),            &
              ScreenSite(myatom)%Neighbors%Position(3,nrsclu),          &
              ScreenSite(myatom)%Neighbors%IndexIJK(3,nrsclu),          &
              ScreenSite(myatom)%R_MuNu(3,max_NrNu,NumAtoms) )
   ScreenSite(myatom)%R_MuNu(:,:,:) = R_munu(:,:,:)
   do m = 1,nrsclu
      n = AtomType(m)
      ScreenSite(myatom)%Neighbors%LmaxClu(m) = LmaxClu(m)
      ScreenSite(myatom)%Neighbors%AtomType(m) = n
      ScreenSite(myatom)%Neighbors%Position(1:3,m) = rsclu(1:3,m)
      ScreenSite(myatom)%Neighbors%IndexIJK(:,m) = AtomIndexIJK(:,m)
   enddo
!
   deallocate(LmaxClu)
!
   if ( AtomMode(1) /= myatom ) then
      call ErrorHandler( 'neighbors','AtomMode(1), myatom ',           &
                         AtomMode(1),myatom )
   endif
!
   end subroutine neighbors
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sort_rsclu( posv, atdsqr, n, AtomType, AtomMode,         &
                      AtomIndexIJK, IndexIJK, rsclu, rsclusq, nrsclu )
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(inout) ::  nrsclu
   integer (kind=IntKind), intent(inout) ::  AtomType(nrsclu+1)
   integer (kind=IntKind), intent(inout) ::  AtomMode(nrsclu+1)
   integer (kind=IntKind), intent(inout) ::  AtomIndexIJK(3,nrsclu+1)
   integer (kind=IntKind), intent(in) :: IndexIJK(3)
   integer (kind=IntKind), intent(in) ::  n
!
   integer (kind=IntKind) ::  ic
   integer (kind=IntKind) ::  nv
   integer (kind=IntKind) ::  nvm
!
   real (kind=RealKind), intent(inout) ::  rsclu(3,nrsclu+1)
   real (kind=RealKind), intent(inout) ::  rsclusq(nrsclu+1)
   real (kind=RealKind), intent(in) ::  posv(3)
   real (kind=RealKind), intent(in) ::  atdsqr
!
!  *******************************************************************
!  inserts a vector in a list of vectors such that they are in
!  order of increasing length
!  *******************************************************************
!
!  ===================================================================
   nv = 1
   do while( nv <= nrsclu .and. rsclusq(nv) < atdsqr )
      nv = nv+1
   enddo
   do nvm = nrsclu,nv,-1
      rsclusq(nvm+1) = rsclusq(nvm)
      AtomMode(nvm+1) = AtomMode(nvm)
      AtomType(nvm+1) = AtomType(nvm)
      AtomIndexIJK(:,nvm+1) = AtomIndexIJK(:,nvm)
      rsclu(1:3,nvm+1) = rsclu(1:3,nvm)
   enddo
   rsclusq(nv) = atdsqr
   AtomMode(nv) = n
   AtomType(nv) = ScreenSite(n)%AtomType
   AtomIndexIJK(:,nv) = IndexIJK(:)
   rsclu(1:3,nv) = posv(1:3)
   nrsclu = nrsclu+1
!  ===================================================================
!
   return
   end subroutine sort_rsclu
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine size_clu_mat( rcirclu, rsclusq, nrsclu,                  &
                            nrmat, lmaxclu, lmax)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nrsclu
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind), intent(inout) :: nrmat
   integer (kind=IntKind), intent(inout) :: lmaxclu(nrsclu)
!
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: lkeep
!
   real (kind=RealKind), intent(in) :: rcirclu
   real (kind=RealKind), intent(in) :: rsclusq(nrsclu)
   real (kind=RealKind) :: rdis
   real (kind=RealKind) :: rsteps(lmax+1)
!
!  ===================================================================
!  Set array that controls stepping down of lmax LIZ shells
!
   rsteps(1) = 11.500
   rsteps(2) = 13.200
   rsteps(3) = 14.000
   rsteps(4) = 99.900
   do i = 5,lmax+1
      rsteps(i) = rcirclu
   enddo
   if ( lmax > 4 ) then
      do i = lmax,lmax-3,-1
         rsteps(i) = rsteps(i-lmax+4)
      enddo
      do i = 1,lmax-4
         rsteps(i) = ten2m6
      enddo
   endif
!
!  ===================================================================
!  set up the size of the cluster matrix
!
   nrmat = 0
   do i = 1,nrsclu
      rdis = sqrt(rsclusq(i))
      lkeep = lmax
      do n1 = 1,lmax+1
         if ( rdis > rsteps(n1) ) then
            lkeep = lkeep-1
         endif
      enddo
      lmaxclu(i) = lkeep
      nrmat = nrmat + (lkeep+1)*(lkeep+1)
   enddo
!
   end subroutine size_clu_mat
!
!  ===================================================================
!
!  *******************************************************************
!  matrix operation realted subroutines
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getTAU0j_RowBlk( myatom, kkrsz, nrmat_swiss, nrs_swiss,  &
                               WBig, Tau_j0 )
!  ===================================================================
!
   implicit none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
!
   integer (kind=IntKind), intent(in) :: kkrsz
   integer (kind=IntKind), intent(in) :: nrmat_swiss
   integer (kind=IntKind), intent(in) :: nrs_swiss
   integer (kind=IntKind), intent(in) :: myatom
!
   complex (kind=CmplxKind), intent(inout) :: WBig(nrmat_swiss,nrmat_swiss)
   complex (kind=CmplxKind), intent(inout) :: Tau_j0(kkrsz,kkrsz,nrs_swiss)
!
!  ===================================================================
!  Local variables
!  ===================================================================
!
   integer (kind=IntKind) :: ir1
   integer (kind=IntKind) :: atom_type
!
!  ===================================================================
!  invert [1-t*G] : obtain the kmax_max*kmax_max block only
!  ===================================================================
!
!  -------------------------------------------------------------------
   call block_inv_t0j( WBig, kkrsz, nrmat_swiss )
!  -------------------------------------------------------------------
!  ===================================================================
!  Calculate Tau_j0: For this case all the single-site t-matrices
!  in the local interaction zone(LIZ) is needed to calculate Tau_j0
!  ===================================================================
   do ir1 = 1,nrs_swiss
!     ----------------------------------------------------------------
      atom_type = ScreenSite(myatom)%Neighbors%AtomType(ir1)
      call zgemm( 'n','n',kkrsz,kkrsz,kkrsz,cone,                      &
                  WBig((ir1-1)*kkrsz+1,1), nrmat_swiss,                &
                  ScreenM(UniqPots(atom_type))%Tswiss(:,:), kkrsz,     &
                  czero, Tau_j0(:,:,ir1), kkrsz)
!     ----------------------------------------------------------------
!      call writeMatrix('Tau0j::', Tau_j0(:,:,ir1),kkrsz,kkrsz)
   enddo
!
!  ===================================================================
!
   return
!
   end subroutine getTAU0j_RowBlk
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine block_inv_t0j( c_matrix, dim_block, dim_matrix )
!  ===================================================================
!
   implicit   none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
!
   integer (kind=IntKind) :: dim_block
   integer (kind=IntKind) :: dim_matrix
!
   complex (kind=CmplxKind) :: c_matrix(dim_matrix,dim_matrix)
!
!  ===================================================================
!  Local variables
!  ===================================================================
!
   integer (kind=IntKind) :: ipiv(dim_matrix)
   integer (kind=IntKind) :: info
!
   complex (kind=CmplxKind) :: vecs(dim_block*dim_matrix)
!
!  *******************************************************************
!  PURPOSE:   inverse the first block of a complex matrix: c_matrix
!             using LU method
!
!  INPUT:     c_matrix,   the complex matrix to be inverted
!             dim_block,  the dimension of the block to be inverted
!             dim_matrix, the dimension of the c_matrix
!
!  OUTPUT:    c_matrix,   contains the first block of the inverted
!                         matrix
!  *******************************************************************
!
!  ===================================================================
!   Invert the KKR-Matrix using the LU method
!  ===================================================================
!  -------------------------------------------------------------------
   call zgetrf( dim_matrix, dim_matrix, c_matrix, dim_matrix, ipiv,info )
!  -------------------------------------------------------------------
   if (info /= 0) then
      call ErrorHandler( 'BLOCK_INV_T0J:: zgetrf', 'bad LU, info', info)
   endif
!  -------------------------------------------------------------------
   call zgetri( dim_matrix,c_matrix,dim_matrix,ipiv,vecs,dim_matrix,info )
!  -------------------------------------------------------------------
   if (info /= 0) then
      call ErrorHandler( 'BLOCK_INV_T0J:: zgetri', 'bad LU, info', info )
   endif
!  -------------------------------------------------------------------
!
   return
!
   end subroutine block_inv_t0j
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cmtrins(ma,mb,kr,kc,rst,cst,fac,a,lda,b,ldb)
!  ===================================================================
!
   implicit   none
!
   integer (kind=IntKind) :: ma
   integer (kind=IntKind) :: mb
   integer (kind=IntKind) :: kr
   integer (kind=IntKind) :: kc
   integer (kind=IntKind) :: rst
   integer (kind=IntKind) :: cst
   integer (kind=IntKind) :: lda
   integer (kind=IntKind) :: ldb
   integer (kind=IntKind) :: i,j
!
   complex (kind=CmplxKind) :: fac
   complex (kind=CmplxKind) :: a(lda,ma)
   complex (kind=CmplxKind) :: b(ldb,mb)
!
!  *******************************************************************
!  a :     a matrix with dimension (lda x ma)
!  =
!
!  b :     a matrix assigned within dimension (kr x kc)
!  =
!  construct a new matrix b with dimension (ldb x mb):
!                         =
!
!          | b       |
!          | =       |
!  b   =   |         |
!  =       |   fac*a |
!          |       = |
!  *******************************************************************
!
!  testing
!     write(*,*)'enter cmtrins'
!
   do j = 1,kc
!     ----------------------------------------------------------------
      call zcopy( kr,a(1,j),1,b(rst+1,cst+j),1 )
      call zscal( kr,fac,b(rst+1,cst+j),1 )
!     ----------------------------------------------------------------
   enddo
!
   return
!
   end subroutine cmtrins
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine lu_inverse( matrix, mat_size )
!  ===================================================================
!
   implicit none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
   integer (kind=IntKind), intent(in) :: mat_size
!
   complex (kind=CmplxKind), intent(inout) :: matrix(mat_size,mat_size)
!
!  ===================================================================
!  Local variables
!  ===================================================================
!
   integer (kind=IntKind)   :: info
   integer (kind=IntKind)   :: ipiv(mat_size)
!
   complex (kind=CmplxKind) :: work(mat_size)
!
!  ===================================================================
!   Find LU factorization of matrix....................................
!  ===================================================================
!  -------------------------------------------------------------------
   call zgetrf( mat_size, mat_size, matrix, mat_size, ipiv, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the LU factorization worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler( 'getDeltaM :: lu_inverse',                    &
                         'problem with zgetrf: info ',                 &
                         info )
      call StopHandler( 'lu_inverse', 'info not equal 0', info)
   endif
!
!  ===================================================================
!  Find inverse of matrix
!  ===================================================================
!  -------------------------------------------------------------------
   call zgetri(mat_size,matrix,mat_size,ipiv,work,mat_size,info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the Inverse worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler( 'getDeltaM :: lu_inverse',                    &
                         'problem with zgetri: info ',                 &
                         info )
      call StopHandler( 'lu_inverse', 'info not equal 0', info)
   endif
!
   end subroutine lu_inverse
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTauBlockJ( mu, nrs )                     result(pTauJ)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: mu,nrs
!
   complex (kind=CmplxKind), pointer :: pTauJ(:,:)
!
   pTauJ => ScreenM(UniqPots(mu))%tauS0jp(:,:,nrs)
!
   end function getTauBlockJ
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTau00K()                     result(pTau00K)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind) :: mdim
!
   complex (kind=CmplxKind), pointer :: pTau00k(:,:)
!
   if ( isSpinSpace ) then
      mdim = 2*kmax_max
   else
      mdim = kmax_max
   endif
!
   if ( .not.allocated( TauK_Big) ) then
      call ErrorHandler( 'getTau00K',                                  &
                         'calScreenTauK should be called first')
   endif
!
   pTau00K => TauK_Big(1:mdim,1:mdim)
!
   end function getTau00K
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTauK_munu(mu,nu)                     result(pTauK_mn)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: mu, nu
!
   integer (kind=IntKind) :: mdim
!
   complex (kind=CmplxKind), pointer :: pTauK_mn(:,:)
!
   if ( mu > NumAtoms .or. nu > NumAtoms ) then
      call ErrorHandler('getTauK_munu','Index too large', mu, nu)
   endif
   if ( .not.allocated( TauK_Big) ) then
      call ErrorHandler( 'getTau00K',                                 &
                         'calScreenTauK should be called first')
   endif
!
   if ( isSpinSpace ) then
      mdim = 2*kmax_max
   else
      mdim = kmax_max
   endif
!
   pTauK_mn => TauK_Big( (mu-1)*mdim+1:(mu-1)*mdim,                   &
                            (nu-1)*mdim+1:(nu-1)*mdim )
!
   end function getTauK_munu
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTauK()                     result(pTauK)
!  ===================================================================
!
   implicit none
!
   complex (kind=CmplxKind), pointer :: pTauK(:,:)
!
   if ( .not.allocated( TauK_Big) ) then
      call ErrorHandler( 'getTau00K',                                 &
                         'calScreenTauK should be called first')
   endif
!
   pTauK => TauK_Big
!
   end function getTauK
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDeltaM( myAtom )                     result(pdelta)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: myAtom
!
   complex (kind=CmplxKind), pointer :: pdelta(:,:)
!
   if ( .not.associated(ScreenM(UniqPots(myAtom))%deltaMp) .and.       &
        myAtom <= NumAtoms ) then
      call calDeltaM( myAtom )
   endif
!
   pdelta => ScreenM(UniqPots(myAtom))%deltaMp(:,:)
!
   end function getDeltaM
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getScreenLmax()                     result(lmax)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind) :: lmax
!
   lmax = lmax_max
!
   end function getScreenLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printScreenTau0J(ns,nb)
!  ===================================================================
!
   implicit   none
!
   integer (kind=IntKind), intent(in) :: ns,nb
!
   integer (kind=IntKind) :: i, j
   integer (kind=IntKind) :: n, nrs
!
   write(6,'(/,a)')'********************************'
   write(6,'( a )')'*  Output from printScreenStr  *'
   write(6,'(a,/)')'********************************'
   write(6,'(80(''=''))')
!
   write(6,*)'energy   ::',energy
   write(6,'(80(''=''))')
!
   do i = 1,NumAtoms
      write(6,'(/,/,a,i5,/)')'     MyAtom   ::', i
      write(6,'(a,3f9.5)')' Position ::',ScreenSite(i)%Position(1:3)
      write(6,*)'Lmax     :: ',ScreenSite(i)%Lmax
      write(6,*)'NrMat    :: ',ScreenSite(i)%NrMat
      write(6,*)'Rmt      :: ',ScreenSite(i)%Rmt
      write(6,*)'Rcut     :: ',ScreenSite(i)%Rcut
      write(6,*)'NrsClu   :: ',ScreenSite(i)%Neighbors%NumLocalSites
      write(6,'(9x,''i'',16x,''rs_swiss'')')
      write(6,'(9x,50(''=''))')
      nrs = ScreenSite(i)%Neighbors%NumLocalSites
      do n = 1,nrs
         write(6,'(5x,i5,3f12.5)')n,                                   &
              (ScreenSite(i)%Neighbors%Position(j,n),j=1,3)
!           (n,(ScreenSite(i)%Neighbors%Position(j,n),j=1,3),n=1,nrs)
      enddo
   enddo
!
   write(6,'(/,a,/)')'ScreenTau0J block matrices(non-zero elements) ::'
!
   do i = 1,NumUniqPots
      write(6,'(/,a,i2,/)')"    Site ::",i
      nrs = ScreenSite(i)%Neighbors%NumLocalSites
      if ( associated( ScreenM(i)%tauS0jp ) ) then
         do n = 1,nrs
            if (n <= nb .and. n <= nrs ) then
               write(6,'(/,a,i5.3,/)')'TauOJ :: block', n
               call writeMatrix( 'Tau0j::',ScreenM(i)%tauS0jp(:,:,n),   &
                                 kmax_max,kmax_max )
            endif
         enddo
      else
         write(6,*)'tauS0jp not allocated at this point'
      endif
!
      if ( associated( ScreenM(i)%tauS00p ) ) then
         write(6,'(/,a,/)')'TauSO0 :: '
         call writeMatrix( 'Tau00', ScreenM(i)%tauS00p(:,:), kmax_max, &
                           kmax_max)
      else
         write(6,*)'TauS00 not allocated at this point'
      endif
!
      if ( associated( ScreenM(i)%Tswiss ) ) then
         write(6,'(/,a,i5.3/)')'Tswiss :: block', i
         call writeMatrix( 'Tswiss::', ScreenM(i)%Tswiss(:,:),         &
                           kmax_max, kmax_max )
      else
         write(6,*)'Tswiss not allocated at this point'
      endif
      if ( associated( ScreenM(i)%deltaMp ) ) then
         write(6,'(/,a,i5.3/)')'DeltaM :: block', i
         call writeMatrix( 'DeltaM::', ScreenM(i)%deltaMp(:,:),         &
                           kmax_max*ns, kmax_max*ns )
      else
         write(6,*)'DeltaM not allocated at this point'
      endif
   enddo
!
   write(6,'(80(''-''))')
!
   end subroutine printScreenTau0J
!  ===================================================================
!
end module ScreenStrConstModule


