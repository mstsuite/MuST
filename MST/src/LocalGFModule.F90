module LocalGFModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, THREE, PI, PI2, PI4, &
                               CZERO, CONE, TWO, HALF, SQRTm1, Y0, TEN2m6
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
   use IntegerFactorsModule, only : lofk, lofj, mofj, m1m, mofk, jofk
   use TimerModule, only : getTime, getRoutineCallTiming
   use MPPModule, only : MyPE, syncAllPEs
!
public :: initLocalGF,          &
          endLocalGF,           &
          computeLocalGF,       &
          getLocalGF,           &
          getNumOrbitals,       &
          getEnergyMesh,        &
          printLocalGF
!
   interface getLocalGF
      module procedure getLocalGF_v0, getLocalGF_v2, getLocalGF_v4, &
                       getLocalGF_v5
   end interface getLocalGF
!
private
   integer (kind=IntKind) :: lmax_kkr
   integer (kind=IntKind) :: kmax_kkr
   integer (kind=IntKind) :: lmax_gf
   integer (kind=IntKind) :: kmax_gf
   integer (kind=IntKind) :: NumLocalAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: Relativity
   integer (kind=IntKind) :: rspace_size
   integer (kind=IntKind) :: fspace_size
   integer (kind=IntKind) :: pspace_size
   integer (kind=IntKind) :: pint_size
   integer (kind=IntKind) :: MaxNumOrbitals
   integer (kind=IntKind), parameter :: BasisType = 0 ! The default local orbitals
                                                      ! are the local atomic orbitals 
!
   complex (kind=CmplxKind), target, allocatable :: rspace(:)
   complex (kind=CmplxKind), target, allocatable :: fspace(:)
   complex (kind=CmplxKind), target, allocatable :: pspace(:)
   complex (kind=CmplxKind), target, allocatable :: pint(:)
!
   real (kind=RealKind), allocatable :: sqrt_r(:)
   complex (kind=CmplxKind), allocatable :: fgr2(:)
!
   type LocalGFStruct
      integer (kind=IntKind) :: NumEs
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: NumOrbitals
      complex (kind=CmplxKind), pointer :: energy(:)
      complex (kind=CmplxKind), pointer :: gf(:)
   end type LocalGFStruct
!
   type (LocalGFStruct), allocatable :: local_gf(:)
!
   logical :: Initialized = .false.
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initLocalGF(na,lmax,ne,pola,cant,rel)
!  ===================================================================
   use RadialGridModule, only : getMaxNumRmesh
!
   use PotentialTypeModule, only : isFullPotential
!
   use AtomModule, only : getLocalNumSpecies
!
   use OrbitalBasisModule, only : initOrbitalBasis
   use OrbitalBasisModule, only : computeOrbitalBasis
!
   use InputModule, only : getKeyValue
! 
   use StringModule, only : initString, endString
   use StringModule, only : setString, getNumTokens, readToken
!
   use SortModule, only : QuickSort
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, lmax
   integer (kind=IntKind), intent(in) :: ne
   integer (kind=IntKind), intent(in) :: pola, cant, rel
   integer (kind=IntKind) :: nsize, nr, id, n, i, rstatus, is, ia
   integer (kind=IntKind), allocatable :: num_orbs(:), num_ls(:), orb_l(:,:)
!
   character (len=20) :: local_orbitals, s
   character (len=1) :: orb
!
   real (kind=RealKind) :: e_orb
!
   if (isFullPotential()) then
      call ErrorHandler('initLocalGF',                                &
                        'LocalGF module is not working for full-potential yet')
   endif
!
   lmax_kkr = lmax
   lmax_gf = 2*lmax
   kmax_kkr = (lmax_kkr+1)**2
   kmax_gf = (lmax_gf+1)**2
!
   NumLocalAtoms = na
   n_spin_pola = pola
   n_spin_cant = cant
   Relativity = rel
!
   allocate(local_gf(na))
   allocate(num_orbs(na), num_ls(na), orb_l(4,na))
!
!  ===================================================================
!  Extract the orbital information from the input.
!     so far, only works for all sites having the same local orbitals ...
!  -------------------------------------------------------------------
   rstatus = getKeyValue(1,'Default Local Orbital Labels',local_orbitals)
!  -------------------------------------------------------------------
!  rstatus = getKeyIndexValue(1,'Local Orbital Labels',local_orbitals)
!  -------------------------------------------------------------------
   call initString(80)
   call setString(adjustl(local_orbitals))
!  -------------------------------------------------------------------
   num_orbs = 0
   n = 0
   do i = 1, getNumTokens()
!     ----------------------------------------------------------------
      call readToken(i,s)
!     ----------------------------------------------------------------
      if (len_trim(s) > 1) then
         call ErrorHandler('initLocalGF','Invalid local orbital label',trim(s))
      else
         n = n + 1
         orb = trim(adjustl(s))
         if (orb == 's') then
            do id = 1, na
               orb_l(n,id) = 0
               num_orbs(id) = num_orbs(id) + 1
            enddo
         else if (orb == 'p') then
            do id = 1, na
               orb_l(n,id) = 1
               num_orbs(id) = num_orbs(id) + 3
            enddo
         else if (orb == 'd') then
            do id = 1, na
               orb_l(n,id) = 2
               num_orbs(id) = num_orbs(id) + 5
            enddo
         else if (orb == 'f') then
            do id = 1, na
               orb_l(n,id) = 3
               num_orbs(id) = num_orbs(id) + 7
            enddo
         else
            call ErrorHandler('initLocalGF','Invalid local orbital label',trim(s))
         endif
      endif
   enddo
!  -------------------------------------------------------------------
   call endString()
!  -------------------------------------------------------------------
   do id = 1, na
      num_ls(id) = n
!     ----------------------------------------------------------------
      call QuickSort(n,orb_l(1:n,id))
!     ----------------------------------------------------------------
   enddo
!
   if (BasisType == 0) then
!     ----------------------------------------------------------------
      rstatus = getKeyValue(1,'Default Local Orbital Energy',e_orb)
!     ----------------------------------------------------------------
!     rstatus = getKeyIndexValue(1,'Local Orbital Energy',e_orb)
!     ----------------------------------------------------------------
   endif
!  ===================================================================
!
!  -------------------------------------------------------------------
   call initOrbitalBasis(na,lmax,pola,cant,rel,btype=BasisType,num_ls=num_ls,orb_l=orb_l)
!  -------------------------------------------------------------------
!
   MaxNumOrbitals = 0
   do id = 1,NumLocalAtoms
      local_gf(id)%NumEs = ne
      local_gf(id)%NumSpecies = getLocalNumSpecies(id)
      local_gf(id)%NumOrbitals = num_orbs(id)
      nsize = local_gf(id)%NumOrbitals*local_gf(id)%NumOrbitals
      nsize = nsize*n_spin_cant*n_spin_cant*local_gf(id)%NumEs*local_gf(id)%NumSpecies
      allocate(local_gf(id)%energy(ne), local_gf(id)%gf(nsize))
      local_gf(id)%energy = CZERO; local_gf(id)%gf = CZERO
      MaxNumOrbitals = max(MaxNumOrbitals, num_orbs(id))
   enddo
!
   do is = 1, n_spin_pola
      do id = 1, NumLocalAtoms
!        =============================================================
!        Note: For the orbital basis type = 0, the single site solver has been
!              called inside subroutine computeOrbitalBasis.
!        =============================================================
         if (BasisType == 0) then
!           ==========================================================
!           The orbital is calculated at enery = e_orb, read from the input
!           ==========================================================
            do ia = 1, local_gf(id)%NumSpecies
!              -------------------------------------------------------
               call computeOrbitalBasis(spin=is,site=id,atom=ia,e0=e_orb)
!              -------------------------------------------------------
            enddo
         else
            do ia = 1, local_gf(id)%NumSpecies
!              -------------------------------------------------------
               call computeOrbitalBasis(spin=is,site=id,atom=ia)
!              -------------------------------------------------------
            enddo
         endif
      enddo
   enddo
!
   nr = getMaxNumRmesh()
!
   rspace_size = nr*kmax_gf
   allocate(rspace(rspace_size))
   allocate(sqrt_r(0:nr), fgr2(0:nr))
!
   fspace_size = 2*(nr+1)*kmax_gf*MaxNumOrbitals*n_spin_cant
   allocate(fspace(fspace_size))
!
   pspace_size = 2*nr*kmax_gf*MaxNumOrbitals*n_spin_cant
   allocate(pspace(pspace_size))
!
   pint_size = nr
   allocate(pint(pint_size))
!
   deallocate(num_orbs, num_ls, orb_l)
!
   Initialized = .true.
!
   end subroutine initLocalGF
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endLocalGF()
!  ===================================================================
   use OrbitalBasisModule, only : endOrbitalBasis
!
   implicit none
!
   integer (kind=IntKind) :: id
!
!  -------------------------------------------------------------------
   call endOrbitalBasis()
!  -------------------------------------------------------------------
!
   do id = 1,NumLocalAtoms
      deallocate(local_gf(id)%energy, local_gf(id)%gf)
      nullify(local_gf(id)%energy, local_gf(id)%gf)
   enddo
!
   deallocate(local_gf, rspace)
   deallocate(sqrt_r, fgr2)
   deallocate(fspace, pspace, pint)
!
   Initialized = .false.
!
   end subroutine endLocalGF
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeLocalGF(spin,ie,e)
!  ===================================================================
   use PotentialTypeModule, only : isFullPotential
!
   use AtomModule, only : getLocalNumSpecies
!
   use ScfDataModule, only : isKKR, isLSMS, isKKRCPA
!
   use MSSolverModule, only : computeMSTMatrix
!
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : getRegSolution, getIrrSolution
   use SSSolverModule, only : getJostInvMatrix
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RadialGridModule, only : getRmesh, getRadialIntegration, getNumRmesh
!
   use ClusterMatrixModule, only : getClusterKau => getKau
!
   use CrystalMatrixModule, only : getCrystalKau => getKau
!
   use CPAMediumModule, only : getImpurityMatrix
!
   use OrbitalBasisModule, only : getOrbitalBasis
   use OrbitalBasisModule, only : getLFlag
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, ie
   integer (kind=IntKind) :: site, atom
   integer (kind=IntKind) :: js, jsr, jsl, ns, iend, nsize
   integer (kind=IntKind) :: iorb, jorb, kl, klp, ir, rpow
   integer (kind=IntKind) :: m, kl_bar, norbs
   integer (kind=IntKind), pointer :: i_Lflag(:), j_Lflag(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: prod(:), prod_fp(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:), HLr(:,:,:)
   complex (kind=CmplxKind), pointer :: bf(:,:), bf_bar(:,:)
   complex (kind=CmplxKind), pointer :: kau00(:,:,:), p_kau00(:,:)
   complex (kind=CmplxKind), pointer :: gfe(:,:,:), gf(:,:,:)
   complex (kind=CmplxKind), pointer :: JostInv(:,:)
   complex (kind=CmplxKind), pointer :: xfun(:), sfun(:), pfun(:)
   complex (kind=CmplxKind), pointer :: PhiFunc(:,:,:,:)
   complex (kind=CmplxKind), pointer :: LambdaFunc(:,:,:,:)
   complex (kind=CmplxKind), pointer :: XiPhiLr(:,:,:,:)
   complex (kind=CmplxKind), pointer :: XiHLr(:,:,:,:)
   complex (kind=CmplxKind) :: Amat, Bmat, Dmat, LambdaRc, PhiRc
   complex (kind=CmplxKind) :: cfac, bfac, bfac_bar
!
   interface adjustEnergy
      function adjustEnergy_r(is,e) result(energy)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: is
         real (kind=RealKind), intent(in) :: e
         real (kind=RealKind) :: energy
      end function adjustEnergy_r
!
      function adjustEnergy_c(is,e) result(energy)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: is 
         complex (kind=CmplxKind), intent(in) :: e 
         complex (kind=CmplxKind) :: energy 
      end function adjustEnergy_c 
   end interface adjustEnergy
!
   interface
      subroutine computeCmplxProdExpan(n,lf,f,act,lg,g,lh,h)
        use KindParamModule, only : IntKind, CmplxKind
        character (len=1), intent(in) :: act
        integer (kind=IntKind), intent(in) :: n, lf, lg, lh
        complex (kind=CmplxKind), intent(in) :: f(:,:), g(:,:)
        complex (kind=CmplxKind), intent(out) :: h(:,:)
      end subroutine computeCmplxProdExpan
   end interface
!
   if (spin < 1 .or. spin > n_spin_pola) then
      call ErrorHandler('computeLocalGF','The spin index is out of range',spin)
   endif
   do site = 1, NumLocalAtoms
      if (ie < 1 .or. ie > local_gf(site)%NumEs) then
         call ErrorHandler('computeLocalGF','The energy index is out of range',ie)
      endif
   enddo
!
   do site = 1, NumLocalAtoms
      local_gf(site)%energy(ie) = e
      do js = 1, n_spin_cant
         ns = max(js,spin)
         if (isFullPotential()) then
!           ----------------------------------------------------------
            call solveSingleScattering(spin=ns,site=site,e=adjustEnergy(spin,e), &
                                       vshift=CZERO)
!           ----------------------------------------------------------
         else
!           ---------------------------------------------------------- 
            call solveSingleScattering(spin=ns,site=site,e=adjustEnergy(spin,e), &
                                       vshift=CZERO,isSphSolver=.true.,useIrrSol='H')
!           ----------------------------------------------------------
         endif
      enddo
   enddo
!  -------------------------------------------------------------------
   call computeMSTMatrix(spin,e)
!  -------------------------------------------------------------------
!
   do site = 1, NumLocalAtoms
      norbs = local_gf(site)%NumOrbitals
      r_mesh => getRmesh(site)
!
      nsize = norbs*norbs*n_spin_cant*n_spin_cant
      gfe => aliasArray3_c(local_gf(site)%gf,                         &
                           nsize,                                     &
                           local_gf(site)%NumEs,                      &
                           local_gf(site)%NumSpecies)
!
      if (isFullPotential()) then  ! The following code for the full-potential case needs some work ...
         iend = getNumRmesh(site)
      else
         iend = getNumRmesh(site,end_point='MT')
      endif
      prod_fp => aliasArray2_c(rspace,iend,kmax_gf)
!
      nsize = (iend+1)*kmax_kkr*norbs*n_spin_cant
      PhiFunc => aliasArray4_c(fspace(1:nsize),iend+1,kmax_kkr,norbs,n_spin_cant)
      LambdaFunc => aliasArray4_c(fspace(nsize+1:2*nsize),iend+1,kmax_kkr,norbs,n_spin_cant)
!
      nsize = iend*kmax_kkr*norbs*n_spin_cant
      XiPhiLr => aliasArray4_c(pspace(1:nsize),iend,kmax_kkr,norbs,n_spin_cant)
      XiHLr => aliasArray4_c(pspace(nsize+1:2*nsize),iend,kmax_kkr,norbs,n_spin_cant)
!
      do atom = 1, local_gf(site)%NumSpecies
         gf => aliasArray3_c(gfe(:,ie,atom),norbs,norbs,              &
                             n_spin_cant*n_spin_cant)
         do js = 1, n_spin_cant
            ns = max(js,spin)
            PhiLr => getRegSolution(spin=js,site=site,atom=atom)
            HLr => getIrrSolution(spin=js,site=site,atom=atom)
            do iorb = 1, norbs
!              bf => getOrbitalBasis(spin=js,site=site,atom=atom,orb=iorb,rpow=rpow,star=.false.,bfac=bfac)
               bf => getOrbitalBasis(spin=js,site=site,atom=atom,orb=iorb,rpow=rpow,bfac=bfac)
               i_Lflag => getLFlag(js,site,atom,iorb)
               if (rpow < 0 .or. rpow > 1) then
                  call ErrorHandler('computeLocalGF','rpow < 0 or rpow > 1',rpow)
               endif
               if (isFullPotential()) then  ! The following code for the full-potential case needs some work ...
!                 ==========================! ******************************************************************
!                 bf_bar => getOrbitalBasis(spin=js,site=site,atom=atom,orb=iorb,rpow=rpow,star=.true.,bfac=bfac_bar)
!!!               do kl = 1, kmax_kkr
!                    -------------------------------------------------
!!!                  call computeCmplxProdExpan(site,kmax_kkr,bf_bar,'c',kmax_kkr,PhiLr(:,:,kl),kmax_gf,prod_fp)
!                    -------------------------------------------------
!!!                  bsp(kl,iorb,js) = bfac_bar*getVolumeIntegration(site,iend,r_mesh,kmax_gf,rpow+1,prod_fp)
!
!!!                  m = mofk(kl)
!!!                  kl_bar = kl-2*m
!!!                  cfac = m1m(m)*bfac
!                    -------------------------------------------------
!!!                  call computeCmplxProdExpan(site,kmax_kkr,PhiLr(:,:,kl_bar),'s',kmax_kkr,bf,kmax_gf,prod_fp)
!                    -------------------------------------------------
!!!                  psb(kl,iorb,js) = cfac*getVolumeIntegration(site,iend,r_mesh,kmax_gf,rpow+1,prod_fp)
!
!                    -------------------------------------------------
!!!                  call computeCmplxProdExpan(site,kmax_kkr,HLr(:,:,kl_bar),'s',kmax_kkr,bf,kmax_gf,prod_fp)
!                    -------------------------------------------------
!!!                  hsb(kl,iorb,js) = cfac*getVolumeIntegration(site,iend,r_mesh,kmax_gf,rpow+1,prod_fp)
!!!               enddo
               else
!                 ====================================================
!                 In the muffin-tin case
!                 ====================================================
                  do kl = 1, kmax_kkr
                     if (i_Lflag(kl) > 0) then
                        xfun => bf(:,kl)
!
!                       ==============================================
!                       calculate Phi function
!                       ==============================================
                        sfun => PhiLr(:,kl,kl)
                        pfun => PhiFunc(:,kl,iorb,js)
!                       ----------------------------------------------
                        call integrateFuncProduct(nr=iend,r=r_mesh,f=xfun,g=sfun,rpow=1-rpow,p=pfun)
!                       ----------------------------------------------
!
!                       ==============================================
!                       calculate Lambda function
!                       ==============================================
                        sfun => HLr(:,kl,kl)
                        pfun => LambdaFunc(:,kl,iorb,js)
!                       ----------------------------------------------
                        call integrateFuncProduct(nr=iend,r=r_mesh,f=xfun,g=sfun,rpow=1-rpow,p=pfun)
!                       ----------------------------------------------
!
!                       ==============================================
!                       calculate Orbital function times Regular solution
!                       ==============================================
                        prod => XiPhiLr(:,kl,iorb,js)
                        do ir = 1, iend 
                           prod(ir) = PhiLr(ir,kl,kl)*xfun(ir)
                        enddo
!
!                       ==============================================
!                       calculate Orbital function times Irregular solution
!                       ==============================================
                        prod => XiHLr(:,kl,iorb,js)
                        do ir = 1, iend 
                           prod(ir) = HLr(ir,kl,kl)*xfun(ir)
                        enddo
                     endif
                  enddo
               endif
            enddo
         enddo
!
         if (isLSMS()) then
            kau00 => getClusterKau(local_id=site) ! Kau00 = energy * S^{-1} * [Tau00 - t_matrix] * S^{-1*}
         else if (isKKR()) then
            kau00 => getCrystalKau(local_id=site) ! Kau00 = energy * S^{-1} * [Tau00 - t_matrix] * S^{-1*}
         else if ( isKKRCPA()) then
            kau00 => getImpurityMatrix('Kau_a',site=site,atom=atom)
         endif
!
         ns = 0
         gf = CZERO
         do jsr = 1, n_spin_cant
            JostInv => getJostInvMatrix(spin=jsr,site=site,atom=atom)
            do jsl = 1, n_spin_cant
               ns = ns + 1
               p_kau00 => kau00(:,:,ns)
               if (jsr == jsl) then
                  cfac = -SQRTm1*sqrt(e)
               else
                  cfac = CZERO
               endif
               if (isFullPotential()) then  ! The following code for the full-potential case needs some work ...
                  call ErrorHandler('computeLocalGF','The full-potential local Green function is under construction')
               else
                  do jorb = 1, norbs
                     j_Lflag => getLFlag(jsl,site,atom,jorb)
                     do klp = 1, kmax_kkr
                        if (j_Lflag(klp) > 0) then
                           sfun => PhiFunc(1:iend+1,klp,jorb,jsr)
                           PhiRc = sfun(iend+1)
                           xfun => LambdaFunc(1:iend+1,klp,jorb,jsr)
                           LambdaRc = xfun(iend+1)
                           do iorb = 1, norbs
                              i_Lflag => getLFlag(jsr,site,atom,iorb)
                              do kl = 1, kmax_kkr
                                 if (i_Lflag(kl) > 0) then
!                                   ==================================
!                                   calculate A and D matrices.
!                                   ==================================
                                    prod => XiPhiLr(1:iend,kl,iorb,jsl)
!                                   ----------------------------------
                                    call integrateFuncProduct(nr=iend,r=r_mesh,f=prod,g=sfun(2:iend+1),rpow=1-rpow,p=pint)
!                                   ----------------------------------
                                    Amat = pint(iend+1) ! The index of p starts from 0 while pint start from 1
!                                   ----------------------------------
                                    call integrateFuncProduct(nr=iend,r=r_mesh,f=prod,g=xfun(2:iend+1),rpow=1-rpow,p=pint)
!                                   ----------------------------------
                                    Dmat = pint(iend+1) ! The index of p starts from 0 while pint starts from 1
!                                   ==================================
!                                   calculate B matrix.
!                                   ==================================
                                    prod => XiHLr(1:iend,kl,iorb,jsl)
!                                   ----------------------------------
                                    call integrateFuncProduct(nr=iend,r=r_mesh,f=prod,g=sfun,rpow=1-rpow,p=pint)
!                                   ----------------------------------
                                    Bmat = pint(iend+1) ! The index of p starts from 0 while pint starts from 1
!
!                                   ==================================
!                                   calculate local Green function
!                                   ==================================
                                    gf(iorb,jorb,ns) = gf(iorb,jorb,ns)                            &
                                             + e*Amat*(p_kau00(klp,kl)-p_kau00(kl,klp))            &
                                             + e*PhiFunc(iend+1,kl,iorb,jsl)*p_kau00(kl,klp)*PhiRc &
                                             + cfac*(Bmat-Dmat+PhiFunc(iend+1,kl,iorb,jsl)*LambdaRc)*JostInv(kl,klp)
                                 endif
                              enddo
                           enddo
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
   enddo
!
   nullify(prod, prod_fp, PhiLr, HLr, bf, bf_bar, i_Lflag, j_Lflag)
   nullify(kau00, p_kau00, gfe, gf, sfun, xfun, XiPhiLr, XiHLr, JostInv)
!
   end subroutine computeLocalGF
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalGF_v0(site,atom,ie,jsr,jsl,jorb,iorb) result (gf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, ie, jsr, jsl, jorb, iorb
   integer (kind=IntKind) :: nsize, n
!
   complex (kind=CmplxKind) :: gf
   complex (kind=CmplxKind), pointer :: gfe(:,:,:,:,:)
!
   nsize = local_gf(site)%NumOrbitals*local_gf(site)%NumOrbitals
   gfe => aliasArray5_c(local_gf(site)%gf,                            &
                        nsize,                                        &
                        n_spin_cant,                                  &
                        n_spin_cant,                                  &
                        local_gf(site)%NumEs,                         &
                        local_gf(site)%NumSpecies)
!
   n = (jorb-1)*local_gf(site)%NumOrbitals + iorb
   gf = gfe(n,jsl,jsr,ie,atom)
!
   nullify(gfe)
!
   end function getLocalGF_v0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalGF_v2(site,atom,ie,jsr,jsl) result (gf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, ie, jsr, jsl
   integer (kind=IntKind) :: nsize
!
   complex (kind=CmplxKind), pointer :: gf(:,:)
   complex (kind=CmplxKind), pointer :: gfe(:,:,:,:,:)
!
   nsize = local_gf(site)%NumOrbitals*local_gf(site)%NumOrbitals
   gfe => aliasArray5_c(local_gf(site)%gf,                            &
                        nsize,                                        &
                        n_spin_cant,                                  &
                        n_spin_cant,                                  &
                        local_gf(site)%NumEs,                         &
                        local_gf(site)%NumSpecies)
!
   gf => aliasArray2_c(gfe(:,jsl,jsr,ie,atom),                        &
                       local_gf(site)%NumOrbitals,                    &
                       local_gf(site)%NumOrbitals)
!
   nullify(gfe)
!
   end function getLocalGF_v2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalGF_v4(site, atom, ie) result (gf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, ie
   integer (kind=IntKind) :: nsize
!
   complex (kind=CmplxKind), pointer :: gf(:,:,:,:)
   complex (kind=CmplxKind), pointer :: gfe(:,:,:)
!
   nsize = local_gf(site)%NumOrbitals*local_gf(site)%NumOrbitals
   gfe => aliasArray3_c(local_gf(site)%gf,                            &
                        nsize*n_spin_cant*n_spin_cant,                &
                        local_gf(site)%NumEs,                         &
                        local_gf(site)%NumSpecies)
!   
   gf => aliasArray4_c(gfe(:,ie,atom),                                &
                       local_gf(site)%NumOrbitals,                    &
                       local_gf(site)%NumOrbitals,                    &
                       n_spin_cant, n_spin_cant)
!
   nullify(gfe)
!
   end function getLocalGF_v4
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalGF_v5(site, atom) result (gf)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom
   integer (kind=IntKind) :: nsize
!
   complex (kind=CmplxKind), pointer :: gf(:,:,:,:,:)
   complex (kind=CmplxKind), pointer :: gfe(:,:)
!
   nsize = local_gf(site)%NumOrbitals*local_gf(site)%NumOrbitals
   nsize = nsize*n_spin_cant*n_spin_cant*local_gf(site)%NumEs
   gfe => aliasArray2_c(local_gf(site)%gf,                            &
                        nsize, local_gf(site)%NumSpecies)
!
   gf => aliasArray5_c(gfe(:,atom),                                   &
                       local_gf(site)%NumOrbitals,                    &
                       local_gf(site)%NumOrbitals,                    &
                       n_spin_cant, n_spin_cant, local_gf(site)%NumEs)
!
   nullify(gfe)
!
   end function getLocalGF_v5
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumOrbitals(site) result (n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind) :: n
!
   n = local_gf(site)%NumOrbitals
!
   end function getNumOrbitals
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEnergyMesh(site, NumEs) result (e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site
   integer (kind=IntKind), intent(out) :: NumEs
!
   complex (kind=CmplxKind), pointer :: e(:)
!
   NumEs = local_gf(site)%NumEs
   e => local_gf(site)%energy
!
   end function getEnergyMesh
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printLocalGF(site, atom, ie, jsr, jsl)
!  ===================================================================
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, ie, jsr, jsl
   integer (kind=IntKind) :: nsize, norbs
!
   complex (kind=CmplxKind), pointer :: gf(:,:)
   complex (kind=CmplxKind), pointer :: gfe(:,:,:,:,:)
!
   nsize = local_gf(site)%NumOrbitals*local_gf(site)%NumOrbitals
   gfe => aliasArray5_c(local_gf(site)%gf,                            &
                        nsize,                                        &
                        n_spin_cant,                                  &
                        n_spin_cant,                                  &
                        local_gf(site)%NumEs,                         &
                        local_gf(site)%NumSpecies)
!
   norbs = local_gf(site)%NumOrbitals
   gf => aliasArray2_c(gfe(:,jsl,jsr,ie,atom),norbs,norbs)
!
!  -------------------------------------------------------------------
   call writeMatrix('Local GF',gf,norbs,norbs,TEN2m6)
!  -------------------------------------------------------------------
!
   nullify(gfe, gf)
!
   end subroutine printLocalGF
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine integrateFuncProduct(nr, r, f, g, rpow, p)
!  ===================================================================
!  Given functions f(r) and g(r), calculate the following integration
!      p(r) = int^r_0 f(r')*g(r')*r'^rpow dr'
!  *******************************************************************
   use InterpolationModule, only : FitInterp
!
   use IntegrationModule, only : calIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(nr)
!
   complex (kind=CmplxKind), intent(in) :: f(nr), g(nr)
   complex (kind=CmplxKind), intent(out) :: p(0:nr) ! Note: p array index starts from 0
   complex (kind=CmplxKind) :: dps
!
   sqrt_r(0) = ZERO
   do ir = 1, nr
      sqrt_r(ir) = sqrt(r(ir))
   enddo
!
   do ir = 1, nr
      fgr2(ir) = TWO*f(ir)*g(ir)*r(ir)**rpow
   enddo
!
!  -------------------------------------------------------------------
   call FitInterp( 4, sqrt_r(1:4), fgr2(1:4), ZERO, fgr2(0), dps )
!  -------------------------------------------------------------------
   call calIntegration(nr+1,sqrt_r(0:),fgr2(0:),p(0:),1)
!  -------------------------------------------------------------------
!
   end subroutine integrateFuncProduct
!  ===================================================================
end module LocalGFModule
