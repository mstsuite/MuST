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
   integer (kind=IntKind) :: pspace_size
!
   complex (kind=CmplxKind), target, allocatable :: rspace(:)
   complex (kind=CmplxKind), target, allocatable :: pspace(:)
   complex (kind=CmplxKind), pointer :: bsp(:,:,:)
   complex (kind=CmplxKind), pointer :: psb(:,:,:)
   complex (kind=CmplxKind), pointer :: hsb(:,:,:)
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
   subroutine initLocalGF(na,lmax,ne,norb,pola,cant,rel)
!  ===================================================================
   use RadialGridModule, only : getMaxNumRmesh
!
   use PotentialTypeModule, only : isFullPotential
!
   use AtomModule, only : getLocalNumSpecies
!
   use OrbitalBasisModule, only : initOrbitalBasis
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, lmax
   integer (kind=IntKind), intent(in) :: ne, norb(na)
   integer (kind=IntKind), intent(in) :: pola, cant, rel
   integer (kind=IntKind) :: nsize, nr, id, MaxNumOrbitals
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
!
   MaxNumOrbitals = 0
   do id = 1,NumLocalAtoms
      local_gf(id)%NumEs = ne
      local_gf(id)%NumSpecies = getLocalNumSpecies(id)
      local_gf(id)%NumOrbitals = norb(id)
      nsize = local_gf(id)%NumOrbitals*local_gf(id)%NumOrbitals
      nsize = nsize*n_spin_cant*n_spin_cant*local_gf(id)%NumEs*local_gf(id)%NumSpecies
      allocate(local_gf(id)%energy(ne), local_gf(id)%gf(nsize))
      local_gf(id)%energy = CZERO; local_gf(id)%gf = CZERO
      MaxNumOrbitals = max(MaxNumOrbitals, norb(id))
   enddo
!
   nsize = kmax_kkr*MaxNumOrbitals*n_spin_cant
   pspace_size = 3*nsize
   allocate(pspace(pspace_size))
   psb => aliasArray3_c(pspace,kmax_kkr,MaxNumOrbitals,n_spin_cant)
   hsb => aliasArray3_c(pspace(nsize+1:2*nsize),kmax_kkr,MaxNumOrbitals,n_spin_cant)
   if (isFullPotential()) then
      bsp => aliasArray3_c(pspace(2*nsize+1:3*nsize),kmax_kkr,MaxNumOrbitals,n_spin_cant)
   else
      bsp => psb
   endif
!
   rspace_size = getMaxNumRmesh()*kmax_gf
   allocate(rspace(rspace_size))
!
!  -------------------------------------------------------------------
   call initOrbitalBasis(na,lmax,pola,cant,rel,btype=0)
!  -------------------------------------------------------------------
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
   nullify(psb, hsb, bsp)
   deallocate(local_gf, rspace, pspace)
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
   use OrbitalBasisModule, only : computeOrbitalBasis, getOrbitalBasis
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, ie
   integer (kind=IntKind) :: site, atom
   integer (kind=IntKind) :: js, jsr, jsl, ns, iend, nsize
   integer (kind=IntKind) :: iorb, jorb, kl, klp, ir, rpow
   integer (kind=IntKind) :: m, kl_bar
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), pointer :: prod(:), prod_fp(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:), HLr(:,:,:)
   complex (kind=CmplxKind), pointer :: bf(:,:), bf_bar(:,:)
   complex (kind=CmplxKind), pointer :: kau00(:,:,:), p_kau00(:,:)
   complex (kind=CmplxKind), pointer :: gfe(:,:,:), gf(:,:,:)
   complex (kind=CmplxKind), pointer :: pk(:,:), pt(:,:), JostInv(:,:)
   complex (kind=CmplxKind) :: cfac
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
!     ================================================================
!     Note: For the orbital basis type = 0, the single site solver has been
!           called inside subroutine computeOrbitalBasis.
!     ================================================================
      do atom = 1, local_gf(site)%NumSpecies
!        -------------------------------------------------------------
         call computeOrbitalBasis(spin=spin,site=site,atom=atom,      &
                                  norbs=local_gf(site)%NumOrbitals,e=e)
!        -------------------------------------------------------------
      enddo
   enddo
!  -------------------------------------------------------------------
   call computeMSTMatrix(spin,e)
!  -------------------------------------------------------------------
!
   do site = 1, NumLocalAtoms
      nsize = local_gf(site)%NumOrbitals*local_gf(site)%NumOrbitals*n_spin_cant*n_spin_cant
      gfe => aliasArray3_c(local_gf(site)%gf,                         &
                           nsize,                                     &
                           local_gf(site)%NumEs,                      &
                           local_gf(site)%NumSpecies)
      r_mesh => getRmesh(site)
      if (isFullPotential()) then  ! The following code for the full-potential case needs some work ...
         iend = getNumRmesh(site)
      else
         iend = getNumRmesh(site,end_point='MT')
      endif
      prod_fp => aliasArray2_c(rspace,iend,kmax_gf)
      prod => prod_fp(:,1)
!
      nsize = kmax_kkr*local_gf(site)%NumOrbitals
      pk => aliasArray2_c(rspace(1:nsize),kmax_kkr,local_gf(site)%NumOrbitals)
      pt => aliasArray2_c(rspace(nsize+1:2*nsize),kmax_kkr,local_gf(site)%NumOrbitals)
!
      do atom = 1, local_gf(site)%NumSpecies
         gf => aliasArray3_c(gfe(:,ie,atom),                          &
                             local_gf(site)%NumOrbitals,              &
                             local_gf(site)%NumOrbitals,              &
                             n_spin_cant*n_spin_cant)
         do js = 1, n_spin_cant
            ns = max(js,spin)
            PhiLr => getRegSolution(spin=js,site=site,atom=atom)
            HLr => getIrrSolution(spin=js,site=site,atom=atom)
            do iorb = 1, local_gf(site)%NumOrbitals
               bf => getOrbitalBasis(spin=js,site=site,atom=atom,orb=iorb,rpow=rpow)
               if (rpow < 0 .or. rpow > 1) then
                  call ErrorHandler('computeLocalGF','rpow < 0 or rpow > 1',rpow)
               endif
               if (isFullPotential()) then  ! The following code for the full-potential case needs some work ...
                  bf_bar => getOrbitalBasis(spin=js,site=site,atom=atom,orb=iorb,rpow=rpow,star=.true.)
                  do kl = 1, kmax_kkr
!                    -------------------------------------------------
                     call computeCmplxProdExpan(site,kmax_kkr,bf_bar,'c',kmax_kkr,PhiLr(:,:,kl),kmax_gf,prod_fp)
!                    -------------------------------------------------
                     bsp(kl,iorb,js) = getVolumeIntegration(site,iend,r_mesh,kmax_gf,rpow+1,prod_fp)
!
                     m = mofk(kl)
                     kl_bar = kl-2*m
                     cfac = m1m(m)
!                    -------------------------------------------------
                     call computeCmplxProdExpan(site,kmax_kkr,PhiLr(:,:,kl_bar),'s',kmax_kkr,bf,kmax_gf,prod_fp)
!                    -------------------------------------------------
                     psb(kl,iorb,js) = cfac*getVolumeIntegration(site,iend,r_mesh,kmax_gf,rpow+1,prod_fp)
!
!                    -------------------------------------------------
                     call computeCmplxProdExpan(site,kmax_kkr,HLr(:,:,kl_bar),'s',kmax_kkr,bf,kmax_gf,prod_fp)
!                    -------------------------------------------------
                     hsb(kl,iorb,js) = cfac*getVolumeIntegration(site,iend,r_mesh,kmax_gf,rpow+1,prod_fp)
                  enddo
               else
                  do kl = 1, kmax_kkr
                     do ir = 1, iend 
                        prod(ir) = PhiLr(ir,kl,kl)*bf(ir,kl)
                     enddo
                     if (rpow == 0) then
                        psb(kl,iorb,js) = getRadialIntegration(site,iend,prod,rpow=1)
                     else
                        psb(kl,iorb,js) = getRadialIntegration(site,iend,prod)
                     endif
                     do ir = 1, iend 
                        prod(ir) = HLr(ir,kl,kl)*bf(ir,kl)
                     enddo
                     if (rpow == 0) then
                        hsb(kl,iorb,js) = getRadialIntegration(site,iend,prod,rpow=1)
                     else
                        hsb(kl,iorb,js) = getRadialIntegration(site,iend,prod)
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
            pt = CZERO
            do iorb = 1, local_gf(site)%NumOrbitals
               do kl = 1, kmax_kkr
                  do klp = 1, kmax_kkr
                     pt(kl,iorb) = pt(kl,iorb) + bsp(klp,iorb,jsr)*JostInv(klp,kl)
                  enddo
               enddo
            enddo
            do jsl = 1, n_spin_cant
               ns = ns + 1
               p_kau00 => kau00(:,:,ns)
               pk = CZERO
               if (jsr == jsl) then
                  cfac = -SQRTm1*sqrt(e)
               else
                  cfac = CZERO
               endif
               do iorb = 1, local_gf(site)%NumOrbitals
                  do kl = 1, kmax_kkr
                     do klp = 1, kmax_kkr
                        pk(kl,iorb) = pk(kl,iorb) + bsp(klp,iorb,jsl)*p_kau00(klp,kl)
                     enddo
                  enddo
               enddo
               do jorb = 1, local_gf(site)%NumOrbitals
                  do iorb = 1, local_gf(site)%NumOrbitals
                     do kl = 1, kmax_kkr
                        gf(iorb,jorb,ns) = gf(iorb,jorb,ns) +               &
                                           pk(kl,iorb)*psb(kl,jorb,jsr) +    &
                                           cfac*pt(kl,iorb)*hsb(kl,jorb,jsr)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
!
   nullify(prod, prod_fp, PhiLr, HLr, bf, bf_bar)
   nullify(kau00, p_kau00, gfe, gf, pk, pt, JostInv)
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
!  -------------------------------------------------------------------
   call writeMatrix('Local GF',gf,kmax_kkr,kmax_kkr,TEN2m6)
!  -------------------------------------------------------------------
!
   nullify(gfe, gf)
!
   end subroutine printLocalGF
!  ===================================================================
end module LocalGFModule
