module OrbitalBasisModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
   use MPPModule, only : MyPE
!
public :: initOrbitalBasis,    &
          endOrbitalBasis,     &
          computeOrbitalBasis, &
          getOrbitalBasis,     &
          getLFlag,            &
          getBasisType
!
   interface getOrbitalBasis
      module procedure getOrbitalBasis_mt, getOrbitalBasis_fp
   end interface getOrbitalBasis
!
private
   integer (kind=IntKind) :: NumLocalAtoms
   integer (kind=IntKind) :: lmax_kkr
   integer (kind=IntKind) :: kmax_kkr
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: Relativity
   integer (kind=IntKind) :: MaxNumSpecies
   integer (kind=IntKind) :: MaxNumRs
   integer (kind=IntKind) :: BasisType
!
   complex (kind=CmplxKind) :: BasisFactor   ! This is a prefactor for the basis function
!
   type LocalOrbitalStruct
      integer (kind=IntKind) :: NumOrbitals
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind) :: r_power
      integer (kind=IntKind), pointer :: lm_flag(:,:) ! For a given orbital (2nd dimension index), it is
                                                      ! 1 if the kl index (1st dimension index) of the 
                                                      ! corresponding orbital in spherical harmonics
                                                      ! expansion has non-zero component, and is 0 if
                                                      ! the expansion has zero component due to symmetry
      integer (kind=IntKind), pointer :: orb2kl(:)    ! It relates each orbital to a kl index
!     complex (kind=CmplxKind), pointer :: orbital_star(:,:,:) ! for real energy, orbital_star is just
                                                               ! the complex conjugate of orbital
      complex (kind=CmplxKind), pointer :: orbital_mt(:,:)     ! If the orbital has spherical symmetry, obtained
                                                               ! from the single site solution of a muffin-tin
                                                               ! potential
      complex (kind=CmplxKind), pointer :: orbital_fp(:,:,:)   ! If the orbital is based on the single site
                                                               ! solution obtained from a full-potential
   end type LocalOrbitalStruct
!
   type (LocalOrbitalStruct), allocatable :: LocalOrbitals(:,:,:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initOrbitalBasis(na,lmax,pola,cant,rel,btype,num_ls,orb_l)
!  ===================================================================
   use OrthoNormModule, only : initOrthoNorm
   use AtomModule, only : getLocalNumSpecies
   use RadialGridModule, only : getNumRmesh
   use PotentialTypeModule, only : isFullPotential
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, lmax, pola, cant, rel, btype
   integer (kind=IntKind), intent(in) :: num_ls(na), orb_l(4,na)
   integer (kind=IntKind) :: id, ia, norbs, l, m, n, js, NumRs
   integer (kind=IntKind) :: orb2kl(25)
!
   NumLocalAtoms = na
   lmax_kkr = lmax
   kmax_kkr = (lmax+1)**2
   n_spin_pola = pola
   n_spin_cant = cant
   Relativity = rel
   BasisType = btype ! = 0, Using the local regular solution as the basis function.
                     ! = else, to be determined.
!
   if (btype /= 0) then ! Using the local regular solution as the basis function.
      call ErrorHandler('initOrbitalBasis','Invalid basis function type',btype)
   endif
!
   MaxNumSpecies = 0
   do id = 1, NumLocalAtoms
      MaxNumSpecies = max(MaxNumSpecies,getLocalNumSpecies(id))
   enddo
!
   allocate(LocalOrbitals(n_spin_cant,MaxNumSpecies,NumLocalAtoms))
   MaxNumRs = 0
   do id = 1, NumLocalAtoms
      if (isFullPotential()) then
         NumRs = getNumRmesh(id)
      else
         NumRs = getNumRmesh(id,end_point='MT')
      endif
      MaxNumRs = max(MaxNumRs,NumRs)
      norbs = 0
      do n = 1, num_ls(id)
         l = orb_l(n,id)
         do m = -l, l
            norbs = norbs + 1
            orb2kl(norbs) = (l+1)**2-l+m
         enddo
      enddo
      do ia = 1, getLocalNumSpecies(id)
         do js = 1, n_spin_cant
            LocalOrbitals(js,ia,id)%NumOrbitals = norbs
            LocalOrbitals(js,ia,id)%NumRs = NumRs
            if (BasisType == 0) then
               LocalOrbitals(js,ia,id)%r_power = 1
            else ! This needs to be re-examined..........
               LocalOrbitals(js,ia,id)%r_power = 0
            endif
            allocate(LocalOrbitals(js,ia,id)%lm_flag(kmax_kkr,norbs))
            LocalOrbitals(js,ia,id)%lm_flag = 0
            if (isFullPotential()) then
               allocate(LocalOrbitals(js,ia,id)%orbital_fp(NumRs,kmax_kkr,norbs))
               nullify(LocalOrbitals(js,ia,id)%orbital_mt)
            else
               allocate(LocalOrbitals(js,ia,id)%orbital_mt(NumRs,norbs))
               nullify(LocalOrbitals(js,ia,id)%orbital_fp)
            endif
!           allocate(LocalOrbitals(js,ia,id)%orbital_star(NumRs,kmax_kkr,norbs))
            if (btype == 0) then
               allocate(LocalOrbitals(js,ia,id)%orb2kl(norbs))
               LocalOrbitals(js,ia,id)%orb2kl(1:norbs) = orb2kl(1:norbs)
            else
               nullify(LocalOrbitals(js,ia,id)%orb2kl)
            endif
         enddo
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call initOrthoNorm(MaxNumRs,lmax_kkr,isFullPotential())
!  -------------------------------------------------------------------
!
   end subroutine initOrbitalBasis
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endOrbitalBasis()
!  ===================================================================
   use OrthoNormModule, only : endOrthoNorm
   use PotentialTypeModule, only : isFullPotential
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   integer (kind=IntKind) :: id, ia, js
!
   do id = 1, NumLocalAtoms
      do ia = 1, getLocalNumSpecies(id)
         do js = 1, n_spin_cant
            if (associated(LocalOrbitals(js,ia,id)%lm_flag)) then
               deallocate(LocalOrbitals(js,ia,id)%lm_flag)
            endif
            nullify(LocalOrbitals(js,ia,id)%lm_flag)
            if (associated(LocalOrbitals(js,ia,id)%orb2kl)) then
               deallocate(LocalOrbitals(js,ia,id)%orb2kl)
            endif
            nullify(LocalOrbitals(js,ia,id)%orb2kl)
            if (isFullPotential()) then
               if (associated(LocalOrbitals(js,ia,id)%orbital_fp)) then
                  deallocate(LocalOrbitals(js,ia,id)%orbital_fp)
               endif
               nullify(LocalOrbitals(js,ia,id)%orbital_fp)
            else
               if (associated(LocalOrbitals(js,ia,id)%orbital_mt)) then
                  deallocate(LocalOrbitals(js,ia,id)%orbital_mt)
               endif
               nullify(LocalOrbitals(js,ia,id)%orbital_mt)
            endif
!           if (associated(LocalOrbitals(js,ia,id)%orbital_star)) then
!              deallocate(LocalOrbitals(js,ia,id)%orbital_star)
!           endif
!           nullify(LocalOrbitals(js,ia,id)%orbital_star)
         enddo
      enddo
   enddo
!
   if (isFullPotential()) then
      deallocate(LocalOrbitals)
   else
      deallocate(LocalOrbitals)
   endif
!
!  -------------------------------------------------------------------
   call endOrthoNorm()
!  -------------------------------------------------------------------
!
   end subroutine endOrbitalBasis
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeOrbitalBasis(spin,site,atom,e0)
!  ===================================================================
!
!  Note: spin = 1 or n_spin_pola/n_spin_cant
!
!  *******************************************************************
   use MathParamModule, only : CZERO, ZERO, ONE, TWO, PI, CONE, TEN2m8
!
   use IntegerFactorsModule, only : m1m, mofk
!
   use OrthoNormModule, only : orthogonalize, normalize
!
   use PotentialTypeModule, only : isFullPotential
!
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RadialGridModule, only : getRmesh
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin,site,atom
   integer (kind=IntKind) :: kl, klp, ir
!  integer (kind=IntKind) :: m, mp, kl_bar, klp_bar
   integer (kind=IntKind) :: js, ns, NumRs, iorb, norbs
   integer (kind=IntKind) :: info(2)
!
   real (kind=RealKind), intent(in), optional :: e0
   real (kind=RealKind) :: fnorm, fnorm_ws, fnorm_mt
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: f(:,:,:)
   complex (kind=CmplxKind) :: e
!  complex (kind=CmplxKind) :: cfac
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
   if (BasisType == 0) then ! Using the local regular solution as the basis function.
      if (.not.present(e0)) then
         call ErrorHandler('computeOrbitalBasis',                        &
                           'Energy parameter is required for this basis type.')
      else
         e = e0
      endif
   endif
!
   do js = 1, n_spin_cant
      if (LocalOrbitals(js,atom,site)%NumOrbitals == 0) then
         call ErrorHandler('computeOrbitalBasis','NumOrbitals = 0')
      endif
   enddo
!
   r_mesh => getRmesh(site)
!
   if (BasisType == 0) then ! Using the local regular solution as the basis function.
!     BasisFactor = sqrt(TWO*e/PI) ! This needs to be checked...
      BasisFactor = CONE
      if (isFullPotential()) then
         do js = 1, n_spin_cant
            ns = max(js,spin)
!           ----------------------------------------------------------
            call solveSingleScattering(spin=ns,site=site,e=adjustEnergy(spin,e), &
                                       vshift=CZERO)
!           ----------------------------------------------------------
         enddo
!        Needs to determine LocalOrbitals(js,atom,site)%lm_flags here ...
      else
         do js = 1, n_spin_cant
            norbs = LocalOrbitals(js,atom,site)%NumOrbitals
            ns = max(js,spin)
!           ----------------------------------------------------------
            call solveSingleScattering(spin=ns,site=site,e=adjustEnergy(spin,e), &
                                       vshift=CZERO,isSphSolver=.true.,useIrrSol='H')
!           ----------------------------------------------------------
            do iorb = 1, norbs
               kl = LocalOrbitals(js,atom,site)%orb2kl(iorb)
               LocalOrbitals(js,atom,site)%lm_flag(kl,iorb) = 1
            enddo
         enddo
      endif
      do js = 1, n_spin_cant
         norbs = LocalOrbitals(js,atom,site)%NumOrbitals
         NumRs = LocalOrbitals(js,atom,site)%NumRs
         f => getRegSolution(spin=js,site=site,atom=atom)
!
         if (isFullPotential()) then
            LocalOrbitals(js,atom,site)%orbital_fp = CZERO
         else
            LocalOrbitals(js,atom,site)%orbital_mt = CZERO
         endif
!        LocalOrbitals(js,atom,site)%orbital_star = CZERO
         do iorb = 1, norbs
            kl = LocalOrbitals(js,atom,site)%orb2kl(iorb)
            if (isFullPotential()) then
               do klp = 1, kmax_kkr
                  if (LocalOrbitals(js,atom,site)%lm_flag(klp,iorb) > 0) then
                     do ir = 1, NumRs
                         LocalOrbitals(js,atom,site)%orbital_fp(ir,klp,iorb) = f(ir,klp,kl)
                     enddo
                  endif
               enddo
            else
               do ir = 1, NumRs
                  LocalOrbitals(js,atom,site)%orbital_mt(ir,iorb) = f(ir,kl,kl)
               enddo
            endif
!
!           m = mofk(kl)
!           kl_bar = kl - 2*m
!           do klp = 1, kmax_kkr
!              if (LocalOrbitals(js,atom,site)%lm_flag(klp,iorb) > 0) then
!                 mp = mofk(klp)
!                 cfac = m1m(mp+m)
!                 klp_bar = klp - 2*mp
!                 do ir = 1, NumRs
!                    LocalOrbitals(js,atom,site)%orbital_star(ir,klp,iorb) = cfac*f(ir,klp_bar,kl_bar)
!                 enddo
!              endif
!           enddo
         enddo
!        =============================================================
!        Orthogonalize and normalize the orbitals
!        =============================================================
         info(1) = site; info(2) = 2*LocalOrbitals(js,atom,site)%r_power
         if (isFullPotential()) then
!           ----------------------------------------------------------
            call orthogonalize(info,NumRs,lmax_kkr,norbs,r_mesh, &
                               LocalOrbitals(js,atom,site)%orbital_fp)
!           ----------------------------------------------------------
            do iorb = 1, norbs
!              -------------------------------------------------------
               call normalize(info,NumRs,lmax_kkr,r_mesh,LocalOrbitals(js,atom,site)%orbital_fp(:,:,iorb))
!              -------------------------------------------------------
            enddo
!           ----------------------------------------------------------
         else ! In the muffin-tin, case, the orbital is already orthogonalized due to spherical harmonics.
              ! We only need to normalize the orbitals
            do iorb = 1, norbs
!              -------------------------------------------------------
               call normalize(info,NumRs,r_mesh,LocalOrbitals(js,atom,site)%orbital_mt(:,iorb))
!              -------------------------------------------------------
            enddo
         endif
      enddo
   else
!     ----------------------------------------------------------------
      call ErrorHandler('computeOrbitalBasis','Invalid basis function type',BasisType)
!     ----------------------------------------------------------------
      BasisFactor = CONE
   endif
!
   nullify(f)
!
   end subroutine computeOrbitalBasis
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  function getOrbitalBasis(spin,site,atom,orb,rpow,star,bfac) result(f)
   function getOrbitalBasis_fp(spin,site,atom,orb,rpow,bfac) result(f)
!  ===================================================================
!
!  Note: spin = 1 or n_spin_cant
!
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, site, atom, orb
   integer (kind=IntKind), intent(out) :: rpow
!
!  logical, intent(in) :: star
!
   complex (kind=CmplxKind), intent(out) :: bfac
   complex (kind=CmplxKind), pointer :: f(:,:) ! returns the basis function multiplied by r^rpow
!
   if (LocalOrbitals(spin,atom,site)%NumOrbitals == 0) then
      call ErrorHandler('getOrbitalBasis','The orbital basis functions are not available')
   endif
!
   f => LocalOrbitals(spin,atom,site)%orbital_fp(:,:,orb)
!
!  if (star) then
!     f => LocalOrbitals(spin,atom,site)%orbital_star(:,:,orb)
!  endif
!
   rpow = LocalOrbitals(spin,atom,site)%r_power
!
   bfac = BasisFactor
!
   end function getOrbitalBasis_fp
!  ===================================================================
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  function getOrbitalBasis(spin,site,atom,orb,rpow,star,bfac) result(f)
   function getOrbitalBasis_mt(spin,site,atom,orb,kl,rpow,bfac) result(f)
!  ===================================================================
!
!  Note: spin = 1 or n_spin_cant
!
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, site, atom, orb
   integer (kind=IntKind), intent(out) :: kl, rpow
!
!  logical, intent(in) :: star
!
   complex (kind=CmplxKind), intent(out) :: bfac
   complex (kind=CmplxKind), pointer :: f(:) ! returns the basis function multiplied by r^rpow
!
   if (LocalOrbitals(spin,atom,site)%NumOrbitals == 0) then
      call ErrorHandler('getOrbitalBasis','The orbital basis functions are not available')
   endif
!
   f => LocalOrbitals(spin,atom,site)%orbital_mt(:,orb)
!
!  if (star) then
!     f => LocalOrbitals(spin,atom,site)%orbital_star(:,orb)
!  endif
!
   kl = LocalOrbitals(spin,atom,site)%orb2kl(orb)
   rpow = LocalOrbitals(spin,atom,site)%r_power
!
   bfac = BasisFactor
!
   end function getOrbitalBasis_mt
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLFlag(spin,site,atom,orb) result(lm_flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, site, atom, orb
   integer (kind=IntKind), pointer :: lm_flag(:)
!
   if (LocalOrbitals(spin,atom,site)%NumOrbitals == 0) then
      call ErrorHandler('getOrbitalBasis','The orbital basis functions are not available')
   endif
!
   lm_flag => LocalOrbitals(spin,atom,site)%lm_flag(:,orb)
!
   end function getLFlag
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBasisType() result(t)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: t
!
   t = BasisType
!
   end function getBasisType
!  ===================================================================
end module OrbitalBasisModule
