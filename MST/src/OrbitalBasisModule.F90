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
private
   integer (kind=IntKind) :: NumLocalAtoms
   integer (kind=IntKind) :: kmax_kkr
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: Relativity
   integer (kind=IntKind) :: MaxNumSpecies
   integer (kind=IntKind) :: BasisType
!
   complex (kind=CmplxKind) :: BasisFactor   ! This is a prefactor for the basis function
!
   type LocalOrbitalStruct
      integer (kind=IntKind) :: NumOrbitals
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind) :: r_power
      integer (kind=IntKind), pointer :: lm_flag(:,:)
      integer (kind=IntKind), pointer :: orb2kl(:)
      complex (kind=CmplxKind), pointer :: orbital_star(:)
      complex (kind=CmplxKind), pointer :: orbital(:,:,:)
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
   kmax_kkr = (lmax+1)**2
   n_spin_pola = pola
   n_spin_cant = cant
   Relativity = rel
   BasisType = btype
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
   do id = 1, NumLocalAtoms
      if (isFullPotential()) then
         NumRs = getNumRmesh(id)
      else
         NumRs = getNumRmesh(id,end_point='MT')
      endif
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
            allocate(LocalOrbitals(js,ia,id)%orbital_star(NumRs*kmax_kkr*norbs))
            nullify(LocalOrbitals(js,ia,id)%orbital)
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
   end subroutine initOrbitalBasis
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endOrbitalBasis()
!  ===================================================================
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
            if (associated(LocalOrbitals(js,ia,id)%orbital_star)) then
               deallocate(LocalOrbitals(js,ia,id)%orbital_star)
            endif
            nullify(LocalOrbitals(js,ia,id)%orbital_star)
            nullify(LocalOrbitals(js,ia,id)%orbital)
         enddo
      enddo
   enddo
   deallocate(LocalOrbitals)
!
   end subroutine endOrbitalBasis
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeOrbitalBasis(spin,site,atom,e)
!  ===================================================================
!
!  Note: spin = 1 or n_spin_pola/n_spin_cant
!
!  *******************************************************************
   use MathParamModule, only : CZERO, TWO, PI, CONE
!
   use IntegerFactorsModule, only : m1m, mofk
!
   use PotentialTypeModule, only : isFullPotential
!
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin,site,atom
   integer (kind=IntKind) :: kl, klp, ir, m, mp, kl_bar, klp_bar
   integer (kind=IntKind) :: js, ns, osize, NumRs, iorb, norbs
!
   complex (kind=CmplxKind), intent(in), optional :: e
   complex (kind=CmplxKind), pointer :: f_star(:,:,:)
   complex (kind=CmplxKind) :: cfac
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
      if (.not.present(e)) then
         call ErrorHandler('computeOrbitalBasis',                        &
                           'Energy parameter is required for this basis type.')
      endif
   endif
!
   do js = 1, n_spin_cant
      if (LocalOrbitals(js,atom,site)%NumOrbitals == 0) then
         call ErrorHandler('computeOrbitalBasis','NumOrbitals = 0')
      endif
   enddo
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
         LocalOrbitals(js,atom,site)%orbital => getRegSolution(spin=js,site=site,atom=atom)
         NumRs = LocalOrbitals(js,atom,site)%NumRs
         f_star => aliasArray3_c(LocalOrbitals(js,atom,site)%orbital_star,NumRs,kmax_kkr,norbs)
         f_star = CZERO
         do iorb = 1, norbs
            kl = LocalOrbitals(js,atom,site)%orb2kl(iorb)
            m = mofk(kl)
            kl_bar = kl - 2*m
            do klp = 1, kmax_kkr
               if (LocalOrbitals(js,atom,site)%lm_flag(klp,iorb) > 0) then
                  mp = mofk(klp)
                  cfac = m1m(mp+m)
                  klp_bar = klp - 2*mp
                  do ir = 1, NumRs
                     f_star(ir,klp,iorb) = cfac*LocalOrbitals(js,atom,site)%orbital(ir,klp_bar,kl_bar)
                  enddo
               endif
            enddo
         enddo
      enddo
   else
      call ErrorHandler('computeOrbitalBasis','Invalid basis function type',BasisType)
      BasisFactor = CONE
   endif
!
   nullify(f_star)
!
   end subroutine computeOrbitalBasis
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOrbitalBasis(spin,site,atom,orb,rpow,star,bfac) result(f)
!  ===================================================================
!
!  Note: spin = 1 or n_spin_cant
!
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, site, atom, orb
   integer (kind=IntKind), intent(out) :: rpow
   integer (kind=IntKind) :: NumRs, NumOrbitals, kl
!
   logical, intent(in) :: star
!
   complex (kind=CmplxKind), intent(out) :: bfac
   complex (kind=CmplxKind), pointer :: f(:,:) ! returns the basis function multiplied by r^rpow
   complex (kind=CmplxKind), pointer :: f_star(:,:,:)
!
   if (LocalOrbitals(spin,atom,site)%NumOrbitals == 0) then
      call ErrorHandler('getOrbitalBasis','The orbital basis functions are not available')
   endif
!
   kl = LocalOrbitals(spin,atom,site)%orb2kl(orb)
   if (star) then
      NumOrbitals = LocalOrbitals(spin,atom,site)%NumOrbitals
      NumRs = LocalOrbitals(spin,atom,site)%NumRs
      f_star => aliasArray3_c(LocalOrbitals(spin,atom,site)%orbital_star,NumRs,kmax_kkr,NumOrbitals)
      f => f_star(:,:,kl)
   else
      f => LocalOrbitals(spin,atom,site)%orbital(:,:,kl)
   endif
   rpow = LocalOrbitals(spin,atom,site)%r_power
!
   bfac = BasisFactor
!
   end function getOrbitalBasis
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
