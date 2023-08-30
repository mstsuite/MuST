module OrbitalBasisModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initOrbitalBasis,    &
          endOrbitalBasis,     &
          computeOrbitalBasis, &
          getOrbitalBasis
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
   type LocalOrbitalStruct
      integer (kind=IntKind) :: NumOrbitals
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind) :: r_power
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
   subroutine initOrbitalBasis(na,lmax,pola,cant,rel,btype)
!  ===================================================================
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na, lmax, pola, cant, rel, btype
   integer (kind=IntKind) :: id, ia, js
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
      do ia = 1, getLocalNumSpecies(id)
         do js = 1, n_spin_cant
            LocalOrbitals(js,ia,id)%NumOrbitals = 0
            LocalOrbitals(js,ia,id)%NumRs = 0
            if (BasisType == 0) then
               LocalOrbitals(js,ia,id)%r_power = 1
            else ! This needs to be re-examined..........
               LocalOrbitals(js,ia,id)%r_power = 0
            endif
            nullify(LocalOrbitals(js,ia,id)%orbital_star)
            nullify(LocalOrbitals(js,ia,id)%orbital)
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
   subroutine computeOrbitalBasis(spin,site,atom,norbs,e)
!  ===================================================================
!
!  Note: spin = 1 or n_spin_pola/n_spin_cant
!
!  *******************************************************************
   use MathParamModule, only : CZERO
!
   use IntegerFactorsModule, only : m1m, mofk
!
   use PotentialTypeModule, only : isFullPotential
!
   use RadialGridModule, only : getNumRmesh
!
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin,site,atom,norbs
   integer (kind=IntKind) :: kl, klp, ir, m, mp, kl_bar, klp_bar
   integer (kind=IntKind) :: js, ns, osize, NumRs
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
   do js = 1, n_spin_cant
      if (LocalOrbitals(js,atom,site)%NumOrbitals == 0) then
         LocalOrbitals(js,atom,site)%NumOrbitals = norbs
         if (BasisType == 0) then ! Using the local regular solution as the basis function.
            LocalOrbitals(js,atom,site)%NumRs = getSolutionRmeshSize(site=site)
         endif
         osize = getNumRmesh(site,'Max')*kmax_kkr*norbs
         allocate(LocalOrbitals(js,atom,site)%orbital_star(osize))
      else if (LocalOrbitals(js,atom,site)%NumOrbitals /= norbs) then
         call ErrorHandler('computeOrbitalBasis',                        &
                           'NumOrbitals <> norbs',LocalOrbitals(js,atom,site)%NumOrbitals,norbs)
      endif
   enddo
!
   if (BasisType == 0) then ! Using the local regular solution as the basis function.
      if (.not.present(e)) then
         call ErrorHandler('computeOrbitalBasis',                     &
                           'Energy parameter is required for this basis type.')
      endif
      if (isFullPotential()) then
         do js = 1, n_spin_cant
            ns = max(js,spin)
!           ----------------------------------------------------------
            call solveSingleScattering(spin=ns,site=site,e=adjustEnergy(spin,e), &
                                       vshift=CZERO)
!           ----------------------------------------------------------
         enddo
      else
         do js = 1, n_spin_cant
            ns = max(js,spin)
!           ----------------------------------------------------------
            call solveSingleScattering(spin=ns,site=site,e=adjustEnergy(spin,e), &
                                       vshift=CZERO,isSphSolver=.true.,useIrrSol='H')
!           ----------------------------------------------------------
         enddo
      endif
      do js = 1, n_spin_cant
         LocalOrbitals(js,atom,site)%orbital => getRegSolution(spin=js,site=site,atom=atom)
         NumRs = LocalOrbitals(js,atom,site)%NumRs
         f_star => aliasArray3_c(LocalOrbitals(js,atom,site)%orbital_star,NumRs,kmax_kkr,norbs)
         do kl = 1, norbs
            m = mofk(kl)
            kl_bar = kl - 2*m
            do klp = 1, kmax_kkr
               mp = mofk(klp)
               cfac = m1m(mp+m)
               klp_bar = klp - 2*mp
               do ir = 1, NumRs
                  f_star(ir,klp,kl) = cfac*LocalOrbitals(js,atom,site)%orbital(ir,klp_bar,kl_bar)
               enddo
            enddo
         enddo
      enddo
   else
      call ErrorHandler('computeOrbitalBasis','Invalid basis function type',BasisType)
   endif
!
   nullify(f_star)
!
   end subroutine computeOrbitalBasis
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOrbitalBasis(spin,site,atom,orb,rpow,star) result(f)
!  ===================================================================
!
!  Note: spin = 1 or n_spin_cant
!
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, site, atom, orb
   integer (kind=IntKind), intent(out) :: rpow
   integer (kind=IntKind) :: NumRs, NumOrbitals
!
   logical, intent(in), optional :: star
   logical :: apply_star = .false.
!
   complex (kind=CmplxKind), pointer :: f(:,:) ! returns the basis function multiplied by r^rpow
   complex (kind=CmplxKind), pointer :: f_star(:,:,:)
!
   if (LocalOrbitals(spin,atom,site)%NumOrbitals == 0) then
      call ErrorHandler('getOrbitalBasis','The orbital basis functions are not available')
   endif
!
   if (present(star)) then
      apply_star = star
   else
      apply_star = .false.
   endif
!
   if (apply_star) then
      NumOrbitals = LocalOrbitals(spin,atom,site)%NumOrbitals
      NumRs = LocalOrbitals(spin,atom,site)%NumRs
      f_star => aliasArray3_c(LocalOrbitals(spin,atom,site)%orbital_star,NumRs,kmax_kkr,NumOrbitals)
      f => f_star(:,:,orb)
   else
      f => LocalOrbitals(spin,atom,site)%orbital(:,:,orb)
   endif
   rpow = LocalOrbitals(spin,atom,site)%r_power
!
   end function getOrbitalBasis
!  ===================================================================
end module OrbitalBasisModule
