!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOrbitalBasis(orb,spin,site,atom,btype,star,rpow) result(f)
!  ===================================================================
   use KindParamModule, only : IntKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getRegSolution
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: orb,spin,site,atom
   integer (kind=IntKind), intent(in) :: btype
   integer (kind=IntKind), intent(out) :: rpow
!
   logical, intent(in), optional :: star
   logical :: apply_star = .false.
!
   complex (kind=CmplxKind), pointer :: f(:,:) ! returns the basis function multiplied by r^rpow
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:)
!
   if (btype == 0) then ! Using the local regular solution as the basis function.
      PhiLr => getRegSolution(spin=spin,site=site,atom=atom)
      if (present(star)) then
         apply_star = star
      else
         apply_star = .false.
      endif
      if (apply_star) then
      else
         f => PhiLr(:,:,orb)
      endif
   else
      call ErrorHandler('getOrbitalBasis','Invalid basis function type',btype)
   endif
!
   end function getOrbitalBasis
