module SMatrixPolesModule
!  *******************************************************************
!  * Purpose: determine the poles of single scattering S-matrix      *
!  *          poles, and identigy the bound states and resonance     *
!  *          states from those poles.                               *
!  *******************************************************************
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE,    &
                               Ten2m8, FOURTH, CZERO, TEN2m7, TEN2m12, TEN2m10
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : MyPE, syncAllPEs
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEInGroup
   use GroupCommModule, only : GlobalSumInGroup, bcastMessageInGroup
   use GroupCommModule, only : syncAllPEsInGroup
!
public :: initSMatrixPoles,          &
          endSMatrixPoles,           &
          clearSMatrixPoles,         &
          getNumBoundStates,         &
          getNumBoundStateDegen,     &
          getBoundStateEnergy,       &
          getNumResonanceStates,     &
          getNumResonanceStateDegen, &
          getResonanceStateEnergy,   &
          isEnergyInResonanceRange,  &
          findSMatrixPoles,          &
          computeBoundStateDensity,  &
          getBoundStateDensity,      &
          getBoundStateChargeInCell, &
          computeResonanceStateDensity, & ! compute the density arising from the
                                          ! Lorenzian state
          getIntegratedResStateDensity, & ! Returns the Lorenzian state density
          getResidualResStateDensity,   & ! Returns the residual density at
                                          ! an energy in the resonance region
          printSMatrixPoleInfo,      &
          printBoundStateDensity,    &
          isSMatrixPolesInitialized
!
   interface computeBoundStateDensity
      module procedure computeBSD_residual, computeBSD_coutour
   end interface
!
private
   logical :: isInitialized = .false.
   logical :: rad_derivative = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola, n_spin_cant
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: eGID, NumPEsInEGroup, MyPEInEGroup
   integer (kind=IntKind) :: lmax_kkr_max, kmax_kkr_max, kmax_kkr_save
   integer (kind=IntKind) :: lmax_rho_max, kmax_rho_max, jmax_rho_max
   integer (kind=IntKind) :: MaxNumRs
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   integer (kind=IntKind), parameter :: MaxNumBoundStates = 10
   integer (kind=IntKind), parameter :: MaxNumResonanceStates = 5
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
   real (kind=RealKind), allocatable :: gaunt(:,:,:)
!
   real (kind=RealKind), parameter :: degen_tol = TEN2m6
!
   complex (kind=CmplxKind), allocatable, target :: wspace0(:), wspace1(:), wspace2(:), wspace3(:)
   complex (kind=CmplxKind), allocatable, target :: wspace4(:), wspace5(:), wspace6(:)
!
   type PoleDensityStruct
      integer (kind=IntKind) :: NumDegens
      real (kind=RealKind) :: PoleE
      real (kind=RealKind) :: PoleWidth
      real (kind=RealKind) :: Qvp
      real (kind=RealKind) :: Qmt
      complex (kind=CmplxKind), allocatable :: ResidualMat(:)
      complex (kind=CmplxKind), allocatable :: AuxiliaryMat(:)
      complex (kind=CmplxKind), allocatable :: Density(:)
      complex (kind=CmplxKind), allocatable :: Deriv_Density(:)
   end type PoleDensityStruct
!
   type PoleStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind) :: kmax_kkr
      integer (kind=IntKind) :: jmax_rho
      integer (kind=IntKind), allocatable :: NumBoundPoles(:,:)
      integer (kind=IntKind), allocatable :: NumResPoles(:,:)
      integer (kind=IntKind), allocatable :: BoundSI(:,:,:)
      integer (kind=IntKind), allocatable :: ResSI(:,:,:)
!
      real (kind=RealKind), allocatable :: ebot(:,:)
      real (kind=RealKind), allocatable :: etop(:,:)
!
      type (PoleDensityStruct), allocatable :: BoundState(:,:,:)
      type (PoleDensityStruct), allocatable :: ResState(:,:,:)
   end type PoleStruct
!
   type (PoleStruct), allocatable :: Pole(:)
!
   real (kind=RealKind) :: ResWidth
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSMatrixPoles(nla,npola,ncant,nspecies,lmax_kkr,lmax_rho,iprint)
!  ===================================================================
   use GauntFactorsModule, only : getK3, getNumK3, getGauntFactor
!
   use RadialGridModule, only : getMaxNumRmesh
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, &
                                     isQuadraticMatrixInitialized
!
   use IntegerFactorsModule, only : initIntegerFactors,               &
                                    isIntegerFactorsInitialized
!
   use ScfDataModule, only : Harris
!
   use ExchCorrFunctionalModule, only : isGGAFunctional
!
   use SingleScatteringDOSModule, only : setSScatteringDOSParam
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nla, npola, ncant, iprint
   integer (kind=IntKind), intent(in) :: nspecies(nla)
   integer (kind=IntKind), intent(in) :: lmax_kkr(nla)
   integer (kind=IntKind), intent(in) :: lmax_rho(nla)
   integer (kind=IntKind) :: id, kl1, kl2, kl3, lmax_max, i2, n
!
   logical, parameter :: isGeneral = .false.
!
   LocalNumAtoms = nla
   n_spin_pola = npola
   n_spin_cant = ncant
   print_level = iprint
!
   allocate(Pole(LocalNumAtoms))
!
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
!
!  NumPEsInEGroup = 1; MyPEinEGroup = 0
!
   if (iprint >= 0) then
      write(6,'(a,3i5)')'MyPE, MyPEinEGroup, NumPEsInEGroup = ',MyPE,MyPEinEGroup,NumPEsInEGroup
   endif
!
   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
!
   lmax_kkr_max = 0
   lmax_rho_max = 0
   do id = 1,LocalNumAtoms
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(id))
      lmax_rho_max = max(lmax_rho_max,lmax_rho(id))
      Pole(id)%NumSpecies = nspecies(id)
      Pole(id)%jmax_rho = (lmax_rho(id)+1)*(lmax_rho(id)+2)/2
      Pole(id)%kmax_kkr = (lmax_kkr(id)+1)**2
      Pole(id)%NumRs = 0
   enddo
!
   kmax_kkr_max = (lmax_kkr_max+1)**2
   kmax_kkr_save = kmax_kkr_max
   kmax_rho_max = (lmax_rho_max+1)**2
   jmax_rho_max = (lmax_rho_max+1)*(lmax_rho_max+2)/2
   lmax_max = max(lmax_kkr_max, lmax_rho_max)
!
   allocate( gaunt(kmax_kkr_max,kmax_kkr_max,kmax_rho_max) )
   gaunt = ZERO
   do kl3 = 1, kmax_rho_max
      do kl1 = 1, kmax_kkr_max
         do i2 = 1, nj3(kl1,kl3)
            kl2 = kj3(i2,kl1,kl3)
            if (kl2 <= kmax_kkr_max) then
               gaunt(kl2,kl1,kl3) = cgnt(i2,kl1,kl3)
            endif
         enddo
      enddo
   enddo
!
!  ===================================================================
!  calculate the charge density associated with each bound state
!  ===================================================================
   MaxNumRs = getMaxNumRmesh()
   allocate( wspace0(kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace1(kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace2(MaxNumRs*kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace3(MaxNumRs*kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace4(MaxNumRs*kmax_kkr_max*kmax_kkr_max) )
   allocate( wspace5(MaxNumRs*jmax_rho_max) )
   allocate( wspace6(MaxNumRs*jmax_rho_max) )
!
   do id = 1, LocalNumAtoms
      n = Pole(id)%NumSpecies
      allocate( Pole(id)%NumBoundPoles(n_spin_pola,n) ); Pole(id)%NumBoundPoles = 0
      allocate( Pole(id)%NumResPoles(n_spin_pola,n) );   Pole(id)%NumResPoles = 0
!
      allocate( Pole(id)%BoundState(MaxNumBoundStates,n_spin_pola,n) )
      allocate( Pole(id)%ResState(MaxNumResonanceStates,n_spin_pola,n) )
!
      allocate( Pole(id)%BoundSI(MaxNumBoundStates,n_spin_pola,n) )
      allocate( Pole(id)%ResSI(MaxNumResonanceStates,n_spin_pola,n) )
!
      allocate( Pole(id)%ebot(n_spin_pola,n), Pole(id)%etop(n_spin_pola,n) )
      Pole(id)%ebot = ZERO
      Pole(id)%etop = ZERO
   enddo
!
   if (.not.isQuadraticMatrixInitialized()) then
!     ----------------------------------------------------------------
      call initQuadraticMatrix(kmax_kkr_max,isGeneral)
!     ----------------------------------------------------------------
   endif
   if (.not.isIntegerFactorsInitialized()) then
!     ----------------------------------------------------------------
      call initIntegerFactors(lmax_max)
!     ----------------------------------------------------------------
   endif
!
   ResWidth = 0.001d0
!
   rad_derivative = isGGAFUnctional()
!  -------------------------------------------------------------------
   call setSScatteringDOSParam(n_spin_pola,n_spin_cant,rad_derivative,Harris,print_level)
!  -------------------------------------------------------------------
!
   isInitialized = .true.
!
   end subroutine initSMatrixPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSMatrixPoles()
!  ===================================================================
   use QuadraticMatrixModule, only : endQuadraticMatrix, &
                                     isQuadraticMatrixInitialized
!
   implicit none
!
   integer (kind=IntKind) :: id, ib, ia, is
!
   deallocate( gaunt )
   deallocate( wspace0 )
   deallocate( wspace1 )
   deallocate( wspace2 )
   deallocate( wspace3 )
   deallocate( wspace4 )
   deallocate( wspace5 )
   deallocate( wspace6 )
!
   do id = 1, LocalNumAtoms
      do ia = 1, Pole(id)%NumSpecies
         do is = 1, n_spin_pola
            do ib = 1, Pole(id)%NumBoundPoles(is,ia)
               call finalizePoleDensity(Pole(id)%BoundState(ib,is,ia))
            enddo
            do ib = 1, Pole(id)%NumResPoles(is,ia)
               call finalizePoleDensity(Pole(id)%ResState(ib,is,ia))
            enddo
         enddo
      enddo
      deallocate( Pole(id)%ebot, Pole(id)%etop )
      deallocate( Pole(id)%NumBoundPoles )
      deallocate( Pole(id)%NumResPoles )
      deallocate( Pole(id)%BoundState )
      deallocate( Pole(id)%ResState )
      deallocate( Pole(id)%BoundSI )
      deallocate( Pole(id)%ResSI )
   enddo
   deallocate(Pole)
!
   if (isQuadraticMatrixInitialized()) then
!     ----------------------------------------------------------------
      call endQuadraticMatrix()
!     ----------------------------------------------------------------
   endif
!
   isInitialized = .false.
!
   end subroutine endSMatrixPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine clearSMatrixPoles()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id, ib, ia, is
!
   do id = 1, LocalNumAtoms
      do ia = 1, Pole(id)%NumSpecies
         do is = 1, n_spin_pola
            do ib = 1, Pole(id)%NumBoundPoles(is,ia)
               call finalizePoleDensity(Pole(id)%BoundState(ib,is,ia))
            enddo
            do ib = 1, Pole(id)%NumResPoles(is,ia)
               call finalizePoleDensity(Pole(id)%ResState(ib,is,ia))
            enddo
         enddo
      enddo
      Pole(id)%NumBoundPoles = 0
      Pole(id)%NumResPoles = 0
!
      Pole(id)%BoundSI = 0
      Pole(id)%ResSI = 0
!
      Pole(id)%ebot = ZERO
      Pole(id)%etop = ZERO
   enddo
!
   end subroutine clearSMatrixPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSMatrixPolesInitialized() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   y = isInitialized
!
   end function isSMatrixPolesInitialized
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumBoundStates(id,ia,is) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
   integer (kind=IntKind) :: n
!
   n = Pole(id)%NumBoundPoles(is,ia)
!
   end function getNumBoundStates
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumBoundStateDegen(id,ia,is,ibs,sorted) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ibs, is, id, ia
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: n, ib, NumBPs
!
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
   if (NumBPs < 1) then
      n = 0
      return
   else if (ibs < 1 .or. ibs > NumBPs) then
      call ErrorHandler('getNumBoundStateDegen','Invalid bound state index',ibs)
   endif
!
   ib = ibs
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(id)%BoundSI(ibs,is,ia)
      endif
   endif
!
   n = Pole(id)%BoundState(ib,is,ia)%NumDegens
!
   end function getNumBoundStateDegen
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBoundStateEnergy(id,ia,is,ibs,sorted) result(e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ibs, ia
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: ib, NumBPs
!
   real (kind=RealKind) :: e
!
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
   if (NumBPs < 1) then
      e = ZERO
      return
   else if (ibs < 1 .or. ibs > NumBPs) then
      call ErrorHandler('getBoundStateEnergy','Invalid bound state index',ibs)
   endif
!
   ib = ibs
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(id)%BoundSI(ibs,is,ia)
      endif
   endif
!
   e = Pole(id)%BoundState(ib,is,ia)%PoleE
!
   end function getBoundStateEnergy
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBoundStateChargeInCell(id,ia,is,ibs,qmt,sorted) result(qvp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ibs, ia
   logical, optional, intent(in) :: sorted
   real (kind=RealKind), intent(out), optional :: qmt
!
   integer (kind=IntKind) :: ib, NumBPs
   real (kind=RealKind) :: qvp
!
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
   if (NumBPs < 1) then
      qvp = ZERO
      return
   else if (ibs < 1 .or. ibs > NumBPs) then
      call ErrorHandler('getBoundStateChargeInCell','Invalid bound state index',ibs)
   endif
!
   ib = ibs
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(id)%BoundSI(ibs,is,ia)
      endif
   endif
!
   qvp = Pole(id)%BoundState(ib,is,ia)%Qvp
!
   if (present(qmt)) then
      qmt = Pole(id)%BoundState(ib,is,ia)%Qmt
   endif
!
   end function getBoundStateChargeInCell
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBoundStateDensity(id,ia,is,ibs,NumRs,jmax_rho,derivative,sorted) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is, ibs
   integer (kind=IntKind), intent(out), optional :: NumRs, jmax_rho
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: NumBPs
   integer (kind=IntKind) :: ib
!
   complex (kind=CmplxKind), intent(out), optional, pointer :: derivative(:,:)
   complex (kind=CmplxKind), pointer :: p(:,:)
!
   if (present(NumRs)) then
       NumRs = Pole(id)%NumRs
   endif
!
   if (present(jmax_rho)) then
      jmax_rho = Pole(id)%jmax_rho
   endif
!
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
   if (NumBPs < 1) then
      nullify(p)
      return
   else if (ibs < 1 .or. ibs > NumBPs) then
      call ErrorHandler('getBoundStateDensity','Invalid bound state index',ibs)
   endif
!
   ib = ibs
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(id)%BoundSI(ibs,is,ia)
      endif
   endif
!
   if (present(derivative)) then
      if (rad_derivative) then
         derivative => aliasArray2_c(Pole(id)%BoundState(ib,is,ia)%Deriv_Density, &
                                     Pole(id)%NumRs,Pole(id)%jmax_rho)
      else
         nullify(derivative)
      endif
   endif
!
   p => aliasArray2_c(Pole(id)%BoundState(ib,is,ia)%Density,          &
                      Pole(id)%NumRs,Pole(id)%jmax_rho)
!
   end function getBoundStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumResonanceStates(site,atom,spin) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: spin, site, atom
   integer (kind=IntKind) :: n
!
   n = Pole(site)%NumResPoles(spin,atom)
!
   end function getNumResonanceStates
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumResonanceStateDegen(id,ia,is,ibs,sorted) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ibs, is, id, ia
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: n, ib, NumRPs
!
   NumRPs = Pole(id)%NumResPoles(is,ia)
   if (NumRPs < 1) then
      n = 0
      return
   else if (ibs < 1 .or. ibs > NumRPs) then
      call ErrorHandler('getNumResonanceStateDegen','Invalid resonance state index',ibs)
   endif
!
   ib = ibs
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(id)%ResSI(ibs,is,ia)
      endif
   endif
!
   n = Pole(id)%ResState(ib,is,ia)%NumDegens
!
   end function getNumResonanceStateDegen
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getResonanceStateEnergy(id,ia,is,ibs,hw,sorted) result(e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ibs, is, id, ia
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: ib, NumRPs
!
   real (kind=RealKind), optional, intent(out) :: hw ! Half width
   real (kind=RealKind) :: e
!
   NumRPs = Pole(id)%NumResPoles(is,ia)
   if (NumRPs < 1) then
      e = ZERO
      return
   else if (ibs < 1 .or. ibs > NumRPs) then
      call ErrorHandler('getResonanceStateEnergy','Invalid resonance state index',ibs)
   endif
!
   ib = ibs
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(id)%ResSI(ibs,is,ia)
      endif
   endif
!
   e = Pole(id)%ResState(ib,is,ia)%PoleE
   if (present(hw)) then
      hw = Pole(id)%ResState(ib,is,ia)%PoleWidth
   endif
!
   end function getResonanceStateEnergy
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isEnergyInResonanceRange(e,site,atom,spin,res_id,sorted) result(y)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, spin
   integer (kind=IntKind), intent(out), optional :: res_id
   logical, optional, intent(in) :: sorted
   integer (kind=IntKind) :: i, j
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind) :: el, er
!
   logical :: y
!
   y = .false.
   j = 0
   LOOP_i: do i = 1, Pole(site)%NumResPoles(spin,atom)
!     el = Pole(site)%ResState(i,spin,atom)%PoleE - HALF*Pole(site)%ResState(i,spin,atom)%PoleWidth
      el = Pole(site)%ResState(i,spin,atom)%PoleE - HALF*ResWidth
!     er = Pole(site)%ResState(i,spin,atom)%PoleE + HALF*Pole(site)%ResState(i,spin,atom)%PoleWidth
      er = Pole(site)%ResState(i,spin,atom)%PoleE + HALF*ResWidth
      if (e >= el .and. e <= er) then
         j = i
         y = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
!
   if (present(res_id)) then
      res_id = j
      if (present(sorted)) then
         if (sorted) then
            res_id = Pole(site)%ResSI(j,spin,atom)
         endif
      endif
   endif
!
   end function isEnergyInResonanceRange
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIntegratedResStateDensity(site,atom,spin,rstate,       &
                                         NumRs,jmax_rho,derivative,sorted) result(p)
!  ===================================================================
!
!  Note: This function returns the integrated Lorenzian density, which is
!        the density of the resonance compotent (usually for l > 1)
!        integrated over the resonance energy.
!
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, spin, rstate
   integer (kind=IntKind), intent(out), optional :: NumRs, jmax_rho
   integer (kind=IntKind) :: ib
   logical, optional, intent(in) :: sorted
!
   complex (kind=CmplxKind), intent(out), optional, pointer :: derivative(:,:)
   complex (kind=CmplxKind), pointer :: p(:,:)
!
   if (present(NumRs)) then
       NumRs = Pole(site)%NumRs
   endif
!
   if (present(jmax_rho)) then
      jmax_rho = Pole(site)%jmax_rho
   endif
!
   if (Pole(site)%NumResPoles(spin,atom) < 1) then
      nullify(p)
      return
   else if (rstate < 1 .or. rstate > Pole(site)%NumResPoles(spin,atom)) then
      call ErrorHandler('getIntegratedResStateDensity',               &
                        'Invalid resonance state index',rstate)
   endif
!
   ib = rstate
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(site)%ResSI(rstate,spin,atom)
      endif
   endif
!
   if (present(derivative)) then
      if (rad_derivative) then
         derivative => aliasArray2_c(Pole(site)%ResState(ib,spin,atom)%Deriv_Density,   &
                                     Pole(site)%NumRs,Pole(site)%jmax_rho)
      else
         nullify(derivative)
      endif
   endif
!
   p => aliasArray2_c(Pole(site)%ResState(ib,spin,atom)%Density,      &
                      Pole(site)%NumRs,Pole(site)%jmax_rho)
!
   end function getIntegratedResStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getResidualResStateDensity(site,atom,spin,rstate,e,       &
                                       NumRs,jmax_rho,derivative,sorted) result(p)
!  ===================================================================
!
!  Note: This function returns the residual density, which is the density
!        excluding the resonance component, for an energy e in the resonance
!        region.
!
!  *******************************************************************
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, spin, rstate
   integer (kind=IntKind), intent(out), optional :: NumRs, jmax_rho
   integer (kind=IntKind) :: ib
   logical, optional, intent(in) :: sorted
!
   complex (kind=CmplxKind), intent(out), optional, pointer :: derivative(:,:)
   complex (kind=CmplxKind), pointer :: p(:,:)
!
   real (kind=RealKind), intent(in) :: e
!
   if (present(NumRs)) then
       NumRs = Pole(site)%NumRs
   endif
!
   if (present(jmax_rho)) then
      jmax_rho = Pole(site)%jmax_rho
   endif
!
   if (Pole(site)%NumResPoles(spin,atom) < 1) then
      nullify(p)
      return
   else if (rstate < 1 .or. rstate > Pole(site)%NumResPoles(spin,atom)) then
      call ErrorHandler('getResidualResStateDensity',                 &
                        'Invalid resonance state index',rstate)
   endif
!
   ib = rstate
   if (present(sorted)) then
      if (sorted) then
         ib = Pole(site)%ResSI(rstate,spin,atom)
      endif
   endif
!
!  -------------------------------------------------------------------
   call computeResidualRSDensity(site,atom,spin,ib,e)
!  -------------------------------------------------------------------
!
   if (present(derivative)) then
      derivative => aliasArray2_c(wspace6,Pole(site)%NumRs,Pole(site)%jmax_rho)
   endif
!
   p => aliasArray2_c(wspace5,Pole(site)%NumRs,Pole(site)%jmax_rho)
!
   end function getResidualResStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSMatrixPoleInfo(id,ia,is)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
   integer (kind=IntKind) :: ib, ibs
!
   write(6,'(/,3(a,i3))')'Spin index: ',is,',  Site index: ',id,',  Species index: ',ia
   if (Pole(id)%ebot(is,ia) < -TEN2m8) then
      write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of bound states found within (',  &
                                       Pole(id)%ebot(is,ia),', ',ZERO,'): ',Pole(id)%NumBoundPoles(is,ia)
      do ibs = 1, Pole(id)%NumBoundPoles(is,ia)
         ib = Pole(id)%BoundSI(ibs,is,ia)
         write(6,'(a,i2,5x,a,f20.12)')'Degeneracy = ',Pole(id)%BoundState(ib,is,ia)%NumDegens, &
                                      ', Bound state energy = ',Pole(id)%BoundState(ib,is,ia)%PoleE
      enddo
   endif
   if (Pole(id)%etop(is,ia) > TEN2m8) then
      write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of resonance states found within (',      &
                                       Pole(id)%ebot(is,ia),', ',Pole(id)%etop(is,ia),'): ', &
                                       Pole(id)%NumResPoles(is,ia)
      write(6,'(a,f13.8)')'The width threshold for a resonance state = ',ResWidth
      do ibs = 1, Pole(id)%NumResPoles(is,ia)
         ib = Pole(id)%ResSI(ibs,is,ia)
         write(6,'(a,i2,5x,a,f20.12,a,f20.12)')'Degeneracy = ',Pole(id)%ResState(ib,is,ia)%NumDegens,           &
                                               ', Resonance state energy = ',Pole(id)%ResState(ib,is,ia)%PoleE, &
                                               ', Width = ',Pole(id)%ResState(ib,is,ia)%PoleWidth
      enddo
   endif
!
   end subroutine printSMatrixPoleInfo
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printBoundStateDensity(id,ia,is)
!  ===================================================================
   use MathParamModule, only : Y0, TEN2m6
!
   use RadialGridModule, only : getGrid
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use WriteFunctionModule, only : writeFunction
!
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
   integer (kind=IntKind) :: ib, ibs, ir, jl, jmax_rho, kmax_rho, NumRs, occ, NumBPs
!
   real (kind=RealKind) :: rfac, q_VP, q_mt, deg
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable :: rho0(:)
!
   type (GridStruct), pointer :: Grid
!
   logical :: non_zero
!
   character (len=6) :: state_string, jl_string
   character (len=8) :: app_string
   character (len=50) :: file_name
!
   complex (kind=CmplxKind), pointer :: Bdensity(:,:)
!
   Grid => getGrid(id)
   r_mesh => Grid%r_mesh
   NumRs = Pole(id)%NumRs
   jmax_rho = Pole(id)%jmax_rho
   kmax_rho = kofj(jmax_rho)
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
   if (NumBPs < 1) then
      return
   endif
   if (jmax_rho < 1) then
      call ErrorHandler('printBoundStateDensity','Invalid jmax_rho',jmax_rho)
   else if (NumRs < 1) then
      call ErrorHandler('printBoundStateDensity','Invalid NumRs',NumRs)
   endif
   allocate(rho0(NumRs))
!
   write(6,'(/,a)')'*********************************************'
   write(6,  '(a)')'*   Print out from printBoundStateDensity   *'
   write(6,'(a,/)')'*********************************************'
   write(6,'(3(a,i2))')'Local site index = ',id,', species index = ',ia,', spin index = ',is
   write(6,'(a,/,a,/,a)')                                             &
      '=============================================================',&
      '     energy       Degeneracy    Occupancy     q_VP      q_MT ',&
      '-------------------------------------------------------------'
   do ibs = 1, NumBPs
      ib = Pole(id)%BoundSI(ibs,is,ia)
!     ================================================================
!     Note: ResidualMat already contains a summation of the
!           contributions from the degenerate poles.
!     ================================================================
      occ = (3-n_spin_pola)*Pole(id)%BoundState(ib,is,ia)%NumDegens
      deg = real(Pole(id)%BoundState(ib,is,ia)%NumDegens,kind=RealKind)
      Bdensity => aliasArray2_c(Pole(id)%BoundState(ib,is,ia)%Density,NumRs,jmax_rho)
      q_VP=getVolumeIntegration(id,NumRs,r_mesh,kmax_rho,jmax_rho,0,Bdensity,q_MT)
      write(6,'(f15.8,4x,i5,8x,i5,4x,2f10.5)')                        &
         Pole(id)%BoundState(ib,is,ia)%PoleE,Pole(id)%BoundState(ib,is,ia)%NumDegens,occ,q_VP/deg,q_MT/deg
      if (id < 10 .and. ib < 10) then
         write(app_string,'(a,i1,a,i1,a,i1)')'a',id,'s',is,'e',ib
      else if (id >= 10 .and. ib < 10) then
         write(app_string,'(a,i2,a,i1,a,i1)')'a',id,'s',is,'e',ib
      else if (id < 10 .and. ib >= 10) then
         write(app_string,'(a,i1,a,i1,a,i2)')'a',id,'s',is,'e',ib
      else
         write(app_string,'(a,i2,a,i1,a,i2)')'a',id,'s',is,'e',ib
      endif
      if (occ < 10) then
         write(state_string,'(a,i1)')'Occ',occ
      else
         write(state_string,'(a,i2)')'Occ',occ
      endif
!
!     ================================================================
!     Notes: The spherical density data in the file includes r^2, but
!            does NOT include a factor for the number of degeneracies.
!     ================================================================
      rfac = Y0/deg
      do ir = 1, NumRs
         rho0(ir) = real(Bdensity(ir,1),kind=RealKind)*rfac
      enddo
      file_name = 'BndState_'//trim(state_string)//trim(app_string)//'Sph'
!     ----------------------------------------------------------------
      call writeFunction(file_name,NumRs,r_mesh,rho0,2)
!     ----------------------------------------------------------------
!
!     ================================================================
!     Notes: The full density data in the file does not include r^2
!     ================================================================
      do jl = 1, jmax_rho
         non_zero = .false.
         LOOP_ir: do ir = 1, NumRs
            if (abs(Bdensity(ir,jl)) > TEN2m6) then
               non_zero = .true.
               exit LOOP_ir
            endif
         enddo LOOP_ir
         if (non_zero) then
            if (lofj(jl) < 10) then
               write(jl_string,'(a,i1,a,i1)')'l',lofj(jl),'m',mofj(jl)
            else if (mofj(jl) < 10) then
               write(jl_string,'(a,i2,a,i1)')'l',lofj(jl),'m',mofj(jl)
            else
               write(jl_string,'(a,i2,a,i2)')'l',lofj(jl),'m',mofj(jl)
            endif
            file_name = 'BndState_'//trim(state_string)//trim(app_string)//jl_string
!           ----------------------------------------------------------
            call writeFunction(file_name,NumRs,r_mesh,Bdensity(:,jl))
!           ----------------------------------------------------------
         endif
      enddo
   enddo
   write(6,'(a)')                                                     &
      '============================================================='
   deallocate(rho0)
!
   end subroutine printBoundStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine findSMatrixPoles(id,ia,is,eb,et,Delta,MaxResWidth,CheckPoles,PanelOnZero)
!  ===================================================================
   use SSSolverModule, only : solveSingleScattering, getJostMatrix
   use SSSolverModule, only : getSolutionRmeshSize
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix
   use QuadraticMatrixModule, only : solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation, getEigenVector,     &
                                     getResidualMatrix, getAuxiliaryMatrix
!
   use WriteMatrixModule,  only : writeMatrix
!
   use ScfDataModule, only : CurrentScfIteration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
!
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind), intent(in), optional :: Delta
   real (kind=RealKind), intent(in), optional :: MaxResWidth
!
   logical, optional, intent(in) :: CheckPoles
   logical, optional, intent(in) :: PanelOnZero
!
   integer (kind=IntKind) :: ie, iw, kmax_kkr, NumWindows, info
   integer (kind=IntKind) :: i, j, l, m, kl, lp, mp, klp, nv, nb, nr, je, ip
   integer (kind=IntKind) :: MyNumWindows, nb0, nr0, ib, ir
   integer (kind=IntKind) :: bpdeg(kmax_kkr_max), rpdeg(kmax_kkr_max), degens(kmax_kkr_max)
   integer (kind=IntKind), allocatable :: nbr(:,:)
!
   logical :: isZeroInterval = .false.
   logical :: chkpole = .false.
   logical :: found
!
   real (kind=RealKind) :: WindowWidth
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0, err, w
!
   real (kind=RealKind) :: bpe(kmax_kkr_max), ebr(kmax_kkr_max), bpe_prev, rpe_prev
!
   complex (kind=CmplxKind) :: e, a2l, a2lp, c, cde, ce0, det
   complex (kind=CmplxKind), pointer :: jost_mat(:,:)
   complex (kind=CmplxKind), pointer :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), pointer :: sm(:,:)
   complex (kind=CmplxKind), pointer :: pv(:), evr(:), evl(:), em(:,:), am(:,:)
   complex (kind=CmplxKind) :: vt(kmax_kkr_max), diag(kmax_kkr_max)
   complex (kind=CmplxKind) :: rpe(kmax_kkr_max)
   complex (kind=CmplxKind) :: erc(kmax_kkr_max)
   complex (kind=CmplxKind) :: bmat(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max)
   complex (kind=CmplxKind) :: brmat(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max)
   complex (kind=CMplxKind) :: rmat(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max)
   complex (kind=CMplxKind) :: amat(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max)
   complex (kind=CMplxKind) :: camat(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max)
!
   if (eb > et) then
      call ErrorHandler('findSMatrixPoles','eb > et',eb,et)
   else
!     ----------------------------------------------------------------
      call clearSMatrixPoles()
!     ----------------------------------------------------------------
   endif
!
   if (present(Delta)) then
      WindowWidth = 4.0d0*Delta
   else
      WindowWidth = 0.01d0
   endif
!
   if (present(MaxResWidth)) then
      ResWidth = abs(MaxResWidth)
   endif
!
   kmax_kkr = Pole(id)%kmax_kkr
!
   if (.not.present(CheckPoles)) then
      chkpole = .false.
   else
      chkpole = CheckPoles
   endif
!
   if (.not.present(PanelOnZero)) then
      isZeroInterval = .false.
   else
      isZeroInterval = PanelOnZero
   endif
!
   Pole(id)%ebot(is,ia) = eb
   Pole(id)%etop(is,ia) = et
!
   if (isZeroInterval) then
      NumWindows = 1
      MyNumWindows = 1
!     de = FOURTH*(et-eb)
      de = HALF*(et-eb)
   else
      NumWindows = ceiling((et-eb)/WindowWidth)
!     NumWindows = NumWindows - mod(NumWindows,4)
      if (NumWindows < NumPEsInEGroup) then
         NumWindows = NumPEsInEGroup
         MyNumWindows = 1
      else
         NumWindows = ceiling(NumWindows/real(NumPEsInEGroup))*NumPEsInEGroup
         MyNumWindows = NumWindows/NumPEsInEGroup
      endif
!     WindowWidth = (et-eb)/real(NumWindows,kind=RealKind)
      if (present(Delta)) then
         de = Delta
      else
         de = WindowWidth/4.0d0
      endif
   endif
!
   de2 = de*TWO; dede2 = de*de*TWO
   cde = de
!  write(6,'(a,2d15.8,i5)')'de,win,nw = ',de,WindowWidth,NumWindows
!
   s0 => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   s1 => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
   s2 => aliasArray2_c(wspace2,kmax_kkr,kmax_kkr)
   sm => aliasArray2_c(wspace3,kmax_kkr,kmax_kkr)
!
!  ===================================================================
!  Hopefully, QuadraticMatrixModule could be updated in the future
!  so that kmax_kkr_max is made to be the leading matrix dimension and
!  the following 4 lines of the code become unnecessary
!  ===================================================================
   if (kmax_kkr /= kmax_kkr_save) then
      call endQuadraticMatrix()
      call initQuadraticMatrix(kmax_kkr)
      kmax_kkr_save = kmax_kkr
   endif
!
   nr = 0; nb = 0
   bpe = ZERO; rpe = CZERO
   bpe_prev = ZERO; rpe_prev = ZERO
   bpdeg = 0; rpdeg = 0
   bmat = CZERO; rmat = CZERO
   do iw = 1, MyNumWindows
      w0 = eb + (iw+MyPEInEGroup*MyNumWindows-1)*WindowWidth
      e0 = w0 + (HALF)*WindowWidth
!     write(6,'(a,i3,a,f6.3,a,f6.3,a)')'Window:',iw,'  (',w0,',',w0+WindowWidth,')'
      if (isZeroInterval) then
         e0 = ZERO
      else if ((abs(e0) < Ten2m6 .or. abs(e0-de) < Ten2m6 .or. abs(e0+de) < Ten2m6)) then
         if (e0 < ZERO) then
            e0 = e0 - HALF*de
         else
            e0 = e0 + HALF*de
         endif
      endif
      ce0 = e0
!
      if (isZeroInterval) then
         s0(1:kmax_kkr,1:kmax_kkr) = CZERO
      else
         e = cmplx(e0,ZERO,kind=CmplxKind)
!        -------------------------------------------------------------
         call solveSingleScattering(is, id, ce0, CZERO, atom=ia)
!        -------------------------------------------------------------
         jost_mat => getJostMatrix(spin=is,site=id,atom=ia)
!        -------------------------------------------------------------
         call zcopy(kmax_kkr*kmax_kkr,jost_mat,1,s0,1)
!        -------------------------------------------------------------
      endif
!
      e = cmplx(e0+de,ZERO,kind=CmplxKind)
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     call solveSingleScattering(is, id, ce0, -cde, atom=ia)
!     ----------------------------------------------------------------
      jost_mat => getJostMatrix(spin=is,site=id,atom=ia)
!     ----------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,jost_mat,1,s2,1)
!     ----------------------------------------------------------------
!
      e = cmplx(e0-de,ZERO,kind=CmplxKind)
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     call solveSingleScattering(is, id, ce0, cde, atom=ia)
!     ----------------------------------------------------------------
      jost_mat => getJostMatrix(spin=is,site=id,atom=ia)
!     ----------------------------------------------------------------
      call zcopy(kmax_kkr*kmax_kkr,jost_mat,1,s1,1)
!     ----------------------------------------------------------------
!
      s1 = (s2 - jost_mat)/de2
      s2 = (s2 + jost_mat - TWO*s0)/dede2
!
      if (isZeroInterval) then
!        -------------------------------------------------------------
         call solveLinearEquation(s1,s2,info)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call solveQuadraticEquation(s0,s1,s2,info)
!        -------------------------------------------------------------
      endif
!
      if (info /= 0) then
         stop 'Error in s0, s1, s2'
      endif
!
!     ----------------------------------------------------------------
      pv => getEigenValue(nv)
!     ----------------------------------------------------------------
      do ie = 1, nv
         if (.false.) then ! change to .false. to turn off self-checking
!           ==========================================================
!           Check eigenvalues and eigenvectors
!           ==========================================================
            write(6,'(a,i5,2d15.8)')'Index, Eigenvalue = ',ie,pv(ie)
            sm = s0 + s1*pv(ie) + s2*(pv(ie)*pv(ie))
            evr => getEigenVector('R',ie)
            vt = CZERO
            do j = 1, kmax_kkr
               do i = 1, kmax_kkr
                  vt(i) = vt(i) + sm(i,j)*evr(j)
               enddo
            enddo
            do i = 1, kmax_kkr
               err = abs(vt(i))
               if (err > ten2m7) then
                  call ErrorHandler('findSMatrixPoles','Right-side eigenvector error > 10^7',err)
               endif
            enddo
            write(6,'(a)')'Right-side eigenvector passed!'
!
            evl => getEigenVector('L',ie)
            vt = CZERO
            do j = 1, kmax_kkr
               do i = 1, kmax_kkr
                  vt(j) = vt(j) + evl(i)*sm(i,j)
               enddo
            enddo
            do i = 1, kmax_kkr
               err = abs(vt(i))
               if (err > ten2m7) then
                  call WarningHandler('findSMatrixPoles','Left-side eigenvector error > 10^7',err)
               endif
            enddo
            write(6,'(a)')'Left-side eigenvector passed!'
         endif
!        =============================================================
!        if (aimag(pv(ie)) > ZERO .and. real(pv(ie),kind=RealKind) + e0 > ZERO) then
!        if (abs(aimag(pv(ie))) < Ten2m8 .and. real(pv(ie),kind=RealKind)+e0 < ZERO) then   ! Bound states
         if (abs(aimag(pv(ie))) < Ten2m8 .and. real(pv(ie),kind=RealKind)+e0 < TEN2m6) then   ! Bound states
            pe = real(pv(ie),kind=RealKind) + e0
            if (pe >= w0 .and. pe <= w0+WindowWidth) then
!              -------------------------------------------------------
               em => getResidualMatrix(ie) ! em is the residule matrix of
                                           ! integrating sm^{-1} around its eigenvalue
!              -------------------------------------------------------
               if (size(em,1) /= kmax_kkr) then
                  call ErrorHandler('findSMatrixPoles','inconsistent matrix size',size(em,1),kmax_kkr)
               endif
               found = .false.
               do i = 1, nb
                  if (abs(pe-bpe(i)) < degen_tol) then
                     found = .true.
                     bpdeg(i) = bpdeg(i) + 1
!                    -------------------------------------------------
                     call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,bmat(1,i),1)
!                    -------------------------------------------------
                     exit
                  endif
               enddo
               if (.not.found) then
                  nb = nb + 1
                  bpe(nb) = pe
                  bpdeg(nb) = 1
!                 ----------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,em,1,bmat(1,nb),1)
!                 ----------------------------------------------------
               endif
    !          if (abs(pe-bpe_prev) > degen_tol) then
    !             nb = nb + 1
    !             bpe(nb) = pe
    !             bpdeg(nb) = 1
!                 ----------------------------------------------------
    !             call zcopy(kmax_kkr*kmax_kkr,em,1,bmat(1,nb),1)
!                 ----------------------------------------------------
    !             bpe_prev = pe
    !          else if (nb == 0) then
    !             call ErrorHandler('findSMatrixPoles','bound state pe = ZERO',pe)
    !          else ! In degeneracy case, em is added to bmat of the same energy
    !             bpdeg(nb) = bpdeg(nb) + 1
!                 ----------------------------------------------------
    !             call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,bmat(1,nb),1)
!                 ----------------------------------------------------
    !          endif
            endif
         else if (aimag(sqrt(pv(ie)+e0)) < ZERO) then  ! Resonance states
            pe = real(pv(ie),kind=RealKind) + e0
            w = aimag(sqrt(pv(ie)))**2
            if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > ZERO .and. TWO*w <= ResWidth) then
!              -------------------------------------------------------
               em => getResidualMatrix(ie) ! em is the residule matrix of
                                           ! integrating sm^{-1} around its eigenvalue
!              -------------------------------------------------------
               found = .false.
               do i = 1, nr
                  if (abs(pe-real(rpe(i),kind=RealKind)) < degen_tol) then
                     found = .true.
                     rpdeg(nr) = rpdeg(nr) + 1
!                    -------------------------------------------------
                     call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,rmat(1,i),1)
!                    -------------------------------------------------
                     exit
                  endif
               enddo
               if (.not.found) then
                  nr = nr + 1
                  rpe(nr) = cmplx(pe,aimag(sqrt(pv(ie)))**2,kind=CmplxKind)
                  rpdeg(nr) = 1
!                 ----------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,em,1,rmat(1,nr),1)
!                 ----------------------------------------------------
                  am => getAuxiliaryMatrix(ie,degen_tol)
!                 ----------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,am,1,amat(1,nr),1)
!                 ----------------------------------------------------
               endif
     !         if (abs(pe-rpe_prev) > degen_tol) then
     !            nr = nr + 1
     !            rpe(nr) = cmplx(pe,aimag(sqrt(pv(ie)))**2,kind=CmplxKind)
     !            rpdeg(nr) = 1
!                 ----------------------------------------------------
     !            call zcopy(kmax_kkr*kmax_kkr,em,1,rmat(1,nr),1)
!                 ----------------------------------------------------
     !            rpe_prev = pe
!
     !            am => getAuxiliaryMatrix(ie,degen_tol)
!                 ----------------------------------------------------
     !            call zcopy(kmax_kkr*kmax_kkr,am,1,amat(1,nr),1)
!                 ----------------------------------------------------
     !         else if (nr == 0) then
     !            call ErrorHandler('findSMatrixPoles','resonance state pe = ZERO',pe)
     !         else ! Add contribution from the degenerate state......
     !            rpdeg(nr) = rpdeg(nr) + 1
!    !            ----------------------------------------------------
     !            call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,rmat(1,nr),1)
!    !            ----------------------------------------------------
     !         endif
            endif
         endif
      enddo
   enddo
!
   if (chkpole) then
      do je = 1, nb
         e0 = bpe(je)
         do ie = -10, 10
            e = e0 + ie*0.001d0
            call solveSingleScattering(is, id, e, CZERO, atom=ia)
            jost_mat => getJostMatrix()
            call calcDet(jost_mat,kmax_kkr,det,diag)
            write(6,'(a,f12.5,2x,2d16.8)')'e,det = ',real(e),det
         enddo
         write(6,'(/)')
      enddo
      do je = 1, nr
         e0 = rpe(je)
         do ie = -10, 10
            e = e0 + ie*0.001d0
            call solveSingleScattering(is, id, e, CZERO, atom=ia)
            jost_mat => getJostMatrix()
            call calcDet(jost_mat,kmax_kkr,det,diag)
            write(6,'(a,f12.5,2x,2d16.8)')'e,det = ',real(e),det
         enddo
         write(6,'(/)')
      enddo
   endif
!
   allocate(nbr(2,NumPEsInEGroup))
   nbr = 0
   nbr(1,MyPEinEGroup+1) = nb
   nbr(2,MyPEinEGroup+1) = nr
   if (NumPEsInEGroup > 1) then
!     ----------------------------------------------------------------
      call GlobalSumInGroup(eGID,nbr,2,NumPEsInEGroup)
!     ----------------------------------------------------------------
   endif
   Pole(id)%NumBoundPoles(is,ia) = 0
   Pole(id)%NumResPoles(is,ia) = 0
   do ip = 1, NumPEsInEGroup
      Pole(id)%NumBoundPoles(is,ia) = Pole(id)%NumBoundPoles(is,ia) + nbr(1,ip)
      Pole(id)%NumResPoles(is,ia) = Pole(id)%NumResPoles(is,ia) + nbr(2,ip)
   enddo
   if (Pole(id)%NumBoundPoles(is,ia) > MaxNumBoundStates) then
      call ErrorHandler('findSMatrixPoles','NumBoundPoles > MaxNumBoundStates', &
                        Pole(id)%NumBoundPoles(is,ia), MaxNumBoundStates)
   endif
   if (Pole(id)%NumResPoles(is,ia) > MaxNumResonanceStates) then
      call ErrorHandler('findSMatrixPoles','NumResPoles > MaxNumResonanceStates', &
                        Pole(id)%NumResPoles(is,ia), MaxNumResonanceStates)
   endif
!
   do ib = 1, Pole(id)%NumBoundPoles(is,ia)
      call initializePoleDensity(Pole(id)%BoundState(ib,is,ia),aux=.false.)
   enddo
!
   do ib = 1, Pole(id)%NumResPoles(is,ia)
      call initializePoleDensity(Pole(id)%ResState(ib,is,ia),aux=.true.)
   enddo
!
   ebr = ZERO; erc = CZERO
   nb0 = 0; nr0 = 0
   do ip = 1, NumPEsInEGroup
      if (nbr(1,ip) > 0) then
         if (MyPEinEGroup == ip-1) then
            ebr(1:nbr(1,ip)) = bpe(1:nbr(1,ip))
!           ----------------------------------------------------------
            call zcopy(kmax_kkr_max*kmax_kkr_max*nbr(1,ip),bmat,1,brmat,1)
!           ----------------------------------------------------------
            degens(1:nbr(1,ip)) = bpdeg(1:nbr(1,ip))
         endif
         if (NumPEsInEGroup > 1) then
!           ----------------------------------------------------------
            call bcastMessageInGroup(eGID,ebr,nbr(1,ip),ip-1)
            call bcastMessageInGroup(eGID,degens,nbr(1,ip),ip-1)
            call bcastMessageInGroup(eGID,brmat,kmax_kkr_max*kmax_kkr_max,nbr(1,ip),ip-1)
!           ----------------------------------------------------------
         endif
         do ib = 1, nbr(1,ip)
            Pole(id)%BoundState(nb0+ib,is,ia)%PoleE = ebr(ib)
            Pole(id)%BoundState(nb0+ib,is,ia)%NumDegens = degens(ib)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,brmat(1,ib),1,               &
                       Pole(id)%BoundState(nb0+ib,is,ia)%ResidualMat,1)
!           ----------------------------------------------------------
         enddo
         nb0 = nb0 + nbr(1,ip)
      endif
!
      if (nbr(2,ip) > 0) then
         if (MyPEinEGroup == ip-1) then
            erc(1:nbr(2,ip)) = rpe(1:nbr(2,ip))
!           ----------------------------------------------------------
            call zcopy(kmax_kkr_max*kmax_kkr_max*nbr(2,ip),rmat,1,brmat,1)
!           ----------------------------------------------------------
            degens(1:nbr(2,ip)) = rpdeg(1:nbr(2,ip))
!           ----------------------------------------------------------
            call zcopy(kmax_kkr_max*kmax_kkr_max*nbr(2,ip),amat,1,camat,1)
!           ----------------------------------------------------------
         endif
         if (NumPEsInEGroup > 1) then
!           ----------------------------------------------------------
            call bcastMessageInGroup(eGID,erc,nbr(2,ip),ip-1)
            call bcastMessageInGroup(eGID,degens,nbr(2,ip),ip-1)
            call bcastMessageInGroup(eGID,brmat,kmax_kkr_max*kmax_kkr_max,nbr(2,ip),ip-1)
            call bcastMessageInGroup(eGID,camat,kmax_kkr_max*kmax_kkr_max,nbr(2,ip),ip-1)
!           ----------------------------------------------------------
         endif
         do ir = 1, nbr(2,ip)
            Pole(id)%ResState(nr0+ir,is,ia)%PoleE = real(erc(ir),kind=RealKind)
            Pole(id)%ResState(nr0+ir,is,ia)%PoleWidth = TWO*aimag(erc(ir))
            Pole(id)%ResState(nr0+ir,is,ia)%NumDegens = degens(ir)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,brmat(1,ir),1,               &
                       Pole(id)%ResState(nr0+ir,is,ia)%ResidualMat,1)
!           ----------------------------------------------------------
            call zcopy(kmax_kkr*kmax_kkr,camat(1,ir),1,               &
                       Pole(id)%ResState(nr0+ir,is,ia)%AuxiliaryMat,1)
!           ----------------------------------------------------------
         enddo
         nr0 = nr0 + nbr(2,ip)
      endif
   enddo
!
   Pole(id)%NumRs = getSolutionRmeshSize(id,isCSRadius=.true.)
!
!  -------------------------------------------------------------------
   call sortPoleStates(Pole(id)%NumBoundPoles(is,ia),                 &
                       Pole(id)%BoundState(:,is,ia),                  &
                       Pole(id)%BoundSI(:,is,ia)) 
   call sortPoleStates(Pole(id)%NumResPoles(is,ia),                   &
                       Pole(id)%ResState(:,is,ia),                    &
                       Pole(id)%ResSI(:,is,ia)) 
!  -------------------------------------------------------------------
!
   deallocate(nbr)
   nullify(pv, evl, evr, em, am)
   nullify(s0, s1, s2, sm)
!
   end subroutine findSMatrixPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initializePoleDensity(pd,aux)
!  ===================================================================
   implicit none
!
   type (PoleDensityStruct), intent(inout) :: pd
!
   logical, intent(in) :: aux
!
   pd%NumDegens = 0
   pd%PoleE = ZERO
   pd%PoleWidth = ZERO
   pd%Qvp = ZERO
   pd%Qmt = ZERO
!
   if ( .not.allocated(pd%ResidualMat) ) then
      allocate( pd%ResidualMat(kmax_kkr_max*kmax_kkr_max) )
   endif
   pd%ResidualMat = CZERO
!
   if (aux) then
      if ( .not.allocated(pd%AuxiliaryMat) ) then
         allocate( pd%AuxiliaryMat(kmax_kkr_max*kmax_kkr_max) )
      endif
      pd%AuxiliaryMat = CZERO
   endif 
!
   if ( .not.allocated(pd%Density) ) then
      allocate( pd%Density(MaxNumRs*jmax_rho_max) )
   endif
   pd%Density = CZERO
!
   if (rad_derivative) then
      if ( .not.allocated(pd%Deriv_Density) ) then
         allocate( pd%Deriv_Density(MaxNumRs*jmax_rho_max) )
      endif
      pd%Deriv_Density = CZERO
   endif
!
   end subroutine initializePoleDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine finalizePoleDensity(pd)
!  ===================================================================
   implicit none
!
   type (PoleDensityStruct), intent(inout) :: pd
!
   if ( allocated(pd%ResidualMat) ) then
      deallocate( pd%ResidualMat )
   endif
!
   if ( allocated(pd%AuxiliaryMat) ) then
      deallocate( pd%AuxiliaryMat )
   endif
!
   if ( allocated(pd%Density) ) then
      deallocate( pd%Density )
   endif
!
   if (rad_derivative) then
      if ( allocated(pd%Deriv_Density) ) then
         deallocate( pd%Deriv_Density )
      endif
   endif
!
   end subroutine finalizePoleDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeBSD_residual(id,ia,is)
!  ===================================================================
!  Note: The density associated with the bound state pole includes the 
!        degeneracy of the pole, since the residual matrix has already
!        been multiplied by the number of degeneracies.
!  *******************************************************************
   use SSSolverModule, only : solveSingleScattering, getSineMatrix
   use SSSolverModule, only : getRegSolution, getRegSolutionDerivative
!  use SSSolverModule, only : getJostInvMatrix, getOmegaHatMatrix
!  use SSSolverModule, only : computeDOS, getDOS
!
   use RadialGridModule, only : getGrid
!
   use MatrixModule, only : computeAStarTInv
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use ScfDataModule, only : CurrentScfIteration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ia
!
   integer (kind=IntKind) :: ie, ib, info, ip
   integer (kind=IntKind) :: kmax_kkr, jmax_rho, kmax_rho, NumRs, NumBPs
   integer (kind=IntKind) :: kl, klp, klp_bar, kl1, kl2, kl3, kl3_bar, m3, mp, ir, jl3
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: Deriv_Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: sine_mat(:,:), smat_inv(:,:), BSinv(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:), DerPhiLr(:,:,:)
   complex (kind=CmplxKind), pointer :: BPhiLr(:,:,:), DerBPhiLr(:,:,:), PPr(:,:,:)
!  complex (kind=CmplxKind), pointer :: jost_inv(:,:)
!  complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind) :: e, cfac, cfac0, cfac1, cfac2, kappa
!
   type (GridStruct), pointer :: Grid
!
   logical :: non_zero
!
   kmax_kkr = Pole(id)%kmax_kkr
   jmax_rho = Pole(id)%jmax_rho
   kmax_rho = kofj(jmax_rho)
   NumRs = Pole(id)%NumRs
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
!
   if (NumBPs < 1) then
      return
   endif
!
!  ===================================================================
!  Note: Bdensity stores the density multiplied by the number of degeneracies
!        of the bound state.
!  ===================================================================
   smat_inv => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   BSinv => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
!
   Grid => getGrid(id)
   r_mesh => Grid%r_mesh
   BPhiLr => aliasArray3_c(wspace2,NumRs,kmax_kkr,kmax_kkr)
   PPr => aliasArray3_c(wspace3,NumRs,kmax_kkr,kmax_kkr)
   DerBPhiLr => aliasArray3_c(wspace4,NumRs,kmax_kkr,kmax_kkr)
!
   do ib = 1, NumBPs
      Pole(id)%BoundState(ib,is,ia)%Density = CZERO
      if (rad_derivative) then
         Pole(id)%BoundState(ib,is,ia)%Deriv_Density = CZERO
      endif
   enddo
!
   do ib = MyPEinEGroup+1, NumBPs, NumPEsInEGroup
      Bdensity => aliasArray2_c(Pole(id)%BoundState(ib,is,ia)%Density,NumRs,jmax_rho)
      if (rad_derivative) then
         Deriv_Bdensity => aliasArray2_c(Pole(id)%BoundState(ib,is,ia)%Deriv_Density,NumRs,jmax_rho)
      else
         nullify(Deriv_Bdensity)
      endif
      e = Pole(id)%BoundState(ib,is,ia)%PoleE
      kappa = sqrt(e)
      cfac0 = HALF*kappa
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     ----------------------------------------------------------------
      sine_mat => getSineMatrix()
!     ================================================================
!     calculate sine_mat^(-T*) and store the result in smat_inv
!     ----------------------------------------------------------------
      call computeAStarTInv(sine_mat,kmax_kkr,kmax_kkr,smat_inv)
!     ----------------------------------------------------------------
!
!     ===================================================================
!     calculate ResidualMat*sine_mat^{-T*} and store the result in BSinv
!     ----------------------------------------------------------------
      call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
                 Pole(id)%BoundState(ib,is,ia)%ResidualMat,kmax_kkr,smat_inv,kmax_kkr, &
                 CZERO,BSinv,kmax_kkr)
!     ----------------------------------------------------------------
!
      PhiLr => getRegSolution()
!     ----------------------------------------------------------------
      call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
                 PhiLr,NumRs*kmax_kkr,BSinv,kmax_kkr,                    &
                 CZERO,BPhiLr,NumRs*kmax_kkr)
!     ----------------------------------------------------------------
      PPr = CZERO
      do klp = 1, kmax_kkr
         mp = mofk(klp)
         klp_bar = bofk(klp)
         cfac = m1m(mp)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               do ir = 1, NumRs
                  PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*BPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar)
               enddo
            enddo
         enddo
      enddo
      do jl3 = 1, jmax_rho
         m3 = mofj(jl3)
         kl3 = kofj(jl3)
         kl3_bar = bofk(kl3)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               cfac1 = cfac0*gaunt(kl1,kl2,kl3)
               cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
               do ir = 1, NumRs
                  Bdensity(ir,jl3) = Bdensity(ir,jl3)           &
                                      + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
!+ cfac1*PPr(ir,kl1,kl2) - conjg(cfac2*PPr(ir,kl1,kl2))
               enddo
            enddo
         enddo
      enddo
!
      if (rad_derivative) then
         DerPhiLr => getRegSolutionDerivative()
!        -------------------------------------------------------------
         call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                    DerPhiLr,NumRs*kmax_kkr,BSinv,kmax_kkr,              &
                    CZERO,DerBPhiLr,NumRs*kmax_kkr)
!        -------------------------------------------------------------
         PPr = CZERO
         do klp = 1, kmax_kkr
            mp = mofk(klp)
            klp_bar = bofk(klp)
            cfac = m1m(mp)
            do kl2 = 1, kmax_kkr
               do kl1 = 1, kmax_kkr
                  do ir = 1, NumRs
                     PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*DerBPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar) &
                                                       + cfac*BPhiLr(ir,kl1,klp)*DerPhiLr(ir,kl2,klp_bar)
                  enddo
               enddo
            enddo
         enddo
         do jl3 = 1, jmax_rho
            m3 = mofj(jl3)
            kl3 = kofj(jl3)
            kl3_bar = bofk(kl3)
            do kl2 = 1, kmax_kkr
               do kl1 = 1, kmax_kkr
                  cfac1 = cfac0*gaunt(kl1,kl2,kl3)
                  cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
                  do ir = 1, NumRs
                     Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)  &
                                               + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
!+ cfac1*PPr(ir,kl1,kl2) - conjg(cfac2*PPr(ir,kl1,kl2))
                  enddo
               enddo
            enddo
         enddo
      endif
!
!     ================================================================
!     Get rid of r^2 from Bdensity and Deriv_Bdensity
!     ================================================================
      do jl3 = 1, jmax_rho
         do ir = 1, NumRs
            Bdensity(ir,jl3) = Bdensity(ir,jl3)/r_mesh(ir)**2
         enddo
      enddo
      if (rad_derivative) then
         do jl3 = 1, jmax_rho
            do ir = 1, NumRs
               Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)/r_mesh(ir)**2
            enddo
         enddo
      endif
   enddo
   if (NumPEsInEGroup > 1) then
!     ---------------------------------------------------------------
      call syncAllPEsInGroup(eGID)
!     ---------------------------------------------------------------
   endif
!
   do ib = 1, NumBPs
      ip = mod(ib-1,NumPEsInEGroup)
      if (NumPEsInEGroup > 1) then
!        ------------------------------------------------------------
         call bcastMessageInGroup(eGID,Pole(id)%BoundState(ib,is,ia)%Density,NumRs*jmax_rho,ip)
!        ------------------------------------------------------------
         if (rad_derivative) then
!           ---------------------------------------------------------
            call bcastMessageInGroup(eGID,Pole(id)%BoundState(ib,is,ia)%Deriv_Density,NumRs*jmax_rho,ip)
!           ---------------------------------------------------------
         endif
      endif
      Bdensity => aliasArray2_c(Pole(id)%BoundState(ib,is,ia)%Density,NumRs,jmax_rho)
      Pole(id)%BoundState(ib,is,ia)%Qvp = getVolumeIntegration(id,NumRs,r_mesh,kmax_rho,   &
                                                               jmax_rho,0,Bdensity,        &
                                                               Pole(id)%BoundState(ib,is,ia)%Qmt)
      if (print_level >= 0) then
         write(6,'(a,f18.13)')'In computeBoundStateDensity: Qmt per spin = ',Pole(id)%BoundState(ib,is,ia)%Qmt
         write(6,'(a,f18.13)')'In computeBoundStateDensity: Qvp per spin = ',Pole(id)%BoundState(ib,is,ia)%Qvp
      endif
   enddo
!
   nullify(sine_mat, BSinv, smat_inv, Grid, r_mesh)
   nullify(BPhiLr, PhiLr, DerBPhiLr, DerPhiLr, PPr)
   nullify(Bdensity, Deriv_Bdensity)
!
   end subroutine computeBSD_residual
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeBSD_coutour(id,ia,is,bse_rad,bse_index,chempot,wk_dos)
!  ===================================================================
!  Note: This routine uses contour integration method to compute the
!        density due to the S-matrix poles, i.e. the bound states. 
!        The computed density is a sum of the contribution from all the poles
!        within the contour radius e_rad. Therefore, the density includes
!        the degeneracy of the poles.
!  *******************************************************************
   use SSSolverModule, only : solveSingleScattering, getSineMatrix
   use SSSolverModule, only : getRegSolution, getRegSolutionDerivative
!  use SSSolverModule, only : getJostInvMatrix, getOmegaHatMatrix
!  use SSSolverModule, only : computeDOS, getDOS
!
   use RadialGridModule, only : getGrid
!
   use MatrixModule, only : computeAStarTInv
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use ScfDataModule, only : CurrentScfIteration, Temperature
!
   use InputModule, only : getKeyValue
!
   use SingleScatteringDOSModule, only : setSScatteringDOSParam
   use SingleScatteringDOSModule, only : getSScatteringDOSofCmplxE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ia
   integer (kind=IntKind), intent(in) :: bse_index
!
   integer (kind=IntKind) :: ie, info(6), rstatus, nsize, jl, ir
   integer (kind=IntKind) :: jmax_rho, kmax_rho, NumRs, NumBPs
   integer (kind=IntKind) :: NumGQPs
   integer (kind=IntKind), parameter :: MaxGQPs = 30
!
   real (kind=RealKind), intent(in) :: bse_rad, chempot
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind) :: ssDOS, eb
!
   complex (kind=CmplxKind), intent(inout) :: wk_dos(:)
   complex (kind=CmplxKind), pointer :: Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: Deriv_Bdensity(:,:)
   complex (kind=CmplxKind) :: eg(MaxGQPs), ew(MaxGQPs)
   complex (kind=CmplxKind) :: ec
   complex (kind=CmplxKind), allocatable :: wk_tmp(:)
!
   type (GridStruct), pointer :: Grid
!
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
   if (NumBPs < 1) then
      return
   else if (bse_index < 1) then
      call ErrorHandler('computeBoundStateDensity','bse_index < 1',bse_index)
   else if (bse_index > NumBPs) then
      call ErrorHandler('computeBoundStateDensity','bse_index > NumBPs',bse_index,NumBPs)
   endif
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('computeBoundStateDensity','local site index is out of bound',id)
   endif
!
   if (ia < 1 .or. ia > Pole(id)%NumSpecies) then
      call ErrorHandler('computeBoundStateDensity','local atom species index is out of bound',ia)
   endif
!
   jmax_rho = Pole(id)%jmax_rho
   kmax_rho = kofj(jmax_rho)
   NumRs = Pole(id)%NumRs
!
   Pole(id)%BoundState(bse_index,is,ia)%Density = CZERO
   if (rad_derivative) then
      Pole(id)%BoundState(bse_index,is,ia)%Deriv_Density = CZERO
   endif
!
!  ===================================================================
!  Note: Bdensity stores the density multiplied by the number of degeneracies
!        of the bound state.
!  ===================================================================
   Grid => getGrid(id)
   r_mesh => Grid%r_mesh

!  ===================================================================
!  Setup Gaussian quadrature on a semi-circle contour with radius = bse_rad
!  ===================================================================
   rstatus = getKeyValue(1,'No. Gauss Pts. along Bound State Contour',NumGQPs)
   if (bse_rad > 0.01d0 .and. NumGQPs <= 5) then
      NumGQPs = 5*int(bse_rad/0.01d0+HALF)
   endif
   if (NumGQPs > MaxGQPs) then
      call ErrorHandler('computeBSD_coutour','NumGQPs > MaxGQPs',NumGQPs,MaxGQPs)
   endif
!  -------------------------------------------------------------------
   call setupSemiCircleContour(NumGQPs,bse_rad,eg,ew)
!  -------------------------------------------------------------------
   eb = Pole(id)%BoundState(bse_index,is,ia)%PoleE
!
!  ===================================================================
!  Perform contour integration with radius = HALF*(e2-e1)
!  ===================================================================
   if ( print_level >= 0) then
      write(6,'(a,f8.5,a)')'Performing contour integration, with radius = ',bse_rad,   &
                             ', around the bound state over the'
      write(6,'(a,f8.5,a,f8.5,a,i5,a)')'energy domain: (',eb-bse_rad,',',eb+bse_rad,'), with ', &
                                       NumGQPs,' Gaussian quadrature points.'
      write(6,'(a)') &
         '=========================================================================================='
   endif
!
   info = 0
   nsize = size(wk_dos)
   allocate(wk_tmp(nsize))
!
!  -------------------------------------------------------------------
   call setSScatteringDOSParam(id,NumRs,jmax_rho)
   call setSScatteringDOSParam(chempot,Temperature)
!  -------------------------------------------------------------------
!
   info(1) = is; info(2) = id; info(3) = ia; info(4) = -1; info(5) = lofk(Pole(id)%kmax_kkr)
   ssDOS = ZERO; wk_dos = CZERO; wk_tmp = CZERO
   do ie = MyPEinEGroup+1, NumGQPs, NumPEsInEGroup
      ec = eb + eg(ie)
!     ================================================================
!     ew is the Gaussian quadrature weight and included in the DOS result
!     ================================================================
      ssDOS = ssDOS + getSScatteringDOSofCmplxE(info,ec,wk_tmp,ew(ie))
      wk_dos = wk_dos + wk_tmp
   enddo
   if (mod(NumGQPs,MyPEinEGroup+1) > 0) then
!     ================================================================
!     This is needed to make sure that all processes
!     in the group are properly synchronized...
!     ================================================================
      ec = eb + eg(NumGQPs)
      ssDOS = ssDOS + ZERO*getSScatteringDOSofCmplxE(info,ec,wk_tmp,ew(NumGQPs))
   endif
!  ===================================================================
!  Sum over the processors to get the integrated value
!  -------------------------------------------------------------------
   call GlobalSumInGroup(eGID,ssDOS)
   call GlobalSumInGroup(eGID,wk_dos,nsize)
!  -------------------------------------------------------------------
   deallocate(wk_tmp)
!
!  ===================================================================
!  Decode wk_dos and copy the results into 
!     Pole(id)%BoundState(bse_index,is,ia)%Density
!     Pole(id)%BoundState(bse_index,is,ia)%Deriv_Density
!  -------------------------------------------------------------------
   call decodeDOSDATA(NumRs,jmax_rho,wk_dos,Pole(id)%BoundState(bse_index,is,ia))
!  -------------------------------------------------------------------

   if (NumPEsInEGroup > 1) then
!     ---------------------------------------------------------------
      call syncAllPEsInGroup(eGID)
!     ---------------------------------------------------------------
   endif
!
   Bdensity => aliasArray2_c(Pole(id)%BoundState(bse_index,is,ia)%Density,NumRs,jmax_rho)
!  ===================================================================
!  Get rid of r^2 from Bdensity and Deriv_Bdensity
!  ===================================================================
   do jl = 1, jmax_rho
      do ir = 1, NumRs
         Bdensity(ir,jl) = Bdensity(ir,jl)/r_mesh(ir)**2
      enddo
   enddo
   if (rad_derivative) then
      Deriv_Bdensity => aliasArray2_c(Pole(id)%BoundState(bse_index,is,ia)%Deriv_Density,NumRs,jmax_rho)
      do jl = 1, jmax_rho
         do ir = 1, NumRs
            Deriv_Bdensity(ir,jl) = Deriv_Bdensity(ir,jl)/r_mesh(ir)**2
         enddo
      enddo
   endif
   Pole(id)%BoundState(bse_index,is,ia)%Qvp = getVolumeIntegration(id,NumRs,r_mesh,kmax_rho,   &
                                                                   jmax_rho,0,Bdensity,        &
                                                                   Pole(id)%BoundState(bse_index,is,ia)%Qmt)
   if (print_level >= 0) then
      write(6,'(a,f18.13)')'In computeBoundStateDensity: Qmt per spin = ',Pole(id)%BoundState(bse_index,is,ia)%Qmt
      write(6,'(a,f18.13)')'In computeBoundStateDensity: Qvp per spin = ',Pole(id)%BoundState(bse_index,is,ia)%Qvp
   endif
!
   nullify(Grid, r_mesh)
   nullify(Bdensity, Deriv_Bdensity)
!
   end subroutine computeBSD_coutour
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeResonanceStateDensity(id,ia,is)
!  ===================================================================
!  Note: The density associated with the resonance state pole includes the 
!        degeneracy of the pole, since the residual matrix has already
!        been multiplied by the number of degeneracies.
!  *******************************************************************
   use SSSolverModule, only : solveSingleScattering, getSineMatrix
   use SSSolverModule, only : getRegSolution, getRegSolutionDerivative
!  use SSSolverModule, only : getJostInvMatrix, getOmegaHatMatrix
!  use SSSolverModule, only : computeDOS, getDOS
!
   use RadialGridModule, only : getGrid
!
   use MatrixModule, only : computeAStarTInv
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ia
!
   integer (kind=IntKind) :: ie, ib, info, ip
   integer (kind=IntKind) :: kmax_kkr, jmax_rho, kmax_rho, NumRs, NumResPs
   integer (kind=IntKind) :: kl, klp, klp_bar, kl1, kl2, kl3, kl3_bar, m3, mp, ir, jl3
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: Deriv_Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: sine_mat(:,:), smat_inv(:,:), BSinv(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:), DerPhiLr(:,:,:)
   complex (kind=CmplxKind), pointer :: BPhiLr(:,:,:), DerBPhiLr(:,:,:), PPr(:,:,:)
!  complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind) :: e, cfac, cfac0, cfac1, cfac2, kappa
!
   type (GridStruct), pointer :: Grid
!
   logical :: non_zero
!
   kmax_kkr = Pole(id)%kmax_kkr
   jmax_rho = Pole(id)%jmax_rho
   kmax_rho = kofj(jmax_rho)
   NumRs = Pole(id)%NumRs
   NumResPs = Pole(id)%NumResPoles(is,ia)
!
   if (NumResPs < 1) then
      return
   endif
!
!  ===================================================================
!  Note: Bdensity stores the density multiplied by the number of degeneracies
!        of the resonance state.
!  ===================================================================
   smat_inv => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   BSinv => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
!
   Grid => getGrid(id)
   r_mesh => Grid%r_mesh
   BPhiLr => aliasArray3_c(wspace2,NumRs,kmax_kkr,kmax_kkr)
   PPr => aliasArray3_c(wspace3,NumRs,kmax_kkr,kmax_kkr)
   DerBPhiLr => aliasArray3_c(wspace4,NumRs,kmax_kkr,kmax_kkr)
!
   do ib = 1, NumResPs
      Pole(id)%ResState(ib,is,ia)%Density = CZERO
   enddo
   if (rad_derivative) then
      do ib = 1, NumResPs
         Pole(id)%ResState(ib,is,ia)%Deriv_Density = CZERO
      enddo
   endif
!
   do ib = MyPEinEGroup+1, NumResPs, NumPEsInEGroup
      Bdensity => aliasArray2_c(Pole(id)%ResState(ib,is,ia)%Density,NumRs,jmax_rho)
      if (rad_derivative) then
         Deriv_Bdensity => aliasArray2_c(Pole(id)%ResState(ib,is,ia)%Deriv_Density,NumRs,jmax_rho)
      else
         nullify(Deriv_Bdensity)
      endif
      e = Pole(id)%ResState(ib,is,ia)%PoleE
      kappa = sqrt(e)
      cfac0 = HALF*kappa
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     ----------------------------------------------------------------
      sine_mat => getSineMatrix()
!     ================================================================
!     calculate sine_mat^(-T*) and store the result in smat_inv
!     ----------------------------------------------------------------
      call computeAStarTInv(sine_mat,kmax_kkr,kmax_kkr,smat_inv)
!     ----------------------------------------------------------------
!
!     ===================================================================
!     calculate ResidualMat*sine_mat^{-T*} and store the result in BSinv
!     ----------------------------------------------------------------
      call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
                 Pole(id)%ResState(ib,is,ia)%ResidualMat,kmax_kkr,smat_inv,kmax_kkr, &
                 CZERO,BSinv,kmax_kkr)
!     ----------------------------------------------------------------
!
      PhiLr => getRegSolution()
!     ----------------------------------------------------------------
      call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
                 PhiLr,NumRs*kmax_kkr,BSinv,kmax_kkr,                    &
                 CZERO,BPhiLr,NumRs*kmax_kkr)
!     ----------------------------------------------------------------
      PPr = CZERO
      do klp = 1, kmax_kkr
         mp = mofk(klp)
         klp_bar = bofk(klp)
         cfac = m1m(mp)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               do ir = 1, NumRs
                  PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*BPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar)
               enddo
            enddo
         enddo
      enddo
      do jl3 = 1, jmax_rho
         m3 = mofj(jl3)
         kl3 = kofj(jl3)
         kl3_bar = bofk(kl3)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               cfac1 = cfac0*gaunt(kl1,kl2,kl3)
               cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
               do ir = 1, NumRs
                  Bdensity(ir,jl3) = Bdensity(ir,jl3)           &
                                      + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
               enddo
            enddo
         enddo
      enddo
!
      if (rad_derivative) then
         DerPhiLr => getRegSolutionDerivative()
!        -------------------------------------------------------------
         call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                    DerPhiLr,NumRs*kmax_kkr,BSinv,kmax_kkr,              &
                    CZERO,DerBPhiLr,NumRs*kmax_kkr)
!        -------------------------------------------------------------
         PPr = CZERO
         do klp = 1, kmax_kkr
            mp = mofk(klp)
            klp_bar = bofk(klp)
            cfac = m1m(mp)
            do kl2 = 1, kmax_kkr
               do kl1 = 1, kmax_kkr
                  do ir = 1, NumRs
                     PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*DerBPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar) &
                                                       + cfac*BPhiLr(ir,kl1,klp)*DerPhiLr(ir,kl2,klp_bar)
                  enddo
               enddo
            enddo
         enddo
         do jl3 = 1, jmax_rho
            m3 = mofj(jl3)
            kl3 = kofj(jl3)
            kl3_bar = bofk(kl3)
            do kl2 = 1, kmax_kkr
               do kl1 = 1, kmax_kkr
                  cfac1 = cfac0*gaunt(kl1,kl2,kl3)
                  cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
                  do ir = 1, NumRs
                     Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)  &
                                               + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
                  enddo
               enddo
            enddo
         enddo
      endif
!
!     ================================================================
!     Get rid of r^2 from Bdensity and Deriv_Bdensity
!     ================================================================
   !  do jl3 = 1, jmax_rho
   !     do ir = 1, NumRs
   !        Bdensity(ir,jl3) = Bdensity(ir,jl3)/r_mesh(ir)**2
   !     enddo
   !  enddo
   !  if (rad_derivative) then
   !     do jl3 = 1, jmax_rho
   !        do ir = 1, NumRs
   !           Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)/r_mesh(ir)**2
   !        enddo
   !     enddo
   !  endif
   enddo
   if (NumPEsInEGroup > 1) then
!     ---------------------------------------------------------------
      call syncAllPEsInGroup(eGID)
!     ---------------------------------------------------------------
   endif
!
   do ib = 1, NumResPs
      ip = mod(ib-1,NumPEsInEGroup)
      if (NumPEsInEGroup > 1) then
!        ------------------------------------------------------------
         call bcastMessageInGroup(eGID,Pole(id)%ResState(ib,is,ia)%Density,NumRs*jmax_rho,ip)
!        ------------------------------------------------------------
         if (rad_derivative) then
!           ---------------------------------------------------------
            call bcastMessageInGroup(eGID,Pole(id)%ResState(ib,is,ia)%Deriv_Density,NumRs*jmax_rho,ip)
!           ---------------------------------------------------------
         endif
      endif
      Bdensity => aliasArray2_c(Pole(id)%ResState(ib,is,ia)%Density,NumRs,jmax_rho)
      Pole(id)%ResState(ib,is,ia)%Qvp = getVolumeIntegration(id,NumRs,r_mesh,kmax_rho,   &
                                                             jmax_rho,2,Bdensity,        &
                                                             Pole(id)%ResState(ib,is,ia)%Qmt)
      if (print_level >= 0) then
         write(6,'(a,f12.6)')'In computeResonanceStateDensity: Qvp = ',Pole(id)%ResState(ib,is,ia)%Qvp
      endif
   enddo
!
   nullify(sine_mat, BSinv, smat_inv, Grid, r_mesh)
   nullify(BPhiLr, PhiLr, DerBPhiLr, DerPhiLr, PPr)
   nullify(Bdensity, Deriv_Bdensity)
!
   end subroutine computeResonanceStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeResidualRSDensity(id,ia,is,ib,e)
!  ===================================================================
   use MathParamModule, only : SQRTm1, PI
!
   use SSSolverModule, only : solveSingleScattering, getSineMatrix
   use SSSolverModule, only : getRegSolution, getRegSolutionDerivative
!  use SSSolverModule, only : getJostInvMatrix, getOmegaHatMatrix
!  use SSSolverModule, only : computeDOS, getDOS
!
   use RadialGridModule, only : getGrid
!
   use MatrixModule, only : computeAStarTInv
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ia, ib
!
   integer (kind=IntKind) :: ie, info, ip
   integer (kind=IntKind) :: kmax_kkr, jmax_rho, kmax_rho, NumRs, NumResPs
   integer (kind=IntKind) :: kl, klp, klp_bar, kl1, kl2, kl3, kl3_bar, m3, mp, ir, jl3
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   real (kind=RealKind), intent(in) :: e
!
   complex (kind=CmplxKind), pointer :: Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: Deriv_Bdensity(:,:)
   complex (kind=CmplxKind), pointer :: sine_mat(:,:), smat_inv(:,:), BSinv(:,:)
   complex (kind=CmplxKind), pointer :: PhiLr(:,:,:), DerPhiLr(:,:,:)
   complex (kind=CmplxKind), pointer :: BPhiLr(:,:,:), DerBPhiLr(:,:,:), PPr(:,:,:)
!  complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind) :: cfac, cfac0, cfac1, cfac2, kappa, ec
!
   type (GridStruct), pointer :: Grid
!
   logical :: non_zero
!
   if (ib < 1 .or. ib > Pole(id)%NumResPoles(is,ia)) then
      call ErrorHandler('computeResidualRSDensity','ib is out of range',ib, &
                        Pole(id)%NumResPoles(is,ia))
   endif
!
   kmax_kkr = Pole(id)%kmax_kkr
   jmax_rho = Pole(id)%jmax_rho
   kmax_rho = kofj(jmax_rho)
   NumRs = Pole(id)%NumRs
!
   wspace5 = CZERO; wspace6 = CZERO
   Bdensity => aliasArray2_c(wspace5,NumRs,jmax_rho)
   if (rad_derivative) then
      Deriv_Bdensity => aliasArray2_c(wspace6,NumRs,jmax_rho)
   else
      nullify(Deriv_Bdensity)
   endif
!
!  ===================================================================
!  Note: Bdensity stores the density multiplied by the number of degeneracies
!        of the bound state.
!  ===================================================================
   smat_inv => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   BSinv => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
!
   Grid => getGrid(id)
   r_mesh => Grid%r_mesh
   BPhiLr => aliasArray3_c(wspace2,NumRs,kmax_kkr,kmax_kkr)
   PPr => aliasArray3_c(wspace3,NumRs,kmax_kkr,kmax_kkr)
   DerBPhiLr => aliasArray3_c(wspace4,NumRs,kmax_kkr,kmax_kkr)
!
   ec = e
   kappa = sqrt(ec)
   cfac0 = HALF*kappa*SQRTm1/PI
!  -------------------------------------------------------------------
   call solveSingleScattering(is, id, ec, CZERO, atom=ia)
!  -------------------------------------------------------------------
   sine_mat => getSineMatrix()
!  ===================================================================
!  calculate sine_mat^(-T*) and store the result in smat_inv
!  -------------------------------------------------------------------
   call computeAStarTInv(sine_mat,kmax_kkr,kmax_kkr,smat_inv)
!  -------------------------------------------------------------------
!
!  ======================================================================
!  calculate AuxiliaryMat*sine_mat^{-T*} and store the result in BSinv
!  -------------------------------------------------------------------
   call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
              Pole(id)%ResState(ib,is,ia)%AuxiliaryMat,kmax_kkr,smat_inv,kmax_kkr, &
              CZERO,BSinv,kmax_kkr)
!  -------------------------------------------------------------------
!
   PhiLr => getRegSolution()
!  -------------------------------------------------------------------
   call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,          &
              PhiLr,NumRs*kmax_kkr,BSinv,kmax_kkr,                    &
              CZERO,BPhiLr,NumRs*kmax_kkr)
!  -------------------------------------------------------------------
   PPr = CZERO
   do klp = 1, kmax_kkr
      mp = mofk(klp)
      klp_bar = bofk(klp)
      cfac = m1m(mp)
      do kl2 = 1, kmax_kkr
         do kl1 = 1, kmax_kkr
            do ir = 1, NumRs
               PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*BPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar)
            enddo
         enddo
      enddo
   enddo
   do jl3 = 1, jmax_rho
      m3 = mofj(jl3)
      kl3 = kofj(jl3)
      kl3_bar = bofk(kl3)
      do kl2 = 1, kmax_kkr
         do kl1 = 1, kmax_kkr
            cfac1 = cfac0*gaunt(kl1,kl2,kl3)
            cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
            do ir = 1, NumRs
               Bdensity(ir,jl3) = Bdensity(ir,jl3)           &
                                   + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
            enddo
         enddo
      enddo
   enddo
! 
   if (rad_derivative) then
      DerPhiLr => getRegSolutionDerivative()
!     ----------------------------------------------------------------
      call zgemm('n','n',NumRs*kmax_kkr,kmax_kkr,kmax_kkr,CONE,       &
                 DerPhiLr,NumRs*kmax_kkr,BSinv,kmax_kkr,              &
                 CZERO,DerBPhiLr,NumRs*kmax_kkr)
!     ----------------------------------------------------------------
      PPr = CZERO
      do klp = 1, kmax_kkr
         mp = mofk(klp)
         klp_bar = bofk(klp)
         cfac = m1m(mp)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               do ir = 1, NumRs
                  PPr(ir,kl1,kl2) = PPr(ir,kl1,kl2) + cfac*DerBPhiLr(ir,kl1,klp)*PhiLr(ir,kl2,klp_bar) &
                                                    + cfac*BPhiLr(ir,kl1,klp)*DerPhiLr(ir,kl2,klp_bar)
               enddo
            enddo
         enddo
      enddo
      do jl3 = 1, jmax_rho
         m3 = mofj(jl3)
         kl3 = kofj(jl3)
         kl3_bar = bofk(kl3)
         do kl2 = 1, kmax_kkr
            do kl1 = 1, kmax_kkr
               cfac1 = cfac0*gaunt(kl1,kl2,kl3)
               cfac2 = cfac0*m1m(m3)*gaunt(kl1,kl2,kl3_bar)
               do ir = 1, NumRs
                  Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)  &
                                            + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
               enddo
            enddo
         enddo
      enddo
   endif
!
!  ===================================================================
!  Get rid of r^2 from Bdensity and Deriv_Bdensity
!  ===================================================================
!  do jl3 = 1, jmax_rho
!     do ir = 1, NumRs
!        Bdensity(ir,jl3) = Bdensity(ir,jl3)/r_mesh(ir)**2
!     enddo
!  enddo
!  if (rad_derivative) then
!     do jl3 = 1, jmax_rho
!        do ir = 1, NumRs
!           Deriv_Bdensity(ir,jl3) = Deriv_Bdensity(ir,jl3)/r_mesh(ir)**2
!        enddo
!     enddo
!  endif
!
   ip = mod(ib-1,NumPEsInEGroup)
   if (NumPEsInEGroup > 1) then
!     ---------------------------------------------------------------
      call bcastMessageInGroup(eGID,Bdensity,NumRs,jmax_rho,ip)
!     ---------------------------------------------------------------
      if (rad_derivative) then
!        ------------------------------------------------------------
         call bcastMessageInGroup(eGID,Deriv_Bdensity,NumRs,jmax_rho,ip)
!        ------------------------------------------------------------
      endif
   endif
!
   nullify(sine_mat, BSinv, smat_inv, Grid, r_mesh)
   nullify(BPhiLr, PhiLr, DerBPhiLr, DerPhiLr, PPr)
   nullify(Bdensity, Deriv_Bdensity)
!
   end subroutine computeResidualRSDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calcDet(mat,kmax,det,diag)
!  ===================================================================
   use KindParamModule, only : RealKind, CmplxKind, IntKind
!
   use MathParamModule, only : CONE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax
   integer (kind=IntKind) :: kl
!
   complex (kind=CmplxKind), intent(out) :: det
   complex (kind=CmplxKind), intent(out) :: diag(kmax)
   complex (kind=CmplxKind), intent(in) :: mat(kmax,kmax)
   complex (kind=CmplxKind) :: matU(kmax,kmax)
!
!  ----------------------------------------------------------
   call zcopy(kmax*kmax,mat,1,matU,1)
!  ----------------------------------------------------------
   call GaussianElim(matU,kmax)
!  ----------------------------------------------------------
   det = CONE
   do kl = 1,kmax
      det = det*matU(kl,kl)
      diag(kl) = matU(kl,kl)
   enddo
!
   end subroutine calcDet
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sortPoleStates(np,ps,idx)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: np
   integer (kind=IntKind), intent(out) :: idx(np)
   integer (kind=IntKind) :: i, j, ki, kj
!
   type (PoleDensityStruct), intent(in) :: ps(np)
!
   do i = 1, np
      idx(i) = i
   enddo
!
   do j = 1, np-1
      kj = idx(j)
      do i = j+1, np
         ki = idx(i)
         if (ps(ki)%PoleE < ps(kj)%PoleE) then
            idx(i) = kj
            kj = ki
         endif
      enddo
      idx(j) = kj
   enddo
!
   end subroutine sortPoleStates
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine swapPoleStates(p1,p2)
!  ===================================================================
   implicit none
!
   type (PoleDensityStruct), intent(inout) :: p1, p2
!
   integer (kind=IntKind) :: n
!
   real (kind=RealKind) :: r
   real (kind=RealKind) :: mat(kmax_kkr_max*kmax_kkr_max)
   real (kind=RealKind) :: dm(MaxNumRs*jmax_rho_max)
!
   n = p1%NumDegens
   p1%NumDegens = p2%NumDegens
   p2%NumDegens = n
!
   r = p1%PoleE
   p1%PoleE = p2%PoleE
   p2%PoleE = r
!
   r = p1%PoleWidth
   p1%PoleWidth = p2%PoleWidth
   p2%PoleWidth = r
!
   r = p1%Qvp
   p1%Qvp = p2%Qvp
   p2%Qvp = r
!
   r = p1%Qmt
   p1%Qmt = p2%Qmt
   p2%Qmt = r
!
   mat = p1%ResidualMat
   p1%ResidualMat = p2%ResidualMat
   p2%ResidualMat = mat
!
   mat = p1%AuxiliaryMat
   p1%AuxiliaryMat = p2%AuxiliaryMat
   p2%AuxiliaryMat = mat
!
   dm = p1%Density
   p1%Density = p2%Density
   p2%Density = dm
!
   if (rad_derivative) then
      dm = p1%Deriv_Density
      p1%Deriv_Density = p2%Deriv_Density
      p2%Deriv_Density = dm
   endif
!
   end subroutine swapPoleStates
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine decodeDOSDATA(NumRs,jmax_rho,dos_data,PoleState)
!  ===================================================================
   use IntegerFactorsModule, only : kofj
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: jmax_rho, NumRs
!
   integer (kind=IntKind) :: n, n0
! 
   complex (kind=CmplxKind), intent(in) :: dos_data(:)
!
   type (PoleDensityStruct), intent(inout) :: PoleState
!
   n0 = 0
   n = NumRs*jmax_rho
!  -------------------------------------------------------------------
   call zcopy(n,dos_data(n0+1:n0+n),1,PoleState%Density,1)
!  -------------------------------------------------------------------
   n0 = n0 + n
   if (rad_derivative) then
!     ----------------------------------------------------------------
      call zcopy(n,dos_data(n0+1:n0+n),1,PoleState%Deriv_Density,1)
!     ----------------------------------------------------------------
      n0 = n0 + n
   endif
   n0 = n0 + 4
!
   end subroutine decodeDOSDATA
!  ===================================================================
end module SMatrixPolesModule
