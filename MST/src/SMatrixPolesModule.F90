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
                               Ten2m8, FOURTH, CZERO, TEN2m7
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MPPModule, only : syncAllPEs, MyPE
!
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEInGroup
   use GroupCommModule, only : GlobalSumInGroup, bcastMessageInGroup
!
public :: initSMatrixPoles,          &
          endSMatrixPoles,           &
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
          computeResonanceStateDensity,  &
          getResonanceStateDensity,      &
          printSMatrixPoleInfo,      &
          printBoundStateDensity,    &
          isSMatrixPolesInitialized
!
private
   logical :: isInitialized = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: eGID, NumPEsInEGroup, MyPEInEGroup
   integer (kind=IntKind) :: lmax_kkr_max, kmax_kkr_max, kmax_kkr_save
   integer (kind=IntKind) :: lmax_rho_max, kmax_rho_max, jmax_rho_max
   integer (kind=IntKind) :: MaxNumRs
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)
!
   integer (kind=IntKind), parameter :: MaxNumBoundStates = 10
!
   real (kind=RealKind), pointer :: cgnt(:,:,:)
   real (kind=RealKind), allocatable :: gaunt(:,:,:)
!
   complex (kind=CmplxKind), allocatable, target :: wspace0(:), wspace1(:), wspace2(:), wspace3(:)
   complex (kind=CmplxKind), allocatable, target :: wspace4(:)
!
   type PoleStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind) :: NumRs
      integer (kind=IntKind) :: kmax_kkr
      integer (kind=IntKind) :: jmax_rho
      integer (kind=IntKind), allocatable :: NumBoundPoles(:,:)
      integer (kind=IntKind), allocatable :: NumBPDegens(:,:,:)
      integer (kind=IntKind), allocatable :: NumResPoles(:,:)
      integer (kind=IntKind), allocatable :: NumRPDegens(:,:,:)
!
      real (kind=RealKind), allocatable :: ebot(:,:)
      real (kind=RealKind), allocatable :: etop(:,:)
      real (kind=RealKind), allocatable :: BoundPoles(:,:,:)
      real (kind=RealKind), allocatable :: ResPoles(:,:,:)
      real (kind=RealKind), allocatable :: ResWidth(:,:,:)
      real (kind=RealKind), allocatable :: Qvp(:,:,:)
      real (kind=RealKind), allocatable :: Qmt(:,:,:)
!
      complex (kind=CmplxKind), allocatable :: BmatResidual(:,:,:,:)
      complex (kind=CmplxKind), allocatable :: RmatResidual(:,:,:,:)
      complex (kind=CmplxKind), allocatable :: Density(:,:,:)
      complex (kind=CmplxKind), allocatable :: Deriv_Density(:,:,:)
   end type PoleStruct
!
   type (PoleStruct), allocatable :: Pole(:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSMatrixPoles(nla,npola,nspecies,lmax_kkr,lmax_rho,iprint)
!  ===================================================================
   use GauntFactorsModule, only : getK3, getNumK3, getGauntFactor
!
   use RadialGridModule, only : getMaxNumRmesh
!
   use QuadraticMatrixModule, only : initQuadraticMatrix
!
   use IntegerFactorsModule, only : initIntegerFactors,               &
                                    isIntegerFactorsInitialized
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nla, npola, iprint
   integer (kind=IntKind), intent(in) :: nspecies(nla)
   integer (kind=IntKind), intent(in) :: lmax_kkr(nla)
   integer (kind=IntKind), intent(in) :: lmax_rho(nla)
   integer (kind=IntKind) :: id, kl1, kl2, kl3, lmax_max, i2, n
!
   logical, parameter :: isGeneral = .false.
!
   LocalNumAtoms = nla
   n_spin_pola = npola
   print_level = iprint
!
   allocate(Pole(LocalNumAtoms))
!
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
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
!
   do id = 1, LocalNumAtoms
      n = Pole(id)%NumSpecies
      allocate( Pole(id)%NumBoundPoles(n_spin_pola,n) );            Pole(id)%NumBoundPoles = 0
      allocate( Pole(id)%NumBPDegens(kmax_kkr_max,n_spin_pola,n) ); Pole(id)%NumBPDegens = 0
      allocate( Pole(id)%NumResPoles(n_spin_pola,n) );              Pole(id)%NumResPoles = 0
      allocate( Pole(id)%NumRPDegens(kmax_kkr_max,n_spin_pola,n) ); Pole(id)%NumRPDegens = 0
!
      allocate( Pole(id)%BoundPoles(kmax_kkr_max,n_spin_pola,n) )
      allocate( Pole(id)%ResPoles(kmax_kkr_max,n_spin_pola,n) )
      allocate( Pole(id)%ResWidth(kmax_kkr_max,n_spin_pola,n) )
      allocate( Pole(id)%BmatResidual(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max,n_spin_pola,n) )
      allocate( Pole(id)%RmatResidual(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max,n_spin_pola,n) )
      allocate( Pole(id)%Density(MaxNumRs*jmax_rho_max*MaxNumBoundStates,n_spin_pola,n) )
      allocate( Pole(id)%Deriv_Density(MaxNumRs*jmax_rho_max*MaxNumBoundStates,n_spin_pola,n) )
      Pole(id)%Density = CZERO
      Pole(id)%Deriv_Density = CZERO
!
      allocate( Pole(id)%ebot(n_spin_pola,n), Pole(id)%etop(n_spin_pola,n) )
      allocate( Pole(id)%Qvp(MaxNumBoundStates,n_spin_pola,n) )
      allocate( Pole(id)%Qmt(MaxNumBoundStates,n_spin_pola,n) )
      Pole(id)%ebot = ZERO
      Pole(id)%etop = ZERO
      Pole(id)%Qvp = ZERO
      Pole(id)%Qmt = ZERO
   enddo
!
!  -------------------------------------------------------------------
   call initQuadraticMatrix(kmax_kkr_max,isGeneral)
!  -------------------------------------------------------------------
   if (.not.isIntegerFactorsInitialized()) then
!     ----------------------------------------------------------------
      call initIntegerFactors(lmax_max)
!     ----------------------------------------------------------------
   endif
!
   isInitialized = .true.
!
   end subroutine initSMatrixPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSMatrixPoles()
!  ===================================================================
   use QuadraticMatrixModule, only : endQuadraticMatrix
!
   implicit none
!
   integer (kind=IntKind) :: id
!
   deallocate( gaunt )
   deallocate( wspace0 )
   deallocate( wspace1 )
   deallocate( wspace2 )
   deallocate( wspace3 )
   deallocate( wspace4 )
!
   do id = 1, LocalNumAtoms
      deallocate( Pole(id)%NumBoundPoles )
      deallocate( Pole(id)%NumBPDegens )
      deallocate( Pole(id)%NumResPoles )
      deallocate( Pole(id)%NumRPDegens )
      deallocate( Pole(id)%BoundPoles )
      deallocate( Pole(id)%ResPoles )
      deallocate( Pole(id)%ResWidth )
      deallocate( Pole(id)%BmatResidual )
      deallocate( Pole(id)%RmatResidual )
      deallocate( Pole(id)%Density )
      deallocate( Pole(id)%Deriv_Density )
      deallocate( Pole(id)%ebot, Pole(id)%etop )
      deallocate( Pole(id)%Qvp, Pole(id)%Qmt )
   enddo
   deallocate(Pole)
!
!  -------------------------------------------------------------------
   call endQuadraticMatrix()
!  -------------------------------------------------------------------
!
   isInitialized = .false.
!
   end subroutine endSMatrixPoles
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
   function getNumBoundStateDegen(id,ia,is,ib) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ib, is, id, ia
   integer (kind=IntKind) :: n
!
   n = Pole(id)%NumBPDegens(ib,is,ia)
!
   end function getNumBoundStateDegen
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBoundStateEnergy(id,ia,is,ib) result(e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ib, ia
!
   real (kind=RealKind) :: e
!
   e = Pole(id)%BoundPoles(ib,is,ia)
!
   end function getBoundStateEnergy
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBoundStateChargeInCell(id,ia,is,ib,qmt) result(qvp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ib, ia
   real (kind=RealKind), intent(out), optional :: qmt
!
   real (kind=RealKind) :: qvp
!
   qvp = Pole(id)%Qvp(ib,is,ia)
!
   if (present(qmt)) then
      qmt = Pole(id)%Qmt(ib,is,ia)
   endif
!
   end function getBoundStateChargeInCell
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBoundStateDensity(id,ia,is,ib,NumRs,jmax_rho,derivative) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia, is, ib
   integer (kind=IntKind), intent(out), optional :: NumRs, jmax_rho
   integer (kind=IntKind) :: NumBPs
!
   logical, intent(in), optional :: derivative
   logical :: deriv_den
!
   complex (kind=CmplxKind), pointer :: Bdensity(:,:,:)
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
   else if (ib < 1 .or. ib > NumBPs) then
      call ErrorHandler('getBoundStateDensity','Invalid bound state index',ib)
   endif
!
   if (present(derivative)) then
      deriv_den = derivative
   else
      deriv_den = .false.
   endif
!
   if (deriv_den) then
      Bdensity => aliasArray3_c(Pole(id)%Deriv_Density(:,is,ia),   &
                                Pole(id)%NumRs,Pole(id)%jmax_rho,NumBPs)
   else
      Bdensity => aliasArray3_c(Pole(id)%Density(:,is,ia),         &
                                Pole(id)%NumRs,Pole(id)%jmax_rho,NumBPs)
   endif
   p => Bdensity(:,:,ib)
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
   function getNumResonanceStateDegen(id,ia,is,ib) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ib, is, id, ia
   integer (kind=IntKind) :: n
!
   n = Pole(id)%NumRPDegens(ib,is,ia)
!
   end function getNumResonanceStateDegen
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getResonanceStateEnergy(id,ia,is,ib,hw) result(e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ib, is, id, ia
!
   real (kind=RealKind), optional :: hw ! Half width
   real (kind=RealKind) :: e
!
   e = Pole(id)%ResPoles(ib,is,ia)
   hw = Pole(id)%ResWidth(ib,is,ia)
!
   end function getResonanceStateEnergy
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isEnergyInResonanceRange(e,site,atom,spin,res_id) result(y)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, spin
   integer (kind=IntKind), intent(out), optional :: res_id
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
      el = Pole(site)%ResPoles(i,spin,atom) - HALF*Pole(site)%ResWidth(i,spin,atom)
      er = Pole(site)%ResPoles(i,spin,atom) + HALF*Pole(site)%ResWidth(i,spin,atom)
      if (e >= el .and. e <= er) then
         j = i
         y = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
!
   if (present(res_id)) then
      res_id = j
   endif
!
   end function isEnergyInResonanceRange
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeResonanceStateDensity(site,atom,spin)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, spin, atom
!
!
   end subroutine computeResonanceStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getResonanceStateDensity(site,atom,spin,rstate,peak_term,e, &
                                     NumRs,jmax_rho,derivative) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: site, atom, spin, rstate
   integer (kind=IntKind), intent(out), optional :: NumRs, jmax_rho
   integer (kind=IntKind) :: NumRPs
!
   logical, intent(in) :: peak_term
   logical, intent(in), optional :: derivative
   logical :: deriv_den
!
   real (kind=RealKind), intent(in), optional :: e
!
   complex (kind=CmplxKind), pointer :: Rdensity(:,:,:)
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
   NumRPs = Pole(site)%NumResPoles(spin,atom)
   if (NumRPs < 1) then
      nullify(p)
      return
   else if (rstate < 1 .or. rstate > NumRPs) then
      call ErrorHandler('getBoundStateDensity','Invalid resonance state index',rstate)
   endif
!
   if (present(derivative)) then
      deriv_den = derivative
   else
      deriv_den = .false.
   endif
!
   if (deriv_den) then
      Rdensity => aliasArray3_c(Pole(site)%Deriv_Density(:,spin,atom),   &
                                Pole(site)%NumRs,Pole(site)%jmax_rho,NumRPs)
   else
      Rdensity => aliasArray3_c(Pole(site)%Density(:,spin,atom),         &
                                Pole(site)%NumRs,Pole(site)%jmax_rho,NumRPs)
   endif
   p => Rdensity(:,:,rstate)
!
   end function getResonanceStateDensity
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSMatrixPoleInfo(id,ia,is)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, id, ia
   integer (kind=IntKind) :: ib
!
   write(6,'(3(a,i3))')'Spin index: ',is,',  Site index: ',id,',  Species index: ',ia
   if (Pole(id)%ebot(is,ia) < -TEN2m8) then
      write(6,'(a,f13.8,a,f13.8,a,i5)')'The Number of bound states found within (',  &
                                       Pole(id)%ebot(is,ia),', ',ZERO,'): ',Pole(id)%NumBoundPoles(is,ia)
      do ib = 1, Pole(id)%NumBoundPoles(is,ia)
         write(6,'(a,i2,5x,a,f20.12)')'Degeneracy = ',Pole(id)%NumBPDegens(ib,is,ia), &
                                      ', Bound state energy = ',Pole(id)%BoundPoles(ib,is,ia)
      enddo
   endif
   if (Pole(id)%etop(is,ia) > TEN2m8) then
      write(6,'(/,a,f13.8,a,f13.8,a,i5)')'The Number of resonance states found within (',  &
                                         Pole(id)%ebot(is,ia),', ',Pole(id)%etop(is,ia),'): ', &
                                         Pole(id)%NumResPoles(is,ia)
      do ib = 1, Pole(id)%NumResPoles(is,ia)
         write(6,'(a,i2,5x,a,f20.12,a,f20.12)')'Degeneracy = ',Pole(id)%NumRPDegens(ib,is,ia), &
                                               ', Resonance state energy = ',Pole(id)%ResPoles(ib,is,ia), &
                                               ', Width = ',Pole(id)%ResWidth(ib,is,ia)
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
   integer (kind=IntKind) :: ib, ir, jl, jmax_rho, kmax_rho, NumRs, occ, NumBPs
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
   complex (kind=CmplxKind), pointer :: Bdensity(:,:,:)
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
   Bdensity => aliasArray3_c(Pole(id)%Density(:,is,ia),NumRs,jmax_rho,NumBPs)
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
   do ib = 1, NumBPs
!     ================================================================
!     Note: BmatResidual already contains a summation of the
!           contributions from the degenerate poles.
!     ================================================================
      occ = (3-n_spin_pola)*Pole(id)%NumBPDegens(ib,is,ia)
      deg = real(Pole(id)%NumBPDegens(ib,is,ia),kind=RealKind)
      q_VP=getVolumeIntegration(id,NumRs,r_mesh,kmax_rho,jmax_rho,0,Bdensity(:,:,ib),q_MT)
      write(6,'(f15.8,4x,i5,8x,i5,4x,2f10.5)')                        &
            Pole(id)%BoundPoles(ib,is,ia),Pole(id)%NumBPDegens(ib,is,ia),occ,q_VP/deg,q_MT/deg
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
         rho0(ir) = real(Bdensity(ir,1,ib),kind=RealKind)*rfac
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
            if (abs(Bdensity(ir,jl,ib)) > TEN2m6) then
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
            call writeFunction(file_name,NumRs,r_mesh,Bdensity(:,jl,ib))
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
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix
   use QuadraticMatrixModule, only : solveQuadraticEquation, getEigenValue,   &
                                     solveLinearEquation, getEigenVector,     &
                                     getEigenMatrix
!
   use WriteMatrixModule,  only : writeMatrix
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
   integer (kind=IntKind) :: MyNumWindows, nbr(2), nb0, nr0, ib, ir
   integer (kind=IntKind) :: bpdeg(kmax_kkr_max), rpdeg(kmax_kkr_max), degens(kmax_kkr_max)
!
   logical :: isZeroInterval = .false.
   logical :: chkpole = .false.
!
   real (kind=RealKind) :: WindowWidth, ResWidth
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0, err, w
!
   real (kind=RealKind) :: bpe(kmax_kkr_max), ebr(kmax_kkr_max), bpe_prev, rpe_prev
!
   complex (kind=CmplxKind) :: e, a2l, a2lp, c, cde, ce0, det
   complex (kind=CmplxKind), pointer :: jost_mat(:,:)
   complex (kind=CmplxKind), pointer :: s0(:,:), s1(:,:), s2(:,:)
   complex (kind=CmplxKind), pointer :: sm(:,:)
   complex (kind=CmplxKind), pointer :: pv(:), evr(:), evl(:), em(:,:)
   complex (kind=CmplxKind) :: vt(kmax_kkr_max), diag(kmax_kkr_max)
   complex (kind=CmplxKind) :: rpe(kmax_kkr_max)
   complex (kind=CmplxKind) :: erc(kmax_kkr_max)
   complex (kind=CmplxKind) :: bmat(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max)
   complex (kind=CMplxKind) :: rmat(kmax_kkr_max*kmax_kkr_max,kmax_kkr_max)
!
   if (eb > et) then
      call ErrorHandler('findSMatrixPoles','eb > et',eb,et)
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
   else
      ResWidth = 1.0d20
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
      NumWindows = int((et-eb)/WindowWidth)
      NumWindows = NumWindows - mod(NumWindows,4)
      if (NumWindows < NumPEsInEGroup) then
         NumWindows = NumPEsInEGroup
         MyNumWindows = 1
      else
         NumWindows = ceiling(NumWindows/real(NumPEsInEGroup))*NumPEsInEGroup
         MyNumWindows = NumWindows/NumPEsInEGroup
      endif
      WindowWidth = (et-eb)/real(NumWindows,kind=RealKind)
!     de = Delta
      de = WindowWidth/4.0d0
   endif
!
   de2 = de*TWO; dede2 = de*de*TWO
   cde = de
!
   s0 => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   s1 => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
   s2 => aliasArray2_c(wspace2,kmax_kkr,kmax_kkr)
   sm => aliasArray2_c(wspace3,kmax_kkr,kmax_kkr)
!
   Pole(id)%NumBoundPoles(is,ia) = 0; Pole(id)%NumBPDegens(:,is,ia) = 0
   Pole(id)%NumResPoles(is,ia) = 0; Pole(id)%NumRPDegens(:,is,ia) = 0
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
         if (abs(aimag(pv(ie))) < Ten2m8 .and. real(pv(ie),kind=RealKind)+e0 < ZERO) then   ! Bound states
            pe = real(pv(ie),kind=RealKind) + e0
            if (pe >= w0 .and. pe <= w0+WindowWidth) then
!              -------------------------------------------------------
               em => getEigenMatrix(ie) ! em is the residule matrix of
                                        ! integrating sm^{-1} around its eigenvalue
!              -------------------------------------------------------
               if (size(em,1) /= kmax_kkr) then
                  call ErrorHandler('findSMatrixPoles','inconsistent matrix size',size(em,1),kmax_kkr)
               endif
               if (abs(pe-bpe_prev) > TEN2m6) then
                  nb = nb + 1
!                 Pole(id)%BoundPoles(nb,is,ia) = pe
                  bpe(nb) = pe
                  bpdeg(nb) = 1
!write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
!                 ----------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,em,1,bmat(1,nb),1)
!                 ----------------------------------------------------
                  bpe_prev = pe
               else if (nb == 0) then
                  call ErrorHandler('findSMatrixPoles','bound state pe = ZERO',pe)
               else ! In degeneracy case, em is added to bmat of the same energy
                  bpdeg(nb) = bpdeg(nb) + 1
!                 ----------------------------------------------------
                  call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,bmat(1,nb),1)
!                 ----------------------------------------------------
               endif
!write(6,'(a,3i4,a,f12.8)')'MyPEinEGroup,nb,deg,pe = ',MyPEinEGroup,nb,bpdeg(nb),', ',pe
            endif
         else if (aimag(sqrt(pv(ie)+e0)) < ZERO) then  ! Resonance states
            pe = real(pv(ie),kind=RealKind) + e0
            w = aimag(sqrt(pv(ie)))**2
            if (pe >= w0 .and. pe <= w0+WindowWidth .and. pe > ZERO .and. TWO*w <= ResWidth) then
!              -------------------------------------------------------
               em => getEigenMatrix(ie) ! em is the residule matrix of
                                        ! integrating sm^{-1} around its eigenvalue
!              -------------------------------------------------------
               if (abs(pe-rpe_prev) > TEN2m6) then
                  nr = nr + 1
!                 Pole(id)%ResPoles(nr,is,ia) = pe
!                 Pole(id)%ResWidth(nr,is,ia) = w
                  rpe(nr) = cmplx(pe,aimag(sqrt(pv(ie)))**2)
                  rpdeg(nr) = 1
!write(6,'(a,2d15.8,a,2d15.8)')'Pole = ',pv(ie)+e0,', kappa = ',sqrt(pv(ie)+e0)
!                 ----------------------------------------------------
                  call zcopy(kmax_kkr*kmax_kkr,em,1,rmat(1,nr),1)
!                 ----------------------------------------------------
                  rpe_prev = pe
               else if (nr == 0) then
                  call ErrorHandler('findSMatrixPoles','resonance state pe = ZERO',pe)
               else
                  rpdeg(nr) = rpdeg(nr) + 1
!                 ----------------------------------------------------
                  call zaxpy(kmax_kkr*kmax_kkr,CONE,em,1,rmat(1,nr),1)
!                 ----------------------------------------------------
               endif
            endif
         endif
      enddo
   enddo
!
   if (chkpole) then
      do je = 1, nb
!        e0 = Pole(id)%BoundPoles(je,is,ia)
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
!        e0 = Pole(id)%ResPoles(je,is,ia)
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
   ebr = ZERO; erc = CZERO
   do ip = 1, NumPEsInEGroup
      if (MyPEinEGroup == ip-1) then
         nbr(1) = nb
         nbr(2) = nr
      endif
      call bcastMessageInGroup(eGID,nbr,2,ip-1)
      if (nbr(1) > 0) then
         nb0 = Pole(id)%NumBoundPoles(is,ia)
         if (MyPEinEGroup == ip-1) then
            ebr(1:nbr(1)) = bpe(1:nbr(1))
!           ----------------------------------------------------------
!           call zcopy(kmax_kkr_max*kmax_kkr_max*kmax_kkr_max,bmat,1,brmat,1)
!           do ib = 1, nbr(1)
!              -------------------------------------------------------
!              call zcopy(kmax_kkr_max*kmax_kkr_max,bmat(1,ib),1,Pole(id)%BmatResidual(1,nb0+ib,is,ia),1)
!              -------------------------------------------------------
!           enddo
!           ----------------------------------------------------------
            call zcopy(kmax_kkr_max*kmax_kkr_max*nbr(1),bmat,1,Pole(id)%BmatResidual(1,nb0+1,is,ia),1)
!           ----------------------------------------------------------
            degens(1:nbr(1)) = bpdeg(1:nbr(1))
         endif
!        -------------------------------------------------------------
         call bcastMessageInGroup(eGID,ebr,nbr(1),ip-1)
!        call bcastMessageInGroup(eGID,brmat,kmax_kkr_max*kmax_kkr_max,kmax_kkr_max,ip-1)
         call bcastMessageInGroup(eGID,degens,nbr(1),ip-1)
!        -------------------------------------------------------------
         do ib = 1, nbr(1)
            Pole(id)%BoundPoles(nb0+ib,is,ia) = ebr(ib)
            Pole(id)%NumBPDegens(nb0+ib,is,ia) = degens(ib)
!           ----------------------------------------------------------
!           call zcopy(kmax_kkr*kmax_kkr,brmat(1,ib),1,          &
!                      Pole(id)%BmatResidual(1,nb0+ib,is,ia),1)
!           ----------------------------------------------------------
         enddo
!        -------------------------------------------------------------
         call bcastMessageInGroup(eGID,Pole(id)%BmatResidual(1:,nb0+1:,is,ia),&
                                  kmax_kkr_max*kmax_kkr_max,nbr(1),ip-1)
!        -------------------------------------------------------------
      endif
      if (nbr(2) > 0) then
         nr0 = Pole(id)%NumResPoles(is,ia)
         if (MyPEinEGroup == ip-1) then
            erc(1:nbr(2)) = rpe(1:nbr(2))
!           ----------------------------------------------------------
!           call zcopy(kmax_kkr_max*kmax_kkr_max*kmax_kkr_max,rmat,1,brmat,1)
            call zcopy(kmax_kkr_max*kmax_kkr_max*nbr(2),rmat,1,Pole(id)%RmatResidual(1,nr0+1,is,ia),1)
!           ----------------------------------------------------------
            degens(1:nbr(2)) = rpdeg(1:nbr(2))
         endif
!        -------------------------------------------------------------
         call bcastMessageInGroup(eGID,erc,nbr(2),ip-1)
!        call bcastMessageInGroup(eGID,brmat,kmax_kkr_max*kmax_kkr_max,nbr(2),ip-1)
         call bcastMessageInGroup(eGID,degens,nbr(2),ip-1)
!        -------------------------------------------------------------
         do ir = 1, nbr(2)
            Pole(id)%ResPoles(nr0+ir,is,ia) = real(erc(ir),kind=RealKind)
            Pole(id)%ResWidth(nr0+ir,is,ia) = TWO*aimag(erc(ir))
            Pole(id)%NumRPDegens(nr0+ir,is,ia) = degens(ir)
!           ----------------------------------------------------------
!           call zcopy(kmax_kkr*kmax_kkr,brmat(1,ir),1,          &
!                      Pole(id)%RmatResidual(1,nr0+ir,is,ia),1)
!           ----------------------------------------------------------
         enddo
!        -------------------------------------------------------------
         call bcastMessageInGroup(eGID,Pole(id)%RmatResidual(1:,nr0+1:,is,ia),&
                                  kmax_kkr_max*kmax_kkr_max,nbr(2),ip-1)
!        -------------------------------------------------------------
      endif
      Pole(id)%NumBoundPoles(is,ia) = Pole(id)%NumBoundPoles(is,ia) + nbr(1)
      Pole(id)%NumResPoles(is,ia) = Pole(id)%NumResPoles(is,ia) + nbr(2)
   enddo
!
   nullify(pv, evl, evr, em)
   nullify(s0, s1, s2, sm)
!
   end subroutine findSMatrixPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeBoundStateDensity(id,ia,is)
!  ===================================================================
   use SSSolverModule, only : solveSingleScattering, getSineMatrix
   use SSSolverModule, only : getRegSolution, getSolutionRmeshSize
   use SSSolverModule, only : getRegSolutionDerivative
!  use SSSolverModule, only : getJostInvMatrix, getOmegaHatMatrix
!  use SSSolverModule, only : computeDOS, getDOS
!
   use RadialGridModule, only : getGrid
!
   use MatrixModule, only : computeAStarTInv
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use IBZRotationModule, only : symmetrizeMatrix
!
   use IntegerFactorsModule, only : lofj, mofj, kofj, lofk, mofk, bofk, m1m
!
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ia
!
   integer (kind=IntKind) :: ie, ib, info
   integer (kind=IntKind) :: kmax_kkr, jmax_rho, kmax_rho, NumRs, NumBPs
   integer (kind=IntKind) :: kl, klp, klp_bar, kl1, kl2, kl3, kl3_bar, m3, mp, ir, jl3
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: Bdensity(:,:,:)
   complex (kind=CmplxKind), pointer :: Deriv_Bdensity(:,:,:)
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
   NumRs = getSolutionRmeshSize()
   Pole(id)%NumRs = NumRs
   NumBPs = Pole(id)%NumBoundPoles(is,ia)
!
   if (NumBPs > 0) then
      Bdensity => aliasArray3_c(Pole(id)%Density(:,is,ia),NumRs,jmax_rho,NumBPs)
      Deriv_Bdensity => aliasArray3_c(Pole(id)%Deriv_Density(:,is,ia),NumRs,jmax_rho,NumBPs)
   else
      return
   endif
!
!  ===================================================================
!  Note: Bdensity stores the density multiplied the number of degeneracies
!        of the bound state.
!  ===================================================================
   Bdensity = CZERO
   Deriv_Bdensity = CZERO
   smat_inv => aliasArray2_c(wspace0,kmax_kkr,kmax_kkr)
   BSinv => aliasArray2_c(wspace1,kmax_kkr,kmax_kkr)
!
   Grid => getGrid(id)
   r_mesh => Grid%r_mesh
   BPhiLr => aliasArray3_c(wspace2,NumRs,kmax_kkr,kmax_kkr)
   PPr => aliasArray3_c(wspace3,NumRs,kmax_kkr,kmax_kkr)
   DerBPhiLr => aliasArray3_c(wspace4,NumRs,kmax_kkr,kmax_kkr)
!
   do ib = MyPEinEGroup+1, NumBPs, NumPEsInEGroup
      e = Pole(id)%BoundPoles(ib,is,ia)
      kappa = sqrt(e)
      cfac0 = HALF*kappa
!     ----------------------------------------------------------------
      call solveSingleScattering(is, id, e, CZERO, atom=ia)
!     ----------------------------------------------------------------
      sine_mat => getSineMatrix()
!jost_inv => getJostInvMatrix()
!     ================================================================
!     calculate sine_mat^(-T*) and store the result in smat_inv
!     ----------------------------------------------------------------
      call computeAStarTInv(sine_mat,kmax_kkr,kmax_kkr,smat_inv)
!     ----------------------------------------------------------------
!
!     ===================================================================
!     calculate BmatResidual*sine_mat^{-T*} and store the result in BSinv
!     ----------------------------------------------------------------
      call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,             &
                 Pole(id)%BmatResidual(1,ib,is,ia),kmax_kkr,smat_inv,kmax_kkr, &
                 CZERO,BSinv,kmax_kkr)
!     ----------------------------------------------------------------
      call symmetrizeMatrix(BSinv,kmax_kkr)
!     ----------------------------------------------------------------
!call zgemm('n','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,jost_inv,kmax_kkr,smat_inv,kmax_kkr,CZERO,BSinv,kmax_kkr)
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
                  Bdensity(ir,jl3,ib) = Bdensity(ir,jl3,ib)           &
                                      + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
!+ cfac1*PPr(ir,kl1,kl2) - conjg(cfac2*PPr(ir,kl1,kl2))
               enddo
            enddo
         enddo
      enddo
!
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
                  Deriv_Bdensity(ir,jl3,ib) = Deriv_Bdensity(ir,jl3,ib)  &
                                            + cfac1*PPr(ir,kl1,kl2) + conjg(cfac2*PPr(ir,kl1,kl2))
!+ cfac1*PPr(ir,kl1,kl2) - conjg(cfac2*PPr(ir,kl1,kl2))
               enddo
            enddo
         enddo
      enddo
!
!     ================================================================
!     Get rid of r^2 from Bdensity and Deriv_Bdensity
!     ================================================================
      do jl3 = 1, jmax_rho
         do ir = 1, NumRs
            Bdensity(ir,jl3,ib) = Bdensity(ir,jl3,ib)/r_mesh(ir)**2
         enddo
         do ir = 1, NumRs
            Deriv_Bdensity(ir,jl3,ib) = Deriv_Bdensity(ir,jl3,ib)/r_mesh(ir)**2
         enddo
      enddo
!
!call computeDOS()
!dos_r_jl => getDOS()
!do jl3 = 1, jmax_rho
!non_zero = .false.
!LOOP_ir_3: do ir = 1, Grid%jend
!if (abs(dos_r_jl(ir,jl3)) > TEN2m6) then
!non_zero = .true.
!exit LOOP_ir_3
!endif
!enddo LOOP_ir_3
!if (non_zero) then
!write(6,'(a,3i5)')'None zero dos_r_jl component: ib,lp,mp = ',ib,lofj(jl3),mofj(jl3)
!write(6,'(a,2d15.8,2x,2d15.8)')'den,dos = ',Bdensity(1000,jl3,ib),dos_r_jl(1000,jl3)
!endif
!enddo
   enddo
!  ------------------------------------------------------------------
   call GlobalSumInGroup(eGID,Bdensity,NumRs,jmax_rho,NumBPs)
   call GlobalSumInGroup(eGID,Deriv_Bdensity,NumRs,jmax_rho,NumBPs)
!  ------------------------------------------------------------------
!
   do ib = 1, NumBPs
      Pole(id)%Qvp(ib,is,ia) = getVolumeIntegration(id,NumRs,r_mesh,kmax_rho,    &
                                                    jmax_rho,0,Bdensity(:,:,ib), &
                                                    Pole(id)%Qmt(ib,is,ia))
   enddo
!
   nullify(sine_mat, BSinv, smat_inv, Grid, r_mesh)
   nullify(BPhiLr, PhiLr, DerBPhiLr, DerPhiLr, PPr)
   nullify(Bdensity, Deriv_Bdensity)
!
   end subroutine computeBoundStateDensity
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
end module SMatrixPolesModule
