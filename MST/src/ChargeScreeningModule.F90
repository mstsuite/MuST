module ChargeScreeningModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, TWO, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct

!
public :: initChargeScreeningModule,     &
          calPotentialCorrection,        &
          calCPAEnergyShift,             &
!         calSROEnergyShift,             &
          calChargeCorrection,           &
          getSpeciesPotentialCorrection, &
          getEnergyCorrectionTerm,       &
          endChargeScreeningModule      
!

private
  integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
  integer (kind=IntKind) :: NumSpecies
  integer (kind=IntKind) :: NumSROShells = 2 ! Assume upto 2 neighboring shells

  type ChargeCorrectionData
      real (kind=RealKind) :: fs_radius
      real (kind=RealKind) :: ss_radius
      real (kind=RealKind) :: echarge
      real (kind=RealKind), allocatable :: vmt1_corr(:)
   end type ChargeCorrectionData

  type (ChargeCorrectionData), allocatable :: scr(:)
  real (kind=RealKind), allocatable :: w_ab(:,:,:)

contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initChargeScreeningModule (nlocal, num_atoms)
!  ===================================================================
   use SystemModule, only : getLatticeConstant
!
   use Atom2ProcModule, only : getGlobalIndex
   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies,   &
           getLocalSpeciesContent
!
   use NeighborModule, only : getShellRadius
!
   use ScfDataModule, only : retrieveSROParams, isKKRCPASRO,  &
                   getSpeciesSlope, getSpeciesIntercept

   integer (kind=IntKind), intent(in) :: nlocal, num_atoms
   integer (kind=IntKind) :: i, j, sro_param_num, temp, na
   real (kind=RealKind) :: spec_i, spec_j
   real (kind=RealKind), allocatable :: sro_params(:)
!
   logical :: isWarrenCowley = .false.

   GlobalNumSites = num_atoms
   LocalNumSites  = nlocal

   allocate(scr(LocalNumSites))
   do i = 1, LocalNumSites
      NumSpecies = getLocalNumSpecies(i)
      allocate(scr(i)%vmt1_corr(getLocalNumSpecies(i)))
      scr(i)%fs_radius = getShellRadius(i,1)
!     scr(i)%ss_radius = getLatticeConstant()
      scr(i)%ss_radius = getShellRadius(i,2)
      scr(i)%echarge = ZERO
      do j = 1, getLocalNumSpecies(i)
         scr(i)%vmt1_corr(j) = ZERO
      enddo
   enddo
   allocate(w_ab(NumSpecies,NumSpecies,NumSROShells))

   if (isKKRCPASRO() .eqv. .true.) then
   !  --------------------------------------------------------
      call retrieveSROParams(sro_param_list=sro_params, param_num=sro_param_num, &
                             isWC=isWarrenCowley)
   !  --------------------------------------------------------
!     In the following code, we only consider 1st neighboring shell charge
!     correlation case, will add the 2nd neighboring shell charge correlation
!     later.
!     ================================================================
      w_ab = ZERO
      do na = 1, LocalNumSites
        do i = 1, getLocalNumSpecies(na)
          spec_i = getLocalSpeciesContent(na, i)
          do j = 1, getLocalNumSpecies(na)
            spec_j = getLocalSpeciesContent(na, j)
            if (j < i) then
               w_ab(i,j,1) = (spec_j/spec_i)*w_ab(j,i,1)
            else
               temp = (i - 1)*getLocalNumSpecies(na) - (i - 1)*(i - 2)/2
               if (isWarrenCowley) then
                  w_ab(i,j,1) = spec_j*(ONE - sro_params(temp + j - i + 1))
               else
                  w_ab(i,j,1) = sro_params(temp + j - i + 1)
               endif
            endif
          enddo
        enddo
      enddo
   endif

   end subroutine initChargeScreeningModule
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calPotentialCorrection ()
!  ===================================================================

   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies,   &
                          getLocalSpeciesContent
!
   use ScfDataModule, only : isLinRel, isKKRCPASRO,  getSpeciesSlope, getSpeciesIntercept
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use NeighborModule, only : getShellRadius, getNumAtomsOnShell
!
   use ChargeDistributionModule, only : getGlobalOnSiteElectronTableOld, &
                                     getGlobalTableLine
                                        
   integer (kind=IntKind) :: i, j, ia, na, lig, ja, mig
   real (kind=RealKind) :: qtemp, slope, intercept, dq, avq
   real (kind=RealKind), pointer :: Q_Table(:)
   integer (kind=IntKind), pointer :: global_table_line(:)

   Q_Table => getGlobalOnSiteElectronTableOld()
   global_table_line => getGlobalTableLine()


   do na = 1, LocalNumSites
      j = getGlobalIndex(na)
      do ia = 1, getLocalNumSpecies(na)
        lig = global_table_line(j) + ia
        qtemp = getLocalAtomicNumber(j, ia) - Q_Table(lig)
        if (isKKRCPASRO()) then
           scr(na)%vmt1_corr(ia) = ZERO
           do ns = 1, NumSROShells
              avq = ZERO
              do ja = 1, getLocalNumSpecies(na)
                 mig = global_table_line(j) + ja
                 dq = getLocalAtomicNumber(j,ja) - Q_Table(mig)
                 avq = avq + w_ab(ia,ja,ns)*dq
              enddo
              scr(na)%vmt1_corr(ia) = scr(na)%vmt1_corr(ia)     &
                 -TWO*avq*getNumAtomsOnShell(na,ns)/getShellRadius(na,ns)
           enddo
        else if (isLinRel()) then
           slope = getSpeciesSlope(ia)
           intercept = getSpeciesIntercept(ia)
!          ===========================================================
!          Note: <vmad> will be subtracted from this term when function
!          getSpeciesPotentialCorrection is called from
!          PotentialGenerationModule
!          ===========================================================
           scr(na)%vmt1_corr(ia) = -slope*qtemp + intercept
        else
           scr(na)%vmt1_corr(ia) = TWO*(qtemp/scr(na)%fs_radius)
        endif
     enddo
   enddo

   end subroutine calPotentialCorrection
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calCPAEnergyShift ()
!  ===================================================================
   use MPPModule, only : MyPE
!
   use MathParamModule, only : PI4, THREE, HALF, THIRD, ZERO
!
   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies,  &
                                        getLocalSpeciesContent

   use Atom2ProcModule, only : getGlobalIndex
!
   use ScfDataModule, only : isLinRel, isKKRCPASRO 
!
   use ChargeDistributionModule, only : getGlobalOnSiteElectronTable, &
                                        getGlobalTableLine
!
   use SystemModule, only : getAtomicNumber, getNumAlloyElements, getAlloyElementContent
!
   use SystemVolumeModule, only : getAtomicVPVolume
!
   use NeighborModule, only : getShellRadius, getNumAtomsOnShell
!
   integer (kind=IntKind) :: i, j, ia, na, lig, ja, mig, ns
   real (kind=RealKind) :: qtemp, uc
!  real (kind=RealKind) :: uc, avq
   real (kind=RealKind), pointer :: Q_Table(:)
   integer (kind=IntKind), pointer :: global_table_line(:)
!
   Q_Table => getGlobalOnSiteElectronTable()
   global_table_line => getGlobalTableLine()
!   
   do na = 1, LocalNumSites
      scr(na)%echarge = ZERO
      j = getGlobalIndex(na)
      do ia = 1, getLocalNumSpecies(na)
        lig = global_table_line(j) + ia
        qtemp = getLocalAtomicNumber(j, ia) - Q_Table(lig)
        if (isKKRCPASRO()) then
           uc = ZERO
           do ns = 1, NumSROShells
              avq = ZERO
              do ja = 1, getLocalNumSpecies(na)
                 mig = global_table_line(j) + ja
                 dq = getLocalAtomicNumber(j,ja) - Q_Table(mig)
                 avq = avq + w_ab(ia,ja,ns)*dq
              enddo
              uc = uc + avq*getNumAtomsOnShell(na,ns)/getShellRadius(na,ns)
           enddo
           scr(na)%echarge = scr(na)%echarge + getLocalSpeciesContent(na,ia)*qtemp*uc
        else if (isLinRel()) then
           uc = HALF*qtemp*getSpeciesPotentialCorrection(na,ia,ZERO)
           scr(na)%echarge = scr(na)%echarge - getLocalSpeciesContent(na,ia)*uc
        else
           scr(na)%echarge = scr(na)%echarge -  &
              getLocalSpeciesContent(na, ia)*((qtemp*qtemp)/scr(na)%fs_radius)
        endif
      enddo
   enddo

   end subroutine calCPAEnergyShift
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSROEnergyShift()
!  ===================================================================

   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies,   &
                              getLocalSpeciesContent
   use Atom2ProcModule, only : getGlobalIndex
   use ChargeDistributionModule, only : getGlobalOnSiteElectronTable, &
                                     getGlobalTableLine

   integer (kind=IntKind) :: i, j, ig, ia, na, lig, lig1
   real (kind=RealKind) :: dq_a, dq_b
   real (kind=RealKind), pointer :: Q_Table(:)
   integer (kind=IntKind), pointer :: global_table_line(:)

   Q_Table => getGlobalOnSiteElectronTable()
   global_table_line => getGlobalTableLine()

   do na = 1, LocalNumSites
      scr(na)%echarge = ZERO
      ig = getGlobalIndex(na)
      do i = 1, getLocalNumSpecies(na)
         lig = global_table_line(ig) + i
         dq_a = getLocalAtomicNumber(na, i) - Q_Table(lig)
         do j = 1, getLocalNumSpecies(na)
           lig1 = global_table_line(ig) + j
           dq_b = getLocalAtomicNumber(na, j) - Q_Table(lig1)
           scr(na)%echarge = scr(na)%echarge + &
                   getLocalSpeciesContent(na, i)*w_ab(i,j,1)* &
                  ((dq_b*(dq_b - 2*dq_a))/scr(na)%fs_radius & 
                  -(dq_b - dq_a)**2/scr(na)%ss_radius)
         enddo
      enddo
   enddo

   end subroutine calSROEnergyShift
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calChargeCorrection()
!  ===================================================================

   use ScfDataModule, only : isKKRCPASRO, isChargeCorr

   if (.not.isChargeCorr()) then
      return
!  else if (isKKRCPASRO()) then
!    ! -----------------------------
!      call calPotentialCorrection()
!      call calSROEnergyShift()
!    ! -----------------------------
   else
     ! -----------------------------
       call calPotentialCorrection()
       call calCPAEnergyShift()
     ! -----------------------------
   endif

   end subroutine calChargeCorrection
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSpeciesPotentialCorrection(na, ia, vmad)  result(corr)
!  ===================================================================
   use ScfDataModule, only : isLinRel

   use AtomModule, only : getLocalNumSpecies

   integer (kind=IntKind), intent(in) :: na, ia
   real (kind=RealKind), intent(in) :: vmad
   real (kind=RealKind) :: corr

   if (na < 1 .or. na > LocalNumSites) then
      call ErrorHandler('getSpeciesPotentialCorrection', &
              'Invalid Local Index', na)
   else if (ia < 1 .or. ia > getLocalNumSpecies(na)) then
      call ErrorHandler('getSpeciesPotentialCorrection', &
              'Invalid Species Index', ia)
   endif

   if (isLinRel()) then
      corr = scr(na)%vmt1_corr(ia) - vmad
   else
      corr = scr(na)%vmt1_corr(ia)
   endif

   end function getSpeciesPotentialCorrection
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEnergyCorrectionTerm(na)  result(corr)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: na

   real (kind=RealKind) :: corr

   if (na < 1 .or. na > LocalNumSites) then
      call ErrorHandler('getEnergyCorrection', &
              'Invalid Local Index', na)
   endif

   corr = scr(na)%echarge

   end function getEnergyCorrectionTerm
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endChargeScreeningModule()
!  ===================================================================
   
   integer (kind=IntKind) :: i

   do i = 1, LocalNumSites
      deallocate(scr(i)%vmt1_corr)
   enddo

   deallocate(scr)
   deallocate(w_ab)

   end subroutine endChargeScreeningModule
!  ===================================================================
end module ChargeScreeningModule
