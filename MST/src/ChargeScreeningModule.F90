module ChargeScreeningModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, TWO, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
  use ErrorHandlerModule, only : ErrorHandler, WarningHandler
  use PublicTypeDefinitionsModule, only : NeighborStruct

!
public :: initChargeScreeningModule,     &
          calPotentialCorrection,        &
          calCPAEnergyShift,             &
          calSROEnergyShift,             &
          calChargeCorrection,           &
          getSpeciesPotentialCorrection, &
          getEnergyCorrectionTerm,       &
          endChargeScreeningModule      
!

private
  integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
  integer (kind=IntKind) :: NumSpecies

  type ChargeCorrectionData
      real (kind=RealKind) :: fs_radius
      real (kind=RealKind) :: ss_radius
      real (kind=RealKind) :: echarge
      real (kind=RealKind), allocatable :: vmt1_corr(:)
   end type ChargeCorrectionData

  type (ChargeCorrectionData), allocatable :: scr(:)
  real (kind=RealKind), allocatable :: w_ab(:,:)

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
   use PolyhedraModule, only : getNeighborDistance
!
   use ScfDataModule, only : retrieveSROParams, isKKRCPASRO,  &
                   getSpeciesSlope, getSpeciesIntercept

   integer (kind=IntKind), intent(in) :: nlocal, num_atoms
   integer (kind=IntKind) :: i, j, sro_param_num, temp, na
   real (kind=RealKind) :: spec_i, spec_j
   real (kind=RealKind), allocatable :: sro_params(:)

   GlobalNumSites = num_atoms
   LocalNumSites  = nlocal

   allocate(scr(LocalNumSites))
   do i = 1, LocalNumSites
      NumSpecies = getLocalNumSpecies(i)
      allocate(scr(i)%vmt1_corr(getLocalNumSpecies(i)))
      scr(i)%fs_radius = getNeighborDistance(i, dmin=.true.)
      scr(i)%ss_radius = getLatticeConstant()
      scr(i)%echarge = ZERO
      do j = 1, getLocalNumSpecies(i)
         scr(i)%vmt1_corr(j) = ZERO
      enddo
   enddo
   allocate(w_ab(NumSpecies, NumSpecies))

   if (isKKRCPASRO() .eqv. .true.) then
   !  --------------------------------------------------------
      call retrieveSROParams(sro_param_list=sro_params, param_num=sro_param_num)
   !  --------------------------------------------------------
      w_ab = ZERO
      do na = 1, LocalNumSites
        do i = 1, getLocalNumSpecies(na)
          spec_i = getLocalSpeciesContent(na, i)
          do j = 1, getLocalNumSpecies(na)
            spec_j = getLocalSpeciesContent(na, j)
            if (j < i) then
               w_ab(i, j) = (spec_j/spec_i)*w_ab(j, i)
            else
               temp = (i - 1)*getLocalNumSpecies(na) - (i - 1)*(i - 2)/2
               w_ab(i, j) = sro_params(temp + j - i + 1)
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
   use ScfDataModule, only : isLinRel, getSpeciesSlope, getSpeciesIntercept
   use Atom2ProcModule, only : getGlobalIndex
   use ChargeDistributionModule, only : getGlobalOnSiteElectronTableOld, &
                                     getGlobalTableLine
                                        
   integer (kind=IntKind) :: i, j, ia, na, lig
   real (kind=RealKind) :: qtemp, slope, intercept
   real (kind=RealKind), pointer :: Q_Table(:)
   integer (kind=IntKind), pointer :: global_table_line(:)

   Q_Table => getGlobalOnSiteElectronTableOld()
   global_table_line => getGlobalTableLine()


   do na = 1, LocalNumSites
      j = getGlobalIndex(na)
      do ia = 1, getLocalNumSpecies(na)
        lig = global_table_line(j) + ia
        qtemp = getLocalAtomicNumber(j, ia) - Q_Table(lig)
        if (isLinRel()) then
           slope = getSpeciesSlope(ia)
           intercept = getSpeciesIntercept(ia)
           scr(na)%vmt1_corr(ia) = slope*qtemp + intercept
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

   use AtomModule, only : getLocalAtomicNumber, getLocalNumSpecies,  &
                                        getLocalSpeciesContent

   use Atom2ProcModule, only : getGlobalIndex

   use ChargeDistributionModule, only : getGlobalOnSiteElectronTable, &
                                     getGlobalTableLine

   integer (kind=IntKind) :: i, j, ia, na, lig
   real (kind=RealKind) :: qtemp
   real (kind=RealKind), pointer :: Q_Table(:)
   integer (kind=IntKind), pointer :: global_table_line(:)

   Q_Table => getGlobalOnSiteElectronTable()
   global_table_line => getGlobalTableLine()
   
   do na = 1, LocalNumSites
      scr(na)%echarge = ZERO
      j = getGlobalIndex(na)
      do ia = 1, getLocalNumSpecies(na)
        lig = global_table_line(j) + ia
        qtemp = getLocalAtomicNumber(j, ia) - Q_Table(lig)
        scr(na)%echarge = scr(na)%echarge -  &
           getLocalSpeciesContent(na, ia)*((qtemp*qtemp)/scr(na)%fs_radius)
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
                   getLocalSpeciesContent(na, i)*w_ab(i, j)* &
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

   use ScfDataModule, only : isKKRCPASRO, isKKRCPA, isLinRel

   if (isLinRel() .and. isKKRCPA()) then
     ! -----------------------------
       call calPotentialCorrection()
     ! -----------------------------
   else
     if (isKKRCPA()) then
       ! -----------------------------
         call calPotentialCorrection()
         call calCPAEnergyShift()
       ! -----------------------------
     else if (isKKRCPASRO()) then
       ! -----------------------------
         call calPotentialCorrection()
         call calSROEnergyShift()
       ! -----------------------------
     endif
   endif

   end subroutine calChargeCorrection
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSpeciesPotentialCorrection(na, ia)  result(corr)
!  ===================================================================

   use AtomModule, only : getLocalNumSpecies

   integer (kind=IntKind), intent(in) :: na, ia
   real (kind=RealKind) :: corr

   if (na < 1 .or. na > LocalNumSites) then
      call ErrorHandler('getSpeciesPotentialCorrection', &
              'Invalid Local Index', na)
   else if (ia < 1 .or. ia > getLocalNumSpecies(na)) then
      call ErrorHandler('getSpeciesPotentialCorrection', &
              'Invalid Species Index', ia)
   endif

   corr = scr(na)%vmt1_corr(ia)

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
