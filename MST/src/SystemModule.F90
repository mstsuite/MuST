module SystemModule
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : ZERO, ONE
   use ChemElementModule, only : MaxLenOfAtomName
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initSystem,             &
          endSystem,              &
          getSystemID,            &
          getSystemTitle,         &
          getNumAtoms,            &
          getNumVacancies,        &
          getNumCoherentPots,     &
          getNumAtomTypes,        &
          getAtomType,            &
          getAtomTypeName,        &
          getNumAtomsOfType,      &
          getAlloyTableSize,      &
          getNumAlloyElements,    &
          getAlloyElementName,    &
          getAlloyElementContent, &
          getAlloySpeciesIndex,   &
          setBravaisLattice,      &
          getBravaisLattice,      &
          setSystemTable,         &
          getAtomName,            &
          getAtomicNumber,        &
          getAtomPosition,        &
          setAtomPosition,        &
          getSiteLIZ,             &
          setSiteLIZ,             &
          resetSiteLIZ,           &
          getScalingFactor,       &
          getMomentDirection,     &
          setMomentDirection,     &
          getMomentDirectionOld,  &
          setMomentDirectionOld,  &
          getConstrainField,      &
          setConstrainField,      &
          getForce,               &
          setForce,               &
          getRMSInfo,             &
          getAtomEnergy,          &
          setAtomEnergy,          &
          setRMSInfo,             &
          setLatticeConstant,     &
          getLatticeConstant,     &
          getRadicalPlaneRatio,   &
          getLmaxKKR,             &
          getLmaxPhi,             &
          getLmaxRho,             &
          getLmaxPot,             &
          getLmaxMax,             &
          getUniformGridParam,    &
          printSystem,            &
          updateSystem,           &
          getMomentDirectionMixingParam, &
          writeAtomPositionData,  &
          writeAtomPositionMovie, &
          writeMomentDirectionData, &
          writeMomentMovie,       &
          getAdditionalElectrons, &
          updateSystemMovie,      &
          resetSystemMovie
!
   interface getLmaxKKR
      module procedure getLmaxKKR0, getLmaxKKR1
   end interface getLmaxKKR
!
   interface getLmaxPhi
      module procedure getLmaxPhi0, getLmaxPhi1
   end interface getLmaxPhi
!
   interface getLmaxPot
      module procedure getLmaxPot0, getLmaxPot1
   end interface getLmaxPot
!
   interface getLmaxRho
      module procedure getLmaxRho0, getLmaxRho1
   end interface getLmaxRho
!
   interface getAtomType
      module procedure getAtomType_0, getAtomType_1
   end interface getAtomType
!
   interface getSiteLIZ
      module procedure getSiteLIZ_one, getSiteLIZ_all
   end interface getSiteLIZ
!
   interface getAtomName
      module procedure getAtomName_one, getAtomName_all
   end interface getAtomName
!
   interface getAtomicNumber
      module procedure getAtomicNumber_one, getAtomicNumber_all
   end interface getAtomicNumber
!
   interface getAtomPosition
      module procedure getAtomPosition_one, getAtomPosition_all
   end interface getAtomPosition
!
   interface getForce
      module procedure getForce_one, getForce_all
   end interface getForce
!
   interface getMomentDirection
      module procedure getMomentDirection_one, getMomentDirection_all
   end interface getMomentDirection
!
   interface getMomentDirectionOld
      module procedure getMomentDirectionOld_one, getMomentDirectionOld_all
   end interface getMomentDirectionOld
!
   interface getConstrainField
      module procedure getConstrainField_one, getConstrainField_all
   end interface getConstrainField
!
   integer (kind=IntKind), public :: printSysMovie
private
   integer (kind=IntKind), parameter :: MaxSpinIndex = 2
!
   character (len=50) :: SystemID
   character (len=70) :: SystemTitle
   character (len=45) :: fposi_in
   character (len=45) :: fposi_out
   character (len=45) :: fevec_in
   character (len=45) :: fevec_out
   character (len=50) :: file_path
   character (len=1), pointer :: AtomNameCharacter(:)
   character (len=MaxLenOfAtomName), target, allocatable :: AtomName(:)
   character (len=MaxLenOfAtomName), allocatable :: AtomTypeName(:)
!
   integer (kind=IntKind) :: table_size
   integer (kind=IntKind), allocatable :: table_line(:)
   character (len=MaxLenOfAtomName), allocatable :: AlloyElementName(:)
   integer (kind=IntKind), allocatable :: AlloyElementAN(:)
   integer (kind=IntKind), allocatable :: NumAlloyElements(:)
   integer (kind=IntKind) :: MaxAlloyElements
   real (kind=RealKind), allocatable :: AlloyElementContent(:)
!
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: NumVacancies
   integer (kind=IntKind) :: NumCoherentPots
   integer (kind=IntKind) :: NumAtomTypes
   integer (kind=IntKind) :: uniform_grid(3)
   integer (kind=IntKind) :: nspin
   integer (kind=IntKind), pointer :: AtomicNum(:)
   integer (kind=IntKind), allocatable, target :: AtomType(:)
   integer (kind=IntKind), allocatable :: NumAtomsOfType(:)
   integer (kind=IntKind), allocatable :: LmaxKKR(:)
   integer (kind=IntKind), allocatable :: LmaxPhi(:)
   integer (kind=IntKind), allocatable :: LmaxRho(:)
   integer (kind=IntKind), allocatable :: LmaxPot(:)
   integer (kind=IntKind), allocatable, target :: SiteLIZ(:)
   integer (kind=IntKind) :: LmaxMax
   integer (kind=IntKind) :: LmaxRhoMax
   integer (kind=IntKind) :: LmaxPotMax
   integer (kind=IntKind) :: LmaxKKRMax
   integer (kind=IntKind) :: LmaxPhiMax
!
   real (kind=RealKind) :: scaling
   real (kind=RealKind) :: alat
   real (kind=RealKind), pointer :: Bravais(:,:)
   real (kind=RealKind), pointer :: AtomPosition(:,:)
   real (kind=RealKind), pointer :: EnPres(:,:)
   real (kind=RealKind), pointer :: Force(:,:)
   real (kind=RealKind), pointer :: RMS(:,:)
   real (kind=RealKind), pointer :: Evec(:,:)
   real (kind=RealKind), pointer :: Evec_Old(:,:)
   real (kind=RealKind), pointer :: ConstrainField(:,:)
   real (kind=RealKind), pointer :: moment(:)
   real (kind=RealKind), pointer :: emix(:)
   real (kind=RealKind), allocatable :: RadicalPlaneRatio(:)
!
   real (kind=RealKind) :: AdditionalElectrons = ZERO ! additional electrons artificially added to the system
!
   integer (kind=IntKind) :: MessageID = 9000
   integer (kind=IntKind) :: NumPositionUpdates = 0
   integer (kind=IntKind) :: NumEvecUpdates = 0
   integer (kind=IntKind) :: NumEvecOldUpdates = 0
   integer (kind=IntKind) :: NumConstrFieldUpdates = 0
   integer (kind=IntKind) :: NumFoceVectUpdates = 0
   integer (kind=IntKind) :: NumEnergyUpdates =0
   integer (kind=IntKind) :: NumRMSUpdates = 0
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSystem(tbl_id)
!  ===================================================================
   use MPPModule, only : MyPE
   use ChemElementModule, only : getZtot, getZval, getName
!
   use MathParamModule, only : TEN2m8, THIRD, ONE
!
   use DataServiceCenterModule, only : isDataStorageExisting,         &
                                       createDataStorage,             &
                                       getDataStorageSize,            &
                                       getDataStorageLDA,             &
                                       getDataStorage,                &
                                       RealType, IntegerType, CharacterType, &
                                       IntegerMark, RealMark, CharacterMark
!
   use InputModule, only : getKeyValue, getTableIndex, getKeyIndexValue
   use InputModule, only : isKeyExisting
!
   use PublicParamDefinitionsModule, only : ASA, MuffinTin, MuffinTinASA
!
   implicit none
!
   character (len=50) :: info_table
   character (len=80) :: svalue
   character (len=80), allocatable :: value(:)
!
   integer (kind=IntKind), intent(in) :: tbl_id
   integer (kind=IntKind) :: i, info_id, ig, j, reset_lmax, pot_type
   integer (kind=IntKind) :: lig
   integer (kind=IntKind), allocatable :: lmax_kkr(:), ind_lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:), ind_lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:), ind_lmax_rho(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:), ind_lmax_pot(:)
   integer (kind=IntKind), allocatable :: ind_radplane(:)
!
   integer (kind=IntKind) :: NumAlloySubLatts, MaxComponents, rstatus
   integer (kind=IntKind), pointer :: p_AlloySublattIndex(:)
   integer (kind=IntKind), pointer :: p_NumComponents(:)
   integer (kind=IntKind), pointer :: p_AlloyElement(:,:)
   real (kind=RealKind), pointer :: p_AlloyContent(:,:)
   real (kind=RealKind), allocatable :: radplane(:)
!
   real (kind=RealKind) :: volume
   real (kind=RealKind) :: evt(3), em, alpev, cn
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   interface
      subroutine readPositionData(fname,NumAtomsIn,NumAtomsOut)
         use KindParamModule, only : IntKind
         character (len=*), intent(in) :: fname
         integer (kind=IntKind), intent(in), optional :: NumAtomsIn
         integer (kind=IntKind), intent(out), optional :: NumAtomsOut
      end subroutine readPositionData
   end interface
!
   interface
      subroutine copyCharArray2String(a,b,n)
         use KindParamModule, only : IntKind
         integer (kind=IntKind), intent(in) :: n
         character (len=1), intent(in) :: a(:)
         character (len=*), intent(out) :: b
      end subroutine copyCharArray2String
   end interface
!
   alat = 0.0d0
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Text Identification',SystemID)
   rstatus = getKeyValue(tbl_id,'Alloy System Description',SystemTitle)
   rstatus = getKeyValue(tbl_id,'No. Atoms in System (> 0)',NumAtoms)
   rstatus = getKeyValue(tbl_id,'Spin Index Param (>= 1)',nspin)
   rstatus = getKeyValue(tbl_id,'Additional Electrons',AdditionalElectrons)
   rstatus = getKeyValue(tbl_id,'Generate System Movie',printSysMovie)
   if ( rstatus == 1 ) then
      printSysMovie = 0
   endif
   rstatus = getKeyValue(tbl_id,'Uniform Grid Parameters',svalue)
   read(svalue,*) uniform_grid(1:3)
!  -------------------------------------------------------------------
!
   allocate( LmaxKKR(NumAtoms), lmax_kkr(0:NumAtoms), ind_lmax_kkr(NumAtoms) )
   allocate( LmaxPhi(NumAtoms), lmax_phi(0:NumAtoms), ind_lmax_phi(NumAtoms) )
   allocate( LmaxRho(NumAtoms), lmax_rho(0:NumAtoms), ind_lmax_rho(NumAtoms) )
   allocate( LmaxPot(NumAtoms), lmax_pot(0:NumAtoms), ind_lmax_pot(NumAtoms) )
   allocate( RadicalPlaneRatio(NumAtoms), radplane(0:NumAtoms), &
             ind_radplane(NumAtoms) )
   allocate( SiteLIZ(NumAtoms) )
!
   SiteLIZ(1:NumAtoms) = 0
   RadicalPlaneRatio(1:NumAtoms) = ONE
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Current File Path',file_path)
!  -------------------------------------------------------------------
   file_path = adjustl(file_path)
   if ( getKeyValue(tbl_id,'Info Table File Name',info_table) == 0) then
      info_id=getTableIndex(trim(file_path)//adjustl(info_table))
   else 
      info_id = tbl_id
   endif
!  -------------------------------------------------------------------
!
   if ( isKeyExisting(tbl_id,'Atomic Position File Name') ) then
!     ----------------------------------------------------------------
      rstatus = getKeyValue(tbl_id,'Atomic Position File Name',fposi_in)
!     ----------------------------------------------------------------
   else if ( getKeyValue(info_id,'Atomic Position File Name',fposi_in) /= 0) then
      fposi_in = 'position.dat'
   else
!     ----------------------------------------------------------------
!     call ErrorHandler('initSystem','Atom Position Data file is not specified')
!     ----------------------------------------------------------------
      fposi_in = 'position.dat'
   endif
!
   fposi_out = trim(fposi_in)//'_out'
!
   if ( len_trim(fposi_in) > 0 .and. .not.nocaseCompare(fposi_in,'None') ) then
!     ----------------------------------------------------------------
      call readPositionData(trim(file_path)//trim(fposi_in),NumAtomsIn=NumAtoms)
!     ----------------------------------------------------------------
   endif
!
   if ( nspin > 2 ) then
      if ( isKeyExisting(tbl_id,'Moment Direction File Name') ) then
!        -------------------------------------------------------------
         rstatus = getKeyValue(tbl_id,'Moment Direction File Name',fevec_in)
!        -------------------------------------------------------------
      else if ( getKeyValue(info_id,'Moment Direction File Name',fevec_in) /= 0) then
            fevec_in ='None'
      else
         fevec_in ='None'
      endif
      if ( len_trim(fevec_in) > 0 .and. .not.nocaseCompare(fevec_in,'None') ) then
!        -------------------------------------------------------------
         call readMomentDirectionData(trim(file_path)//trim(fevec_in), &
                                      NumAtoms)
!        -------------------------------------------------------------
         fevec_out = trim(fevec_in)//'_out'
      else
         fevec_out = 'Evec_out.dat'
      endif
      allocate(moment(NumAtoms))
   endif
!
   if ( isKeyExisting(tbl_id,'Additional Electrons') ) then
!     ----------------------------------------------------------------
      rstatus = getKeyValue(tbl_id,'Additional Electrons',AdditionalElectrons)
!     ----------------------------------------------------------------
   else
      AdditionalElectrons = ZERO
   endif
!
   allocate( AtomName(NumAtoms) )
   if (isDataStorageExisting('Atomic Number')) then
      AtomicNum => getDataStorage('Atomic Number',NumAtoms,IntegerMark)
      AtomNameCharacter => getDataStorage('Atomic Name',MaxLenOfAtomName*NumAtoms,CharacterMark)
      do i=1,NumAtoms
!        AtomName(i)=AtomNameCharacter((i-1)*MaxLenOfAtomName+1:i*MaxLenOfAtomName)
!        -------------------------------------------------------------
         call copyCharArray2String(AtomNameCharacter((i-1)*MaxLenOfAtomName+1:i*MaxLenOfAtomName), &
                                   AtomName(i),MaxLenOfAtomName)
!        -------------------------------------------------------------
      enddo
   else
!     ----------------------------------------------------------------
      call createDataStorage('Atomic Number',NumAtoms,IntegerType)
!     ----------------------------------------------------------------
      AtomicNum => getDataStorage('Atomic Number',NumAtoms,IntegerMark)
!     ----------------------------------------------------------------
      rstatus = getKeyValue(info_id,'Atom Name',AtomName,NumAtoms)
!     ----------------------------------------------------------------
      do i=1,NumAtoms
         if (len_trim(AtomName(i)) < 1) then
            call ErrorHandler('initSystem','Atom name is not found')
         endif
         AtomicNum(i)=getZtot(AtomName(i))
      enddo
   endif
!
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      Bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call createDataStorage('Bravais Vector',9,RealType)
!     ----------------------------------------------------------------
      Bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
      rstatus = getKeyValue(info_id,'Bravais Vector',3,Bravais,3)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(info_id,'Default Lmax-T matrix',lmax_kkr(0))
   rstatus = getKeyIndexValue(info_id,'Lmax-T matrix',ind_lmax_kkr,lmax_kkr(1:NumAtoms),NumAtoms)
!
!  rstatus = getKeyValue(info_id,'Default Lmax-Wave Func',lmax_phi(0))
   if (getKeyValue(info_id,'Default Lmax-Wave Func',lmax_phi(0),default_param=.false.) /= 0) then
      lmax_phi(0) = lmax_kkr(0)
   endif
   rstatus = getKeyIndexValue(info_id,'Lmax-Wave Func',ind_lmax_phi,lmax_phi(1:NumAtoms),NumAtoms)
!
!  rstatus = getKeyValue(info_id,'Default Lmax-Potential',lmax_pot(0))
   if (getKeyValue(info_id,'Default Lmax-Potential',lmax_pot(0),default_param=.false.) /= 0) then
      rstatus = getKeyValue(tbl_id,'Potential Type (>= 0)',pot_type)
      if (pot_type == ASA .or. pot_type == MuffinTin .or.  pot_type == MuffinTinASA) then
         lmax_pot(0) = 0
      else
         lmax_pot(0) = 2*lmax_kkr(0)
      endif
   endif
   rstatus = getKeyIndexValue(info_id,'Lmax-Potential',ind_lmax_pot,lmax_pot(1:NumAtoms),NumAtoms)
!
!  rstatus = getKeyValue(info_id,'Default Lmax-Charge Den',lmax_rho(0))
   if (getKeyValue(info_id,'Default Lmax-Charge Den',lmax_rho(0),default_param=.false.) /= 0) then
      lmax_rho(0) = lmax_pot(0)
   endif
   rstatus = getKeyIndexValue(info_id,'Lmax-Charge Den',ind_lmax_rho,lmax_rho(1:NumAtoms),NumAtoms)
!
   rstatus = getKeyValue(info_id,'Default Radical Plane Ratio',radplane(0))
   rstatus = getKeyIndexValue(info_id,'Radical Plane Ratio',ind_radplane,radplane(1:NumAtoms),NumAtoms)
!  -------------------------------------------------------------------
   reset_lmax = 0
   do i = 1, NumAtoms
      LmaxKKR(i) = lmax_kkr(ind_lmax_kkr(i))
      LmaxPhi(i) = lmax_phi(ind_lmax_phi(i))
      if (LmaxPhi(i) < LmaxKKR(i) .or. LmaxPhi(i) < 0) then
         LmaxPhi(i) =  LmaxKKR(i)
         reset_lmax = 1
      endif
      LmaxRho(i) = lmax_rho(ind_lmax_rho(i))
      if (LmaxRho(i) > 2*LmaxPhi(i) .or. LmaxRho(i) < 0) then
         LmaxRho(i) = 2*LmaxPhi(i)
         reset_lmax = 1
      endif
      LmaxPot(i) = lmax_pot(ind_lmax_pot(i))
      RadicalPlaneRatio(i) = radplane(ind_radplane(i))
   enddo
   if (MyPE == 0 .and. reset_lmax > 0) then
      call WarningHandler('initSystem','LmaxPhi and LmaxRho are set to new values')
   endif
   deallocate( lmax_kkr, ind_lmax_kkr )
   deallocate( lmax_phi, ind_lmax_phi )
   deallocate( lmax_rho, ind_lmax_rho )
   deallocate( lmax_pot, ind_lmax_pot )
   deallocate( radplane, ind_radplane )
!
   NumVacancies = 0
   NumCoherentPots = 0
   do i=1,NumAtoms
      if (AtomicNum(i) == 0) then
         NumVacancies = NumVacancies + 1
      else if (AtomicNum(i) == -1) then
         NumCoherentPots = NumCoherentPots + 1
      endif
   enddo
!
   allocate(AtomTypeName(NumAtoms),AtomType(NumAtoms))
   NumAtomTypes = 1
   AtomTypeName(1) = AtomName(1)
   AtomType(1) = NumAtomTypes
   LOOP_i: do i=2,NumAtoms
      do ig=1,NumAtomTypes
         if (nocaseCompare(AtomName(i),AtomTypeName(ig)) .and. .not.nocaseCompare(AtomName(i),'CPA')) then
            AtomType(i) = ig
            cycle LOOP_i
         endif
      enddo
      NumAtomTypes = NumAtomTypes + 1
      AtomTypeName(NumAtomTypes) = AtomName(i)
      AtomType(i) = NumAtomTypes
   enddo LOOP_i
   allocate(NumAtomsOfType(NumAtomTypes))
   do ig=1,NumAtomTypes
      NumAtomsOfType(ig) = 0
   enddo
   do i=1,NumAtoms
      ig = AtomType(i)
      NumAtomsOfType(ig) = NumAtomsOfType(ig) + 1
   enddo
!
   allocate(NumAlloyElements(NumAtoms))
   do ig=1,NumAtoms
      NumAlloyElements(ig) = 1
   enddo
   MaxAlloyElements = 1
   if (isDataStorageExisting('Alloy Sublattice Index')) then
      NumAlloySubLatts = getDataStorageSize('Alloy Sublattice Index')
      p_AlloySublattIndex => getDataStorage('Alloy Sublattice Index',NumAlloySubLatts,IntegerMark)
      p_NumComponents => getDataStorage('Number of Components on Sublattice',NumAlloySubLatts,IntegerMark)
      MaxComponents = getDataStorageLDA('Alloy Element')
      p_AlloyElement => getDataStorage('Alloy Element',MaxComponents,NumAlloySubLatts,IntegerMark)
      p_AlloyContent => getDataStorage('Alloy Content',MaxComponents,NumAlloySubLatts,RealMark)
      do i = 1, NumAlloySubLatts
         ig = p_AlloySublattIndex(i)
         NumAlloyElements(ig) = p_NumComponents(i)
         MaxAlloyElements = max(MaxAlloyElements,NumAlloyElements(ig))
      enddo
   else
      NumAlloySubLatts = 0
   endif
!
   allocate(table_line(NumAtoms))
   table_line = 0
   table_size = 0
   do ig = 1, NumAtoms
      table_line(ig) = table_size
      table_size = table_size + NumAlloyElements(ig)
   enddo
   allocate(AlloyElementName(table_size))
   allocate(AlloyElementAN(table_size))
   allocate(AlloyElementContent(table_size))
   if (MaxAlloyElements > 1) then
      do i = 1, NumAlloySubLatts
         ig = p_AlloySublattIndex(i)
         do j = 1, NumAlloyElements(ig)
            lig = table_line(ig) + j
            AlloyElementName(lig) = getName(p_AlloyElement(j,i))
            AlloyElementAN(lig) = p_AlloyElement(j,i)
            AlloyElementContent(lig) = p_AlloyContent(j,i)
         enddo
      enddo
   else
      do ig = 1, NumAtoms
         AlloyElementName(ig) = AtomName(ig)
         AlloyElementAN(ig) = AtomicNum(ig)
         AlloyElementContent(ig) = ONE
      enddo
   endif
!
   allocate( value(NumAtoms) )
!
   if (isDataStorageExisting('Atomic Position')) then
      scaling = getDataStorage('Position Scaling Factor',RealMark)
      AtomPosition => getDataStorage('Atomic Position',3,NumAtoms,RealMark)
   else
      scaling = ONE
!     ----------------------------------------------------------------
      call createDataStorage('Atomic Position',3*NumAtoms,RealType)
!     ----------------------------------------------------------------
      AtomPosition => getDataStorage('Atomic Position',3,NumAtoms,RealMark)
      if ( isKeyExisting(info_id,'Atomic Position') ) then
!        -------------------------------------------------------------
         rstatus = getKeyValue(info_id,'Atomic Position',value,NumAtoms)
!        -------------------------------------------------------------
         do i=1,NumAtoms
            read(value(i),*) AtomPosition(1:3,i)
         enddo
      else
!        -------------------------------------------------------------
         call ErrorHandler('initSystem',                              &
                           'Not able to identify Atom Position Data')
!        -------------------------------------------------------------
      endif
   endif
!
   if (isDataStorageExisting('Moment Direction')) then
      Evec => getDataStorage('Moment Direction',3,NumAtoms,RealMark)
   else if (nspin > 2) then
!     ----------------------------------------------------------------
      call createDataStorage('Moment Direction',3*NumAtoms,RealType)
!     ----------------------------------------------------------------
      Evec => getDataStorage('Moment Direction',3,NumAtoms,RealMark)
      if ( isKeyExisting(info_id,'Moment Direction') ) then
!        -------------------------------------------------------------
         rstatus = getKeyValue(info_id,'Moment Direction',value,NumAtoms)
!        -------------------------------------------------------------
         do i=1,NumAtoms
            read(value(i),*) Evec(1:3,i)
         enddo
      else if ( isKeyExisting(info_id,'Default Moment Direction') ) then
!        -------------------------------------------------------------
         rstatus = getKeyValue(info_id,'Default Moment Direction',3,evt(1:3))
!        -------------------------------------------------------------
         do i=1,NumAtoms
            Evec(1:3,i) = evt(1:3)
         enddo
      else
!        -------------------------------------------------------------
         call ErrorHandler('initSystem',                              &
                           'Not able to identify Moment Direction Data')
!        -------------------------------------------------------------
      endif
   endif
!
!  ===================================================================
!  Normalize Evec
!  ===================================================================
   if (nspin > 2) then
      do i=1,NumAtoms
         em = sqrt(Evec(1,i)*Evec(1,i)+Evec(2,i)*Evec(2,i)+Evec(3,i)*Evec(3,i))
         if (abs(em-ONE) > TEN2m8) then
            Evec(1,i)=Evec(1,i)/em
            Evec(2,i)=Evec(2,i)/em
            Evec(3,i)=Evec(3,i)/em
         else if (em < TEN2m8) then
            Evec(1,i)=ZERO
            Evec(2,i)=ZERO
            Evec(3,i)=ONE
         endif
      enddo
   endif
!
   if (isDataStorageExisting('Old Moment Direction')) then
      Evec_Old => getDataStorage('Old Moment Direction',3,NumAtoms,RealMark)
   else if (nspin > 2) then
!     ----------------------------------------------------------------
      call createDataStorage('Old Moment Direction',3*NumAtoms,RealType)
!     ----------------------------------------------------------------
      Evec_Old => getDataStorage('Old Moment Direction',3,NumAtoms,RealMark)
      do i=1,NumAtoms
         Evec_Old(1:3,i) = Evec(1:3,i)
      enddo
   endif
!
   if (isDataStorageExisting('Constrain Field')) then
      ConstrainField => getDataStorage('Constrain Field',3,NumAtoms,RealMark)
   else if (nspin > 2) then
!     ----------------------------------------------------------------
      call createDataStorage('Constrain Field',3*NumAtoms,RealType)
!     ----------------------------------------------------------------
      ConstrainField => getDataStorage('Constrain Field',3,NumAtoms,RealMark)
      if ( isKeyExisting(info_id,'Constrain Field') ) then
!        -------------------------------------------------------------
         rstatus = getKeyValue(info_id,'Constrain Field',value,NumAtoms)
!        -------------------------------------------------------------
         do i=1,NumAtoms
            read(value(i),*) ConstrainField(1:3,i)
         enddo
      else if ( isKeyExisting(info_id,'Default Constrain Field') ) then
!        -------------------------------------------------------------
         rstatus = getKeyValue(info_id,'Default Constrain Field',3,evt(1:3))
!        -------------------------------------------------------------
         do i=1,NumAtoms
            ConstrainField(1:3,i) = evt(1:3)
         enddo
      else
!        -------------------------------------------------------------
         call ErrorHandler('initSystem',                              &
                           'Not able to identify Constrain Field Data')
!        -------------------------------------------------------------
      endif
   endif
!
   if (isDataStorageExisting('Evec Mixing Parameters')) then
      emix => getDataStorage('Evec Mixing Parameters',NumAtoms,RealMark)
   else if (nspin > 2) then
!     ----------------------------------------------------------------
      call createDataStorage('Evec Mixing Parameters',NumAtoms,RealType)
!     ----------------------------------------------------------------
      emix => getDataStorage('Evec Mixing Parameters',NumAtoms,RealMark)
      if ( isKeyExisting(info_id,'Evec Mix Param.') ) then
!        -------------------------------------------------------------
         rstatus = getKeyValue(info_id,'Evec Mix Param.',value,NumAtoms)
!        -------------------------------------------------------------
         do i=1,NumAtoms
            read(value(i),*)emix(i)
         enddo
      else
         emix(1:NumAtoms) = -ONE
      endif
   endif
!
   if (nspin > 2) then
!     if ( isKeyExisting(info_id,'Default Evec Mix Param.') ) then
      if ( getKeyValue(info_id,'Default Evec Mix Param.',alpev) == 0 ) then
!        -------------------------------------------------------------
!        rstatus = getKeyValue(info_id,'Default Evec Mix Param.',alpev)
!        -------------------------------------------------------------
         do i=1,NumAtoms
            if (emix(i) < ZERO) then
               emix(i) = alpev
            endif
         enddo
      else
!        -------------------------------------------------------------
         call ErrorHandler('initSystem',                              &
                           'Not able to identify Evec Mixing Parameters')
!        -------------------------------------------------------------
      endif
   endif
!
   if (isDataStorageExisting('Force Data')) then
      Force => getDataStorage('Force Data',3,NumAtoms,RealMark)
   else
!     ----------------------------------------------------------------
      call createDataStorage('Force Data',3*NumAtoms,RealType)
!     ----------------------------------------------------------------
      Force => getDataStorage('Force Data',3,NumAtoms,RealMark)
      Force = ZERO
   endif
!
   if (isDataStorageExisting('Energy Pressure')) then
      EnPres => getDataStorage('Energy Pressure',2,NumAtoms,RealMark)
   else
!     ----------------------------------------------------------------
      call createDataStorage('Energy Pressure',2*NumAtoms,RealType)
!     ----------------------------------------------------------------
      EnPres => getDataStorage('Energy Pressure',2,NumAtoms,RealMark)
   endif
!
   if (isDataStorageExisting('RMS Info')) then
      if (nspin <= 2) then
         RMS => getDataStorage('RMS Info',2,NumAtoms,RealMark)
      else
         RMS => getDataStorage('RMS Info',4,NumAtoms,RealMark)
      endif
   else if (nspin <= 2) then
!     ----------------------------------------------------------------
      call createDataStorage('RMS Info',2*NumAtoms,RealType)
!     ----------------------------------------------------------------
      RMS => getDataStorage('RMS Info',2,NumAtoms,RealMark)
   else if (nspin > 2) then
!     ----------------------------------------------------------------
      call createDataStorage('RMS Info',4*NumAtoms,RealType)
!     ----------------------------------------------------------------
      RMS => getDataStorage('RMS Info',4,NumAtoms,RealMark)
   endif
!
   deallocate(value)
!
   volume = ( Bravais(2,1)*Bravais(3,2)-Bravais(3,1)*Bravais(2,2)) &
             *Bravais(1,3)+                                        &
            ( Bravais(3,1)*Bravais(1,2)-Bravais(1,1)*Bravais(3,2)) &
             *Bravais(2,3)+                                        &
            ( Bravais(1,1)*Bravais(2,2)-Bravais(2,1)*Bravais(1,2)) &
             *Bravais(3,3)
   volume=abs(volume)
   alat = volume**THIRD
!
   LmaxPotMax = maxval(LmaxPot(1:NumAtoms))
   LmaxKKRMax = maxval(LmaxKKR(1:NumAtoms))
   LmaxRhoMax = maxval(LmaxRho(1:NumAtoms))
   LmaxPhiMax = maxval(LmaxPhi(1:NumAtoms))
   LmaxMax = max(LmaxPotMax, LmaxKKRMax, LmaxRhoMax, LmaxPhiMax)
!
   end subroutine initSystem
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSystem()
!  ===================================================================
   implicit none
!
   nullify( AtomicNum )
   deallocate( AtomName, AtomTypeName, AtomType, NumAtomsOfType )
   deallocate( NumAlloyElements, SiteLIZ )
   nullify( AtomPosition )
   nullify( Evec, emix )
   nullify( ConstrainField )
   nullify( Force, EnPres )
!
   deallocate(table_line)
   deallocate(AlloyElementName, AlloyElementAN, AlloyElementContent)
   MaxAlloyElements = 0
!
   deallocate(LmaxKKR, LmaxRho, LmaxPot)
   deallocate(RadicalPlaneRatio)
   scaling = ONE
!
   if (nspin>2) then
      deallocate(moment)
   endif
!
   end subroutine endSystem
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSystemID() result (sid)
!  ===================================================================
   implicit none
   character (len=50) :: sid
!
   sid = SystemID
!
   end function getSystemID
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSystemTitle() result (st)
!  ===================================================================
   implicit none
   character (len=70) :: st
!
   st = SystemTitle
!
   end function getSystemTitle
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumAtoms() result (na)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: na
!
   na = NumAtoms
!
   end function getNumAtoms
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSiteLIZ(i,n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i,n
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setSiteLIZ','Invalid atom index',i)
   endif
   SiteLIZ(i) = n
!
   end subroutine setSiteLIZ
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetSiteLIZ()
!  ===================================================================
   implicit none
!
   SiteLIZ(1:NumAtoms) = 0
!
   end subroutine resetSiteLIZ
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSiteLIZ_one(i)                                 result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: n
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getSiteLIZ','Invalid atom index',i)
   endif
   n = SiteLIZ(i)
!
   end function getSiteLIZ_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSiteLIZ_all()                              result(pLIZ)
!  ===================================================================
   implicit none
   integer (kind=IntKind), pointer :: pLIZ(:)
!
   pLIZ => SiteLIZ
!
   end function getSiteLIZ_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumAtomTypes() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   n = NumAtomTypes
!
   end function getNumAtomTypes
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomType_0(i) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
!
   integer (kind=IntKind) :: n
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getAtomType','Invalid atom index',i)
   endif
   n = AtomType(i)
!
   end function getAtomType_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomType_1() result(p_ntype)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: p_ntype(:)
!
   p_ntype => AtomType
!
   end function getAtomType_1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomTypeName(i) result(a)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
!
   character (len=MaxLenOfAtomName) :: a
!
   if (i<1 .or. i>NumAtomTypes) then
      call ErrorHandler('getAtomTypeName','Invalid type index',i)
   endif
   a = AtomTypeName(i)
!
   end function getAtomTypeName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumAtomsOfType(i) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
!
   integer (kind=IntKind) :: n
!
   if (i<1 .or. i>NumAtomTypes) then
      call ErrorHandler('getNumAtomsOfType','Invalid type index',i)
   endif
   n = NumAtomsOfType(i)
!
   end function getNumAtomsOfType
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumAlloyElements(i) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
!
   integer (kind=IntKind) :: n
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getNumAlloyElements','Invalid global atom index',i)
   endif
   n = NumAlloyElements(i)
!
   end function getNumAlloyElements
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAlloyTableSize() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   n = table_size
!
   end function getAlloyTableSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAlloyElementName(ig,ia) result(a)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ia,ig
   integer (kind=IntKind) :: lig
!
   character (len=MaxLenOfAtomName) :: a
!
   if (ig < 1 .or. ig > NumAtoms) then
      call ErrorHandler('getAlloyElementName','Invalid global atom index',ig)
   else if (ia < 1 .or. ia > NumAlloyElements(ig)) then
      call ErrorHandler('getAlloyElementName','Invalid alloy element index',ia)
   endif
   lig = table_line(ig) + ia
   a = AlloyElementName(lig)
!
   end function getAlloyElementName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAlloyElementContent(ig,ia) result(c)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ia,ig
   integer (kind=IntKind) :: lig
!
   real (kind=RealKind) :: c
!
   if (ig < 1 .or. ig > NumAtoms) then
      call ErrorHandler('getAlloyElementContent','Invalid global atom index',ig)
   else if (ia < 1 .or. ia > NumAlloyElements(ig)) then
      call ErrorHandler('getAlloyElementContent','Invalid alloy element index',ia)
   endif
   lig = table_line(ig) + ia
   c = AlloyElementContent(lig)
!
   end function getAlloyElementContent
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAlloySpeciesIndex(ig,ia) result(lig)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ia,ig
   integer (kind=IntKind) :: lig
!
   lig = table_line(ig) + ia
!
   end function getAlloySpeciesIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumVacancies() result (na)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: na
!
   na = NumVacancies
!
   end function getNumVacancies
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCoherentPots() result (na)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: na
!
   na = NumCoherentPots
!
   end function getNumCoherentPots
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setBravaisLattice(br)
!  ===================================================================
   implicit none
   real (kind=RealKind), intent(in) :: br(3,3)
!
   Bravais(1:3, 1:3) = br(1:3, 1:3)
!
   end subroutine setBravaisLattice
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getBravaisLattice() result (br)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: br(3,3)
!
   br(1:3, 1:3) = Bravais(1:3, 1:3)
!
   end function getBravaisLattice
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSystemTable(anam, anum, apos)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: anam(*)
   integer (kind=IntKind), intent(in) :: anum(*)
   real (kind=RealKind), intent(in) :: apos(3,*)
!
   integer (kind=IntKind) :: i
!
   if (len_trim(anam(1)) > MaxLenOfAtomName) then
      call ErrorHandler('setSystemTable','length of anam > MaxLenOfAtomName', &
                        anam(1))
   endif
!
   AtomName(1:NumAtoms) = anam(1:NumAtoms)
   AtomicNum(1:NumAtoms) = anum(1:NumAtoms)
   do i = 1,NumAtoms
      AtomPosition(1:3,i) = apos(1:3,i)
   enddo
!
   end subroutine setSystemTable
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomName_one(ig,ia) result (nm)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ig
   integer (kind=IntKind), intent(in), optional :: ia
   integer (kind=IntKind) :: lig
   character (len=MaxLenOfAtomName) :: nm
!
   if (ig<1 .or. ig>NumAtoms) then
      call ErrorHandler('getAtomName','Invalid atom/site index',ig)
   endif
   if (present(ia)) then
      if (ia < 1 .or. ia > NumAlloyElements(ig)) then
         call ErrorHandler('getAtomName','Invalid alloy element index at site',ia,ig)
      else
         lig = table_line(ig) + ia
         nm = AlloyElementName(lig)
      endif
   else
      nm = AtomName(ig)
   endif
!
   end function getAtomName_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomName_all() result (nm)
!  ===================================================================
   implicit none
   character (len=MaxLenOfAtomName), pointer :: nm(:)
!
   nm => AtomName
!
   end function getAtomName_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicNumber_one(ig,ia) result (nm)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: ig
   integer (kind=IntKind), intent(in), optional :: ia
   integer (kind=IntKind) :: nm, lig
!
   if (ig<1 .or. ig>NumAtoms) then
      call ErrorHandler('getAtomicNumber','Invalid atom index',ig)
   endif
   if (present(ia)) then
      if (ia < 1 .or. ia > NumAlloyElements(ig)) then
         call ErrorHandler('getAtomicNumber','Invalid alloy element index at site',ia,ig)
      else
         lig = table_line(ig) + ia
         nm = AlloyElementAN(lig)
      endif
   else
      nm = AtomicNum(ig)
   endif
!
   end function getAtomicNumber_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicNumber_all() result (p_nm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), pointer :: p_nm(:)
!
   p_nm => AtomicNum
!
   end function getAtomicNumber_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomPosition_one(i) result (po)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: po(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getAtomPosition','Invalid atom index',i)
   endif
   po(1:3) = AtomPosition(1:3,i)
!
   end function getAtomPosition_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomPosition_all() result (p_po)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: p_po(:,:)
!
   p_po => AtomPosition(1:3,1:NumAtoms)
!
   end function getAtomPosition_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getScalingFactor() result (s)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: s
!
   s = scaling
!
   end function getScalingFactor
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setAtomPosition(i,po)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: po(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setAtomPosition','Invalid atom index',i)
   endif
   AtomPosition(1:3,i) = po(1:3)
!
   end subroutine setAtomPosition
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomentDirectionMixingParam(i) result (alpev)
!  ===================================================================
   use MathParamModule, only : ONE
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: alpev
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getMomentDirectionMixingParam','Invalid atom index',i)
   else if (nspin <= 2) then
      alpev = ZERO
   else
      alpev = emix(i)
   endif

   end function getMomentDirectionMixingParam
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomentDirection_one(i) result (evec_out)
!  ===================================================================
   use MathParamModule, only : ONE
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: evec_out(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getMomentDirection','Invalid atom index',i)
   else if (nspin < 3) then
      evec_out(1) = ZERO
      evec_out(2) = ZERO
      evec_out(3) = ONE
   else
      evec_out(1:3) = Evec(1:3,i)
   endif

   end function getMomentDirection_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomentDirection_all() result (evec_l)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: evec_l(:,:)
!
   evec_l => Evec

   end function getMomentDirection_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setMomentDirection(i,evec_l)
!  ===================================================================
   use MathParamModule, only : TEN2m8, ONE
!
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: evec_l(3)
   real (kind=RealKind) :: em
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setMomentDirection','Invalid atom index',i)
   else if (nspin < 3) then
      return
   endif
!
   em = sqrt(evec_l(1)*evec_l(1)+evec_l(2)*evec_l(2)+evec_l(3)*evec_l(3))
   if (em < TEN2m8) then
      call ErrorHandler('setMomentDirection','evec is zero')
   else if (abs(em - ONE) > TEN2m8) then
      Evec(1:3,i) = evec_l(1:3)/em
   else
      Evec(1:3,i) = evec_l(1:3)
   endif

   end subroutine setMomentDirection
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomentDirectionOld_one(i) result (evec_l)
!  ===================================================================
   use MathParamModule, only : ONE
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: evec_l(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getMomentDirectionOld','Invalid atom index',i)
   else if (nspin < 3) then
      evec_l(1) = ZERO
      evec_l(2) = ZERO
      evec_l(3) = ONE
   else
      evec_l(1:3) = Evec_Old(1:3,i)
   endif

   end function getMomentDirectionOld_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMomentDirectionOld_all() result (evec_l)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: evec_l(:,:)
!
   evec_l => Evec_Old

   end function getMomentDirectionOld_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setMomentDirectionOld(i,evec_in)
!  ===================================================================
   use MathParamModule, only : TEN2m8, ONE
!
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: evec_in(3)
   real (kind=RealKind) :: em
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setMomentDirectionOld','Invalid atom index',i)
   else if (nspin < 3) then
      return
   endif
!
   em = sqrt(evec_in(1)*evec_in(1)+evec_in(2)*evec_in(2)+evec_in(3)*evec_in(3))
   if (em < TEN2m8) then
      call ErrorHandler('setMomentDirectioOldn','evec is zero')
   else if (abs(em - ONE) > TEN2m8) then
      Evec_Old(1:3,i) = evec_in(1:3)/em
   else
      Evec_Old(1:3,i) = evec_in(1:3)
   endif

   end subroutine setMomentDirectionOld
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getConstrainField_one(i) result (b_con)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: b_con(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getConstrainField','Invalid atom index',i)
   else if (nspin <= 2) then
      b_con(1:3) = ZERO
   else
      b_con(1:3) = ConstrainField(1:3,i)
   endif

   end function getConstrainField_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getConstrainField_all() result (b_con)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: b_con(:,:)
!
   b_con => ConstrainField
!
   end function getConstrainField_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setConstrainField(i,b_con)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: b_con(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setConstrainField','Invalid atom index',i)
   else if (nspin <= 2) then
      return
   endif
   ConstrainField(1:3,i) = b_con(1:3)

   end subroutine setConstrainField
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getForce_one(i)                                     result(pf)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: pf(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getForce','Invalid atom index',i)
   endif
   pf(1:3) = Force(1:3,i)
!
   end function getForce_one
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getForce_all()                                 result(pf)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: pf(:,:)
!
   pf => Force
!
   end function getForce_all
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setForce(i,pf)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: pf(3)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setForce','Invalid atom index',i)
   endif
   Force(1:3,i) = pf(1:3)
!
   end subroutine setForce
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setEnPres(i,pen)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: pen(2)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setForce','Invalid atom index',i)
   endif
   EnPres(1:2,i) = pen(1:2)
!
   end subroutine setEnPres
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEnPres()                                    result(pf)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: pf(:,:)
!
   pf => EnPres
!
   end function getEnPres
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setRMSInfo(i,rms_in)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in), target :: rms_in(1:4)
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('setForce','Invalid atom index',i)
   endif
   RMS(1:2,i) = rms_in(1:2)  ! charge, potential
   if ( nspin>2 ) then
      RMS(3:4,i) = rms_in(3:4) ! evec,constrain field
   endif
!
   end subroutine setRMSInfo
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRMSInfo()                                    result(pf)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: pf(:,:)
!
   pf => RMS
!
   end function getRMSInfo
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomEnergy()                                 result(pen)
!  ===================================================================
   implicit none
   real (kind=RealKind), pointer :: pen(:,:)
!
   pen => EnPres
!
   end function getAtomEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setAtomEnergy(en)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: en(:,:)
!
   EnPres = en
!
   end subroutine setAtomEnergy
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLatticeConstant(a0)
!  ===================================================================
   implicit none
   real (kind=RealKind), intent(in) :: a0
!
   alat=a0
   end subroutine setLatticeConstant
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLatticeConstant() result(a0)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: a0
!
   a0=alat
   end function getLatticeConstant
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAdditionalElectrons() result(ae)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: ae
!
   ae=AdditionalElectrons
   end function getAdditionalElectrons
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getUniformGridParam() result(ngu)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: ngu(3)
!
   ngu(1:3) = uniform_grid(1:3)
   end function getUniformGridParam
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSystem(fu_in)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), optional, intent(in) :: fu_in
   integer(kind=IntKind) :: fu
!
   if ( present(fu_in) ) then
      fu = fu_in
   else
      fu = 6
      write(fu,'(/,80(''-''))')
      write(fu,'(/,24x,a)')'*******************************'
      write(fu,'( 24x,a )')'*   Output from printSystem   *'
      write(fu,'(24x,a,/)')'*******************************'
      write(fu,'(/,80(''=''))')
   endif
   write(fu,'(''# Number of Atoms         :'',i7)')NumAtoms
   write(fu,'(''# Bravais Lattice Vector  :'')')
   write(fu,'(''#      '',3f13.8)')Bravais(1:3,1)
   write(fu,'(''#      '',3f13.8)')Bravais(1:3,2)
   write(fu,'(''#      '',3f13.8)')Bravais(1:3,3)
   if ( fu==6 ) then
      write(fu,'(80(''=''))')
   endif
   end subroutine printSystem
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxKKR0() result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: lmax
!
   lmax = LmaxKKRMax
!
   end function getLmaxKKR0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxKKR1(i) result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: lmax
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getLmaxKKR','Invalid atom index',i)
   endif
   lmax = LmaxKKR(i)
!
   end function getLmaxKKR1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxPhi0() result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: lmax
!
   lmax = LmaxPhiMax
!
   end function getLmaxPhi0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxPhi1(i) result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: lmax
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getLmaxPhi','Invalid atom index',i)
   endif
   lmax = LmaxPhi(i)
!
   end function getLmaxPhi1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxRho0() result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: lmax
!
   lmax = LmaxRhoMax
!
   end function getLmaxRho0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxRho1(i) result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: lmax
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getLmaxKKR','Invalid atom index',i)
   endif
   lmax = LmaxRho(i)
!
   end function getLmaxRho1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxPot0() result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: lmax
!
   lmax = LmaxPotMax
!
   end function getLmaxPot0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxPot1(i) result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: lmax
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getLmaxKKR','Invalid atom index',i)
   endif
   lmax = LmaxPot(i)
!
   end function getLmaxPot1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLmaxMax() result (lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: lmax
!
   lmax = LmaxMax
!
   end function getLmaxMax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRadicalPlaneRatio(i) result(radpln)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: radpln
!
   if (i<1 .or. i>NumAtoms) then
      call ErrorHandler('getRadicalPlaneRatio','Invalid atom index',i)
   endif
   radpln = RadicalPlaneRatio(i)
!
   end function getRadicalPlaneRatio
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateSystem(vname_in,getAtom2Proc)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getMyPEinGroup, GlobalSumInGroup
!
   use DataServiceCenterModule, only : isDataStorageExisting,         &
                                       getDataStorage,                &
                                       RealType, RealMark
!
   implicit none
!
   character (len=*), intent(in) :: vname_in
   character (len=80) :: vname
!
   integer (kind=IntKind) :: i, id, MyPEinGroup, n
!
   real (kind=RealKind), pointer :: vec(:,:)
   real (kind=RealKind), pointer :: vecp(:)
!
   interface
      function getAtom2Proc(a) result(p)
         use KindParamModule, only : IntKind
         implicit none
!
         integer (kind=IntKind), intent(in) :: a
         integer (kind=IntKind) :: p
      end function getAtom2Proc
   end interface
!
   id = getGroupID('Unit Cell')
   MyPEinGroup = getMyPEinGroup(id)
!
   vname = vname_in
   vname = adjustl(vname)
!
   if ( trim(vname)=="LIZ Size" ) then
      call GlobalSumInGroup(id,SiteLIZ,NumAtoms)
      return
   else if ( .not.isDataStorageExisting(trim(vname)) ) then
      call ErrorHandler('updateSystem','The data does not exist',trim(vname))
   else if (trim(vname) == 'Atomic Position') then
      NumPositionUpdates = NumPositionUpdates + 1
      n = 3
   else if (trim(vname) == 'Moment Direction') then
      NumEvecUpdates = NumEvecUpdates + 1
      n = 3
   else if (trim(vname) == 'Old Moment Direction') then
      NumEvecOldUpdates = NumEvecOldUpdates + 1
      n = 3
   else if (trim(vname) == 'Constrain Field') then
      NumConstrFieldUpdates = NumConstrFieldUpdates + 1
      n = 3
   else if (trim(vname) == 'Force Data') then
      NumFoceVectUpdates = NumFoceVectUpdates + 1
      n = 3
   else if (trim(vname) == 'Energy Pressure') then
      NumEnergyUpdates = NumEnergyUpdates + 1
      n = 2
   else if (trim(vname) == 'RMS Info') then
      NumEnergyUpdates = NumRMSUpdates + 1
      if (nspin>2) then
         n = 4
      else
         n=2
      endif
   endif
!
   vec => getDataStorage(trim(vname),n,NumAtoms,RealMark)
!   if (NumAtoms <= 42666) then
      do i = 1, NumAtoms
         if (getAtom2Proc(i) /= MyPEinGroup) then
            vec(1:n,i) = ZERO
         endif
      enddo
!   endif
   nullify(vec)
!
   vecp => getDataStorage(trim(vname),NumAtoms*n,RealMark)
!  -------------------------------------------------------------------
   call GlobalSumInGroup(id,vecp,NumAtoms*n)
!  -------------------------------------------------------------------
   nullify(vecp)
!
   end subroutine updateSystem
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeMomentDirectionData()
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
   implicit none
!
   integer (kind=IntKind) :: i
   integer (kind=IntKind), parameter :: fu = 210
!
   if (nspin <= 2) then
      return
   endif
!
   if (MyPE == 0) then
      open(unit=fu,file=trim(file_path)//trim(fevec_out),form='formatted', &
           status='unknown')
      write(fu,'(a)')   '! ****************************************************'
      write(fu,'(a)')   '! '
      write(fu,'(a,i4)')'! This is Moment Direction update:',NumEvecUpdates
      write(fu,'(a,i4)')'! This is Constrain Fields update:',NumConstrFieldUpdates
      write(fu,'(a)')   '! '
      write(fu,'(a)')   '! ****************************************************'
      do i = 1, NumAtoms
         write(fu,'(a10,2x,6f22.13,1x,f10.5)') AtomName(i), Evec(1:3,i),  &
                                               ConstrainField(1:3,i),emix(i)
      enddo
      close(fu)
   endif
!
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
!
   end subroutine writeMomentDirectionData
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeMomentMovie(iscf,isd)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
   implicit none
!
   logical :: file_exist
!
   character (len=160) :: fname
!
   integer (kind=IntKind) :: i
   integer (kind=IntKind), intent(in) :: iscf,isd
   integer (kind=IntKind), parameter :: fu = 209
!
   if (nspin <= 2) then
      return
   endif
!
   if (MyPE == 0) then
      fname = trim(file_path)//trim(fevec_out)//'_movie'
      inquire(file=fname,exist=file_exist)
      if (file_exist) then
         open(unit=fu,file=fname,form='formatted',status='old',position='append')
      else
         open(unit=fu,file=fname,form='formatted',status='unknown')
      endif
      write(fu,'(a)')   '! ****************************************************'
      write(fu,'(a,i4)')'! SCF Iteration: ',iscf
      write(fu,'(a,i4)')'! S-D Time Step: ',isd
      write(fu,'(a)')   '! ****************************************************'
      do i = 1, NumAtoms
         write(fu,'(a10,1x,7f18.14)')AtomName(i), moment(i),Evec(1:3,i), &
                                     ConstrainField(1:3,i)
      enddo
      close(fu)
   endif
!
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
!
   end subroutine writeMomentMovie
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeAtomPositionData()
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
   implicit none
!
   integer (kind=IntKind) :: i
   integer (kind=IntKind), parameter :: fu = 211
!
   if (MyPE == 0) then
      open(unit=fu,file=trim(file_path)//trim(fposi_out),form='formatted', &
           status='unknown')
      write(fu,'(a)')   '! ****************************************************'
      write(fu,'(a)')   '! '
      write(fu,'(a,i4)')'! This is Atomic Positions update:',NumPositionUpdates
      write(fu,'(a)')   '! '
      write(fu,'(a)')   '! ****************************************************'
      write(fu,'(3f22.13)')Bravais(1:3,1)
      write(fu,'(3f22.13)')Bravais(1:3,2)
      write(fu,'(3f22.13)')Bravais(1:3,3)
      write(fu,'(a)')   '! '
      do i = 1, NumAtoms
         write(fu,'(a10,2x,3f22.13)') AtomName(i), AtomPosition(1:3,i)
      enddo
      close(fu)
   endif
!
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
!
   end subroutine writeAtomPositionData
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine writeAtomPositionMovie(iscf)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
   implicit none
!
   logical :: file_exist
!
   character (len=160) :: fname
!
   integer (kind=IntKind), intent(in) :: iscf
   integer (kind=IntKind) :: i
   integer (kind=IntKind), parameter :: fu = 212
!
   if (MyPE == 0) then
      fname = trim(file_path)//trim(fposi_out)//'_movie'
      inquire(file=fname,exist=file_exist)
      if (file_exist) then
         open(unit=fu,file=fname,form='formatted',status='old',position='append')
      else
         open(unit=fu,file=fname,form='formatted',status='unknown')
      endif
      write(fu,'(a)')   '! ****************************************************'
      write(fu,'(a,i4)')'! Iteration: ',iscf
      write(fu,'(a)')   '! ****************************************************'
      write(fu,'(3f22.13)')Bravais(1:3,1)
      write(fu,'(3f22.13)')Bravais(1:3,2)
      write(fu,'(3f22.13)')Bravais(1:3,3)
      do i = 1, NumAtoms
         write(fu,'(a10,2x,3f22.13)')AtomName(i),                      &
                    AtomPosition(1:3,i)
      enddo
      close(fu)
   endif
!
!  -------------------------------------------------------------------
   call SyncAllPEs()
!  -------------------------------------------------------------------
!
   end subroutine writeAtomPositionMovie
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine updateSystemMovie(moment_in)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: moment_in(1:NumAtoms)
!
   integer (kind=IntKind) :: GroupID
!
   if (nspin <= 2) then
      return
   endif
!
   moment(1:NumAtoms) = moment_in(1:NumAtoms)
!
   end subroutine updateSystemMovie
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine resetSystemMovie()
!  ===================================================================
   implicit none
!
   if (nspin <= 2) then
      return
   endif
!
   Force(1:3,1:NumAtoms)  = ZERO
   EnPres(1:2,1:NumAtoms) = ZERO
   RMS(1:2,1:NumAtoms)    = ZERO
   if ( nspin > 2 ) then
      RMS(3:4,1:NumAtoms) = ZERO
      Evec(1:3,1:NumAtoms) = ZERO
      ConstrainField(1:3,1:NumAtoms) = ZERO
   endif
!
   end subroutine resetSystemMovie
!  ===================================================================
end module SystemModule
