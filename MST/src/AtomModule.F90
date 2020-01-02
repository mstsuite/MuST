module AtomModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : LizLmaxStruct
   use PublicParamDefinitionsModule, only : MaxLenFileName
   use ChemElementModule, only : MaxLenOfAtomName
!
public :: initAtom,    &
          endAtom,     &
          setLocalMoment, &
          getLocalMoment, &
          setLocalEvecOld,   &
          getLocalEvecOld,   &
          setLocalEvecNew,   &
          getLocalEvecNew,   &
          setLocalEvecOut,   &
          getLocalEvecOut,   &
          setLocalConstrainField, &
          getLocalConstrainField, &
          getLocalAtomicNumber, &
          getLocalAtomName,  &
          getLocalAtomNickName,  &
          getLocalAtomPosition, &
          getLocalNumSpecies,  &
          getLocalSpeciesContent,  &
          getLizLmax, &
          getPotLmax, &
          setPotLmax, &
          getTruncPotLmax, &
          setTruncPotLmax, &
          getRhoLmax, &
          setRhoLmax, &
          getPhiLmax, &
          setPhiLmax, &
          getKKRLmax, &
          getMaxLmax, &
          getStepFuncLmax, &
          getInPotFileName,  &
          getInPotFileForm,  &
          getOutPotFileName,  &
          getOutPotFileForm,  &
          getInValDenFileName,  &
          getInValDenFileForm,  &
          getOutValDenFileName,  &
          getOutValDenFileForm,  &
          getGridData, &
          getMixingParam4Rho,    &
          getMixingParam4Pot,    &
          getMixingParam4Mom,    &
          getMixingParam4Charge, &
          getMixingParam4Evec,   &
          getScreenPotential,    &
          getScreenRcut,         &
          getRcutPseudo,         &
          getMaxRad,             &
          getAtomCoreRad,        &
          setAtomCoreRad,        &
          getAtomMuffinTinRad,   &
          setAtomMuffinTinRad,   &
          printAtom,             &
          printAtomMomentInfo
!
   interface setLocalMoment
      module procedure setLocalMoment_s, setLocalMoment_v
   end interface

   interface getInPotFileName
      module procedure getIPFN0, getIPFN1, getIPFN2
   end interface
!
   interface getInPotFileForm
      module procedure getIPFF0, getIPFF1
   end interface
!
   interface getOutPotFileName
      module procedure getOPFN0, getOPFN1, getOPFN2
   end interface
!
   interface getOutPotFileForm
      module procedure getOPFF0, getOPFF1
   end interface
!
   interface getInValDenFileName
      module procedure getIVDFN0, getIVDFN1
   end interface
!
   interface getOutValDenFileName
      module procedure getOVDFN0, getOVDFN1
   end interface
!
private
   type LmaxStruct
      integer (kind=IntKind) :: lmax_kkr
      integer (kind=IntKind) :: lmax_step
      integer (kind=IntKind) :: lmax_phi
      integer (kind=IntKind) :: lmax_rho
      integer (kind=IntKind) :: lmax_pot
      integer (kind=IntKind) :: lmax_pot_trunc
   end type LmaxStruct
!
   type MixParamStruct
      real (kind=RealKind) :: alpha_rho
      real (kind=RealKind) :: alpha_pot
      real (kind=RealKind) :: alpha_chg
      real (kind=RealKind) :: alpha_mom
      real (kind=RealKind) :: alpha_evec
   end type MixParamStruct
!
   logical :: Initialized = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: print_level
!
!  ===================================================================
!  LizLmax Data.
!  -------------------------------------------------------------------
   type (LizLmaxStruct), allocatable, target :: LizLmax(:)
!  ===================================================================
!
!  ===================================================================
!  Lmax, Jmax, Kmax data
!  -------------------------------------------------------------------
   type (LmaxStruct), allocatable, target :: Lmax(:)
   integer (kind=IntKind) :: lmax_max
   real (kind=RealKind) :: max_rad
!  ===================================================================
!
!  ===================================================================
!  Mixing Parameter Data
!  -------------------------------------------------------------------
   type (MixParamStruct), allocatable, target :: MixParam(:)
!  ===================================================================
!
   integer (kind=IntKind), parameter :: MessageID = 50000
   real (kind=RealKind), parameter   :: Rcut_Pseudo_def = 1.7d0
!
!  ===================================================================
!  Atomic property data.
!  -------------------------------------------------------------------
   type AtomPropertyStruct
      character (len=2), allocatable :: AtomName(:)
      character (len=MaxLenOfAtomName), allocatable :: NickName(:)
      character (len=MaxLenFileName), allocatable :: potinfile(:)
      character (len=11) :: potinform
      character (len=MaxLenFileName), allocatable :: potoutfile(:)
      character (len=11) :: potoutform
      integer (kind=IntKind), allocatable :: AtomicNumber(:)
      integer (kind=IntKind) :: GlobalIndex
      integer (kind=IntKind) :: NumSpecies
      logical :: Rmt_fixed
      real (kind=RealKind), allocatable :: AtomContent(:)
      real (kind=RealKind) :: Rmt
      real (kind=RealKind) :: Rcore
      real (kind=RealKind) :: potScreen
      real (kind=RealKind) :: Rcut_pseudo
      real (kind=RealKind) :: Position(3)
      real (kind=RealKind) :: VolumeMT
      real (kind=RealKind) :: VolumeINSC
      real (kind=RealKind) :: VolumeVP
      real (kind=RealKind) :: VolumeWS
   end type AtomPropertyStruct
   type (AtomPropertyStruct), allocatable :: AtomProperty(:)
!
   character (len=MaxLenFileName) :: InPotFileName
   character (len=MaxLenFileName) :: OutPotFileName
   character (len=11) :: InPotFileForm
   character (len=11) :: OutPotFileForm
!  ===================================================================
!
!  ===================================================================
!  Local Moment Data
!  -------------------------------------------------------------------
   type LocalMomentStruct
      real (kind=RealKind) :: moment
      real (kind=RealKind) :: evec_old(3)
      real (kind=RealKind) :: evec_new(3)
      real (kind=RealKind) :: evec_out(3)
      real (kind=RealKind) :: b_con(3)
   end type LocalMomentStruct
   type (LocalMomentStruct), allocatable :: LocalMoment(:)
!  ===================================================================
!
!  ===================================================================
!  Local Grid Data
!  -------------------------------------------------------------------
   type LocalGridStruct
      integer (kind=IntKind) :: ndivin
      integer (kind=IntKind) :: ndivout
      integer (kind=IntKind) :: nmult
   end type LocalGridStruct
   type (LocalGridStruct), allocatable :: GridData(:)
!  ===================================================================
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initAtom(info_id,istop,iprint)
!  ===================================================================
   use MathParamModule, only : ZERO, ONE, TEN2m3
   use ChemElementModule, only : getZtot
   use StringModule, only : initString, endString, getNumTokens, readToken
   use MPPModule, only : MyPE
   use PublicParamDefinitionsModule, only : ASA, MuffinTin, MuffinTinASA
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
   use InputModule, only : getKeyValue, getKeyIndexValue
   use ScfDataModule, only : inputpath, isKKRCPA,  getPotentialTypeParam
   use SystemModule, only : getNumAtoms, getNumAlloyElements, getAlloyElementContent
   use SystemModule, only : getAlloyElementName
   use SystemModule, only : getAtomPosition
   use SystemModule, only : getMomentDirection, getConstrainField
   use SystemModule, only : getMomentDirectionMixingParam
!
   implicit none
!
   character (len=*) :: istop
   character (len=1) :: dummy
   character (len=2) :: s2
   character (len=50) :: fname
   character (len=160) :: path_fname
   character (len=150), allocatable :: lmax_shell(:)
   character (len=150), allocatable :: potinname(:)
   character (len=50), allocatable :: potoutname(:)
!
   integer (kind=IntKind), intent(in) :: info_id, iprint
   integer (kind=IntKind) :: i, j, ig, n, nt, GlobalNumAtoms
   integer (kind=IntKind) :: anum, funit, flen, lmax_n, ia, rstatus
   integer (kind=IntKind), allocatable :: potinform(:)
   integer (kind=IntKind), allocatable :: potoutform(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_step(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:)
   integer (kind=IntKind), allocatable :: lmax_pot_trunc(:)
   integer (kind=IntKind), allocatable :: lmax_rho(:)
   integer (kind=IntKind), allocatable :: num_shells(:)
   integer (kind=IntKind), allocatable :: nmax_liz(:)
   integer (kind=IntKind), allocatable :: ndivin(:)
   integer (kind=IntKind), allocatable :: ndivout(:)
   integer (kind=IntKind), allocatable :: nmult(:)
!
   integer (kind=IntKind), allocatable :: ind_lmax_shell(:)
   integer (kind=IntKind), allocatable :: ind_potinname(:)
   integer (kind=IntKind), allocatable :: ind_potinform(:)
   integer (kind=IntKind), allocatable :: ind_potoutname(:)
   integer (kind=IntKind), allocatable :: ind_potoutform(:)
   integer (kind=IntKind), allocatable :: ind_lmax_kkr(:)
   integer (kind=IntKind), allocatable :: ind_lmax_step(:)
   integer (kind=IntKind), allocatable :: ind_lmax_phi(:)
   integer (kind=IntKind), allocatable :: ind_lmax_pot(:)
   integer (kind=IntKind), allocatable :: ind_lmax_pot_trunc(:)
   integer (kind=IntKind), allocatable :: ind_lmax_rho(:)
   integer (kind=IntKind), allocatable :: ind_num_shells(:)
   integer (kind=IntKind), allocatable :: ind_nmax_liz(:)
   integer (kind=IntKind), allocatable :: ind_alpha_rho(:)
   integer (kind=IntKind), allocatable :: ind_alpha_chg(:)
   integer (kind=IntKind), allocatable :: ind_alpha_pot(:)
   integer (kind=IntKind), allocatable :: ind_alpha_mom(:)
   integer (kind=IntKind), allocatable :: ind_ndivin(:)
   integer (kind=IntKind), allocatable :: ind_ndivout(:)
   integer (kind=IntKind), allocatable :: ind_nmult(:)
   integer (kind=IntKind), allocatable :: ind_cutoff_r(:)
   integer (kind=IntKind), allocatable :: ind_cutoff_r_s(:)
   integer (kind=IntKind), allocatable :: ind_pseudo_r(:)
   integer (kind=IntKind), allocatable :: ind_potscreen(:)
   integer (kind=IntKind), allocatable :: ind_rmt_desire(:)
   integer (kind=IntKind), allocatable :: ind_rcr_desire(:)
!
   real (kind=RealKind), allocatable :: alpha_rho(:)
   real (kind=RealKind), allocatable :: alpha_pot(:)
   real (kind=RealKind), allocatable :: alpha_mom(:)
   real (kind=RealKind), allocatable :: alpha_chg(:)
   real (kind=RealKind), allocatable :: cutoff_r(:)
   real (kind=RealKind), allocatable :: pseudo_r(:)
   real (kind=RealKind), allocatable :: potScreen(:)
   real (kind=RealKind), allocatable :: cutoff_r_s(:)
   real (kind=RealKind), allocatable :: rmt_desire(:)
   real (kind=RealKind), allocatable :: rcr_desire(:)
   real (kind=RealKind) :: Za
!
   if (.not.Initialized) then
      Initialized=.true.
   else
      return
   endif
!
   stop_routine = istop
   print_level = iprint
   GlobalNumAtoms = getNumAtoms()
   LocalNumAtoms = getLocalNumAtoms()
!
!  ===================================================================
   allocate( AtomProperty(LocalNumAtoms) )
   allocate( LizLmax(LocalNumAtoms) )
   allocate( Lmax(LocalNumAtoms) )
   allocate( MixParam(LocalNumAtoms) )
   allocate( LocalMoment(LocalNumAtoms) )
   allocate( GridData(LocalNumAtoms) )
!  ===================================================================
   do i=1,LocalNumAtoms
      LocalMoment(i)%moment=ZERO
      AtomProperty(i)%Rcut_pseudo = Rcut_pseudo_def
   enddo
!  ===================================================================
!  Setup parameters:
!  kmax_kkr ->   l = 0, 1, ..., lmax_kkr, m = -l, -l+1, ..., l
!
!  jmax_kkr ->   l = 0, 1, ..., lmax_kkr, m =  0,   1,  ..., l
!
!  kmax_phi ->   l = 0, 1, ..., lmax_phi, m = -l, -l+1, ..., l
!
!  jmax_phi ->   l = 0, 1, ..., lmax_phi, m =  0,   1,  ..., l
!
!  kmax_pot ->   l = 0, 1, ..., lmax_pot, m = -l, -l+1, ..., l
!
!  jmax_pot ->   l = 0, 1, ..., lmax_pot, m =  0,   1,  ..., l
!
!  kmax_rho ->   l = 0, 1, ..., lmax_rho, m = -l, -l+1, ..., l
!
!  jmax_rho ->   l = 0, 1, ..., lmax_rho, m =  0,   1,  ..., l
!  ===================================================================
!  allocate(atom_index(GlobalNumAtoms))
   allocate(lmax_kkr(0:GlobalNumAtoms), lmax_phi(0:GlobalNumAtoms),    &
            lmax_pot(0:GlobalNumAtoms), lmax_rho(0:GlobalNumAtoms),    &
            num_shells(0:GlobalNumAtoms), lmax_shell(0:GlobalNumAtoms),&
            potinname(0:GlobalNumAtoms), potinform(0:GlobalNumAtoms),  &
            potoutname(0:GlobalNumAtoms), potoutform(0:GlobalNumAtoms),&
            alpha_rho(0:GlobalNumAtoms), alpha_pot(0:GlobalNumAtoms),  &
            alpha_mom(0:GlobalNumAtoms), nmax_liz(0:GlobalNumAtoms),   &
            alpha_chg(0:GlobalNumAtoms), ndivin(0:GlobalNumAtoms),     & 
            ndivout(0:GlobalNumAtoms), nmult(0:GlobalNumAtoms))
   allocate(ind_lmax_kkr(GlobalNumAtoms), ind_lmax_phi(GlobalNumAtoms),    &
            ind_lmax_pot(GlobalNumAtoms), ind_lmax_rho(GlobalNumAtoms),    &
            ind_num_shells(GlobalNumAtoms), ind_lmax_shell(GlobalNumAtoms),&
            ind_potinname(GlobalNumAtoms),  ind_potinform(GlobalNumAtoms), &
            ind_potoutname(GlobalNumAtoms), ind_potoutform(GlobalNumAtoms),&
            ind_alpha_rho(GlobalNumAtoms), ind_alpha_pot(GlobalNumAtoms),  &
            ind_alpha_chg(GlobalNumAtoms), ind_alpha_mom(GlobalNumAtoms),  &
            ind_nmax_liz(GlobalNumAtoms), ind_nmult(GlobalNumAtoms),       &
            ind_ndivin(GlobalNumAtoms), ind_ndivout(GlobalNumAtoms))
   allocate(ind_lmax_step(GlobalNumAtoms),lmax_step(0:GlobalNumAtoms))
   allocate(ind_lmax_pot_trunc(GlobalNumAtoms), lmax_pot_trunc(0:GlobalNumAtoms))
   allocate(ind_cutoff_r(GlobalNumAtoms), cutoff_r(0:GlobalNumAtoms))
   allocate(ind_potScreen(GlobalNumAtoms), potScreen(0:GlobalNumAtoms))
   allocate(ind_pseudo_r(GlobalNumAtoms), pseudo_r(0:GlobalNumAtoms))
   allocate(ind_cutoff_r_s(GlobalNumAtoms), cutoff_r_s(0:GlobalNumAtoms))
   allocate(ind_rmt_desire(GlobalNumAtoms), rmt_desire(0:GlobalNumAtoms))
   allocate(ind_rcr_desire(GlobalNumAtoms), rcr_desire(0:GlobalNumAtoms))
!
!  -------------------------------------------------------------------
!  rstatus = getKeyValue(info_id,'Atom Index',atom_index,GlobalNumAtoms)
!  -------------------------------------------------------------------
   if (getKeyValue(info_id,'Default Potential Input File Name',potinname(0)) /= 0) then
      call ErrorHandler('initAtom','Input potential file name is missing from input')
   endif
   rstatus = getKeyIndexValue(info_id,'Potential Input File Name',    &
                              ind_potinname,potinname(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Potential Input File Form',potinform(0)) /= 0) then
      call ErrorHandler('initAtom','Input potential file form is missing from input')
   endif
   rstatus = getKeyIndexValue(info_id,'Potential Input File Form',    &
                              ind_potinform,potinform(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Potential Output File Name',potoutname(0)) /= 0) then
      call ErrorHandler('initAtom','Output potential file name is missing from input')
   endif
   rstatus = getKeyIndexValue(info_id,'Potential Output File Name',   &
                              ind_potoutname,potoutname(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Potential Output File Form',potoutform(0)) /= 0) then
      call ErrorHandler('initAtom','Output potential file form is missing from input')
   endif
   rstatus = getKeyIndexValue(info_id,'Potential Output File Form',   &
                              ind_potoutform,potoutform(1:GlobalNumAtoms), GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Lmax-T matrix',lmax_kkr(0)) /= 0) then
      call ErrorHandler('initAtom','Lmax for T-matrix is missing from input')
   endif
   ind_lmax_kkr = 0
   rstatus = getKeyIndexValue(info_id,'Lmax-T matrix',                &
                              ind_lmax_kkr,lmax_kkr(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Lmax-Potential',lmax_pot(0),default_param=.false.) /= 0) then
      if (getPotentialTypeParam() == ASA .or. getPotentialTypeParam() == MuffinTin .or. &
          getPotentialTypeParam() == MuffinTinASA) then
         lmax_pot(0) = 0
      else
         lmax_pot(0) = 2*lmax_kkr(0)
      endif
   endif
   ind_lmax_pot = 0
   rstatus = getKeyIndexValue(info_id,'Lmax-Potential',               &
                              ind_lmax_pot,lmax_pot(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Lmax-Wave Func',lmax_phi(0),default_param=.false.) /= 0) then
      lmax_phi(0) = lmax_kkr(0)
   endif
   ind_lmax_phi = 0
   rstatus = getKeyIndexValue(info_id,'Lmax-Wave Func',               &
                              ind_lmax_phi,lmax_phi(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Lmax-Step Func',lmax_step(0),default_param=.false.) /= 0) then
!     call ErrorHandler('initAtom','Lmax for step function is missing from input')
      lmax_step(0) = 4*lmax_kkr(0)
   endif
   ind_lmax_step = 0
   rstatus = getKeyIndexValue(info_id,'Lmax-Step Func',               &
                              ind_lmax_step,lmax_step(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Lmax-Trunc Pot',lmax_pot_trunc(0),default_param=.false.) /= 0) then
      lmax_pot_trunc(0) = lmax_pot(0)
   endif
   ind_lmax_pot_trunc = 0
   rstatus = getKeyIndexValue(info_id,'Lmax-Trunc Pot',               &
                              ind_lmax_pot_trunc,lmax_pot_trunc(1:GlobalNumAtoms),GlobalNumAtoms)
   if (rstatus /= 0) then
      lmax_pot_trunc(1:GlobalNumAtoms) = lmax_pot(1:GlobalNumAtoms) + 4
   else 
      do ig = 1, GlobalNumAtoms
         if (lmax_pot_trunc(ig) < 0) then
            lmax_pot_trunc(ig) = lmax_pot(ig)
         else if (lmax_pot_trunc(ig) > lmax_pot(ig) + lmax_step(ig)) then
            lmax_pot_trunc(ig) = lmax_pot(ig)
!           lmax_pot_trunc(ig) = lmax_pot(ig) + lmax_step(ig)
         endif
      enddo
   endif
!
   if (getKeyValue(info_id,'Default Lmax-Charge Den',lmax_rho(0),default_param=.false.) /= 0) then
      lmax_rho(0) = lmax_pot(0)
   endif
   ind_lmax_rho = 0
   rstatus = getKeyIndexValue(info_id,'Lmax-Charge Den',              &
                              ind_lmax_rho,lmax_rho(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default LIZ # Neighbors',nmax_liz(0)) /= 0) then
      call ErrorHandler('initAtom','Default LIZ # Neighbors is missing from input')
   endif
   ind_nmax_liz = 0
   rstatus = getKeyIndexValue(info_id,'LIZ # Neighbors',              &
                              ind_nmax_liz,nmax_liz(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default LIZ # NN Shells',num_shells(0),default_param=.false.) /= 0) then
      num_shells(0) = 8
   endif
   ind_num_shells = 0
   rstatus = getKeyIndexValue(info_id,'LIZ # NN Shells',              &
                              ind_num_shells,num_shells(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default LIZ Shell Lmax',lmax_shell(0)) /= 0) then
      write(s2,'(i2)')lmax_kkr(0)
      lmax_shell(0) = s2
      do i = 2, num_shells(0)
         lmax_shell(0) = trim(lmax_shell(0))//' '//s2
      enddo
   endif
   ind_lmax_shell = 0
   rstatus = getKeyIndexValue(info_id,'LIZ Shell Lmax',               &
                              ind_lmax_shell,lmax_shell(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default LIZ Cutoff Radius',cutoff_r(0)) /= 0) then
      call ErrorHandler('initAtom','Default LIZ Cutoff Radius is missing from input')
   endif
   ind_cutoff_r = 0
   rstatus = getKeyIndexValue(info_id,'LIZ Cutoff Radius',            &
                              ind_cutoff_r,cutoff_r(1:GlobalNumAtoms),GlobalNumAtoms)
!
   rstatus = getKeyValue(info_id,'Default Rcut-Screen',cutoff_r_s(0))
   ind_cutoff_r_s = 0
   rstatus = getKeyIndexValue(info_id,'Rcut-Screen',                  &
                              ind_cutoff_r_s,cutoff_r_s(1:GlobalNumAtoms),GlobalNumAtoms)
!
   rstatus = getKeyValue(info_id,'Default Pseudo Charge Radius',pseudo_r(0))
   ind_pseudo_r = 0
   rstatus = getKeyIndexValue(info_id,'Pseudo Charge Radius',         &
                              ind_pseudo_r,pseudo_r(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Rho  Mix Param.',alpha_rho(0),default_param=.false.) /= 0) then
      if (getKeyValue(info_id,'Default Mixing Parameter',alpha_rho(0)) /= 0) then
         call ErrorHandler('initAtom','Rho  Mix Param. is missing from input')
      endif
   endif
   ind_alpha_rho = 0
   rstatus = getKeyIndexValue(info_id,'Rho  Mix Param.',              &
                              ind_alpha_rho,alpha_rho(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Pot  Mix Param.',alpha_pot(0),default_param=.false.) /= 0) then
      if (getKeyValue(info_id,'Default Mixing Parameter',alpha_pot(0)) /= 0) then
         call ErrorHandler('initAtom','Pot  Mix Param. is missing from input')
      endif
   endif
   ind_alpha_pot = 0
   rstatus = getKeyIndexValue(info_id,'Pot  Mix Param.',              &
                              ind_alpha_pot,alpha_pot(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Mom  Mix Param.',alpha_mom(0),default_param=.false.) /= 0) then
      if (getKeyValue(info_id,'Default Mixing Parameter',alpha_mom(0)) /= 0) then
         call ErrorHandler('initAtom','Mom  Mix Param. is missing from input')
      endif
   endif
   ind_alpha_mom = 0
   rstatus = getKeyIndexValue(info_id,'Mom  Mix Param.',              &
                              ind_alpha_mom,alpha_mom(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Chg  Mix Param.',alpha_chg(0)) /= 0) then
      call ErrorHandler('initAtom','Default Chg  Mix Param is missing from input')
   endif
   ind_alpha_chg = 0
   rstatus = getKeyIndexValue(info_id,'Chg  Mix Param.',              &
                              ind_alpha_chg,alpha_chg(1:GlobalNumAtoms),GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default No. Rad Points ndivin',ndivin(0)) /= 0) then
      call ErrorHandler('initAtom','Default No. Rad Points ndivin is missing from input')
   endif
   ind_ndivin = 0
   rstatus = getKeyIndexValue(info_id,'No. Rad Points ndivin',        &
                              ind_ndivin,ndivin(1:GlobalNumAtoms), GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default No. Rad Points ndivout',ndivout(0)) /= 0) then
      call ErrorHandler('initAtom','No. Rad Points ndivout is missing from input')
   endif
   ind_ndivout = 0
   rstatus = getKeyIndexValue(info_id,'No. Rad Points ndivout',       &
                              ind_ndivout,ndivout(1:GlobalNumAtoms), GlobalNumAtoms)
!
   if (getKeyValue(info_id,'Default Integer Factor nmult',nmult(0)) /= 0) then
      call ErrorHandler('initAtom','Integer Factor nmult is missing from input')
   endif
   ind_nmult = 0
   rstatus = getKeyIndexValue(info_id,'Integer Factor nmult',         &
                              ind_nmult,nmult(1:GlobalNumAtoms), GlobalNumAtoms)
!
   rstatus = getKeyValue(info_id,'Default Screen Potential',potScreen(0))
   ind_potScreen = 0
   rstatus = getKeyIndexValue(info_id,'Screen Potential',             &
                              ind_potScreen,potScreen(1:GlobalNumAtoms),GlobalNumAtoms)
!
   rstatus = getKeyValue(info_id,'Desired Muffin-tin Radius',rmt_desire(0))
   ind_rmt_desire = 0
   rstatus = getKeyIndexValue(info_id,'Desired Muffin-tin Radius',    &
                              ind_rmt_desire,rmt_desire(1:GlobalNumAtoms),GlobalNumAtoms)
!
   rstatus = getKeyValue(info_id,'Desired Core Radius',rcr_desire(0))
   ind_rcr_desire = 0
   rstatus = getKeyIndexValue(info_id,'Desired Core Radius',          &
                              ind_rcr_desire,rcr_desire(1:GlobalNumAtoms),GlobalNumAtoms)
!
!  ===================================================================
!  Process ind_potinname(:) to taking care of the fact that potinname(0)
!  might contain a list of formatted potential file names
!  ===================================================================
   do n = 1, LocalNumAtoms
      ig = getGlobalIndex(n)
      AtomProperty(n)%GlobalIndex = ig
      AtomProperty(n)%NumSpecies = getNumAlloyElements(ig)
      if (AtomProperty(n)%NumSpecies > 1 .and. .not. isKKRCPA()) then
!        -------------------------------------------------------------
         call ErrorHandler('initAtomModule','Number of Species > 1 for global index',ig)
!        -------------------------------------------------------------
      endif
      allocate( AtomProperty(n)%potinfile(AtomProperty(n)%NumSpecies) )
      allocate( AtomProperty(n)%potoutfile(AtomProperty(n)%NumSpecies) )
      allocate( AtomProperty(n)%NickName(AtomProperty(n)%NumSpecies) )
      allocate( AtomProperty(n)%AtomName(AtomProperty(n)%NumSpecies) )
      allocate( AtomProperty(n)%AtomicNumber(AtomProperty(n)%NumSpecies) )
      allocate( AtomProperty(n)%AtomContent(AtomProperty(n)%NumSpecies) )
   enddo
!
   if (potinform(0) == 0) then
      call initString(potinname(0))
      nt = getNumTokens()
      if (nt >= 1) then
         do j = 1, nt
!           ==========================================================
            call readToken(j,fname,flen)
            path_fname = trim(inputpath)//fname
            funit = 300+j
            open(unit=funit,file=path_fname,status='old',form='formatted')
            read(funit,'(a)')dummy
            read(funit,'(a)')dummy
            read(funit,'(a)')dummy
            read(funit,*)Za
            close(unit=funit)
!           ==========================================================
            anum = int(Za)
            do n = 1, LocalNumAtoms
               ig = getGlobalIndex(n)
               if ( ind_potinname(ig) == 0 ) then
                  do ia = 1, AtomProperty(n)%NumSpecies
                     if ( anum == getZtot(getAlloyElementName(ig,ia)) ) then
                        AtomProperty(n)%potinfile(ia) = path_fname
                     endif
                  enddo
               endif
            enddo
         enddo
      endif
      call endString()
!
      do n = 1, LocalNumAtoms
         ig = getGlobalIndex(n)
         if ( ind_potinname(ig) > 0 ) then
            call initString(potinname(ind_potinname(ig)))
            nt = getNumTokens()
            if (nt < AtomProperty(n)%NumSpecies) then
!              -------------------------------------------------------
               call ErrorHandler('initAtom','The potential files are not completely provided', &
                                 potinname(ig))
!              -------------------------------------------------------
            else if (nt == 1) then
               AtomProperty(n)%potinfile(1) = trim(inputpath)//adjustl(potinname(ind_potinname(ig)))
            else
               do j = 1, nt
!                 ====================================================
                  call readToken(j,fname,flen)
                  path_fname = trim(inputpath)//fname
                  funit = 300+j
                  open(unit=funit,file=path_fname,status='old',form='formatted')
                  read(funit,'(a)')dummy
                  read(funit,'(a)')dummy
                  read(funit,'(a)')dummy
                  read(funit,*)Za
                  close(unit=funit)
!                 ====================================================
                  anum = int(Za)
                  do ia = 1, AtomProperty(n)%NumSpecies
                     if ( anum == getZtot(getAlloyElementName(ig,ia)) ) then
                        AtomProperty(n)%potinfile(ia) = path_fname
                     endif
                  enddo
               enddo
            endif
            call endString()
         else
            do ia = 1, AtomProperty(n)%NumSpecies
               if (len(AtomProperty(n)%potinfile(ia)) <= 1) then
!                  ---------------------------------------------------
                   call ErrorHandler('initAtom','The potential files are not completely provided')
!                  ---------------------------------------------------
               endif
            enddo
         endif
      enddo
   else
      do n = 1, LocalNumAtoms
         ig = getGlobalIndex(n)
         do ia = 1, AtomProperty(n)%NumSpecies
            AtomProperty(n)%potinfile(ia) = trim(inputpath)//adjustl(potinname(ind_potinname(ig)))
         enddo
      enddo
   endif
!  ===================================================================
!
   lmax_max=0
   do n = 1, LocalNumAtoms
      ig = getGlobalIndex(n)
      do ia = 1, AtomProperty(n)%NumSpecies
         AtomProperty(n)%potoutfile(ia)=trim(inputpath)//adjustl(potoutname(ind_potoutname(ig)))
         AtomProperty(n)%NickName(ia)=getAlloyElementName(ig,ia)
         AtomProperty(n)%AtomName(ia)=AtomProperty(n)%NickName(ia)
         AtomProperty(n)%AtomicNumber(ia)=getZtot(AtomProperty(n)%AtomName(ia))
         AtomProperty(n)%AtomContent(ia)=getAlloyElementContent(ig,ia)
      enddo
!
      if (potinform(ind_potinform(ig)) == 0) then
         AtomProperty(n)%potinform = 'FORMATTED'
      else if (potinform(ind_potinform(ig)) == 1) then
         AtomProperty(n)%potinform = 'XDR'
      else if (potinform(ind_potinform(ig)) == 2) then
         AtomProperty(n)%potinform = 'HDF'
      else if (potinform(ind_potinform(ig)) == 3) then
         AtomProperty(n)%potinform = 'UNFORMATTED'
      else
!        -------------------------------------------------------------
         call ErrorHandler('initAtom','unknown potential form',       &
                           potinform(ind_potinform(ig)))
!        -------------------------------------------------------------
      endif
      if (potoutform(ind_potoutform(ig)) == 0) then
         AtomProperty(n)%potoutform = 'FORMATTED'
      else if (potoutform(ind_potoutform(ig)) == 1) then
         AtomProperty(n)%potoutform = 'XDR'
      else if (potoutform(ind_potoutform(ig)) == 2) then
         AtomProperty(n)%potoutform = 'HDF'
      else if (potoutform(ind_potoutform(ig)) == 3) then
         AtomProperty(n)%potoutform = 'UNFORMATTED'
      else
!        -------------------------------------------------------------
         call ErrorHandler('initAtom','unknown potential form',       &
                           potoutform(ind_potoutform(ig)))
!        -------------------------------------------------------------
      endif
!
      AtomProperty(n)%Position(1:3)=getAtomPosition(ig)
      AtomProperty(n)%potScreen=potScreen(ind_potScreen(ig))
      if ( rmt_desire(ind_rmt_desire(ig)) > 0.10d0 ) then
         AtomProperty(n)%Rmt=rmt_desire(ind_rmt_desire(ig))
      else
         AtomProperty(n)%Rmt=ZERO
      endif
      if ( rcr_desire(ind_rcr_desire(ig)) > 0.10d0 ) then
         AtomProperty(n)%Rcore=rcr_desire(ind_rcr_desire(ig))
      else
         AtomProperty(n)%Rcore=ZERO
      endif
      if ( pseudo_r(ind_pseudo_r(ig)) >= 0.5d0 ) then
         AtomProperty(n)%Rcut_pseudo=pseudo_r(ind_pseudo_r(ig))
      else
         AtomProperty(n)%Rcut_pseudo=ZERO
      endif
      Lmax(n)%lmax_kkr=lmax_kkr(ind_lmax_kkr(ig))
      Lmax(n)%lmax_step=lmax_step(ind_lmax_step(ig))
      Lmax(n)%lmax_phi=lmax_phi(ind_lmax_phi(ig))
      if (Lmax(n)%lmax_phi < Lmax(n)%lmax_kkr .and. Lmax(n)%lmax_phi < 0) then
         Lmax(n)%lmax_phi = Lmax(n)%lmax_kkr
      endif
      Lmax(n)%lmax_pot=lmax_pot(ind_lmax_pot(ig))
      Lmax(n)%lmax_pot_trunc=lmax_pot_trunc(ind_lmax_pot_trunc(ig))
      Lmax(n)%lmax_rho=lmax_rho(ind_lmax_rho(ig))
      if (Lmax(n)%lmax_rho > 2*Lmax(n)%lmax_phi .or. Lmax(n)%lmax_rho < 0) then
         Lmax(n)%lmax_rho =  2*Lmax(n)%lmax_phi
      endif
      lmax_max=max(lmax_max,Lmax(n)%lmax_phi, Lmax(n)%lmax_pot, &
                   Lmax(n)%lmax_kkr,Lmax(n)%lmax_rho, Lmax(n)%lmax_step)
      lmax_n = max( Lmax(n)%lmax_pot,Lmax(n)%lmax_rho )
      LizLmax(n)%nmax=nmax_liz(ind_nmax_liz(ig))
      LizLmax(n)%NumShells=num_shells(ind_num_shells(ig))
      LizLmax(n)%rad = cutoff_r(ind_cutoff_r(ig))
      LizLmax(n)%rad_s = cutoff_r_s(ind_cutoff_r_s(ig))
      allocate( LizLmax(n)%lmax_shell( 1:LizLmax(n)%NumShells ) )
      read(lmax_shell(ind_lmax_shell(ig)),*) &
          LizLmax(n)%lmax_shell(1:LizLmax(n)%NumShells)
      MixParam(n)%alpha_rho=alpha_rho(ind_alpha_rho(ig))
      MixParam(n)%alpha_pot=alpha_pot(ind_alpha_pot(ig))
      MixParam(n)%alpha_mom=alpha_mom(ind_alpha_mom(ig))
      MixParam(n)%alpha_chg=alpha_chg(ind_alpha_chg(ig))
      GridData(n)%ndivin=ndivin(ind_ndivin(ig))
      GridData(n)%ndivout=ndivout(ind_ndivout(ig))
      GridData(n)%nmult=nmult(ind_nmult(ig))
   enddo
!
   do n = 1, LocalNumAtoms
      ig = getGlobalIndex(n)
      LocalMoment(n)%evec_old(1:3) = getMomentDirection(ig)
      LocalMoment(n)%evec_new(1:3) = LocalMoment(n)%evec_old(1:3)
      LocalMoment(n)%evec_out(1:3) = LocalMoment(n)%evec_old(1:3)
      LocalMoment(n)%b_con(1:3) = getConstrainField(ig)
      MixParam(n)%alpha_evec = getMomentDirectionMixingParam(ig)
   enddo
!
   InPotFileName = trim(inputpath)//adjustl(potinname(0))
   OutPotFileName = trim(inputpath)//adjustl(potoutname(0))
   if (potinform(0) == 0) then
      InPotFileForm = 'FORMATTED'
   else if (potinform(0) == 1) then
      InPotFileForm = 'XDR'
   else if (potinform(0) == 2) then
      InPotFileForm = 'HDF'
   else if (potinform(0) == 3) then
      InPotFileForm = 'UNFORMATTED'
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initAtom','unknown potential form',potinform(0))
!     ----------------------------------------------------------------
   endif
   if (potoutform(0) == 0) then
      OutPotFileForm = 'FORMATTED'
   else if (potoutform(0) == 1) then
      OutPotFileForm = 'XDR'
   else if (potoutform(0) == 2) then
      OutPotFileForm = 'HDF'
   else if (potoutform(0) == 3) then
      OutPotFileForm = 'UNFORMATTED'
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initAtom','unknown potential form',potoutform(0))
!     ----------------------------------------------------------------
   endif
!
   max_rad = maxval(cutoff_r)
!
   deallocate(lmax_kkr,lmax_phi,lmax_pot,lmax_rho, lmax_step, num_shells, &
              lmax_shell, alpha_rho, alpha_pot, nmax_liz,                 &
              potinname, potinform, potoutname, potoutform,               &
              alpha_mom, alpha_chg, ndivin, ndivout, nmult)
   deallocate(ind_lmax_pot_trunc, lmax_pot_trunc)
   deallocate(ind_lmax_kkr, ind_lmax_phi, ind_lmax_pot, ind_lmax_rho,     &
              ind_num_shells, ind_lmax_shell, ind_lmax_step,              &
              ind_alpha_rho, ind_alpha_pot, ind_alpha_mom, ind_nmax_liz,  &
              ind_potinname, ind_potinform, ind_potoutname, ind_potoutform, &
              ind_ndivin, ind_ndivout,ind_nmult, ind_alpha_chg)
!
   deallocate(ind_cutoff_r, cutoff_r, ind_potScreen, potScreen )
   deallocate(ind_cutoff_r_s, cutoff_r_s)
   deallocate(ind_pseudo_r, pseudo_r)
   deallocate(ind_rmt_desire,rmt_desire)
   deallocate(ind_rcr_desire,rcr_desire)
!
   end subroutine initAtom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endAtom()
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: i
!
   do i=1,LocalNumAtoms
      deallocate( LizLmax(i)%lmax_shell )
   enddo
   deallocate( AtomProperty, LizLmax, Lmax, LocalMoment, MixParam, GridData )
!
   end subroutine endAtom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLocalMoment_s(id,moment)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), intent(in) :: moment
!
   LocalMoment(id)%moment=moment
!
   end subroutine setLocalMoment_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLocalMoment_v(id,moment)
!  ===================================================================
   use MathParamModule, only : TEN2m4, ZERO, ONE
!
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), intent(in) :: moment(3)
   real (kind=RealKind) :: m
!
   m=sqrt( moment(1)*moment(1)+moment(2)*moment(2)+moment(3)*moment(3))
   LocalMoment(id)%moment=m
   if (m < TEN2m4) then
      LocalMoment(id)%evec_out(1)=LocalMoment(id)%evec_old(1)
      LocalMoment(id)%evec_out(2)=LocalMoment(id)%evec_old(2)
      LocalMoment(id)%evec_out(3)=LocalMoment(id)%evec_old(3)
   else
      LocalMoment(id)%evec_out(1)=moment(1)/m
      LocalMoment(id)%evec_out(2)=moment(2)/m
      LocalMoment(id)%evec_out(3)=moment(3)/m
   endif
!
   end subroutine setLocalMoment_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalMoment(id) result(moment)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: moment
!
   moment=LocalMoment(id)%moment
!
   end function getLocalMoment
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLocalEvecOld(id,evec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), intent(in) :: evec(3)
!
   LocalMoment(id)%evec_old(1)=evec(1)
   LocalMoment(id)%evec_old(2)=evec(2)
   LocalMoment(id)%evec_old(3)=evec(3)
!
   end subroutine setLocalEvecOld
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalEvecOld(id) result(evec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: evec(3)
!
   evec(1)=LocalMoment(id)%evec_old(1)
   evec(2)=LocalMoment(id)%evec_old(2)
   evec(3)=LocalMoment(id)%evec_old(3)
!
   end function getLocalEvecOld
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLocalEvecNew(id,evec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), intent(in) :: evec(3)
!
   LocalMoment(id)%evec_new(1)=evec(1)
   LocalMoment(id)%evec_new(2)=evec(2)
   LocalMoment(id)%evec_new(3)=evec(3)
!
   end subroutine setLocalEvecNew
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalEvecNew(id) result(evec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: evec(3)
!
   evec(1)=LocalMoment(id)%evec_new(1)
   evec(2)=LocalMoment(id)%evec_new(2)
   evec(3)=LocalMoment(id)%evec_new(3)
!
   end function getLocalEvecNew
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLocalEvecOut(id,evec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), intent(in) :: evec(3)
!
   LocalMoment(id)%evec_out(1)=evec(1)
   LocalMoment(id)%evec_out(2)=evec(2)
   LocalMoment(id)%evec_out(3)=evec(3)
!
   end subroutine setLocalEvecOut
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalEvecOut(id) result(evec)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: evec(3)
!
   evec(1)=LocalMoment(id)%evec_out(1)
   evec(2)=LocalMoment(id)%evec_out(2)
   evec(3)=LocalMoment(id)%evec_out(3)
!
   end function getLocalEvecOut
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLocalConstrainField(id,b_con)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind), intent(in) :: b_con(3)
!
   LocalMoment(id)%b_con(1)=b_con(1)
   LocalMoment(id)%b_con(2)=b_con(2)
   LocalMoment(id)%b_con(3)=b_con(3)
!
   end subroutine setLocalConstrainField
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalConstrainField(id) result(b_con)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: b_con(3)
!
   b_con(1)=LocalMoment(id)%b_con(1)
   b_con(2)=LocalMoment(id)%b_con(2)
   b_con(3)=LocalMoment(id)%b_con(3)
!
   end function getLocalConstrainField
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalAtomicNumber(id,ia) result(z)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in), optional :: ia
!
   integer (kind=IntKind) :: z
!
   if (present(ia)) then
      z=AtomProperty(id)%AtomicNumber(ia)
   else
      z=AtomProperty(id)%AtomicNumber(1)
   endif
!
   end function getLocalAtomicNumber
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalAtomName(id,ia) result(name)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in), optional :: ia
!
   character (len=2) :: name
!
   if (present(ia)) then
      name=AtomProperty(id)%AtomName(ia)
   else
      name=AtomProperty(id)%AtomName(1)
   endif
!
   end function getLocalAtomName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalAtomNickName(id,ia) result(name)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in), optional :: ia
!
   character (len=MaxLenOfAtomName) :: name
!
   if (present(ia)) then
      name=AtomProperty(id)%NickName(ia)
   else
      name=AtomProperty(id)%NickName(1)
   endif
!
   end function getLocalAtomNickName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalNumSpecies(id) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: n
!
   n = AtomProperty(id)%NumSpecies
!
   end function getLocalNumSpecies
!  ===================================================================
!
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalSpeciesContent(id,ia) result(c)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, ia
!
   real (kind=RealKind) :: c
!
   c = AtomProperty(id)%AtomContent(ia)
!
   end function getLocalSpeciesContent
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalAtomPosition(id) result(position)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: position(3)
!
   position(1)=AtomProperty(id)%Position(1)
   position(2)=AtomProperty(id)%Position(2)
   position(3)=AtomProperty(id)%Position(3)
!
   end function getLocalAtomPosition
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKKRLmax(id) result(l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getKKRLmax','Invalid atom index',id)
   endif
   l = Lmax(id)%lmax_kkr
   end function getKKRLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStepFuncLmax(id) result(l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getPhiLmax','Invalid atom index',id)
   endif
   l = Lmax(id)%lmax_step
   end function getStepFuncLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPhiLmax(id) result(l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getPhiLmax','Invalid atom index',id)
   endif
   l = Lmax(id)%lmax_phi
   end function getPhiLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPhiLmax(id,l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('setPhiLmax','Invalid atom index',id)
   else if (l < 0) then
      call ErrorHandler('setPhiLmax','Invalid l value',l)
   endif
   Lmax(id)%lmax_phi = l
   end subroutine setPhiLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRhoLmax(id) result(l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getRhoLmax','Invalid atom index',id)
   endif
   l = Lmax(id)%lmax_rho
   end function getRhoLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setRhoLmax(id,l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('setRhoLmax','Invalid atom index',id)
   else if (l < 0) then
      call ErrorHandler('setRhoLmax','Invalid l value',l)
   endif
   Lmax(id)%lmax_rho = l
   end subroutine setRhoLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotLmax(id) result(l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getPotLmax','Invalid atom index',id)
   endif
   l = Lmax(id)%lmax_pot
   end function getPotLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPotLmax(id,l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('setPotLmax','Invalid atom index',id)
   else if (l < 0) then
      call ErrorHandler('setPotLmax','Invalid lmax value',l)
   endif
   Lmax(id)%lmax_pot = l
   end subroutine setPotLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTruncPotLmax(id) result(l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getPotLmax','Invalid atom index',id)
   endif
   l = Lmax(id)%lmax_pot_trunc
   end function getTruncPotLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setTruncPotLmax(id,l)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, l
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('setTruncPotLmax','Invalid atom index',id)
   else if (l < 0) then
      call ErrorHandler('setTruncPotLmax','Invalid lmax value',l)
   endif
   Lmax(id)%lmax_pot_trunc = l
   end subroutine setTruncPotLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLizLmax(id) result(s)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   type (LizLmaxStruct), pointer :: s
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getLizLmax','Invalid atom index',id)
   endif
   s => LizLmax(id)
   end function getLizLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMaxLmax() result(l)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: l
   l=lmax_max
!
   end function getMaxLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIPFN0() result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
!
   s=InPotFileName
!
   end function getIPFN0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIPFN1(id) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getInPotFileName','Invalid atom index',id)
   else if (AtomProperty(id)%NumSpecies /= 1) then
      call ErrorHandler('getInPotFileName','Alloy element is not specified')
   endif
   s=AtomProperty(id)%potinfile(1)
!
   end function getIPFN1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIPFN2(id,ia) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id, ia
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getInPotFileName','Invalid atom index',id)
   else if (ia < 1 .or. ia > AtomProperty(id)%NumSpecies) then
      call ErrorHandler('getInPotFileName','Invalid alloy element index',ia)
   endif
   s=AtomProperty(id)%potinfile(ia)
!
   end function getIPFN2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIPFF0() result(s)
!  ===================================================================
   implicit none
   character (len=11) :: s
!
   s=InPotFileForm
!
   end function getIPFF0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIPFF1(id) result(s)
!  ===================================================================
   implicit none
   character (len=11) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getInPotFileForm','Invalid atom index',id)
   endif
   s=AtomProperty(id)%potinform
!
   end function getIPFF1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOPFN0() result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
!
   s=OutPotFileName
!
   end function getOPFN0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOPFN1(id) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getOutPotFileName','Invalid atom index',id)
   else if (AtomProperty(id)%NumSpecies /= 1) then
      call ErrorHandler('getInPotFileName','Alloy element is not specified')
   endif
   s=AtomProperty(id)%potoutfile(1)
!
   end function getOPFN1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOPFN2(id,ia) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id, ia
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getOutPotFileName','Invalid atom index',id)
   else if (ia < 1 .or. ia > AtomProperty(id)%NumSpecies) then
      call ErrorHandler('getInPotFileName','Invalid alloy element index',ia)
   endif
   s=AtomProperty(id)%potoutfile(ia)
!
   end function getOPFN2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOPFF0() result(s)
!  ===================================================================
   implicit none
   character (len=11) :: s
!
   s=OutPotFileForm
!
   end function getOPFF0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOPFF1(id) result(s)
!  ===================================================================
   implicit none
   character (len=11) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getOutPotFileName','Invalid atom index',id)
   endif
   s=AtomProperty(id)%potoutform
!
   end function getOPFF1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIVDFN0(id) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getInValDenFileName','Invalid atom index',id)
   else if (AtomProperty(id)%NumSpecies /= 1) then
      call ErrorHandler('getInPotFileName','Alloy element is not specified')
   endif
   s=AtomProperty(id)%potinfile(1)
!
   end function getIVDFN0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getIVDFN1(id,ia) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id, ia
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getInValDenFileName','Invalid atom index',id)
   else if (ia < 1 .or. ia > AtomProperty(id)%NumSpecies) then
      call ErrorHandler('getInPotFileName','Invalid alloy element index',ia)
   endif
   s=AtomProperty(id)%potinfile(ia)
!
   end function getIVDFN1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInValDenFileForm(id) result(s)
!  ===================================================================
   implicit none
   character (len=11) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getInValDenFileName','Invalid atom index',id)
   endif
   s=AtomProperty(id)%potinform
!
   end function getInValDenFileForm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOVDFN0(id) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getOutValDenFileName','Invalid atom index',id)
   else if (AtomProperty(id)%NumSpecies /= 1) then
      call ErrorHandler('getInPotFileName','Alloy element is not specified')
   endif
   s=AtomProperty(id)%potoutfile(1)
!
   end function getOVDFN0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOVDFN1(id,ia) result(s)
!  ===================================================================
   implicit none
   character (len=MaxLenFileName) :: s
   integer (kind=IntKind), intent(in) :: id, ia
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getOutValDenFileName','Invalid atom index',id)
   else if (ia < 1 .or. ia > AtomProperty(id)%NumSpecies) then
      call ErrorHandler('getInPotFileName','Invalid alloy element index',ia)
   endif
   s=AtomProperty(id)%potoutfile(ia)
!
   end function getOVDFN1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getOutValDenFileForm(id) result(s)
!  ===================================================================
   implicit none
   character (len=11) :: s
   integer (kind=IntKind), intent(in) :: id
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getOutValDenFileName','Invalid atom index',id)
   endif
   s=AtomProperty(id)%potoutform
!
   end function getOutValDenFileForm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getGridData(id,nin,nout,nmul)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: nin
   integer (kind=IntKind), intent(out) :: nout
   integer (kind=IntKind), intent(out) :: nmul
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getGridData','Invalid atom index',id)
   endif
!
   nin=GridData(id)%ndivin
   nout=GridData(id)%ndivout
   nmul=GridData(id)%nmult
!
   end subroutine getGridData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMixingParam4Rho(id) result(alpha)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: alpha
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getMixingParam4Rho','Invalid atom index',id)
   endif
!
   alpha=MixParam(id)%alpha_rho
!
   end function getMixingParam4Rho
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMixingParam4Pot(id) result(alpha)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: alpha
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getMixingParam4Pot','Invalid atom index',id)
   endif
!
   alpha=MixParam(id)%alpha_pot
!
   end function getMixingParam4Pot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMixingParam4Charge(id) result(alpha)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: alpha
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getMixingParam4Charge','Invalid atom index',id)
   endif
!
   alpha=MixParam(id)%alpha_chg
!
   end function getMixingParam4Charge
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMixingParam4Mom(id) result(alpha)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: alpha
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getMixingParam4Mom','Invalid atom index',id)
   endif
!
   alpha=MixParam(id)%alpha_mom
!
   end function getMixingParam4Mom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMixingParam4Evec(id) result(alpha)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   real (kind=RealKind) :: alpha
!
   if (id<1 .or. id>LocalNumAtoms) then
      call ErrorHandler('getMixingParam4Evec','Invalid atom index',id)
   endif
!
   alpha=MixParam(id)%alpha_evec
!
   end function getMixingParam4Evec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getScreenPotential(i) result (vswiss)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: vswiss
!
   if (i<1 .or. i>LocalNumAtoms) then
      call ErrorHandler('getScreenPotential','Invalid atom index',i)
   endif
   vswiss = AtomProperty(i)%potScreen
!
   end function getScreenPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomMuffinTinRad(i) result (rmt)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: rmt
!
   if (i<1 .or. i>LocalNumAtoms) then
      call ErrorHandler('getAtomMuffinTinRad','Invalid atom index',i)
   endif
   rmt = AtomProperty(i)%Rmt
!
   end function getAtomMuffinTinRad
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setAtomMuffinTinRad(i,rmt)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: rmt
!
   if (i<1 .or. i>LocalNumAtoms) then
      call ErrorHandler('setAtomMuffinTinRad','Invalid atom index',i)
   endif
   AtomProperty(i)%Rmt = rmt
!
   end subroutine setAtomMuffinTinRad
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomCoreRad(i) result (rcr)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: rcr
!
   if (i<1 .or. i>LocalNumAtoms) then
      call ErrorHandler('getAtomCoreRad','Invalid atom index',i)
   endif
   rcr = AtomProperty(i)%Rcore
!
   end function getAtomCoreRad
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setAtomCoreRad(i,rcr)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind), intent(in) :: rcr
!
   if (i<1 .or. i>LocalNumAtoms) then
      call ErrorHandler('setAtomCoreRad','Invalid atom index',i)
   endif
   AtomProperty(i)%Rcore = rcr
!
   end subroutine setAtomCoreRad
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getScreenRcut(i) result (rcut)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: rcut
!
   if (i<1 .or. i>LocalNumAtoms) then
      call ErrorHandler('getScreenRcut','Invalid atom index',i)
   endif
   rcut = LizLmax(i)%rad_s
!
   end function getScreenRcut
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRcutPseudo(i) result (rcut)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: rcut
!
   if (i<1 .or. i>LocalNumAtoms) then
      call ErrorHandler('getScreenPotential','Invalid atom index',i)
   endif
   rcut = AtomProperty(i)%Rcut_pseudo
!
   end function getRcutPseudo
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMaxRad() result (mr)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: mr
!
   mr = max_rad
!
   end function getMaxRad
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printAtom()
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: i, ia
!
   write(6,'(/,80(''-''))')
   write(6,'(/,27x,a)')'*************************'
   write(6,'( 27x,a )')'* Output from printAtom *'
   write(6,'(27x,a,/)')'*************************'
   do i=1,LocalNumAtoms
      write(6,'(80(''=''))')
      do ia = 1, AtomProperty(i)%NumSpecies
         write(6,'(''Atom Name       : '',a)')AtomProperty(i)%AtomName(ia)
         write(6,'(''Input  Potential: '',a50)') AtomProperty(i)%potinfile(ia)
         write(6,'(''Output Potential: '',a50)') AtomProperty(i)%potoutfile(ia)
      enddo
      write(6,'(''Global Index    : '',i6)')AtomProperty(i)%GlobalIndex
      write(6,'(''Atom Position   : '',3f10.5)')AtomProperty(i)%Position(1:3)
      write(6,'(''Local Moment    : '',f10.5)')LocalMoment(i)%moment
      write(6,'(''Local evec_old  : '',3f10.5)')LocalMoment(i)%evec_old(1:3)
      write(6,'(''Constrain Fld   : '',3f10.5)')LocalMoment(i)%b_con(1:3)
      write(6,'(''Potential Mixing Parameter: '',f8.5)')MixParam(i)%alpha_pot
      write(6,'(''Density   Mixing Parameter: '',f8.5)')MixParam(i)%alpha_rho
      write(6,'(''Charge    Mixing Parameter: '',f8.5)')MixParam(i)%alpha_chg
      write(6,'(''Moment    Mixing Parameter: '',f8.5)')MixParam(i)%alpha_mom
      write(6,'(''Evec      Mixing Algorithm: '',f8.5)')MixParam(i)%alpha_evec
      write(6,'(''Number of Shells in LIZ: '',i4)')LizLmax(i)%NumShells
      write(6,'(''Lmax for each LIZ shell: '',10i4)') &
            LizLmax(i)%lmax_shell(1:LizLmax(i)%NumShells)
      write(6,'(''Number of r-mesh less or equal than Rmt: '',i5)')   &
            GridData(i)%ndivin
      write(6,'(''Number of r-mesh greater than Rmt      : '',i5)')   &
            GridData(i)%ndivout
      write(6,'(''Ratio of r-mesh steps, hin and hout    : '',i5)')   &
            GridData(i)%nmult
      if (AtomProperty(i)%Rmt > 0.000001d0) then
         write(6,'(''Desired Rmt         : '',f10.5)') AtomProperty(i)%Rmt
      endif
      if (AtomProperty(i)%Rcore > 0.000001d0) then
         write(6,'(''Desired Core Radius : '',f10.5)') AtomProperty(i)%Rcore
      endif
   enddo
   write(6,'(80(''=''))')
!
   end subroutine printAtom
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printAtomMomentInfo()
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: i, ig
!
   write(6,'(/,80(''-''))')
   write(6,'(/,21x,a)')'*************************************'
   write(6,'( 21x,a )')'* Moment Info Output from printAtom *'
   write(6,'(21x,a,/)')'*************************************'
   do i=1,LocalNumAtoms
      ig = AtomProperty(i)%GlobalIndex
      write(6,'(80(''=''))')
      write(6,'(''Global Index    : '',i6)')AtomProperty(i)%GlobalIndex
      write(6,'(''Local Moment    : '',f12.8)') LocalMoment(i)%moment
      write(6,'(''Local evec_old  : '',3f12.8)')LocalMoment(i)%evec_old(1:3)
      write(6,'(''Local evec_new  : '',3f12.8)')LocalMoment(i)%evec_new(1:3)
      write(6,'(''Local evec_out  : '',3f12.8)')LocalMoment(i)%evec_out(1:3)
      write(6,'(''Constrain Fld   : '',3f12.8)')LocalMoment(i)%b_con(1:3)
      write(6,'(''Moment    Mixing Parameter: '',f8.5)')MixParam(i)%alpha_mom
      write(6,'(''Evec      Mixing Algorithm: '',f8.5)')MixParam(i)%alpha_evec
   enddo
   write(6,'(80(''=''))')
!
   end subroutine printAtomMomentInfo
!  ===================================================================
end module AtomModule
