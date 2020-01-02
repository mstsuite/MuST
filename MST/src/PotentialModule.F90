module PotentialModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : GridStruct
   use MathParamModule, only : ZERO, CZERO, CONE, ONE, TWO, HALF, Y0, &
                               SQRT_PI, TEN2m7, TEN2m8, TEN2m9, PI2, PI4, TEN2m6, TEN, SQRTm1
!
public :: initPotential,         &
          endPotential,          &
          readPotential,         &
          writePotential,        &       ! unfinished
          printPot_L,            &
          setHeader,             &
          getHeader,             &
          setTitle,              &
          getTitle,              &
          getPotential,          &
          getTruncatedPotential, &
          isPotComponentZero,    &
          getPotComponentFlag,   &
          getTruncatedPotComponentFlag,   &
          setPotential,          &
          getSphPotr,            &
          setSphPotr,            &
          setV0,                 &
          getV0,                 &
          setVdif,               &
          getVdif,               &
          getPotLmax,            &
          getTruncPotLmax,       &
          getPotEf,              &
          setPotEf,              &
          setPotentialAccess,    &
          setPotentialOutsideMT, &
          isSphericalInputFile,  &
          pushPotentialToAccel,  &
          pushTruncatedPotToAccel, &
          deletePotentialOnAccel,  &
          deleteTruncatedPotOnAccel
!
   interface getPotential
      module procedure getPoten_l0, getPoten_l1, getPoten_l2,         &
                       getPoten_l3, getPoten_r
   end interface
   interface getTruncatedPotential
      module procedure getTruncPoten_l0, getTruncPoten_l1,            &
                       getTruncPoten_r
   end interface
!
private
!
   logical :: Initialized = .false.
   logical :: SphericalInputFile = .true.
   logical :: isTruncatedPotential = .true.
   logical :: isChargeSymmOn = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind), parameter :: n_inter = 5 ! order of polynomial
                                                    ! interpolation
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: GlobalNumAtoms
!
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind), allocatable :: print_level(:)
   integer (kind=IntKind), allocatable :: SiteMaxNumc(:)
   integer (kind=IntKind) :: jmtmax
   integer (kind=IntKind) :: jwsmax
   integer (kind=IntKind) :: numcmax
   integer (kind=IntKind), allocatable :: lofj(:), mofj(:), kofj(:)
   integer (kind=IntKind), allocatable :: lofk(:), mofk(:), jofk(:)
   integer (kind=IntKind), allocatable :: m1m(:)
   integer (kind=IntKind) :: GroupID
   integer (kind=IntKind) :: ClusterIndexOfGroup
!
!  ===================================================================
!  header:   Statement header of the potential file
!  jtitle:   sub-header of the potential data for each spin index
!  lmax:     angular momentum l cutoff
!  jmax:     = (lmax+1)*(lmax+2)/2
!  potr_sph: the spherical part of potential data times r-mesh
!            Note: pot = potr_sph/r + sum_{lm} pot_l * Y_{lm}
!  v0:       constant shift of potential data for each spin
!  Grid:     pointer to radial grid data
!  pot_l:    potential data for 0 <= l <= lmax
!  ===================================================================
   type PotentialStruct
      character (len=80) :: header
      character (len=80), allocatable :: jtitle(:,:)
      logical, allocatable :: isTruncDone(:)
      integer (kind=IntKind) :: lmax
      integer (kind=IntKind) :: jmax
      integer (kind=IntKind) :: lmax_trunc
      integer (kind=IntKind) :: jmax_trunc
      integer (kind=IntKind) :: iend_trunc
      integer (kind=IntKind) :: NumSpecies
      real (kind=RealKind), allocatable :: ztss(:)
      integer (kind=IntKind), pointer :: PotCompFlag(:)
      integer (kind=IntKind), pointer :: PotCompFlag_trunc(:)
      real (kind=RealKind), pointer :: potr_sph(:,:,:)
      real (kind=RealKind), pointer :: potr_sph_trunc(:,:,:)
      type (GridStruct), pointer :: Grid
      complex (kind=CmplxKind), pointer :: pot_l(:,:,:,:)
      complex (kind=CmplxKind), pointer :: pot_l_trunc(:,:,:,:)
   end type PotentialStruct
!
   type (PotentialStruct), allocatable, target :: Potential(:)
!
   integer (kind=IntKind) :: RECORD_LENGTH
!
   real (kind=RealKind) :: efermi
   real (kind=RealKind), allocatable :: efermi_in(:,:)
   real (kind=RealKind), target :: vdif(1)
   real (kind=RealKind) :: v0(2)
   real (kind=RealKind), parameter :: pot_tol = TEN2m7
!
   integer (kind=IntKind) :: MaxNumSpecies
   integer (kind=IntKind), allocatable :: LdaPlusU_DataAccum(:)
   integer (kind=IntKind), allocatable :: NonSphPot_DataAccum(:)
   integer (kind=IntKind), allocatable :: NumSpecies(:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPotential(na,lmax,lmax_step,npola,ncant,istop,iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, GlobalMaxInGroup
   use GroupCommModule, only : GlobalSumInGroup, getMyClusterIndex
   use ChemElementModule, only : MaxLenOfAtomName, getNumCoreStates
   use SystemModule, only : getAtomName, getNumAtoms, getNumAlloyElements, &
                            getAlloyElementName
   use Atom2ProcModule, only : getGlobalIndex, getMaxLocalNumAtoms,  &
                               getLocalNumAtoms
   use AtomModule, only : getTPL => getTruncPotLmax
   use RadialGridModule, only : getGrid
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       setDataStorageLDA,     &
                                       setDataStorage2Value,  &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark, &
                                       IntegerType, IntegerMark
!
   use LdaCorrectionModule, only : checkLdaCorrection, getDataPackSize
   use PotentialTypeModule, only : isFullPotential
   use ScfDataModule, only : getSingleSiteSolverMethod, isChargeSymm
!
   use RadialGridModule, only : getNumRmesh
!
   implicit none
!
   character (len=*) :: istop
   character (len=MaxLenOfAtomName) :: atname
!
   integer (kind=IntKind), intent(in) :: na, npola, ncant, iprint(na)
   integer (kind=IntKind), intent(in) :: lmax(na), lmax_step(na)
   integer (kind=IntKind) :: id, iend, jl, kl, ia, present_atom, iend_diff
!
   integer (kind=IntKind) :: lmax_max_trunc, lmax_pot_trunc, jmax_pot_trunc
   integer (kind=IntKind) :: jmax_trunc, kmax_trunc
   integer (kind=IntKind) :: lmax_max_step
   integer (kind=IntKind) :: lmax_pot, jmax_pot, jend, jmt
   integer (kind=IntKind) :: lmax_max, jmax_max, l, m, ip, glb_lmax_max
   integer (kind=IntKind) :: numc, ig, msg(5), MaxLocalNumAtoms
!
   integer (kind=IntKind), allocatable :: DataSize(:)
   integer (kind=IntKind), allocatable :: DataSize_trunc(:)
   integer (kind=IntKind), allocatable :: DataSizeComm(:,:)
!
   real (kind=RealKind) :: a = ZERO
!
   type (GridStruct), pointer :: Grid
!
   n_spin_pola = npola
   n_spin_cant = ncant
   LocalNumAtoms=na
   stop_routine = istop
   GlobalNumAtoms = getNumAtoms()
!
   GroupID = getGroupID('Unit Cell')
   ClusterIndexOfGroup = getMyClusterIndex(GroupID)
!
!  if ( getSingleSiteSolverMethod() > 0 ) then
   if ( getSingleSiteSolverMethod() > -2 ) then
      isTruncatedPotential = .true.
   else
      isTruncatedPotential = .false.
   endif
!
   allocate( DataSize(LocalNumAtoms), NumSpecies(LocalNumAtoms) )
   if ( isTruncatedPotential ) then
      allocate( DataSize_trunc(LocalNumAtoms) )
   endif
!
   allocate(print_level(LocalNumAtoms))
   print_level(:) = iprint(:)
!
   isChargeSymmOn = isChargeSymm()
!
   lmax_max = 0
   lmax_max_trunc = 0
   lmax_max_step  = 0
   MaxNumSpecies = 0
   do id = 1, LocalNumAtoms
      lmax_max=max(lmax_max, lmax(id))
      if ( isTruncatedPotential) then
!        lmax_max_trunc = max(lmax_max_trunc,lmax_step(id)+lmax(id))
         lmax_max_trunc = max(lmax_max_trunc,getTPL(id))
         lmax_max_step  = max(lmax_max_step,lmax_step(id))
      endif
      ig = getGlobalIndex(id)
      NumSpecies(id) = getNumAlloyElements(ig)
      MaxNumSpecies = max(MaxNumSpecies, NumSpecies(id))
   enddo
   jmax_max = (lmax_max+1)*(lmax_max+2)/2
   if ( .not.isTruncatedPotential ) then
     lmax_max_trunc = lmax_max
   endif
!
   jmax_trunc = ((lmax_max_trunc+1)*(lmax_max_trunc+2))/2
   kmax_trunc = (lmax_max_trunc+1)*(lmax_max_trunc+1)
   l = max(lmax_max_trunc,lmax_max_step,lmax_max)
   jl = ((l+1)*(l+2))/2
   kl = (l+1)*(l+1)
   allocate( lofj(jl), mofj(jl), kofj(jl) )
   allocate( lofk(kl), mofk(kl), jofk(kl) )
   allocate( m1m(-l:l) )
!
   jl=0
   kl=0
   do l = 0, max(lmax_max_trunc,lmax_max_step,lmax_max)
      do m = -l, l
         kl = kl + 1
         lofk(kl) = l
         mofk(kl) = m
         if (m >= 0) then
            jl = jl + 1
            lofj(jl) = l
            mofj(jl) = m
            kofj(jl) = kl
            jofk(kl) = jl
            jofk(kl-2*m) = jl
         endif
      enddo
   enddo
!
   l = max(lmax_max_trunc,lmax_max_step,lmax_max)
   m1m(-l)=1-2*mod(l,2)
   do m=-l+1,l
      m1m(m)=-m1m(m-1)
   enddo
!
   inquire(iolength=RECORD_LENGTH) a
!
   jmtmax = 0
   jwsmax = 0
   iend = 0
!
   do id = 1,LocalNumAtoms
      Grid => getGrid(id)
      DataSize(id) = Grid%jend*NumSpecies(id)*n_spin_pola
      if ( isTruncatedPotential ) then
         iend_diff = Grid%jend-Grid%jmt+1
         DataSize_trunc(id) = iend_diff*n_spin_pola*NumSpecies(id)
      endif
   enddo
!
   if ( .not.isDataStorageExisting('OldLDAPotential_r') ) then
!     ---------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'OldLDAPotential_r', &
                             DataSize,RealType)
!     ---------------------------------------------------------------
      call setDataStorage2Value('OldLDAPotential_r',zero)
!     ---------------------------------------------------------------
      do id = 1,LocalNumAtoms
         Grid => getGrid(id)
         jend = Grid%jend
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'OldLDAPotential_r',jend)
!        -------------------------------------------------------------
      enddo
      if ( isTruncatedPotential ) then
!        ------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'OldLDATruncPotential_r', &
                                DataSize_trunc,RealType)
!        ------------------------------------------------------------
         call setDataStorage2Value('OldLDATruncPotential_r',zero)
!        ------------------------------------------------------------
         do id = 1,LocalNumAtoms
            Grid => getGrid(id)
            iend_diff = Grid%jend-Grid%jmt+1
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'OldLDATruncPotential_r',iend_diff)
!           ----------------------------------------------------------
         enddo
      endif
   endif
!
   do id = 1,LocalNumAtoms
      l = lmax(id)
      jl= (l+1)*(l+2)/2
      Grid => getGrid(id)
      DataSize(id) = Grid%jend*jl*n_spin_pola*NumSpecies(id)
      if ( isTruncatedPotential ) then
!        l = lmax_step(id)+lmax(id)
         l = getTPL(id)
         jl = ((l+1)*(l+2))/2
         iend_diff = Grid%jend - Grid%jmt + 1
         DataSize_trunc(id) = iend_diff*jl*n_spin_pola*NumSpecies(id)
      endif
   enddo
!
   if ( .not.isDataStorageExisting('OldLDAPotential_l') ) then
!     ---------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'OldLDAPotential_l', &
                             DataSize,ComplexType)
!     ---------------------------------------------------------------
      call setDataStorage2Value('OldLDAPotential_l',czero)
!     ---------------------------------------------------------------
      do id = 1,LocalNumAtoms
         Grid => getGrid(id)
         jend = Grid%jend
!        -------------------------------------------------------------
         call setDataStorageLDA(id,'OldLDAPotential_l',jend)
!        -------------------------------------------------------------
      enddo
      if ( isTruncatedPotential ) then
!        -------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,'OldLDATruncPotential_l', &
                                DataSize_trunc,ComplexType)
!        ------------------------------------------------------------
         call setDataStorage2Value('OldLDATruncPotential_l',czero)
!        ------------------------------------------------------------
         do id = 1,LocalNumAtoms
            Grid => getGrid(id)
            iend_diff = Grid%jend - Grid%jmt + 1
!           ----------------------------------------------------------
            call setDataStorageLDA(id,'OldLDATruncPotential_l',iend_diff)
!           ----------------------------------------------------------
         enddo
      endif
   endif
!
   do id = 1,LocalNumAtoms
      DataSize(id) = (getNumRmesh(id)+1)*NumSpecies(id)
   enddo
!
   if ( .not.isDataStorageExisting('OldSphericalElectronDensity') ) then
!     ---------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'OldSphericalElectronDensity', &
                             DataSize,RealType)
!     ---------------------------------------------------------------
      call setDataStorage2Value('OldSphericalElectronDensity',zero)
!     ---------------------------------------------------------------
   endif
   if (n_spin_pola == 2 .and. &
       .not.isDataStorageExisting('OldSphericalMomentDensity') ) then
!     ---------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,                          &
                             'OldSphericalMomentDensity', DataSize,RealType)
!     ---------------------------------------------------------------
      call setDataStorage2Value('OldSphericalMomentDensity',zero)
!     ---------------------------------------------------------------
   endif
!
   allocate( Potential(LocalNumAtoms), efermi_in(MaxNumSpecies, LocalNumAtoms) )
   allocate( SiteMaxNumc(LocalNumAtoms) )
   SiteMaxNumc(1:LocalNumAtoms) = 0
   numcmax = 0
   do id = 1, LocalNumAtoms
      lmax_pot=lmax(id)
      Grid => getGrid(id)
      jmax_pot=((lmax_pot+1)*(lmax_pot+2))/2
      jend = Grid%jend
      jmt = Grid%jmt
      allocate(Potential(id)%jtitle(n_spin_pola,NumSpecies(id)))
      allocate(Potential(id)%ztss(NumSpecies(id)))
      allocate(Potential(id)%isTruncDone(NumSpecies(id)))
      allocate(Potential(id)%PotCompFlag(1:jmax_pot))
      Potential(id)%NumSpecies = NumSpecies(id)
      Potential(id)%PotCompFlag(1:jmax_pot)=0
      Potential(id)%pot_l=>getDataStorage(id,"OldLDAPotential_l",    &
                                          jend,jmax_pot,n_spin_pola,NumSpecies(id),ComplexMark)
      Potential(id)%potr_sph=>getDataStorage(id,"OldLDAPotential_r", &
                                             jend,n_spin_pola,NumSpecies(id),RealMark)
      Potential(id)%Grid=>Grid
      Potential(id)%lmax=lmax_pot
      Potential(id)%jmax=jmax_pot
      if ( isTruncatedPotential ) then
         iend_diff = jend - jmt + 1
!        lmax_pot_trunc = (lmax(id)+lmax_step(id))
         lmax_pot_trunc = getTPL(id)
         jmax_pot_trunc = ((lmax_pot_trunc+1)*(lmax_pot_trunc+2))/2
         Potential(id)%lmax_trunc = lmax_pot_trunc
         Potential(id)%jmax_trunc = jmax_pot_trunc
         Potential(id)%iend_trunc = iend_diff
         allocate(Potential(id)%PotCompFlag_trunc(1:jmax_pot_trunc))
         Potential(id)%pot_l_trunc=>getDataStorage(id,"OldLDATruncPotential_l", &
                            iend_diff,jmax_pot_trunc,n_spin_pola,NumSpecies(id),ComplexMark)
         Potential(id)%potr_sph_trunc=>getDataStorage(id,"OldLDATruncPotential_r", &
                                        iend_diff,n_spin_pola,NumSpecies(id),RealMark)
         Potential(id)%isTruncDone(1:NumSpecies(id)) = .false.
         Potential(id)%PotCompFlag_trunc(1:jmax_pot_trunc)=0
      else
         nullify(Potential(id)%PotCompFlag_trunc)
         nullify(Potential(id)%pot_l_trunc)
         nullify(Potential(id)%potr_sph_trunc)
         Potential(id)%isTruncDone = .false.
         Potential(id)%lmax_trunc=-1
         Potential(id)%jmax_trunc=-1
      endif
!
      jmtmax = max(jmtmax,Potential(id)%Grid%jmt)
      jwsmax = max(jwsmax,Potential(id)%Grid%jend)
      ig = getGlobalIndex(id)
      do ia = 1, Potential(id)%NumSpecies
         atname = getAlloyElementName(ig,ia)
         numc = getNumCoreStates(atname)
         SiteMaxNumc(id) = max(SiteMaxNumc(id),numc)
      enddo
      iend=max(iend,jend)
      numcmax = max(numcmax, SiteMaxNumc(id))
   enddo
!
   if ( .not.isDataStorageExisting('OldEstimatedCoreEnergy') ) then
      do id = 1,LocalNumAtoms
         DataSize(id) = SiteMaxNumc(id)*n_spin_pola*Potential(id)%NumSpecies
         if (DataSize(id) < 1) then      ! in case of empty cell
            DataSize(id) = n_spin_pola
         endif
      enddo
!     ---------------------------------------------------------------
      call createDataStorage(LocalNumAtoms,'OldEstimatedCoreEnergy', &
                             DataSize,RealType)
      call createDataStorage(LocalNumAtoms,'OldEstimatedCoreStates', &
                             DataSize,IntegerType)
!     ---------------------------------------------------------------
   endif
!
   msg(1) = jmtmax
   msg(2) = jwsmax
   msg(3) = numcmax
   msg(4) = lmax_max
   msg(5) = MaxNumSpecies
!  ------------------------------------------------------------------
   call GlobalMaxInGroup(GroupID,msg,5)
!  ------------------------------------------------------------------
   jmtmax = msg(1)
   jwsmax = msg(2)
   numcmax = msg(3)
   glb_lmax_max = msg(4)
   MaxNumSpecies = msg(5)
!
!  ==================================================================
!  determine the data size to be written for orbital dependent potentials
!  from LDA+U calculation.
!  ==================================================================
!
   allocate(LdaPlusU_DataAccum(GlobalNumAtoms))
   LdaPlusU_DataAccum(1:GlobalNumAtoms) = 0
   if ( checkLdaCorrection() ) then
      do id = 1, LocalNumAtoms
         ig = getGlobalIndex(id)
         if (checkLdaCorrection(id,1)) then ! temporary fix.....
            LdaPlusU_DataAccum(ig) = getDataPackSize(id)*NumSpecies(id)
         endif
      enddo
!     ---------------------------------------------------------------
      call GlobalSumInGroup(GroupID,LdaPlusU_DataAccum,GlobalNumAtoms)
!     ---------------------------------------------------------------
      do id = 2, GlobalNumAtoms
         LdaPlusU_DataAccum(id) = LdaPlusU_DataAccum(id) +           &
                                  LdaPlusU_DataAccum(id-1)
      enddo
   endif
!
   allocate(NonSphPot_DataAccum(GlobalNumAtoms))
   NonSphPot_DataAccum(1:GlobalNumAtoms) = 0
   if ( isFullPotential() ) then
      do id = 1, LocalNumAtoms
         Grid => getGrid(id)
         ig = getGlobalIndex(id)
!!!!!!      Here factor 2 is for complex type taking twice as much space
!!!!!!      as the real type
         NonSphPot_DataAccum(ig) = 2*Grid%jend*Potential(id)%jmax*   &
                                   NumSpecies(id)*n_spin_pola
      enddo
!     ---------------------------------------------------------------
      call GlobalSumInGroup(GroupID,NonSphPot_DataAccum,GlobalNumAtoms)
!     ---------------------------------------------------------------
      do id = 2, GlobalNumAtoms
         NonSphPot_DataAccum(id) = NonSphPot_DataAccum(id) +         &
                                   NonSphPot_DataAccum(id-1)
      enddo
   endif
!
   deallocate( DataSize )
!
   Initialized = .true.
   efermi=100.0d0
   vdif(1) = ZERO
!
   end subroutine initPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endPotential()
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: id
!
   do id = 1, LocalNumAtoms
      deallocate(Potential(id)%PotCompFlag)
      deallocate(Potential(id)%jtitle)
      deallocate(Potential(id)%ztss)
      nullify(Potential(id)%PotCompFlag)
      nullify(Potential(id)%pot_l)
      nullify(Potential(id)%potr_sph)
      if ( isTruncatedPotential ) then
         deallocate(Potential(id)%PotCompFlag_trunc)
         nullify(Potential(id)%pot_l_trunc)
         nullify(Potential(id)%potr_sph_trunc)
         Potential(id)%isTruncDone = .false.
      endif
      nullify(Potential(id)%Grid)
   enddo
   deallocate(Potential, print_level, SiteMaxNumc)
   deallocate( lofj, mofj, kofj, lofk, mofk, jofk, m1m )
   deallocate( efermi_in )
   deallocate( NonSphPot_DataAccum )
   deallocate( LdaPlusU_DataAccum)
   deallocate( NumSpecies )
   isChargeSymmOn = .false.
!
   end subroutine endPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSphPotr(id,ia,is) result(potr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), pointer :: potr(:)
!
   if ( .not.Initialized ) then
      call ErrorHandler("getSphrPotr",                                &
                        "The PotentialModule has to be initialized first")
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getSphPotr','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getSphPotr','invalid species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getSphPotr','invalid spin index',is)
   endif
!
   iend = Potential(id)%Grid%jend
   potr => Potential(id)%potr_sph(1:iend,is,ia)
!
   end function getSphPotr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSphPotr(id,ia,is,potr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: ir, iend
!
   real (kind=RealKind), intent(in) :: potr(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   iend = size(potr)
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setSphPotr','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('setSphPotr','invalid species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('setSphPotr','invalid spin index',is)
   else if ( iend > Potential(id)%Grid%jend .or.               &
             iend < Potential(id)%Grid%jmt ) then
      if (print_level(id) >= 0) then
         call WarningHandler('setSphPotr',                                 &
              'Input potential size differs from the internal array size', &
              iend,Potential(id)%Grid%jend)
      endif
      iend = min(iend,Potential(id)%Grid%jend)
   endif
   r_mesh => Potential(id)%Grid%r_mesh(1:iend)
   Potential(id)%potr_sph(1:iend,is,ia)=potr(1:iend)
   do ir=1,Potential(id)%Grid%jend
      Potential(id)%pot_l(ir,1,is,ia)=cmplx( potr(ir)/(Y0*r_mesh(ir)),&
                                          ZERO, kind=CmplxKind )
   enddo
!
   Potential(id)%PotCompFlag(1)=0
   LOOP_ir: do ir = 1,iend
      if (abs(potr(ir)) > TEN2m8) then
         Potential(id)%PotCompFlag(1)=1
         exit LOOP_ir
      endif
   enddo LOOP_ir
!
   if ( isTruncatedPotential ) then
      Potential(id)%isTruncDone = .false.
   endif
!
   end subroutine setSphPotr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isPotComponentZero(id,jl) result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: jl
!
   logical :: flag
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('isPotComponentZero','invalid id',id)
   else if (jl < 1 .or. jl > Potential(id)%jmax) then
      flag = .true.
   else if (Potential(id)%PotCompFlag(jl) == 0) then
      flag = .true.
   else
      flag = .false.
   endif
!
   end function isPotComponentZero
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSphericalInputFile()                       result(isSIF)
!  ===================================================================
   implicit none
!
   logical :: isSIF
!
   if ( SphericalInputFile ) then
      isSIF = .true.
   else
      isSIF = .false.
   endif
!
   end function isSphericalInputFile
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotComponentFlag(id)                    result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   integer, pointer :: flag(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPotComponentFlag','invalid id',id)
   endif
!
   flag => Potential(id)%PotCompFlag
!
   end function getPotComponentFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTruncatedPotComponentFlag(id)           result(flag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: ia
!
   integer, pointer :: flag(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTruncatedPotComponentFlag','invalid id',id)
   else if ( .not.isTruncatedPotential ) then
      call ErrorHandler('getTruncatedPotComponentFlag',               &
                        'Truncation have not been specified!')
   endif
!
   do ia = 1, NumSpecies(id)
      if ( .not.Potential(id)%isTruncDone(ia) ) then
         call truncatePotential(id,ia)
         Potential(id)%isTruncDone(ia) = .true.
      endif
   enddo
   flag => Potential(id)%PotCompFlag_trunc
!
   end function getTruncatedPotComponentFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPoten_l0(id,ia,is) result(pot_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: iend, jmax
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPoten_l0','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getPoten_l0','invalid species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getPoten_l0','invalid is',is)
   endif
!
   jmax = Potential(id)%jmax
   iend = Potential(id)%Grid%jend
   pot_l => Potential(id)%pot_l(1:iend,1:jmax,is,ia)
!
   end function getPoten_l0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPoten_l1(id,ia,is,jl) result(pot_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind) :: iend
!
   complex (kind=CmplxKind), pointer :: pot_l(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPoten_l1','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getPoten_l1','invalid species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getPoten_l1','invalid is',is)
   else if (jl < 1 .or. jl > Potential(id)%jmax) then
      call ErrorHandler('getPoten_l1','invalid jl',jl)
   endif
!
   iend = Potential(id)%Grid%jend
   pot_l => Potential(id)%pot_l(1:iend,jl,is,ia)
!
   end function getPoten_l1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPoten_l2(id,ia,posi,is) result(pot)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   use PotentialTypeModule, only : isMuffinTinPotential
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
   integer (kind=IntKind), optional :: is
!
   real (kind=RealKind), intent(in) :: posi(3)
!
   integer (kind=IntKind) :: iend, irp, ir, l, jl, kl, ns, jmt
!
   real (kind=RealKind) :: r, x, pot, err
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: pot_in
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
   complex (kind=CmplxKind), allocatable :: ylm(:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPoten_l','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getPoten_l2','invalid species index',ia)
   else if ( present(is) .and. ( is < 1 .or. is > n_spin_pola ) ) then
      call ErrorHandler('getPoten_l','invalid is',is)
   endif
!
   r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
   if (r > Potential(id)%Grid%rend) then
      pot = ZERO
      return
   endif
   x = log(r)
!
   iend = Potential(id)%Grid%jend
   r_mesh => Potential(id)%Grid%r_mesh(1:iend)
!  -------------------------------------------------------------------
   call hunt(iend,r_mesh(1:iend),r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   if ( isMuffinTinPotential() .and. Potential(id)%jmax==1 ) then
      jmt = Potential(id)%Grid%jmt
      if ( ir > jmt ) then
         pot = ZERO
         return
      else if (ir>jmt-(n_inter-1)/2 ) then
         irp = jmt-n_inter+1
      endif
   endif
!
   if ( r<=TEN2m6 ) then
      pot = -TWO*Potential(id)%ztss(ia)
      return
   endif
!
   allocate( ylm((Potential(id)%lmax+1)**2))
!  -------------------------------------------------------------------
   call calYlm(posi,Potential(id)%lmax,ylm)
!  -------------------------------------------------------------------
!
   pot = ZERO
   if ( present(is) ) then
      pot_l => Potential(id)%pot_l(1:iend,1:Potential(id)%jmax,is,ia)
!
      do jl = 1,Potential(id)%jmax
         l = lofj(jl)
!        -------------------------------------------------------------
         call PolyInterp(n_inter, r_mesh(irp:irp+n_inter-1),           &
                         pot_l(irp:irp+n_inter-1,jl), r, pot_in, err)
!        -------------------------------------------------------------
         kl = (l+1)*(l+1)-l+mofj(jl)
         if (mofj(jl) == 0) then
            pot = pot + real(pot_in*ylm(kl),RealKind)*r
         else
            pot = pot + TWO*real(pot_in*ylm(kl),RealKind)*r
         endif
      enddo
   else
      do ns = 1,n_spin_pola
         pot_l => Potential(id)%pot_l(1:iend,1:Potential(id)%jmax,ns,ia)
!
         do jl = 1,Potential(id)%jmax
            l = lofj(jl)
!           ----------------------------------------------------------
            call PolyInterp(n_inter, r_mesh(irp:irp+n_inter-1),           &
                         pot_l(irp:irp+n_inter-1,jl), r, pot_in, err)
!           ----------------------------------------------------------
            kl = (l+1)*(l+1)-l+mofj(jl)
            if (mofj(jl) == 0) then
               pot = pot + real(pot_in*ylm(kl),RealKind)*r
            else
               pot = pot + TWO*real(pot_in*ylm(kl),RealKind)*r
            endif
         enddo
      enddo
      pot = pot/(ONE+real(n_spin_pola-1,kind=RealKind))
   endif
!
   deallocate( ylm )
!
   end function getPoten_l2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPoten_l3(id,ia,jmax_in,posi,is) result(pot)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   use InterpolationModule, only : PolyInterp
!
   use PotentialTypeModule, only : isMuffinTinPotential
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), optional :: is
   integer (kind=IntKind), intent(in) :: jmax_in
!
   real (kind=RealKind), intent(in) :: posi(3)
!
   integer (kind=IntKind) :: iend, irp, ir, l, kl, jl, ns, jmt, jmax
!
   real (kind=RealKind) :: r, x, pot, err
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: pot_in
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
   complex (kind=CmplxKind), allocatable :: ylm(:)
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPoten_l3','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getPoten_l3','invalid species index',ia)
   else if ( jmax_in > Potential(id)%jmax .or. jmax_in<1 ) then
      call ErrorHandler('getPoten_l3','Invalid jmax',jmax_in)
   else if ( present(is) .and. (is < 1 .or. is > n_spin_pola) ) then
      call ErrorHandler('getPoten_l3','invalid is',is)
   endif
!
   if ( jmax_in<1 ) then
      jmax = Potential(id)%jmax
   else
      jmax = jmax_in
   endif
   pot = ZERO
   r = sqrt( posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
   if (r > Potential(id)%Grid%rend) then
      return
   endif
   x = log(r)
!
   iend = Potential(id)%Grid%jend
   r_mesh => Potential(id)%Grid%r_mesh(1:iend)
!  -------------------------------------------------------------------
   call hunt(iend,r_mesh(1:iend),r,ir)
!  -------------------------------------------------------------------
   if (ir > iend-(n_inter-1)/2) then
      irp=iend-n_inter+1
   else if (2*ir+1 > n_inter) then
      irp=ir-(n_inter-1)/2
   else
      irp=1
   endif
!
   if ( isMuffinTinPotential() .and. Potential(id)%jmax==1 ) then
      jmt = Potential(id)%Grid%jmt
      if ( ir > jmt ) then
         return
      else if ( ir>jmt-(n_inter-1)/2 ) then
         irp = jmt-n_inter+1
      endif
   endif
!
   if ( r<=TEN2m6 ) then
      pot = -TWO*Potential(id)%ztss(ia)
      return
   endif
!
   allocate( ylm((Potential(id)%lmax+1)**2))
!  -------------------------------------------------------------------
   call calYlm(posi,Potential(id)%lmax,ylm)
!  -------------------------------------------------------------------
!
   if ( present(is) ) then
      pot_l => Potential(id)%pot_l(1:iend,1:jmax,is,ia)
!
      do jl = 1,jmax
         l = lofj(jl)
!        -------------------------------------------------------------
         call PolyInterp(n_inter, r_mesh(irp:irp+n_inter-1),           &
                         pot_l(irp:irp+n_inter-1,jl), r, pot_in, err)
!        -------------------------------------------------------------
         kl = (l+1)*(l+1)-l+mofj(jl)
         if (mofj(jl) == 0) then
            pot = pot + real(pot_in*ylm(kl),RealKind)*r
         else
            pot = pot + TWO*real(pot_in*ylm(kl),RealKind)*r
         endif
      enddo
   else
      do ns = 1,n_spin_pola
         pot_l => Potential(id)%pot_l(1:iend,1:jmax,ns,ia)
!
         do jl = 1,jmax
            l = lofj(jl)
!           ----------------------------------------------------------
            call PolyInterp(n_inter, r_mesh(irp:irp+n_inter-1),           &
                         pot_l(irp:irp+n_inter-1,jl), r, pot_in, err)
!           ----------------------------------------------------------
            kl = (l+1)*(l+1)-l+mofj(jl)
            if (mofj(jl) == 0) then
               pot = pot + real(pot_in*ylm(kl),RealKind)*r
            else
               pot = pot + TWO*real(pot_in*ylm(kl),RealKind)*r
            endif
         enddo
      enddo
      pot = pot/(ONE+real(n_spin_pola-1,kind=RealKind))
   endif
!
   deallocate( ylm )
!
   end function getPoten_l3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPoten_r(id,ia,is,jl,ir) result(pot_r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind), intent(in) :: ir
   integer (kind=IntKind) :: iend
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: pot_r
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPoten_r','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getPoten_r','invalid species index',ia)
   else if (ir < 1) then
      call ErrorHandler('getPoten_r','ir < 1',ir)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getPoten_r','invalid is',is)
   else if (jl < 1 .or. jl > Potential(id)%jmax) then
      pot_r = CZERO
      return
   else if (Potential(id)%PotCompFlag(jl) == 0) then
      pot_r = CZERO
      return
   endif
!
   iend = Potential(id)%Grid%jend
   if (ir > iend) then
      call ErrorHandler('getPoten_r','ir < iend',ir)
   endif
!
   r_mesh=>Potential(id)%Grid%r_mesh(1:iend)
   if (jl == 1) then
      pot_r=Potential(id)%potr_sph(ir,is,ia)/(Y0*r_mesh(ir))
   else
      pot_r=Potential(id)%pot_l(ir,jl,is,ia)
   endif
!
   end function getPoten_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTruncPoten_l0(id,ia,is) result(pot_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: iend, jmax
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTruncPoten_l0','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getTruncPoten_l0','invalid species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getTruncPoten_l0','invalid is',is)
   else if ( .not.isTruncatedPotential ) then
      call ErrorHandler('getTruncPoten_l0','Truncation have not been specified!')
   endif
!
   if ( .not.Potential(id)%isTruncDone(ia) ) then
      call truncatePotential(id,ia)
      Potential(id)%isTruncDone(ia) = .true.
   endif
!
   jmax = Potential(id)%jmax_trunc
   iend = Potential(id)%iend_trunc
   pot_l => Potential(id)%pot_l_trunc(1:iend,1:jmax,is,ia)
!
   end function getTruncPoten_l0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTruncPoten_l1(id,ia,is,jl) result(pot_l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind) :: iend
!
   complex (kind=CmplxKind), pointer :: pot_l(:)
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTruncPoten_l1','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getTruncPoten_l1','invalid species index',ia)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getTruncPoten_l1','invalid is',is)
   else if ( .not.isTruncatedPotential ) then
      call ErrorHandler('getTruncPoten_l1','Truncation have not been specified!')
   else if (jl < 1 .or. jl > Potential(id)%jmax_trunc) then
      call ErrorHandler('getTruncPoten_l1','invalid jl',jl)
   endif
!
   if ( .not.Potential(id)%isTruncDone(ia) ) then
      call truncatePotential(id,ia)
      Potential(id)%isTruncDone(ia) = .true.
   endif
!
   iend = Potential(id)%iend_trunc
   pot_l => Potential(id)%pot_l_trunc(1:iend,jl,is,ia)
!
   end function getTruncPoten_l1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTruncPoten_r(id,ia,is,jl,ir) result(pot_r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind), intent(in) :: ir
   integer (kind=IntKind) :: iend, jmt
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: pot_r
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTruncPoten_r','invalid id',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getTruncPoten_r','invalid species index',ia)
   else if (ir < 1) then
      call ErrorHandler('getTruncPoten_r','ir < 1',ir)
   else if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('getTruncPoten_r','invalid is',is)
   else if ( .not.isTruncatedPotential ) then
      call ErrorHandler('getTruncPoten_r','Truncation have not been specified!')
   else if (jl < 1 .or. jl > Potential(id)%jmax_trunc) then
      pot_r = CZERO
      return
   else if (Potential(id)%PotCompFlag_trunc(jl) == 0) then
      pot_r = CZERO
      return
   endif
!
   if ( .not.Potential(id)%isTruncDone(ia) ) then
      call truncatePotential(id,ia)
      Potential(id)%isTruncDone(ia) = .true.
   endif
!
   iend = Potential(id)%iend_trunc
   if (ir > iend) then
      call ErrorHandler('getTruncPoten_r','ir < iend',ir)
   endif
!
   jmt = Potential(id)%Grid%jmt
   r_mesh=>Potential(id)%Grid%r_mesh(jmt:iend)
   if ( ir>jmt ) then
      if (jl == 1) then
         pot_r=Potential(id)%potr_sph_trunc(ir-jmt+1,is,ia)/(Y0*r_mesh(ir))
      else
         pot_r=Potential(id)%pot_l_trunc(ir-jmt+1,jl,is,ia)
      endif
   else
      if (jl == 1) then
         pot_r=Potential(id)%potr_sph(ir,is,ia)/(Y0*r_mesh(ir))
      else
         pot_r=Potential(id)%pot_l(ir,jl,is,ia)
      endif
   endif
!
   end function getTruncPoten_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPotential(id,ia,is,jl,pot_l,fl)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: jl
   integer (kind=IntKind), intent(in) :: fl
   integer (kind=IntKind) :: ir, iend
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), intent(in) :: pot_l(:)
!
   logical :: nonzero
!
   iend = Potential(id)%Grid%jend
!
   if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('setPotential','invalid is',is)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('setPotential','invalid species index',ia)
   else if (jl < 1 .or. jl > Potential(id)%jmax) then
      call ErrorHandler('setPotential','jl out of bound',jl)
   else if (iend > size(pot_l)) then
      call ErrorHandler('setPotential','iend is out of range',iend,size(pot_l))
   endif
!
   if (jl == 1) then
      r_mesh=>Potential(id)%Grid%r_mesh(1:iend)
      do ir = 1,iend
         Potential(id)%potr_sph(ir,is,ia)=Y0*pot_l(ir)*r_mesh(ir)
      enddo
   endif
!
!  check the integraty of the potential ...
!
   nonzero = .false.
   LOOP_ir: do ir = 1, iend
      if (abs(pot_l(ir)) > pot_tol) then
         nonzero = .true.
         exit LOOP_ir
      endif
   enddo LOOP_ir
   if (nonzero .and. fl == 0) then
      call WarningHandler('setPotential',                          &
                          'The potential component is forced to be ZERO',jl)
      Potential(id)%PotCompFlag(jl) = 0
      do ir = 1,iend
         Potential(id)%pot_l(ir,jl,is,ia) = CZERO
      enddo
   else if (.not.nonzero .and. fl > 0) then
      call WarningHandler('setPotential',                          &
                          'The potential component is found to be ZERO',jl)
      Potential(id)%PotCompFlag(jl) = 0
      do ir = 1,iend
         Potential(id)%pot_l(ir,jl,is,ia) = CZERO
      enddo
   else
      Potential(id)%PotCompFlag(jl) = fl
      do ir = 1,iend
         Potential(id)%pot_l(ir,jl,is,ia) = pot_l(ir)
      enddo
   endif
!
   if ( isTruncatedPotential ) then
      Potential(id)%isTruncDone = .false.
   endif
!
   end subroutine setPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPotentialOutsideMT(id,ia,is,v0,nfit)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in), optional :: nfit
!
   real (kind=RealKind), intent(in), optional :: v0
!
   integer (kind=IntKind) :: ir, iend, jmt, jl, jl_max, iex
!
   real (kind=RealKind), pointer :: r_mesh(:), x_mesh(:)
!
   interface
      function ylag(xi,x,y,ind,n,imax,iex) result(lag)
         use KindParamModule, only : IntKind, RealKind
         integer (kind=IntKind), intent(in) :: ind, n, imax
         integer (kind=IntKind), intent(out) :: iex
         real (kind=RealKind), intent(in) :: xi, x(*), y(*)
         real (kind=RealKind) :: lag
      end function
   end interface
!
   interface
      function cylag(xi,x,y,ind,n,imax,iex) result(lag)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         integer (kind=IntKind), intent(in) :: ind, n, imax
         integer (kind=IntKind), intent(out) :: iex
         real (kind=RealKind), intent(in) :: xi, x(*)
         complex (kind=CmplxKind), intent(in) :: y(*)
         complex (kind=CmplxKind) :: lag
      end function
   end interface
!
   if (is < 1 .or. is > n_spin_pola) then
      call ErrorHandler('setPotential','invalid is',is)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('setPotential','invalid species index',ia)
   endif
!
   jl_max = Potential(id)%jmax
   iend = Potential(id)%Grid%jend
   jmt = Potential(id)%Grid%jmt
   if ( present(nfit) ) then
      x_mesh=>Potential(id)%Grid%x_mesh(1:iend)
      do ir = jmt+1,iend
         Potential(id)%potr_sph(ir,is,ia) =                           &
             ylag(x_mesh(ir),x_mesh,Potential(id)%potr_sph(:,is,ia),0,nfit,jmt,iex)
         Potential(id)%pot_l(ir,1,is,ia) =                            &
             cylag(x_mesh(ir),x_mesh,Potential(id)%pot_l(:,1,is,ia),0,nfit,jmt,iex)
      enddo
   else if ( present(v0) ) then
      r_mesh=>Potential(id)%Grid%r_mesh(1:iend)
      do ir = jmt+1,iend
         Potential(id)%potr_sph(ir,is,ia) = v0*r_mesh(ir)
         Potential(id)%pot_l(ir,1,is,ia) = v0/Y0
      enddo
   else
      do ir = jmt+1,iend
         Potential(id)%potr_sph(ir,is,ia) = Potential(id)%potr_sph(jmt,is,ia)
         Potential(id)%pot_l(ir,1,is,ia) = Potential(id)%pot_l(jmt,1,is,ia)
      enddo
   endif
!  do jl = 2,jl_max
!     Potential(id)%pot_l(1:iend,jl,is,ia) = ZERO
!  enddo
!  -------------------------------------------------------------------
   call setPotCompFlag(id,is,ia)
!  -------------------------------------------------------------------
   if ( isTruncatedPotential ) then
      Potential(id)%isTruncDone = .false.
   endif
!
   end subroutine setPotentialOutsideMT
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getHeader(id) result(hd)
!  ===================================================================
   implicit none
!
   character (len=80) :: hd
   integer (kind=IntKind), intent(in) :: id
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getHeader','Invalid atom index',id)
   endif
   hd = Potential(id)%header
   end function getHeader
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setHeader(id,hd)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: hd
   integer (kind=IntKind), intent(in) :: id
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getHeader','Invalid atom index',id)
   endif
   Potential(id)%header = hd
   end subroutine setHeader
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTitle(id,ia,is) result(tl)
!  ===================================================================
   implicit none
!
   character (len=80) :: tl
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTitle','Invalid atom index',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('getTitle','invalid species index',ia)
   else if (is > n_spin_pola .or. is < 1) then
      call ErrorHandler('getTitle','Invalid spin index',is)
   endif
   tl = Potential(id)%jtitle(is,ia)
   end function getTitle
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setTitle(id,ia,is,tl)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: tl
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ia
   integer (kind=IntKind), intent(in) :: is
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('setTitle','Invalid atom index',id)
   else if (ia < 1 .or. ia > NumSpecies(id)) then
      call ErrorHandler('setTitle','invalid species index',ia)
   else if (is > n_spin_pola .or. is < 1) then
      call ErrorHandler('setTitle','Invalid spin index',is)
   endif
   Potential(id)%jtitle(is,ia) = tl
!
   end subroutine setTitle
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setV0(is,vzero)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: iend, id, jmt, ia
!
   real (kind=RealKind), intent(in) :: vzero
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   if (is > n_spin_pola .or. is < 1) then
      call ErrorHandler('setV0','Invalid spin index',is)
   endif
!
   do id = 1, LocalNumAtoms
      iend = Potential(id)%Grid%jend
      if (abs(vzero) > TEN2m8) then
         r_mesh=>Potential(id)%Grid%r_mesh(1:iend)
         jmt = Potential(id)%Grid%jmt
         do ia = 1, NumSpecies(id)
            Potential(id)%potr_sph(1:jmt,is,ia) = Potential(id)%potr_sph(1:jmt,is,ia) &
                                                 -vzero*r_mesh(1:jmt)
            Potential(id)%pot_l(1:jmt,1,is,ia) = Potential(id)%pot_l(1:jmt,1,is,ia)   &
                                                 -vzero/Y0
          enddo
      endif
   enddo
   v0(is) = ZERO
!
   end subroutine setV0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getV0(is) result(vzero)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind) :: vzero
!
   if (is > n_spin_pola .or. is < 1) then
      call ErrorHandler('getV0','Invalid spin index',is)
   endif
   vzero = v0(is)
!
   end function getV0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setVdif(vd)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: vd
!
   vdif(1) = vd
!
   end subroutine setVdif
!  ===================================================================
!
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVdif() result(vd)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), pointer :: vd(:)
!
   vd => vdif
!
   end function getVdif
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotEf() result(ef)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: ef
!
   ef = efermi
!
   end function getPotEf
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPotEf(ef)
!  ===================================================================
   implicit none
   real (kind=RealKind), intent(in) :: ef
!
   efermi = ef
!
   end subroutine setPotEf
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotLmax(id) result(lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: lmax
!
   lmax = Potential(id)%lmax
   end function getPotLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTruncPotLmax(id) result(lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: lmax
!
   lmax = Potential(id)%lmax_trunc
   end function getTruncPotLmax
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readPotential()
!  ===================================================================
   use GroupCommModule, only : GlobalMinInGroup
!
   use ChemElementModule, only : getZtot
!
   use SystemModule, only : getNumAlloyElements
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use AtomModule, only : getInPotFileName, getInPotFileForm
!
   use ParallelIOModule, only : isInputProc
!
   use PotentialTypeModule, only : isFullPotential
!
   implicit none
!
   logical :: isChgSymm_save
!  character (len=10) :: acc
   character (len=13), parameter :: sname='readPotential'
   character (len=10), parameter :: file_status = 'READ'
   character (len=11) :: file_form
!
   integer (kind=IntKind) :: id, ia, is, vunit, str_len, iend, ir
   integer (kind=IntKind) :: n, ig
!
   real (kind=RealKind) :: pot_shift
!
   file_form = adjustl(getInPotFileForm())
!
   if (trim(file_form) == 'FORMATTED') then
!     if (getInPotFileForm(id) == 'FORMATTED') then
!        acc = 'SEQUENTIAL'
!     else
!        acc = 'DIRECT'
!     endif
!     open(unit=funit,file=filename,form=fileform,access=acc)
!     isChgSymm_save = isChargeSymmOn
!     isChargeSymmOn = .false.
      do id =1, LocalNumAtoms
         Potential(id)%pot_l = CZERO
!        -------------------------------------------------------------
         call readFormattedData(id)
!        -------------------------------------------------------------
         do ia = 1, NumSpecies(id)
            do is = 1, n_spin_pola
!              -------------------------------------------------------
               call setPotCompFlag(id,is,ia)
!              -------------------------------------------------------
            enddo
         enddo
      enddo
!     isChargeSymmOn = isChgSymm_save
   else if (trim(file_form) == 'UNFORMATTED' .or. trim(file_form) == 'XDR') then
!     isChgSymm_save = isChargeSymmOn
!     isChargeSymmOn = .false.
      if( isInputProc() ) then
         str_len = len_trim(adjustl(getInPotFileName()))
!        -------------------------------------------------------------
         call c_gopen(vunit,trim(adjustl(getInPotFileName())),str_len, &
                      trim(file_status),len_trim(file_status),         &
                      trim(file_form),len_trim(file_form))
!        -------------------------------------------------------------
      endif
      do id =1, LocalNumAtoms
         Potential(id)%pot_l = CZERO
!        -------------------------------------------------------------
         call readUnformattedData(vunit,id)
!        -------------------------------------------------------------
         do ia = 1, NumSpecies(id)
            do is = 1, n_spin_pola
!              -------------------------------------------------------
               call setPotCompFlag(id,is,ia)
!              -------------------------------------------------------
            enddo
         enddo
      enddo
      if( isInputProc() ) then
         call c_close(vunit)
      endif
!     isChargeSymmOn = isChgSymm_save
   else if (trim(file_form) == 'HDF') then
      print *,'HDF file I/O needs to be implemented'
      stop 'incomplete'
   else
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'Unknown potential data format',        &
                        getInPotFileForm())
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   call GlobalMinInGroup(GroupID,efermi)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  The starting potential may come with different fermi energy.
!  Shift the potential accordingly.
!  ===================================================================
#ifdef POT_DIPOL
!   do id = 1, LocalNumAtoms
!      jmt = Potential(id)%Grid%jmt
!      do is = 1, n_spin_pola
!         if ( Potential(id)%jmax >= 2 ) then
!            do ia = 1, NumSpecies(id)
!               Potential(id)%pot_l(1:Potential(id)%Grid%jend, &
!                                   2:Potential(id)%jmax,is,ia) = CZERO
!               do ir = 1,Potential(id)%Grid%jend
!                  Potential(id)%pot_l(ir,2,is,ia) =  &
!                                   0.05d0*(Potential(id)%Grid%r_mesh(ir))
!                  Potential(id)%pot_l(1:Potential(id)%Grid%jend, &
!                                      11,is,ia) = 0.5d0*(Potential(id)%Grid%r_mesh(ir))**4
!               enddo
!            enddo
!            Potential(id)%PotCompFlag(2,1) = 1
!            Potential(id)%PotCompFlag(2,n_spin_pola) = 1
!         endif
!      enddo
!   enddo
#endif
   do id = 1, LocalNumAtoms
      do ia = 1, NumSpecies(id)
         pot_shift = efermi_in(ia,id) - efermi
         if (abs(pot_shift) > TEN2m6) then
            if ( isFullPotential() ) then
               iend = Potential(id)%Grid%jend
            else
               iend = Potential(id)%Grid%jmt
            endif
            do is = 1, n_spin_pola
               do ir = 1,iend
                  Potential(id)%pot_l(ir,1,is,ia)  = Potential(id)%pot_l(ir,1,is,ia) -   &
                                                     pot_shift/Y0
                  Potential(id)%potr_sph(ir,is,ia) = Potential(id)%pot_l(ir,1,is,ia)*Y0* &
                                                     Potential(id)%Grid%r_mesh(ir)
               enddo
            enddo
         endif
      enddo
   enddo
!
#ifdef TMP_ACCEL
      call pushPotentialToAccel()
      call pushTruncatedPotToAccel()
#endif
   end subroutine readPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readFormattedData(id)
!  ===================================================================
   use MathParamModule, only : ZERO
!
   use PublicParamDefinitionsModule, only : MaxLenFileName
!
   use InterpolationModule, only : PolyInterp, getInterpolation
!
   use MPPModule, only : MyPE
!
   use SystemModule, only : setLatticeConstant, getAtomName, getAlloyElementName
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : getZcor, getZsem, getZtot
   use ChemElementModule, only : getCoreStateN, getCoreStateL
   use ChemElementModule, only : getCoreStateKappa
   use ChemElementModule, only : getNumCoreStates, setNumCoreStates
   use ChemElementModule, only : getCoreStateIndex
!
   use RadialGridModule, only : getNumRmesh
!
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       setDataStorage2Value,  &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark, &
                                       IntegerType, IntegerMark
!
   use AtomModule, only : getInPotFileName
!
   implicit none
!
   logical :: nmd
!
   character (len=MaxLenFileName) :: filename
!
   character (len=1) :: dummy
   character (len=MaxLenOfAtomName) :: atname
   character (len=MaxLenFileName) :: inp
   character (len=17), parameter :: sname='readFormattedData'
   character (len=10), parameter :: acc = 'SEQUENTIAL'
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind), parameter :: funit=91
!
   integer (kind=IntKind) :: is, ig, js, n, ia
   integer (kind=IntKind) :: j_inter, irp
   integer (kind=IntKind) :: ns,j,jmt,nrrho,nrcor,jmax
   integer (kind=IntKind) :: numc,jz,jc,ios,lmax,jend,jl,nr,ir,jm,jlr
   integer (kind=IntKind) :: nc, lc, kc
   integer (kind=IntKind) :: DataSize(LocalNumAtoms)
   integer (kind=IntKind), pointer :: pc0(:,:,:)
!
   real (kind=RealKind) :: ztss,alat,zcss,edum
   real (kind=RealKind) :: xstart,xmt,vzero(2),xvalws(2)
   real (kind=RealKind) :: err,h,hin,hout,x
   real (kind=RealKind), allocatable :: vr(:,:), x_mesh(:), r_mesh(:)
   real (kind=RealKind), allocatable :: rhoin(:,:), rhotot(:)
   real (kind=RealKind), pointer :: rho0(:,:)
   real (kind=RealKind), pointer :: mom0(:,:)
   real (kind=RealKind), pointer :: ec(:,:,:)
!
   complex (kind=CmplxKind) :: vc
   complex (kind=CmplxKind), allocatable :: pot_l(:)
!
   type (PotentialStruct), pointer :: ThisP
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   ThisP=>Potential(id)
!
   if (SiteMaxNumc(id) >= 1) then
      ec => getDataStorage(id,'OldEstimatedCoreEnergy', SiteMaxNumc(id), n_spin_pola, &
                           NumSpecies(id), RealMark)
      pc0 => getDataStorage( id, 'OldEstimatedCoreStates',            &
                             SiteMaxNumc(id), n_spin_pola, NumSpecies(id), IntegerMark)
   endif
!
   nr = getNumRmesh(id)
   allocate(rhotot(nr))
   rho0 => getDataStorage( id, 'OldSphericalElectronDensity',         &
                           nr+1, NumSpecies(id), RealMark )
   if (n_spin_pola == 2) then
      mom0 => getDataStorage(id,'OldSphericalMomentDensity',          &
                             nr+1, NumSpecies(id), RealMark )
   endif
!
   efermi = 100.0d0
   do ia = 1, NumSpecies(id)
      filename = getInPotFileName(id,ia)
!     ----------------------------------------------------------------
      open(unit=funit,file=filename,form='formatted',access=acc)
!     ----------------------------------------------------------------
!
      read(funit,'(a)') ThisP%header
      read(funit,'(i5,3x,d20.13)') ns,vdif(1)
      if (ns /= n_spin_pola .and. MyPE == 0) then
         inquire(unit=funit,named=nmd,name=inp)
         write(6,'(''File Name: '',a)')trim(inp)
         call WarningHandler(sname,'Spin in file <> n_spin_pola',ns,n_spin_pola)
      endif
      do is=1,ns
         js = min(is,n_spin_pola)
         read(funit,'(a)')ThisP%jtitle(js,ia)
!        read(funit,'(f5.0,17x,f12.5,f5.0,e20.13)')ztss,alat,zcss,efermi_in(ia,id)
         read(funit,*)ztss,alat,zcss,efermi_in(ia,id)
         jz=ztss
         ThisP%ztss(ia) = ztss
         jc=zcss
         ig = getGlobalIndex(id)
         atname = getAlloyElementName(ig,ia)
         if(jz /= getZtot(atname)) then
            call ErrorHandler(sname,'Inconsistent atom type',ztss)
         else if (jc /= getZcor(atname)) then
            zcss = getZcor(atname)
            jc=zcss
            call WarningHandler(sname,'Deep core charge is changed to',zcss)
         endif
         read(funit,'(17x,2e20.13,i5)') xstart,xmt,jmt
!
         if (is == 1) then
            allocate(vr(jmt,max(ns,n_spin_pola)))
         endif
!        =============================================================
!        read in formatted one-electron LDA potential..............
!        =============================================================
         read(funit,'(4e20.13)') (vr(j,is),j=1,jmt)
         read(funit,'(35x,e20.13)') vzero(is)
!
!        =============================================================
         read(funit,'(i5,e20.13)') nrrho,xvalws(is)
!        =============================================================
!
         if (is == 1) then
            allocate(rhoin(nrrho,max(ns,n_spin_pola)))
         endif
!        =============================================================
!        read in formatted total charge density....................
!        =============================================================
         read(funit,'(4e20.13)') (rhoin(j,is),j=1,nrrho)
!
!        =============================================================
!        read in formatted core state information..................
!        =============================================================
         read(funit,'(2i5)') numc,nrcor
         do j=1,numc
            if (j <= getNumCoreStates(atname)) then
               read(funit,*)nc,lc,kc,ec(j,js,ia)
               pc0(j,js,ia) = getCoreStateIndex(nc,lc,kc)
!              if (nc /= getCoreStateN(atname,j)) then
!                 call ErrorHandler(sname,'nc <> getCoreStateN()',nc, &
!                                   getCoreStateN(atname,j))
!              else if (lc /= getCoreStateL(atname,j)) then
!                 call ErrorHandler(sname,'lc <> getCoreStateL()',lc, &
!                                   getCoreStateL(atname,j))
!              else if (kc /= getCoreStateKappa(atname,j)) then
!                 call ErrorHandler(sname,'kc <> getCoreStateKappa()',kc, &
!                                   getCoreStateKappa(atname,j))
!              endif
            else
               read(funit,*)nc,lc,kc,edum
            endif
         enddo
!
         if (numc > getNumCoreStates(atname)) then
            if (print_level(id) >= 0) then
               call WarningHandler(sname,'numc > getNumCoreStates()',numc,&
                                   getNumCoreStates(atname))
            endif
            call setNumCoreStates(atname,numc)
         else if (numc < getNumCoreStates(atname) ) then
            if (print_level(id) >= 0) then
               call WarningHandler(sname,'numc < getNumCoreStates', numc, &
                                   getNumCoreStates(atname))
            endif
            call setNumCoreStates(atname,numc)
         endif
!
         if (nrcor > 0) then
            nrcor=ceiling(nrcor/4.0)       ! data is stored 4 numbers per line
            read(funit,'(a)')(dummy,j=1,nrcor)
         endif
      enddo
!
      if (ns < n_spin_pola) then
         xvalws(1) = HALF*xvalws(1)
         xvalws(2) = xvalws(1)
         do j = 1, jmt
            vr(j,2) = vr(j,1)
         enddo
         vzero(2) = vzero(1)
         do j=1,numc
            if (j <= getNumCoreStates(atname)) then
               ec(j,2,ia) = ec(j,1,ia)
               pc0(j,2,ia)=pc0(j,1,ia)
            endif
         enddo
         do j = 1, nrrho
            rhoin(j,1) = HALF*rhoin(j,1)
         enddo
         do j = 1, nrrho
            rhoin(j,2) = rhoin(j,1)
         enddo
      else if (ns > n_spin_pola) then
         xvalws(1) = xvalws(1)+ xvalws(2)
         do j = 1, jmt
            vr(j,1) = HALF*(vr(j,1)+vr(j,2))
         enddo
         do j = 1, nrrho
            rhoin(j,1) = rhoin(j,1)+rhoin(j,2)
         enddo
      endif
!
!     ================================================================
!     generate grid for potential read-in ......................
!     ================================================================
      nrrho=max(nrrho,jmt)
!
      allocate(x_mesh(nrrho), r_mesh(nrrho))
!
      h=(xmt-xstart)/dble(jmt-1)
      do j=1,nrrho
         x_mesh(j)=xstart+(j-1)*h
      enddo
      do j=1,nrrho
         r_mesh(j)=exp(x_mesh(j))
      enddo
!
      do is = 1, n_spin_pola
         do j=1,jmt
            vr(j,is)=vr(j,is)-vzero(is)*r_mesh(j)
         enddo
!
!        =============================================================
!        interpolate, potential and chg density onto GRID mesh
!        =============================================================
         if (ThisP%Grid%jmt < n_inter) then
!           ----------------------------------------------------------
            call ErrorHandler(sname,'jmt < n_inter',ThisP%Grid%jmt,n_inter)
!           ----------------------------------------------------------
         endif
         j_inter=1
         do j=1,ThisP%Grid%jmt
!           ----------------------------------------------------------
   !!       call hunt(jmt,x_mesh,ThisP%Grid%x_mesh(j),j_inter)
!           ----------------------------------------------------------
   !!       if (j_inter > jmt-(n_inter-1)/2) then
   !!          irp=jmt-n_inter+1
   !!       else if (2*j_inter+1 > n_inter) then
   !!          irp=j_inter-(n_inter-1)/2
   !!       else
   !!          irp=1
   !!       endif
!           ----------------------------------------------------------
   !!       call PolyInterp(n_inter,x_mesh(irp:irp+n_inter-1),           &
   !!                       vr(irp:irp+n_inter-1,is),                    &
   !!                       ThisP%Grid%x_mesh(j),ThisP%potr_sph(j,is,ia),err)
   !!       call interp(r_mesh(1:jmt),vr(1:jmt,is),jmt,ThisP%Grid%r_mesh(j), &
   !!                   ThisP%potr_sph(j,is,ia),err)
            ThisP%potr_sph(j,is,ia) =                                    &
                        getInterpolation(jmt,x_mesh(1:jmt),vr(1:jmt,is), &
                                         ThisP%Grid%x_mesh(j),err)
            ThisP%pot_l(j,1,is,ia)=ThisP%potr_sph(j,is,ia)/(Y0*ThisP%Grid%r_mesh(j))
!           ----------------------------------------------------------
         enddo
!
         do j=1,nr
!           ----------------------------------------------------------
  !!        call hunt(nrrho,x_mesh,ThisP%Grid%x_mesh(j),j_inter)
!           ----------------------------------------------------------
  !!        if (j_inter > nrrho-(n_inter-1)/2) then
  !!           irp=nrrho-n_inter+1
  !!        else if (2*j_inter+1 > n_inter) then
  !!           irp=j_inter-(n_inter-1)/2
  !!        else
  !!           irp=1
  !!        endif
!           ----------------------------------------------------------
  !!        call PolyInterp(n_inter,x_mesh(irp:irp+n_inter-1),        &
  !!                        rhoin(irp:irp+n_inter-1,is),              &
  !!                        ThisP%Grid%x_mesh(j),rhotot(j),err)
!           ----------------------------------------------------------
            rhotot(j) =                                               &
                  getInterpolation(nrrho,x_mesh(1:nrrho),             &
                                   rhoin(1:nrrho,is),ThisP%Grid%x_mesh(j),err)
         enddo
         if (is == 1) then
            do ir = 1,nr
               rho0(ir,ia) = rhotot(ir)/                              &
                            (PI4*ThisP%Grid%r_mesh(ir)*ThisP%Grid%r_mesh(ir))
            enddo
            rho0(nr+1,ia) = xvalws(is)
         else
            do ir = 1,nr
               rho0(ir,ia) = rho0(ir,ia)+rhotot(ir)/                  &
                            (PI4*ThisP%Grid%r_mesh(ir)*ThisP%Grid%r_mesh(ir))
            enddo
            rho0(nr+1,ia) = rho0(nr+1,ia) + xvalws(is)
         endif
         if ( n_spin_pola==2 ) then
            if (is == 1) then
               do ir = 1,nr
                  mom0(ir,ia) = rhotot(ir)/                           &
                            (PI4*ThisP%Grid%r_mesh(ir)*ThisP%Grid%r_mesh(ir))
               enddo
               mom0(nr+1,ia) = xvalws(is)
            else
               do ir = 1,nr
                  mom0(ir,ia) = mom0(ir,ia) - rhotot(ir)/             &
                            (PI4*ThisP%Grid%r_mesh(ir)*ThisP%Grid%r_mesh(ir))
               enddo
               mom0(nr+1,ia) = mom0(nr+1,ia) - xvalws(is)
            endif
         endif
      enddo
!
      deallocate(vr,rhoin,x_mesh,r_mesh)
!
      v0(1:2) = ZERO
      call setLatticeConstant(ThisP%Grid%rmt)
!
      read(funit,'(5i8)',iostat=ios)lmax,jm,jmt,jend
      if (ios == 0 .and. lmax > 0) then
         jmax=min(jm,ThisP%jmax)
         allocate(x_mesh(jend), r_mesh(jend), pot_l(jend))
         pot_l = CZERO
         read(funit,'(4d20.13)')xstart,xmt,hin,hout
         do ir=1,jmt
            x_mesh(ir)=xstart+(ir-1)*hin
         enddo
         do ir=jmt+1,jend
            x_mesh(ir)=xmt+(ir-jmt)*hout
         enddo
         do ir=1,jend
            r_mesh(ir)=exp(x_mesh(ir))
         enddo
         do is=1,n_spin_pola
            do jl = 1, jmax
               read(funit,'(2i5)')jlr,j
               if (jl /= jlr) then
                  call ErrorHandler('readFormattedData','jl <> jlr',jl,jlr)
               endif
               if (j > 0) then
                  read(funit,'(4d20.13)')(pot_l(ir),ir=1,jend)
                  j_inter=1
                  do ir=1,ThisP%Grid%jend
                     x=ThisP%Grid%x_mesh(ir)
!                    -------------------------------------------------
                     call hunt(jend,x_mesh,x,j_inter)
!                    -------------------------------------------------
                     if (j_inter > jend-(n_inter-1)/2) then
                        irp=jend-n_inter+1
                     else if (2*j_inter+1 > n_inter) then
                        irp=j_inter-(n_inter-1)/2
                     else
                        irp=1
                     endif
!                    -------------------------------------------------
                     call PolyInterp(n_inter,x_mesh(irp:irp+n_inter-1),     &
                                     pot_l(irp:irp+n_inter-1),x,vc,err)
!                    -------------------------------------------------
                     ThisP%pot_l(ir,jlr,is,ia)=vc
      !!             ThisP%pot_l(ir,jlr,is,ia)=pot_l(ir)
                  enddo
               endif
            enddo
         enddo
!        =============================================================
         deallocate(x_mesh, r_mesh, pot_l)
!        =============================================================
         SphericalInputFile = .false.
      else
         do is=1,n_spin_pola
            do ir=1,ThisP%Grid%jmt
               ThisP%pot_l(ir,1,is,ia)=ThisP%potr_sph(ir,is,ia)/(Y0*ThisP%Grid%r_mesh(ir))
            enddo
         enddo
      endif
      efermi=min(efermi,efermi_in(ia,id))
!
      close(unit=funit)
   enddo
   deallocate(rhotot)
!
   end subroutine readFormattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readUnformattedData(vunit,id)
!  ===================================================================
   use PotentialTypeModule, only : isMuffinTinPotential, isFullPotential
   use SystemModule, only : setLatticeConstant, getNumAtoms, getAtomName, &
                            getAlloyElementName
   use AtomModule, only: getLocalEvecOld
   use Atom2ProcModule, only : getGlobalIndex
   use RadialGridModule, only : getNumRmesh
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : getZtot, getZcor, getZsem
   use ChemElementModule, only : getNumCoreStates, setNumCoreStates
   use ChemElementModule, only : getCoreStateIndex
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark, &
                                       IntegerType, IntegerMark
!
   use MPPModule, only: GlobalMax
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: vunit
   integer (kind=IntKind), intent(in) :: id
!
   character (len=19), parameter :: sname='readUnformattedData'
   character (len=20) :: istop
   character (len=80) :: header
   character (len=5) :: lst(numcmax,n_spin_pola,MaxNumSpecies)
   character (len=MaxLenOfAtomName) :: atname(MaxNumSpecies)
!
   integer (kind=IntKind) :: na, id_g, nspin
   integer (kind=IntKind) :: jmt, jws, i, is, nr, nrmax, n, ia
   integer (kind=IntKind) :: numc(MaxNumSpecies), DataSize(LocalNumAtoms)
   integer (kind=IntKind) :: nc(numcmax,n_spin_pola,MaxNumSpecies)
   integer (kind=IntKind) :: lc(numcmax,n_spin_pola,MaxNumSpecies)
   integer (kind=IntKind) :: kc(numcmax,n_spin_pola,MaxNumSpecies)
   integer (kind=IntKind), pointer :: pc0(:,:,:)
!
   real (kind=RealKind) :: evec(3)
   real (kind=RealKind) :: ztotss(MaxNumSpecies), zcorss(MaxNumSpecies)
   real (kind=RealKind) :: xvalws(n_spin_pola,MaxNumSpecies)
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), pointer :: rho0(:,:)
   real (kind=RealKind), pointer :: mom0(:,:)
   real (kind=RealKind), pointer :: ec0(:,:,:)
   real (kind=RealKind), allocatable :: rhotot(:,:,:)
   real (kind=RealKind) :: ec(numcmax,n_spin_pola,MaxNumSpecies)
!
   type (PotentialStruct), pointer :: ThisP
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('readPotential','invalid id',id)
   endif
!
   if (n_spin_pola == 1) then
      nspin = 1
   else if (n_spin_cant == 1) then
      nspin = 2
   else
      nspin = 3
   endif
!
   ThisP=>Potential(id)
!
   na = getNumAtoms()
   nr = getNumRmesh(id)
   r_mesh => ThisP%Grid%r_mesh(1:nr)
   jmt = ThisP%Grid%jmt
   jws = ThisP%Grid%jend
!
   nrmax = nr+30
!  nrmax = nr
   call GlobalMax(nrmax)
!
   id_g = getGlobalIndex(id)
!
   do ia = 1, NumSpecies(id)
      atname(ia) = getAlloyElementName(id_g,ia)
      ztotss(ia) = getZtot(atname(ia))
      zcorss(ia) = getZcor(atname(ia))
   enddo
!
   allocate( rhotot(nr,n_spin_pola,NumSpecies(id)) )
!
!  -------------------------------------------------------------------
   call getpotg( id, na, NumSpecies(id), MaxNumSpecies,               &
                 nspin, n_spin_pola, efermi_in(1,id), evec,           &
                 nr, nrmax, jmt, jws, r_mesh, ThisP%jmax,             &
                 ThisP%potr_sph, vdif(1), rhotot, xvalws,             &
                 ztotss, zcorss, numc, nc, lc, kc, ec, lst, numcmax,  &
                 ThisP%pot_l, SphericalInputFile, header, vunit,      &
                 print_level(id), istop )
!  -------------------------------------------------------------------
!
   if ( isFullPotential() .and. .not.SphericalInputFile ) then
      do ia = 1, NumSpecies(id)
         do is = 1, n_spin_pola
            do i = 1, nr
               ThisP%potr_sph(i,is,ia) = ThisP%pot_l(i,1,is,ia)*(Y0*r_mesh(i))
            enddo
         enddo
      enddo
   else
      do ia = 1, NumSpecies(id)
         do is = 1, n_spin_pola
            do i = 1, nr
               ThisP%pot_l(i,1,is,ia) = ThisP%potr_sph(i,is,ia)/(Y0*r_mesh(i))
            enddo
         enddo
      enddo
   endif
!  ===================================================================
!
   ThisP%header = header
   v0(1:2) = ZERO
!
   rho0 => getDataStorage( id, 'OldSphericalElectronDensity',         &
                           nr+1, NumSpecies(id), RealMark )
   if (SiteMaxNumc(id) >= 1) then
      ec0  => getDataStorage( id, 'OldEstimatedCoreEnergy',           &
                              SiteMaxNumc(id),n_spin_pola, NumSpecies(id), RealMark )
      pc0 => getDataStorage( id, 'OldEstimatedCoreStates',            &
                             SiteMaxNumc(id), n_spin_pola, NumSpecies(id), IntegerMark)
   endif
   if (n_spin_pola == 2) then
      mom0 => getDataStorage(id,'OldSphericalMomentDensity', &
                             nr+1, NumSpecies(id), RealMark )
   endif
!
   efermi = 100.0d0
   do ia = 1, NumSpecies(id)
      ThisP%ztss(ia) = ztotss(ia)
      if ( numc(ia) > getNumCoreStates(atname(ia)) ) then
         if (print_level(id) >= 0) then
            call WarningHandler('readUnformattedData','numc > getNumCoreStates', &
                                numc(ia),getNumCoreStates(atname(ia)))
         endif
         call setNumCoreStates(atname(ia),numc(ia))
      else if (numc(ia) < getNumCoreStates(atname(ia)) ) then
         if (print_level(id) >= 0) then
            call WarningHandler('readUnformattedData','numc < getNumCoreStates', &
                                 numc(ia),getNumCoreStates(atname(ia)))
         endif
         call setNumCoreStates(atname(ia),numc(ia))
      endif
!
      efermi = min(efermi, efermi_in(ia,id))
!
      do is = 1, n_spin_pola
         do i = 1, numc(ia)
            ec0(i,is,ia) = ec(i,is,ia)
            pc0(i,is,ia) = getCoreStateIndex(nc(i,is,ia),lc(i,is,ia),kc(i,is,ia))
         enddo
      enddo
!
      if (n_spin_pola == 1) then
         do i=1,jws
            rho0(i,ia) = rhotot(i,1,ia)/(PI4*r_mesh(i)*r_mesh(i))
         enddo
         do i=jws+1,nr
            rho0(i,ia) = ZERO
         enddo
         rho0(nr+1,ia) = xvalws(1,ia)
      else
         do i = 1,jws
            rho0(i,ia) = (rhotot(i,1,ia) + rhotot(i,2,ia))/(PI4*r_mesh(i)*r_mesh(i))
         enddo
         do i=jws+1,nr
            rho0(i,ia) = ZERO
         enddo
         do i = 1,jws
            mom0(i,ia) = (rhotot(i,1,ia) - rhotot(i,2,ia))/(PI4*r_mesh(i)*r_mesh(i))
         enddo
         do i=jws+1,nr
            mom0(i,ia) = ZERO
         enddo
         rho0(nr+1,ia) = xvalws(1,ia)+xvalws(2,ia)
         mom0(nr+1,ia) = xvalws(1,ia)-xvalws(2,ia)
      endif
   enddo
!
   nullify(ec0,pc0,rho0,mom0)
   deallocate( rhotot )
!
   end subroutine readUnformattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writePotential()
!  ===================================================================
   use MPPModule, only : syncAllPEs
!
   use AtomModule, only : getOutPotFileName, getOutPotFileForm
!
   use ParallelIOModule, only : isOutputProc, getMyOutputProc
!
   implicit none
!
   character (len=10), parameter :: file_status = 'WRITE'
   character (len=11) :: file_form
!
   integer (kind=IntKind) :: id, str_len, wunit, ia
   integer (kind=IntKind), parameter :: funit=92
!
   file_form = getOutPotFileForm()
!
   if (file_form == 'HDF') then
      print *,'HDF file I/O needs to be implemented'
      stop 'incomplete'
   else if ( .not. (file_form == 'FORMATTED' .or.                    &
                    trim(file_form) == 'UNFORMATTED' .or.            &
                    trim(file_form) == 'XDR') ) then
      call ErrorHandler('writePotential','Unknown potential data format', &
                        getOutPotFileForm())
   endif
!
!  ------------------------------------------------------------------
   call syncAllPEs()
!  ------------------------------------------------------------------
   if ( getMyOutputProc() >= 0 ) then
      if (file_form == 'FORMATTED') then
         do id = 1, LocalNumAtoms
!           ---------------------------------------------------------
            call writeFormattedData(id)
!           ---------------------------------------------------------
         enddo
      else if (trim(file_form) == 'UNFORMATTED' .or. trim(file_form) == 'XDR') then
         if( isOutputProc() ) then
            str_len = len_trim(adjustl(getOutPotFileName()))
!           ----------------------------------------------------------
            call c_gopen(wunit,adjustl(getOutPotFileName()),str_len,     &
                         trim(file_status),len_trim(file_status),        &
                         trim(file_form),len_trim(file_form))
!           ----------------------------------------------------------
         endif
         do id = 1, LocalNumAtoms
!           ---------------------------------------------------------
            call writeUnformattedData(wunit,id)
!           ---------------------------------------------------------
         enddo
         if( isOutputProc() ) then
            call c_close(wunit)
         endif
      endif
   endif
!  ------------------------------------------------------------------
   call syncAllPEs()
!  ------------------------------------------------------------------
!
   end subroutine writePotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFormattedData(id)
!  ===================================================================
   use PublicParamDefinitionsModule, only : MaxLenFileName
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : getZcor, getZsem, getZval, getZtot
   use ChemElementModule, only : getCoreStateN, getCoreStateL, getNumCoreStates
   use ChemElementModule, only : getCoreStateKappa, getCoreStateSymbol
!
   use DataServiceCenterModule, only : getDataStorage,                &
                                       RealType, ComplexType,         &
                                       RealMark, ComplexMark,         &
                                       IntegerType, IntegerMark
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SystemModule, only : getAtomName, getAlloyElementName
!
   use AtomModule, only : getLocalEvecNew
   use AtomModule, only : getOutPotFileName, getOutPotFileForm
!
   use RadialGridModule, only : getNumRmesh
!
   implicit none
!
   character (len=MaxLenFileName) :: filename
   character (len=MaxLenOfAtomName) :: atname
   character (len=7) :: stmp
   character (len=80) :: header , jtitle
   character (len=120) :: ofile
   character (len=5), allocatable  :: lst(:,:)
!
   integer (kind=IntKind), intent(in) :: id
!
   integer (kind=IntKind) :: ns,ir,ic,jl,n,ia
   integer (kind=IntKind) :: nspin
   integer (kind=IntKind) :: nr,nrcor,numc,present_atom,jmt,jws
   integer (kind=IntKind), allocatable :: nc(:,:), lc(:,:), kc(:,:)
   integer (kind=IntKind), pointer :: pc0(:,:,:)
   integer (kind=IntKind), parameter :: ounit = 201
!
   real (kind=RealKind) :: ztotss, zcorss, zsemss, zvalss
   real (kind=RealKind) :: xstart,xmt,rmt,xvalws(2),evec(3)
   real (kind=RealKind), pointer :: rho0(:,:), mom0(:,:), ec0(:,:,:), r_mesh(:)
   real (kind=RealKind), allocatable :: ec(:,:), rhotot(:,:)
!
   type (PotentialStruct), pointer :: ThisP
!
   ThisP=>Potential(id)
!
   present_atom = getGlobalIndex(id)
!
   header = ThisP%header
   jmt = ThisP%Grid%jmt
   jws = ThisP%Grid%jend
   xmt = ThisP%Grid%xmt
   rmt = ThisP%Grid%rmt
   xstart = ThisP%Grid%xstart
   evec = getLocalEvecNew(id)
   nr = getNumRmesh(id)
   r_mesh => ThisP%Grid%r_mesh(1:nr)
!
   if (n_spin_pola == 1) then
      nspin = 1
   else if (n_spin_cant == 1) then
      nspin = 2
   else
      nspin = 3
   endif
!
   allocate( nc(numcmax,n_spin_pola) )
   allocate( lc(numcmax,n_spin_pola) )
   allocate( kc(numcmax,n_spin_pola) )
   allocate( lst(numcmax,n_spin_pola) )
   allocate( ec(numcmax,n_spin_pola) )
   allocate( rhotot(nr,n_spin_pola) )
!
!  -------------------------------------------------------------------
   rho0 => getDataStorage( id, 'NewSphericalElectronDensity',         &
                           nr+1, NumSpecies(id), RealMark)
!  -------------------------------------------------------------------
   if (n_spin_pola == 2) then
!     ----------------------------------------------------------------
      mom0 => getDataStorage(id,'NewSphericalMomentDensity',          &
                             nr+1, NumSpecies(id), RealMark )
!     ----------------------------------------------------------------
   endif
!
   if (SiteMaxNumc(id) >= 1) then
!     ----------------------------------------------------------------
      ec0 => getDataStorage( id, 'NewEstimatedCoreEnergy',            &
                             SiteMaxNumc(id), n_spin_pola, NumSpecies(id), RealMark)
      pc0 => getDataStorage( id, 'NewEstimatedCoreStates',            &
                             SiteMaxNumc(id), n_spin_pola, NumSpecies(id), IntegerMark)
!     ----------------------------------------------------------------
   endif
!
   do ia = 1, NumSpecies(id)
      filename = getOutPotFileName(id,ia)
      atname = getAlloyElementName(present_atom,ia)
      ztotss = getZtot(atname)
      zcorss = getZcor(atname)
      zsemss = getZsem(atname)
      zvalss = getZval(atname)
      numc = getNumCoreStates(atname)
!
      do ns = 1,n_spin_pola
         do ic=1,numc
            n = pc0(ic,ns,ia)
            nc(ic,ns)=getCoreStateN(atname,n)
            lc(ic,ns)=getCoreStateL(atname,n)
            kc(ic,ns)=getCoreStateKappa(atname,n)
            lst(ic,ns)= getCoreStateSymbol(atname,n)
            ec(ic,ns) = ec0(ic,ns,ia)
         enddo
      enddo
!
      if (n_spin_pola == 1) then
         do ir = 1,jws
            rhotot(ir,1) = PI4*rho0(ir,ia)*r_mesh(ir)*r_mesh(ir)
         enddo
         xvalws(1) = rho0(nr+1,ia)
      else
         do ir = 1,jws
            rhotot(ir,1) = PI2*(rho0(ir,ia)+mom0(ir,ia))*r_mesh(ir)*r_mesh(ir)
         enddo
         do ir = 1,jws
            rhotot(ir,2) = PI2*(rho0(ir,ia)-mom0(ir,ia))*r_mesh(ir)*r_mesh(ir)
         enddo
         xvalws(1) = (rho0(nr+1,ia)+mom0(nr+1,ia))/TWO
         xvalws(2) = (rho0(nr+1,ia)-mom0(nr+1,ia))/TWO
      endif
!
      write(stmp,'(i5,a,i1)')10000+present_atom,'x',ia
      stmp(1:1) = '_'
      write(ofile,'(a,a)')trim(adjustl(filename)),stmp
      open(unit=ounit,file=trim(ofile),form='formatted',status='unknown')
!
      write(ounit,'(a)') header
      write(ounit,'(i5,2x,d20.13)') n_spin_pola,vdif(1)
!
      do ns = 1, n_spin_pola
         write(jtitle,'(a,i5,3x,a2,4(a,f3.0),a,f10.5)')               &
            ' PUTPOTG:',present_atom,atname,'  zt=',ztotss,           &
            ' zc=',zcorss,' zs=',zsemss,' zv=',zvalss,                &
            ' xv=',xvalws(ns)
!
         write(ounit,'(a,t18,a2,t25,a,f4.0,t35,a,f10.5)')             &
                  ' LSMS:',atname,'z=',ztotss,'xvalws=',xvalws(ns)
!        write(ounit,'(f5.0,17x,f12.5,1x,f5.0,1x,d20.13)')            &
         write(ounit,'(f5.0,17x,f12.5,f5.0,d20.13)')            &
                    ztotss,rmt,zcorss,efermi
         write(ounit,'(17x,2d20.13,2i5)') xstart,xmt,jmt,jws
!
         write(ounit,'(4d20.13)') (ThisP%potr_sph(ir,ns,ia),ir=1,jmt)
         write(ounit,'(35x,d20.13)') ZERO
!
         write(ounit,'(i5,d20.13)') jws,xvalws(ns)
         write(ounit,'(4d20.13)') (rhotot(ir,ns),ir=1,jws)
!
         nrcor = 0
         write(ounit,'(2i5)') numc,nrcor
         do ic = 1, numc
            write(ounit,'(3i5,f12.5,2x,a5)') nc(ic,ns),lc(ic,ns),     &
                                             kc(ic,ns),ec(ic,ns),lst(ic,ns)
         enddo
      enddo
!
      if (ThisP%lmax > 0) then
         write(ounit,'(5i8)')ThisP%lmax,ThisP%jmax,ThisP%Grid%jmt,ThisP%Grid%jend
         write(ounit,'(4d20.13)')ThisP%Grid%xstart,ThisP%Grid%xmt,ThisP%Grid%hin,ThisP%Grid%hout
         do ns = 1,n_spin_pola
            do jl = 1, ThisP%jmax
               write(ounit,'(2i5)')jl,ThisP%PotCompFlag(jl)
               if (ThisP%PotCompFlag(jl) > 0) then
                  write(ounit,'(4d20.13)')(ThisP%pot_l(ir,jl,ns,ia),ir=1,ThisP%Grid%jend)
               endif
            enddo
         enddo
      endif
      close(unit=ounit)
   enddo
!
   deallocate( nc, lc, kc, lst, ec, rhotot )
   nullify( ec0, pc0, rho0, r_mesh )
   if (n_spin_pola == 2) then
      nullify( mom0 )
   endif
!
   end subroutine writeFormattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeUnformattedData(wunit,id)
!  ===================================================================
   use SystemModule, only : getLatticeConstant, getNumAtoms, getAtomName, &
                            getAlloyElementName
   use AtomModule, only: getLocalEvecNew
   use Atom2ProcModule, only : getGlobalIndex
   use RadialGridModule, only : getGrid, getNumRmesh
   use ChemElementModule, only : MaxLenOfAtomName
   use ChemElementModule, only : getZtot, getZcor, getZsem, getZval
   use ChemElementModule, only : getNumCoreStates
   use ChemElementModule, only : getCoreStateN
   use ChemElementModule, only : getCoreStateL
   use ChemElementModule, only : getCoreStateKappa, getCoreStateSymbol
   use DataServiceCenterModule, only : getDataStorage,        &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark, &
                                       IntegerMark, IntegerType
!
   implicit none
   integer (kind=IntKind), intent(in) :: wunit
   integer (kind=IntKind), intent(in) :: id
!
   character (len=20) :: istop
   character (len=80) :: header
   character (len=MaxLenOfAtomName) :: atname(MaxNumSpecies)
   character (len=5) :: lst(numcmax,n_spin_pola,MaxNumSpecies)
!
   integer (kind=IntKind) :: na, id_g, n, ia
   integer (kind=IntKind) :: jmt, jws
   integer (kind=IntKind) :: nspin, i, is, nr, ic
   integer (kind=IntKind), pointer :: pc0(:,:,:)
   integer (kind=IntKind) :: nc(numcmax,n_spin_pola,MaxNumSpecies)
   integer (kind=IntKind) :: lc(numcmax,n_spin_pola,MaxNumSpecies)
   integer (kind=IntKind) :: kc(numcmax,n_spin_pola,MaxNumSpecies)
   integer (kind=IntKind) :: numc(MaxNumSpecies)
!
   real (kind=RealKind) :: xstart, xmt, evec(3)
   real (kind=RealKind) :: ztotss(MaxNumSpecies)
   real (kind=RealKind) :: zcorss(MaxNumSpecies)
   real (kind=RealKind) :: zsemss(MaxNumSpecies)
   real (kind=RealKind) :: zvalss(MaxNumSpecies)
   real (kind=RealKind) :: xvalws(n_spin_pola,MaxNumSpecies)
   real (kind=RealKind) :: ec(numcmax,n_spin_pola,MaxNumSpecies)
   real (kind=RealKind), allocatable :: rhotot(:,:,:)
   real (kind=RealKind), allocatable :: r_pot_l(:)
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), pointer :: rho0(:,:)
   real (kind=RealKind), pointer :: mom0(:,:)
   real (kind=RealKind), pointer :: ec0(:,:,:)
!
   type (GridStruct),pointer :: Grid
!
   if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('writePotential','invalid id',id)
   endif
!
   id_g = getGlobalIndex(id)
!
   header = Potential(id)%header
   na = getNumAtoms()
   Grid => getGrid(id)
   jmt = Grid%jmt
   jws = Grid%jend
   xmt = Grid%xmt
   xstart = Grid%xstart
   evec = getLocalEvecNew(id)
   nr = getNumRmesh(id)
   r_mesh => Grid%r_mesh(1:nr)
!
   do ia = 1, NumSpecies(id)
      atname(ia) = getAlloyElementName(id_g,ia)
      ztotss(ia) = getZtot(atname(ia))
      zcorss(ia) = getZcor(atname(ia))
      zsemss(ia) = getZsem(atname(ia))
      zvalss(ia) = getZval(atname(ia))
      numc(ia) = getNumCoreStates(atname(ia))
   enddo
!
   if (n_spin_pola == 1) then
      nspin = 1
   else if (n_spin_cant == 1) then
      nspin = 2
   else
      nspin = 3
   endif
!
   if (SiteMaxNumc(id) >= 1) then
!     ----------------------------------------------------------------
      ec0 => getDataStorage( id, 'NewEstimatedCoreEnergy',            &
                             SiteMaxNumc(id), n_spin_pola, NumSpecies(id), RealMark)
      pc0 => getDataStorage( id, 'NewEstimatedCoreStates',            &
                             SiteMaxNumc(id), n_spin_pola, NumSpecies(id), IntegerMark)
!     ----------------------------------------------------------------
   endif
!
   do ia = 1, NumSpecies(id)
      do is = 1,n_spin_pola
         do ic = 1, numc(ia)
            n = pc0(ic,is,ia)
            nc(ic,is,ia)=getCoreStateN(atname(ia),n)
            lc(ic,is,ia)=getCoreStateL(atname(ia),n)
            kc(ic,is,ia)=getCoreStateKappa(atname(ia),n)
            lst(ic,is,ia)=getCoreStateSymbol(atname(ia),n)
            ec(ic,is,ia)=ec0(ic,is,ia)
         enddo
      enddo
   enddo
!
   allocate(rhotot(nr,n_spin_pola,NumSpecies(id)))
!  -------------------------------------------------------------------
   rho0 => getDataStorage( id, 'NewSphericalElectronDensity',         &
                           nr+1, NumSpecies(id), RealMark)
!  -------------------------------------------------------------------
!
   if ( n_spin_pola == 1 ) then
      do ia = 1, NumSpecies(id) 
         do i = 1,jws
            rhotot(i,1,ia) = PI4*rho0(i,ia)*r_mesh(i)*r_mesh(i)
         enddo
         xvalws(1,ia) = rho0(nr+1,ia)
      enddo
   else
!     ----------------------------------------------------------------
      mom0 => getDataStorage(id,'NewSphericalMomentDensity',          &
                             nr+1, NumSpecies(id), RealMark )
!     ----------------------------------------------------------------
      do ia = 1, NumSpecies(id)
         do i = 1,jws
            rhotot(i,1,ia) = PI2*(rho0(i,ia)+mom0(i,ia))*r_mesh(i)*r_mesh(i)
         enddo
         do i = 1,jws
            rhotot(i,2,ia) = PI2*(rho0(i,ia)-mom0(i,ia))*r_mesh(i)*r_mesh(i)
         enddo
         xvalws(1,ia) = (rho0(nr+1,ia)+mom0(nr+1,ia))/TWO
         xvalws(2,ia) = (rho0(nr+1,ia)-mom0(nr+1,ia))/TWO
      enddo
   endif
!
!  n = nr*Potential(id)%jmax*n_spin_pola
   n = jwsmax*Potential(id)%jmax*n_spin_pola*NumSpecies(id)
   allocate( r_pot_l(2*n) )
   r_pot_l = transfer(Potential(id)%pot_l,r_pot_l)
!  -------------------------------------------------------------------
   call putpotg(id, na, NumSpecies(id),                               &
                nspin, n_spin_pola, efermi, evec,                     &
                jmt, xmt, jws, xstart, Potential(id)%jmax,            &
                Potential(id)%potr_sph, vdif(1), rhotot, xvalws,      &
                atname, ztotss, zcorss, zsemss, zvalss,               &
                numc, nc, lc, kc, ec, lst,                            &
                jmtmax, jwsmax, nr, numcmax,                          &
                LdaPlusU_DataAccum, NonSphPot_DataAccum,              &
                r_pot_l, header, wunit, print_level(id), istop )
!  -------------------------------------------------------------------
!
   nullify(ec0, pc0, rho0, r_mesh)
   if (n_spin_pola == 2) then
      nullify( mom0 )
   endif
   deallocate(rhotot, r_pot_l)
!
   end subroutine writeUnformattedData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPotCompFlag(id,is,ia)
!  ===================================================================
   use SystemSymmetryModule, only : getSymmetryFlags
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: is, ia
!
   integer (kind=IntKind) :: ir, jl
   integer (kind=IntKind), pointer :: flags(:)
!
   Potential(id)%PotCompFlag(1) = 1
!
   if ( .not.isChargeSymmOn ) then
      do jl=2,Potential(id)%jmax
         Potential(id)%PotCompFlag(jl) = 0
         LOOP_ir: do ir = 1,Potential(id)%Grid%jend
            if (abs(Potential(id)%pot_l(ir,jl,is,ia)) > pot_tol) then
               Potential(id)%PotCompFlag(jl) = 1
               exit LOOP_ir
            endif
         enddo LOOP_ir
      enddo
   else
      flags => getSymmetryFlags(id)
      do jl=2,Potential(id)%jmax
         Potential(id)%PotCompFlag(jl) = flags(jl)
      enddo
   endif
   nullify(flags)
!
   end subroutine setPotCompFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setPotentialAccess(key_l, key_r)
!  ===================================================================
   use DataServiceCenterModule, only : createDataStorage,     &
                                       getDataStorage,        &
                                       isDataStorageExisting, &
                                       setDataStorage2Value,  &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark
   implicit none
!
   character (len=*), intent(in) :: key_l, key_r
!
   integer (kind=IntKind) :: i, jmax_pot, jend
   integer (kind=IntKind) :: DataSize(LocalNumAtoms)
!
   if ( isDataStorageExisting(key_l) .and. &
        isDataStorageExisting(key_r) ) then
!
      do i = 1,LocalNumAtoms
         jmax_pot = Potential(i)%jmax
         jend = Potential(i)%Grid%jend
         Potential(i)%pot_l=>getDataStorage( i,key_l,jend,    &
                                    jmax_pot,n_spin_pola,NumSpecies(i),ComplexMark)
         Potential(i)%potr_sph=>getDataStorage( i,key_r,jend, &
                                    n_spin_pola,NumSpecies(i),RealMark)
      enddo
   else
      do i = 1,LocalNumAtoms
         DataSize(i) = Potential(i)%Grid%jend*NumSpecies(i)*n_spin_pola
      enddo
!
      if ( .not.isDataStorageExisting(key_r) ) then
!        ------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,key_r,DataSize,RealType)
!        ------------------------------------------------------------
         call setDataStorage2Value(key_r,zero)
!        ------------------------------------------------------------
      endif
!
      do i = 1,LocalNumAtoms
         jmax_pot = Potential(i)%jmax
         if ( jmax_pot < 1 ) then
            call ErrorHandler('setPotentialAccess',                  &
                              "invalid jmax_pot",jmax_pot)
         endif
         DataSize(i) = Potential(i)%Grid%jend*jmax_pot*NumSpecies(i)*n_spin_pola
      enddo
!
      if ( .not.isDataStorageExisting(key_l) ) then
!        ------------------------------------------------------------
         call createDataStorage(LocalNumAtoms,key_l,DataSize,ComplexType)
!        ------------------------------------------------------------
         call setDataStorage2Value(key_l,czero)
!        ------------------------------------------------------------
      endif
!
      do i = 1,LocalNumAtoms
         jmax_pot = Potential(i)%jmax
         jend = Potential(i)%Grid%jend
         Potential(i)%pot_l=>getDataStorage( i,key_l,jend,    &
                                    jmax_pot,n_spin_pola,NumSpecies(i),ComplexMark)
         Potential(i)%potr_sph=>getDataStorage( i,key_r,jend, &
                                    n_spin_pola,NumSpecies(i),RealMark)
      enddo
!
   endif
!
   end subroutine setPotentialAccess
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine truncatePotential(id,ia)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs
!
   use StepFunctionModule, only : truncate
!
   use SystemSymmetryModule, only : getSymmetryFlags
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ia
!
   character (len=9*Potential(id)%jmax) :: aflag
!
   integer (kind=IntKind) :: jmax_pot, lmax_pot, kmax_pot
   integer (kind=IntKind) :: jmax_trunc, lmax_trunc, jmax_min
   integer (kind=IntKind) :: ir, iend, jmt, iend_diff, is, jl, i, n
   integer (kind=IntKind) :: izamax, ModifiedTruncation
   integer (kind=IntKind), pointer :: flags_jl(:), p_flags(:)
!
   real (kind=RealKind) :: pot_r, pot_i
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind), pointer :: pot(:,:), potl(:,:)
!
   iend = Potential(id)%Grid%jend
   jmt  = Potential(id)%Grid%jmt
   iend_diff = iend-jmt+1
   lmax_pot = Potential(id)%lmax
   jmax_pot = Potential(id)%jmax
   kmax_pot = (lmax_pot+1)*(lmax_pot+1)
   lmax_trunc = Potential(id)%lmax_trunc
!
   jmax_trunc = ((lmax_trunc+1)*(lmax_trunc+2))/2
   r_mesh => Potential(id)%Grid%r_mesh
!
!  do is = 1,n_spin_pola
!     potl =>Potential(id)%pot_l_trunc(1:iend_diff,1:jmax_trunc,is,ia)
!     potl = CZERO
!  enddo
!
   do is = 1,n_spin_pola
      pot => Potential(id)%pot_l(1:iend,1:jmax_pot,is,ia)
      flags_jl => Potential(id)%PotCompFlag
!
      potl =>Potential(id)%pot_l_trunc(1:iend_diff,1:jmax_trunc,is,ia)
      potl = CZERO
!
!     ================================================================
!     Print out non-zero potential channels...........................
!     ================================================================
      if (MyPE == 0) then
         ir = 0
         do jl = 1, jmax_pot
            if (flags_jl(jl) /= 0) then
               write(aflag(ir+1:ir+9),'(a,i2,a,i2,a)')'(',lofj(jl),',',mofj(jl),'), '
               ir = ir + 9
            endif
         enddo
         write(6,'(/,80(''=''))')
         write(6,'(a)')'Potential non-zero channels:'
         n = ceiling((ir-2)/90.0)
         do i = 1, n
            if (i*90 <= ir-2) then
               write(6,'(a)')aflag((i-1)*90+1:i*90)
            else
               write(6,'(a)')aflag((i-1)*90+1:ir-2)
            endif
!           write(6,'(a)')aflag(1:ir-2)
         enddo
      endif
!     ---------------------------------------------------------------
      call truncate( id, jmt, iend, r_mesh, pot, flags_jl, jmax_pot, &
                     iend_diff, potl, jmax_trunc )
!     ---------------------------------------------------------------
!
      flags_jl => Potential(id)%PotCompFlag_trunc
!
      if ( is==1 ) then
         if ( .not.isChargeSymmOn ) then
            do jl = 1,jmax_trunc
!
               ir = izamax(iend_diff,potl(1:iend_diff,jl),1)
!              if ( abs(potl(ir,jl)) >= TEN2m6*(TEN**(-2*lofj(jl))) ) then
               if ( abs(potl(ir,jl)) > pot_tol) then
                  flags_jl(jl) = 1
               else 
                  flags_jl(jl) = 0
                  potl(1:iend_diff,jl) = CZERO
               endif
!
if(.false.) then
               if ( flags_jl(jl)/=0 ) then
                  if ( mofj(jl) /= 0 ) then
                     do ir = 1,iend_diff
                        pot_r = real( potl(ir,jl),kind=RealKind )
                        pot_i = real( -SQRTm1*potl(ir,jl),kind=RealKind )
                        if ( (abs(pot_r) /= ZERO) .and.               &
                             (abs(pot_i)/abs(pot_r) < TEN2m8) ) then
                           potl(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                        else if ( (abs(pot_i) /= ZERO) .and.          &
                             (abs(pot_r)/abs(pot_i) < TEN2m8) ) then
                           potl(ir,jl) = cmplx(ZERO,pot_i,kind=CmplxKind)
                           flags_jl(jl) = 2
                        else
                           flags_jl(jl) = 3
                        endif
                     enddo
                  else
                     do ir = 1,iend_diff
                        pot_r = real( potl(ir,jl),kind=RealKind)
                        potl(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                     enddo
                  endif
               else
                  potl(1:iend_diff,jl) = CZERO
               endif
endif
!
            enddo
         else
!           ==========================================================
!           The following statement overwrites PotCompFlag_trunc with
!           symmetry flags of the system
!           ==========================================================
            p_flags => getSymmetryFlags(id)
            jmax_min = min( size(p_flags,1), jmax_trunc )
            do jl = 1, jmax_min
               flags_jl(jl) = p_flags(jl)
            enddo
!
            do jl = 1,jmax_trunc
               if ( flags_jl(jl) == 0 ) then
                  potl(1:iend_diff,jl) = CZERO
               endif
if(.false.) then
               if ( flags_jl(jl) == 0 ) then
                  potl(1:iend_diff,jl) = ZERO
               else if ( flags_jl(jl)==1 ) then
                  do ir = 1,iend_diff
                     pot_r = real( potl(ir,jl),kind=RealKind)
                     potl(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
                  enddo
               else if ( flags_jl(jl)==2 ) then
                  do ir = 1,iend_diff
                     pot_i = real( -SQRTm1*potl(ir,jl),kind=RealKind )
                     potl(ir,jl) = cmplx(ZERO,pot_i,kind=CmplxKind)
                  enddo
               endif
endif
            enddo
!
         endif
!
!        =============================================================
!        Print out non-zero truncated potential channels..............
!        =============================================================
         if (MyPE == 0) then
            ir = 0
            do jl = 1, jmax_trunc
               if (flags_jl(jl) /= 0) then
                  write(aflag(ir+1:ir+9),'(a,i2,a,i2,a,$)')'(',lofj(jl),',',mofj(jl),'), '
                  ir = ir + 9
               endif
            enddo
            write(6,'(/,a)')'Truncated potential non-zero channels:'
            n = ceiling((ir-2)/90.0)
            do i = 1, n
               if (i*90 <= ir-2) then
                  write(6,'(a)')aflag((i-1)*90+1:i*90)
               else
                  write(6,'(a)')aflag((i-1)*90+1:ir-2)
               endif
            enddo
!           write(6,'(a)')aflag(1:ir-2)
            write(6,'(80(''=''),/)')
         endif
      else
         do jl = 1,jmax_trunc
            if ( flags_jl(jl) == 0 ) then
               potl(1:iend_diff,jl) = CZERO
            endif
if(.false.) then
            if ( flags_jl(jl) == 0 ) then
               potl(1:iend_diff,jl) = ZERO
            else if ( flags_jl(jl)==1 ) then
               do ir = 1,iend_diff
                  pot_r = real( potl(ir,jl),kind=RealKind)
                  potl(ir,jl) = cmplx(pot_r,ZERO,kind=CmplxKind)
               enddo
            else if ( flags_jl(jl)==2 ) then
               do ir = 1,iend_diff
                  pot_i = real( -SQRTm1*potl(ir,jl),kind=RealKind )
                  potl(ir,jl) = cmplx(ZERO,pot_i,kind=CmplxKind)
               enddo
            endif
endif
         enddo
      endif
!
   enddo
!
   nullify( flags_jl, p_flags, r_mesh, pot, potl )
!
   end subroutine truncatePotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPot_L(id, pot_type, aux_name)
!  ===================================================================
   use MPPModule, only : MyPE
   use MathParamModule, only: TEN2m8
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   character (len=*), intent(in), optional :: pot_type, aux_name
!
   character(len=50) :: file_potl
   integer (kind=IntKind) :: is, l, m, jl, jmax, ir, NumRs, funit, count_flag
   integer (kind=IntKind) :: nc, jmt, len_potname, ir_start, ia
   integer (kind=IntKind) :: offset = 100000
   integer (kind=IntKind) :: offset_at = 1000
   integer (kind=IntKind), pointer :: flags(:)
!
   real (kind=RealKind), pointer :: r_mesh(:)
!
   complex (kind=CmplxKind) :: pot
   complex (kind=CmplxKind), pointer :: potl(:,:,:)
!
   do ia = 1, NumSpecies(id)
      file_potl = " "
      if ( .not.present(pot_type) ) then
         nc = 5
         file_potl(1:nc)  = "PotL_"
      else
         len_potname = len(pot_type)
         nc = 5+len_potname
         file_potl(4+len_potname:nc)  = "L_"
         file_potl(1:3)  = "Pot"
         file_potl(4:3+len_potname)  = pot_type(1:len_potname)
      endif
      if ( present(aux_name) ) then
         len_potname = len(trim(aux_name))
         file_potl(nc+1:nc+len_potname) = aux_name(1:len_potname)
         nc = nc + len_potname
      endif
!
      write(file_potl(nc+1:nc+6),'(i6)') offset+MyPE
      file_potl(nc+1:nc+1) = 'n'
      file_potl(nc+7:nc+7) = '_'
      write(file_potl(nc+8:nc+11),'(i4)') offset_at+id
      file_potl(nc+8:nc+8) = 'a'
      if (ia < 10) then
         write(file_potl(nc+9:nc+9),'(i1)')ia
      else
         write(file_potl(nc+9:nc+10),'(i2)')ia
      endif
      funit = 55+MyPE+id
!
      ir_start = 1
      if (present(pot_type) ) then
         if ( pot_type(1:3)=="Old" ) then
            NumRs = Potential(id)%Grid%jend
            jmax =  Potential(id)%jmax
            r_mesh => Potential(id)%Grid%r_mesh(1:NumRs)
            potl => Potential(id)%pot_l(1:NumRs,1:jmax,1:n_spin_pola,ia)
            flags => Potential(id)%PotCompFlag
         else if ( pot_type(1:5)=="Trunc" ) then
            if ( .not.isTruncatedPotential ) then
               call WarningHandler("printPot_L", &
                    "Truncated potential have not been initialized")
!               return
            endif
            if ( .not.Potential(id)%isTruncDone(ia) ) then
               call truncatePotential(id,ia)
               Potential(id)%isTruncDone(ia) = .true.
            endif
            NumRs = Potential(id)%iend_trunc
            jmax  = Potential(id)%jmax_trunc
            jmt   = Potential(id)%Grid%jmt
            r_mesh => Potential(id)%Grid%r_mesh(jmt:jmt+NumRs-1)
            potl => Potential(id)%pot_l_trunc(1:NumRs,1:jmax,1:n_spin_pola,ia)
            flags => Potential(id)%PotCompFlag_trunc
            ir_start = jmt
         endif
      else
         NumRs = Potential(id)%Grid%jend
         jmax =  Potential(id)%jmax
         r_mesh => Potential(id)%Grid%r_mesh(1:NumRs)
         potl => Potential(id)%pot_l(1:NumRs,1:jmax,1:n_spin_pola,ia)
         flags => Potential(id)%PotCompFlag
      endif
      open(unit=funit,file=trim(file_potl),status='unknown')
      count_flag = 0
      do jl = 1,jmax
         if ( flags(jl) /= 0 ) then
            count_flag = count_flag+1
         endif
      enddo
      write(funit,'(a,$)') "#Ind     r_mesh   (lm):"
      do jl = 1,jmax
         if ( flags(jl) /= 0 ) then
            l = lofj(jl)
            m = mofj(jl)
            write(funit,'(11x,a2,i2,1x,i2,a2,11x,$)') "( ", l, m," )"
         endif
      enddo
      write(funit,'(a)') "  "
      do ir = 1, NumRs
         write(funit,'(i4,1x,d16.8,$)') id, r_mesh(ir)
         JlLoop: do jl = 1,jmax
            if ( flags(jl) /= 0 ) then
               pot = CZERO
               do is =1, n_spin_pola
                  pot = pot + potl(ir,jl,is)
               enddo
               l = lofj(jl)
               write(funit,'(1x,(2(1x,d16.8)),$)') potl(ir,jl,1)
            endif
         enddo JlLoop
         write(funit,'(a)') " "
      enddo
      close(funit)
   enddo
!
   end subroutine printPot_L
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pushPotentialToAccel()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id, np, nj, jmax, iend, n, jl
   integer (kind=IntKind) :: iend_packed(LocalNumAtoms)
   integer (kind=IntKind) :: jmax_packed(LocalNumAtoms)
   integer (kind=IntKind), allocatable :: flags_packed(:)
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:,:,:)
   complex (kind=CmplxKind), allocatable :: pot_packed(:)
!
   np = 0; nj = 0
   do id = 1, LocalNumAtoms
      jmax=Potential(id)%jmax
      iend = Potential(id)%Grid%jend
      nj = nj + jmax
      np = np + iend*jmax*n_spin_pola*NumSpecies(id)
   enddo
   allocate(pot_packed(np),flags_packed(nj))
!
   np = 0; nj = 0
   do id = 1, LocalNumAtoms
      iend = Potential(id)%Grid%jend
      jmax=Potential(id)%jmax
      do jl = 1, jmax
         flags_packed(nj+jl) = Potential(id)%PotCompFlag(jl)
      enddo
      iend_packed(id) = iend
      jmax_packed(id) = jmax
      pot_l => Potential(id)%pot_l
      n = jmax*iend*n_spin_pola*NumSpecies(id)
!     ----------------------------------------------------------------
      call zcopy(n,pot_l,1,pot_packed(np+1),1)
!     ----------------------------------------------------------------
      nj = nj + jmax
      np = np + n
   enddo
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call initialize_potential(LocalNumAtoms,NumSpecies,iend_packed,jmax_packed,0)
   call push_potential(pot_packed,flags_packed,0)
!  -------------------------------------------------------------------
#endif
   deallocate(pot_packed, flags_packed)
!
   end subroutine pushPotentialToAccel
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pushTruncatedPotToAccel()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id, np, nj, jmax, iend, n, jl, ia
   integer (kind=IntKind) :: iend_packed(LocalNumAtoms)
   integer (kind=IntKind) :: jmax_packed(LocalNumAtoms)
   integer (kind=IntKind), allocatable :: flags_packed(:)
!
   complex (kind=CmplxKind), pointer :: pot_l(:,:,:,:)
   complex (kind=CmplxKind), allocatable :: pot_packed(:)
!
   np = 0; nj = 0
   do id = 1, LocalNumAtoms
      do ia = 1, NumSpecies(id)
         if ( .not.Potential(id)%isTruncDone(ia) ) then
            call truncatePotential(id,ia)
            Potential(id)%isTruncDone(ia) = .true.
         endif
      enddo
      iend = Potential(id)%iend_trunc
      jmax = Potential(id)%jmax_trunc
      nj = nj + jmax
      np = np + iend*jmax*n_spin_pola*NumSpecies(id)
   enddo
   allocate(pot_packed(np), flags_packed(nj))
!
   np = 0; nj = 0
   do id = 1, LocalNumAtoms
      iend = Potential(id)%iend_trunc
      jmax = Potential(id)%jmax_trunc
      do jl = 1, jmax
         flags_packed(nj+jl) = Potential(id)%PotCompFlag_trunc(jl)
      enddo
      iend_packed(id) = iend
      jmax_packed(id) = jmax
      pot_l => Potential(id)%pot_l_trunc
      n = jmax*iend*n_spin_pola*NumSpecies(id)
!     ----------------------------------------------------------------
      call zcopy(n,pot_l,1,pot_packed(np+1),1)
!     ----------------------------------------------------------------
      nj = nj + jmax
      np = np + n
   enddo
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call initialize_potential(LocalNumAtoms,NumSpecies,iend_packed,jmax_packed,1)
   call push_potential(pot_packed,flags_packed,1)
!  -------------------------------------------------------------------
#endif
   deallocate(pot_packed, flags_packed)
!  -------------------------------------------------------------------
!
   end subroutine pushTruncatedPotToAccel
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deletePotentialOnAccel()
!  ===================================================================
   implicit none
!
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call delete_potential(0)
!  -------------------------------------------------------------------
#endif
!
   end subroutine deletePotentialOnAccel
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteTruncatedPotOnAccel()
!  ===================================================================
   implicit none
!
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call delete_potential(1)
!  -------------------------------------------------------------------
#endif
!
   end subroutine deleteTruncatedPotOnAccel
!  ===================================================================
end module PotentialModule
