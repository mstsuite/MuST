! Needs some work and needs to be checked for multi-species case
module LdaCorrectionModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!   use ArrayToolsModule, only : copyArray, setTarget
!
public :: initLdaCorrection,         &
          endLdaCorrection,          &
          insertCorrectionOrbitals,  &
          computeLdaPlusU,           &
          getPotentialCorrection,    &
          getEnergyCorrection,       &
          getTransformMatrixInverse, &
          getTransformMatrix,        &
          transformScatteringMatrix, &
          transformWaveFunction,     &
          checkLdaCorrection,        &
          getNumCorrOrbitals,        &
          getCorrOrbital,            &
          getDataPackSize,           &
          getDataPack4Output,        &
          insertDataPackFromInput
          
!
   interface insertCorrectionOrbitals
      module procedure insert_withfile, insert_withdata0, insert_withdata1
   end interface
!
   interface checkLdaCorrection
      module procedure isLCN, isLCN_a, isLCN_al
   end interface
!
   interface getEnergyCorrection
      module procedure getEnergyCorrection_al, getEnergyCorrection_a
   end interface
!
   interface transformWaveFunction
      module procedure transformWF2, transformWF3
   end interface
!
private
!
   logical :: Initialized = .false.
   logical :: CorrectionNeeded = .false.
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
!
   integer (kind=IntKind), parameter :: MaxOrbitals = 2
!
   type CorrectionStruct
      integer (kind=IntKind) :: NumCorrOrbitals
      integer (kind=IntKind), pointer :: Orbitals(:)
      real (kind=RealKind), pointer :: Uparam(:)
      real (kind=RealKind), pointer :: Jparam(:)
      real (kind=RealKind), pointer :: SlaterIntegral(:,:)
      real (kind=RealKind), pointer :: v_lda_plus_u_full(:,:,:,:)
      complex (kind=CmplxKind), pointer :: v_lda_plus_u_diag(:,:,:)
      complex (kind=CmplxKind), pointer :: transform_inverse(:,:,:,:)
      complex (kind=CmplxKind), pointer :: transform_matrix(:,:,:,:)
      real (kind=RealKind), pointer :: e_lda_plus_u(:)
   end type CorrectionStruct
!
   type SiteStruct
      integer (kind=IntKind) :: NumSpecies
      integer (kind=IntKind), pointer :: AtomicNumber(:)
      type (CorrectionStruct), pointer :: SpeciesList(:)
   end type SiteStruct
!
   type (SiteStruct), allocatable :: AtomList(:)
!
   real (kind=RealKind) :: ak(-3:3,-3:3,-3:3,-3:3,2:3,1:4)
   real (kind=RealKind), allocatable, target :: DataPack(:)
   integer (kind=IntKind), allocatable :: DataPackSize(:)
   integer (kind=IntKind) :: MaxDataPackSize
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initLdaCorrection(na,index,getAtomicNumber,getNumSpecies)
!  ===================================================================
   use MathParamModule, only : ZERO, ONE, PI4
   use GauntFactorsModule, only : isGauntInitialized => isInitialized
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors,  &
                                  getNumK3, getK3, getGauntFactor
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: na
   integer (kind=IntKind), intent(in) :: index(na)
   integer (kind=IntKind) :: i, j, k, q, ic, n
   integer (kind=IntKind) :: l, lc, m, mp, mp2, mp3, klq, kl, klp, kn
!
   real (kind=RealKind) :: fac
   real (kind=RealKind) :: gfact(-6:6,-3:3,-3:3), gft(-6:6)
!
   logical :: gaunt_end
!
   interface
      function getAtomicNumber(ia,ic) result(z)
         use KindParamModule, only : IntKind
         implicit none
         integer (kind=IntKind), intent(in) :: ia
         integer (kind=IntKind), intent(in), optional :: ic
         integer (kind=IntKind) :: z
      end function getAtomicNumber
!
      function getNumSpecies(ia) result(n)
         use KindParamModule, only : IntKind
         implicit none
         integer (kind=IntKind), intent(in) :: ia
         integer (kind=IntKind) :: n
      end function getNumSpecies
   end interface
!
   if (na < 1) then
      call ErrorHandler('initLdaCorrection','Invalid number of local atoms',na)
      return
   endif
!
   allocate( AtomList(na), GlobalIndex(na) )
   allocate( DataPackSize(na) )
   GlobalIndex(1:na) = index(1:na)
!
   do i = 1, na
      n = getNumSpecies(i)
      if (n < 1) then
         call ErrorHandler('initLdaCorrection','Number of species < 1',n)
      endif
      AtomList(i)%NumSpecies = n
      allocate(AtomList(i)%SpeciesList(n))
      allocate(AtomList(i)%AtomicNumber(n))
      do ic = 1, AtomList(i)%NumSpecies
         AtomList(i)%AtomicNumber(ic) = getAtomicNumber(i,ic)
         AtomList(i)%SpeciesList(ic)%NumCorrOrbitals = 0
         nullify( AtomList(i)%SpeciesList(ic)%Orbitals )
         nullify( AtomList(i)%SpeciesList(ic)%Uparam, AtomList(i)%SpeciesList(ic)%Jparam )
         nullify( AtomList(i)%SpeciesList(ic)%v_lda_plus_u_full, AtomList(i)%SpeciesList(ic)%v_lda_plus_u_diag )
         nullify( AtomList(i)%SpeciesList(ic)%transform_inverse, AtomList(i)%SpeciesList(ic)%transform_matrix )
         nullify( AtomList(i)%SpeciesList(ic)%e_lda_plus_u )
      enddo
   enddo
!
   LocalNumAtoms = na
   DataPackSize(:) = 0
   MaxDataPackSize = 0
!
   Initialized = .true.
   CorrectionNeeded = .false.
!
   if (.not.isGauntInitialized()) then
      call initGauntFactors(3,'none',0)   ! l=3 for f-electrons
      gaunt_end = .true.
   else
      gaunt_end = .false.
   endif
!
   do j = 1, 4
      l = j - 1
      k = 2*l
      fac = 2*k+1; fac = PI4/fac
      do lc = 2, 3
!
         kn = lc*lc + lc + 1
         do mp = -lc, lc
            klp = kn + mp
            do m = -lc, lc
               kl = kn + m
               do q = -k, k
                  klq = kn + q
!                 ====================================================
!                 store the following Gaunt factor in gfact:
!
!                 gfact(q,m,mp,lc,j) =
!                        int_{4pi} Y^{*}_{lc,m} * Y_{l,q} * Y_{lc,mp}
!
!                 with kl = {lc, m}, klp = {lc,mp}, klq = {l,q}
!                 ====================================================
                  gfact(q,m,mp) = getGauntFactor(klp,kl,klq)
               enddo
            enddo
         enddo
!
         do mp3 = -lc, lc
            do mp2 = -lc, lc
               do q = -k, k
                  gft(q) = gfact(q,mp3,mp2)
               enddo
               do mp = -lc, lc
                  do m = -lc, lc
                     ak(m,mp,mp3,mp2,lc,j) = ZERO
                     do q = -k, k
                        ak(m,mp,mp2,mp3,lc,j) = ak(m,mp,mp2,mp3,lc,j) &
                                              + gfact(q,m,mp)*gft(q)
                     enddo
!                    =================================================
!                    ak(m,mp,mp2,mp3,lc,k) 
!                        = a_k(m,m',m'',m'''), defined in PRB 52, R5467 (1995)
!                    =================================================
                     ak(m,mp,mp2,mp3,lc,j) = fac*ak(m,mp,mp2,mp3,lc,j)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
!
   if (gaunt_end) then
      call endGauntFactors()
   endif
!
   end subroutine initLdaCorrection
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endLdaCorrection()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, ic
!
   do i = 1, LocalNumAtoms
      do ic = 1, AtomList(i)%NumSpecies
         if (AtomList(i)%SpeciesList(ic)%NumCorrOrbitals > 0) then
            deallocate( AtomList(i)%SpeciesList(ic)%Orbitals )
            deallocate( AtomList(i)%SpeciesList(ic)%Uparam )
            deallocate( AtomList(i)%SpeciesList(ic)%Jparam )
            deallocate( AtomList(i)%SpeciesList(ic)%v_lda_plus_u_full )
            deallocate( AtomList(i)%SpeciesList(ic)%v_lda_plus_u_diag )
            deallocate( AtomList(i)%SpeciesList(ic)%transform_inverse )
            deallocate( AtomList(i)%SpeciesList(ic)%transform_matrix )
            deallocate( AtomList(i)%SpeciesList(ic)%e_lda_plus_u )
            deallocate( AtomList(i)%SpeciesList(ic)%SlaterIntegral )
         endif
      enddo
      deallocate(AtomList(i)%AtomicNumber)
      deallocate(AtomList(i)%SpeciesList)
   enddo
!
   deallocate( AtomList, GlobalIndex, DataPackSize )
   if ( allocated(DataPack) ) then
      deallocate( DataPack )
   endif
   Initialized = .false.
   CorrectionNeeded = .false.
   MaxDataPackSize = 0
!
   end subroutine endLdaCorrection
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLCN() result(icn)
!  ===================================================================
   implicit none
!
   logical :: icn
!
   icn = CorrectionNeeded
!
   end function isLCN
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLCN_a(ia,ic) result(icn)
!  ===================================================================
   implicit none
!
   logical :: icn
   integer (kind=IntKind), intent(in) :: ia, ic
!
   icn = .false.
   if ( CorrectionNeeded ) then
      if (AtomList(ia)%SpeciesList(ic)%NumCorrOrbitals > 0) then
         icn = .true.
      endif
   endif
!
   end function isLCN_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLCN_al(ia,ic,l) result(icn)
!  ===================================================================
   implicit none
!
   logical :: icn
   integer (kind=IntKind), intent(in) :: ia, ic, l
   integer (kind=IntKind) :: i
!
   icn = .false.
   if ( CorrectionNeeded ) then
      LOOP_i: do i = 1, AtomList(ia)%SpeciesList(ic)%NumCorrOrbitals
         if ( l == AtomList(ia)%SpeciesList(ic)%Orbitals(i) ) then
            icn = .true.
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   end function isLCN_al
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insert_withdata0(id,ic,norb,lc,uparam,jparam)
!  ===================================================================
   use MathParamModule, only : ZERO, CZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind), intent(in) :: norb
   integer (kind=IntKind), intent(in) :: lc(norb)
   integer (kind=IntKind) :: i, n
!
   real (kind=RealKind), intent(in) :: uparam(norb)
   real (kind=RealKind), intent(in) :: jparam(norb)
!
   if (.not.Initialized) then
      call ErrorHandler('insertCorrectionOrbitals',                   &
                        'Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('insertCorrectionOrbitals','invalid atom index',id)
   else if (norb < 1) then
      call WarningHandler('insertCorrectionOrbitals',                 &
                          'Number of orbitals < 1',norb)
      return
   endif
!
   AtomList(id)%SpeciesList(ic)%NumCorrOrbitals = norb; n = 1
   allocate( AtomList(id)%SpeciesList(ic)%Orbitals(1:norb) ); n = n + norb
   allocate( AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-3:3,-3:3,1:2,1:norb) ); n = n + 7*7*2*norb
   allocate( AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(-3:3,1:2,1:norb) ); n = n + 7*2*norb*2
   allocate( AtomList(id)%SpeciesList(ic)%transform_inverse(1:7,1:7,1:2,1:norb) ); n = n + 7*7*2*norb*2
   allocate( AtomList(id)%SpeciesList(ic)%transform_matrix(1:7,1:7,1:2,1:norb) ); n = n + 7*7*2*norb*2
   allocate( AtomList(id)%SpeciesList(ic)%e_lda_plus_u(1:norb) ); n = n + norb
   allocate( AtomList(id)%SpeciesList(ic)%Uparam(1:norb), AtomList(id)%SpeciesList(ic)%Jparam(1:norb) ); n = n + 2*norb
   allocate( AtomList(id)%SpeciesList(ic)%SlaterIntegral(1:4,1:norb) ); n = n + 4*norb
   AtomList(id)%SpeciesList(ic)%Orbitals(1:norb) = lc(1:norb)
   DataPackSize(id) = n
!
   do i = 1, norb
      AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(:,:,:,i) = ZERO
      AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(:,:,i) = CZERO
      AtomList(id)%SpeciesList(ic)%transform_inverse(:,:,:,i) = CZERO
      AtomList(id)%SpeciesList(ic)%transform_matrix(:,:,:,i) = CZERO
      AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i) = ZERO
      AtomList(id)%SpeciesList(ic)%Uparam(i) = uparam(i)
      AtomList(id)%SpeciesList(ic)%Jparam(i) = jparam(i)
!     ----------------------------------------------------------------
      call computeSlaterIntegral(id,ic,i,lc(i),uparam(i),jparam(i))
!     ----------------------------------------------------------------
   enddo
!
   CorrectionNeeded = .true.
!
   end subroutine insert_withdata0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insert_withdata1(id,ic,norb,lc,uparam,jparam,slater)
!  ===================================================================
   use MathParamModule, only : ZERO, TEN2m6, CZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind), intent(in) :: norb
   integer (kind=IntKind), intent(in) :: lc(norb)
   integer (kind=IntKind) :: i, n
!
   real (kind=RealKind), intent(in) :: uparam(norb)
   real (kind=RealKind), intent(in) :: jparam(norb)
   real (kind=RealKind), intent(in) :: slater(3,norb)
!
   if (.not.Initialized) then
      call ErrorHandler('insertCorrectionOrbitals',                   &
                        'Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('insertCorrectionOrbitals','invalid atom index',id)
   else if (norb < 1) then
      call WarningHandler('insertCorrectionOrbitals',                 &
                          'Number of orbitals < 1',norb)
      return
   endif
!
   AtomList(id)%SpeciesList(ic)%NumCorrOrbitals = norb; n = 1
   allocate( AtomList(id)%SpeciesList(ic)%Orbitals(1:norb) ); n = n + norb
   allocate( AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-3:3,-3:3,1:2,1:norb) ); n = n + 7*7*2*norb
   allocate( AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(-3:3,1:2,1:norb) ); n = n + 7*2*norb*2
   allocate( AtomList(id)%SpeciesList(ic)%transform_inverse(1:7,1:7,1:2,1:norb) ); n = n + 7*7*2*norb*2
   allocate( AtomList(id)%SpeciesList(ic)%transform_matrix(1:7,1:7,1:2,1:norb) ); n = n + 7*7*2*norb*2
   allocate( AtomList(id)%SpeciesList(ic)%e_lda_plus_u(1:norb) ); n = n + norb
   allocate( AtomList(id)%SpeciesList(ic)%Uparam(1:norb), AtomList(id)%SpeciesList(ic)%Jparam(1:norb) ); n = n + 2*norb
   allocate( AtomList(id)%SpeciesList(ic)%SlaterIntegral(1:4,1:norb) ); n = n + 4*norb
   AtomList(id)%SpeciesList(ic)%Orbitals(1:norb) = lc(1:norb)
   DataPackSize(id) = n
!
   do i = 1, norb
      AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(:,:,:,i) = ZERO
      AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(:,:,i) = CZERO
      AtomList(id)%SpeciesList(ic)%transform_inverse(:,:,:,i) = CZERO
      AtomList(id)%SpeciesList(ic)%transform_matrix(:,:,:,i) = CZERO
      AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i) = ZERO
      AtomList(id)%SpeciesList(ic)%Uparam(i) = uparam(i)
      AtomList(id)%SpeciesList(ic)%Jparam(i) = jparam(i)
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(1,i) = uparam(i)
      if ( abs(slater(1,i)) > TEN2m6 .or. abs(slater(2,i)) > TEN2m6 .or. &
           abs(slater(3,i)) > TEN2m6 ) then
         AtomList(id)%SpeciesList(ic)%SlaterIntegral(2:4,i) = slater(1:3,i)
      else
!        -------------------------------------------------------------
         call computeSlaterIntegral(id,ic,i,lc(i),uparam(i),jparam(i))
!        -------------------------------------------------------------
      endif
   enddo
!
   CorrectionNeeded = .true.
!
   end subroutine insert_withdata1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insert_withfile(ujfile)
!  ===================================================================
   use MPPModule, only : MyPE, syncAllPEs, bcastMessage
!
   use MathParamModule, only : ZERO
!
   use ChemElementModule, only : getZtot
!
   use StringModule, only : initString, setString, endString,         &
                            getNumTokens, readToken
!
   implicit none
!
   logical :: file_exist
!
   character (len=*), intent(in) :: ujfile
   character (len=1) :: dummy
   character (len=25) :: s
   character (len=80) :: text
   character (len=80), allocatable :: ctemp(:)
!
   integer (kind=IntKind), parameter :: fu = 102
   integer (kind=IntKind) :: lc(MaxOrbitals,LocalNumAtoms)
   integer (kind=IntKind) :: atom_id(LocalNumAtoms)
   integer (kind=IntKind) :: num_orb(LocalNumAtoms)
   integer (kind=IntKind) :: i, j, ic, id, jd, nc, norb, nt, nl, n, slen
   integer (kind=IntKind) :: anum, gindex
   integer (kind=IntKind) :: status, StartingLine
!
   real (kind=RealKind) :: uparam(MaxOrbitals,LocalNumAtoms)
   real (kind=RealKind) :: jparam(MaxOrbitals,LocalNumAtoms)
   real (kind=RealKind) :: slater(3,MaxOrbitals,LocalNumAtoms)
!
   interface
      function isNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isNumber
   end interface
!
   interface
      function isRealNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isRealNumber
   end interface
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('insertCorrectionOrbitals',                   &
                        'Needs to call initLdaCorrection first')
   else if ( nocaseCompare(ujfile,'None') ) then
      call ErrorHandler('insertCorrectionOrbitals','File does not exist',ujfile)
   endif
!
   inquire(file=ujfile,exist=file_exist)
   if (.not.file_exist) then
      call ErrorHandler('insertCorrectionOrbitals','File does not exist',ujfile)
   endif
!
   slater(:,:,:) = ZERO
!
!  -------------------------------------------------------------------
   call initString(80)
!  -------------------------------------------------------------------
!
   nt = 0
   StartingLine = 1
   if (MyPE == 0) then
      open(unit=fu,file=ujfile,form='formatted',status='old')
      LOOP_do: do
         read(fu,'(a)',iostat=status)text
         if (status < 0) then
            call WarningHandler('insertCorrectionOrbitals','Invalid data',ujfile)
            exit LOOP_do
         else
            StartingLine = StartingLine + 1
         endif
         text=adjustl(text)
         if (text(1:1) == '#' .or. text(1:1) == '!'                   &
                              .or. len_trim(text) == 0) then
            cycle LOOP_do
         else
!           ----------------------------------------------------------
            call setString(text)
!           ----------------------------------------------------------
            if (getNumTokens() == 1) then   ! the first line contains 1 number
               if (.not.isNumber(text) .or. isRealNumber(text)) then
                  call WarningHandler('readPositionData','Invalid line',text)
                  exit LOOP_do
               else
                  read(text,*)nt
                  exit LOOP_do
               endif
            else
               call WarningHandler('insertCorrectionOrbitals',        &
                                   'Invalid data',text)
            endif
         endif
      enddo LOOP_do
   endif
!  -------------------------------------------------------------------
   call bcastMessage(nt,0)
!  -------------------------------------------------------------------
   if (nt < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('insertCorrectionOrbitals','Fail to read data',ujfile)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
   allocate( ctemp(nt) )
!  -------------------------------------------------------------------
   if (MyPE == 0) then
      rewind(fu)
      do i = 1, StartingLine-1
         read(fu,'(a)')dummy
      enddo
      nl = 0
      LOOP_do2: do
         read(fu,'(a)',iostat=status)text
         if (status < 0) then
            call ErrorHandler('insertCorrectionOrbitals','Fail to read',ujfile)
         endif
         text=adjustl(text)
         if (text(1:1) == '#' .or. text(1:1) == '!'                   &
                              .or. len_trim(text) == 0) then
            cycle LOOP_do2
         else
            nl = nl + 1
            ctemp(nl) = text 
         endif
         if (nl == nt) then
            exit LOOP_do2
         endif
      enddo LOOP_do2
      close(unit=fu)
   endif
!  -------------------------------------------------------------------
   call bcastMessage(ctemp,nt,0)
!  -------------------------------------------------------------------
!
   nc = 0
   num_orb(1:LocalNumAtoms) = 0
   atom_id(1:LocalNumAtoms) = 0
   lc(1:MaxOrbitals,1:LocalNumAtoms) = -1
   uparam(1:MaxOrbitals,1:LocalNumAtoms) = ZERO
   jparam(1:MaxOrbitals,1:LocalNumAtoms) = ZERO
   do i = 1, nt
!     ----------------------------------------------------------------
      call setString(ctemp(i))
!     ----------------------------------------------------------------
      n = getNumTokens()
      if (n < 5 .or. n > 8) then
!        -------------------------------------------------------------
         call ErrorHandler('insertCorrectionOrbitals',                &
                           'Invalid data format',ctemp(i))
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
      call readToken(1,s,slen)
!     ----------------------------------------------------------------
      if (isNumber(s)) then
         read(s,*)gindex
      else
!        -------------------------------------------------------------
         call ErrorHandler('insertCorrectionOrbitals',                &
                           'Invalid data format',ctemp(i))
!        -------------------------------------------------------------
      endif
!     ----------------------------------------------------------------
      call readToken(2,s,slen)
!     ----------------------------------------------------------------
      if (isNumber(s)) then
         read(s,*)anum
      else
         anum = getZtot(s(1:slen))
      endif
      LOOP_id: do id = 1, LocalNumAtoms
         if ( anum == AtomList(id)%AtomicNumber(1) .and.                          &
             (gindex < 1 .or. gindex == GlobalIndex(id)) ) then
            jd = 0
            LOOP_j: do j = 1, nc
               if (atom_id(j) == id) then
                  jd = j
                  exit LOOP_j
               endif
            enddo LOOP_j
            if (jd > 0) then
               num_orb(id) = num_orb(id) + 1
            else
               nc = nc + 1
               atom_id(nc) = id
               num_orb(id) = 1
            endif
!           ----------------------------------------------------------
            call readToken(3,s,slen)
!           ----------------------------------------------------------
            if (isNumber(s)) then
               read(s,*)lc(num_orb(id),id)
            else
!              -------------------------------------------------------
               call ErrorHandler('insertCorrectionOrbitals',          &
                                 'Invalid data format',ctemp(i))
!              -------------------------------------------------------
            endif
            if (lc(num_orb(id),id) < 2 .or. lc(num_orb(id),id) > 3) then
!              -------------------------------------------------------
               call ErrorHandler('insertCorrectionOrbitals',          &
                                 'Invalid l value',lc(num_orb(id),id))
!              -------------------------------------------------------
            else if (lc(num_orb(id),id) == 2 .and. n > 5 .and. n /= 7) then
!              -------------------------------------------------------
               call ErrorHandler('insertCorrectionOrbitals',          &
                                 'Invalid data format',ctemp(i))
!              -------------------------------------------------------
            else if (lc(num_orb(id),id) == 3 .and. n > 5 .and. n /= 8) then
!              -------------------------------------------------------
               call ErrorHandler('insertCorrectionOrbitals',          &
                                 'Invalid data format',ctemp(i))
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            call readToken(4,s,slen)
!           ----------------------------------------------------------
            if (isNumber(s)) then
               read(s,*)uparam(num_orb(id),id)
            else
!              -------------------------------------------------------
               call ErrorHandler('insertCorrectionOrbitals',          &
                                 'Invalid data format',ctemp(i))
!              -------------------------------------------------------
            endif
!           ----------------------------------------------------------
            call readToken(5,s,slen)
!           ----------------------------------------------------------
            if (isNumber(s)) then
               read(s,*)jparam(num_orb(id),id)
            else
!              -------------------------------------------------------
               call ErrorHandler('insertCorrectionOrbitals',          &
                                 'Invalid data format',ctemp(i))
!              -------------------------------------------------------
            endif
!
            do j = 6, n
!              -------------------------------------------------------
               call readToken(j,s,slen)
!              -------------------------------------------------------
               if (isNumber(s)) then
                  read(s,*)slater(j-5,num_orb(id),id)
               else
!                 ----------------------------------------------------
                  call ErrorHandler('insertCorrectionOrbitals',          &
                                    'Invalid data format',ctemp(i))
!                 ----------------------------------------------------
               endif
            enddo
         endif
      enddo LOOP_id
   enddo
!  -------------------------------------------------------------------
   call endString()
!  -------------------------------------------------------------------
   deallocate( ctemp )
!
!  This piece of codes need to be fixed for the multi-component case
   do ic = 1, nc
      id = atom_id(ic)
      norb = num_orb(id)
!     ----------------------------------------------------------------
      call insert_withdata1(id,1,norb,lc(1:norb,id),                 &
                            uparam(1:norb,id),jparam(1:norb,id),      &
                            slater(1:3,1:norb,id))
!     ----------------------------------------------------------------
   enddo
!
!  -------------------------------------------------------------------
   call syncAllPEs()
!  -------------------------------------------------------------------
!
   end subroutine insert_withfile
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeLdaPlusU(id,ic,kmax,dm)
!  ===================================================================
   use MathParamModule, only : ZERO, HALF, ONE, TEN2m8
   use MatrixInverseModule, only : MtxInv_GE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic, kmax
   integer (kind=IntKind) :: i, k, n, is, isp
   integer (kind=IntKind) :: lc, klp2, klp3, kn, kl, klp
   integer (kind=IntKind) :: m, mp, mp2, mp3, j, jp
!
   real (kind=RealKind), intent(in) :: dm(1:kmax,1:kmax,2)
   real (kind=RealKind) :: vee(-3:3,-3:3,-3:3,-3:3)
   real (kind=RealKind) :: vee_sym1(-3:3,-3:3,-3:3,-3:3)
   real (kind=RealKind) :: vee_sym2(-3:3,-3:3,-3:3,-3:3)
   real (kind=RealKind) :: occ(2), tocc, nmm
!
   integer (kind=IntKind), parameter :: lda = 7
   integer (kind=IntKind) :: info
   real (kind=RealKind) :: vm(lda,lda), wr(lda), wi(lda)
   real (kind=RealKind) :: vl(lda,lda), vr(lda,lda), work(4*lda)
   complex (kind=CmplxKind) :: tr2(5,5), tr2i(5,5)
   complex (kind=CmplxKind) :: tr3(7,7), tr3i(7,7)
!
   do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      lc = AtomList(id)%SpeciesList(ic)%Orbitals(i)
      kn = lc*lc + lc + 1
!
      if (lc < 2 .or. lc > 3) then
!        -------------------------------------------------------------
         call ErrorHandler('computeLdaPlusU','l is not of d- or f-orbital',lc)
!        -------------------------------------------------------------
      endif
!
      vee(:,:,:,:) = ZERO
      do k = 1, lc+1
!        =============================================================
!        ak(m,mp,mp2,mp3,lc,k) 
!                = a_k(m,m',m'',m'''), defined in PRB 52, R5467 (1995)
!        =============================================================
         vee(-lc:lc,-lc:lc,-lc:lc,-lc:lc) = vee(-lc:lc,-lc:lc,-lc:lc,-lc:lc) &
                                   + ak(-lc:lc,-lc:lc,-lc:lc,-lc:lc,lc,k)    &
                                     *AtomList(id)%SpeciesList(ic)%SlaterIntegral(k,i)
      enddo
!
!     ================================================================
!     Symmetrize vee so that the LDA+U potential is symmetric
!     ================================================================
      do mp3 = -lc, lc
         do mp2 = -lc, lc
            do mp = -lc, lc
               do m = -lc, lc
                  vee_sym1(m,mp,mp2,mp3) = HALF*( vee(m,mp,mp2,mp3)   &
                                                 +vee(mp,m,mp2,mp3))
                  vee_sym2(m,mp,mp2,mp3) = HALF*( vee(m,mp,mp2,mp3)   &
                                                 +vee(mp3,mp,mp2,m))
               enddo
            enddo
         enddo
      enddo
!
      do is = 1, 2
         occ(is) = ZERO
         do m = -lc, lc
            kl = kn + m
            occ(is) = occ(is) + dm(kl,kl,is)
         enddo
      enddo
      tocc = occ(1) + occ(2)
!
      AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i) =                                  &
           - HALF*AtomList(id)%SpeciesList(ic)%Uparam(i)*(tocc-ONE)*tocc              &
           + HALF*AtomList(id)%SpeciesList(ic)%Jparam(i)*(  occ(1)*(occ(1)-ONE)       &
                                          + occ(2)*(occ(2)-ONE) )
!
      do is = 1, 2
         isp = 3-is
         AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-lc:lc,-lc:lc,is,i) =         &
                               -AtomList(id)%SpeciesList(ic)%Uparam(i)*(tocc-HALF)    &
                               +AtomList(id)%SpeciesList(ic)%Jparam(i)*(occ(is)-HALF)
         do mp3 = -lc, lc
            klp3 = kn + mp3
            do mp2 = -lc, lc
               klp2 = kn + mp2
!              =======================================================
!              vee(m,mp,mp2,mp3) = <m,m''|V_ee|m',m'''>, defined in 
!                                                        PRB 52, R5467 (1995)
!              =======================================================
               nmm = dm(klp2,klp3,is)+dm(klp2,klp3,isp)
               AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-lc:lc,-lc:lc,is,i) =   &
                  AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-lc:lc,-lc:lc,is,i)  &
                                     + vee_sym1(-lc:lc,-lc:lc,mp2,mp3)*nmm
!                                    + vee(-lc:lc,-lc:lc,mp2,mp3)*nmm
               do mp = -lc, lc
                  klp = kn + mp
                  do m = -lc, lc
                     kl = kn + m
                     AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i) =                   &
                           AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i)               &
                                + vee_sym1(m,mp,mp2,mp3)*dm(kl,klp,is)*nmm
!                               + vee(m,mp,mp2,mp3)*dm(kl,klp,is)*nmm
                  enddo
               enddo
            enddo
         enddo
         do mp = -lc, lc
            klp = kn + mp
            do mp2 = -lc, lc
               klp2 = kn + mp2
               do mp3 = -lc, lc
                  klp3 = kn + mp3
                  AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-lc:lc,mp,is,i) =    &
                     AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-lc:lc,mp,is,i)   &
                              - vee_sym2(-lc:lc,mp3,mp2,mp)*dm(klp3,klp2,is)
!                                  - vee(-lc:lc,mp3,mp2,mp)*dm(klp3,klp2,is)
                  do m = -lc, lc
                     kl = kn + m
                     AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i) =                   &
                           AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i)               &
                                   - vee_sym2(m,mp3,mp2,mp)           &
!                                  - vee(m,mp3,mp2,mp)                &
                                         *dm(klp3,klp2,is)*dm(kl,klp,is)
                  enddo
               enddo
            enddo
         enddo
!
!        =============================================================
!        Diagonalize the LDA+U potential
!        =============================================================
         do mp = -lc, lc
            jp = mp + lc + 1
            do m = -lc, lc
               j = m + lc + 1
               vm(j,jp) = AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(m,mp,is,i)
            enddo
         enddo
         n = 2*lc + 1
!        -------------------------------------------------------------
         call dgeev('N','V',n,vm,lda,wr,wi,vl,lda,vr,lda,work,4*lda,info)
!        -------------------------------------------------------------
         do m = -lc, lc
            j = m + lc + 1
            AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(m,is,i) =                  &
                                          cmplx(wr(j),wi(j),kind=CmplxKind)
         enddo
         jp = 1    
         do while (jp <= n) 
            if (abs(wi(jp)) < TEN2m8) then 
               do j = 1, n  
                  AtomList(id)%SpeciesList(ic)%transform_matrix(j,jp,is,i) =           &
                                          cmplx(vr(j,jp),ZERO,kind=CmplxKind)
               enddo
               jp = jp + 1
            else
               do j = 1, n
                  AtomList(id)%SpeciesList(ic)%transform_matrix(j,jp,is,i) =           &
                                     cmplx(vr(j,jp), vr(j,jp+1),kind=CmplxKind)
                  AtomList(id)%SpeciesList(ic)%transform_matrix(j,jp+1,is,i) =         &
                                     cmplx(vr(j,jp),-vr(j,jp+1),kind=CmplxKind)
               enddo
               jp = jp + 2
            endif 
         enddo
         if (lc == 2) then
            do j = 1, 5
               tr2(1:5,j) = AtomList(id)%SpeciesList(ic)%transform_matrix(1:5,j,is,i)
            enddo
!           ----------------------------------------------------------
            call MtxInv_GE(5,tr2,tr2i)
!           ----------------------------------------------------------
            do j = 1, 5
               AtomList(id)%SpeciesList(ic)%transform_inverse(1:5,j,is,i) = tr2i(1:5,j)
            enddo
         else if (lc == 3) then
            do j = 1, 7
               tr3(1:7,j) = AtomList(id)%SpeciesList(ic)%transform_matrix(1:7,j,is,i)
 write(6,'(a,i5,14d15.8)')'i,j,tr3 = ',j,tr3(1:7,j)
            enddo
!           ----------------------------------------------------------
            call MtxInv_GE(7,tr3,tr3i)
!           ----------------------------------------------------------
            do j = 1, 7
               AtomList(id)%SpeciesList(ic)%transform_inverse(1:7,j,is,i) = tr3i(1:7,j)
            enddo
         else
!           ----------------------------------------------------------
            call ErrorHandler('computeLdaPlusU','the orbital is not d or f',lc)
!           ----------------------------------------------------------
         endif
!
      enddo
   enddo
!
   end subroutine computeLdaPlusU
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeSlaterIntegral(id,ic,n,lc,uparam,jparam)
!  ===================================================================
   use MathParamModule, only : ZERO, ONE
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: lc
!
   real (kind=RealKind), intent(in) :: uparam
   real (kind=RealKind), intent(in) :: jparam
!
   real (kind=RealKind) :: F0, F2, F4, F6
!  real (kind=RealKind), parameter :: F4_over_F2 = 0.625d0
   real (kind=RealKind), parameter :: F4_over_F2 = 451.0d0/675.0d0
   real (kind=RealKind), parameter :: F6_over_F2 = 1001.0d0/2025.0d0
!
   AtomList(id)%SpeciesList(ic)%SlaterIntegral(1:4,n) = ZERO
   if (lc == 2) then
      F0 = uparam
      F2 = 14.0d0*jparam/(ONE+F4_over_F2)
      F4 = 14.0d0*jparam*F4_over_F2/(ONE+F4_over_F2)
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(1,n) = F0
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(2,n) = F2
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(3,n) = F4
   else if (lc == 3) then
      F0 = uparam
      F2 = 6435.0d0*jparam/(286.0d0+195.0d0*F4_over_F2+250.0d0*F6_over_F2)
      F4 = F2*F4_over_F2
      F6 = F2*F6_over_F2
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(1,n) = F0
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(2,n) = F2
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(3,n) = F4
      AtomList(id)%SpeciesList(ic)%SlaterIntegral(4,n) = F6
   else
      call ErrorHandler('computeSlaterIntegral',                      &
                        'This l value is not implemented',lc)
   endif
!
   end subroutine computeSlaterIntegral
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotentialCorrection(id,ic,is,lc,m) result(vcorr)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic, is, lc, m
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind) :: vcorr
!
   if (.not.Initialized) then
      call ErrorHandler('getPotentialCorrection',                     &
                        'Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getPotentialCorrection','invalid atom index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('getPotentialCorrection','invalid spin index',is)
   else if (lc < 0) then
      call ErrorHandler('getPotentialCorrection','invalid orbital index',lc)
   else if (m < -lc .or. m > lc) then
      call ErrorHandler('getPotentialCorrection','invalid m-index',m)
   endif
!
   LOOP_i: do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      if (AtomList(id)%SpeciesList(ic)%Orbitals(i) == lc) then
         vcorr = AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(m,is,i)
         exit LOOP_i
      endif
   enddo LOOP_i
!
   end function getPotentialCorrection
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTransformMatrixInverse(id,ic,is,lc) result(ul)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic, is, lc
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), pointer :: ul(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getTransformMatrixInverse',                     &
                        'Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTransformMatrixInverse','invalid atom index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('getTransformMatrixInverse','invalid spin index',is)
   else if (lc < 0) then
      call ErrorHandler('getTransformMatrixInverse','invalid orbital index',lc)
   endif
!
   nullify( ul )
   LOOP_i: do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      if (AtomList(id)%SpeciesList(ic)%Orbitals(i) == lc) then
         ul => AtomList(id)%SpeciesList(ic)%transform_inverse(:,:,is,i)
         exit LOOP_i
      endif
   enddo LOOP_i
!
   end function getTransformMatrixInverse
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTransformMatrix(id,ic,is,lc) result(ur)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic, is, lc
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), pointer :: ur(:,:)
!
   if (.not.Initialized) then
      call ErrorHandler('getTransformMatrix',                    &
                        'Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getTransformMatrix','invalid atom index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('getTransformMatrix','invalid spin index',is)
   else if (lc < 0) then
      call ErrorHandler('getTransformMatrix','invalid orbital index',lc)
   endif
!
   nullify( ur )
   LOOP_i: do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      if (AtomList(id)%SpeciesList(ic)%Orbitals(i) == lc) then
         ur => AtomList(id)%SpeciesList(ic)%transform_matrix(:,:,is,i)
         exit LOOP_i
      endif
   enddo LOOP_i
!
   end function getTransformMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEnergyCorrection_al(id,ic,lc) result(ecorr)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic, lc
   integer (kind=IntKind) :: i
!
   real (kind=RealKind) :: ecorr
!
   if (.not.Initialized) then
      call ErrorHandler('getEnergyCorrection','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getEnergyCorrection','invalid atom index',id)
   else if (lc < 0) then
      call ErrorHandler('getEnergyCorrection','invalid orbital index',lc)
   endif
!
   ecorr = ZERO
   LOOP_i: do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      if (AtomList(id)%SpeciesList(ic)%Orbitals(i) == lc) then
         ecorr = AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i)
         exit LOOP_i
      endif
   enddo LOOP_i
!
   end function getEnergyCorrection_al
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEnergyCorrection_a(id,ic) result(ecorr)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind) :: i
!
   real (kind=RealKind) :: ecorr
!
   if (.not.Initialized) then
      call ErrorHandler('getEnergyCorrection','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getEnergyCorrection','invalid atom index',id)
   endif
!
   ecorr = ZERO
   LOOP_i: do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      ecorr = ecorr + AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i)
   enddo LOOP_i
!
   end function getEnergyCorrection_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine transformScatteringMatrix(id,ic,is,mtx)
!  ===================================================================
   use MathParamModule, only : CZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind) :: i, l, kl, msize
   integer (kind=IntKind) :: kl1, kl2, kl3
!
   complex (kind=CmplxKind), intent(inout) :: mtx(:,:)
   complex (kind=CmplxKind) :: utm(7), mut(7,7)
!
   if (.not.Initialized) then
      call ErrorHandler('transformScatteringMatrix','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('transformScatteringMatrix','invalid atom index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('transformScatteringMatrix','invalid spin index',is)
   endif
!
   do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      l = AtomList(id)%SpeciesList(ic)%Orbitals(i)
      msize = 2*l+1
      kl = l*l
      mut(:,:) = CZERO
      do kl1 = 1, msize
         utm(:) = CZERO
         do kl2 = 1, msize
            do kl3 = 1, msize
               utm(kl2) = utm(kl2) +                                  &
                          AtomList(id)%SpeciesList(ic)%transform_inverse(kl3,kl2,is,i)*mtx(kl+kl3,kl+kl1)
            enddo
         enddo
         do kl2 = 1, msize
            do kl3 = 1, msize
               mut(kl3,kl2) = mut(kl3,kl2) +                          &
                          utm(kl3)*AtomList(id)%SpeciesList(ic)%transform_matrix(kl2,kl1,is,i)
            enddo
         enddo
      enddo
      do kl1 = 1, msize
         do kl2 = 1, msize
            mtx(kl+kl2,kl+kl1) = mut(kl2,kl1)
         enddo
      enddo
   enddo
!
   end subroutine transformScatteringMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine transformWF2(id,ic,is,nr,wave)
!  ===================================================================
   use MathParamModule, only : CZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind) :: i, l, kl, msize
   integer (kind=IntKind) :: kl1, kl2
!
   complex (kind=CmplxKind), intent(inout) :: wave(:,:)
   complex (kind=CmplxKind) :: wt(nr,7)
!
   if (.not.Initialized) then
      call ErrorHandler('transformWaveFunction','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('transformWaveFunction','invalid atom index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('transformWaveFunction','invalid spin index',is)
   endif
!
   do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      l = AtomList(id)%SpeciesList(ic)%Orbitals(i)
      msize = 2*l+1
      kl = l*l
      wt(:,:) = CZERO
      do kl1 = 1, msize
         do kl2 = 1, msize
            wt(1:nr,kl2) = wt(1:nr,kl2) +               &
                          AtomList(id)%SpeciesList(ic)%transform_matrix(kl2,kl1,is,i)*wave(1:nr,kl+kl1)
         enddo
      enddo
      do kl2 = 1, msize
         wave(1:nr,kl+kl2) = wt(1:nr,kl2)
      enddo
   enddo
!
   end subroutine transformWF2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine transformWF3(id,ic,is,nr,kmax,wave)
!  ===================================================================
   use MathParamModule, only : CZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind), intent(in) :: is
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: kmax
   integer (kind=IntKind) :: i, l, kl, msize
   integer (kind=IntKind) :: kl1, kl2
!
   complex (kind=CmplxKind), intent(inout) :: wave(:,:,:)
   complex (kind=CmplxKind) :: wt(nr,kmax,7)
!
   if (.not.Initialized) then
      call ErrorHandler('transformWaveFunction','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('transformWaveFunction','invalid atom index',id)
   else if (is < 1 .or. is > 2) then
      call ErrorHandler('transformWaveFunction','invalid spin index',is)
   endif
!
   do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      l = AtomList(id)%SpeciesList(ic)%Orbitals(i)
      msize = 2*l+1
      kl = l*l
      wt(:,:,:) = CZERO
      do kl1 = 1, msize
         do kl2 = 1, msize
            wt(1:nr,1:kmax,kl2) = wt(1:nr,1:kmax,kl2) +               &
                          AtomList(id)%SpeciesList(ic)%transform_matrix(kl2,kl1,is,i)*wave(1:nr,1:kmax,kl+kl1)
         enddo
      enddo
      do kl2 = 1, msize
         wave(1:nr,1:kmax,kl+kl2) = wt(1:nr,1:kmax,kl2)
      enddo
   enddo
!
   end subroutine transformWF3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCorrOrbitals(id,ic) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumCorrOrbitals','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getNumCorrOrbitals','invalid site index',id)
   else if (ic < 1 .or. ic > AtomList(id)%NumSpecies) then
      call ErrorHandler('getNumCorrOrbitals','invalid species index',ic)
   endif
!
   n = AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
!
   end function getNumCorrOrbitals
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCorrOrbital(id,ic,i) result(l)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic, i
   integer (kind=IntKind) :: l
!
   if (.not.Initialized) then
      call ErrorHandler('getCorrOrbital','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getCorrOrbital','invalid atom index',id)
   else if (i < 0 .or. i > AtomList(id)%SpeciesList(ic)%NumCorrOrbitals) then
      call ErrorHandler('getCorrOrbital','invalid orbital index',i)
   endif
!
   if (i > 0) then
      l = AtomList(id)%SpeciesList(ic)%Orbitals(i)
   else
      l = -1 
   endif
!
   end function getCorrOrbital
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDataPackSize(id) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getDataPackSize','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDataPackSize','invalid atom index',id)
   endif
!
   n = DataPackSize(id)
!
   end function getDataPackSize
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDataPack4Output(id,ic,n) result(pk)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic
   integer (kind=IntKind), intent(out) :: n
   integer (kind=IntKind) :: i, j
!
   real (kind=RealKind), pointer :: pk(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getDataPack4Output','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('getDataPack4Output','invalid atom index',id)
   else if (AtomList(id)%SpeciesList(ic)%NumCorrOrbitals < 1) then
      nullify(pk)
      return
   endif
!
   if (MaxDataPackSize < 1) then
      MaxDataPackSize = maxval( DataPackSize(1:LocalNumAtoms) )
   endif
!
   if (MaxDataPackSize < 1) then
      nullify(pk)
      return
   else if (DataPackSize(id) < 1) then
      call ErrorHandler('getDataPack4Output','invalid data pack size',DataPackSize(id))
   else if (.not.allocated(DataPack)) then
      allocate( DataPack(1:MaxDataPackSize) )
   endif
!
   n = DataPackSize(id)
!
   if (n /= 640*AtomList(id)%SpeciesList(ic)%NumCorrOrbitals) then
      call ErrorHandler('getDataPack4Output','data pack size <> 640*num_orbitals',n)
   endif
!
   j = 1
   DataPack(j) = AtomList(id)%SpeciesList(ic)%NumCorrOrbitals + 0.0001
   do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      j = j + 1
      DataPack(j) = AtomList(id)%SpeciesList(ic)%Orbitals(i) + 0.0001
      j = j + 1
      DataPack(j) = AtomList(id)%SpeciesList(ic)%Uparam(i)
      j = j + 1
      DataPack(j) = AtomList(id)%SpeciesList(ic)%Jparam(i)
      DataPack(j+1:j+4) = AtomList(id)%SpeciesList(ic)%SlaterIntegral(1:4,i)
      j = j + 5
      DataPack(j:j+97) = aliasArray1_r(AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-3:3,-3:3,1:2,i),98)
      j = j + 98
      call copyCmplx2Real(AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(-3:3,1:2,i),DataPack(j:j+27),28)
      j = j + 28
      call copyCmplx2Real(AtomList(id)%SpeciesList(ic)%transform_inverse(1:7,1:7,1:2,i),DataPack(j:j+195),196)
      j = j + 196
      call copyCmplx2Real(AtomList(id)%SpeciesList(ic)%transform_matrix(1:7,1:7,1:2,i),DataPack(j:j+195),196)
      j = j + 196
      DataPack(j) = AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i)
   enddo
   pk => DataPack(1:n)
!
   end function getDataPack4Output
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insertDataPackFromInput(id,ic,n,pk)
!  ===================================================================
   use MathParamModule, only : TEN2m6
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, ic, n
   integer (kind=IntKind) :: i, j
   real (kind=RealKind), intent(in) :: pk(:)
   real (kind=RealKind) :: diff(4)
!
   if (.not.Initialized) then
      call ErrorHandler('insertDataPackFromInput','Needs to call initLdaCorrection first')
   else if (id < 1 .or. id > LocalNumAtoms) then
      call ErrorHandler('insertDataPackFromInput','invalid atom index',id)
   else if (n < 1) then
      call ErrorHandler('insertDataPackFromInput','invalid pack size',n)
   endif
!
   AtomList(id)%SpeciesList(ic)%NumCorrOrbitals = floor(pk(1),kind=IntKind)
!
   if (n < 640*AtomList(id)%SpeciesList(ic)%NumCorrOrbitals) then
      call ErrorHandler('insertDataPackFromInput','data pack size < 640*num_orbitals',n)
   endif
!
   j = 1
   do i = 1, AtomList(id)%SpeciesList(ic)%NumCorrOrbitals
      j = j + 1
      if (AtomList(id)%SpeciesList(ic)%Orbitals(i) /= floor(pk(j),kind=IntKind)) then
         call ErrorHandler('insertDataPackFromInput','inconsistent orbitals', &
                           AtomList(id)%SpeciesList(ic)%Orbitals(i),floor(pk(j),kind=IntKind))
      endif
      j = j + 1
      if (abs(AtomList(id)%SpeciesList(ic)%Uparam(i)-pk(j)) > TEN2m6) then
         call WarningHandler('insertDataPackFromInput','U parameter has changed', &
                           AtomList(id)%SpeciesList(ic)%Uparam(i),pk(j))
      endif
      j = j + 1
      if (abs(AtomList(id)%SpeciesList(ic)%Jparam(i)-pk(j)) > TEN2m6) then
         call WarningHandler('insertDataPackFromInput','J parameter has changed', &
                           AtomList(id)%SpeciesList(ic)%Jparam(i),pk(j))
      endif
      diff(1:4) = AtomList(id)%SpeciesList(ic)%SlaterIntegral(1:4,i) - pk(j+1:j+4)
      if (abs(diff(1)) > TEN2m6 .or. abs(diff(2)) > TEN2m6 .or.                   &
          abs(diff(3)) > TEN2m6 .or. abs(diff(4)) > TEN2m6) then
         call WarningHandler('insertDataPackFromInput','Slater Integrals are changed')
      endif
      j = j + 5
      AtomList(id)%SpeciesList(ic)%v_lda_plus_u_full(-3:3,-3:3,1:2,i) = aliasArray3_r(DataPack(j:j+97),7,7,2)
      j = j + 98
      call copyReal2Cmplx(pk(j:j+27),AtomList(id)%SpeciesList(ic)%v_lda_plus_u_diag(-3:3,1:2,i),28)
      j = j + 28
      call copyReal2Cmplx(pk(j:j+195),AtomList(id)%SpeciesList(ic)%transform_inverse(1:7,1:7,1:2,i),196)
      j = j + 196
      call copyReal2Cmplx(pk(j:j+195),AtomList(id)%SpeciesList(ic)%transform_matrix(1:7,1:7,1:2,i),196)
      j = j + 196
      AtomList(id)%SpeciesList(ic)%e_lda_plus_u(i) = pk(j)
   enddo
!
   end subroutine insertDataPackFromInput
!  ===================================================================
end module LdaCorrectionModule
