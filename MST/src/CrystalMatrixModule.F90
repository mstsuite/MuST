!  *******************************************************************
!  * Module for calculating distributed Tau matrix                   *
!  * External Modules used:                                          *
!  *   KindParamModule                                               *
!  *   MathParamModule                                               *
!  *   WriteMatrixModule                                             *
!  *   ErrorHandlerModule                                            *
!  *   PhysParamModule                                               *
!  *   SingleScatteringModule                                        *
!  *                                                                 *
!  *Public functions:                                                *
!  *                                                                 *
!  *  initCrystalMatrix(nla, cant, lmax_kkr, rel, istop, ipr)        *
!  *  Purpose: Initialize the module for calculating Tau matrix      *
!  *  Input:   nla   = number of atoms allocated to this process     *
!  *           cant  = the interger defines whether it is a spin     *
!  *                   canted calculation                            *
!  *           lmax_kkr = largest possible l-cut off for single      *
!  *                      scattering t-matrix                        *
!  *           rel   = the interger defines whether it is a          *
!  *                   relativistic calculation                      *
!  *                    an 1-d array of LocalNumAtoms size           *
!  *           istop = routine name to stop (character string)       *
!  *           iprint= print level (integer)                         *
!  *                                                                 *
!  *  endCrystalMatrix()                                             *
!  *  Purpose: deallocate the internal allocated arrays and clear    *
!  *           the storage                                           *
!  *                                                                 *
!  *  calCrystalMatrix(e,getSM,tau,use_tmat,configuration)           *
!  *  Purpose: calculate k-space integration of the Tau matrix       *
!  *  Input:   e  = the energy in Rydberg units                      *
!  *           getSM (interface) = function interface to get single  *
!  *                               scattering matrix.                *
!  *           tau (logic) = whether to calculate tau-matrix....     *
!  *           use_tmat (logic) = whether to use t-matrix formula to *
!  *                              calculate tau                      *
!  *           configuration = if case of using effective medium,    *
!  *                           an atomic site may be probable to     *
!  *                           have different atomic species. This   *
!  *                           integer array specifies a particular  *
!  *                           configuration of atoms on the atomic  *
!  *                           sites mapped on the current process.  *
!  *           Note: for atomic site i, configuration(i) = 0 is the  *
!  *                 effective medium "atom"; configuration(i) = j   *
!  *                 means that the site has the j'th atomic species *
!  *                                                                 *
!  *  calCrystalMatrix(e,kpts,getSM,tauk,use_tmat,configuration)     *
!  *  Purpose: calculate the Tau matrix for the given k-point        *
!  *  Input:   e  = the energy in Rydberg units                      *
!  *           kpts = the array defines the k-points                 *
!  *           getSM (interface) = function interface to get single  *
!  !                               scattering matrix.                *
!  *           tau (logic) = whether to calculate tau-matrix....     *
!  *           use_tmat (logic) = whether to use t-matrix formula to *
!  *                              calculate tau                      *
!  *           configuration = if case of using effective medium,    *
!  *                           an atomic site may be probable to     *
!  *                           have different atomic species. This   *
!  *                           integer array specifies a particular  *
!  *                           configuration of atoms on the atomic  *
!  *                           sites mapped on the current process.  *
!  *           Note: for atomic site i, configuration(i) = 0 is the  *
!  *                 effective medium "atom"; configuration(i) = j   *
!  *                 means that the site has the j'th atomic species *
!  *                                                                 *
!  *  getTau(id,i,j)                                                 *
!  *  Purpose: retrieve the {i,j} block of the big Tau matrix        *
!  *  Input:   id  = the local  index of atom                        *
!  *           i   = the global index of atom i                      *
!  *           j   = the global index of atom j                      *
!  *  Result:  pointer tau => Tau^{i,j}, in the local frame of       *
!  *                            reference in spin space              *
!  *  Note:  the j-atom needs to be on the local process             *
!  *                                                                 *
!  *  getKau(id,i,j)                                                 *
!  *  Purpose: retrieve the {i,j} block of the big Kau matrix        *
!  *  Input:   id  = the local  index of atom                        *
!  *           i   = the global index of atom i                      *
!  *           j   = the global index of atom j                      *
!  *  Result:  pointer kau => Kau^{i,j}, in the local frame of       *
!  *                            reference in spin space              *
!  *  Note:  the j-atom needs to be on the local process             *
!  *                                                                 *
!  *  printCrystalMatrix(i,j)                                        *
!  *  Purpose: Prints out Kau^{i,j} matrix                           *
!  *  Input:   i   = the global index of atom i                      *
!  *           i   = the global index of atom j                      *
!  *  Note:  the j-atom needs to be on the local process             *
!  *                                                                 *
!  *******************************************************************
module CrystalMatrixModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1
   use MathParamModule, only : two, half, third
   use MathParamModule, only : ten2m8, ten2m6, ten2m14, ten2m9, ten2m10
   use MathParamModule, only : sqrtm1, pi2, pi
   use TimerModule, only: getTime
   use WriteMatrixModule, only  : writeMatrix
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler,        &
                                  StopHandler
   use PhysParamModule, only : LightSpeed !xianglin
   use PublicTypeDefinitionsModule, only : MatrixBandStruct, ScmBlockStruct
!
public :: initCrystalMatrix, &
          endCrystalMatrix,  &
          calCrystalMatrix,  &
          getTau,            &
          getKau,            &
          printCrystalMatrix
!
   interface calCrystalMatrix
      module procedure calCrystalMatrix_sumk , calCrystalMatrix_k
   end interface ! calCrystalMatrix
!
private
!
   type (ScmBlockStruct), allocatable :: sc_blocks(:,:)
!
   type (MatrixBandStruct), allocatable :: MatrixBand(:)  ! Matrix column band
!
   complex (kind=CmplxKind), allocatable, target :: Kau_MatrixDiag(:)
   complex (kind=CmplxKind), allocatable, target :: Tau_MatrixDiag(:)
   complex (kind=CmplxKind), allocatable, target :: TMP_MatrixBand(:)
   complex (kind=CmplxKind), allocatable, target :: KKR_MatrixBand(:)
   complex (kind=CmplxKind), allocatable, target :: SCM_MatrixBand(:)
   complex (kind=CmplxKind), allocatable, target :: strconrel(:,:)
!
   logical :: isRelativistic = .false.
   real (kind=RealKind), parameter :: Me = 0.5d0 !xianglin
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: max_print_level
!
   integer (kind=IntKind) :: GlobalNumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: nSpinCant
   integer (kind=IntKind) :: lmax_kkr_max, kmax_kkr_max, tsize
   integer (kind=IntKind) :: KKRMatrixSize
   integer (kind=IntKind) :: BandSize
   integer (kind=IntKind) :: KKRMatrixSizeCant
   integer (kind=IntKind) :: BandSizeCant
!
   integer (kind=IntKind), allocatable :: print_level(:)
   integer (kind=IntKind), allocatable :: ip_array(:) ! Relates a matrix block row index to the mapped processor index (0, 1, ...)
   integer (kind=IntKind), allocatable :: il_array(:) ! Relates a matrix block row index to the local atom index on the mapped processor
   integer (kind=IntKind), allocatable :: id_array(:) ! Relates a matrix block row index to the global index of the corresponding atom
   integer (kind=IntKind), allocatable :: jd_array(:) ! Relates a local atom index to the global index of the atom
   integer (kind=IntKind), allocatable :: gid_array(:)! Relates a global index to the corresponding matrix block row index
   integer (kind=IntKind), allocatable :: lmaxi_array(:)
   integer (kind=IntKind), allocatable :: lmaxj_array(:)
!
   integer (kind=IntKind) :: LWORK, LIWORK
   integer (kind=IntKind), allocatable :: IWORK(:)
   complex (kind=CmplxKInd), allocatable, target :: WORK(:)
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
!
   complex (kind=CmplxKind) :: energy, kappa
!
   complex (kind=CmplxKind), allocatable, target :: jinv_g(:,:) ! Jost matrix inverse
   complex (kind=CmplxKind), allocatable, target :: sine_g(:,:) ! Sine matrix in global frame
   complex (kind=CmplxKind), allocatable, target :: stcm_g(:,:) ! kappa*S^{T*}*C matrix in global frame
   complex (kind=CmplxKind), allocatable, target :: tmat_g(:,:) ! t-matrix in global frame
   complex (kind=CmplxKind), pointer :: cosine_g(:,:)  ! Cosine matrix in global frame
!
#ifdef USE_SCALAPACK
!  ===================================================================
!  *****      ScaLAPACK parameters
!  ===================================================================
   integer (kind=IntKind), parameter :: DLEN_ = 9
   integer (kind=IntKind) :: ICTXT
   integer (kind=IntKind) :: NPROW
   integer (kind=IntKind) :: NPCOL
   integer (kind=IntKind) :: MYROW
   integer (kind=IntKind) :: MYCOL
   integer (kind=IntKind) :: DESC_A( DLEN_ )
#endif
   integer (kind=IntKind) :: INFO
   integer (kind=IntKind), allocatable :: IPVT(:)
!  ===================================================================
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCrystalMatrix( nla, cant, lmax_kkr, rel, istop, iprint)
!  ===================================================================
   use MPPModule, only : MyPE
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
   use GroupCommModule, only : syncAllPEsInGroup
   use SystemModule, only  : getNumAtoms, getLmaxKKR, getAtomPosition, &
                             getLmaxPhi, getBravaisLattice
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms
   use StrConstModule, only : initStrConst
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: nla, cant, rel
   integer (kind=IntKind), intent(in) :: lmax_kkr(nla)
   integer (kind=IntKind), intent(in) :: iprint(nla)
!
   character (len=20) :: sname = "initCrystalMatrix"
!
   integer (kind=IntKind) :: i, j, ig, na, n, k, nk
   integer (kind=IntKind) :: lmaxi, kmaxi, kmaxj, t0size, nsize
   integer (kind=IntKind) :: status
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: global_posi(:,:)
!
   complex (kind=CmplxKind), pointer :: p1(:)
!
   stop_routine = istop
!
   GlobalNumAtoms = getNumAtoms()
   LocalNumAtoms = nla
   nSpinCant = cant
!
   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
!
   if (rel>1) then
      isRelativistic = .true.
      nSpinCant = 2  !fully relativistic override spin cant Xianglin
   endif
!
!  ===================================================================
!  Initialize the structure constant module.
!  ===================================================================
   allocate( global_posi(1:3,1:GlobalNumAtoms) )
   do i = 1, GlobalNumAtoms
      global_posi(1:3,i) = getAtomPosition(i)
   enddo
   bravais(1:3,1:3) = getBravaisLattice()
!  -------------------------------------------------------------------
   call initStrConst(getLmaxPhi(), GlobalNumAtoms, global_posi, bravais, &
                     istop, maxval(iprint))
!  -------------------------------------------------------------------
   deallocate(global_posi)
!  ===================================================================
!
   allocate ( print_level(LocalNumAtoms), STAT=status )
   allocate ( MatrixBand(LocalNumAtoms), STAT=status )
!
   allocate( id_array(GlobalNumAtoms), jd_array(LocalNumAtoms) )
   allocate( ip_array(GlobalNumAtoms), il_array(GlobalNumAtoms) )
   allocate( lmaxi_array(GlobalNumAtoms), lmaxj_array(LocalNumAtoms) )
   allocate( gid_array(GlobalNumAtoms), sc_blocks(GlobalNumAtoms,LocalNumAtoms) )
!
   max_print_level = -1
   BandSize = 0
   n = 0
   do j = 1, LocalNumAtoms
      print_level(j) = iprint(j)
      max_print_level = max(max_print_level,iprint(j))
      MatrixBand(j)%lmax_kkr = lmax_kkr(j)
      MatrixBand(j)%kmax_kkr = (lmax_kkr(j)+1)**2
      MatrixBand(j)%global_index = getGlobalIndex(j)
      BandSize = BandSize + MatrixBand(j)%kmax_kkr
      allocate( MatrixBand(j)%MatrixBlock(GlobalNumAtoms) )
      MatrixBand(j)%column_index = n + 1
      n = n + MatrixBand(j)%kmax_kkr*nSpinCant
      jd_array(j) = MatrixBand(j)%global_index
      lmaxj_array(j) = lmax_kkr(j)
   enddo
!
   lmax_kkr_max = 0
   KKRMatrixSize = 0
   do i = 1, GlobalNumAtoms
      lmaxi = getLmaxKKR(i)
      kmaxi = (lmaxi+1)**2
      lmax_kkr_max = max( lmax_kkr_max,lmaxi )
      KKRMatrixSize = KKRMatrixSize + kmaxi
   enddo
   kmax_kkr_max = (lmax_kkr_max+1)**2
!
   KKRMatrixSizeCant = KKRMatrixSize*nSpinCant
   BandSizeCant = BandSize*nSpinCant
   tsize = kmax_kkr_max*kmax_kkr_max*nSpinCant*nSpinCant
!
   allocate ( Kau_MatrixDiag(tsize*LocalNumAtoms) )
   allocate ( Tau_MatrixDiag(tsize*LocalNumAtoms) )
   allocate ( TMP_MatrixBand(KKRMatrixSizeCant*BandSizeCant) )
   allocate ( KKR_MatrixBand(KKRMatrixSizeCant*BandSizeCant) )
   allocate ( SCM_MatrixBand(KKRMatrixSize*BandSize) )
   if (isRelativistic) then
      allocate( strconrel(kmax_kkr_max*nSpinCant,kmax_kkr_max*nSpinCant) )
   endif
!
   LWORK = KKRMatrixSizeCant*BandSizeCant
   LIWORK = 2*KKRMatrixSizeCant
   allocate( WORK(1:max(LWORK,2*BandSizeCant*BandSizeCant)), IWORK(1:LIWORK) )
   WORK = CZERO; IWORK = 0
!
!  ===================================================================
!  Set up the matrix block row/column index in such way that each
!  a band of columns of the matrix is stored on a processor 
!      n = the matrix block row/column index in the band matrix
!     ig = the global index of the atom that corresponds to the matrix block
!  ===================================================================
   n = 0
   nk = 0
   nsize = 0
   do k = 1, NumPEsInGroup
      na = getLocalNumAtoms(k-1)
      do i = 1, na
         n = n + 1
         ig = getGlobalIndex(i,k-1)
         lmaxi = getLmaxKKR(ig)
         kmaxi = (lmaxi+1)**2
         t0size = kmaxi*kmaxi*nSpinCant*nSpinCant
         do j = 1, LocalNumAtoms
            MatrixBand(j)%MatrixBlock(n)%lmax_kkr = lmaxi
            MatrixBand(j)%MatrixBlock(n)%kmax_kkr = kmaxi
            MatrixBand(j)%MatrixBlock(n)%global_index = ig
            MatrixBand(j)%MatrixBlock(n)%row_index = nk + 1
            if (ig == getGlobalIndex(j)) then   ! For now, we only store the diagonal blocks of the big Kau matrix.
               if (isRelativistic) then !by xianglin
!                 MatrixBand(j)%MatrixBlock(n)%kau_l => aliasArray3_c(Kau_MatrixDiag(nsize+1:nsize+t0size), &
!                                                                     kmaxi*nSpinCant,kmaxi*nSpinCant,1)
                  p1 => Kau_MatrixDiag(nsize+1:nsize+t0size)
                  MatrixBand(j)%MatrixBlock(n)%kau_l => aliasArray3_c(p1,kmaxi*nSpinCant,kmaxi*nSpinCant,1)
!                 MatrixBand(j)%MatrixBlock(n)%tau_l => aliasArray3_c(Tau_MatrixDiag(nsize+1:nsize+t0size), &
!                                                                     kmaxi*nSpinCant,kmaxi*nSpinCant,1)
                  p1 => Tau_MatrixDiag(nsize+1:nsize+t0size)
                  MatrixBand(j)%MatrixBlock(n)%tau_l => aliasArray3_c(p1,kmaxi*nSpinCant,kmaxi*nSpinCant,1)
               else
!                 MatrixBand(j)%MatrixBlock(n)%kau_l => aliasArray3_c(Kau_MatrixDiag(nsize+1:nsize+t0size), &
!                                                                     kmaxi,kmaxi,nSpinCant*nSpinCant)
                  p1 => Kau_MatrixDiag(nsize+1:nsize+t0size)
                  MatrixBand(j)%MatrixBlock(n)%kau_l => aliasArray3_c(p1,kmaxi,kmaxi,nSpinCant*nSpinCant)
!                 MatrixBand(j)%MatrixBlock(n)%tau_l => aliasArray3_c(Tau_MatrixDiag(nsize+1:nsize+t0size), &
!                                                                     kmaxi,kmaxi,nSpinCant*nSpinCant)
                  p1 => Tau_MatrixDiag(nsize+1:nsize+t0size)
                  MatrixBand(j)%MatrixBlock(n)%tau_l => aliasArray3_c(p1,kmaxi,kmaxi,nSpinCant*nSpinCant)
               endif
               nsize = nsize + kmaxi*kmaxi*nSpinCant*nSpinCant
            endif
         enddo
         nk = nk + kmaxi*nSpinCant
         id_array(n) = ig
         gid_array(ig) = n
         ip_array(n) = k-1
         il_array(n) = i
         lmaxi_array(n) = lmaxi
      enddo
   enddo
!
   nsize = 0
   do j = 1, LocalNumAtoms
      kmaxj = (lmax_kkr(j)+1)**2
      do n = 1, GlobalNumAtoms
         kmaxi = MatrixBand(j)%MatrixBlock(n)%kmax_kkr
!        sc_blocks(n,j)%strcon_matrix => aliasArray2_c(SCM_MatrixBand(nsize+1:nsize+kmaxi*kmaxj), &
!                                                      kmaxi,kmaxj)
         p1 => SCM_MatrixBand(nsize+1:nsize+kmaxi*kmaxj)
         sc_blocks(n,j)%strcon_matrix => aliasArray2_c(p1,kmaxi,kmaxj)
         nsize = nsize + kmaxi*kmaxj
      enddo
   enddo
!
!  ===================================================================
   if ( max_print_level >0 .or. MyPE==0 ) then
      write(6,*) "===================================================="
      write(6,*) "init CrystalMatrix Module ::  "
      write(6,'(a,i5)') "   LocalNumAtoms   : ", LocalNumAtoms
      write(6,'(a,i5)') "   lmax            : ", lmax_kkr_max
      write(6,'(a,i5)') "   n_spin_cant     : ", nSpinCant
      write(6,'(a,i5)') "   KKR_row_dim     : ", KKRMatrixSize
      write(6,'(a,i5)') "   KKR_col_dim     : ", BandSize
      do j = 1,LocalNumAtoms
         write(6,'(a,i5)') "   Local Atom Index : ", j
         write(6,'(a,i5)') "       LocalCol_Ind : ", MatrixBand(j)%column_index
         write(6,'(a,i5)') "       kmax         : ", MatrixBand(j)%kmax_kkr
         write(6,'(a,i5)') "       lmax         : ", MatrixBand(j)%lmax_kkr
         write(6,'(a,i5)') "       global index : ", MatrixBand(j)%global_index
         if (GlobalNumAtoms < 101) then
            do i = 1, GlobalNumAtoms
               write(6,'(a,4i5)')"   block index, starting row, kmax, global Index : ", &
                                 i, MatrixBand(j)%MatrixBlock(i)%row_index,             &
                                 MatrixBand(j)%MatrixBlock(i)%kmax_kkr,                 &
                                 MatrixBand(j)%MatrixBlock(i)%global_index
            enddo
         endif
      enddo
      write(6,*) "===================================================="
   endif
!
!  ===================================================================
!  check the distribution of the matrix blocks. If lmax_kkr is the same
!  for all atoms and each process has the same number of atoms, BandSize
!  will be the same on all the processes
!  ===================================================================
   n = BandSizeCant
   call GlobalMaxInGroup(GroupID,n)
   if (n /= BandSizeCant) then
!     ----------------------------------------------------------------
      call ErrorHandler('initCrystalMatrix',                                       &
                        'ScaLAPACK will fail as BandSize is unevenly distrubuted', &
                         BandSizeCant)
!     ----------------------------------------------------------------
   endif
!
   allocate( jinv_g(tsize,GlobalNumAtoms) )
   allocate( sine_g(tsize,LocalNumAtoms) )
   allocate( tmat_g(tsize,LocalNumAtoms) )
   allocate( stcm_g(tsize,LocalNumAtoms) )
   jinv_g = CZERO ! It stores the single site Jinv matrix in the global frame for all the atoms in the unit cell
   sine_g = CZERO ! It stores the sine matrix in the global frame
   tmat_g = CZERO ! It stores the t-matrix in the global frame
   stcm_g = CZERO ! It stores the S^{T*}*C in the global frame
!
!  sould be dimension of the local matrix+blocksize
   allocate( IPVT(1:KKRMatrixSizeCant+BandSizeCant) )
!
   if (NumPEsInGroup == 1) then  ! ScaLapack will not be used for 1 process case
      return
   endif
!
#ifdef USE_SCALAPACK
!  ===================================================================
!  Initialize ScaLAPACK and set up matrix distribution
!  ===================================================================
   ICTXT = getGroupCommunicator(GroupID)
   call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEsInGroup )
!  ===================================================================
   call BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!  -------------------------------------------------------------------
   if (NPROW /= 1 .or. NPCOL /= NumPEsInGroup .or. MYROW /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('initCrystalMatrix',                                 &
              'Failed: NPROW /= 1 || NPCOL /= NumPEsInGroup || MYROW /= 0',  &
              NPROW, NPCOL, MYROW)
!     ----------------------------------------------------------------
   else if (MYCOL /= MyPEinGroup) then
!     ----------------------------------------------------------------
      call ErrorHandler('initCrystalMatrix','MYCOL /= MyPEinGroup',MYCOL,MyPEinGroup)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call syncAllPEsInGroup(GroupID)
!  -------------------------------------------------------------------
   call DESCINIT( DESC_A, KKRMatrixSizeCant, KKRMatrixSizeCant,       &
                  BandSizeCant, BandSizeCant, 0, 0, ICTXT,            &
                  KKRMatrixSizeCant, INFO )
!  -------------------------------------------------------------------
#endif
!
   end subroutine initCrystalMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCrystalMatrix()
!  ===================================================================
   use StrConstModule, only :  endStrConst
!
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   do j = 1, LocalNumAtoms
      do i = 1, GlobalNumAtoms
         nullify( MatrixBand(j)%MatrixBlock(i)%tau_l )
         nullify( MatrixBand(j)%MatrixBlock(i)%kau_l )
         nullify( sc_blocks(i,j)%strcon_matrix )
      enddo
      deallocate( MatrixBand(j)%MatrixBlock )
   enddo
   nullify(cosine_g)
   deallocate( MatrixBand, sc_blocks )
   deallocate( Kau_MatrixDiag, Tau_MatrixDiag )
   deallocate( TMP_MatrixBand )
   deallocate( KKR_MatrixBand )
   deallocate( SCM_MatrixBand )
   deallocate( id_array, jd_array, ip_array, il_array )
   deallocate( lmaxi_array, lmaxj_array, gid_array )
   deallocate( jinv_g, sine_g, tmat_g, stcm_g )
!
   if (isRelativistic) then
      deallocate( strconrel )
   endif
!
   deallocate( WORK, IWORK )
!
   isRelativistic = .false.
!
   deallocate( print_level )
   deallocate( IPVT )
!
!  -------------------------------------------------------------------
   call endStrConst()
!  -------------------------------------------------------------------
!
   if (NumPEsInGroup == 1) then  ! ScaLapack will not be used for 1 process case
      return
   endif
!
#ifdef USE_SCALAPACK
!  -------------------------------------------------------------------
   call BLACS_GRIDEXIT(ICTXT)
!  -------------------------------------------------------------------
!  call BLACS_EXIT(1)
!
#endif
!
   end subroutine endCrystalMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupSSMatrixBuf(getSingleScatteringMatrix,method,      &
                               site_config)
!  ===================================================================
   use MPPModule, only : AnyPE, MyPE
   use MPPModule, only : nbsendMessage, nbrecvMessage
   use MPPModule, only : waitMessage, setCommunicator, resetCommunicator
!
   use WriteMatrixModule,  only : writeMatrix
!
   use MatrixModule, only : computeAStar
!
   use GroupCommModule, only : getGroupCommunicator
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SpinRotationModule, only : rotateLtoG
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: method, site_config(:)
   integer (kind=IntKind) :: comm, t0size, kkri_ns, prev_pe, next_pe
   integer (kind=IntKind) :: i, n, nr, send_msgid, recv_msgid, kmax_kkr
!
   complex (kind=CmplxKind), pointer :: sm1(:,:), sm2(:,:)
   complex (kind=CmplxKind), pointer :: tm1(:,:), tm2(:,:)
   complex (kind=CmplxKind), pointer :: jm1(:,:), jm2(:,:)
   complex (kind=CmplxKind), pointer :: cm1(:,:), cm2(:,:)
   complex (kind=CmplxKind), pointer :: pm(:), gmat(:,:), p1(:)
   complex (kind=CmplxKind), pointer :: send_buf(:,:), recv_buf(:,:)
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   if (method == 2) then ! Only the t-matrix is needed
      do i = 1, LocalNumAtoms
!        =============================================================
!        Obtain the t-matrix in Global frame
!        =============================================================
         if (.not.isRelativistic .and. nSpinCant == 2) then
            tm1 => getSingleScatteringMatrix('T-Matrix',site=i,       &
                                             atom=site_config(i),spin=1,dsize=kmax_kkr)
            tm2 => getSingleScatteringMatrix('T-Matrix',site=i,       &
                                             atom=site_config(i),spin=2,dsize=kmax_kkr)
            kkri_ns = kmax_kkr*nSpinCant
            pm => tmat_g(:,i)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           ----------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, tm1, tm2, gmat)
!           ----------------------------------------------------------
         else
            tm1 => getSingleScatteringMatrix('T-Matrix',site=i,       &
                                             atom=site_config(i),dsize=kkri_ns)
         endif
         t0size = kkri_ns*kkri_ns
!        -------------------------------------------------------------
         call zcopy(t0size,tm1(1,1),1,tmat_g(1,i),1)
!        -------------------------------------------------------------
      enddo
      nullify(tm1,tm2,pm,gmat)
      return
   endif
!
!  -------------------------------------------------------------------
   comm = getGroupCommunicator(GroupID)
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
!
   send_buf => aliasArray2_c(WORK,tsize,LocalNumAtoms)
!
   if (isRelativistic) then !xianglin in between
      do i = 1, LocalNumAtoms
   !     ================================================================
   !     Obtain the Jinv-matrix, Sine-Matrix, and t-matrix in Global frame
   !     ================================================================
         kmax_kkr = MatrixBand(i)%kmax_kkr
         kkri_ns =  kmax_kkr*nSpinCant
         t0size = kmax_kkr*kmax_kkr*nSpinCant*nSpinCant
         n  = gid_array(getGlobalIndex(i)) 
            if (method == 0) then
               jm1 => getSingleScatteringMatrix('JostInv-Matrix',site=i, &
                                                atom=site_config(i))
            else
               sm1 => getSingleScatteringMatrix('Sine-Matrix',site=i,    &
                                                atom=site_config(i))
!              jm1 => aliasArray2_c(WORK(1:t0size),kkri_ns,kkri_ns)
               p1 => WORK(1:t0size)
               jm1 => aliasArray2_c(p1,kkri_ns,kkri_ns)
               jm1 = sm1
   !           ----------------------------------------------------------
   !            call computeAStar(sm1,kkri_ns,kkri_ns,jm1)                !not needed unless there is B_y
   !           ----------------------------------------------------------
            endif
            sm1 => getSingleScatteringMatrix('Sine-Matrix',site=i,atom=site_config(i))
   !        -------------------------------------------------------------
            call zcopy( t0size, jm1, 1, jinv_g(1,n), 1 )   !save jinv(Method 0) or s^t(Method 1) into jinv_g 
            call zcopy( t0size, sm1, 1, sine_g(1,i), 1 )   !save s into sine_g
   !        -------------------------------------------------------------
            tm1 => getSingleScatteringMatrix('T-Matrix',site=i,atom=site_config(i))
            call zcopy( t0size, tm1, 1, tmat_g(1,i), 1 )
   !        -------------------------------------------------------------
   !        call writeMatrix('Jinvt',jm1,kkri_ns,kkri_ns)
   !        call writeMatrix('t-mat',tm1,kkri_ns,kkri_ns)
            if (method /= 0) then
               cm1 => getSingleScatteringMatrix('Cosine-Matrix',site=i,  &
                                                atom=site_config(i))
   !           ----------------------------------------------------------
               call zgemm('t','n',kkri_ns,kkri_ns,kkri_ns,CONE,      &
                          jinv_g(1,n),kkri_ns,cm1,kkri_ns,CZERO,stcm_g(1,i),kkri_ns) !stcm_g=S^t*C
   !           ----------------------------------------------------------
            endif
   !
            if (NumPEsInGroup > 1) then
      !        =============================================================
      !        Store the jinv-matrix into the local communication buffer
      !        -------------------------------------------------------------
               call zcopy(tsize,jinv_g(1,n),1,send_buf(1,i),1)
      !        -------------------------------------------------------------
            endif
      enddo
   else
      do i = 1, LocalNumAtoms
!     ================================================================
!     Obtain the Jinv-matrix, Sine-Matrix, and t-matrix in Global frame
!     ================================================================
         kmax_kkr = MatrixBand(i)%kmax_kkr
         t0size = kmax_kkr*kmax_kkr
         n  = gid_array(getGlobalIndex(i)) ! n is the row index of the matrix block
         if ( nSpinCant == 2 ) then
            kkri_ns =  kmax_kkr*nSpinCant
!           =============================================================
!           calculate jinv_g and sine_g in global frame of reference.
!           =============================================================
            if (method == 0) then ! jinv_g contains (iS-C)^{-1}
               jm1 => getSingleScatteringMatrix('JostInv-Matrix',spin=1, &
                                                site=i,atom=site_config(i))
               jm2 => getSingleScatteringMatrix('JostInv-Matrix',spin=2, &
                                                site=i,atom=site_config(i))
            else                  ! jinv_g contains S^{*}
               sm1 => getSingleScatteringMatrix('Sine-Matrix',spin=1,    &
                                                site=i,atom=site_config(i))
               sm2 => getSingleScatteringMatrix('Sine-Matrix',spin=2,    &
                                                site=i,atom=site_config(i))
!              jm1 => aliasArray2_c(WORK(1:t0size),kmax_kkr,kmax_kkr)
               p1 => WORK(1:t0size)
               jm1 => aliasArray2_c(p1,kmax_kkr,kmax_kkr)
!              jm2 => aliasArray2_c(WORK(t0size+1:2*t0size),kmax_kkr,kmax_kkr)
               p1 => WORK(t0size+1:2*t0size)
               jm2 => aliasArray2_c(p1,kmax_kkr,kmax_kkr)
!              ----------------------------------------------------------
               call computeAStar(sm1,kmax_kkr,kmax_kkr,jm1)
               call computeAStar(sm2,kmax_kkr,kmax_kkr,jm2)
!              ----------------------------------------------------------
            endif
            pm => jinv_g(:,n)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           -------------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, jm1, jm2, gmat)
!           -------------------------------------------------------------
            sm1 => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=i,&
                                             atom=site_config(i))
            sm2 => getSingleScatteringMatrix('Sine-Matrix',spin=2,site=i,&
                                             atom=site_config(i))
            pm => sine_g(:,i)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           -------------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, sm1, sm2, gmat)
!           -------------------------------------------------------------
            tm1 => getSingleScatteringMatrix('T-Matrix',spin=1,site=i,   &
                                             atom=site_config(i))
            tm2 => getSingleScatteringMatrix('T-Matrix',spin=2,site=i,   &
                                             atom=site_config(i))
            pm => tmat_g(:,i)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           -------------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, tm1, tm2, gmat)
!           -------------------------------------------------------------
            if (method /= 0) then
               cm1 => getSingleScatteringMatrix('Cosine-Matrix',spin=1,  &
                                                site=i,atom=site_config(i))
               cm2 => getSingleScatteringMatrix('Cosine-Matrix',spin=2,  &
                                                site=i,atom=site_config(i))
               gmat => aliasArray2_c(WORK,kkri_ns,kkri_ns)
!              ----------------------------------------------------------
               call rotateLtoG(i, kmax_kkr, kmax_kkr, cm1, cm2, gmat)
               call zgemm('t','n',kkri_ns,kkri_ns,kkri_ns,CONE,jinv_g(1,n),&
                          kkri_ns,gmat,kkri_ns,CZERO,stcm_g(1,i),kkri_ns)
!              ----------------------------------------------------------
            endif
         else
            if (method == 0) then
               jm1 => getSingleScatteringMatrix('JostInv-Matrix',spin=1, &
                                                site=i,atom=site_config(i))
            else
               sm1 => getSingleScatteringMatrix('Sine-Matrix',spin=1,    &
                                                site=i,atom=site_config(i))
!              jm1 => aliasArray2_c(WORK(1:t0size),kmax_kkr,kmax_kkr)
               p1 => WORK(1:t0size)
               jm1 => aliasArray2_c(p1,kmax_kkr,kmax_kkr)
!              ----------------------------------------------------------
               call computeAStar(sm1,kmax_kkr,kmax_kkr,jm1)
!              ----------------------------------------------------------
            endif
            sm1 => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=i,&
                                             atom=site_config(i))
!           -------------------------------------------------------------
            call zcopy( t0size, jm1, 1, jinv_g(1,n), 1 )
            call zcopy( t0size, sm1, 1, sine_g(1,i), 1 )
!           -------------------------------------------------------------
            tm1 => getSingleScatteringMatrix('T-Matrix',spin=1,site=i,   &
                                             atom=site_config(i))
            call zcopy( t0size, tm1, 1, tmat_g(1,i), 1 )
!           -------------------------------------------------------------
!           call writeMatrix('Jinvt',jm1,kmax_kkr,kmax_kkr)
!           call writeMatrix('t-mat',tm1,kmax_kkr,kmax_kkr)
            if (method /= 0) then
               cm1 => getSingleScatteringMatrix('Cosine-Matrix',spin=1,  &
                                                site=i,atom=site_config(i))
!              ----------------------------------------------------------
               call zgemm('t','n',kmax_kkr,kmax_kkr,kmax_kkr,CONE,      &
                          jinv_g(1,n),kmax_kkr,cm1,kmax_kkr,CZERO,stcm_g(1,i),kmax_kkr)
!              ----------------------------------------------------------
            endif
         endif
!
         if (NumPEsInGroup > 1) then
!           =============================================================
!           Store the jinv-matrix into the local communication buffer
!           -------------------------------------------------------------
            call zcopy(tsize,jinv_g(1,n),1,send_buf(1,i),1)
!           -------------------------------------------------------------
         endif
      enddo
!
      nullify(sm1,sm2,jm1,jm2,tm1,tm2,pm,gmat)
!
      recv_buf => aliasArray2_c(Tau_MatrixDiag,tsize,LocalNumAtoms) ! Use Tau_MatrixBand as temporary receiving buffer space
!
      nr = 1
      do while (nr <=  NumPEsInGroup-1)
         next_pe = mod(MyPEInGroup+nr,NumPEsInGroup)
         prev_pe = mod(MyPEInGroup-nr+NumPEsInGroup,NumPEsInGroup)
!        ----------------------------------------------------------------
         recv_msgid = nbrecvMessage(recv_buf,tsize,LocalNumAtoms,1010,prev_pe)
!        ----------------------------------------------------------------
         send_msgid = nbsendMessage(send_buf,tsize,LocalNumAtoms,1010,next_pe)
!        ----------------------------------------------------------------
         call waitMessage(recv_msgid)
!        ----------------------------------------------------------------
         do i = 1, LocalNumAtoms
            n  = gid_array(getGlobalIndex(i,prev_pe)) ! n is the row index of the matrix block
!           -------------------------------------------------------------
            call zcopy(tsize,recv_buf(1,i),1,jinv_g(1,n),1)
!           -------------------------------------------------------------
         enddo
!        ----------------------------------------------------------------
         call waitMessage(send_msgid)
!        ----------------------------------------------------------------
         nr = nr + 1
      enddo
   endif
!
!  -------------------------------------------------------------------
   call resetCommunicator(sync=.true.)
!  -------------------------------------------------------------------
   nullify(recv_buf, send_buf)
!
   end subroutine setupSSMatrixBuf
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calCrystalMatrix_sumk(e,getSingleScatteringMatrix,      &
                                    tau_needed,use_tmat,configuration)
!  ===================================================================
   use MPPModule, only : MyPE
   use MatrixModule, only : computeAStarT
   use BZoneModule, only : getNumKs, getAllWeights, getAllKPoints,    &
                           getWeightSum
   use IBZRotationModule, only : getNumIBZRotations
   use ProcMappingModule, only : isKPointOnMyProc, getNumKsOnMyProc,  &
                                 getKPointIndex, getNumRedundantKsOnMyProc
   use GroupCommModule, only : getGroupID, GlobalSumInGroup, getMyPEinGroup
   use StrConstModule, only : getStrConstMatrix
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   logical, intent(in), optional :: tau_needed
   logical, intent(in), optional :: use_tmat
   logical :: calculate_tau = .false.
!
   character (len=20) :: sname = "calCrystalMatrix"
!
   integer (kind=IntKind), intent(in), optional :: configuration(:)
   integer (kind=IntKind) :: k_loc, k, i, row, col, MyPEinKGroup, method
   integer (kind=IntKind) :: NumKs, kGID, aGID, NumKsOnMyProc, NumRedunKs
   integer (kind=IntKind) :: site_config(LocalNumAtoms)
!
   real (kind=RealKind), pointer :: kpts(:,:), weight(:)
   real (kind=RealKind) :: kfac, kaij
   real (kind=RealKind) :: weightSum, kvec(1:3), aij(3)
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind) :: wfac, efac
   complex (kind=CmplxKind), pointer :: rotmat(:,:), scm(:,:)
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   energy = e
   if (isRelativistic) then !xianglin
      kappa = sqrt(2.d0*Me*e + e**2/LightSpeed**2)
   else
      kappa = sqrt(e)
   endif
!
!  ===================================================================
!  There are different ways to construct the KKR matrix or its modified
!  form. A (modified) KKR Matrix (or the M-matrix) can be constructed 
!  as follows.
!  Method 0:
!    p_MatrixBand = 
!      [1 - OmegaHat(e) * S(e)^{T*} * (B(e,k)+i*kappa) * S(e)/kappa]^{-1} - 1
!  where, OmegaHat(e) * S(e)^{T*} = Jost(e)^{-1}.
!
!  Method 1:
!    p_MatrixBand =
!      -[S(e)^{T*} * C(e) + S(e)^{T*} * B(k,e) * S(e)/kappa]^{-1}
!
!  Method 2:
!    p_MatrixBand = [ 1 - (B(k,e)+i*kappa) * t(e) ]^{-1}
!
!  By default, set method = 0 so that both sine and Jost matrices are
!  involved.
!  Using method = 1, the numerical noise seems breaking the symmetry
!  for lmax > 3 and gives rise to a non-zero force, even though
!  very small, acting on each atom in a crystal.
!  Using method = 2, calculting tau-matrix will only involve t-matrix.
!
!  KKR_MatrixBand is the integration of p_MatrixBand over IBZ
!  ===================================================================
   method = 0
   if (present(use_tmat)) then
      if (use_tmat) then
         method = 2 ! When constructing tau-matrix, only t-matrix used
      endif
   endif
!
   if (present(configuration)) then
      site_config(1:LocalNumAtoms) = configuration(1:LocalNumAtoms)
   else
      site_config = 1
   endif
!
!  ===================================================================
!  Exchange single site scattering matrix among processes.
!  -------------------------------------------------------------------
   call setupSSMatrixBuf(getSingleScatteringMatrix,method,site_config)
!  -------------------------------------------------------------------
!
   if (method == 2) then
      calculate_tau = .true.
   else if ( present(tau_needed) ) then
      calculate_tau = tau_needed
   else 
      calculate_tau = .false.
   endif
!
   kGID = getGroupID('K-Mesh')
   aGID = getGroupID('Unit Cell')
   NumKsOnMyProc = getNumKsOnMyProc()
   NumRedunKs = getNumRedundantKsOnMyProc()
   MyPEinKGroup = getMyPEinGroup(kGID)
!
   NumKs = getNumKs()
   kpts => getAllKPoints(kfac)
   weight => getAllWeights()
   weightSum = getWeightSum()
!
!  ===================================================================
!  Loop over k-points on mesh
!  ===================================================================
   KKR_MatrixBand = CZERO
   do k_loc = 1,NumKsOnMyProc
!     ================================================================
!     Normorlize BZ integration weights
!     ================================================================
      k = getKPointIndex(k_loc)
      kvec(1:3) = kpts(1:3,k)*kfac
!     ================================================================
!     get structure constant matrix for the k-point and energy
      do col = 1, LocalNumAtoms
         do row = 1, GlobalNumAtoms
!           ----------------------------------------------------------
            scm => getStrConstMatrix(kvec,kappa,id_array(row),jd_array(col), &
                                     lmaxi_array(row),lmaxj_array(col),aij)
!           ----------------------------------------------------------
!           kaij = kvec(1)*aij(1)+kvec(2)*aij(2)+kvec(3)*aij(3)
!           efac = exp(sqrtm1*kaij)
!           sc_blocks(row,col)%strcon_matrix = scm*efac
            sc_blocks(row,col)%strcon_matrix = scm
         enddo
      enddo
!
      wfac = weight(k)/weightSum
!     write(6,'(a,i4,a,3f12.5,a,2d12.5,a,2d12.5)')'k-ind = ',k,       &
!             ', kvec = ',kvec,', wfac = ',wfac,', weightsum = ',weightSum
!
!     ================================================================
!     Compute the modified KKR matrix, which is stored in TMP_MatrixBand
!     ----------------------------------------------------------------
      call computeMatrixBand(TMP_MatrixBand,method)
!     ----------------------------------------------------------------
!     call checkMatrixBandRotation(kvec,TMP_MatrixBand)
!     ----------------------------------------------------------------
!
!     ================================================================
!     Sum over the matrix over the processors that take care different
!     k-points. The loop will continue for the redundant k-points, if 
!     there are any. The redundant k-points are needed for the load 
!     balance purpose.
!     ================================================================
      if (k_loc <= NumKsOnMyProc - NumRedunKs .or. MyPEinKGroup == 0) then
!        -------------------------------------------------------------
         call zaxpy(KKRMatrixSizeCant*BandSizeCant,wfac,TMP_MatrixBand,1, &
                    KKR_MatrixBand,1)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call zaxpy(KKRMatrixSizeCant*BandSizeCant,CZERO,TMP_MatrixBand,1, &
                    KKR_MatrixBand,1)
!        -------------------------------------------------------------
      endif
!     call writeMatrix('KKR_MatrixBand',KKR_MatrixBand,KKRMatrixSizeCant, &
!                      BandSizeCant,TEN2m8,.true.)
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(kGID,KKR_MatrixBand,KKRMatrixSizeCant*BandSizeCant)
!  -------------------------------------------------------------------
!
   if (method == 0 .or. method == 1) then
!     ----------------------------------------------------------------
      call computeKauMatrix(getSingleScatteringMatrix,method,calculate_tau,site_config)
!     ----------------------------------------------------------------
   else if (method == 2) then
!     ----------------------------------------------------------------
      call computeTauMatrix(getSingleScatteringMatrix,site_config)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('calCrystalMatrix_sumk','The method is invalid',method)
!     ----------------------------------------------------------------
   endif
!
   if (getNumIBZRotations() > 1) then
      if (isRelativistic) then
         call ErrorHandler('calCrystalMatrix_sumk','Relativistic sumIBZRotation not implemented',INFO)
      else
!        ----------------------------------------------------------------
         call sumIBZRotation(calculate_tau)
!        ----------------------------------------------------------------
      endif
   endif
!
   nullify( weight, kpts )
!
   if (trim(stop_routine) ==trim(sname)) then
      stop 'At calCrystalMatrix'
   endif
!
   end subroutine calCrystalMatrix_sumk
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calCrystalMatrix_k(e,kpts,getSingleScatteringMatrix,    &
                                 tau_needed,use_tmat,configuration)
!  ===================================================================
   use StrConstModule, only : getStrConstMatrix
!
   implicit none
!
   logical, intent(in), optional :: tau_needed
   logical, intent(in), optional :: use_tmat
   logical :: calculate_tau
!
   character (len=20) :: sname = "calCrystalMatrix"
!
   integer (kind=IntKind), intent(in), optional :: configuration(:)
   integer (kind=IntKind) :: row, col, method
   integer (kind=IntKind) :: site_config(LocalNumAtoms)
!
   real (kind=RealKind), intent(in) :: kpts(1:3)
!
   complex (kind=CmplxKind), intent(in) :: e
!
   complex (kind=CmplxKind), pointer :: scm(:,:)
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   energy = e
   if (isRelativistic) then !xianglin
      kappa = sqrt(2.d0*Me*e + e**2/LightSpeed**2)
   else
      kappa = sqrt(e)
   endif
!
   if ( present(tau_needed) ) then
      calculate_tau = tau_needed
   else 
      calculate_tau = .false.
   endif
!
!  ===================================================================
!  There are different ways to construct the KKR matrix or its modified
!  form. A (modified) KKR Matrix (or the M-matrix) can be constructed 
!  as follows.
!  Method 0:
!    p_MatrixBand = 
!      [1 - OmegaHat(e) * S(e)^{T*} * (B(e,k)+i*kappa) * S(e)/kappa]^{-1} - 1
!  where, OmegaHat(e) * S(e)^{T*} = Jost(e)^{-1}.
!
!  Method 1:
!    p_MatrixBand =
!      -[S(e)^{T*} * C(e) + S(e)^{T*} * B(k,e) * S(e)/kappa]^{-1}
!
!  Method 2:
!    p_MatrixBand = [ 1 - (B(k,e)+i*kappa) * t(e) ]^{-1}
!
!  By default, set method = 0 so that both sine and Jost matrices are
!  involved.
!  Using method = 1, the numerical noise seems breaking the symmetry
!  for lmax > 3 and gives rise to a non-zero force, even though
!  very small, acting on each atom in a crystal.
!  Using method = 2, calculting tau-matrix will only involve t-matrix.
!  ===================================================================
   method = 0
   if (present(use_tmat)) then
      if (use_tmat) then
         method = 2
      endif
   endif
!
   if (present(configuration)) then
      site_config(1:LocalNumAtoms) = configuration(1:LocalNumAtoms)
   else
      site_config = 1
   endif
!
!  ===================================================================
!  Exchange single site scattering matrix among processes.
!  -------------------------------------------------------------------
   call setupSSMatrixBuf(getSingleScatteringMatrix,method,site_config)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  get structure constant matrix for the k-point and energy
!  ===================================================================
   do col = 1, LocalNumAtoms
      do row = 1, GlobalNumAtoms
!        -------------------------------------------------------------
         scm => getStrConstMatrix(kpts,kappa,id_array(row),jd_array(col), &
                                  lmaxi_array(row),lmaxj_array(col))
!        -------------------------------------------------------------
         sc_blocks(row,col)%strcon_matrix = scm
      enddo
   enddo
!
!  ===================================================================
!  Compute the modified KKR matrix, which is stored in KKR_MatrixBand
!  -------------------------------------------------------------------
   call computeMatrixBand(KKR_MatrixBand,method)
!  -------------------------------------------------------------------
   if (method == 0. .or. method == 1) then
!     ----------------------------------------------------------------
      call computeKauMatrix(getSingleScatteringMatrix,method,         &
                         calculate_tau,site_config)
!     ----------------------------------------------------------------
   else if (method == 2) then
!     ----------------------------------------------------------------
      call computeTauMatrix(getSingleScatteringMatrix,site_config)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('calCrystalMatrix_k','The method is invalid',method)
!     ----------------------------------------------------------------
   endif
!
   end subroutine calCrystalMatrix_k
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeKauMatrix(getSingleScatteringMatrix,method,      &
                               calculate_tau,site_config)
!  ===================================================================
   use MatrixModule, only : computeUAUts
!
   use WriteMatrixModule,  only : writeMatrix
!
   use SpinRotationModule, only : rotateGtoL
!
   use StrConstModule, only : getStrConstMatrix
!
   implicit none
!
   logical, intent(in) :: calculate_tau
!
   integer (kind=IntKind), intent(in) :: method, site_config(:)
   integer (kind=IntKind) :: i, j, ig, np, ni, nj, kl, klp, is, js, ns
   integer (kind=IntKind) :: kmaxj, kmaxj_ns, t0size
!
   real (kind=RealKind) :: kvec(3)
!
   complex (kind=CmplxKind), pointer :: kau_l(:,:), tau_l(:,:)
   complex (kind=CmplxKind), pointer :: sm1(:,:), sm2(:,:), pw(:,:), poi(:,:)
   complex (kind=CmplxKind), pointer :: tm(:,:), om(:,:), oim(:,:)
   complex (kind=CmplxKind), pointer :: wau_g(:,:), wau_l(:,:,:), p1(:)
   complex (kind=CmplxKind),allocatable :: mat_tmp(:,:) , mat_tmp2(:,:)!xianglin
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
!  ===================================================================
!  calculate Kau, the modified KKR Matrix as follows:
!    Kau = e * [S(e)]^{-1} * [tau(k,e)-t(e)] * [S(e)]^{-T*} 
!        = ( [ 1 - (iS(e)-C(e))^{-1} * (B(k,e)+i*kappa) * S(e)/kappa ]^{-1}
!           - 1) * OmegaHat(e) * kappa
!    Note:
!
!    Method = 0 =>  KKR_MatrixBand = 
!        [ 1 - (iS(e)-C(e))^{-1} * (B(k,e)+i*kappa) * S(e)/kappa ]^{-1} - 1
!
!    Method = 1 =>  KKR_MatrixBand = 
!       - [ S(e)^{T*} * C(e) + S(e)^{T*} * B(k,e) * S(e)/kappa ]^{-1}
!  ===================================================================
   do j = 1, LocalNumAtoms
      ig = MatrixBand(j)%global_index
      i = gid_array(ig)  ! Here, we only calculate wau_g for the diagonal components of the Kau matrices
      kmaxj = MatrixBand(j)%MatrixBlock(i)%kmax_kkr
      kmaxj_ns = kmaxj*nSpinCant
      wau_g => aliasArray2_c(WORK,kmaxj_ns,kmaxj_ns)
      nj = MatrixBand(j)%column_index-1
      ni = MatrixBand(j)%MatrixBlock(i)%row_index-1
      do kl = 1, kmaxj_ns
         np = KKRMatrixSizeCant*(nj+kl-1)+ni
         do klp = 1, kmaxj_ns
            wau_g(klp,kl) = KKR_MatrixBand(np+klp)
         enddo
      enddo
      if (isRelativistic) then !xianglin
         wau_l => aliasArray3_c(WORK,kmaxj_ns,kmaxj_ns,1)
      else
         if ( nSpinCant == 2 ) then
            t0size = kmaxj_ns*kmaxj_ns
            p1 => WORK(t0size+1:t0size*2)
            wau_l => aliasArray3_c(p1,kmaxj,kmaxj,nSpinCant*nSpinCant)
   !        -------------------------------------------------------------
            call rotateGtoL(j, kmaxj, kmaxj, wau_g, wau_l)
   !        -------------------------------------------------------------
         else
            wau_l => aliasArray3_c(WORK,kmaxj,kmaxj,1) ! wau_l = wau_g
         endif
      endif
!
      if (isRelativistic) then !xianglin
         om => getSingleScatteringMatrix('OmegaHat-Matrix',site=j,       &
                                         atom=site_config(j))
         kau_l => MatrixBand(j)%MatrixBlock(i)%kau_l(:,:,1)
         pw => wau_l(:,:,1)
!   open(131,file="pw",action="write")
!   write(131,*) pw
!   close(131)

!   open(131,file="om",action="write")
!   write(131,*) om
!   close(131)
         if (method == 0) then
   !              -------------------------------------------------------
                  call zgemm('n', 'n', kmaxj_ns, kmaxj_ns, kmaxj_ns, kappa,       &
                             pw, kmaxj_ns, om, kmaxj_ns, CZERO, kau_l, kmaxj_ns)
   !              -------------------------------------------------------
         else
            kau_l = kappa*(pw-om)
         endif 
      else
         ns = 0
         do js = 1, nSpinCant
            om => getSingleScatteringMatrix('OmegaHat-Matrix',spin=js,&
                                            site=j,atom=site_config(j))
   !        call writeMatrix('Ohat-mat',om,kmaxj,kmaxj)
            do is = 1, nSpinCant
               ns = ns + 1
               kau_l => MatrixBand(j)%MatrixBlock(i)%kau_l(:,:,ns)
               pw => wau_l(:,:,ns)
               if (method == 0) then
   !              -------------------------------------------------------
                  call zgemm('n', 'n', kmaxj, kmaxj, kmaxj, kappa,       &
                             pw, kmaxj, om, kmaxj, CZERO, kau_l, kmaxj)
   !              -------------------------------------------------------
   !              call writeMatrix('kau_l',kau_l,kmaxj,kmaxj)
               else
                  oim => getSingleScatteringMatrix('OmegaHatInv-Matrix', &
                                                   spin=js,site=j,atom=site_config(j))
                  if (is == js) then
                     kau_l = kappa*(pw-om)  ! Warning: For large l component, 
                                            !          both pw and om are big numbers
   !ywg              poi => MatrixBand(j)%MatrixBlock(i)%tau_l(:,:,ns) ! used as temporary space
   !                 ----------------------------------------------------
   !ywg              call zgemm('n', 'n', kmaxj, kmaxj, kmaxj, CONE,     &
   !ywg                         pw, kmaxj, oim, kmaxj, CZERO, poi, kmaxj)
   !                 ----------------------------------------------------
   !ywg              do kl = 1, kmaxj
   !ywg                 poi(kl,kl) = poi(kl,kl)-CONE
   !ywg              enddo
   !                 ----------------------------------------------------
   !ywg              call zgemm('n', 'n', kmaxj, kmaxj, kmaxj, kappa,    &
   !ywg                         poi, kmaxj, om, kmaxj, CZERO, kau_l, kmaxj)
   !                 ----------------------------------------------------
                  else
                     kau_l = kappa*pw
                  endif
               endif
            enddo
         enddo
      endif
!
      if (calculate_tau) then
         if (isRelativistic) then
!           by xianglin, need to be revised to use s_d when include B_y
            sm2 => getSingleScatteringMatrix('Sine-Matrix',site=j,atom=site_config(j))
            sm1 => getSingleScatteringMatrix('Sine-Matrix',site=j,atom=site_config(j))
            tm => getSingleScatteringMatrix('T-Matrix',site=j,atom=site_config(j))
            kau_l => MatrixBand(j)%MatrixBlock(i)%kau_l(:,:,1)
            tau_l => MatrixBand(j)%MatrixBlock(i)%tau_l(:,:,1)
            allocate(mat_tmp(kmaxj_ns,kmaxj_ns))
            allocate(mat_tmp2(kmaxj_ns,kmaxj_ns))
            mat_tmp2=kau_l !note in the following CONE/kappa**2 used instead of divide by energy
            call zgemm( 'n', 't', kmaxj_ns, kmaxj_ns, kmaxj_ns, CONE/kappa**2,&
            mat_tmp2, kmaxj_ns, sm2, kmaxj_ns, CZERO, mat_tmp, kmaxj_ns)
            call zgemm( 'n', 'n', kmaxj_ns, kmaxj_ns, kmaxj_ns, CONE,&
            sm1, kmaxj_ns, mat_tmp, kmaxj_ns, CZERO, tau_l, kmaxj_ns)
            tau_l = tau_l + tm
            deallocate(mat_tmp, mat_tmp2)
         else
            ns = 0
            do js = 1, nSpinCant
               sm2 => getSingleScatteringMatrix('Sine-Matrix',spin=js,&
                                                site=j,atom=site_config(j))
               tm => getSingleScatteringMatrix('T-Matrix',spin=js,    &
                                                site=j,atom=site_config(j))
               do is = 1, nSpinCant
                  sm1 => getSingleScatteringMatrix('Sine-Matrix',     &
                                                   spin=is,site=j,atom=site_config(j))
                  ns = ns + 1
                  kau_l => MatrixBand(j)%MatrixBlock(i)%kau_l(:,:,ns)
                  tau_l => MatrixBand(j)%MatrixBlock(i)%tau_l(:,:,ns)
   !              -------------------------------------------------------
                  call computeUAUts(sm1,kmaxj,kmaxj,sm2,kmaxj,CONE/energy,&
                                    kau_l,kmaxj,CZERO,tau_l,kmaxj,WORK)
   !              -------------------------------------------------------
                  if (is == js) then
                     tau_l = tau_l + tm
                  endif
               enddo
            enddo
         endif
      endif
   enddo
!  *********************
!  test for L = 0 case *
!  *********************
!  sm2 => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=1,atom=site_config(1))
!  tm => getSingleScatteringMatrix('T-Matrix',spin=1,site=1,atom=site_config(1))
!  kvec  = ZERO
!  sm1 => getStrConstMatrix(kvec,kappa,1,1,0,0)
!  kau_l => MatrixBand(1)%MatrixBlock(1)%kau_l(:,:,1)
!  kau_l(1,1) = energy * (CONE/(CONE/tm(1,1) - sm1(1,1) - sqrtm1*kappa)-tm(1,1))/sm2(1,1)**2
!  kau_l(1,1) = energy * CONE/(CONE/tm(1,1) - sm1(1,1) - sqrtm1*kappa)/sm2(1,1)**2
!  =====================
!
   nullify( wau_g, wau_l, kau_l, tau_l, sm1, sm2, tm, om )
!
   end subroutine computeKauMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeTauMatrix(getSingleScatteringMatrix,site_config)
!  ===================================================================
   use WriteMatrixModule,  only : writeMatrix
!
   use SpinRotationModule, only : rotateGtoL
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: site_config(:)
   integer (kind=IntKind) :: i, j, ig, np, ni, nj, kl, klp, is, js, ns
   integer (kind=IntKind) :: kmaxj, kmaxj_ns, t0size
!
   complex (kind=CmplxKind), pointer :: tau_l(:,:), pw(:,:), tm(:,:)
   complex (kind=CmplxKind), pointer :: wau_g(:,:), wau_l(:,:,:), p1(:)
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   do j = 1, LocalNumAtoms
      ig = MatrixBand(j)%global_index
      i = gid_array(ig)  ! Here, we only calculate wau_g for the diagonal components of the Tau matrices
      kmaxj = MatrixBand(j)%MatrixBlock(i)%kmax_kkr
      kmaxj_ns = kmaxj*nSpinCant
      wau_g => aliasArray2_c(WORK,kmaxj_ns,kmaxj_ns)
      nj = MatrixBand(j)%column_index-1
      ni = MatrixBand(j)%MatrixBlock(i)%row_index-1
      do kl = 1, kmaxj_ns
         np = KKRMatrixSizeCant*(nj+kl-1)+ni
         do klp = 1, kmaxj_ns
            wau_g(klp,kl) = KKR_MatrixBand(np+klp)
         enddo
      enddo
!
!     ================================================================
!     Determine wau_l which is wau_g in local frame.
!     ================================================================
      if (isRelativistic) then !xianglin
         wau_l => aliasArray3_c(WORK,kmaxj_ns,kmaxj_ns,1)
      else if ( nSpinCant == 2 ) then
         t0size = kmaxj_ns*kmaxj_ns
         p1 => WORK(t0size+1:t0size*2)
         wau_l => aliasArray3_c(p1,kmaxj,kmaxj,nSpinCant*nSpinCant)
!        -------------------------------------------------------------
         call rotateGtoL(j, kmaxj, kmaxj, wau_g, wau_l)
!        -------------------------------------------------------------
      else
         wau_l => aliasArray3_c(WORK,kmaxj,kmaxj,1) ! wau_l = wau_g
      endif
!
      if (isRelativistic) then !xianglin
         tm => getSingleScatteringMatrix('T-Matrix',site=j,atom=site_config(j))
         tau_l => MatrixBand(j)%MatrixBlock(i)%tau_l(:,:,1)
         pw => wau_l(:,:,1)
!        -------------------------------------------------------------
         call zgemm('n', 'n', kmaxj_ns, kmaxj_ns, kmaxj_ns, CONE,     &
                    tm, kmaxj_ns, pw, kmaxj_ns, CZERO, tau_l, kmaxj_ns)
!        -------------------------------------------------------------
      else
         ns = 0
         do js = 1, nSpinCant
            tm => getSingleScatteringMatrix('T-Matrix',spin=js,       &
                                            site=j,atom=site_config(j))
!           call writeMatrix('tm',tm,kmaxj,kmaxj,TEN2m8)
!           ==========================================================
!           Note: wau_l is the integration of [ 1 - (B(k,e)+i*kappa) * t(e) ]^{-1}
!                 over IBZ.
!           ==========================================================
            do is = 1, nSpinCant
               ns = ns + 1
               tau_l => MatrixBand(j)%MatrixBlock(i)%tau_l(:,:,ns)
               pw => wau_l(:,:,ns)
!              -------------------------------------------------------
               call zgemm('n', 'n', kmaxj, kmaxj, kmaxj, CONE,        &
                          tm, kmaxj, pw, kmaxj, CZERO, tau_l, kmaxj)
!              -------------------------------------------------------
!              call writeMatrix('tau',tau_l,kmaxj,kmaxj,TEN2m8)
            enddo
         enddo
      endif
   enddo
!
   end subroutine computeTauMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeMatrixBand(p_MatrixBand,method)
!  ===================================================================
   use MatrixModule, only : setupUnitMatrix
!
   use MatrixInverseModule, only : MtxInv_LU
!
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: method
   integer (kind=IntKind) :: kmaxj, kmaxj_ns, kmaxi, kmaxi_ns, t0size
   integer (kind=IntKind) :: j, nj, ni, i, is, ig, nc, kl, klp, n, nr
!
   complex (kind=CmplxKind), pointer :: strcon(:,:)
   complex (kind=CmplxKind), pointer :: jinvB(:)
   complex (kind=CmplxKind), pointer :: p_jinvi(:)
   complex (kind=CmplxKind), pointer :: p_sinej(:), w2(:,:)
   complex (kind=CmplxKind), intent(out), target :: p_MatrixBand(:)
   complex (kind=CmplxKInd) :: cfac, wfac
!
   interface
      subroutine convertGijToRel(gij, bgij, kkr1, kkr2, ce)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: kkr1, kkr2
         complex (kind=CmplxKind), intent(in) :: gij(:,:)
         complex (kind=CmplxKind), intent(out) :: bgij(:,:)
         complex (kind=CmplxKind), intent(in) :: ce
      end subroutine convertGijToRel
   end interface
!
!  ===================================================================
!  calculate the following modified KKR Matrix (or the M-matrix).
!  Method 0:
!    p_MatrixBand = 
!      [1 - OmegaHat(e) * S(e)^{T*} * (B(e,k)+i*kappa) * S(e)/kappa]^{-1} - 1
!  where, OmegaHat(e) * S(e)^{T*} = Jost(e)^{-1}.
!
!  Method 1:
!    p_MatrixBand =
!      -[S(e)^{T*} * C(e) + S(e)^{T*} * B(k,e) * S(e)/kappa]^{-1}
!
!  Method 2:
!    p_MatrixBand = [ 1 - (B(k,e)+i*kappa) * t(e) ]^{-1}
!  ===================================================================
   cfac = SQRTm1*kappa
   do j = 1, LocalNumAtoms
      kmaxj = MatrixBand(j)%kmax_kkr
      kmaxj_ns = kmaxj*nSpinCant
      nj = MatrixBand(j)%column_index-1
      ig = MatrixBand(j)%global_index ! "ig" is the global index of the corresponding atom
      nc = gid_array(ig)              ! "nc" is the column index of the block in the big matrix
      do i = 1, GlobalNumAtoms        ! "i" is the row index of the matrix block
         kmaxi = MatrixBand(j)%MatrixBlock(i)%kmax_kkr
         kmaxi_ns = kmaxi*nSpinCant
         t0size = kmaxi_ns*kmaxi_ns
!        =============================================================
!        Method = 0: Obtain p_jinvi = OmegaHat(e) * S(e)^{T*}, 
!        which is the inverse of the Jost matrix: [iS(e) - C(e)],
!                      i.e., p_jinvi = [iS(e) - C(e)]^{-1}
!
!        Method = 1: Obtain p_jinvi = S(e)^{T*}
!        =============================================================
         jinvB => p_MatrixBand(1:kmaxi_ns*kmaxj_ns)  ! Use p_MatrixBand as temporary space...
         p_jinvi => jinv_g(1:kmaxi_ns*kmaxj_ns,i)
!
!        =============================================================
!        Obtain,
!           method 0: jinvB = i*kappa*OmegaHat(e) * S(e)^{T*}
!                                    + OmegaHat(e) * S(e)^{T*} * B(e,k)
!           method 1: jinvB = S(e)^{T*} * B(e,k)
!           method 2: jinvB = i*kappa + B(e,k)
!        =============================================================
         strcon => sc_blocks(i,j)%strcon_matrix(:,:)
!        call writeMatrix('strcon',strcon,kmaxi,kmaxj)
         if (isRelativistic) then
!           ----------------------------------------------------------
            call convertGijToRel(strcon, strconrel, kmaxi, kmaxj, energy)
!           ----------------------------------------------------------
            strcon => strconrel
            if (method == 0) then
               if (i == nc) then
                  jinvB = cfac*p_jinvi
!                 ----------------------------------------------------
                  call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxi_ns, CONE, &
                             p_jinvi, kmaxi_ns, strcon, kmaxi_ns, CONE, jinvB, kmaxi_ns)
!                 ----------------------------------------------------
               else
!                 ----------------------------------------------------
                  call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxi_ns, CONE, &
                             p_jinvi, kmaxi_ns, strcon, kmaxi_ns, CZERO, jinvB, kmaxi_ns)
!                 ----------------------------------------------------
               endif
            else if (method == 1) then
!              -------------------------------------------------------
               call zgemm('t', 'n', kmaxi_ns, kmaxj_ns, kmaxi_ns, CONE, &
                          p_jinvi, kmaxi_ns, strcon, kmaxi_ns, CZERO, jinvB, kmaxi_ns)
!              -------------------------------------------------------
            else if (method == 2) then
               w2 => aliasArray2_c(jinvB,kmaxi_ns,kmaxi_ns)
               if (i == nc) then
!                 ----------------------------------------------------
                  call setupUnitMatrix(kmaxi_ns,w2,cfac)
!                 ----------------------------------------------------
                  w2 = w2 + strcon
               else
                  w2 = strcon
               endif
            endif
         else
            if (method == 0) then
               if (i == nc) then
                  jinvB = cfac*p_jinvi
                  do is = 1, nSpinCant
!                    -------------------------------------------------
                     call zgemm('n', 'n', kmaxi_ns, kmaxj, kmaxi, CONE,    &
                                p_jinvi((is-1)*kmaxi_ns*kmaxi+1), kmaxi_ns,&
                                strcon, kmaxi, CONE, jinvB((is-1)*kmaxi_ns*kmaxj+1), kmaxi_ns)
!                    -------------------------------------------------
                  enddo
               else
                  do is = 1, nSpinCant
!                    -------------------------------------------------
                     call zgemm('n', 'n', kmaxi_ns, kmaxj, kmaxi, CONE,    &
                                p_jinvi((is-1)*kmaxi_ns*kmaxi+1), kmaxi_ns,&
                                strcon, kmaxi, CZERO, jinvB((is-1)*kmaxi_ns*kmaxj+1), kmaxi_ns)
!                    -------------------------------------------------
                  enddo
               endif
            else if (method == 1) then
               do is = 1, nSpinCant
!                 ----------------------------------------------------
                  call zgemm('t', 'n', kmaxi_ns, kmaxj, kmaxi, CONE,    &
                             p_jinvi((is-1)*kmaxi_ns*kmaxi+1), kmaxi_ns,&
                             strcon, kmaxi, CZERO, jinvB((is-1)*kmaxi_ns*kmaxj+1), kmaxi_ns)
!                 ----------------------------------------------------
               enddo
            else if (method == 2) then
               w2 => aliasArray2_c(jinvB,kmaxi_ns,kmaxi_ns)
               if (i == nc) then
!                 ----------------------------------------------------
                  call setupUnitMatrix(kmaxi_ns,w2,cfac)
!                 ----------------------------------------------------
                  w2 = w2 + strcon
               else
                  w2 = strcon
               endif
            endif
         endif
!
!        =============================================================
!        Obtain,
!          method 0: WORK = OmegaHat(e) * S(e)^{T*}
!                                  * (B(e,k)+i*kappa) * S(e)/kappa
!          method 1: WORK = S(e)^{T*} * B(e,k) * S(e)/kappa
!          method 2: WORK = [B(e,k)+i*kappa] * t(e)
!        =============================================================
         t0size = kmaxj_ns*kmaxj_ns
         if (method == 0 .or. method == 1) then
            p_sinej => sine_g(:,j)
            wfac = CONE/kappa
         else
            p_sinej => tmat_g(:,j)
            wfac = CONE
         endif
         ni = MatrixBand(j)%MatrixBlock(i)%row_index-1
!        -------------------------------------------------------------
         call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxj_ns, wfac,     &
                    jinvB, kmaxi_ns, p_sinej, kmaxj_ns, CZERO,        &
                    WORK(KKRMatrixSizeCant*nj+ni+1), KKRMatrixSizeCant)
!        -------------------------------------------------------------
      enddo
   enddo
!  call writeMatrix('WORK',WORK,KKRMatrixSizeCant,KKRMatrixSizeCant,TEN2m8,.true.)
!
   nullify(jinvB)
!
   if (method == 0 .or. method == 2) then
!     ================================================================
!     For method = 0: 
!        WORK = (iS(e)-C(e))^{-1} * (B(k,e)+i*kappa) * S(e)/kappa
!
!        p_MatrixBand = [1 - WORK]^{-1} - 1
!        = [ 1 - (iS(e)-C(e))^{-1} * (B(k,e)+i*kappa) * S(e)/kappa ]^{-1} - 1
!
!        Will call ZGETRF and ZGETRI to calculate WORK^{-1} and then to
!        solve p_MatrixBand.
!
!     For method = 2:
!        WORK = (B(k,e)+i*kappa) * t(e)
!
!        p_MatrixBand = [1 - WORK]^{-1} = [ 1 - (B(k,e)+i*kappa) * t(e) ]^{-1}
!
!        Will call ZGETRF and ZGETRI to solve WORK * p_MatrixBand = 1
!     ================================================================
!     -WORK => p_MatrixBand
!     p_MatrixBand = -WORK
      n = size(p_MatrixBand)
      call zcopy(n,WORK,1,p_MatrixBand,1)
      p_MatrixBand = -p_MatrixBand
!     ----------------------------------------------------------------
      do j = 1, LocalNumAtoms
         nc = MatrixBand(j)%column_index
         LOOP_n2: do n = 1, GlobalNumAtoms
            if ( MatrixBand(j)%global_index ==                        &
                 MatrixBand(j)%MatrixBlock(n)%global_index ) then
               nr = MatrixBand(j)%MatrixBlock(n)%row_index
               do i = 1, MatrixBand(j)%kmax_kkr*nSpinCant
                  p_MatrixBand(i+(nc-1+i-1)*KKRMatrixSizeCant+nr-1) = &
                     CONE+p_MatrixBand(i+(nc-1+i-1)*KKRMatrixSizeCant+nr-1)
               enddo
               exit LOOP_n2
            endif
         enddo LOOP_n2
      enddo
!     call writeMatrix('p_MatrixBand',p_MatrixBand,KKRMatrixSizeCant, &
!                      KKRMatrixSizeCant,TEN2m8)
   else if (method == 1) then
!     ================================================================
!     WORK = S(e)^{T*} * B(k,e) * S(e)/kappa
!
!     p_MatrixBand = -[S(e)^{T*}*C(e) + WORK]^{-1}
!  
!     -[S(e)^{T*}*C(e) + WORK] => WORK
!
!     Will call ZGETRF and ZGETRI to solve p_MatrixBand in the following equation:
!     WORK * p_MatrixBand = 1
!     ================================================================
      do j = 1, LocalNumAtoms
         kmaxj = MatrixBand(j)%kmax_kkr
         kmaxj_ns = kmaxj*nSpinCant
         nj = MatrixBand(j)%column_index-1
         ig = MatrixBand(j)%global_index ! "ig" is the global index of the corresponding atom
         nc = gid_array(ig)              ! "nc" is the column index of the block in the big matrix
         ni = MatrixBand(j)%MatrixBlock(nc)%row_index-1
         do kl = 1, kmaxj_ns
            i = KKRMatrixSizeCant*(nj+kl-1)+ni
            n = kmaxj_ns*(kl-1)
            do klp = 1, kmaxj_ns
               WORK(i+klp) = -WORK(i+klp) - stcm_g(n+klp,j)
            enddo
         enddo
      enddo
!     p_MatrixBand = WORK
      call zcopy(KKRMatrixSizeCant*KKRMatrixSizeCant,WORK,1,p_MatrixBand,1)
   endif
!
   if (NumPEsInGroup == 1) then  ! BandSizeCant = KKRMatrixSizeCant
!     ----------------------------------------------------------------
      call ZGETRF(KKRMatrixSizeCant, KKRMatrixSizeCant,               &
                  p_MatrixBand, KKRMatrixSizeCant, IPVT, INFO)
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calCrystalMatrix','Failed in ZGETRF',INFO)
!        -------------------------------------------------------------
      endif
   else
#ifdef USE_SCALAPACK
!     ----------------------------------------------------------------
      call PZGETRF(KKRMatrixSizeCant, KKRMatrixSizeCant,              &
!???  call PZGETRF(KKRMatrixSizeCant, BandSizeCant,                   &
                   p_MatrixBand, 1, 1, DESC_A, IPVT, INFO)
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calCrystalMatrix','Failed in PZGETRF',INFO)
!        -------------------------------------------------------------
      endif
#else
!     ----------------------------------------------------------------
      call ErrorHandler('calCrystalMatrix','Compiling with -DUSE_SCALAPACK is needed!')
!     ----------------------------------------------------------------
      call ZGETRF(KKRMatrixSizeCant, KKRMatrixSizeCant,               &
                  p_MatrixBand, KKRMatrixSizeCant, IPVT, INFO)
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calCrystalMatrix','Failed in ZGETRF',INFO)
!        -------------------------------------------------------------
      endif
#endif
   endif
!
   if (NumPEsInGroup == 1) then  ! BandSizeCant = KKRMatrixSizeCant
!     ----------------------------------------------------------------
      call ZGETRI( KKRMatrixSizeCant, p_MatrixBand, KKRMatrixSizeCant, IPVT, WORK, LWORK, INFO )
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calCrystalMatrix','Failed in ZGETRS',INFO)
!        -------------------------------------------------------------
      endif
   else
#ifdef USE_SCALAPACK
!     ----------------------------------------------------------------
      call PZGETRI(KKRMatrixSizeCant, p_MatrixBand, 1, 1,             &
                   DESC_A, IPVT, WORK, LWORK, IWORK, LIWORK, INFO )
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calCrystalMatrix','Failed in PZGETRS',INFO)
!        -------------------------------------------------------------
      endif
#else
!     ----------------------------------------------------------------
      call ErrorHandler('calCrystalMatrix','Compiling with -DUSE_SCALAPACK is needed!')
!     ----------------------------------------------------------------
      call ZGETRI( KKRMatrixSizeCant, p_MatrixBand, KKRMatrixSizeCant, IPVT, WORK, LWORK, INFO )
!     ----------------------------------------------------------------
      if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calCrystalMatrix','Failed in ZGETRS',INFO)
!        -------------------------------------------------------------
      endif
#endif
   endif
!
   if (method == 0) then
      do j = 1, LocalNumAtoms
         nc = MatrixBand(j)%column_index
         LOOP_n1: do n = 1, GlobalNumAtoms
            if ( MatrixBand(j)%global_index ==                        &
                 MatrixBand(j)%MatrixBlock(n)%global_index ) then
               nr = MatrixBand(j)%MatrixBlock(n)%row_index
               do i = 1, MatrixBand(j)%kmax_kkr*nSpinCant
                  p_MatrixBand(i+(nc-1+i-1)*KKRMatrixSizeCant+nr-1) = &
                     p_MatrixBand(i+(nc-1+i-1)*KKRMatrixSizeCant+nr-1) - CONE
               enddo
               exit LOOP_n1
            endif
         enddo LOOP_n1
      enddo
   endif
!
   end subroutine computeMatrixBand
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine sumIBZRotation(calculate_tau)
!  ===================================================================
   use MatrixModule, only : computeUAUtc
   use GroupCommModule, only : getGroupID, GlobalSumInGroup
   use IBZRotationModule, only : getNumIBZRotations, getIBZRotationMatrix
   use IBZRotationModule, only : getBasisRotationTable
!
   implicit none
!
   logical, intent(in) :: calculate_tau
!
   integer (kind=IntKind) :: id, jd, ig, irot, nrot, ns, is, js, kkrsz
   integer (kind=IntKind) :: ni, nj, np, kl, klp, klpp, aGID, ig_rot
   integer (kind=IntKind), pointer :: rotation_table(:,:)
!
   complex (kind=CmplxKind), target :: wtmp(kmax_kkr_max*kmax_kkr_max)
   complex (kind=CmplxKind), pointer :: w0(:,:), rotmat(:,:)
   complex (kind=CmplxKind), pointer :: kmb(:,:), tmb(:,:)
   complex (kind=CmplxKind), pointer :: matrix_diag(:,:), pm(:)
   complex (kind=CmplxKind) :: cfac
!
   nrot = getNumIBZRotations()
   cfac = CONE/real(nrot,RealKind)
   aGID = getGroupID('Unit Cell')
!
   rotation_table => getBasisRotationTable()
   if (KKRMatrixSizeCant*BandSizeCant >= tsize*GlobalNumAtoms) then
      matrix_diag => aliasArray2_c(TMP_MatrixBand,tsize,GlobalNumAtoms)
   else
      call ErrorHandler('sumIBZRotation','KKR lmax needs to be the same for all atoms')
   endif
!
   matrix_diag = CZERO
   do jd = 1, LocalNumAtoms
      ig = MatrixBand(jd)%global_index
      id = gid_array(ig)  ! Here, we only consider&rotate the diagonal blocks of the band matrix
      kkrsz = MatrixBand(jd)%MatrixBlock(id)%kmax_kkr
      ns = 0
      do js = 1, nSpinCant
         do is = 1, nSpinCant
            ns = ns + 1
!           ----------------------------------------------------------
            call zcopy(kkrsz*kkrsz,MatrixBand(jd)%MatrixBlock(id)%kau_l(1,1,ns),1, &
                       matrix_diag((ns-1)*kkrsz*kkrsz+1,ig),1)
!           ----------------------------------------------------------
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call GlobalSumInGroup(aGID,matrix_diag,tsize,GlobalNumAtoms)
!  -------------------------------------------------------------------
!
   do jd = 1, LocalNumAtoms
      ig = MatrixBand(jd)%global_index
      id = gid_array(ig)  ! Here, we only consider&rotate the diagonal blocks of the band matrix
                          ! For non-diagonal matrix blocks, there is factor of
                          ! exp(i*k_vector*(Rnm_vector)) needs to be applied to the
                          ! transformation
      kkrsz = MatrixBand(jd)%MatrixBlock(id)%kmax_kkr
      w0 => aliasArray2_c(wtmp,kkrsz,kkrsz)

!ywg  nj = MatrixBand(jd)%column_index-1
!ywg  ni = MatrixBand(jd)%MatrixBlock(id)%row_index-1
      ns = 0
      do js = 1, nSpinCant
         do is = 1, nSpinCant
            ns = ns + 1
!ywg        do kl = 1, kkrsz
!ywg           np = KKRMatrixSizeCant*(nj+kl+(js-1)*kkrsz-1)+ni+(is-1)*kkrsz
!ywg           do klp = 1, kkrsz
!ywg              kmb(klp,kl) = KKR_MatrixBand(np+klp)
!ywg           enddo
!ywg        enddo
!           kmb => MatrixBand(jd)%MatrixBlock(id)%kau_l(:,:,ns)
!           call writeMatrix('kmb before sum rot',kmb,kkrsz,kkrsz,TEN2m8)
            w0 = CZERO
            do irot = 1, nrot
               ig_rot = rotation_table(ig,irot)
               pm => matrix_diag((ns-1)*kkrsz*kkrsz+1:ns*kkrsz*kkrsz,ig_rot)
               kmb => aliasArray2_c(pm,kkrsz,kkrsz)
               rotmat => getIBZRotationMatrix('c',irot)
!              -------------------------------------------------------
!              call checkScatteringMatrixSymmetry(jd,rotmat,kmb,kkrsz, &
!                                                 getSingleScatteringMatrix, &
!                                                 site_config)
!              -------------------------------------------------------
               call computeUAUtc(rotmat,kkrsz,kkrsz,rotmat,kkrsz,cfac, &
                                 kmb,kkrsz,CONE,w0,kkrsz,WORK)
!              -------------------------------------------------------
            enddo
!ywg        do kl = 1, kkrsz
!ywg           np = KKRMatrixSizeCant*(nj+kl+(js-1)*kkrsz-1)+ni+(is-1)*kkrsz
!ywg           do klp = 1, kkrsz
!ywg              KKR_MatrixBand(np+klp) = w0(klp,kl)
!ywg           enddo
!ywg        enddo
            kmb => MatrixBand(jd)%MatrixBlock(id)%kau_l(:,:,ns)
!           ----------------------------------------------------------
            call zcopy(kkrsz*kkrsz,w0,1,kmb,1)
!           ----------------------------------------------------------
!           call writeMatrix('kmb after sum rot',kmb,kkrsz,kkrsz,TEN2m8)
!
         enddo
      enddo
   enddo
!
   if (calculate_tau) then
      matrix_diag = CZERO
      do jd = 1, LocalNumAtoms
         ig = MatrixBand(jd)%global_index
         id = gid_array(ig)  ! Here, we only consider&rotate the diagonal blocks of the band matrix
         kkrsz = MatrixBand(jd)%MatrixBlock(id)%kmax_kkr
         ns = 0
         do js = 1, nSpinCant
            do is = 1, nSpinCant
               ns = ns + 1
!              -------------------------------------------------------
               call zcopy(kkrsz*kkrsz,MatrixBand(jd)%MatrixBlock(id)%tau_l(1,1,ns),1, &
                          matrix_diag((ns-1)*kkrsz*kkrsz+1,ig),1)
!              -------------------------------------------------------
            enddo
         enddo
      enddo
!     ----------------------------------------------------------------
      call GlobalSumInGroup(aGID,matrix_diag,tsize,GlobalNumAtoms)
!     ----------------------------------------------------------------
      do jd = 1, LocalNumAtoms
         ig = MatrixBand(jd)%global_index
         id = gid_array(ig)  ! Here, we only consider&rotate the diagonal blocks of the band matrix
                             ! For non-diagonal matrix blocks, there is factor of
                             ! exp(i*k_vector*(Rnm_vector)) needs to be applied to the
                             ! transformation
         kkrsz = MatrixBand(jd)%MatrixBlock(id)%kmax_kkr
         w0 => aliasArray2_c(wtmp,kkrsz,kkrsz)
         ns = 0
         do js = 1, nSpinCant
            do is = 1, nSpinCant
               ns = ns + 1
!              tmb => MatrixBand(jd)%MatrixBlock(id)%tau_l(:,:,ns)
               w0 = CZERO
               do irot = 1, nrot
                  ig_rot = rotation_table(MatrixBand(jd)%global_index,irot)
                  pm => matrix_diag((ns-1)*kkrsz*kkrsz+1:ns*kkrsz*kkrsz,ig_rot)
                  tmb => aliasArray2_c(pm,kkrsz,kkrsz)
                  rotmat => getIBZRotationMatrix('c',irot)
!                 ----------------------------------------------------
                  call computeUAUtc(rotmat,kkrsz,kkrsz,rotmat,kkrsz,cfac, &
                                    tmb,kkrsz,CONE,w0,kkrsz,WORK)
!                 ----------------------------------------------------
               enddo
               tmb => MatrixBand(jd)%MatrixBlock(id)%tau_l(:,:,ns)
!              -------------------------------------------------------
               call zcopy(kkrsz*kkrsz,w0,1,tmb,1)
!              -------------------------------------------------------
!              call writeMatrix('tau',tmb,kkrsz,kkrsz,TEN2m8)
            enddo
         enddo
      enddo
   endif
!
   nullify(w0, rotmat, kmb, tmb, matrix_diag)
!
   end subroutine sumIBZRotation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkScatteringMatrixSymmetry(id,rotmatc,kmax_phi,      &
                                            getSingleScatteringMatrix,&
                                            site_config)
!  ===================================================================
   use MatrixModule, only : computeUAUtc
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, kmax_phi
   integer (kind=IntKind), intent(in) :: site_config(:)
   integer (kind=IntKind) :: kl, klp
!
   real (kind=RealKind), parameter :: tol = ten2m10
!
   complex (kind=CmplxKind), intent(in) :: rotmatc(:,:)
   complex (kind=CmplxKind), pointer :: tmat(:,:), OmegaHat(:,:)
   complex (kind=CmplxKind), pointer :: smat(:,:), cmat(:,:)
   complex (kind=CmplxKind), allocatable :: tmat_rot(:,:), OmegaHat_rot(:,:)
   complex (kind=CmplxKind), allocatable :: smat_rot(:,:), cmat_rot(:,:)
!
   interface
      function getSingleScatteringMatrix(smt,spin,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in), optional :: spin, site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   allocate(tmat_rot(kmax_phi,kmax_phi),OmegaHat_rot(kmax_phi,kmax_phi))
   allocate(smat_rot(kmax_phi,kmax_phi),cmat_rot(kmax_phi,kmax_phi))
!
   tmat => getSingleScatteringMatrix('T-Matrix',spin=1,site=id,       &
                                     atom=site_config(id))
   OmegaHat => getSingleScatteringMatrix('OmegaHatInv-Matrix',spin=1, &
                                         site=id,atom=site_config(id))
   smat => getSingleScatteringMatrix('Sine-Matrix',spin=1,site=id,    &
                                     atom=site_config(id))
   cmat => getSingleScatteringMatrix('Cosine-Matrix',spin=1,site=id,  &
                                     atom=site_config(id))
   call computeUAUtc(rotmatc,kmax_phi,kmax_phi,rotmatc,kmax_phi,CONE, &
                     tmat,kmax_phi,CZERO,tmat_rot,kmax_phi,WORK)
   call computeUAUtc(rotmatc,kmax_phi,kmax_phi,rotmatc,kmax_phi,CONE, &
                     OmegaHat,kmax_phi,CZERO,OmegaHat_rot,kmax_phi,WORK)
   call computeUAUtc(rotmatc,kmax_phi,kmax_phi,rotmatc,kmax_phi,CONE, &
                     smat,kmax_phi,CZERO,smat_rot,kmax_phi,WORK)
   call computeUAUtc(rotmatc,kmax_phi,kmax_phi,rotmatc,kmax_phi,CONE, &
                     cmat,kmax_phi,CZERO,cmat_rot,kmax_phi,WORK)
!
   do kl = 1, kmax_phi
      do klp = 1, kmax_phi
         if (abs(tmat_rot(klp,kl) - tmat(klp,kl)) > tol) then
            write(6,'(a,3i5)')'id,kl,klp = ',id,kl,klp
            call ErrorHandler('checkScatteringMatrixSymmetry','tmat_rot <> tmat', &
                              tmat_rot(klp,kl),tmat(klp,kl))
!        else if (abs(smat_rot(klp,kl) - smat(klp,kl)) > tol) then
!           write(6,'(a,3i5)')'id,kl,klp = ',id,kl,klp
!           call ErrorHandler('checkScatteringMatrixSymmetry','smat_rot <> smat', &
!                             smat_rot(klp,kl),smat(klp,kl))
!        else if (abs(cmat_rot(klp,kl) - cmat(klp,kl)) > tol) then
!           write(6,'(a,3i5)')'id,kl,klp = ',id,kl,klp
!           call ErrorHandler('checkScatteringMatrixSymmetry','cmat_rot <> cmat', &
!                             cmat_rot(klp,kl),cmat(klp,kl))
!        else if (abs(OmegaHat_rot(klp,kl) - OmegaHat(klp,kl)) > 1000.0d0*tol) then
!           write(6,'(a,3i5)')'id,kl,klp = ',id,kl,klp
!           call ErrorHandler('checkScatteringMatrixSymmetry','OmegaHat_rot <> OmegaHat', &
!                             OmegaHat_rot(klp,kl),OmegaHat(klp,kl))
         endif
      enddo
   enddo
!
   deallocate(tmat_rot,OmegaHat_rot)
   deallocate(smat_rot,cmat_rot)
!
   end subroutine checkScatteringMatrixSymmetry
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTau(local_id,global_i,global_j,dsize) result(tau)
!  ===================================================================
!  local_i  = local site index
!  global_i = global site index
!  global_j = global site index
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: local_id
   integer (kind=IntKind), intent(in), optional :: global_i, global_j
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: ni, ig, jd, id
!
   complex (kind=CmplxKind), pointer :: tau(:,:,:)
!
   if (present(local_id)) then
      if (local_id < 1 .or. local_id > LocalNumAtoms) then
         call ErrorHandler('getTau','invalid local index i',local_id)
      else if (present(global_j)) then
         call ErrorHandler('getTau','The function call is confusing')
      endif
!
      jd = local_id
      if (present(global_i)) then
         ig = global_i
      else
!        take the diagonal block associated with the atom with local_id
         ig = MatrixBand(jd)%global_index
      endif
   else if (present(global_i) .and. present(global_j)) then
      if (global_i < 1 .or. global_i > GlobalNumAtoms) then
         call ErrorHandler('getTau','invalid global index i',global_i)
      else if (global_j < 1 .or. global_j > GlobalNumAtoms) then
         call ErrorHandler('getTau','invalid global index j',global_j)
      endif
!     take the appropriate matrix block associated with the atom specified
!     by the global index pair: (global_i, global_j)
      jd = 0
      LOOP_id: do id = 1, LocalNumAtoms
         if (MatrixBand(id)%global_index == global_j) then
            jd = id
            exit LOOP_id
         endif
      enddo LOOP_id
      if (jd == 0) then
         call ErrorHandler('getTau',                                  &
              'The global j index is not associated with a local atom',global_j)
      endif
      ig = global_i
   else if (present(global_i) .or. present(global_j)) then
      call ErrorHandler('getTau','The function call is confusing')
   else 
!     take the upper-left corner block
      jd = 1
      ig = MatrixBand(jd)%global_index
   endif
!
   ni = gid_array(ig)
!
   tau => MatrixBand(jd)%MatrixBlock(ni)%tau_l
!
   if (present(dsize)) then
      dsize = MatrixBand(jd)%MatrixBlock(ni)%kmax_kkr*nSpinCant
   endif
!
   end function getTau
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKau(local_id,global_i,global_j,dsize) result(kau)
!  ===================================================================
!  local_i  = local site index
!  global_i = global site index
!  global_j = global site index
!  *******************************************************************
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: local_id
   integer (kind=IntKind), intent(in), optional :: global_i, global_j
   integer (kind=IntKind), intent(out), optional :: dsize
   integer (kind=IntKind) :: ni, ig, id, jd
!
   complex (kind=CmplxKind), pointer :: kau(:,:,:)
!
   if (present(local_id)) then
      if (local_id < 1 .or. local_id > LocalNumAtoms) then
         call ErrorHandler('getKau','invalid local index i',local_id)
      else if (present(global_j)) then
         call ErrorHandler('getKau','The function call is confusing')
      endif
!
      jd = local_id
      if (present(global_i)) then
         ig = global_i
      else
!        take the diagonal block associated with the atom with local_id
         ig = MatrixBand(jd)%global_index
      endif
   else if (present(global_i) .and. present(global_j)) then
      if (global_i < 1 .or. global_i > GlobalNumAtoms) then
         call ErrorHandler('getKau','invalid global index i',global_i)
      else if (global_j < 1 .or. global_j > GlobalNumAtoms) then
         call ErrorHandler('getKau','invalid global index j',global_j)
      endif
!     take the appropriate matrix block associated with the atom specified
!     by the global index pair: (global_i, global_j)
      jd = 0
      LOOP_id: do id = 1, LocalNumAtoms
         if (MatrixBand(id)%global_index == global_j) then
            jd = id
            exit LOOP_id
         endif
      enddo LOOP_id
      if (jd == 0) then
         call ErrorHandler('getKau',                                  &
              'The global j index is not associated with a local atom',global_j)
      endif
      ig = global_i
   else if (present(global_i) .or. present(global_j)) then
      call ErrorHandler('getKau','The function call is confusing')
   else 
!     take the upper-left corner block
      jd = 1
      ig = MatrixBand(jd)%global_index
   endif
!
   ni = gid_array(ig)
!
   kau => MatrixBand(jd)%MatrixBlock(ni)%kau_l
!
   if (present(dsize)) then
      dsize = MatrixBand(jd)%MatrixBlock(ni)%kmax_kkr*nSpinCant
   endif
!
   end function getKau
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printCrystalMatrix(i)  ! "i" is the global index
!  ===================================================================
   use Atom2ProcModule, only : getAtom2ProcInGroup, getLocalIndex
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: i
   integer (kind=IntKind) :: ni, ig, id, jd, js
   integer (kind=IntKind) :: kmaxi, kmaxj, kmaxi_ns, kmaxj_ns
!
   complex (kind=CmplxKind), pointer :: kau(:,:,:)
!
   if (present(i)) then
      if (i < 1 .or. i > LocalNumAtoms) then
         call ErrorHandler('printCrystalMatrix','invalid local index i',i)
      else
         id = i
      endif
   else 
      id = 1
   endif
!
   if (getAtom2ProcInGroup(id) /= MyPEinGroup) then
      call ErrorHandler('printCrystalMatrix','getAtom2ProcInGroup(j) /= MyPEinGroup',&
                         getAtom2ProcInGroup(jd), MyPEinGroup)
   endif
!
   ig = MatrixBand(id)%global_index
   ni = gid_array(ig)
   kmaxj = MatrixBand(id)%kmax_kkr
   kmaxj_ns = kmaxj*nSpinCant
!
   kau => MatrixBand(id)%MatrixBlock(ni)%kau_l
!
   write(6,'(/,a)')'****************************************'
   write(6,'( a )')'*    Output from printCrystalMatrix    *'
   write(6,'(a,/)')'****************************************'
   write(6,'(80(''=''))')
!  
   write(6,*)'energy   ::',energy
   write(6,'(80(''=''))')
!
   write(6,'(/,a,i4,2x,i4,/)')"    Sites i :: ",id
   do js = 1, nSpinCant*nSpinCant
      call writeMatrix( 'Kau-matrix in Local Frame of Reference in Spin Space ::', &
                        kau(1:kmaxj_ns,1:kmaxj_ns,js),kmaxj_ns,kmaxj_ns )
   enddo
!
   end subroutine printCrystalMatrix
!  ===================================================================
end module CrystalMatrixModule
