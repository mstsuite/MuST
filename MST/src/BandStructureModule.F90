module BandStructureModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1, TEN2m6, HALF
   use PhysParamModule, only : LightSpeed
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
   use PublicTypeDefinitionsModule, only : MatrixBandStruct, ScmBlockStruct
!
public :: initBandStructure,   &
          endBandStructure,    &
          calBandStructure
!
private
!
   logical :: isRelativistic = .false.
   real (kind=RealKind), parameter :: Me = 0.5d0 !xianglin
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: GlobalNumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: nSpinCant
   integer (kind=IntKind) :: lmax_kkr_max, kmax_kkr_max, tsize
   integer (kind=IntKind) :: KKRMatrixSize
   integer (kind=IntKind) :: BandSize
   integer (kind=IntKind) :: KKRMatrixSizeCant
   integer (kind=IntKind) :: BandSizeCant
!
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
!
   integer (kind=IntKind), allocatable :: ip_array(:) ! Relates a matrix block
                                                      ! row index to the mapped processor index (0, 1, ...)
   integer (kind=IntKind), allocatable :: il_array(:) ! Relates a matrix block
                                                      ! row index to the local atom index on the mapped processor
   integer (kind=IntKind), allocatable :: id_array(:) ! Relates a matrix block
                                                      ! row index to the global index of the corresponding atom
   integer (kind=IntKind), allocatable :: jd_array(:) ! Relates a local atom
                                                      ! index to the global index of the atom
   integer (kind=IntKind), allocatable :: gid_array(:)! Relates a global index
                                                      ! to the corresponding matrix block row index
   integer (kind=IntKind), allocatable :: lmaxi_array(:)
   integer (kind=IntKind), allocatable :: lmaxj_array(:)
!
   type (ScmBlockStruct), allocatable :: sc_blocks(:,:)
!
   type (MatrixBandStruct), allocatable :: MatrixBand(:)  ! Matrix column band
!
   complex (kind=CmplxKind), allocatable, target :: KKR_MatrixBand(:)
   complex (kind=CmplxKind), allocatable, target :: strconrel(:,:)
!
   complex (kind=CmplxKind), allocatable, target :: sine_g(:,:) ! Sine matrix in global frame
   complex (kind=CmplxKind), allocatable, target :: stcm_g(:,:) ! kappa*S^{T*}*C matrix in global frame
   complex (kind=CmplxKind), allocatable, target :: tmat_g(:,:) ! t-matrix in global frame
   complex (kind=CmplxKind), pointer :: cosine_g(:,:)  ! Cosine matrix in global frame
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initBandStructure(nla,cant,lmax_kkr,rel,istop,iprint)
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
   integer (kind=IntKind) :: lmaxi, kmaxi, status, n, i, j
!
   character (len=20) :: sname = "initBandStructure"
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: global_posi(:,:)
!
   call StopHandler(sname,'To be implemented')
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
   allocate ( MatrixBand(LocalNumAtoms), STAT=status )
!
   allocate( id_array(GlobalNumAtoms), jd_array(LocalNumAtoms) )
   allocate( ip_array(GlobalNumAtoms), il_array(GlobalNumAtoms) )
   allocate( lmaxi_array(GlobalNumAtoms), lmaxj_array(LocalNumAtoms) )
   allocate( gid_array(GlobalNumAtoms), sc_blocks(GlobalNumAtoms,LocalNumAtoms) )
!
   BandSize = 0
   n = 0
   do j = 1, LocalNumAtoms
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
   allocate ( KKR_MatrixBand(KKRMatrixSizeCant*BandSizeCant) )
   if (isRelativistic) then
      allocate( strconrel(kmax_kkr_max*nSpinCant,kmax_kkr_max*nSpinCant) )
   endif
!
   allocate( sine_g(tsize,LocalNumAtoms) )
   allocate( tmat_g(tsize,LocalNumAtoms) )
   allocate( stcm_g(tsize,LocalNumAtoms) )
   sine_g = CZERO ! It stores the sine matrix in the global frame
   tmat_g = CZERO ! It stores the t-matrix in the global frame
   stcm_g = CZERO ! It stores the S^{T*}*C in the global frame
!
   print_level = maxval(iprint)
!
   end subroutine initBandStructure
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endBandStructure()
!  ===================================================================
   implicit none
!
   call StopHandler('endBandStructure','To be implemented')
!
   deallocate ( MatrixBand )
!
   deallocate( id_array, jd_array )
   deallocate( ip_array, il_array )
   deallocate( lmaxi_array, lmaxj_array )
   deallocate( gid_array, sc_blocks )
!
   deallocate ( KKR_MatrixBand )
   if (isRelativistic) then
      deallocate( strconrel )
   endif
!
   deallocate( sine_g )
   deallocate( tmat_g )
   deallocate( stcm_g )
!
   end subroutine endBandStructure
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calBandStructure(is,eb,et,ne,kpts,nk)
!  ===================================================================
   use StrConstModule, only : getStrConstMatrix
!
   use SSSolverModule, only : solveSingleScattering
!
   implicit none
!
   character (len=20) :: sname = "calBandStructure"
!
   integer (kind=IntKind), intent(in) :: is, ne, nk
   integer (kind=IntKind) :: ik, col, row, ie, id, js, ns
!
   real (kind=RealKind), intent(in) :: kpts(3,nk)
   real (kind=RealKind) :: kvec(3), de
!
   real (kind=RealKind), intent(in) :: eb, et ! Find energy eigenvalues
                                              ! within (eb, et)
!
   complex (kind=CmplxKind), pointer :: scm(:,:)
   complex (kind=CmplxKind) :: energy, kappa
!
   do ik = 1, nk    ! Loop over the k-point mesh
      kvec(1:3) = kpts(1:3,ik)
      write(6,'(a,3f12.6)')'kvec(1:3) = ',kvec(1:3)
!
      de = (et-eb)/real(ne-1,kind=RealKind)
      do ie = 1, ne ! Loopp over the energy mesh
         energy = eb + (ie-1)*de
         if (abs(energy) < TEN2m6) then
            energy = energy + de*HALF
         endif
         write(6,'(a,f12.6)')'energy = ',real(energy,kind=RealKind)
!
         if (isRelativistic) then !xianglin
            kappa = sqrt(2.d0*Me*energy + energy**2/LightSpeed**2)
         else
            kappa = sqrt(energy)
         endif
!
         do col = 1, LocalNumAtoms
            do row = 1, GlobalNumAtoms
!              -------------------------------------------------------
               scm => getStrConstMatrix(kvec,kappa,id_array(row),jd_array(col), &
                                        lmaxi_array(row),lmaxj_array(col))
!              -------------------------------------------------------
               sc_blocks(row,col)%strcon_matrix = scm
            enddo
         enddo
!
         do id = 1, LocalNumAtoms
            do js = 1, nSpinCant
               ns = max(js,is)
!              -------------------------------------------------------
               call solveSingleScattering(ns,id,energy,CZERO)
!              -------------------------------------------------------
            enddo
         enddo
!
!        =============================================================
!        setup S- and C- matrix in global frame
!        -------------------------------------------------------------
         call setupSCinGlobalFrame()
!        -------------------------------------------------------------
!
         KKR_MatrixBand = CZERO
!        =============================================================
!        Compute the KKR matrix (kappa*C+B*S), which is stored in KKR_MatrixBand
!        -------------------------------------------------------------
         call computeKKRMatrix(energy)
!        -------------------------------------------------------------
!
!        =============================================================
!        At this point, KKR_MatrixBand is calculated. Note:
!           det[ KKR_Matrix ] = CZERO gives rise to the band structure
!           KKR_Matrix is a KKRMatrixSizeCant x KKRMatrixSizeCant 
!           matrix, and is divided into multiple Bands of columns so
!           that it is disibuted on multiple processors. In other words,
!           each band of columns, called KKR_MatrixBand, is allocated 
!           on a processor as follows:
!               allocate ( KKR_MatrixBand(KKRMatrixSizeCant*BandSizeCant) )
!           Physically, KKR_MatrixBand is a matrix (2-D array) with 
!           row size = KKRMatrixSizeCant:
!
!                         [ ( 1, 1 )                 ( 1, 2 )               ...      ( 1, BandSizeCant )          ]
!                         | ( 2, 1 )                 ( 2, 2 )               ...      ( 2, BandSizeCant )          |
!                         |    .                        .                   ...         .                         |
!        KKR_MatrixBand = |    .                        .                   ...         .                         |
!                         |    .                        .                   ...         .                         |
!                         [ (KKRMatrixSizeCant,1)   (KKRMatrixSizeCant,2)   ...  (KKRMatrixSizeCant,BandSizeCant) ]
!
!           Of course, for thye number of processors = 1, BandSizeCant = KKRMatrixSizeCant
!        =============================================================
!
!        Leo: The following part needs to be worked on...
!        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 

!        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        End of the change.
!        =============================================================
      enddo
   enddo
!
   end subroutine calBandStructure
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupSCinGlobalFrame()
!  ===================================================================
   use WriteMatrixModule,  only : writeMatrix
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SpinRotationModule, only : rotateLtoG
!
   use SSSolverModule, only : getScatteringMatrix
!
   use RelSSSolverModule, only : getRelScatteringMatrix
!
   implicit none
!
   integer (kind=IntKind) :: t0size, kkri_ns, i, kmax_kkr
!
   complex (kind=CmplxKind), pointer :: sm1(:,:), sm2(:,:)
   complex (kind=CmplxKind), pointer :: cm1(:,:), cm2(:,:)
   complex (kind=CmplxKind), pointer :: pm(:), gmat(:,:)
!
   cosine_g => stcm_g   ! Use stcm_g as the space for storing the cosine matrix
!
   if (isRelativistic) then !xianglin in between
      do i = 1, LocalNumAtoms
!        ================================================================
!        Obtain the Jinv-matrix, Sine-Matrix, and t-matrix in Global frame
!        ================================================================
         kmax_kkr = MatrixBand(i)%kmax_kkr
         kkri_ns =  kmax_kkr*nSpinCant
         t0size = kmax_kkr*kmax_kkr*nSpinCant*nSpinCant
         sm1 => getRelScatteringMatrix('Sine-Matrix',i)
         cm1 => getRelScatteringMatrix('Cosine-Matrix',i)
!        ----------------------------------------------------------------
         call zcopy( t0size, sm1, 1, sine_g(1,i), 1 )     !save s into sine_g
         call zcopy( t0size, cm1, 1, cosine_g(1,i), 1 )   !save jinv(Method 0) or s^t(Method 1) into jinv_g 
!        ----------------------------------------------------------------
      enddo
   else
      do i = 1, LocalNumAtoms
!        =============================================================
!        Obtain the Cosine-matrix, and Sine-Matrix t-matrix in Global frame
!        =============================================================
         kmax_kkr = MatrixBand(i)%kmax_kkr
         t0size = kmax_kkr*kmax_kkr
         if ( nSpinCant == 2 ) then
            kkri_ns =  kmax_kkr*nSpinCant
!           =============================================================
!           calculate sine_g and cosine_g in global frame of reference.
!           =============================================================
            sm1 => getScatteringMatrix('Sine-Matrix',1,i)
            sm2 => getScatteringMatrix('Sine-Matrix',2,i)
            pm => sine_g(:,i)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           -------------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, sm1, sm2, gmat)
!           -------------------------------------------------------------
            cm1 => getScatteringMatrix('Cosine-Matrix',1,i)
            cm2 => getScatteringMatrix('Cosine-Matrix',2,i)
            pm => cosine_g(:,i)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           -------------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, cm1, cm2, gmat)
!           -------------------------------------------------------------
         else
            sm1 => getScatteringMatrix('Sine-Matrix',1,i)
            cm1 => getScatteringMatrix('Cosine-Matrix',1,i)
!           -------------------------------------------------------------
            call zcopy( t0size, sm1, 1, sine_g(1,i), 1 )
            call zcopy( t0size, cm1, 1, cosine_g(1,i), 1 )
!           -------------------------------------------------------------
         endif
      enddo
   endif
!
   nullify(sm1,sm2,cm1,cm2,pm,gmat)
!
   end subroutine setupSCinGlobalFrame
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeKKRMatrix(energy)
!  ===================================================================
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind) :: kmaxj, kmaxj_ns, kmaxi, kmaxi_ns, t0size
   integer (kind=IntKind) :: j, nj, ni, i, is, ig, nc, kl
!
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind) :: kappa
!
   complex (kind=CmplxKind), pointer :: strcon(:,:)
   complex (kind=CmplxKind), pointer :: p_cosinej(:)
   complex (kind=CmplxKind), pointer :: p_sinej(:)
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
!    KKR_MatrixBand = [kappa*C(e) + B(e,k) * S(e)]
!  ===================================================================
   kappa = sqrt(energy)
   KKR_MatrixBand = CZERO
   do j = 1, LocalNumAtoms
      p_sinej => sine_g(:,j)
      p_cosinej => cosine_g(:,j)
!
      kmaxj = MatrixBand(j)%kmax_kkr
      kmaxj_ns = kmaxj*nSpinCant
      nj = MatrixBand(j)%column_index-1
      ig = MatrixBand(j)%global_index ! "ig" is the global index of the corresponding atom
      nc = gid_array(ig)              ! "nc" is the column index of the block in the big matrix
!
      ni = MatrixBand(j)%MatrixBlock(nc)%row_index-1
      do kl = 1, kmaxj_ns
!        -------------------------------------------------------------
         call zaxpy(kmaxj_ns,kappa,p_cosinej(kmaxj_ns*(kl-1)+1),1,    &
                    KKR_MatrixBand(KKRMatrixSizeCant*(nj+kl-1)+ni+1),1)
!        -------------------------------------------------------------
      enddo
!
      do i = 1, GlobalNumAtoms        ! "i" is the row index of the matrix block
         kmaxi = MatrixBand(j)%MatrixBlock(i)%kmax_kkr
         kmaxi_ns = kmaxi*nSpinCant
         t0size = kmaxi_ns*kmaxi_ns
         ni = MatrixBand(j)%MatrixBlock(i)%row_index-1
!
         strcon => sc_blocks(i,j)%strcon_matrix(:,:)
!
         if (isRelativistic) then
!           ----------------------------------------------------------
            call convertGijToRel(strcon, strconrel, kmaxi, kmaxj, energy)
!           ----------------------------------------------------------
            strcon => strconrel
!           ----------------------------------------------------------
            call zgemm('n', 'n', kmaxi_ns, kmaxj_ns, kmaxj_ns, CONE,  &
                       strcon, kmaxi_ns, p_sinej, kmaxj_ns, CONE,     &
                       KKR_MatrixBand, KKRMatrixSizeCant)
!           ----------------------------------------------------------
         else
            do is = 1, nSpinCant
!              -------------------------------------------------------
               call zgemm('n', 'n', kmaxi, kmaxj, kmaxj, CONE,        &
                          strcon, kmaxi,                              &
                          p_sinej((is-1)*kmaxj_ns*kmaxj+1), kmaxj_ns, &
                          CONE,                                       &
                          KKR_MatrixBand(KKRMatrixSizeCant*nj+ni+1), KKRMatrixSizeCant)
!              -------------------------------------------------------
            enddo
         endif
      enddo
   enddo
!
   nullify(p_sinej, p_cosinej)
!
   end subroutine computeKKRMatrix
!  ===================================================================
end module BandStructureModule
