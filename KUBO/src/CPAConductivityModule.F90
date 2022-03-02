 module CPAConductivityModule
   use KindParamModule, only : IntKind, QuadRealKind, QuadCmplxKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : PI, ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors

public :: initCPAConductivity, &
          calSigmaTilde0,   &
          calSigmaTilde1VC, &
          calOmegaMatrixCPA,   &
          calChiMatrixCPA,     &
          calChiIntegral,      &
          calVertexCorrectionMatrixCPA, &
          endCPAConductivity

private
   logical :: vertex_corr
   integer (kind=IntKind) :: num_atoms, local_index, global_index
   integer (kind=IntKind) :: lmax_kkr, dsize, num_species, spin_pola
   integer (kind=IntKind) :: LWORK, LIWORK, dbg
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind), allocatable :: IWORK(:)
   real (kind=RealKind) :: omega
   real (kind=RealKind), allocatable :: species_content(:)
   complex (kind=CmplxKind) :: eval
   complex (kind=CmplxKind), pointer :: tau_c(:,:), tmat(:,:), tc(:,:)
   complex (kind=CmplxKind), allocatable :: tcc(:,:), tau_cc(:,:)
   complex (kind=CmplxKind), allocatable :: tac(:,:), tau1(:,:), tau2(:,:), tbc(:,:)
   complex (kind=CmplxKind), allocatable :: xa(:,:), xac(:,:)
   complex (kind=CmplxKind), allocatable :: temp1(:,:), temp2(:,:),  &
                                              temp3(:,:), temp4(:,:)
   complex (kind=CmplxKind), allocatable :: J1avg(:), J2avg(:)
   complex (kind=CmplxKind), allocatable :: X(:,:,:), W(:,:,:), A(:,:,:)
   complex (kind=CmplxKind), allocatable, target :: wtmp(:), wtmpsym(:), wint(:)
   complex (kind=CmplxKind), allocatable, target :: TMP_MatrixBand(:), WORK(:)
   complex (kind=CmplxKind), allocatable, target :: tmat_g(:)
   complex (kind=CmplxKind), allocatable :: tmb(:,:), tmbsym(:,:)


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

contains
!
   include '../lib/arrayTools.F90'
!  ==================================================================

!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCPAConductivity(n, is, kmax, efermi, LocalNumAtoms)
!  ==================================================================
 
   use SSSolverModule, only : solveSingleScattering 
   use CPAMediumModule, only : computeCPAMedium, getCPAMatrix, getSingleSiteTmat
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
   use Atom2ProcModule, only : getGlobalIndex
   use SystemVolumeModule, only : getAtomicVPVolume
   use CurrentMatrixModule, only : calCurrentMatrix
   use KuboDataModule, only : isFermiEnergyRealPart, includeVertexCorrections, &
     getFermiEnergyRealPart, useStepFunctionForSigma, getFermiEnergyImagPart
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
   use GroupCommModule, only : syncAllPEsInGroup

   integer (kind=IntKind), intent(in) :: n, is, kmax, LocalNumAtoms
   real (kind=RealKind), intent(in) :: efermi

   integer (kind=IntKind) :: ic, pot_type
   real (kind=RealKind) :: delta

   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)

   local_index = n
   global_index = getGlobalIndex(local_index)
   num_atoms = LocalNumAtoms
   
   spin_pola = is
   dsize = kmax
   global_index = getGlobalIndex(local_index)
   num_species = getLocalNumSpecies(local_index)
   omega = getAtomicVPVolume(local_index)
   
   vertex_corr = includeVertexCorrections()
   delta = getFermiEnergyImagPart()
   pot_type = useStepFunctionForSigma()

   if (isFermiEnergyRealPart()) then
     eval = getFermiEnergyRealPart() + SQRTm1*delta
   else
     eval = efermi + SQRTm1*delta
   endif

!  -------------------------------------------------------------------------------
   call solveSingleScattering(spin=spin_pola,site=local_index,e=eval,vshift=CZERO)
!  -------------------------------------------------------------------------------
   call computeCPAMedium(eval)
!  ------------------------------------------------------------------
   call calCurrentMatrix(local_index,spin_pola,eval,pot_type,3)
!  ------------------------------------------------------------------
   
   
   allocate(species_content(num_species))
   do ic = 1, num_species
     species_content(ic) = getLocalSpeciesContent(local_index, ic)
   enddo
   
   tau_c => getCPAMatrix('Tau',site=local_index,atom=0)
   tc => getSingleSiteTmat('TInv-Matrix', spin=spin_pola, site=local_index, atom=0)
   tmat => getSingleSiteTmat('T-Matrix', spin=spin_pola, site=local_index, atom=0)
   
   allocate(tau_cc(dsize, dsize), &
     tcc(dsize, dsize), tac(dsize, dsize), tau1(dsize, dsize), tau2(dsize, dsize), &
     xa(dsize, dsize), xac(dsize, dsize), temp1(dsize, dsize), temp2(dsize, dsize), &
     temp3(dsize, dsize), temp4(dsize, dsize), J1avg(dsize*dsize), J2avg(dsize*dsize))
   
   tau_cc = conjg(tau_c)
   tcc = conjg(tc)
   tac = CZERO; tau1 = CZERO; tau2 = CZERO; xa = CZERO; xac = CZERO
   temp1 = CZERO; temp2 = CZERO; temp3 = CZERO; temp4 = CZERO
   J1avg = CZERO; J2avg = CZERO
   
   allocate(X(dsize*dsize, dsize*dsize, 4))
   allocate(A(dsize*dsize, dsize*dsize, 4))
   allocate(W(dsize*dsize, dsize*dsize, 4))
   
   X = CZERO; A = CZERO; W = CZERO
   
   allocate(wtmp(dsize*dsize), wtmpsym(dsize*dsize), wint(dsize*dsize))
   allocate(TMP_MatrixBand(dsize*dsize), WORK(dsize*dsize), tmat_g(dsize*dsize))
   allocate(tmb(dsize, dsize), tmbsym(dsize, dsize))
   allocate( IPVT(1:dsize+num_atoms*dsize) )
   LWORK = dsize*dsize
   LIWORK = 2*dsize
   allocate( IWORK(1:LIWORK) )
   WORK = CZERO; IWORK = 0

   wtmp = CZERO; wtmpsym = CZERO; wint = CZERO;
   TMP_MatrixBand = CZERO
   tmb = CZERO; tmbsym = CZERO;
!  -------------------------------------------------------------
   call zcopy(dsize*dsize,tmat(1,1),1,tmat_g(1),1)
!  -------------------------------------------------------------
   if (NumPEsInGroup == 1) then  ! ScaLapack will not be used for 1 process case
      return
   endif
!
#ifdef USE_SCALAPACK
!  ==================================================================
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
   call DESCINIT( DESC_A, dsize, dsize, dsize, dsize, 0, 0, ICTXT, dsize, INFO )
!  -------------------------------------------------------------------
#endif

   end subroutine initCPAConductivity
!  ==================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTilde0(dir1, dir2, caltype) result(sigma0)
!  =================================================================== 
!
!  Calculates sigma_tilde0(z_1, z_2) at z_1 = z_2 = e_F + i*delta
!  Can do calculation for both single site CPA and supercell case
!
   use CurrentMatrixModule, only : getJMatrix

   integer (kind=IntKind), intent(in) :: dir1, dir2, caltype
   integer (kind=IntKind) :: ic, ic1, L
   real (kind=RealKind) :: c_a, c_b, coeff
   complex (kind=CmplxKind), pointer :: Ja(:,:), Jb(:,:), Jc(:,:)
   complex (kind=CmplxKind) :: sigma0

   sigma0 = CZERO

   do ic = 1, num_species
    do ic1 = 1, num_species
      c_a = species_content(ic)
      c_b = species_content(ic1)
      coeff = -(c_a*c_b)/(PI*omega)
      if (caltype == 1) then
        Ja => getJMatrix(local_index, ic, spin_pola, dir1, 1, 1)
        tau1 = tau_c; tau2 = tau_c
      else if (caltype == 2) then
        Ja => getJMatrix(local_index, ic, spin_pola, dir1, 3, 1)
        tau1 = tau_c; tau2 = tau_cc
      else if (caltype == 3) then
        Ja => getJMatrix(local_index, ic, spin_pola, dir1, 2, 1)
        tau1 = tau_cc; tau2 = tau_c
      else if (caltype == 4) then
        Ja => getJMatrix(local_index, ic, spin_pola, dir1, 4, 1)
        tau1 = tau_cc; tau2 = tau_cc
      else
        call ErrorHandler('calSigmaTildeCPA0','Incorrect caltype (1-4)', caltype)
      endif
      Jb => getJMatrix(local_index, ic, spin_pola, dir2, caltype, 0)
      temp4 = Jb
!     ---------------------------------------------------------------------
      call zgemm('N', 'n', dsize, dsize, dsize, CONE, temp4, &
         dsize, tau2, dsize, CZERO, temp1, dsize)
!     ---------------------------------------------------------------------
      call zgemm('N', 'n', dsize, dsize, dsize, CONE, tau1, dsize,  &
         temp1, dsize, CZERO, temp2, dsize)
!     ---------------------------------------------------------------------
      call zgemm('N', 'n', dsize, dsize, dsize, CONE, Ja, &
         dsize, temp2, dsize, CZERO, temp3, dsize)
!     ----------------------------------------------------------------------
      do L = 1, dsize
        sigma0 = sigma0 + coeff*temp3(L,L)
      enddo
      nullify(Ja, Jb, Jc)
    enddo
   enddo

   end function calSigmaTilde0
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTilde1VC(dir1, dir2, caltype) result(sigma1)
!  ===================================================================
   
   use CurrentMatrixModule, only : getJMatrix1D
   
   integer (kind=IntKind), intent(in) :: dir1, dir2, caltype
   integer (kind=IntKind) :: ic, K, L1, L4, K1, K1_t
   real (kind=RealKind) :: c_a, c_b, coeff
   complex (kind=CmplxKind) :: sigma1 
   complex (kind=CmplxKind), pointer :: J1(:), J2(:)
   
   J1avg = CZERO; J2avg = CZERO; sigma1 = CZERO
   coeff = -1.0/(PI*omega)
   
   do ic = 1, num_species
     c_a = species_content(ic)
     if (caltype == 2) then
       J1 => getJMatrix1D(local_index, ic, spin_pola, dir1, 3)
     else if (caltype == 3) then
       J1 => getJMatrix1D(local_index, ic, spin_pola, dir1, 2)
     else
       J1 => getJMatrix1D(local_index, ic, spin_pola, dir1, caltype)
     endif
     J2 => getJMatrix1D(local_index, ic, spin_pola, dir2, caltype)
     J1avg = J1avg + c_a*J1
     J2avg = J2avg + c_a*J2
   enddo
  
   do K = 1, dsize**2
     do L1 = 1, dsize
       do L4 = 1, dsize
         K1 = dsize*(L1 - 1) + L4
         K1_t = dsize*(L4 - 1) + L1
         sigma1 = sigma1 + coeff*J1avg(K1_t)*A(K1, K, caltype)*J2avg(K)
       enddo
     enddo
   enddo

   end function calSigmaTilde1VC
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calOmegaMatrixCPA()
!  ===================================================================

   use SSSolverModule, only : getScatteringMatrix
   use MatrixModule, only : computeAprojB

   integer (kind=IntKind) :: ic, L1, L2, L3, L4, K1, K2
   real (kind=RealKind) :: c_a
   complex (kind=CmplxKind), pointer :: ta(:,:)


   do ic = 1, num_species
     tac = CZERO; temp1 = CZERO; temp2 = CZERO; xa = CZERO; xac = CZERO
     c_a = species_content(ic)
     ta => getScatteringMatrix('TInv-Matrix', spin=spin_pola, site=local_index, atom=ic)
     tac = conjg(ta)
     temp1 = ta - tc
     temp2 = tac - tcc
!    ------------------------------------------------------------
     call computeAprojB('L', dsize, temp1, tau_c, xa)
     call computeAprojB('L', dsize, temp2, tau_cc, xac)
!    ------------------------------------------------------------
     do L4 = 1, dsize
       do L3 = 1, dsize
         do L2 = 1, dsize
           do L1 = 1, dsize
             K1 = dsize*(L1 - 1) + L4
             K2 = dsize*(L2 - 1) + L3
             W(K1,K2,1) = W(K1,K2,1) - c_a*xa(L1,L2)*xa(L3,L4)
             W(K1,K2,2) = W(K1,K2,2) - c_a*xa(L1,L2)*xac(L3,L4)
             W(K1,K2,3) = W(K1,K2,3) - c_a*xac(L1,L2)*xa(L3,L4)
             W(K1,K2,4) = W(K1,K2,4) - c_a*xac(L1,L2)*xac(L3,L4)
           enddo
         enddo
       enddo
     enddo
   enddo

   end subroutine calOmegaMatrixCPA
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calChiMatrixCPA()
!  ===================================================================

   use CrystalMatrixModule, only : calChiIntegralCPA
   use CPAMediumModule, only : getSingleSiteTmat
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind) :: L1,L2,L3,L4,K1,K2
   
 !  X = calChiIntegralCPA(local_index, eval, getSingleSiteTmat)
   call calChiIntegral(getSingleSiteTmat)
 ! call writeMatrix('chi', X(:,:,1), dsize*dsize, dsize*dsize)

    do L4 = 1, dsize
      do L3 = 1, dsize
        do L2 = 1, dsize
         do L1 = 1, dsize
           K1 = dsize*(L1 - 1) + L4
           K2 = dsize*(L2 - 1) + L3
           X(K1,K2,1) = X(K1,K2,1) - tau_c(L1,L2)*tau_c(L3,L4)
           X(K1,K2,2) = X(K1,K2,2) - tau_c(L1,L2)*tau_cc(L3,L4)
           X(K1,K2,3) = X(K1,K2,3) - tau_cc(L1,L2)*tau_c(L3,L4)
           X(K1,K2,4) = X(K1,K2,4) - tau_cc(L1,L2)*tau_cc(L3,L4)
         enddo
       enddo
      enddo
    enddo

   end subroutine calChiMatrixCPA
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calChiIntegral(getSingleScatteringMatrix)
!  ===================================================================
   use MPPModule, only : MyPE
   use BZoneModule, only : getNumKs, getAllWeights, getAllKPoints,    &
                           getWeightSum
   use ProcMappingModule, only : isKPointOnMyProc, getNumKsOnMyProc,  &
                                 getKPointIndex, getNumRedundantKsOnMyProc
   use GroupCommModule, only : getGroupID, GlobalSumInGroup, getMyPEinGroup, getNumPEsInGroup
   use IBZRotationModule, only : getNumIBZRotations, getIBZRotationMatrix
   use StrConstModule, only : getStrConstMatrix, &
                   checkFreeElectronPoles, getFreeElectronPoleFactor
   use WriteMatrixModule,  only : writeMatrix
   use MatrixModule, only : computeUAUtc, setupUnitMatrix
   use SystemModule, only : getLmaxKKR
   use CurrentMatrixModule, only : getLindex

   LOGICAL :: isHost

   integer (kind=IntKind) :: k_loc, k, row, col, MyPEinKGroup, method
   integer (kind=IntKind) :: NumKs, kGID, aGID, NumKsOnMyProc, NumRedunKs
   integer (kind=IntKind) :: itertmp, t0size, nt
   integer (kind=IntKind) :: ig, i, j, kkrsz, L1, L2, L3, L4, K1, K2
   integer (kind=IntKind) :: nrot, irot
   integer (kind=IntKind) :: site_config(num_atoms)

   real (kind=RealKind), pointer :: kpts(:,:), weight(:)
   real (kind=RealKind) :: kfac, kaij, aij(3)
   real (kind=RealKind) :: weightSum, kvec(1:3)

   complex (kind=CmplxKind) :: wfac, cfac, efac, fepf, sint, kappa, sfac
   complex (kind=CmplxKind) :: tauprod, tauprod2, tauprod3, tauprod4
   complex (kind=CmplxKind), pointer :: scm(:,:), tauk(:,:), pm(:)
   complex (kind=CmplxKind), pointer :: w0(:,:), w1(:,:), w2(:,:), rotmat(:,:)
   complex (kind=CmplxKind), pointer :: jinvB(:), p_sinej(:)
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
   isHost = .false.
   kappa = sqrt(eval)

   kkrsz = dsize
   nrot = getNumIBZRotations()
   cfac = CONE/real(nrot, RealKind)
   tmbsym = CZERO
   tmb = CZERO
   wint = CZERO
   site_config = 0
   method = 2
 
!  ===================================================================
!  Exchange single site scattering matrix among processes.
!  -------------------------------------------------------------------
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
   w0 => aliasArray2_c(wtmp, kkrsz, kkrsz)
   w1 => aliasArray2_c(wtmpsym,kkrsz,kkrsz)

   do k_loc = 1,NumKsOnMyProc
     k = getKPointIndex(k_loc)
     kvec(1:3) = kpts(1:3,k)*kfac
!    ================================================================
!    get structure constant matrix for the k-point and energy
!    ----------------------------------------------------------------
     call checkFreeElectronPoles(kvec,kappa)
!    ----------------------------------------------------------------
!    ----------------------------------------------------------
     scm => getStrConstMatrix(kvec,kappa,global_index,local_index, &
                            getLmaxKKR(global_index),getLindex(dsize),aij)
!    ----------------------------------------------------------
!    ----------------------------------------------------------------
     fepf = getFreeElectronPoleFactor() ! Returns the free eletron pole factor
!    ----------------------------------------------------------------
     wfac = weight(k)/weightSum
     sfac = SQRTm1*kappa*fepf 
 ! write(6,'(a,i4,a,3f12.5,a,2d12.5,a,2d12.5)')'k-ind = ',k,       &
!    ', kvec = ',kvec,', wfac = ',wfac,', weightsum = ',weightSum
!  
     jinvB => TMP_MatrixBand(1:kkrsz*kkrsz)  ! Use p_MatrixBand as temporary space...
     w2 => aliasArray2_c(jinvB,kkrsz,kkrsz)
     if (global_index == local_index) then
!      ----------------------------------------------------
       call setupUnitMatrix(kkrsz,w2,sfac)
!      ----------------------------------------------------
       w2 = w2 + scm
     else
       w2 = scm
     endif
     p_sinej => tmat_g(:)
!    -----------------------------------------------------------
     call zgemm('n', 'n', kkrsz, kkrsz, kkrsz, CONE, &
          jinvB, kkrsz, p_sinej, kkrsz, CZERO,    &
          WORK, kkrsz)
!    -----------------------------------------------------------
     call zcopy(size(TMP_MatrixBand),WORK,1,TMP_MatrixBand,1)
     TMP_MatrixBand = -TMP_MatrixBand
     do i = 1, kkrsz
        TMP_MatrixBand(i+(i-1)*kkrsz) = &
              fepf+TMP_MatrixBand(i+(i-1)*kkrsz)
     enddo
     if (NumPEsInGroup == 1) then  ! BandSizeCant = KKRMatrixSizeCant
!        ----------------------------------------------------------------
         call ZGETRF(kkrsz, kkrsz, TMP_MatrixBand, kkrsz, IPVT, INFO)
!        ----------------------------------------------------------------
         if (INFO /= 0) then
!          -------------------------------------------------------------
           call ErrorHandler('calChiIntegral','Failed in ZGETRF',INFO)
!          -------------------------------------------------------------
         endif
     else
#ifdef USE_SCALAPACK
!        ----------------------------------------------------------------
         call PZGETRF(kkrsz, kkrsz, TMP_MatrixBand, 1, 1, DESC_A, IPVT, INFO)
!        ----------------------------------------------------------------
         if (INFO /= 0) then
!          -------------------------------------------------------------
           call ErrorHandler('calChiIntegral','Failed in PZGETRF',INFO)
!          -------------------------------------------------------------
         endif
#else
!        ----------------------------------------------------------------
         call ErrorHandler('calChiIntegral','Compiling with -DUSE_SCALAPACK is needed!')
!        ----------------------------------------------------------------
         call ZGETRF(kkrsz, kkrsz, TMP_MatrixBand, kkrsz, IPVT, INFO)
!        ----------------------------------------------------------------
         if (INFO /= 0) then
!          -------------------------------------------------------------
           call ErrorHandler('calCrystalMatrix','Failed in ZGETRF',INFO)
!          -------------------------------------------------------------
         endif
#endif
     endif
!
     if (NumPEsInGroup == 1) then  ! BandSizeCant = KKRMatrixSizeCant
!      ----------------------------------------------------------------
       call ZGETRI(kkrsz, TMP_MatrixBand, kkrsz, IPVT, WORK, LWORK, INFO )
!      ----------------------------------------------------------------
       if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calChiIntegral','Failed in ZGETRS',INFO)
!        -------------------------------------------------------------
       endif
     else
#ifdef USE_SCALAPACK
!      ----------------------------------------------------------------
       call PZGETRI(kkrsz, TMP_MatrixBand, 1, 1, DESC_A, IPVT, &
                      WORK, LWORK, IWORK, LIWORK, INFO)
!      ----------------------------------------------------------------
       if (INFO /= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('calChiIntegral','Failed in PZGETRS',INFO)
!        -------------------------------------------------------------
       endif
#else
!      ----------------------------------------------------------------
       call ErrorHandler('calChiIntegral','Compiling with -DUSE_SCALAPACK is needed!')
!      ----------------------------------------------------------------
       call ZGETRI(kkrsz, TMP_MatrixBand, kkrsz, IPVT, WORK, LWORK, INFO )
!      ----------------------------------------------------------------
       if (INFO /= 0) then
!         -------------------------------------------------------------
          call ErrorHandler('calCrystalMatrix','Failed in ZGETRS',INFO)
!         -------------------------------------------------------------
       endif
#endif
     endif
     TMP_MatrixBand = fepf*TMP_MatrixBand

     pm => TMP_MatrixBand(1:kkrsz*kkrsz)
     tauk => aliasArray2_c(pm,kkrsz,kkrsz)
 !   tmb = matmul(tmat, tauk)
     call zgemm('n', 'n', kkrsz, kkrsz, kkrsz, CONE, tmat, kkrsz, &
         tauk, kkrsz, CZERO, tmb, kkrsz)

     do j = 1, kkrsz
       do i = 1, kkrsz
         tmbsym(i, j) = ((-1.0)**(getLindex(j) - getLindex(i)))*conjg(tmb(j, i))
       enddo
     enddo

     if (k_loc <= NumKsOnMyProc - NumRedunKs .or. MyPEinKGroup == 0) then
       if (nrot > 1) then
         do irot = 1, nrot
           w0 = CZERO
           w1 = CZERO
           WORK = CZERO
           rotmat => getIBZRotationMatrix('c', irot)
!          ----------------------------------------------------------------
           call computeUAUtc(rotmat,kkrsz,kkrsz,rotmat,kkrsz,CONE, &
                        tmb,kkrsz,CONE,w0,kkrsz,WORK)
!          ----------------------------------------------------------------
           call computeUAUtc(rotmat,kkrsz,kkrsz,rotmat,kkrsz,CONE, &
                        tmbsym,kkrsz,CONE,w1,kkrsz,WORK)
!          ----------------------------------------------------------------
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD MAP(TO:kkrsz,nrot,wfac,w0,w1) MAP(FROM:X) PRIVATE(K1,K2) COLLAPSE(4)  
           do L4 = 1, kkrsz
             do L3 = 1, kkrsz
               do L2 = 1, kkrsz
                 do L1 = 1, kkrsz
                   K1 = (L1 - 1)*kkrsz + L4
                   K2 = (L2 - 1)*kkrsz + L3
   !               isHost = omp_is_initial_device()
                   X(K1,K2,1) = X(K1,K2,1) + (1.0/nrot)*wfac*w0(L1,L2)*w0(L3,L4)
                   X(K1,K2,2) = X(K1,K2,2) + (1.0/nrot)*wfac*w0(L1,L2)*w1(L3,L4)
                   X(K1,K2,3) = X(K1,K2,3) + (1.0/nrot)*wfac*w1(L1,L2)*w0(L3,L4)
                   X(K1,K2,4) = X(K1,K2,4) + (1.0/nrot)*wfac*w1(L1,L2)*w1(L3,L4)
                 enddo
               enddo
             enddo
           enddo
        !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
         enddo
       else
         do L4 = 1, kkrsz
           do L3 = 1, kkrsz
             do L2 = 1, kkrsz
               do L1 = 1, kkrsz
                 K1 = (L1 - 1)*kkrsz + L4
                 K2 = (L2 - 1)*kkrsz + L3
                 X(K1,K2,1) = X(K1,K2,1) + wfac*tmb(L1,L2)*tmb(L3,L4)
                 X(K1,K2,2) = X(K1,K2,2) + wfac*tmb(L1,L2)*tmbsym(L3,L4)
                 X(K1,K2,3) = X(K1,K2,3) + wfac*tmbsym(L1,L2)*tmb(L3,L4)
                 X(K1,K2,4) = X(K1,K2,4) + wfac*tmbsym(L1,L2)*tmbsym(L3,L4)
               enddo
             enddo
           enddo
         enddo
       endif
     endif
   enddo
   call GlobalSumInGroup(kGID,X,kkrsz*kkrsz,kkrsz*kkrsz,4)

   end subroutine calChiIntegral
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calVertexCorrectionMatrixCPA()
!  ===================================================================

   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind) :: i

   call calChiMatrixCPA()
   
   if (vertex_corr) then
     call calOmegaMatrixCPA()
   endif

   do i = 1, 4
!    --------------------------------------------------------------------
     call computeAprojB('L', dsize*dsize, X(:,:,i), W(:,:,i), A(:,:,i))
!    --------------------------------------------------------------------
   enddo

   end subroutine calVertexCorrectionMatrixCPA
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCPAConductivity()
!  ===================================================================

   nullify(tau_c, tc)
   deallocate(tau_cc, tcc, tac, temp1, temp2, temp3, temp4, tau1,  &
    tau2, xa, xac, J1avg, J2avg, X, A, W, species_content)
   deallocate(wtmp, wtmpsym, wint, tmb, tmbsym, TMP_MatrixBand, WORK)  
   deallocate(IPVT, IWORK, tmat_g)

   end subroutine endCPAConductivity
!  ===================================================================
end module CPAConductivityModule
