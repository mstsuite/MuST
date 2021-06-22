module ConductivityModule
   use KindParamModule, only : IntKind, QuadRealKind, QuadCmplxKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : PI, ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors

public :: initConductivity,        &
          endConductivity,   &
          calSigmaTildeCPA0, &
          calSigmaTildeCPA1, &
          computeSROConductivity, &
          computeCPAConductivity, &
          calConductivity
!

private
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: scf, mode, master_size
   logical :: vertex_corr

   integer (kind=IntKind), allocatable :: print_instruction(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: kmax_kkr(:)
   integer (kind=IntKind), allocatable :: kmax_phi(:)
   integer (kind=Intkind), allocatable :: lofk(:), mofk(:), jofk(:), m1m(:)
   integer (kind=Intkind), allocatable :: lofj(:), mofj(:)
   integer (kind=IntKind), allocatable :: num_species(:)
   integer (kind=IntKind) :: lmax_phi_max, kmax_kkr_max, lmax_green_max, &
                             lmax_sigma, lmax_sigma_2, lmax_cg
   integer (kind=IntKind) :: kmax_phi_max, kmax_green_max, iend_max, &
                             kmax_sigma, kmax_sigma_2, kmax_max, kmax_cg
   integer (kind=IntKind) :: NumSpecies, NumClusterSize

   complex (kind=CmplxKind), allocatable :: iden(:,:)
   complex (kind=CmplxKind), allocatable :: wmatSRO(:,:,:)
   complex (kind=CmplxKind), allocatable :: vcSRO(:,:,:)
!  Conductivity data stored here
!  --------------------------------------------------------------------
   complex (kind=CmplxKind), allocatable :: sigmatilde(:,:,:), sigmatilde2(:,:,:), &
                                          sigmatilde3(:,:,:), sigmatilde4(:,:,:)

   complex (kind=CmplxKind), allocatable :: sigma(:,:,:)
!  ---------------------------------------------------------------------

   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, kGID
   logical :: Initialized = .false.
! 
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initConductivity(energy, num_atoms, lmaxkkr, lmaxphi, lmaxgreen,  &
                              pola, cant, rel, istop, iprint, vc)
!  ===================================================================
   use RadialGridModule, only : getNumRmesh, getMaxNumRmesh
   use AtomModule, only : getLocalNumSpecies
   use PolyhedraModule, only : getNumPolyhedra
   use WriteMatrixModule, only : writeMatrix
   use ScfDataModule, only : isKKR, isLSMS, isKKRCPA, isKKRCPASRO, isSROSCF
   use CurrentMatrixModule, only : initCurrentMatrixModule
   use SROModule, only : getNeighSize

   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: lmaxkkr(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxphi(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxgreen(num_atoms)
   integer (kind=IntKind), intent(in) :: pola
   integer (kind=IntKind), intent(in) :: cant
   integer (kind=IntKind), intent(in) :: rel
   character (len=*), intent(in) :: istop
   integer (kind=IntKind), intent(in) :: iprint(num_atoms)
   
   complex (kind=CmplxKind), intent(in) :: energy
   logical, intent(in) :: vc

   integer (kind=IntKind) :: i, lmax_max, jmax, iend, kmax, NumSpecies, NumPolyhedra
   integer (kind=IntKind) :: lmax, kl, jl, m, n, l, j, jsize
   
   complex (kind=CmplxKind) :: tmp

   vertex_corr = vc
   LocalNumAtoms = num_atoms
   n_spin_pola = pola
   n_spin_cant = cant
   NumPolyhedra = getNumPolyhedra()
 
   if (isKKR()) then
      mode = 1
   else if (isLSMS()) then
      mode = 2
   else if (isKKRCPA()) then
      mode = 3
   else if (isKKRCPASRO()) then
      mode = 4
      if (isSROSCF() == 1) then
        scf = 1
      else
        scf = 0
      endif
   endif

   if (mode == 1 .or. mode == 2) then
     call ErrorHandler('initConductivity', &
            'Conductivity for KKR/LSMS not yet implemented')
   endif

   allocate( print_instruction(LocalNumAtoms) ) 
   allocate( lmax_kkr(LocalNumAtoms) )
   allocate( lmax_phi(LocalNumAtoms) )
   allocate( kmax_kkr(LocalNumAtoms) )
   allocate( kmax_phi(LocalNumAtoms) )
   allocate( num_species(LocalNumAtoms) )
!
   kmax_kkr_max = 1
   kmax_phi_max = 1
   kmax_green_max = 1
   lmax_phi_max = 0
   lmax_green_max = 0
   lmax_max = 0
   NumSpecies = 0

   do i=1, LocalNumAtoms
      lmax_kkr(i) = lmaxkkr(i)
      kmax_kkr(i) = (lmaxkkr(i)+1)**2
      lmax_phi(i) = lmaxphi(i)
      kmax_phi(i) = (lmaxphi(i)+1)**2
      print_instruction(i) = iprint(i)
      num_species(i) = getLocalNumSpecies(i)
      kmax_kkr_max = max(kmax_kkr_max, (lmaxkkr(i)+1)**2)
      kmax_phi_max = max(kmax_phi_max, (lmaxphi(i)+1)**2)
      kmax_green_max = max(kmax_green_max, (lmaxgreen(i)+1)**2)
      lmax_phi_max = max(lmax_phi_max, lmaxphi(i))
      lmax_green_max = max(lmax_green_max, lmaxgreen(i))
      lmax_max = max(lmax_max, lmaxgreen(i), lmaxkkr(i), lmaxphi(i), &
                     lmax_sigma, lmax_sigma_2)
      NumSpecies = max(NumSpecies, num_species(i))
   enddo
   
   iend_max = getMaxNumRmesh()
   if (mode == 4) then
     NumClusterSize = getNeighSize(1)
     if (vertex_corr) then
       allocate(wmatSRO(NumClusterSize*NumClusterSize*kmax_kkr_max*kmax_kkr_max, &
         NumClusterSize*NumClusterSize*kmax_kkr_max*kmax_kkr_max, 4))
       wmatSRO = CZERO
     endif
     allocate(vcSRO(NumClusterSize*NumClusterSize*kmax_kkr_max*kmax_kkr_max, &
            NumClusterSize*NumClusterSize*kmax_kkr_max*kmax_kkr_max, 4))
     vcSRO = CZERO
   endif

!  -------------------------------------------------------------------
!  Initialize the arrays which store conductivity data
!  -------------------------------------------------------------------
   jsize = master_size

   allocate(lofk(kmax_kkr_max), mofk(kmax_cg),  &
         jofk(kmax_cg), m1m(-lmax_cg:lmax_cg))
   allocate(lofj((lmax_cg+1)*(lmax_cg+2)/2), mofj((lmax_cg+1)*(lmax_cg+2)/2))

   allocate(iden(jsize, jsize))
   allocate(sigmatilde(3, 3, n_spin_pola), sigmatilde2(3, 3, n_spin_pola), sigmatilde3(3, 3, n_spin_pola),  &
    sigmatilde4(3, 3, n_spin_pola), sigma(3, 3, n_spin_pola))

   iden = CZERO
   sigmatilde = CZERO; sigmatilde2 = CZERO
   sigmatilde3 = CZERO; sigmatilde4 = CZERO
   sigma = CZERO

   do i = 1, jsize
      iden(i, i) = CONE
   enddo

!  -------------------------------------------------------------------
   call initCurrentMatrixModule(energy, num_atoms, lmaxkkr, lmaxphi, &
           lmaxgreen, pola, cant, rel, istop, iprint, mode, vc)
!  -------------------------------------------------------------------

   Initialized = .true.

   end subroutine initConductivity
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endConductivity()
!  ===================================================================

   use CurrentMatrixModule, only : endCurrentMatrixModule

   deallocate(sigmatilde, sigmatilde2, sigmatilde3, sigmatilde4)
   deallocate(sigma)

   
   deallocate(print_instruction, lmax_kkr, lmax_phi, kmax_kkr, &
     kmax_phi, iden)  

!  ----------------------------------------------------
   call endCurrentMatrixModule()
!  ----------------------------------------------------

   end subroutine endConductivity
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeCPA0(n, dir1, dir2, is, caltype) result(sigma0)
!  =================================================================== 
!
!  Calculates sigma_tilde0(z_1, z_2) at z_1 = z_2 = e_F + i*delta
!  Can do calculation for both single site CPA and supercell case
!
   use PolyhedraModule, only : getVolume
   use CPAMediumModule, only : getCPAMatrix
   use SROModule, only : getSROMatrix, getSROParam
   use AtomModule, only : getLocalSpeciesContent, getLocalNumSpecies
   use SystemVolumeModule, only : getAtomicVPVolume
   use WriteMatrixModule, only : writeMatrix
   use CurrentMatrixModule, only : getJMatrix

   integer (kind=IntKind), intent(in) :: n, dir1, dir2, is, caltype
   integer (kind=IntKind) :: ic, ic1, dsize, L, num_species
   real (kind=RealKind) :: Omega, c_a, c_b, wab, coeff
   complex (kind=CmplxKind), pointer :: tau_ctemp(:,:)
   complex (kind=CmplxKind), allocatable :: tau_c(:,:), tau_cc(:,:)
   complex (kind=CmplxKind), allocatable :: temp1(:,:), temp2(:,:), & 
                      temp3(:,:), temp4(:,:), taua(:,:), taub(:,:)
   complex (kind=CmplxKind), pointer :: Ja(:,:), Jb(:,:), Jc(:,:)

   complex (kind=CmplxKind) :: sigma0

   tau_ctemp => getCPAMatrix('Tau',site=n,atom=0)
   dsize = kmax_kkr_max

   allocate(temp1(dsize, dsize), temp2(dsize, dsize), &
    temp3(dsize, dsize), temp4(dsize, dsize), taua(dsize, dsize), &
    tau_cc(dsize, dsize), tau_c(dsize, dsize), taub(dsize, dsize))
   temp1 = CZERO; temp2 = CZERO; temp3 = CZERO; temp4 = CZERO
   taua = CZERO; taub = CZERO

   Omega = getAtomicVPVolume(n)
   tau_c = tau_ctemp(1:dsize, 1:dsize)
   tau_cc = conjg(tau_c)
   num_species = getLocalNumSpecies(n)
   sigma0 = CZERO

   do ic = 1, num_species
    do ic1 = 1, num_species
      c_a = getLocalSpeciesContent(n, ic)
      c_b = getLocalSpeciesContent(n, ic1)
      coeff = -(c_a*c_b)/(PI*Omega)
      if (caltype == 1) then
        Ja => getJMatrix(n, ic, is, dir1, 1, 1)
        taua = tau_c; taub = tau_c
      else if (caltype == 2) then
        Ja => getJMatrix(n, ic, is, dir1, 3, 1)
        taua = tau_c; taub = tau_cc
      else if (caltype == 3) then
        Ja => getJMatrix(n, ic, is, dir1, 2, 1)
        taua = tau_cc; taub = tau_c
      else if (caltype == 4) then
        Ja => getJMatrix(n, ic, is, dir1, 4, 1)
        taua = tau_cc; taub = tau_cc
      else
        call ErrorHandler('calSigmaTildeCPA0','Incorrect caltype (1-4)', caltype)
      endif
      Jb => getJMatrix(n, ic, is, dir2, caltype, 0)
      temp4 = Jb
!     call writeMatrix('Jb - Jc', temp4, kmax_kkr_max, kmax_kkr_max)
!     ---------------------------------------------------------------------
!     call zaxpy(dsize*dsize, -CONE, jtspace(:,:,1,ic1,is,dir2),1,temp4,1)
!     ---------------------------------------------------------------------
      call zgemm('N', 'n', dsize, dsize, dsize, CONE, temp4, &
         dsize, taub, dsize, CZERO, temp1, dsize)
!     ---------------------------------------------------------------------
      call zgemm('N', 'n', dsize, dsize, dsize, CONE, taua, dsize,  &
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

   end function calSigmaTildeCPA0
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeCPA1(n, dir1, dir2, is, eval, caltype) result(sigma1)
!  ===================================================================

   use CPAMediumModule, only : getSingleSiteTmat
   use SystemVolumeModule, only : getAtomicVPVolume
   use AtomModule, only : getLocalSpeciesContent, getLocalNumSpecies
   use CurrentMatrixModule, only : getJMatrix
   use CrystalMatrixModule, only : calSigmaIntegralCPA
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, dir1, dir2, is, caltype
   complex (kind=CmplxKind), intent(in) :: eval 

   integer (kind=IntKind) :: ic1, ic2, num_species
   real (kind=RealKind) :: Omega, c_a, c_b, coeff
   complex (kind=CmplxKind) :: sigma1
   complex (kind=CmplxKind), pointer :: J1(:,:), J2(:,:)

   Omega = getAtomicVPVolume(n)
   num_species = getLocalNumSpecies(n)
   sigma1 = CZERO
    
   do ic1 = 1, num_species
     do ic2 = 1, num_species
       c_a = getLocalSpeciesContent(n, ic1)
       c_b = getLocalSpeciesContent(n, ic2)
       coeff = -(c_a*c_b)/(PI*Omega)
       if (caltype == 1 .or. caltype == 4) then
         J1 => getJMatrix(n, ic1, is, dir1, caltype, 1)
         J2 => getJMatrix(n, ic2, is, dir2, caltype, 1)
       else if (caltype == 2) then
         J1 => getJMatrix(n, ic1, is, dir1, 3, 1)
         J2 => getJMatrix(n, ic2, is, dir2, 2, 1)
       else if (caltype == 3) then
         J1 => getJMatrix(n, ic1, is, dir1, 2, 1)
         J2 => getJMatrix(n, ic2, is, dir2, 3, 1)
       else
         call ErrorHandler('calSigmaTildeCPA1', 'Incorrect caltype (1-4)', caltype)
       endif
         sigma1 = sigma1 + coeff*calSigmaIntegralCPA(n, eval, J1, J2, &
           getSingleSiteTmat, tau_needed=.true.,use_tmat=.true.,caltype=caltype)
       nullify(J1, J2)
     enddo
   enddo
       
   end function calSigmaTildeCPA1
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calOmegaMatrixCPA(n, is) result(wmat)
!  ===================================================================

   use SSSolverModule, only : getScatteringMatrix
   use CPAMediumModule, only : getSingleSiteTmat, getCPAMatrix
   use MatrixModule, only : computeAprojB
   use AtomModule, only : getLocalSpeciesContent, getLocalNumSpecies
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, is
   integer (kind=IntKind) :: num_species, ic, L1, L2, L3, L4, K1, K2
   real (kind=RealKind) :: c_a
   complex (kind=CmplxKind), pointer :: ta(:,:), tc(:,:), tauc(:,:)
   complex (kind=CmplxKind) :: tac(kmax_kkr_max, kmax_kkr_max), &
     tcc(kmax_kkr_max, kmax_kkr_max), taucc(kmax_kkr_max, kmax_kkr_max)
   complex (kind=CmplxKind) :: temp(kmax_kkr_max, kmax_kkr_max), &
     tempc(kmax_kkr_max, kmax_kkr_max)
   complex (kind=CmplxKind) :: xa(kmax_kkr_max, kmax_kkr_max), &
     xac(kmax_kkr_max, kmax_kkr_max)
   complex (kind=CmplxKind) :: wmat(kmax_kkr_max*kmax_kkr_max, &
     kmax_kkr_max*kmax_kkr_max,4)
 
   wmat = CZERO
   num_species = getLocalNumSpecies(n)

   tc => getSingleSiteTmat('TInv-Matrix', spin=is, site=n, atom=0)
   tauc => getCPAMatrix('Tau',site=n,atom=0)
   tcc = conjg(tc)
   taucc = conjg(tauc)

   do ic = 1, num_species
     tac = CZERO; temp = CZERO; tempc = CZERO; xa = CZERO; xac = CZERO
     c_a = getLocalSpeciesContent(n, ic)
     ta => getScatteringMatrix('TInv-Matrix', spin=is, site=n, atom=ic)
     tac = conjg(ta)
     temp = ta - tc
     tempc = tac - tcc
!    ------------------------------------------------------------
     call computeAprojB('L', kmax_kkr_max, temp, tauc, xa)
     call computeAprojB('L', kmax_kkr_max, tempc, taucc, xac)
!    ------------------------------------------------------------
     do L4 = 1, kmax_kkr_max
       do L3 = 1, kmax_kkr_max
         do L2 = 1, kmax_kkr_max
           do L1 = 1, kmax_kkr_max
             K1 = kmax_kkr_max*(L1 - 1) + L4
             K2 = kmax_kkr_max*(L2 - 1) + L3
             wmat(K1,K2,1) = wmat(K1,K2,1) - c_a*xa(L1,L2)*xa(L3,L4)
             wmat(K1,K2,2) = wmat(K1,K2,2) - c_a*xa(L1,L2)*xac(L3,L4)
             wmat(K1,K2,3) = wmat(K1,K2,3) - c_a*xac(L1,L2)*xa(L3,L4)
             wmat(K1,K2,4) = wmat(K1,K2,4) - c_a*xac(L1,L2)*xac(L3,L4)
           enddo
         enddo
       enddo
     enddo
   enddo

   end function calOmegaMatrixCPA
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calChiMatrixCPA(n, is, e) result(chi_mat)
!  ===================================================================

   use CrystalMatrixModule, only : calChiIntegralCPA
   use CPAMediumModule, only : getCPAMatrix, getSingleSiteTmat
   
   integer (kind=IntKind), intent(in) :: n, is
   complex (kind=CmplxKind), intent(in) :: e

   integer (kind=IntKind) :: L1,L2,L3,L4,K1,K2
   complex (kind=CmplxKind), pointer :: tau_c(:,:)
   complex (kind=CmplxKind) :: tau_cc(kmax_kkr_max,kmax_kkr_max)
   complex (kind=CmplxKind) :: chi_mat(kmax_kkr_max*kmax_kkr_max, &
                            kmax_kkr_max*kmax_kkr_max, 4)

   tau_cc = CZERO
   tau_c => getCPAMatrix('Tau',site=n,atom=0)
   tau_cc = conjg(tau_c)

   chi_mat = calChiIntegralCPA(n, e, getSingleSiteTmat)

   do L4 = 1, kmax_kkr_max
     do L3 = 1, kmax_kkr_max
       do L2 = 1, kmax_kkr_max
         do L1 = 1, kmax_kkr_max
           K1 = kmax_kkr_max*(L1 - 1) + L4
           K2 = kmax_kkr_max*(L2 - 1) + L3
           chi_mat(K1,K2,1) = chi_mat(K1,K2,1) - tau_c(L1,L2)*tau_c(L3,L4)
           chi_mat(K1,K2,2) = chi_mat(K1,K2,2) - tau_c(L1,L2)*tau_cc(L3,L4)
           chi_mat(K1,K2,3) = chi_mat(K1,K2,3) - tau_cc(L1,L2)*tau_c(L3,L4)
           chi_mat(K1,K2,4) = chi_mat(K1,K2,4) - tau_cc(L1,L2)*tau_cc(L3,L4)
         enddo
       enddo
     enddo
   enddo

   end function calChiMatrixCPA
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calVertexCorrectionMatrixCPA(n, is, e) result(A)
!  ===================================================================

   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, is
   complex (kind=CmplxKind), intent(in) :: e
 
   integer (kind=IntKind) :: i
   complex (kind=CmplxKind) :: X(kmax_kkr_max*kmax_kkr_max, &
     kmax_kkr_max*kmax_kkr_max, 4), W(kmax_kkr_max*kmax_kkr_max, &
     kmax_kkr_max*kmax_kkr_max, 4), A(kmax_kkr_max*kmax_kkr_max, &
     kmax_kkr_max*kmax_kkr_max, 4)

   X = CZERO; W = CZERO; A = CZERO

   X = calChiMatrixCPA(n, is, e)
   if (vertex_corr) then
     W = calOmegaMatrixCPA(n, is)
   endif
   do i = 1, 4
!    --------------------------------------------------------------------
     call computeAprojB('L', kmax_kkr_max*kmax_kkr_max, X(:,:,i), W(:,:,i), A(:,:,i))
!    --------------------------------------------------------------------
   enddo

   end function calVertexCorrectionMatrixCPA
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeCPA1VC(n, dir1, dir2, is, caltype, chi) result(sigma1)
!  ===================================================================

   use CPAMediumModule, only : getSingleSiteTmat
   use SystemVolumeModule, only : getAtomicVPVolume
   use AtomModule, only : getLocalSpeciesContent, getLocalNumSpecies
   use CurrentMatrixModule, only : getJMatrix1D

   integer (kind=IntKind), intent(in) :: n, dir1, dir2, is, caltype
   complex (kind=CmplxKind), intent(in) :: chi(kmax_kkr_max*kmax_kkr_max, &
                                         kmax_kkr_max*kmax_kkr_max)

   integer (kind=IntKind) :: num_species, ic1, ic2, K, L1, L4, K1, K1_t
   real (kind=RealKind) :: Omega, c_a, c_b, coeff
   complex (kind=CmplxKind) :: sigma1 ! dir1, dir2, caltype
   complex (kind=CmplxKind), pointer :: J1(:), J2(:)
   complex (kind=CmplxKind) :: J1avg(kmax_kkr_max*kmax_kkr_max), &
            J2avg(kmax_kkr_max*kmax_kkr_max)

   J1avg = CZERO; J2avg = CZERO; sigma1 = CZERO
   Omega = getAtomicVPVolume(n)
   num_species = getLocalNumSpecies(n)
   coeff = -1.0/(PI*Omega)

   do ic1 = 1, num_species
     c_a = getLocalSpeciesContent(n, ic1)
     if (caltype == 2) then
       J1 => getJMatrix1D(n, ic1, is, dir1, 3)
     else if (caltype == 3) then
       J1 => getJMatrix1D(n, ic1, is, dir1, 2)
     else
       J1 => getJMatrix1D(n, ic1, is, dir1, caltype)
     endif
     J2 => getJMatrix1D(n, ic1, is, dir2, caltype)
     J1avg = J1avg + c_a*J1
     J2avg = J2avg + c_a*J2
    !call zaxpy(kmax_kkr_max*kmax_kkr_max, c_a, J1, 1, J1avg, 1)
    !call zaxpy(kmax_kkr_max*kmax_kkr_max, c_a, J2, 1, J2avg, 1)
   enddo

   do K = 1, kmax_kkr_max**2
     do L1 = 1, kmax_kkr_max
       do L4 = 1, kmax_kkr_max
         K1 = kmax_kkr_max*(L1 - 1) + L4
         K1_t = kmax_kkr_max*(L4 - 1) + L1
         sigma1 = sigma1 + coeff*J1avg(K1_t)*chi(K1, K)*J2avg(K)
       enddo
     enddo
   enddo


   end function calSigmaTildeCPA1VC
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeCPAConductivity(n, is, delta, pot_type, dirnum, e)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: n, is, pot_type, dirnum
   real (kind=RealKind), intent(in) :: delta
   complex (kind=CmplxKind), intent(in) :: e

   integer (kind=IntKind) :: etype
   integer (kind=IntKind) :: dir, dir1
   complex (kind=CmplxKind) :: int_val(4)
   complex (kind=CmplxKind) :: A(kmax_kkr_max*kmax_kkr_max, &
                        kmax_kkr_max*kmax_kkr_max, 4)

   A = calVertexCorrectionMatrixCPA(n, is, e)

   do dir = 1, dirnum
     do dir1 = 1, dirnum
       int_val = CZERO
       do etype = 1, 4
         int_val(etype) = &
           calSigmaTildeCPA1VC(n, dir, dir1, is, etype, A(:,:,etype)) + &
           calSigmaTildeCPA0(n, dir, dir1, is, etype)
       enddo
       sigmatilde(dir,dir1,is) = int_val(1)
       sigmatilde2(dir,dir1,is) = int_val(2)
       sigmatilde3(dir,dir1,is) = int_val(3)
       sigmatilde4(dir,dir1,is) = int_val(4)
       if (n_spin_pola == 1) then
         sigmatilde(dir,dir1,is) = 2*sigmatilde(dir,dir1,is)
         sigmatilde2(dir,dir1,is) = 2*sigmatilde2(dir,dir1,is)
         sigmatilde3(dir,dir1,is) = 2*sigmatilde3(dir,dir1,is) 
         sigmatilde4(dir,dir1,is) = 2*sigmatilde4(dir,dir1,is)
       endif
       sigma(dir,dir1,is) = 0.25*(sigmatilde(dir,dir1,is) - &
       sigmatilde2(dir,dir1,is) - sigmatilde3(dir,dir1,is) + &
       sigmatilde4(dir,dir1,is))
     enddo
   enddo

   end subroutine computeCPAConductivity
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calOmegaMatrixSRO(n, is)
!  ===================================================================

   use SROModule, only : getSROMatrix
   use AtomModule, only : getLocalSpeciesContent, getLocalNumSpecies
   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, is

   integer (kind=IntKind) :: num_species, ic
   integer (kind=IntKind) :: l,m,z,r,L1,L2,L3,L4,K1,K2,K3,K4,CK1,CK2
   real (kind=RealKind) :: c_a
   complex (kind=CmplxKind), pointer :: tauc(:,:), Ta(:,:), Tc(:,:)
   complex (kind=CmplxKind) :: taucc(kmax_kkr_max*NumClusterSize, &
     kmax_kkr_max*NumClusterSize), Tcc(kmax_kkr_max*NumClusterSize, &
     kmax_kkr_max*NumClusterSize), Tac(kmax_kkr_max*NumClusterSize, &
     kmax_kkr_max*NumClusterSize), xa(kmax_kkr_max*NumClusterSize, &
     kmax_kkr_max*NumClusterSize), xac(kmax_kkr_max*NumClusterSize, &
     kmax_kkr_max*NumClusterSize), tdiff(kmax_kkr_max*NumClusterSize, &
     kmax_kkr_max*NumClusterSize), tdiffc(kmax_kkr_max*NumClusterSize, &
     kmax_kkr_max*NumClusterSize)
   taucc = CZERO; Tcc = CZERO; 

   num_species = getLocalNumSpecies(n)

   tauc => getSROMatrix('blk-tau', n=n, ic=0, is=is)
   Tc => getSROMatrix('blk-tinv', n=n, ic=0, is=is)
   Tcc = conjg(Tc)

   do L4 = 1, kmax_kkr_max
     do L1 = 1, kmax_kkr_max
       do l = 1, NumClusterSize
         do m = 1, NumClusterSize
           Lp = kmax_kkr_max*(l - 1) + L1
           Lpp = kmax_kkr_max*(m - 1) + L4
           taucc(Lp, Lpp) = (-1.0)**(lofk(L1) - lofk(L4)) * conjg(tauc(Lpp, Lp))
         enddo
       enddo
     enddo
   enddo

   do ic = 1, num_species
     Tac = CZERO; xa = CZERO; xac = CZERO; tdiff = CZERO; tdiffc = CZERO
     c_a = getLocalSpeciesContent(n, ic)
     Ta => getSROMatrix('blk-tinv', n=n, ic=ic, is=is)
     Tac = conjg(Ta)
     tdiff = Ta - Tc
     tdiffc = Tac - Tcc
!    -------------------------------------------------------------------------
     call computeAprojB('L', kmax_kkr_max*NumClusterSize, tdiff, tauc, xa)
     call computeAprojB('L', kmax_kkr_max*NumClusterSize, tdiffc, taucc, xac)
!    -------------------------------------------------------------------------
     do z = 1, NumClusterSize
       do L2 = 1, kmax_kkr_max
         K2 = kmax_kkr_max*(z - 1) + L2
         do r = 1, NumClusterSize
           do L3 = 1, kmax_kkr_max
             K3 = kmax_kkr_max*(r - 1) + L3
             do l = 1, NumClusterSize
               do L1 = 1, kmax_kkr_max
                 K1 = kmax_kkr_max*(l - 1) + L1
                 do m = 1, NumClusterSize
                   do L4 = 1, kmax_kkr_max
                     K4 = kmax_kkr_max*(m - 1) + L4
                     CK1 = (kmax_kkr_max**2)*(NumClusterSize*(l-1)+m-1)+kmax_kkr_max*(L1-1)+L4
                     CK2 = (kmax_kkr_max**2)*(NumClusterSize*(z-1)+r-1)+kmax_kkr_max*(L2-1)+L3
                     wmatSRO(CK1,CK2,1) = wmatSRO(CK1,CK2,1) - c_a*xa(K1,K2)*xa(K3,K4)
                     wmatSRO(CK1,CK2,2) = wmatSRO(CK1,CK2,2) - c_a*xa(K1,K2)*xac(K3,K4)
                     wmatSRO(CK1,CK2,3) = wmatSRO(CK1,CK2,3) - c_a*xac(K1,K2)*xa(K3,K4)
                     wmatSRO(CK1,CK2,4) = wmatSRO(CK1,CK2,4) - c_a*xac(K1,K2)*xac(K3,K4)
                   enddo
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo

   end subroutine calOmegaMatrixSRO
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calChiMatrixSRO(n, is, e) result(chi)
!  ===================================================================   

   use CrystalMatrixModule, only : calChiIntegralSRO, getChiIntegralSRO
   use CPAMediumModule, only : getSingleSiteTmat
   use SROModule, only : getSROMatrix

   integer (kind=IntKind), intent(in) :: n, is
   complex (kind=CmplxKind), intent(in) :: e
   
   integer (kind=IntKind) :: Lp, Lpp, l, z, r, m, L1, L2, L3, L4, CK1, CK2
   integer (kind=IntKind) :: K1, K2, K3, K4
   complex (kind=CmplxKind), pointer :: tauc(:,:), chi(:,:,:)
   complex (kind=CmplxKind) :: taucc(kmax_kkr_max*NumClusterSize, &
       kmax_kkr_max*NumClusterSize)

   taucc = CZERO
   tauc => getSROMatrix('blk-tau', n=n, ic=0, is=is)
   
   do L4 = 1, kmax_kkr_max
     do L1 = 1, kmax_kkr_max
       do l = 1, NumClusterSize
         do m = 1, NumClusterSize
           Lp = kmax_kkr_max*(l - 1) + L1 
           Lpp = kmax_kkr_max*(m - 1) + L4 
           taucc(Lp, Lpp) = (-1.0)**(lofk(L1) - lofk(L4)) * conjg(tauc(Lpp, Lp))
         enddo
       enddo
     enddo
   enddo
 
!  -----------------------------------------------------------
   call calChiIntegralSRO(n, e, getSingleSiteTmat)
   chi => getChiIntegralSRO()
!  -----------------------------------------------------------

   do z = 1, NumClusterSize
     do L2 = 1, kmax_kkr_max
       K2 = kmax_kkr_max*(z - 1) + L2
       do r = 1, NumClusterSize
         do L3 = 1, kmax_kkr_max
           K3 = kmax_kkr_max*(r - 1) + L3
           do l = 1, NumClusterSize
             do L1 = 1, kmax_kkr_max
               K1 = kmax_kkr_max*(l - 1) + L1
               do m = 1, NumClusterSize
                 do L4 = 1, kmax_kkr_max
                   K4 = kmax_kkr_max*(m - 1) + L4
                   CK1 = (kmax_kkr_max**2)*(NumClusterSize*(l-1)+m-1)+kmax_kkr_max*(L1-1)+L4
                   CK2 = (kmax_kkr_max**2)*(NumClusterSize*(z-1)+r-1)+kmax_kkr_max*(L2-1)+L3
                   chi(CK1,CK2,1) = chi(CK1,CK2,1) - tauc(K1,K2)*tauc(K3,K4)
                   chi(CK1,Ck2,2) = chi(CK1,CK2,2) - tauc(K1,K2)*taucc(K3,K4)
                   chi(CK1,CK2,3) = chi(CK1,CK2,3) - taucc(K1,K2)*tauc(K3,K4)
                   chi(CK1,CK2,4) = chi(CK1,CK2,4) - taucc(K1,K2)*taucc(K3,K4)
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
     enddo
   enddo

   end function calChiMatrixSRO
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calVertexCorrectionMatrixSRO(n,is,e)
!  ===================================================================

   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, is
   complex (kind=CmplxKind), intent(in) :: e

   integer (kind=IntKind) :: i
   complex (kind=CmplxKind), pointer :: chi(:,:,:)

   chi => calChiMatrixSRO(n, is, e)
   if (vertex_corr) then
     call calOmegaMatrixSRO(n,is)
     do i = 1, 4
!      -----------------------------------------------------------------
       call computeAprojB('L', kmax_kkr_max*kmax_kkr_max*NumClusterSize*NumClusterSize, &
         chi(:,:,i), wmatSRO(:,:,i), vcSRO(:,:,i))
!      -----------------------------------------------------------------
     enddo
   else 
     vcSRO = chi
   endif

   end subroutine calVertexCorrectionMatrixSRO
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeSRO1(n, dir1, dir2, is, caltype) result(sigma1)
!  ===================================================================

   use CurrentMatrixModule, only : getJMatrix1D
   use AtomModule, only : getLocalSpeciesContent, getLocalNumSpecies
   use SystemVolumeModule, only : getAtomicVPVolume

   integer (kind=IntKind), intent(in) :: n, dir1, dir2, is, caltype

   integer (kind=IntKind) :: num_species, ic1, ic2, CK2, CK1, CK1_t, L1,L4,l,m
   real (kind=RealKind) :: Omega, c_a, c_b, coeff
   complex (kind=CmplxKind) :: sigma1 ! dir1, dir2, caltype
   complex (kind=CmplxKind), pointer :: J1(:), J2(:)
   complex (kind=CmplxKind) :: J1avg(kmax_kkr_max*kmax_kkr_max*NumClusterSize*NumClusterSize), &
            J2avg(kmax_kkr_max*kmax_kkr_max*NumClusterSize*NumClusterSize)

   J1avg = CZERO; J2avg = CZERO; sigma1 = CZERO
   Omega = getAtomicVPVolume(n)
   num_species = getLocalNumSpecies(n)
   coeff = -1.0/(PI*Omega)

   do ic1 = 1, num_species
     c_a = getLocalSpeciesContent(n, ic1)
     if (caltype == 2) then
       J1 => getJMatrix1D(n, ic1, is, dir1, 3)
     else if (caltype == 3) then
       J1 => getJMatrix1D(n, ic1, is, dir1, 2)
     else
       J1 => getJMatrix1D(n, ic1, is, dir1, caltype)
     endif
     J2 => getJMatrix1D(n, ic1, is, dir2, caltype)
     J1avg = J1avg + c_a*J1
     J2avg = J2avg + c_a*J2
   enddo

   do CK2 = 1, NumClusterSize*NumClusterSize*kmax_kkr_max*kmax_kkr_max
     do l = 1, NumClusterSize
       do L1 = 1, kmax_kkr_max
         do m = 1, NumClusterSize
           do L4 = 1, kmax_kkr_max
             CK1 = (kmax_kkr_max**2)*(NumClusterSize*(l-1)+m-1)+kmax_kkr_max*(L1-1)+L4
             CK1_t = (kmax_kkr_max**2)*(NumClusterSize*(m-1)+l-1)+kmax_kkr_max*(L4-1)+L1
             sigma1 = sigma1 + coeff*J1avg(CK1_t)*vcSRO(CK1,CK2,caltype)*J2avg(CK2)
           enddo
         enddo
       enddo
     enddo
   enddo

   end function calSigmaTildeSRO1
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeSROConductivity(n, is, delta, pot_type, dirnum, e)
!  ===================================================================   

   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, is, pot_type, dirnum
   real (kind=RealKind), intent(in) :: delta
   complex (kind=CmplxKind), intent(in) :: e

   integer (kind=IntKind) :: etype
   integer (kind=IntKind) :: dir, dir1
   complex (kind=CmplxKind) :: int_val(4)
   complex (kind=CmplxKind), pointer :: chi(:,:,:)

   call calVertexCorrectionMatrixSRO(n,is,e)

   do dir = 1, dirnum
     do dir1 = 1, dirnum
       int_val = CZERO
       do etype = 1, 4
         int_val(etype) = calSigmaTildeSRO1(n, dir, dir1, is, etype)
       enddo
       sigmatilde(dir,dir1,is) = int_val(1)
       sigmatilde2(dir,dir1,is) = int_val(2)
       sigmatilde3(dir,dir1,is) = int_val(3)
       sigmatilde4(dir,dir1,is) = int_val(4)
       if (n_spin_pola == 1) then
         sigmatilde(dir,dir1,is) = 2*sigmatilde(dir,dir1,is)
         sigmatilde2(dir,dir1,is) = 2*sigmatilde2(dir,dir1,is)
         sigmatilde3(dir,dir1,is) = 2*sigmatilde3(dir,dir1,is)
         sigmatilde4(dir,dir1,is) = 2*sigmatilde4(dir,dir1,is)
       endif
       sigma(dir,dir1,is) = 0.25*(sigmatilde(dir,dir1,is) - & 
       sigmatilde2(dir,dir1,is) - sigmatilde3(dir,dir1,is) + & 
       sigmatilde4(dir,dir1,is))
     enddo
   enddo

   end subroutine computeSROConductivity
!  ===================================================================
   
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calConductivity(LocalNumAtoms, n_spin_pola)
!  ===================================================================

   use SSSolverModule, only : solveSingleScattering
   use ScfDataModule, only : getFermiEnergyImagPart, isFermiEnergyRealPart, &
     useStepFunctionForSigma, useCubicSymmetryForSigma, getFermiEnergyRealPart
   use ValenceDensityModule, only : getFermiEnergy
   use WriteMatrixModule, only : writeMatrix
   use CrystalMatrixModule, only : calCrystalMatrix, retrieveTauSRO
   use CPAMediumModule, only : computeCPAMedium, populateBigTCPA, getSingleSiteMatrix
   use SROModule, only : calSpeciesTauMatrix
   use CurrentMatrixModule, only : calCurrentMatrix
   use MPPModule, only : getMyPE

   integer (kind=IntKind), intent(in) :: LocalNumAtoms, n_spin_pola
   integer (kind=IntKind) :: id, is, pot_type, dirnum, MyPE
   real (kind=RealKind) :: delta, efermi, ti, tf
   complex (kind=CmplxKind) :: eval

   call cpu_time(ti)
   delta = getFermiEnergyImagPart()
   pot_type = useStepFunctionForSigma()

   if (useCubicSymmetryForSigma()) then
     dirnum = 1
   else
     dirnum = 3
   endif

   efermi = getFermiEnergy()
   if (isFermiEnergyRealPart()) then
     eval = getFermiEnergyRealPart() + SQRTm1*delta
   else
     eval = efermi + SQRTm1*delta
   endif

   do id = 1, LocalNumAtoms
     do is = 1, n_spin_pola
       if (mode == 3) then
      !  ---------------------------------------------------------------
         call solveSingleScattering(spin=is,site=id,e=eval,vshift=CZERO)
      !  ----------------------------------------------------------------   
         call computeCPAMedium(eval)
      !  --------------------------------------------------------------
         call calCurrentMatrix(id,is,eval,pot_type,mode)
      !  ---------------------------------------------------------------- 
         call computeCPAConductivity(id, is, delta, pot_type, dirnum, eval)
      !  --------------------------------------------------------------
       else if (mode == 4) then
      !  ---------------------------------------------------------------
         call solveSingleScattering(spin=is,site=id,e=eval,vshift=CZERO)
      !  ----------------------------------------------------------------   
      !  Need to investigate properly for multiple sublattices   
         if (scf == 0) then
!          --------------------------------------------------------------
           call computeCPAMedium(eval, do_sro=.true.)
           call populateBigTCPA()
           call calCrystalMatrix(eval, getSingleSiteMatrix ,use_tmat=.true., &
                tau_needed=.true., use_sro=.true.)
           call retrieveTauSRO()
           call calSpeciesTauMatrix()
!          --------------------------------------------------------------
         else if (scf == 1) then
!          --------------------------------------------------------------
           call computeCPAMedium(eval, do_sro = .true.)
!          --------------------------------------------------------------
         endif
!        ----------------------------------------------------------------
         call calCurrentMatrix(id,is,eval,pot_type,mode)
!        ---------------------------------------------------------------- 
         call computeSROConductivity(id, is, delta, pot_type, dirnum, eval)
!        -------------------------------------------------------------- 
       endif
     enddo
   enddo

!  Print calculated results to output file
!  --------------------------------------------------------------------
   if (MyPE == 0) then
     call writeMatrix('sigma', sigma, dirnum, dirnum, n_spin_pola)
     call writeMatrix('sigmatilde', sigmatilde, dirnum, dirnum, n_spin_pola)
     call writeMatrix('sigmatilde2', sigmatilde2, dirnum, dirnum, n_spin_pola)
     call writeMatrix('sigmatilde3', sigmatilde3, dirnum, dirnum, n_spin_pola)
     call writeMatrix('sigmatilde4', sigmatilde4, dirnum, dirnum, n_spin_pola)
!  --------------------------------------------------------------------
     call cpu_time(tf)
     write(*,*) "Time: ", tf-ti, " seconds"
   endif

   end subroutine calConductivity
!  ===================================================================
end module ConductivityModule
