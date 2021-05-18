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
          calSigmaTildeSRO00, &
          computeSROConductivity, &
          computeCPAConductivity, &
          calConductivity
!

private
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: scf, mode, master_size

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
   integer (kind=IntKind) :: NumSpecies

   complex (kind=CmplxKind), allocatable :: iden(:,:)

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
                              pola, cant, rel, istop, iprint)
!  ===================================================================
   use RadialGridModule, only : getNumRmesh, getMaxNumRmesh
   use AtomModule, only : getLocalNumSpecies
   use PolyhedraModule, only : getNumPolyhedra
   use WriteMatrixModule, only : writeMatrix
   use ScfDataModule, only : isKKR, isLSMS, isKKRCPA, isKKRCPASRO, isSROSCF
   use CurrentMatrixModule, only : initCurrentMatrixModule

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

   integer (kind=IntKind) :: i, lmax_max, jmax, iend, kmax, NumSpecies, NumPolyhedra
   integer (kind=IntKind) :: lmax, kl, jl, m, n, l, j, jsize
   
   complex (kind=CmplxKind) :: tmp

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
           lmaxgreen, pola, cant, rel, istop, iprint, mode)
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
      Jc => getJMatrix(n, ic1, is, dir2, caltype, 1)
      temp4 = Jb - Jc
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
   function calSigmaTildeSRO00(n, dir1, dir2, is, caltype) result(sigma00)
!  ===================================================================

   use SROModule, only : getSROMatrix
   use CurrentMatrixModule, only : getJMatrix
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
   use SystemVolumeModule, only : getAtomicVPVolume
 
   integer (kind=IntKind), intent(in) :: n, dir1, dir2, is, caltype
   integer (kind=IntKind) :: ic1, ic2, L
   real (kind=RealKind) :: Omega, c_a, c_b, coeff
   complex (kind=CmplxKind) :: sigma00
   complex (kind=CmplxKind), pointer :: taua(:,:), J1(:,:), J2(:,:), J3(:,:)
   complex (kind=CmplxKind), allocatable :: tauac(:,:), Jdiff(:,:), tau1(:,:), tau2(:,:)
   complex (kind=CmplxKind), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)

   allocate(tauac(kmax_kkr_max, kmax_kkr_max), Jdiff(kmax_kkr_max, kmax_kkr_max), &
      tau1(kmax_kkr_max, kmax_kkr_max), tau2(kmax_kkr_max, kmax_kkr_max), &
      tmp1(kmax_kkr_max, kmax_kkr_max), tmp2(kmax_kkr_max, kmax_kkr_max), &
      tmp3(kmax_kkr_max, kmax_kkr_max))
   tauac = CZERO
   sigma00 = CZERO
   
   Omega = getAtomicVPVolume(n)

   do ic1 = 1, getLocalNumSpecies(n)
     do ic2 = 1, getLocalNumSpecies(n)
       c_a = getLocalSpeciesContent(n, ic1)
       c_b = getLocalSpeciesContent(n, ic2)
       coeff = -(c_a*c_b)/(PI*Omega)
       tauac = CZERO; Jdiff = CZERO; tau1 = CZERO; tau2 = CZERO;
       taua => getSROMatrix('tau11', n, ic1, is)
       tauac = conjg(taua)
       if (caltype == 1) then
         J1 => getJMatrix(n, ic1, is, dir1, caltype ,0)
         tau1 = taua; tau2 = taua;
       else if (caltype == 2) then
         J1 => getJMatrix(n, ic1, is, dir1, 3, 0)
         tau1 = taua; tau2 = tauac;
       else if (caltype == 3) then
         J1 => getJMatrix(n, ic1, is, dir1, 2, 0)
         tau1 = tauac; tau2 = taua;
       else if (caltype == 4) then
         J1 => getJMatrix(n, ic1, is, dir1, caltype, 0)
         tau1 = tauac; tau2 = tauac;
       else
         call ErrorHandler('calSigmaTildeSRO00', 'Incorrect caltype (1-4)', caltype)
       endif
       J2 => getJMatrix(n, ic1, is, dir2, caltype, 0)
       J3 => getJMatrix(n, ic2, is, dir2, caltype, 1)
       Jdiff = J2 - J3
!      ---------------------------------------------------------------
       call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
       CONE, Jdiff, kmax_kkr_max, tau2, kmax_kkr_max, CZERO, tmp1, kmax_kkr_max)
!      ---------------------------------------------------------------
       call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
       CONE, tau1, kmax_kkr_max, tmp1, kmax_kkr_max, CZERO, tmp2, kmax_kkr_max)
!      ---------------------------------------------------------------
       call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
       CONE, J1, kmax_kkr_max, tmp2, kmax_kkr_max, CZERO, tmp3, kmax_kkr_max)
!      ---------------------------------------------------------------
       do L = 1, kmax_kkr_max
         sigma00 = sigma00 + coeff*tmp3(L,L)
       enddo
     enddo
   enddo

   end function calSigmaTildeSRO00
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeSRO010(n, dir1, dir2, is, caltype) result(sigma010)
!  ===================================================================

   use SROModule, only : getDoubleSpeciesTauMatrix, getSROParam, &
                         getNeighSize
   use CurrentMatrixModule, only : getJMatrix
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
   use SystemVolumeModule, only : getAtomicVPVolume
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, dir1, dir2, is, caltype

   integer (kind=IntKind) :: ic1, ic2, L, i, neigh_size
   real (kind=RealKind) :: Omega, c_a, w_ab, coeff
   complex (kind=CmplxKind) :: sigma010
   complex (kind=CmplxKind), pointer :: tauab(:,:), tauabc(:,:), J1(:,:), J2(:,:)
   complex (kind=CmplxKind), allocatable :: tauc(:,:), tau1(:,:), tau2(:,:)
   complex (kind=CmplxKind), allocatable :: taucc(:,:), tmp1(:,:), tmp2(:,:), tmp3(:,:)

   neigh_size = getNeighSize(n)
   Omega = getAtomicVPVolume(n)
   
   allocate(tauc(neigh_size*kmax_kkr_max, neigh_size*kmax_kkr_max), &
      taucc(neigh_size*kmax_kkr_max, neigh_size*kmax_kkr_max), &
      tau1(kmax_kkr_max, kmax_kkr_max), tau2(kmax_kkr_max, kmax_kkr_max), &
      tmp1(kmax_kkr_max, kmax_kkr_max), tmp2(kmax_kkr_max, kmax_kkr_max), &
      tmp3(kmax_kkr_max, kmax_kkr_max))
   sigma010 = CZERO
   
   do ic1 = 1, getLocalNumSpecies(n)
     do ic2 = 1, getLocalNumSpecies(n)
       c_a = getLocalSpeciesContent(n, ic1)
       w_ab = getSROParam(n, ic1, ic2)
       coeff = -(c_a*w_ab)/(PI*Omega)
       tauab => getDoubleSpeciesTauMatrix(n, is, ic1, ic2, 0)
       tauabc => getDoubleSpeciesTauMatrix(n, is, ic1, ic2, 1)
       J2 => getJMatrix(n, ic2, is, dir2, caltype, 0)
      !call writeMatrix('tauab', tauab, kmax_kkr_max*neigh_size, kmax_kkr_max*neigh_size)
       do i = 2, neigh_size
         tmp1 = CZERO; tmp2 = CZERO; tmp3 = CZERO
         tau1 = CZERO; tau2 = CZERO
         if (caltype == 1) then
           tau1 = tauab(1:kmax_kkr_max, &
                (i-1)*kmax_kkr_max + 1:i*kmax_kkr_max); 
           tau2 = tauabc((i-1)*kmax_kkr_max + 1:i*kmax_kkr_max, &
                 1:kmax_kkr_max)
           J1 => getJMatrix(n, ic1, is, dir1, caltype, 0)
       !   call writeMatrix('tau1', tau1, kmax_kkr_max, kmax_kkr_max)
       !   call writeMatrix('tau2', tau2, kmax_kkr_max, kmax_kkr_max)
       !   call ErrorHandler('calSigmaTildeSRO010', 'stop')
         else if (caltype == 2) then
           tau1 = tauab(1:kmax_kkr_max, &
                (i-1)*kmax_kkr_max + 1:i*kmax_kkr_max)
           tau2 = tauabc((i-1)*kmax_kkr_max + 1:i*kmax_kkr_max, &
                 1:kmax_kkr_max)
           J1 => getJMatrix(n, ic1, is, dir1, 3, 0)
         else if (caltype == 3) then
           tau1 = tauabc(1:kmax_kkr_max, &
                (i-1)*kmax_kkr_max+1:i*kmax_kkr_max) 
           tau2 = tauab((i-1)*kmax_kkr_max + 1:i*kmax_kkr_max, &
                  1:kmax_kkr_max)
           J1 => getJMatrix(n, ic1, is, dir1, 2, 0)
         else if (caltype == 4) then
           tau1 = tauabc(1:kmax_kkr_max, &
                (i-1)*kmax_kkr_max+1:i*kmax_kkr_max)
           tau2 = tauabc((i-1)*kmax_kkr_max+1:i*kmax_kkr_max, &
                1:kmax_kkr_max)
           J1 => getJMatrix(n, ic1, is, dir1, caltype, 0)
         endif
!        ---------------------------------------------------------------
         call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
         CONE, J2, kmax_kkr_max, tau2, kmax_kkr_max, CZERO, tmp1, kmax_kkr_max)
!        ---------------------------------------------------------------
         call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
         CONE, tau1, kmax_kkr_max, tmp1, kmax_kkr_max, CZERO, tmp2, kmax_kkr_max)
!        ---------------------------------------------------------------
         call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
         CONE, J1, kmax_kkr_max, tmp2, kmax_kkr_max, CZERO, tmp3, kmax_kkr_max)
!        ---------------------------------------------------------------
         do L = 1, kmax_kkr_max
           sigma010 = sigma010 + coeff*tmp3(L,L)
         enddo
       enddo
     enddo
   enddo

   end function calSigmaTildeSRO010
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeSRO011(n, dir1, dir2, is, caltype) result(sigma011)
!  ===================================================================

   use SROModule, only : getSROMatrix, getNeighSize
   use CurrentMatrixModule, only : getJMatrix
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
   use SystemVolumeModule, only : getAtomicVPVolume

   integer (kind=IntKind), intent(in) :: n, dir1, dir2, is, caltype
   integer (kind=IntKind) :: neigh_size, ic1, ic2, L, i
   real (kind=RealKind) :: Omega, c_a, c_b, coeff
   complex (kind=CmplxKind) :: sigma011
   complex (kind=CmplxKind), pointer :: tauac(:,:), taua(:,:), J1(:,:), J2(:,:)
   complex (kind=CmplxKind), allocatable :: tau1(:,:), tau2(:,:)
   complex (kind=CmplxKind), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)
   
   neigh_size = getNeighSize(n)
   Omega = getAtomicVPVolume(n)

   allocate(tau1(kmax_kkr_max, kmax_kkr_max), tau2(kmax_kkr_max, kmax_kkr_max), &
      tmp1(kmax_kkr_max, kmax_kkr_max), tmp2(kmax_kkr_max, kmax_kkr_max), &
      tmp3(kmax_kkr_max, kmax_kkr_max))
   sigma011 = CZERO

   do ic1 = 1, getLocalNumSpecies(n)
     do ic2 = 1, getLocalNumSpecies(n)
       c_a = getLocalSpeciesContent(n, ic1)
       c_b = getLocalSpeciesContent(n, ic2)
       coeff = (c_a*c_b)/(PI*Omega)
       taua => getSROMatrix('blk-tau', n,ic1, is)
       tauac => getSROMatrix('neg-blk-tau',n,ic1,is)
       J2 => getJMatrix(n, ic2, is, dir2, caltype, 1)
       do i = 2, neigh_size
         if (caltype == 1) then
           J1 => getJMatrix(n, ic1, is, dir1, caltype, 0)
           tau1 = taua(1:kmax_kkr_max, &
             (i - 1)*kmax_kkr_max +1: i*kmax_kkr_max)
           tau2 = taua((i-1)*kmax_kkr_max+1:i*kmax_kkr_max, &
              1:kmax_kkr_max)
         else if (caltype == 2) then
           J1 => getJMatrix(n, ic1, is, dir1, 3, 0)
           tau1 = taua(1:kmax_kkr_max, &
             (i - 1)*kmax_kkr_max + 1:i*kmax_kkr_max)
           tau2 = tauac((i-1)*kmax_kkr_max+1:i*kmax_kkr_max, &
             1:kmax_kkr_max)
         else if (caltype == 3) then
           J1 => getJMatrix(n, ic1, is, dir1, 2, 0)
           tau1 = tauac(1:kmax_kkr_max, &
              (i - 1)*kmax_kkr_max + 1:i*kmax_kkr_max)
           tau2 = taua((i-1)*kmax_kkr_max+1:i*kmax_kkr_max, &
              1:kmax_kkr_max)
         else if (caltype == 4) then
           J1 => getJMatrix(n, ic1, is, dir1, caltype, 0)
           tau1 = tauac(1:kmax_kkr_max, &
              (i - 1)*kmax_kkr_max + 1:i*kmax_kkr_max)
           tau2 = tauac((i-1)*kmax_kkr_max+1:i*kmax_kkr_max, &
              1:kmax_kkr_max)
         endif
!        ---------------------------------------------------------------
         call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
         CONE, J2, kmax_kkr_max, tau2, kmax_kkr_max, CZERO, tmp1, kmax_kkr_max)
!        ---------------------------------------------------------------
         call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
         CONE, tau1, kmax_kkr_max, tmp1, kmax_kkr_max, CZERO, tmp2, kmax_kkr_max)
!        ---------------------------------------------------------------
         call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, &
         CONE, J1, kmax_kkr_max, tmp2, kmax_kkr_max, CZERO, tmp3, kmax_kkr_max)
!        ---------------------------------------------------------------
         do L = 1, kmax_kkr_max
           sigma011 = sigma011 + coeff*tmp3(L,L)
         enddo
       enddo
     enddo
   enddo       
 
   end function calSigmaTildeSRO011
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeSRO10(n, dir, dir1, is, eval, caltype) result(sigma10)
!  ===================================================================

   use CPAMediumModule, only : getSingleSiteTmat
   use SystemVolumeModule, only : getAtomicVPVolume
   use AtomModule, only : getLocalSpeciesContent, getLocalNumSpecies
   use CurrentMatrixModule, only : getJMatrix
   use CrystalMatrixModule, only : calSigmaIntegralSRO
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, dir, dir1, is, caltype
   complex (kind=CmplxKind), intent(in) :: eval 

   integer (kind=IntKind) :: ic1, ic2, num_species
   real (kind=RealKind) :: Omega, c_a, c_b, coeff
   complex (kind=CmplxKind) :: sigma10
   complex (kind=CmplxKind), pointer :: J1(:,:), J2(:,:)

   Omega = getAtomicVPVolume(n)
   num_species = getLocalNumSpecies(n)
   sigma10 = CZERO

   do ic1 = 1, num_species
     do ic2 = 1, num_species
       c_a = getLocalSpeciesContent(n, ic1)
       c_b = getLocalSpeciesContent(n, ic2)
       coeff = -(c_a*c_b)/(PI*Omega)
       if (caltype == 1 .or. caltype == 4) then
         J1 => getJMatrix(n, ic1, is, dir, caltype, 0)
         J2 => getJMatrix(n, ic2, is, dir1, caltype, 1)
       else if (caltype == 2) then
         J1 => getJMatrix(n, ic1, is, dir, 3, 0)
         J2 => getJMatrix(n, ic2, is, dir1, 2, 1)
       else if (caltype == 3) then
         J1 => getJMatrix(n, ic1, is, dir, 2, 0)
         J2 => getJMatrix(n, ic2, is, dir1, 3, 1)
       else
         call ErrorHandler('calSigmaTildeSRO10', 'Incorrect caltype (1-4)', caltype)
       endif
       sigma10 = sigma10 + coeff*calSigmaIntegralSRO(n, eval, ic1, ic2, J1, J2, &
           getSingleSiteTmat, tau_needed=.true.,use_tmat=.true.,caltype=caltype)
       nullify(J1, J2)
     enddo
   enddo  
   
   end function calSigmaTildeSRO10
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeSROConductivity(n, is, delta, pot_type, dirnum, e)
!  ===================================================================   

   integer (kind=IntKind), intent(in) :: n, is, pot_type, dirnum
   real (kind=RealKind), intent(in) :: delta
   complex (kind=CmplxKind), intent(in) :: e

   integer (kind=IntKind) :: etype
   integer (kind=IntKind) :: dir, dir1
   complex(kind=CmplxKind) :: int_val(4)

   do dir = 1, dirnum
     do dir1 = 1, dirnum
       int_val = CZERO
       do etype = 1, 4
         int_val(etype) = calSigmaTildeSRO00(n, dir, dir1, is, etype) &
            + calSigmaTildeSRO010(n, dir, dir1, is, etype) &
            + calSigmaTildeSRO011(n, dir, dir1, is, etype) &
            + calSigmaTildeSRO10(n, dir, dir1, is, e, etype)
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
   subroutine computeCPAConductivity(n, is, delta, pot_type, dirnum, e)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: n, is, pot_type, dirnum
   real (kind=RealKind), intent(in) :: delta
   complex (kind=CmplxKind), intent(in) :: e

   integer (kind=IntKind) :: etype
   integer (kind=IntKind) :: dir, dir1
   complex(kind=CmplxKind) :: int_val(4)

   do dir = 1, dirnum
     do dir1 = 1, dirnum
       int_val = CZERO
       do etype = 1, 4
         int_val(etype) = &
          calSigmaTildeCPA1(n, dir, dir1, is, e, etype) + &
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
   subroutine calConductivity(LocalNumAtoms, n_spin_pola)
!  ===================================================================


   use SSSolverModule, only : solveSingleScattering
   use ScfDataModule, only : getFermiEnergyImagPart, isFermiEnergyRealPart, &
     useStepFunctionForSigma, useCubicSymmetryForSigma, getFermiEnergyRealPart
   use ValenceDensityModule, only : getFermiEnergy
   use WriteMatrixModule, only : writeMatrix
   use CrystalMatrixModule, only : calCrystalMatrix, retrieveTauSRO
   use CPAMediumModule, only : computeCPAMedium, populateBigTCPA, getSingleSiteMatrix
   use SROModule, only : calSpeciesTauMatrix, calNegatives
   use CurrentMatrixModule, only : calCurrentMatrix

   integer (kind=IntKind), intent(in) :: LocalNumAtoms, n_spin_pola
   integer (kind=IntKind) :: id, is, pot_type, dirnum
   real (kind=RealKind) :: delta, efermi
   complex (kind=CmplxKind) :: eval

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
           call calNegatives(id)
!          --------------------------------------------------------------
         else if (scf == 1) then
!          --------------------------------------------------------------
           call computeCPAMedium(eval, do_sro = .true.)
           call calNegatives(id)
!          --------------------------------------------------------------
         endif
      !  ----------------------------------------------------------------
         call calCurrentMatrix(id,is,eval,pot_type,mode)
      !  ---------------------------------------------------------------- 
         call computeSROConductivity(id, is, delta, pot_type, dirnum, eval)
!        -------------------------------------------------------------- 
       endif
     enddo
   enddo

!  Print calculated results to output file
!  --------------------------------------------------------------------
   call writeMatrix('sigma', sigma, dirnum, dirnum, n_spin_pola)
   call writeMatrix('sigmatilde', sigmatilde, dirnum, dirnum, n_spin_pola)
   call writeMatrix('sigmatilde2', sigmatilde2, dirnum, dirnum, n_spin_pola)
   call writeMatrix('sigmatilde3', sigmatilde3, dirnum, dirnum, n_spin_pola)
   call writeMatrix('sigmatilde4', sigmatilde4, dirnum, dirnum, n_spin_pola)
!  --------------------------------------------------------------------

   end subroutine calConductivity
!  ===================================================================
end module ConductivityModule
