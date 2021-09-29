module ConductivityModule
   use KindParamModule, only : IntKind, QuadRealKind, QuadCmplxKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : PI, ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors

public :: initConductivity,        &
          endConductivity,   &
          computeSROConductivity, &
          computeCPAConductivity, &
          calConductivity
!
!  Conductivity data stored here
!  --------------------------------------------------------------------
   complex (kind=CmplxKind), allocatable, public :: sigmatilde(:,:,:), sigmatilde2(:,:,:), &
                                          sigmatilde3(:,:,:), sigmatilde4(:,:,:)

   complex (kind=CmplxKind), allocatable, public :: sigma(:,:,:)
!  ---------------------------------------------------------------------

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

   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, kGID
   logical :: Initialized = .false.
! 
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initConductivity(num_atoms, lmaxkkr, lmaxphi, lmaxgreen,  &
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
   call initCurrentMatrixModule(num_atoms, lmaxkkr, lmaxphi, &
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
   subroutine computeCPAConductivity(is, dirnum)
!  ===================================================================
 
   use CPAConductivityModule, only : calVertexCorrectionMatrixCPA, &
      calSigmaTilde1VC, calSigmaTilde0

   integer (kind=IntKind), intent(in) :: is, dirnum

   integer (kind=IntKind) :: etype
   integer (kind=IntKind) :: dir, dir1
   complex (kind=CmplxKind) :: int_val(4)

   call calVertexCorrectionMatrixCPA()

   do dir = 1, dirnum
     do dir1 = 1, dirnum
       int_val = CZERO
       do etype = 1, 4
         int_val(etype) = &
           calSigmaTilde1VC(dir, dir1, etype) + &
           calSigmaTilde0(dir, dir1, etype)
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
   subroutine computeSROConductivity(is, dirnum)
!  ===================================================================   

   integer (kind=IntKind), intent(in) :: is, dirnum

   integer (kind=IntKind) :: etype
   integer (kind=IntKind) :: dir, dir1
   complex (kind=CmplxKind) :: int_val(4)

   do dir = 1, dirnum
     do dir1 = 1, dirnum
       int_val = CZERO
       do etype = 1, 4
         int_val(etype) = 0 ! template code
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
   subroutine calConductivity(efermi, LocalNumAtoms, n_spin_pola)
!  ===================================================================

   use KuboDataModule, only : useCubicSymmetryForSigma
   use CPAConductivityModule, only : initCPAConductivity, endCPAConductivity

   integer (kind=IntKind), intent(in) :: LocalNumAtoms, n_spin_pola
   integer (kind=IntKind) :: id, is, dirnum
   real (kind=RealKind), intent(in) :: efermi


   if (useCubicSymmetryForSigma()) then
     dirnum = 1
   else
     dirnum = 3
   endif

   do id = 1, LocalNumAtoms
     do is = 1, n_spin_pola
       if (mode == 3) then
      !  ---------------------------------------------------------------
         call initCPAConductivity(id, is, kmax_kkr_max, efermi, LocalNumAtoms)
      !  --------------------------------------------------------------- 
         call computeCPAConductivity(is, dirnum)
      !  ---------------------------------------------------------------
         call endCPAConductivity()
      !  --------------------------------------------------------------
       else if (mode == 4) then
      !  ---------------------------------------------------------------
      !  call initSROConductivity( ------- )   
      !  --------------------------------------------------------------- 
         call computeSROConductivity(is, dirnum)
!        --------------------------------------------------------------
      !  call endSROConductivity(  ------- )
      !  --------------------------------------------------------------
       endif
     enddo
   enddo
   
   if (useCubicSymmetryForSigma()) then
     sigma(2,2,:) = sigma(1,1,:)
     sigma(3,3,:) = sigma(1,1,:)
   endif

   end subroutine calConductivity
!  ===================================================================
end module ConductivityModule
