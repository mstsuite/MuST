module ConductivityModule
   use KindParamModule, only : IntKind, QuadRealKind, QuadCmplxKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : PI, ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors

public :: initConductivity,       &
          endConductivity,        &
          computeSROConductivity, &
          computeCPAConductivity, &
          calConductivity
!
!  Conductivity data stored here
!  --------------------------------------------------------------------
   complex (kind=CmplxKind), public :: sigma(3,3,2)
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
   integer (kind=IntKind) :: kmax_phi_max, kmax_green_max,               &
                             kmax_sigma, kmax_sigma_2, kmax_max, kmax_cg
   integer (kind=IntKind) :: NumSpecies
!
   integer (kind=IntKind) :: NumPEsInEGroup, MyPEinEGroup, eGID
!
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
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use RadialGridModule, only : getNumRmesh, getMaxNumRmesh
   use AtomModule, only : getLocalNumSpecies
   use PolyhedraModule, only : getNumPolyhedra
   use WriteMatrixModule, only : writeMatrix
   use ScfDataModule, only : isKKR, isLSMS, isKKRCPA, isKKRCPASRO, isSROSCF
   use CurrentMatrixModule, only : initCurrentMatrixModule
   use SROModule, only : getNeighSize
!
   implicit none
!
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
!
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
   if (maxval(iprint) >= 0) then
      write(6,'(/,a,i4)')'Number of processes for the energy parallelization: ',NumPEsInEGroup
   endif
! 
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

   if (mode == 1) then
     call ErrorHandler('initConductivity', &
                       'Conductivity for KKR not yet implemented')
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
   
!  -------------------------------------------------------------------
!  Initialize the arrays which store conductivity data
!  -------------------------------------------------------------------
   jsize = master_size

   allocate(lofk(kmax_kkr_max), mofk(kmax_cg),  &
         jofk(kmax_cg), m1m(-lmax_cg:lmax_cg))
   allocate(lofj((lmax_cg+1)*(lmax_cg+2)/2), mofj((lmax_cg+1)*(lmax_cg+2)/2))

   sigma = CZERO

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
!
   implicit none

   deallocate(print_instruction, lmax_kkr, lmax_phi, kmax_kkr, kmax_phi)  

!  ----------------------------------------------------
   call endCurrentMatrixModule(mode)
!  ----------------------------------------------------

   end subroutine endConductivity
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calConductivity(efermi, LocalNumAtoms, n_spin_pola)
!  ===================================================================
   use TimerModule, only : getTime
   use KuboDataModule, only : useCubicSymmetryForSigma
   use CPAConductivityModule, only : initCPAConductivity, endCPAConductivity, &
                                     computeCPAConductivity
   use LSMSConductivityModule, only : initLSMSConductivity, endLSMSConductivity, &
                                      computeLSMSConductivity
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms, n_spin_pola
   integer (kind=IntKind) :: id, is, dirnum, iprint, dir, dir1
   real (kind=RealKind), intent(in) :: efermi
   real (kind=RealKind) :: t0
!
   complex (kind=CmplxKind) :: sigmatilde(3,3,4)

   if (useCubicSymmetryForSigma()) then
      dirnum = 1
   else
      dirnum = 3
   endif

   iprint = maxval(print_instruction)

   if (mode == 2) then
!     Do LSMS Conductivity
!     t0 = getTime()
!     ----------------------------------------------------------------
      call initLSMSConductivity(dirnum, n_spin_pola, kmax_kkr_max, efermi, iprint)
!     ----------------------------------------------------------------
!     print *,'Calling initLSMSConductivity, timing = ',getTime() - t0
      do is = 1, n_spin_pola
!        -------------------------------------------------------------
         call computeLSMSConductivity(is, dirnum, sigmatilde)
!        -------------------------------------------------------------
         if (iprint >= 0) then
            call writeMatrix('sigmatilde1', sigmatilde(:,:,1), 3, 3)
            call writeMatrix('sigmatilde2', sigmatilde(:,:,2), 3, 3)
            call writeMatrix('sigmatilde3', sigmatilde(:,:,3), 3, 3)
            call writeMatrix('sigmatilde4', sigmatilde(:,:,4), 3, 3)
         endif
!
         do dir1 = 1, dirnum
            do dir = 1, dirnum
               sigma(dir,dir1,is) = 0.25d0*(sigmatilde(dir,dir1,1) -        &
                                            sigmatilde(dir,dir1,2) -        &
                                            sigmatilde(dir,dir1,3) +        &
                                            sigmatilde(dir,dir1,4))
            enddo
         enddo
      enddo
      call endLSMSConductivity()
   else
      do id = 1, LocalNumAtoms
         do is = 1, n_spin_pola
            if (mode == 3) then
      !        ----------------------------------------------------------
               call initCPAConductivity(id, is, n_spin_pola, kmax_kkr_max, efermi, LocalNumAtoms, mode)
      !        ---------------------------------------------------------- 
               call computeCPAConductivity(is, dirnum, sigmatilde)
      !        ----------------------------------------------------------
!
               do dir1 = 1, dirnum
                  do dir = 1, dirnum
                     sigma(dir,dir1,is) = 0.25d0*(sigmatilde(dir,dir1,1) -        &
                                                  sigmatilde(dir,dir1,2) -        &
                                                  sigmatilde(dir,dir1,3) +        &
                                                  sigmatilde(dir,dir1,4))
                  enddo
               enddo
               call endCPAConductivity()
      !        ----------------------------------------------------------
            else if (mode == 4) then
      !        ----------------------------------------------------------
      !        Under development
               call initCPAConductivity(id, is, n_spin_pola, kmax_kkr_max, efermi, LocalNumAtoms, mode)
      !        ----------------------------------------------------------
            endif
            if (iprint >= 0) then
               call writeMatrix('sigmatilde1', sigmatilde(:,:,1), 3, 3)
               call writeMatrix('sigmatilde2', sigmatilde(:,:,2), 3, 3)
               call writeMatrix('sigmatilde3', sigmatilde(:,:,3), 3, 3)
               call writeMatrix('sigmatilde4', sigmatilde(:,:,4), 3, 3)
            endif
         enddo
      enddo
   endif
   
   if (useCubicSymmetryForSigma()) then
      sigma(2,2,:) = sigma(1,1,:)
      sigma(3,3,:) = sigma(1,1,:)
   endif

   if (iprint >= 0) then
      call writeMatrix('sigma', sigma, 3, 3, n_spin_pola)
   endif

   end subroutine calConductivity
!  ===================================================================
end module ConductivityModule
