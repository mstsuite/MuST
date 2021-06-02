module SROModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
  use ErrorHandlerModule, only : ErrorHandler, WarningHandler
  use PublicTypeDefinitionsModule, only : NeighborStruct

!
public :: initSROMatrix,             &
          endSROMatrix,              &
          averageSROMatrix,          &
          generateBigTAMatrix,       &
          generateBigTCPAMatrix,     &
          obtainPosition,            &
          obtainNeighborIndex,       &
          populateTau,               &
          assembleTauFromBlocks,     &
          calculateImpurityMatrix,   &
          calSpeciesTauMatrix,       &
          clusterDtilde,             &
          calculateSCFSpeciesTerm,   &
          getKauFromTau,             &
          calculateNewTCPA,          &
          getSROMatrix,              &
          getDoubleSpeciesTauMatrix, &
          getSROParam,               &
          calNegatives,              &
          getNeighSize
!

private
   integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
   integer (kind=IntKind) :: nSpinCant, nSpinPola
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind) :: kmax_kkr_max
   integer (kind=IntKind) :: ndim_Tmat
   integer (kind=IntKind) :: print_instruction
   integer (kind=IntKind) :: sigma
   integer (kind=IntKind), allocatable :: lofk(:), mofk(:), jofk(:) 
!
   complex(kind=CmplxKind), allocatable, target :: WORK0_sro(:), WORK1_sro(:), WORK2_sro(:)
   complex(kind=CmplxKind), pointer :: tm(:,:), tm0(:,:), tm1(:,:), tm2(:,:)
   complex(kind=CmplxKind), allocatable :: z(:,:)
   complex(kind=CmplxKind), allocatable :: y(:,:)
!
   type TauBlockStruct
      complex (kind=CmplxKind), pointer :: tau_neighbor(:,:,:)
   end type TauBlockStruct

!  Ensure data is not overwritten

   type TmatBlockStruct
      complex (kind=CmplxKind), pointer :: tmat(:,:)
      complex (kind=CmplxKind), pointer :: tmat_inv(:,:)
      complex (kind=CmplxKind), pointer :: tmat_tilde_inv(:,:)
      complex (kind=CmplxKind), pointer :: tmat_tilde_inv_nn(:,:)
      complex (kind=CmplxKind), pointer :: T_inv(:,:)
      complex (kind=CmplxKind), pointer :: T_sigma_inv(:,:,:)
      complex (kind=CmplxKind), pointer :: proj_a(:,:)
      complex (kind=CmplxKind), pointer :: proj_b(:,:)
   end type TmatBlockStruct
!
   type SROTMatrixStruct
      real (kind=RealKind), pointer :: sro_param_a(:)
      real (kind=RealKind), pointer :: sro_param_a_nn(:)
      type (TmatBlockStruct), allocatable :: tmat_s(:)
      complex (kind=CmplxKind), pointer :: tau_ab(:,:,:)
      complex (kind=CmplxKind), pointer :: tau_abc(:,:,:) ! if tauab = tau_ab(e+id), then tau_abc = tau_ab(e-id)
      complex (kind=CmplxKind), pointer :: D_ab(:,:,:)
      complex (kind=CmplxKind), pointer :: Dt_ab(:,:,:)
      complex (kind=CmplxKind), pointer :: D_abc(:,:,:)
      complex (kind=CmplxKind), pointer :: Dt_abc(:,:,:)
      complex (kind=CmplxKind), pointer :: tau_sigma(:,:,:,:)
      complex (kind=CmplxKind), pointer :: tau_sigmac(:,:,:,:) 
      complex (kind=CmplxKind), pointer :: kau11(:,:,:) 
   end type SROTMatrixStruct
!
   type SROMediumStruct
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: num_species
      integer (kind=IntKind) :: blk_size
      integer (kind=IntKind) :: neigh_size
      real (kind=RealKind) :: CPosition(3)
      logical :: isCPA
      type(NeighborStruct), pointer :: Neighbor
      type(SROTMatrixStruct), allocatable :: SROTMatrix(:)
      complex (kind=CmplxKind), pointer :: Tcpa(:,:)
      complex (kind=CmplxKind), pointer :: Tcpa_inv(:,:)
      complex (kind=CmplxKind), pointer :: T_CPA(:,:)
      complex (kind=CmplxKind), pointer :: T_CPA_inv(:,:)
      complex (kind=CmplxKind), pointer :: tau_cpa(:,:,:)
      complex (kind=CmplxKind), pointer :: tau_cpac(:,:,:) ! if tau_cpa = tau_cpa(e+id), then tau_cpac = tau_cpa(e-id)
      type(TauBlockStruct), pointer :: tau_c(:)
   end type SROMediumStruct
!
   type(SROMediumStruct), allocatable :: SROMedium(:)
   integer (kind=IntKind) :: next_near_option
   logical :: test_CPA = .false.
   logical :: test_pure = .false.


contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSROMatrix (cant, pola)
!  ===================================================================
   
   use MediumHostModule, only  : getNumSites, getLocalNumSites, &
                     getGlobalSiteIndex, getNumSpecies, getSpeciesContent
   use ScfDataModule, only : retrieveSROParams, isNextNearestSRO, &
                    isSROSCF, isSROCVM, retrieveCVMParams, isConductivity
   use NeighborModule, only : getNeighbor
   use SSSolverModule, only : getScatteringMatrix
   use SystemModule, only : getAtomPosition
   use WriteMatrixModule, only : writeMatrix
   use MatrixInverseModule, only : MtxInv_LU

   integer(kind=IntKind), intent(in) :: cant, pola
   integer(kind=IntKind) :: sro_param_nums, num, il, ic, is, ig, i, j, iter1, iter2, temp
   integer(kind=IntKind) :: in, jn
   integer(kind=IntKind) :: type, kl, jl, m, l, n, lmax_kkr_max
   real (kind=RealKind) :: spec_i, spec_j
   real(kind=RealKind), allocatable :: sro_params(:), sro_params_nn(:)
   real(kind=RealKind) :: cvm_sitei(2)
   real(kind=RealKind), allocatable :: cvm_params(:)
!
   logical :: isWarrenCowley = .false.

   if (isConductivity()) then
     sigma = 1
   else
     sigma = 0
   endif

!  --------------------------------------------------------
   next_near_option = isNextNearestSRO()
!  --------------------------------------------------------

   if (next_near_option == 0) then
   !  --------------------------------------------------------
      call retrieveSROParams(sro_params, sro_param_nums, isWC=isWarrenCowley)
   !  --------------------------------------------------------
   else 
   !  --------------------------------------------------------
      call retrieveSROParams(sro_param_list=sro_params,       &
                             param_num=sro_param_nums,        &
                             sro_param_list_nn=sro_params_nn, &
                             isWC=isWarrenCowley)
   !  --------------------------------------------------------
   endif

   if (isSROCVM()) then
      call retrieveCVMParams(cvm_params)
      cvm_sitei(1) = cvm_params(1)
      cvm_sitei(2) = 1 - cvm_params(1)
   endif

   GlobalNumSites = getNumSites()
   LocalNumSites = getLocalNumSites()
   if (.not. allocated(SROMedium)) then
      allocate(SROMedium(LocalNumSites))
   endif
   nSpinCant = cant
   nSpinPola = pola
 
   do il = 1, LocalNumSites
      SROMedium(il)%Neighbor => getNeighbor(il)
      SROMedium(il)%local_index = il
      ig = getGlobalSiteIndex(il)
      SROMedium(il)%CPosition = getAtomPosition(ig)
      SROMedium(il)%global_index = ig
      num = getNumSpecies(ig)
      SROMedium(il)%num_species = num
      SROMedium(il)%neigh_size = SROMedium(il)%Neighbor%NumAtoms+1

    ! allocate(SROMedium(il)%tau_c((SROMedium(il)%neigh_size)**2))
      
      if (num > 1) then
        SROMedium(il)%isCPA = .true.
      else
        SROMedium(il)%isCPA = .false.
      endif
      
      if (num*(num + 1)/2 /= sro_param_nums) then
         call ErrorHandler('initSROMatrix','Mismatch between number of alloy species and SRO parameters')
      else
         if (.not. allocated(SROMedium(il)%SROTMatrix)) then
            allocate(SROMedium(il)%SROTMatrix(num))
         endif
         do i = 1, num
            allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(nSpinCant**2))
            allocate(SROMedium(il)%SROTMatrix(i)%sro_param_a(num))
            if (next_near_option == 1) then
              allocate(SROMedium(il)%SROTMatrix(i)%sro_param_a_nn(num))
            endif
            if (isSROCVM()) then
               spec_i = cvm_sitei(i)
            else
               spec_i = getSpeciesContent(i, ig)
            endif
            do j = 1, num
              if (isSROCVM()) then
                 spec_j = cvm_sitei(j)
              else
                 spec_j = getSpeciesContent(j, ig)
              endif
              if (j < i) then
                 SROMedium(il)%SROTMatrix(i)%sro_param_a(j) = (spec_j/spec_i)*SROMedium(il)%SROTMatrix(j)%sro_param_a(i)
                 if (next_near_option == 1) then
                    SROMedium(il)%SROTMatrix(i)%sro_param_a_nn(j) =  &
                        (spec_j/spec_i)*SROMedium(il)%SROTMatrix(j)%sro_param_a_nn(i)
                 endif
              else
                 temp = (i - 1)*num - (i - 1)*(i - 2)/2
                 if (isWarrenCowley) then
                    SROMedium(il)%SROTMatrix(i)%sro_param_a(j) = spec_j*(ONE-sro_params(temp + j - i + 1))
                    if (next_near_option == 1) then
                       SROMedium(il)%SROTMatrix(i)%sro_param_a_nn(j) = spec_j*(ONE-sro_params_nn(temp + j - i + 1))
                    endif
                 else
                    SROMedium(il)%SROTMatrix(i)%sro_param_a(j) = sro_params(temp + j - i + 1)
                    if (next_near_option == 1) then
                       SROMedium(il)%SROTMatrix(i)%sro_param_a_nn(j) = sro_params_nn(temp + j - i + 1)
                    endif
                 endif
              endif
            enddo

!           Print *, SROMedium(il)%SROTMatrix(i)%sro_param_a
 
            tm => getScatteringMatrix('T-Matrix',spin=1,site=SROMedium(il)%local_index,atom=i)
            
            SROMedium(il)%blk_size = size(tm, 1)
            kmax_kkr_max = SROMedium(il)%blk_size

            if (nSpinCant == 2) then
               call ErrorHandler('initSROMatrix','SRO is not equipped to deal with spin canting yet')
            endif

            allocate(SROMedium(il)%SROTMatrix(i)%kau11(SROMedium(il)%blk_size, SROMedium(il)%blk_size, nSpinCant**2))
            allocate(SROMedium(il)%SROTMatrix(i)%tau_ab(SROMedium(il)%blk_size*SROMedium(il)%neigh_size,  &
                     SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
            if (sigma == 1) then
              allocate(SROMedium(il)%SROTMatrix(i)%tau_sigma(num, SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                 SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
              allocate(SROMedium(il)%SROTMatrix(i)%tau_abc(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                 SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
              allocate(SROMedium(il)%SROTMatrix(i)%tau_sigmac(num, SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                 SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
              allocate(SROMedium(il)%SROTMatrix(i)%D_ab(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                 SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
              allocate(SROMedium(il)%SROTMatrix(i)%D_abc(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                 SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
              allocate(SROMedium(il)%SROTMatrix(i)%Dt_ab(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                 SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
              allocate(SROMedium(il)%SROTMatrix(i)%Dt_abc(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                 SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
            endif
            
            do is = 1,nSpinCant**2
               allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%tmat_tilde_inv(SROMedium(il)%blk_size, &
                       SROMedium(il)%blk_size))
               allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%tmat_tilde_inv_nn(SROMedium(il)%blk_size, &
                       SROMedium(il)%blk_size))
               allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%T_inv(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                       SROMedium(il)%blk_size*SROMedium(il)%neigh_size))
               if (sigma == 1) then
                  allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%T_sigma_inv(num, &
                   SROMedium(il)%blk_size*SROMedium(il)%neigh_size, SROMedium(il)%blk_size*SROMedium(il)%neigh_size))
               endif
               if (isSROSCF() == 1) then
                  allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%proj_a(SROMedium(il)%blk_size, SROMedium(il)%blk_size))
                  allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%proj_b(SROMedium(il)%blk_size, SROMedium(il)%blk_size))
               endif
            enddo

         enddo
         allocate(SROMedium(il)%Tcpa(SROMedium(il)%blk_size, SROMedium(il)%blk_size))
         allocate(SROMedium(il)%Tcpa_inv(SROMedium(il)%blk_size, SROMedium(il)%blk_size))
         allocate(SROMedium(il)%T_CPA(SROMedium(il)%blk_size*SROMedium(il)%neigh_size,  &
                          SROMedium(il)%blk_size*SROMedium(il)%neigh_size))
         allocate(SROMedium(il)%T_CPA_inv(SROMedium(il)%blk_size*SROMedium(il)%neigh_size,   &
                          SROMedium(il)%blk_size*SROMedium(il)%neigh_size))
         if (sigma == 1) then
           allocate(SROMedium(il)%tau_cpac(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
               SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
         endif
         allocate(SROMedium(il)%tau_cpa(SROMedium(il)%blk_size*SROMedium(il)%neigh_size,  &
                          SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
         allocate(SROMedium(il)%tau_c((SROMedium(il)%neigh_size)**2))
         do in = 1, SROMedium(il)%neigh_size
           do jn = 1, SROMedium(il)%neigh_size
             index = jn + (in - 1)*SROMedium(il)%neigh_size
             allocate(SROMedium(il)%tau_c(index)%tau_neighbor(SROMedium(il)%blk_size,   & 
                     SROMedium(il)%blk_size, nSpinCant**2))
           enddo
         enddo
      endif
   enddo

   if (SROMedium(1)%Neighbor%NumShells < 2 .and. next_near_option == 1) then
      call ErrorHandler('initSROMatrix', 'No next near neighbors - increase LIZ cutoff')
   endif

!  Test to check whether next near neighbors are being detected   
!  do il = 1, SROMedium(1)%neigh_size-1
!      ------------------------------------------
!      call determineNeighborType(1, il, type) 
!      ------------------------------------------
!      if (type == 0) then
!         Print *, "This is a near neighbor"
!      else if (type == 1) then
!         Print *, "This is a next near neighbor"
!      else
!         Print *, "Can't detect"
!      endif
!  enddo

   allocate(z(SROMedium(1)%blk_size*SROMedium(1)%neigh_size, SROMedium(1)%blk_size*SROMedium(1)%neigh_size))
   allocate(y(SROMedium(1)%blk_size*SROMedium(1)%neigh_size, SROMedium(1)%blk_size*SROMedium(1)%neigh_size))
   
   allocate(WORK0_sro(SROMedium(1)%blk_size**2))
   allocate(WORK1_sro(SROMedium(1)%blk_size**2))
   allocate(WORK2_sro(SROMedium(1)%blk_size**2))

   lmax_kkr_max = int(sqrt(1.0*kmax_kkr_max)) - 1
   allocate(lofk(kmax_kkr_max), mofk(kmax_kkr_max),  &
         jofk(kmax_kkr_max))

   kl=0; jl = 0
   do l=0,lmax_kkr_max
      n=(l+1)*(l+2)/2-l
      do m=-l,l
         kl=kl+1
         lofk(kl)=l
         mofk(kl)=m
         jofk(kl)=n+abs(m)
      enddo
   enddo

    
   end subroutine initSROMatrix
!  =================================================================== 

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSROMatrix()
!  ===================================================================

   integer(kind=IntKind) :: n, ic, in
   
   do n = 1, LocalNumSites
     nullify(SROMedium(n)%Tcpa, SROMedium(n)%Neighbor,  & 
       SROMedium(n)%Tcpa_inv, SROMedium(n)%T_CPA &
       ,SROMedium(n)%T_CPA_inv, SROMedium(n)%tau_cpa)
     do in = 1, SROMedium(n)%neigh_size**2
       nullify(SROMedium(n)%tau_c(in)%tau_neighbor)
     enddo
     nullify(SROMedium(n)%tau_c)
     do ic = 1, SROMedium(n)%num_species
       nullify(SROMedium(n)%SROTMatrix(ic)%sro_param_a, & 
        SROMedium(n)%SROTMatrix(ic)%tau_ab, SROMedium(n)%SROTMatrix(ic)%kau11)
       deallocate(SROMedium(n)%SROTMatrix(ic)%tmat_s)
     enddo
     deallocate(SROMedium(n)%SROTMatrix)
   enddo
   deallocate(SROMedium, WORK0_sro, WORK1_sro, WORK2_sro)
   deallocate(z, y)
   deallocate(lofk, mofk, jofk)

   end subroutine endSROMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageSROMatrix (n, ia)
!  ===================================================================

!  NEED TO CHECK IF SITE N IS CPA OR NOT
!  BE CAREFUL WHETHER TO USE LOCAL OR GLOBAL INDEX ANYWHERE
   
   use MatrixInverseModule, only : MtxInv_LU
   integer (kind=IntKind), intent(in) :: n, ia
   integer (kind=IntKind) :: ic, is, nsize
   real (kind=RealKind) :: wab, wab_nn
   
   nsize = SROMedium(n)%blk_size
 
   do is = 1, nSpinCant**2
     SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv = CZERO
     SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv_nn = CZERO
     do ic = 1, SROMedium(n)%num_species
        wab = SROMedium(n)%SROTMatrix(ia)%sro_param_a(ic)
!       -------------------------------------------------------------------------------------
        call zaxpy(nsize*nsize,wab,SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%tmat,1,&
                  SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv, 1)
!       -------------------------------------------------------------------------------------
     enddo
     if (next_near_option == 1) then
        do ic = 1, SROMedium(n)%num_species
           wab_nn = SROMedium(n)%SROTMatrix(ia)%sro_param_a_nn(ic)
!          ----------------------------------------------------------------------------------
           call zaxpy(nsize*nsize,wab_nn,SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%tmat,1, &
                 SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv_nn, 1)
!          ----------------------------------------------------------------------------------
        enddo
     endif

!    Inverting to get tilde_inv
!    ------------------------------------------------------------------------------
     call MtxInv_LU(SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv, nsize)
!    ------------------------------------------------------------------------------
     if (next_near_option == 1) then
!      ---------------------------------------------------------------------------------
       call MtxInv_LU(SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv_nn, nsize)
!      ---------------------------------------------------------------------------------
     endif
   enddo   
   

   end subroutine averageSROMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generateBigTAMatrix (n, ia)
!  ===================================================================

   use AtomModule, only : getLocalNumSpecies

   integer (kind=IntKind), intent(in) :: n, ia
   integer (kind=IntKind) :: nsize, delta, total_size, i, is,iter1,iter2, tmp
   integer (kind=IntKind) :: ttype, ic, num_species

   delta = SROMedium(n)%neigh_size
   nsize = SROMedium(n)%blk_size
   total_size = nsize*delta
   num_species = getLocalNumSpecies(n)   

   do is = 1, nSpinCant
      SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv = CZERO
   
      do iter2 = 1, nsize
         do iter1 = 1, nsize
            SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1,iter2) =  &
              SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_inv(iter1,iter2)
            if (sigma == 1) then
              do ic = 1, num_species
                SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_sigma_inv(ic,iter1,iter2) &
                = SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_inv(iter1,iter2)
              enddo
            endif
         enddo
      enddo
   
      if (next_near_option == 1) then
         do i = 2, delta
            tmp = (i - 1)*nsize
            call determineNeighborType(n, i - 1, ttype)
            if (ttype == 0) then
               do iter2 = 1, nsize
                  do iter1 = 1, nsize
                     SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1+tmp,iter2+tmp) = &
                      SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv(iter1, iter2)
                  enddo
               enddo
            else if (ttype == 1) then
               do iter2 = 1, nsize
                  do iter1 = 1, nsize
                     SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1+tmp,iter2+tmp) = &
                      SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv_nn(iter1, iter2)
                  enddo
               enddo
            endif
         enddo 
      
      else
         do i = 2, delta
            tmp = (i - 1)*nsize
            do iter2 = 1, nsize
               do iter1 = 1, nsize
                 if (test_CPA) then
                    SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1+tmp,iter2+tmp) = &
                     SROMedium(n)%Tcpa_inv(iter1, iter2)
               
                 else if (test_pure) then
                    if (i < 3) then
                       SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1+tmp, iter2+tmp) = &
                         SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv(iter1, iter2)
                    else
                       SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1+tmp, iter2+tmp) = &
                         SROMedium(n)%Tcpa_inv(iter1, iter2)
                    endif
              
                 else
                   SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1+tmp,iter2+tmp) = & 
                     SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv(iter1, iter2)
                 endif

                 if (sigma == 1) then
                   do ic = 1, num_species
                     SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_sigma_inv(ic,iter1+tmp,iter2+tmp) &
                     = SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%tmat_inv(iter1,iter2)
                   enddo
                 endif
               enddo
            enddo
         enddo
    
      endif
   enddo
 
!  call writeMatrix('Big-TA', SROMedium(1)%SROTMatrix(ia)%tmat_s(1)%T_inv, nsize*delta, nsize*delta) 

   end subroutine generateBigTAMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generateBigTCPAMatrix (n, Tcpa_in)
!  ===================================================================

   use MatrixInverseModule, only : MtxInv_LU
!   
   integer (kind=IntKind), intent(in) :: n
   complex (kind=CmplxKind), intent(in), pointer :: Tcpa_in(:, :)
!
   integer (kind=IntKind) :: i, iter1, iter2, tmp, nsize, delta, total_size
!
   SROMedium(n)%Tcpa_inv = Tcpa_in
   nsize = SROMedium(n)%blk_size
   delta = SROMedium(n)%neigh_size
   SROMedium(n)%T_CPA_inv = CZERO
    
   do i = 1, delta
       tmp = (i - 1)*nsize
       do iter2 = 1, nsize
          do iter1 = 1, nsize
             SROMedium(n)%T_CPA_inv(iter1+tmp,iter2+tmp) = SROMedium(n)%Tcpa_inv(iter1,iter2)
          enddo
       enddo
   enddo
   
   end subroutine generateBigTCPAMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine obtainPosition (local_index, p_vec, jn)
!  ===================================================================

!  use Atom2ProcModule, only : getLocalIndex

   integer (kind=IntKind), intent(in) :: local_index, jn
   real (kind=RealKind), intent(out) :: p_vec(3)

   integer (kind=IntKind) :: il

!  il = getLocalIndex(global_index)
   if (jn == 1) then
     p_vec = SROMedium(local_index)%CPosition
   else
     p_vec = SROMedium(local_index)%Neighbor%Position(1:3, jn - 1)
   endif

   end subroutine obtainPosition
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function obtainNeighborIndex (local_index, p_vec)  result(neigh_index)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: local_index
   real (kind=RealKind), intent(in) :: p_vec(3)

   real (kind=RealKind) :: temp(3)
   integer (kind=IntKind) :: in, num

   integer (kind=IntKind) :: neigh_index

   neigh_index = 1
   num = SROMedium(local_index)%neigh_size
   
   if (p_vec(1) == 0 .and. p_vec(2) == 0 .and. p_vec(3) == 0) then
     neigh_index = 1
   else
     do in = 2, num
        temp = SROMedium(local_index)%Neighbor%Position(1:3, in - 1)
        if (p_vec(1) == temp(1) .and. p_vec(2) == temp(2) .and. p_vec(3) == temp(3)) then
           neigh_index = in
           EXIT
        endif
     enddo
   endif

   end function obtainNeighborIndex
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine determineNeighborType(local_index, neigh_index, type)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: local_index, neigh_index
   integer (kind=IntKind), intent(out) :: type
   real (kind=RealKind) :: norm = 0

   norm = SROMedium(local_index)%Neighbor%Position(1, neigh_index)**2 + &
      SROMedium(local_index)%Neighbor%Position(2, neigh_index)**2 + &
      SROMedium(local_index)%Neighbor%Position(3, neigh_index)**2 
   
   if (SQRT(norm) == SROMedium(local_index)%Neighbor%ShellRad(1)) then
      type = 0
   else if (SQRT(norm) == SROMedium(local_index)%Neighbor%ShellRad(2)) then
      type = 1
   else
      call ErrorHandler('determineNeighborType',   &
            'Unknown type of Neighbor - check LIZ cutoff')
   endif

   end subroutine determineNeighborType
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine populateTau (tau, jindex, nindex)
!  ===================================================================

   complex(kind=CmplxKind), intent(in) :: tau(:,:,:)
   integer(kind=IntKind), intent(in) :: jindex, nindex

!
   if (jindex < 1 .or. jindex > LocalNumSites) then
     call ErrorHandler('populateTau', 'Invalid local index', jindex)
   else
     if (nindex < 1 .or. nindex > (SROMedium(jindex)%neigh_size)**2) then
       call ErrorHandler('populateTau', 'Invalid neighbor pair index', nindex)
     else
!      Print *, "Populating nindex", nindex 
       SROMedium(jindex)%tau_c(nindex)%tau_neighbor = tau
     endif
   endif
!
   end subroutine populateTau
!  ==================================================================

!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine assembleTauFromBlocks (n)
!  ==================================================================
   
   integer(kind=IntKind), intent(in) :: n
!
   integer(kind=IntKind) :: in, jn, i, j,index, nsize, dsize, is
 
   nsize = SROMedium(n)%neigh_size
   dsize = SROMedium(n)%blk_size

   SROMedium(n)%tau_cpa = CZERO
   
   do is = 1, nSpinCant
     do jn = 1, nsize
       do in = 1, nsize
         index = in + (jn - 1)*nsize
         do j = 1, dsize
           do i = 1, dsize
             SROMedium(n)%tau_cpa(i + (jn - 1)*dsize, j + (in - 1)*dsize, is)   & 
                = SROMedium(n)%tau_c(index)%tau_neighbor(i, j, is)
           enddo
         enddo
       enddo
     enddo
   enddo 

   end subroutine assembleTauFromBlocks
!  ==================================================================

!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calculateImpurityMatrix(n, ic)
!  ==================================================================

   use MatrixInverseModule, only : MtxInv_LU
   use WriteMatrixModule, only : writeMatrix
   use MatrixModule, only : computeAprojB
   use AtomModule, only : getLocalNumSpecies
!   
   integer(kind=IntKind), intent(in) :: n, ic
!
   integer(kind=IntKind) :: dsize, nsize, is, ic1

   SROMedium(n)%SROTMatrix(ic)%tau_ab = CZERO 

   dsize = SROMedium(n)%blk_size
   nsize = SROMedium(n)%neigh_size

!  do is = 1, nSpinCant**2
   y = CZERO
   z = SROMedium(n)%SROTMatrix(ic)%tmat_s(1)%T_inv - SROMedium(n)%T_CPA_inv

   call computeAprojB('L', dsize*nsize, SROMedium(n)%tau_cpa(:,:, 1), z, y)

   SROMedium(n)%SROTMatrix(ic)%tau_ab(:, :, 1) = y

   if (sigma == 1) then
     call computeAprojB('N', dsize*nsize, SROMedium(n)%tau_cpa(:,:, 1), z, &
                SROMedium(n)%SROTMatrix(ic)%D_ab(:,:,1))
     call computeAprojB('N', dsize*nsize, z, SROMedium(n)%tau_cpa(:,:,1), &
                SROMedium(n)%SROTMatrix(ic)%Dt_ab(:,:,1))
   endif

!  call writeMatrix('tau_a11', SROMedium(n)%SROTMatrix(ic)%tau_ab(1:dsize, 1:dsize, 1), &
!                 dsize, dsize, TEN2m8) 
!  enddo

   if (sigma == 1) then
     do ic1 = 1, getLocalNumSpecies(n)
       y = CZERO
       z = SROMedium(n)%SROTMatrix(ic)%tmat_s(1)%T_sigma_inv(ic1,:,:) - &
           SROMedium(n)%T_CPA_inv
       call computeAprojB('L', dsize*nsize, SROMedium(n)%tau_cpa(:,:,1), z, y)
       SROMedium(n)%SROTMatrix(ic)%tau_sigma(ic1,:,:,1) = y
     enddo
   endif

   end subroutine calculateImpurityMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function clusterDtilde(n, ic, is, kvec, caltype, etype) result(DtildeK)
!  ===================================================================

   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, ic, is, caltype, etype
   real (kind=RealKind), intent(in) :: kvec(3)
   integer (kind=IntKind) :: i, j, L1, L2, iter1, iter2, nsize, dsize, istart, iend
   real (kind=RealKind) :: Rp(3)
   complex (kind=CmplxKind) :: kp, exp_term
   complex (kind=CmplxKind) :: DtildeK(kmax_kkr_max,kmax_kkr_max)
   complex (kind=CmplxKind), allocatable :: tmp1(:,:), tmp2(:,:)
   complex (kind=CmplxKind), allocatable :: iden(:,:), D(:,:), Dp(:,:), Tdiff(:,:), &
                              Tdiffc(:,:), tauc(:,:), taucc(:,:), Dc(:,:)

   Rp = ZERO
   nsize = SROMedium(n)%neigh_size
   dsize = SROMedium(n)%blk_size

   allocate(D(nsize*dsize, nsize*dsize), Tdiff(nsize*dsize, nsize*dsize), &
     Dc(nsize*dsize, nsize*dsize), tauc(nsize*dsize, nsize*dsize), &
     taucc(nsize*dsize, nsize*dsize), Tdiffc(nsize*dsize, nsize*dsize))
   allocate(Dp(dsize, dsize), iden(dsize, dsize), tmp1(dsize, dsize), tmp2(dsize, dsize))
   D = CZERO; Tdiff = CZERO; Dc = CZERO
   Dp = CZERO; DtildeK = CZERO
   iden = CZERO; tauc = CZERO; taucc = CZERO
   tmp1 = CZERO; tmp2 = CZERO
 
   do i = 1, dsize
     iden(i, i) = CONE
   enddo

!  call writeMatrix('iden', iden, dsize, dsize)
   Tdiff = SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%T_inv - &
           SROMedium(n)%T_CPA_inv
   Tdiffc = conjg(Tdiff)
   tauc = SROMedium(n)%tau_cpa(:,:,is)
   taucc = SROMedium(n)%tau_cpac(:,:,is)

   if (caltype == 0) then
     D = SROMedium(n)%SROTMatrix(ic)%D_ab(:,:,is)
     Dc = SROMedium(n)%SROTMatrix(ic)%D_abc(:,:,is)
   ! call computeAprojB('N', dsize*nsize, tauc, Tdiff, D)
   ! call computeAprojB('N', dsize*nsize, taucc, Tdiffc, Dc)
   else if (caltype == 1) then
     D = SROMedium(n)%SROTMatrix(ic)%Dt_ab(:,:,is)
     Dc = SROMedium(n)%SROTMatrix(ic)%Dt_abc(:,:,is)
   ! call computeAprojB('N', dsize*nsize, Tdiff, SROMedium(n)%tau_cpa(:,:,is), D)
   ! call computeAprojB('N', dsize*nsize, Tdiffc, taucc, Dc)
   endif

!  Print *, kvec
   do i = 1, nsize
     tmp1 = CZERO
     istart = (i-1)*dsize + 1
     iend = i*dsize
     call obtainPosition(n, Rp, i)
     kp = kvec(1)*Rp(1) + kvec(2)*Rp(2) + kvec(3)*Rp(3)
!    Print *, Rp
     if (caltype == 0) then
       exp_term = exp(sqrtm1*kp)
       if (etype == 1) then
         tmp1 = D(1:dsize, istart:iend)
       else if (etype == 2) then
         tmp1 = Dc(1:dsize, istart:iend)
       endif
       call zgemm('n', 'n', dsize, dsize, dsize, CONE, tmp1, &
            dsize, iden, dsize, CZERO, Dp, dsize)
     else if (caltype == 1) then
       exp_term = exp(-sqrtm1*kp)
       if (etype == 1) then
         tmp1 = D(istart:iend, 1:dsize)
       else if (etype == 2) then
         tmp1 = Dc(istart:iend, 1:dsize)
       endif
       call zgemm('n', 'n', dsize, dsize, dsize, CONE, tmp1, &
          dsize, iden, dsize, CZERO, Dp, dsize)
     endif
     call zgemm('n', 'n', dsize, dsize, dsize, exp_term, Dp, dsize, &
         iden, dsize, CONE, DtildeK, dsize)
  !  call writeMatrix('Dp', Dp, dsize, dsize)
  !  Print *, exp_term
   enddo

!  call writeMatrix('DtildeK', DtildeK, dsize, dsize) 
!  call ErrorHandler('clusterDtilde', 'stop')

   end function clusterDtilde
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calNegatives(n)
!  ===================================================================

   use MatrixModule, only : computeAprojB
   
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: dsize, nsize, ic1, ic2, i, j, L1, L2
   complex (kind=CmplxKind), allocatable :: taucc(:,:), tauc(:,:)
   complex (kind=CmplxKind), allocatable :: Tdiff(:,:), Tdiffc(:,:)
   complex (kind=CmplxKind), allocatable :: Tdiffab(:,:), Tdiffabc(:,:)
   complex (kind=CmplxKind), allocatable :: tauac(:,:), tausigmac(:,:), tmp1(:,:), &
                                        tmp2(:,:), tmp3(:,:), tmp4(:,:)
   
   nsize = SROMedium(n)%neigh_size
   dsize = SROMedium(n)%blk_size
   allocate(tauc(nsize*dsize, nsize*dsize), taucc(nsize*dsize, nsize*dsize))
   allocate(Tdiff(nsize*dsize, nsize*dsize), Tdiffc(nsize*dsize, nsize*dsize))
   allocate(Tdiffab(nsize*dsize, nsize*dsize), Tdiffabc(nsize*dsize, nsize*dsize))
   allocate(tauac(nsize*dsize, nsize*dsize), tausigmac(nsize*dsize, nsize*dsize))
   allocate(tmp1(nsize*dsize, nsize*dsize), tmp2(nsize*dsize, nsize*dsize), &
      tmp3(nsize*dsize, nsize*dsize), tmp4(nsize*dsize, nsize*dsize))
   taucc = CZERO
   Tdiff = CZERO
   Tdiffc = CZERO
   tausigmac = CZERO
   tauc = SROMedium(n)%tau_cpa(:,:,1)
   
   do i = 1, nsize
     do j = 1, nsize
       do L2 = 1, dsize
         do L1 = 1, dsize
           taucc((i-1)*dsize + L1, (j-1)*dsize + L2)  = &  
            (-1.0)**(lofk(L2) - lofk(L1))*conjg(tauc((j-1)*dsize + L2, (i-1)*dsize + L1)) 
         enddo
       enddo
     enddo
   enddo
   
   SROMedium(n)%tau_cpac(:,:,1) = taucc
   
   do ic1 = 1, SROMedium(n)%num_species
     tmp1 = CZERO; tmp2 = CZERO; tmp3 = CZERO; tmp4 = CZERO
     tauac = CZERO; Tdiff = CZERO; Tdiffc = CZERO
     Tdiff = SROMedium(n)%SROTMatrix(ic1)%tmat_s(1)%T_inv - &
              SROMedium(n)%T_CPA_inv
     Tdiffc = conjg(Tdiff)
     call computeAprojB('L',dsize*nsize, taucc, Tdiffc, tmp1)
     call computeAprojB('N',dsize*nsize, taucc, Tdiffc, tmp2)
     call computeAprojB('N',dsize*nsize, Tdiffc, taucc, tmp3)

     SROMedium(n)%SROTMatrix(ic1)%tau_abc(:,:,1) = tmp1
     SROMedium(n)%SROTMatrix(ic1)%D_abc(:,:,1) = tmp2
     SROMedium(n)%SROTMatrix(ic1)%Dt_abc(:,:,1) = tmp3
     do ic2 = 1, SROMedium(n)%num_species
       tmp4 = CZERO; tausigmac = CZERO; Tdiffab = CZERO; Tdiffabc = CZERO
       Tdiffab = SROMedium(n)%SROTMatrix(ic1)%tmat_s(1)%T_sigma_inv(ic2,:,:) - &
             SROMedium(n)%T_CPA_inv
       Tdiffabc = conjg(Tdiffab)
       call computeAprojB('L',dsize*nsize, taucc, Tdiffabc, tmp4)
       SROMedium(n)%SROTMatrix(ic1)%tau_sigmac(ic2,:,:,1) = tmp4
     enddo
   enddo
   
   end subroutine calNegatives
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calculateSCFSpeciesTerm (n, ic, c_ic)
!  ===================================================================   
   
   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, ic
   real (kind=RealKind), intent(in) :: c_ic
   integer (kind=IntKind) :: dsize, nsize, is

   complex (kind=CmplxKind), allocatable :: tmp(:,:), proj_c(:,:)

   dsize = SROMedium(n)%blk_size
   nsize = SROMedium(n)%neigh_size
   allocate(tmp(dsize*nsize, dsize*nsize))
   allocate(proj_c(dsize*nsize, dsize*nsize))

   do is = 1, nSpinCant**2
      y = CZERO
      z = SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%T_inv - SROMedium(n)%T_CPA_inv

      call zgemm('N', 'n', dsize*nsize, dsize*nsize, dsize*nsize,    &
        c_ic, z, dsize*nsize, SROMedium(n)%SROTMatrix(ic)%tau_ab(1,1,1), dsize*nsize,  &
        CZERO, tmp, dsize*nsize)

!     call writeMatrix('Ta(1 + tau(ta - tcpa))^-1', tmp, dsize*nsize, dsize*nsize, TEN2m8) 
!     ------------------------------------------------------------------------
      call zgemm('N', 'n', dsize*nsize, dsize*nsize, dsize*nsize,    &
        CONE, SROMedium(n)%tau_cpa(1,1,1), dsize*nsize, tmp, dsize*nsize, CZERO, proj_c, dsize*nsize)
!     ------------------------------------------------------------------------
!     call writeMatrix('tauTa(1 + tau(ta - tcpa))^-1', proj_c(1:dsize, 1:dsize), dsize, dsize, TEN2m8)
      SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%proj_a = proj_c(1:dsize, 1:dsize)
   enddo 

!  call writeMatrix('tau_a11', SROMedium(n)%SROTMatrix(ic)%tau_ab(1:dsize,1:dsize,1), dsize, dsize, TEN2m8) 

   end subroutine calculateSCFSpeciesTerm
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calculateNewTCPA(n)  result(total_proj)
!  ===================================================================

   use MatrixInverseModule, only : MtxInv_LU
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: dsize, ic
   
   complex (kind=CmplxKind), allocatable :: tau_inv(:,:), temp(:,:), temp2(:,:)
   complex (kind=CmplxKind) :: total_proj(SROMedium(n)%blk_size,SROMedium(n)%blk_size)

   dsize = SROMedium(n)%blk_size
   allocate(tau_inv(dsize, dsize), temp(dsize, dsize), temp2(dsize, dsize))

   total_proj = CZERO
   temp = CZERO

   do is = 1, nSpinCant**2
     do ic = 1, SROMedium(n)%num_species
!       -------------------------------------------------------------
        call zaxpy(dsize*dsize, CONE, SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%proj_a, &
             1, temp, 1)
!       -------------------------------------------------------------
     enddo
   enddo

!  call writeMatrix('pre-tauinv', temp, dsize, dsize, TEN2m8)

   tau_inv = SROMedium(n)%tau_cpa(1:dsize, 1:dsize, 1)
!  call writeMatrix('tau_cpa11', tau_inv, dsize, dsize, TEN2m8)
!  -----------------------------------------------------------
   call MtxInv_LU(tau_inv, dsize)
!  -----------------------------------------------------------
   call zgemm('N', 'n', dsize, dsize, dsize, CONE, temp, &
        dsize, tau_inv, dsize, CZERO, temp2, dsize)
!  -----------------------------------------------------------
   call zgemm('N', 'n', dsize, dsize, dsize, CONE, tau_inv, &
        dsize, temp2, dsize, CZERO, total_proj, dsize)
!  -----------------------------------------------------------
!  call writeMatrix('new_tcpa_inv', total_proj, dsize, dsize, TEN2m8)'


   end function calculateNewTCPA
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSpeciesTauMatrix()
!  ===================================================================
!  -------------------------------------------------------------------
!  If SCF mode is off, it will calculate tau_a for all species
!  If SCF mode in on, then it will only calculate the average and block matrices (old scheme)
!  ------------------------------------------------------------------

   use ScfDataModule, only : isSROSCF
   use WriteMatrixModule, only : writeMatrix
   use SSSolverModule, only : getScatteringMatrix
!
   integer(kind=IntKind) :: j, is, ic, dsize, nsize

   dsize = SROMedium(1)%blk_size
   nsize = SROMedium(1)%neigh_size

   do j = 1, LocalNumSites
     do ic = 1, SROMedium(j)%num_species
       do is = 1, nSpinCant
        SROMedium(j)%SROTMatrix(ic)%tmat_s(is)%tmat => getScatteringMatrix('T-Matrix', spin=is, site=j, atom=ic)
        SROMedium(j)%SROTMatrix(ic)%tmat_s(is)%tmat_inv => getScatteringMatrix('TInv-Matrix', spin=is, site=j, atom=ic)
!       call writeMatrix('TmatInv', SROMedium(j)%SROTMatrix(ic)%tmat_s(is)%tmat_inv, dsize, dsize, TEN2m8)
       enddo
     enddo

     do ic = 1, SROMedium(j)%num_species
       call averageSROMatrix(j, ic)
!      call writeMatrix('TmatTildeInv', SROMedium(j)%SROTMatrix(ic)%tmat_s(1)%tmat_tilde_inv, dsize, dsize, TEN2m8)
     enddo

     do ic = 1, SROMedium(j)%num_species
       call generateBigTAMatrix(j, ic)
     enddo
   enddo

   do j = 1, LocalNumSites
     do ic = 1, SROMedium(j)%num_species
!      -------------------------------------------
       call calculateImpurityMatrix(j, ic)
!      -------------------------------------------
     enddo
   enddo

!  if (isSROSCF() == 0) then
!    do ic = 1, SROMedium(1)%num_species
!      --------------------------------------------------------------------------------
!      call writeMatrix('tau_a', SROMedium(1)%SROTMatrix(ic)%tau_ab, &
!             dsize*nsize, dsize*nsize, TEN2m8)
!      --------------------------------------------------------------------------------
!    enddo
!  endif
   
   end subroutine calSpeciesTauMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKauFromTau (e, n, ic)  result(kau_a)
!  ===================================================================

   use SSSolverModule, only : getScatteringMatrix
   use MatrixModule, only : computeUAU
   use WriteMatrixModule, only : writeMatrix
!
   complex(kind=CmplxKind), intent(in) :: e
   integer(kind=IntKind), intent(in) :: n, ic
!
   integer(kind=IntKind) :: is, js, dsize, nsize
   complex(kind=CmplxKind) :: kappa
   complex(kind=CmplxKind), pointer :: Jinv(:,:), Tinv(:,:), SinvL(:,:)
   complex(kind=CmplxKind), pointer :: Jost(:,:), OH(:,:), SinvR(:,:)
!
   complex(kind=CmplxKind), pointer :: kau(:,:), mat(:,:)
   complex(kind=CmplxKind), pointer :: kau_a(:,:,:)
!  
   dsize = SROMedium(n)%blk_size
   nsize = SROMedium(n)%neigh_size
   kappa = sqrt(e)
!
   do is = 1, nSpinCant 
      Jinv => getScatteringMatrix('JostInv-Matrix',spin=is,site=n,atom=ic)
      Tinv => getScatteringMatrix('TInv-Matrix',spin=is,site=n,atom=ic)
      SinvL => aliasArray2_c(WORK0_sro,dsize,dsize)
!     ==========================================================
!     S^{-1} = Jost^{-1}*tmat_a^{-1}/kappa
!     ----------------------------------------------------------
      call zgemm( 'n', 'n', dsize, dsize, dsize, CONE/kappa,    &
                 Jinv, dsize, Tinv, dsize, CZERO, SinvL, dsize)
!     ----------------------------------------------------------
      Jost => getScatteringMatrix('Jost-Matrix',spin=is,site=n,atom=ic)
      OH => getScatteringMatrix('OmegaHat-Matrix',spin=is,site=n,atom=ic)
      SinvR => aliasArray2_c(WORK1_sro,dsize,dsize)
!     =======================================================
!     OmegaHat = S^{-1} * tmat_a * S^{-T*}/kappa
!     S^{-T*} = Jost*OmegaHat
!     -------------------------------------------------------
      call zgemm( 'n', 'n', dsize, dsize, dsize, CONE,       &
                  Jost, dsize, OH, dsize, CZERO, SinvR, dsize)
!     -------------------------------------------------------
      kau => SROMedium(n)%SROTMatrix(ic)%kau11(:,:,is)
      mat => SROMedium(n)%SROTMatrix(ic)%tau_ab(1:dsize, 1:dsize, is)
!     =======================================================
!     kau_a = energy * S^{-1} * tau_a * S^{-T*}
!     -------------------------------------------------------
      call computeUAU(SinvL,dsize,dsize,SinvR,dsize,e,   &
           SROMedium(n)%SROTMatrix(ic)%tau_ab(1:dsize, 1:dsize, is),  &
           dsize,CZERO,kau,dsize,WORK2_sro)
!     -------------------------------------------------------
      SROMedium(n)%SROTMatrix(ic)%kau11(:,:,is) = &
          SROMedium(n)%SROTMatrix(ic)%kau11(:,:,is) - kappa*OH
   enddo
 
   kau_a => SROMedium(n)%SROTMatrix(ic)%kau11
!  call writeMatrix('kau_a11', kau_a(:,:,1), dsize, dsize, TEN2m8)
   
   end function getKauFromTau
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSROMatrix(sm_type,n,ic,is,matsize)   result(sro_mat)
!  ===================================================================
   implicit none

   character (len=*), intent(in) :: sm_type
   integer (kind=IntKind), intent(in) :: n, ic, is
   integer (kind=IntKind), intent(out), optional :: matsize
   integer (kind=IntKind) :: dsize, nsize
   logical :: is_size = .false.
   complex (kind=CmplxKind), pointer :: sro_mat(:,:)

   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface

   dsize = SROMedium(n)%blk_size
   nsize = SROMedium(n)%neigh_size

   if (present(matsize)) then
     is_size = .true.
   else
     is_size = .false.
   endif

   if (nocaseCompare(sm_type,'tau11')) then
     if (ic == 0) then
       sro_mat => SROMedium(n)%tau_cpa(1:dsize, 1:dsize, is)
     else 
       sro_mat => SROMedium(n)%SROTMatrix(ic)%tau_ab(1:dsize, 1:dsize, is)
     endif
     if (is_size) then
       matsize = dsize
     endif
   else if (nocaseCompare(sm_type,'neg-tau11')) then
     if (ic == 0) then
       sro_mat => SROMedium(n)%tau_cpac(1:dsize, 1:dsize,is)
     else
       sro_mat => SROMedium(n)%SROTMatrix(ic)%tau_abc(1:dsize, 1:dsize, is)
     endif
     if (is_size) then
       matsize = dsize
     endif
   else if (nocaseCompare(sm_type,'blk-tau')) then
     if (ic == 0) then
       sro_mat => SROMedium(n)%tau_cpa(:,:,is)
     else
       sro_mat => SROMedium(n)%SROTMatrix(ic)%tau_ab(:,:,is)
     endif
     if (is_size) then
       matsize = nsize
     endif
   else if (nocaseCompare(sm_type,'neg-blk-tau')) then
     if (ic == 0) then
       sro_mat => SROMedium(n)%tau_cpac(:,:,is)
     else
       sro_mat => SROMedium(n)%SROTMatrix(ic)%tau_abc(:,:,is)
     endif
     if (is_size) then
       matsize = nsize
     endif
   else if (nocaseCompare(sm_type,'Dmat')) then
     sro_mat => SROMedium(n)%SROTMatrix(ic)%D_ab(:,:,is)
     if (is_size) then
       matsize = nsize
     endif
   else if (nocaseCompare(sm_type,'Dtmat')) then
     sro_mat => SROMedium(n)%SROTMatrix(ic)%Dt_ab(:,:,is)
     if (is_size) then
       matsize = nsize
     endif
   else if (nocaseCompare(sm_type,'neg-Dmat')) then
     sro_mat => SROMedium(n)%SROTMatrix(ic)%D_abc(:,:,is)
     if (is_size) then
       matsize = nsize
     endif
   else if (nocaseCompare(sm_type,'neg-Dtmat')) then
     sro_mat => SROMedium(n)%SROTMatrix(ic)%Dt_abc(:,:,is)
     if (is_size) then
       matsize = nsize
     endif
   else if (nocaseCompare(sm_type,'tcpa-inv')) then
     sro_mat => SROMedium(n)%Tcpa_inv
     if (is_size) then
       matsize = dsize
     endif
   else if (nocaseCompare(sm_type,'blk-tinv')) then
     if (ic == 0) then
       sro_mat => SROMedium(n)%T_CPA_inv
     else
       sro_mat => SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%T_inv
     endif
     if (is_size) then
       matsize = nsize
     endif
   else
     call ErrorHandler('getSROMatrix', 'incorrect control string')
   endif
   
   end function getSROMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDoubleSpeciesTauMatrix(n, is, ic, ic1, neg) result(sro_mat)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: n, is, ic, ic1, neg
   complex (kind=CmplxKind), pointer :: sro_mat(:,:)   

   if (neg == 0) then
     sro_mat => SROMedium(n)%SROTMatrix(ic)%tau_sigma(ic1,:,:,is)
   else if (neg == 1) then
     sro_mat => SROMedium(n)%SROTMatrix(ic)%tau_sigmac(ic1,:,:,is)
   endif

   end function getDoubleSpeciesTauMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSROParam(n, ic1, ic2) result(w12)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: n, ic1, ic2
   real (kind=RealKind) :: w12

   w12 = SROMedium(n)%SROTMatrix(ic1)%sro_param_a(ic2)

   end function getSROParam
!  =================================================================== 

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNeighSize(n) result(neigh_size)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: neigh_size

   neigh_size = SROMedium(n)%neigh_size

   end function getNeighSize
!  ===================================================================
end module SROModule
