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
          populateTau,               &
          assembleTauFromBlocks,     &
          calculateImpurityMatrix,   &
          calSpeciesTauMatrix,       &
          getKauFromTau,             &
!

private
   integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
   integer (kind=IntKind) :: nSpinCant, nSpinPola
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind) :: kmax_kkr_max
   integer (kind=IntKind) :: ndim_Tmat
   integer (kind=IntKind) :: print_instruction
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
      complex  (kind=CmplxKind), pointer :: T_inv(:,:)
   end type TmatBlockStruct
!
   type SROTMatrixStruct
      real (kind=RealKind), pointer :: sro_param_a(:)
      type (TmatBlockStruct), allocatable :: tmat_s(:)
      complex (kind=CmplxKind), pointer :: tau_ab(:,:,:)
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
      type(TauBlockStruct), pointer :: tau_c(:)
   end type SROMediumStruct
!
   type(SROMediumStruct), allocatable :: SROMedium(:)
!  integer (kind=IntKind) :: num     ! Number of Species in CPA Medium
   logical :: test_CPA = .false.
   logical :: test_pure = .false.


contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSROMatrix (cant, pola)
!  ===================================================================
   
   use MediumHostModule, only  : getNumSites, getLocalNumSites, getGlobalSiteIndex, getNumSpecies
   use ScfDataModule, only : retrieveSROParams
   use NeighborModule, only : getNeighbor
   use SSSolverModule, only : getScatteringMatrix
   use SystemModule, only : getAtomPosition
   use WriteMatrixModule, only : writeMatrix
   use MatrixInverseModule, only : MtxInv_LU

   integer(kind=IntKind), intent(in) :: cant, pola
   integer(kind=IntKind) :: sro_param_nums, num, il, ic, is, ig, i, j, iter1, iter2, temp
   integer(kind=IntKind) :: in, jn
   real(kind=RealKind), allocatable :: sro_params(:)
!  --------------------------------------------------------
   call retrieveSROParams(sro_params, sro_param_nums)
!  --------------------------------------------------------
   GlobalNumSites = getNumSites()
   LocalNumSites = getLocalNumSites()
   if (.not. allocated(SROMedium)) then
      allocate(SROMedium(LocalNumSites))
   endif
   nSpinCant = cant
   nSpinPola = pola
   Print *,nSpinPola
 
   do il = 1, LocalNumSites
      SROMedium(il)%Neighbor => getNeighbor(il)
      SROMedium(il)%local_index = il
      ig = getGlobalSiteIndex(il)
      SROMedium(il)%CPosition = getAtomPosition(ig)
      SROMedium(il)%global_index = ig
      num = getNumSpecies(ig)
      SROMedium(il)%num_species = num
      SROMedium(il)%neigh_size = SROMedium(il)%Neighbor%NumAtoms+1
      allocate(SROMedium(il)%tau_c((SROMedium(il)%neigh_size)**2))
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
            
            do j = 1, num
              if (j < i) then
                 SROMedium(il)%SROTMatrix(i)%sro_param_a(j) = SROMedium(il)%SROTMatrix(j)%sro_param_a(i)
              else
                 
                 temp = (i - 1)*n - (i - 1)*(i - 2)/2   !j + (((i - 1)*(2*num - 2 - i))/2)
                 SROMedium(il)%SROTMatrix(i)%sro_param_a(j) = sro_params(temp + j - i + 1)
              endif
            enddo

            tm => getScatteringMatrix('T-Matrix',spin=1,site=SROMedium(il)%local_index,atom=i)
            
            SROMedium(il)%blk_size = size(tm, 1)

            if (nSpinCant == 2) then
               call ErrorHandler('initSROMatrix','SRO is not equipped to deal with spin canting yet')
            endif

            allocate(SROMedium(il)%SROTMatrix(i)%kau11(SROMedium(il)%blk_size, SROMedium(il)%blk_size, nSpinCant**2))
            allocate(SROMedium(il)%SROTMatrix(i)%tau_ab(SROMedium(il)%blk_size*SROMedium(il)%neigh_size,  &
                     SROMedium(il)%blk_size*SROMedium(il)%neigh_size, nSpinCant**2))
 
            do is = 1,nSpinCant**2
               allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%tmat_tilde_inv(SROMedium(il)%blk_size, &
                       SROMedium(il)%blk_size))
               allocate(SROMedium(il)%SROTMatrix(i)%tmat_s(is)%T_inv(SROMedium(il)%blk_size*SROMedium(il)%neigh_size, &
                       SROMedium(il)%blk_size*SROMedium(il)%neigh_size))
            enddo

         enddo
         allocate(SROMedium(il)%Tcpa(SROMedium(il)%blk_size, SROMedium(il)%blk_size))
         allocate(SROMedium(il)%Tcpa_inv(SROMedium(il)%blk_size, SROMedium(il)%blk_size))
         allocate(SROMedium(il)%T_CPA(SROMedium(il)%blk_size*SROMedium(il)%neigh_size,  &
                          SROMedium(il)%blk_size*SROMedium(il)%neigh_size))
         allocate(SROMedium(il)%T_CPA_inv(SROMedium(il)%blk_size*SROMedium(il)%neigh_size,   &
                          SROMedium(il)%blk_size*SROMedium(il)%neigh_size))
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

   allocate(z(SROMedium(1)%blk_size*SROMedium(1)%neigh_size, SROMedium(1)%blk_size*SROMedium(1)%neigh_size))
   allocate(y(SROMedium(1)%blk_size*SROMedium(1)%neigh_size, SROMedium(1)%blk_size*SROMedium(1)%neigh_size))
   
   allocate(WORK0_sro(SROMedium(1)%blk_size**2))
   allocate(WORK1_sro(SROMedium(1)%blk_size**2))
   allocate(WORK2_sro(SROMedium(1)%blk_size**2))
    
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
   real (kind=RealKind) :: wab
   
   nsize = SROMedium(n)%blk_size
 
   do is = 1, nSpinCant**2
     SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv = CZERO
     do ic = 1, SROMedium(n)%num_species
        wab = SROMedium(n)%SROTMatrix(ia)%sro_param_a(ic)
!       -------------------------------------------------------------------------------------
        call zaxpy(nsize*nsize,wab,SROMedium(n)%SROTMatrix(ic)%tmat_s(is)%tmat,1,&
                  SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv, 1)
!       -------------------------------------------------------------------------------------
     enddo

!    Inverting to get tilde_inv
!    ------------------------------------------------------------------------------
     call MtxInv_LU(SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_tilde_inv, nsize)
!    ------------------------------------------------------------------------------
   enddo   
      
   end subroutine averageSROMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generateBigTAMatrix (n, ia)
!  ===================================================================
   
   integer (kind=IntKind), intent(in) :: n, ia
   integer (kind=IntKind) :: nsize, delta, total_size, i, is,iter1,iter2, tmp

   delta = SROMedium(n)%neigh_size
   nsize = SROMedium(n)%blk_size
   total_size = nsize*delta
   
   do is = 1, nSpinCant
      SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv = CZERO

      do iter2 = 1, nsize
         do iter1 = 1, nsize
            SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1,iter2) =  &
              SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%tmat_inv(iter1,iter2)
         enddo
      enddo

      do i = 2, delta
         tmp = (i - 1)*nsize
         do iter2 = 1, nsize
            do iter1 = 1, nsize
              if (test_CPA) then
                 SROMedium(n)%SROTMatrix(ia)%tmat_s(is)%T_inv(iter1+tmp,iter2+tmp) = &
                  SROMedium(n)%Tcpa_inv(iter1, iter2)
               
              else if (test_pure) then
                 if (i < 6) then
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
            enddo
         enddo
      enddo
   enddo
   

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
!   
   integer(kind=IntKind), intent(in) :: n, ic
!
   integer(kind=IntKind) :: dsize, nsize, is, it, k

   SROMedium(n)%SROTMatrix(ic)%tau_ab = CZERO 

   dsize = SROMedium(n)%blk_size
   nsize = SROMedium(n)%neigh_size
   y = CZERO

   z = SROMedium(n)%SROTMatrix(ic)%tmat_s(1)%T_inv - SROMedium(n)%T_CPA_inv

   call computeAprojB('L', dsize*nsize, SROMedium(n)%tau_cpa(:,:, 1), z, y)

   SROMedium(n)%SROTMatrix(ic)%tau_ab(:, :, 1) = y

   end subroutine calculateImpurityMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSpeciesTauMatrix()
!  ===================================================================

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
       enddo
     enddo

     do ic = 1, SROMedium(j)%num_species
       call averageSROMatrix(j, ic)
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
   
   end function getKauFromTau
!  ===================================================================
end module SROModule
