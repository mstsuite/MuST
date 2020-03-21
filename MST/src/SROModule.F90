module SROModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
  use PublicTypeDefinitionsModule, only : NeighborStruct

!
public :: initSROMatrix,            &
          averageSROMatrix,             &
!

private
   integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
   integer (kind=IntKind) :: nSpinCant, nSpinPola
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind) :: kmax_kkr_max
   integer (kind=IntKind) :: ndim_Tmat
   integer (kind=IntKind) :: print_instruction
!
   type SROTMatrixStruct
      real (kind=RealKind), pointer :: sro_param_a(:)
      complex (kind=CmplxKind), pointer :: tmat_a(:,:,:)
      complex (kind=CmplxKind), pointer :: tmat_a_inv(:,:,:)
      complex (kind=CmplxKind), pointer :: tmat_tilde_a(:,:,:) ! In global spin framework
      complex (kind=CmplxKind), pointer :: tmat_tilde_a_inv(:,:,:) ! In global spin framework
      complex (kind=CmplxKind), pointer :: T_a(:,:,:) ! In global spin framework
      complex (kind=CmplxKind), pointer :: T_a_inv(:,:,:)
   end type SROTMatrixStruct
!
   type SROMediumStruct
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: num_species
      integer (kind=IntKind) :: blk_size
      integer (kind=IntKind) :: CPosition(3)
      logical :: isCPA
      type(NeighborStruct), pointer :: Neighbor
      type(SROTMatrixStruct), allocatable :: SROTMatrix(:)
      complex (kind=CmplxKind), pointer :: Tcpa(:,:)
      complex (kind=CmplxKind), pointer :: Tcpa_inv(:,:)
      complex (kind=CmplxKind), pointer :: T_CPA(:,:)
      complex (kind=CmplxKind), pointer :: T_CPA_inv(:,:)
   end type SROMediumStruct
!
   type(SROMediumStruct), allocatable :: SROMedium(:)
!  integer (kind=IntKind) :: num     ! Number of Species in CPA Medium


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

   integer(kind=IntKind), intent(in) :: cant, pola
   integer(kind=IntKind) :: sro_param_nums, num, il, ic, ig, i, j, iter1, iter2, temp
   real(kind=RealKind), allocatable :: sro_params(:)
   complex(kind=CmplxKind), pointer :: tm0(:,:), tm1(:,:)
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
            allocate(SROMedium(il)%SROTMatrix(i)%sro_param_a(num))
            do j = 1, num
               if (j < i) then
                  SROMedium(il)%SROTMatrix(i)%sro_param_a(j) = SROMedium(il)%SROTMatrix(j)%sro_param_a(i)
               else
                  temp = j + (((i - 1)*(2*num - 2 - i))/2)
                  SROMedium(il)%SROTMatrix(i)%sro_param_a(j) = sro_params(temp)
               endif
            enddo

            tm0 => getScatteringMatrix('T-Matrix',spin=1,site=SROMedium(il)%local_index,atom=i)
            SROMedium(il)%blk_size = size(tm0, 1)
            allocate(SROMedium(il)%SROTMatrix(i)%tmat_a(SROMedium(il)%blk_size,SROMedium(il)%blk_size, nSpinPola))
            allocate(SROMedium(il)%SROTMatrix(i)%tmat_a_inv(SROMedium(il)%blk_size,SROMedium(il)%blk_size, nSpinPola))
            if (nSpinCant == 2) then
               call ErrorHandler('initSROMatrix','SRO is not equipped to deal with spin canting yet')
            endif
            do ic = 1,nSpinPola
               tm0 => getScatteringMatrix('T-Matrix',spin=ic,site=SROMedium(il)%local_index,atom=i)
               tm1 => getScatteringMatrix('TInv-Matrix',spin=ic,site=SROMedium(il)%local_index,atom=i)
               do iter2 = 1,SROMedium(il)%blk_size
                  do iter1 = 1,SROMedium(il)%blk_size
                     SROMedium(il)%SROTMatrix(i)%tmat_a(iter1,iter2,ic) = tm0(iter1, iter2)
                     SROMedium(il)%SROTMatrix(i)%tmat_a_inv(iter1,iter2,ic) = tm1(iter1, iter2)
                  enddo
               enddo
            enddo
         enddo
       endif
   enddo
  
   Print *,"Successful init" 
   
!  do i = 1, LocalNumSites
!     do j = 1, SROMedium(i)%num_species
!        Print *, SROMedium(i)%SROTMatrix(j)%sro_param_a
!        Print *, SROMedium(i)%nearest_neighbors
!     enddo
!  enddo
   
   call averageSROMatrix (1, 1)
!  call generateBigTAMatrix(1, 1)
!  call generateBigTCPAMatrix(1)
   
   end subroutine initSROMatrix
!  =================================================================== 

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageSROMatrix (n, ia)
!  ===================================================================

!  NEED TO CHECK IF SITE N IS CPA OR NOT
!  BE CAREFUL WHETHER TO USE LOCAL OR GLOBAL INDEX ANYWHERE

   complex (kind=CmplxKind), pointer :: tm0(:,:)
   integer (kind=IntKind) :: ic, is, nsize
   real (kind=RealKind) :: wab
   
   nsize = SROMedium(n)%blk_size
!  Print *,nsize
   allocate(SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a(nsize, nsize, nSpinPola))
   allocate(SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a_inv(nsize, nsize, nSpinPola))

   do ic = 1, SROMedium(n)%num_species
     wab = SROMedium(n)%SROTMatrix(ia)%sro_param_a(ic)
     do is = 1, nSpinPola
!       -------------------------------------------------------------------------------------
        call zaxpy(nsize*nsize,wab,SROMedium(n)%SROTMatrix(ic)%tmat_a(:,:,is),1,&
                  SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a(:,:,is),1)
        call zaxpy(nsize*nsize,wab,SROMedium(n)%SROTMatrix(ic)%tmat_a_inv(:,:,is),1,&
                  SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a_inv(:,:,is), 1)
!       -------------------------------------------------------------------------------------
     enddo
   enddo
   
!  do ic = 1,SROMedium(n)%num_species
!    Print *,"Species is",ic
!    Print *,SROMedium(n)%SROTMatrix(ic)%tmat_a
!  enddo
!  Print *,"The average is"
!  Print *,SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a
      
   end subroutine averageSROMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generateBigTAMatrix (n, ia)
!  ===================================================================
   
   integer (kind=IntKind) :: nsize, delta, total_size, i, is,iter1,iter2, tmp
   complex (kind=CmplxKind), pointer :: tm0(:,:), tm1(:,:)

   delta = SROMedium(n)%Neighbor%NumAtoms + 1
   nsize = SROMedium(n)%blk_size
   total_size = nsize*delta
!  Print *, delta, nsize, total_size
   allocate(SROMedium(n)%SROTMatrix(ia)%T_a(total_size, total_size, nSpinPola))
   allocate(SROMedium(n)%SROTMatrix(ia)%T_a_inv(total_size, total_size, nSpinPola))
   
   do is = 1, nSpinPola
      tm1 => SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a(:,:,is)
      tm0 => SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a_inv(:,:,is)
      SROMedium(n)%SROTMatrix(ia)%T_a(:,:,is) = CZERO
      SROMedium(n)%SROTMatrix(ia)%T_a_inv(:,:,is) = CZERO
      do i = 1, delta
         tmp = (i - 1)*blk_size
         if (i == 1) then
            tm1 => SROMedium(n)%SROTMatrix(ia)%tmat_a(:,:,is)
            tm0 => SROMedium(n)%SROTMatrix(ia)%tmat_a_inv(:,:,is)
         endif
         do iter2 = 1+tmp,nsize+tmp
            do iter1 = 1+tmp,nsize+tmp
              SROMedium(n)%SROTMatrix(ia)%T_a(iter1,iter2,is) = tm1(iter1, iter2)
              SROMedium(n)%SROTMatrix(ia)%T_a_inv(iter1,iter2,is) = tm0(iter1, iter2)
            enddo
         enddo
      enddo
   enddo
   

   end subroutine generateBigTAMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generateBigTCPAMatrix (n)
!  ===================================================================

   use CPAMediumModule, only : getCPAMatrix
   use MatrixInverseModule, only : MtxInv_LU
   
   integer (kind=IntKind) :: i, iter1, iter2, tmp, nsize, delta, total_size
   
   SROMedium(n)%Tcpa => getCPAMatrix('Tcpa', n)
   SROMedium(n)%Tcpa_inv => getCPAMatrix('Tcpa', n)
   nsize = SROMedium(n)%blk_size
   delta = SROMedium(n)%Neighbor%NumAtoms + 1
   total_size = nsize*delta
!  -------------------------------------------------------------------
   call MtxInv_LU(SROMedium(n)%Tcpa_inv,nsize)
!  -------------------------------------------------------------------
   allocate(SROMedium(n)%T_CPA(total_size, total_size))
   allocate(SROMedium(n)%T_CPA_inv(total_size, total_size))
   SROMedium(n)%T_CPA = CZERO
   SROMedium(n)%T_CPA_inv = CZERO
    
   do i = 1, delta
       tmp = (i - 1)*nsize
       do iter2 = 1+tmp,nsize+tmp
          do iter1 = 1+tmp,nsize+tmp
             SROMedium(n)%T_CPA(iter1,iter2) = SROMedium(n)%Tcpa(iter1,iter2)
             SROMedium(n)%T_CPA_inv(iter1,iter2) = SROMedium(n)%Tcpa_inv(iter1,iter2)
          enddo
       enddo
   enddo
   
   end subroutine generateBigTCPAMatrix
!  ===================================================================
end module SROModule
