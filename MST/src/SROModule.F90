module SROModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1, TEN2m6, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler

!
public :: initSROMatrix,            &
          averageSROMatrix,             &
!

private
   integer (kind=IntKind) :: GlobalNumSites, LocalNumSites
   integer (kind=IntKind) :: nSpinCant
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID
   integer (kind=IntKind) :: kmax_kkr_max
   integer (kind=IntKind) :: ndim_Tmat
   integer (kind=IntKind) :: print_instruction
!
   type SROTMatrixStruct
      real (kind=RealKind), pointer :: sro_param_a(:)
      complex (kind=CmplxKind), pointer :: tmat_a(:,:)
      complex (kind=CmplxKind), pointer :: tmat_a_inv(:,:)
      complex (kind=CmplxKind), pointer :: tmat_tilde_a(:,:) ! In global spin framework
      complex (kind=CmplxKind), pointer :: tmat_tilde_a_inv(:, :) ! In global spin framework
      complex (kind=CmplxKind), pointer :: T_a(:,:) ! In global spin framework
      complex (kind=CmplxKind), pointer :: T_a_inv(:,:)
   end type SROTMatrixStruct
!
   type SROMediumStruct
      integer (kind=IntKind) :: local_index
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: neighbors
      integer (kind=IntKind) :: num_species
      integer (kind=IntKind) :: blk_size
      logical :: isCPA
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
   subroutine initSROMatrix (cant)
!  ===================================================================
   
   use MediumHostModule, only  : getNumSites, getLocalNumSites, getGlobalSiteIndex, getNumSpecies
   use ScfDataModule, only : retrieveSROParams
   use NeighborModule, only : getNumNeighbors
   use SSSolverModule, only : getScatteringMatrix
   use MatrixInverseModule, only : MtxInv_LU

   integer(kind=IntKind), intent(in) :: cant
   integer(kind=IntKind) :: sro_param_nums, num, il, ic, ig, i, j, temp
   real(kind=RealKind), allocatable :: sro_params(:)
!  --------------------------------------------------------
   call retrieveSROParams(sro_params, sro_param_nums)
!  --------------------------------------------------------
   GlobalNumSites = getNumSites()
   LocalNumSites = getLocalNumSites()
   if (.not. allocated(SROMedium)) then
      allocate(SROMedium(LocalNumSites))
   endif
   Print *, LocalNumSites
 
   do il = 1, LocalNumSites
      SROMedium(il)%neighbors = getNumNeighbors(il)
      SROMedium(il)%local_index = il
      ig = getGlobalSiteIndex(il)
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
            SROMedium(il)%SROTMatrix(i)%tmat_a => getScatteringMatrix('T-Matrix',spin=1,site=SROMedium(il)%local_index,atom=i)
            SROMedium(il)%SROTMatrix(i)%tmat_a_inv => getScatteringMatrix('T-Matrix',spin=1,site=SROMedium(il)%local_index,atom=i)
         !  -------------------------------------------------------------------
            call MtxInv_LU(SROMedium(il)%SROTMatrix(i)%tmat_a_inv,size(SROMedium(il)%SROTMatrix(i)%tmat_a_inv, 1))
         !  -------------------------------------------------------------------
         enddo
       endif
      SROMedium(il)%blk_size = size(SROMedium(il)%SROTMatrix(num)%tmat_a, 1)       
   enddo
  
   Print *,"Successful init" 
   
   nSpinCant = cant
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
   
   complex (kind=CmplxKind), pointer :: tm0(:,:)
   integer (kind=IntKind) :: ic, nsize
   real (kind=RealKind) :: wab
   
   nsize = SROMedium(n)%blk_size
!  Print *,nsize
   allocate(SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a(nsize, nsize))
   allocate(SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a_inv(nsize, nsize))

   do ic = 1, SROMedium(n)%num_species
     wab = SROMedium(n)%SROTMatrix(ia)%sro_param_a(ic)
!    -------------------------------------------------------------------------------------
     call zaxpy(nsize*nsize,wab,SROMedium(n)%SROTMatrix(ic)%tmat_a,1,SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a,1)
     call zaxpy(nsize*nsize,wab,SROMedium(n)%SROTMatrix(ic)%tmat_a_inv,1,SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a_inv, 1)
!    -------------------------------------------------------------------------------------
   enddo
      
   end subroutine averageSROMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generateBigTAMatrix (n, ia)
!  ===================================================================
   
   integer (kind=IntKind) :: nsize, delta, total_size, i, tmp
   complex (kind=CmplxKind), pointer :: tm0(:,:), tm1(:,:)

   delta = SROMedium(n)%neighbors + 1
   nsize = SROMedium(n)%blk_size
   total_size = nsize*delta
!  Print *, delta, nsize, total_size
   allocate(SROMedium(n)%SROTMatrix(ia)%T_a(total_size, total_size))
   allocate(SROMedium(n)%SROTMatrix(ia)%T_a_inv(total_size, total_size))
   
   do i = 1, delta
      tmp = (i - 1)*blk_size
      if (i == 1) then
         tm1 => SROMedium(n)%SROTMatrix(ia)%tmat_a
         tm0 => SROMedium(n)%SROTMatrix(ia)%tmat_a_inv
      else
         tm1 => SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a
         tm0 => SROMedium(n)%SROTMatrix(ia)%tmat_tilde_a_inv
      endif
      SROMedium(n)%SROTMatrix(ia)%T_a(1+tmp:nsize+tmp,1+tmp:nsize+tmp) = tm1
      SROMedium(n)%SROTMatrix(ia)%T_a_inv(1+tmp:nsize+tmp,1+tmp:nsize+tmp) = tm0
   enddo
   

   end subroutine generateBigTAMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generateBigTCPAMatrix (n)
!  ===================================================================

   use CPAMediumModule, only : getCPAMatrix
   use MatrixInverseModule, only : MtxInv_LU
   
   integer (kind=IntKind) :: i, tmp, nsize, delta, total_size
   
   SROMedium(n)%Tcpa => getCPAMatrix('Tcpa', n)
   SROMedium(n)%Tcpa_inv => getCPAMatrix('Tcpa', n)
   nsize = SROMedium(n)%blk_size
   delta = SROMedium(n)%neighbors + 1
   total_size = nsize*delta
   Print *, "Okay till here!" 
!  -------------------------------------------------------------------
   call MtxInv_LU(SROMedium(n)%Tcpa_inv,nsize)
!  -------------------------------------------------------------------
   Print *, "If you're reading this, everything's fine!"
   allocate(SROMedium(n)%T_CPA(total_size, total_size))
   allocate(SROMedium(n)%T_CPA_inv(total_size, total_size))
   
   do i = 1, delta
       tmp = (i - 1)*nsize
       SROMedium(n)%T_CPA(1+tmp:nsize+tmp,1+tmp:nsize+tmp) = SROMedium(n)%Tcpa
       SROMedium(n)%T_CPA_inv(1+tmp:nsize+tmp,1+tmp:nsize+tmp) = SROMedium(n)%Tcpa_inv
   enddo
   
   end subroutine generateBigTCPAMatrix
!  ===================================================================
end module SROModule
