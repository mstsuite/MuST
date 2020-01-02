!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  History: Original code written by D.D. Johnson (see PRB 38, 12807)
!                  note:  there are a few typos in that paper but 
!                  the code is working!
!              Rewritten by W. A. Shelton for LSMS code 6/21/94
!                  this version is easy to read (no goto!!!! more comments ...)
!                  and is setup for MPP machines (not tested)
!              Rewritten by T. C. Schulthess, ORNL, March 97
!                  this version should work for any code (see comments below)
!
!     Bug fixes:   TCS, 8/5/97 see comments below 
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     further comments on how to use this subroutine:
!     (if nothing useful stands here I had no time yet to update these
!     comments, please consult usage in lkkr code version 0.6.3 or later,
!     or call Thomas Schulthess (423) 5768067)
!
!      vector(r,i) -> i=1: old vector (input), scratch (ouput)
!                  -> i=2: new vector (input), mixed vector (output)
!      vlen        -> length of vector
!      alpha       -> linear mixing factor
!      rms         -> RMS difference between old and new vector
!      iter        -> iteration number (if 1, linear mixing, broyden reset)
!      broylen     -> number of iterations that are used for mixing
!                     (<=MaxNumBroydenIter)
!      u, vt, f, df, w and vold are working arrays that need to be saved
!                    between call to this subroutine
!      a, b, d, cm, and ipiv are working arrays that need not be saved
!      MaxNumBroydenIter    -> maximum number of iterations that can be saved
!      MaxSizeBroydenVect       -> maximum length of vectors
!
!      See declaration for exact dimentsions and types
!
!      There are two options for matrix inversions, a Gaussian
!      elimination routine called invert1 and calls to lapack routines
!      with pivoting (see comments "using invert1" and "using lapack").
!      Obviously only one can be used, comment out the other one.
!
!      When using this subroutine in a parallel code in which only 
!      parts of the vectors are known on every node, make sure that 
!      the calls to gldsum (global sum) are correct (LKKR and LSMS 
!      codes have different calls).
!
!      In a serial code, either comment out the calls to GlobalSum or
!      provide a dummy subroutine
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module MixingModule
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : zero, one, two, ten2m2, ten2m8, ten2m9, &
                               czero, cone, PI4, THIRD
   use GroupCommModule, only : GlobalSumInGroup
   use PublicTypeDefinitionsModule, only : MixListRealStruct, &
                                           MixListCmplxStruct
!
public :: initMixing,         &
          mixValues,          &
          setBroydenMixing,   &
          setSimpleMixing,    &
          setDGAMixing,       &
          setMixingAlpha,     &
          resetMixing,        &
          endMixing
!
   interface initMixing
      module procedure initMixing_r, initMixing_c
   end interface initMixing
!
   interface mixValues
      module procedure mixValues_r, mixValues_c
   end interface mixValues
!
   interface resetMixing
      module procedure resetMixing_r, resetMixing_c
   end interface resetMixing
!
private
!
!  General parameters
!
   logical :: Initialized = .false.
   logical :: InitializedWkSpace = .false.
!
   integer(kind=IntKind) :: iter_count = 0
   integer(kind=IntKind) :: NumTotalMix
   integer(kind=IntKind) :: NumTypes
   integer(kind=IntKind), allocatable :: NumQuantities(:)
   integer(kind=IntKind) :: MaxVlen
   integer(kind=IntKind) :: MixingMethod = -1 ! -1 - No method have been set yet
                                              !  0 - Simple
                                              !  1 - Broyden 
                                              !  2 - D.G. Anderson
!
!  Simple Mixing Parameters
!
   type SimpleStruct
      logical :: Initialized
!
      integer (kind=IntKind) :: vlen
!
      real(kind=RealKind) :: alpha
!
      real(kind=RealKind), pointer :: vector_old_r(:)
      real(kind=RealKind), pointer :: vector_new_r(:)
!
      complex(kind=CmplxKind), pointer :: vector_old_c(:)
      complex(kind=CmplxKind), pointer :: vector_new_c(:)
   end type SimpleStruct
!
   type (SimpleStruct), allocatable :: SimpleMix(:)
!
!  Broyden Mixing Parameters
!
   integer(kind=IntKind) :: MaxNumBroydenIter = 10
   integer(kind=IntKind) :: MaxSizeBroydenVect
   integer(kind=IntKind) :: BroydenInvMethod = 1
!
   integer(kind=IntKind) :: lastIter = 0
!
   integer(kind=IntKind) :: NumBroydenIter
   integer(kind=IntKind) :: SizeBroydenVect
!
   type BroydenStruct
      logical :: Initialized
!
      integer(kind=IntKind) :: vlen
!
      real(kind=RealKind) :: alpha
!
      real(kind=RealKind) :: rms
      real(kind=RealKind) :: weight
!
      real(kind=RealKind), pointer :: mesh(:)
!
      real(kind=RealKind), pointer :: vector_old_r(:)
      real(kind=RealKind), pointer :: vector_new_r(:)
!
      complex(kind=CmplxKind), pointer :: vector_old_c(:)
      complex(kind=CmplxKind), pointer :: vector_new_c(:)
!
      real(kind=RealKind), pointer :: u_r(:,:)
      real(kind=RealKind), pointer :: vt_r(:,:)
      real(kind=RealKind), pointer :: f_r(:)
      real(kind=RealKind), pointer :: df_r(:)
      real(kind=RealKind), pointer :: w_r(:)
      real(kind=RealKind), pointer :: vold_r(:)
!
      complex(kind=CmplxKind), pointer :: u_c(:,:)
      complex(kind=CmplxKind), pointer :: vt_c(:,:)
      complex(kind=CmplxKind), pointer :: f_c(:)
      complex(kind=CmplxKind), pointer :: df_c(:)
      complex(kind=CmplxKind), pointer :: w_c(:)
      complex(kind=CmplxKind), pointer :: vold_c(:)
   end type BroydenStruct
!
   type (BroydenStruct), allocatable :: BroydenMix(:)
!
!  D.G.Anderson mixing Parameters
!
   integer(kind=IntKind) :: MaxNumDGAIter = 10
   integer(kind=IntKind) :: NumDGAIter
   integer(kind=IntKind) :: GroupID
!
   type DGAStruct
      logical :: Initialized
!
      integer (kind=IntKind) :: vlen
!
      real(kind=RealKind) :: alpha
!
      real(kind=RealKind), pointer :: f_old_r(:,:)
      real(kind=RealKind), pointer :: f_new_r(:,:)
!
      complex(kind=CmplxKind), pointer :: f_old_c(:,:)
      complex(kind=CmplxKind), pointer :: f_new_c(:,:)
!
      real(kind=RealKind), pointer :: vector_old_r(:)
      real(kind=RealKind), pointer :: vector_new_r(:)
!
      complex(kind=CmplxKind), pointer :: vector_old_c(:)
      complex(kind=CmplxKind), pointer :: vector_new_c(:)
   end type DGAStruct
!
   type (DGAStruct), allocatable :: DGAMix(:)
!
!  Working Space
!
   integer (kind=IntKind), allocatable, target :: ipiv(:)
!
   real(kind=RealKind), allocatable, target :: a_r(:,:), b_r(:,:)
   real(kind=RealKind), allocatable, target :: d_r(:,:), cm_r(:)
!
   complex(kind=CmplxKind), allocatable, target :: a_c(:,:), b_c(:,:)
   complex(kind=CmplxKind), allocatable, target :: d_c(:,:), cm_c(:)
!
contains
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMixing_r( nt, nq, RealArrayList )
!  ===================================================================
   use GroupCommModule, only : getGroupID
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: nt, nq(nt) 
   integer(kind=IntKind) :: idt, idq
!
   type (MixListRealStruct), target :: RealArrayList
   type (MixListRealStruct), pointer :: p_RAL
!
   if ( nt<1 ) then
      call ErrorHandler('initMixing','Number of Types is less than 1')
   endif
!
   NumTypes = nt
   allocate( NumQuantities(NumTypes) )
   NumQuantities = nq
   GroupID = getGroupID('Unit Cell')
!
   NumTotalMix = 0
   p_RAL => RealArrayList
   nullify( p_RAL%next )
!
   do idt = 1,NumTypes
      if ( NumQuantities(idt)<1 ) then
         call ErrorHandler('initMixing', &
              'Number of Quntities is less than 1 sor type',idt)
      endif
      do idq = 1,NumQuantities(idt)
         NumTotalMix = NumTotalMix +1
         if ( idq < NumQuantities(idt)) then
            allocate(p_RAL%next)
            p_RAL => p_RAL%next
            nullify( p_RAL%next )
         endif
      enddo
   enddo
!   nullify( p_RAL%next )
!
   Initialized = .true.
!
   end subroutine initMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMixing_c( nt, nq, CmplxArrayList )
!  ===================================================================
   use GroupCommModule, only : getGroupID
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: nt, nq(nt)
   integer(kind=IntKind) :: idt, idq
!
   type (MixListCmplxStruct), target :: CmplxArrayList
   type (MixListCmplxStruct), pointer :: p_CAL
!
   if ( nt<1 ) then
      call ErrorHandler('initMixing','Number of Types is less than 1')
   endif
!
   NumTypes = nt
   allocate( NumQuantities(NumTypes) )
   NumQuantities = nq
   GroupID = getGroupID('Unit Cell')
!
   NumTotalMix = 0
   p_CAL => CmplxArrayList
   nullify( p_CAL%next )
!
   do idt = 1,NumTypes
      if ( NumQuantities(idt)<1 ) then
         call ErrorHandler('initMixing', &
              'Number of Quntities is less than 1 sor type',idt)
      endif
      do idq = 1,NumQuantities(idt)
         NumTotalMix = NumTotalMix +1
         if (idq < NumQuantities(idt)) then
            allocate(p_CAL%next)
            p_CAL => p_CAL%next
            nullify( p_CAL%next )
         endif
      enddo
   enddo
!   nullify( p_CAL%next )
!
   Initialized = .true.
!
   end subroutine initMixing_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMixing()
!  ===================================================================
!
   implicit none
!
   if ( Initialized ) then
      deallocate( NumQuantities )
      call delWorkingSpace()
      call delMixing()
      if ( MixingMethod == 0 ) then
         deallocate( SimpleMix )
      else if ( MixingMethod == 1 ) then
         deallocate( BroydenMix )
      else
         deallocate( DGAMix )
      endif
      Initialized = .false.
   endif
!
   MixingMethod = -1
!
   end subroutine endMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine mixValues_r( list_q )
!  ===================================================================
!
   implicit none
!
   type(MixListRealStruct) :: list_q
!
   iter_count =iter_count+1
!
   if ( MixingMethod == 0 ) then
      call calSimpleMixing_r(list_q)
   else if ( MixingMethod == 1 ) then
      call calBroydenMixing_r(list_q)
   else
      call calDGAMixing_r(list_q)
   endif
!
   end subroutine mixValues_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine mixValues_c( list_q )
!  ===================================================================
!
   implicit none
!
   type(MixListCmplxStruct) :: list_q
!
   iter_count =iter_count+1
!
   if ( MixingMethod == 0 ) then
      call calSimpleMixing_c(list_q)
   else if ( MixingMethod == 1 ) then
      call calBroydenMixing_c(list_q)
   else
      call calDGAMixing_c(list_q)
   endif
!
   end subroutine mixValues_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSimpleMixing( idt, idq, amix )
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: idt, idq
!
   real(kind=RealKind), intent(in) :: amix
!
   integer(kind=IntKind) :: id
!
   if ( .not.Initialized ) then
      call ErrorHandler('setSimpleMixing',                            &
          'The MixingModule needs to be initialized first')
   endif
!
   if ( idt<0 .or. idt>NumTypes ) then
      call ErrorHandler('setSimpleMixing',"Invalid type index",idt)
   endif
   if ( idq<0 .or. idq>NumQuantities(idt) ) then
      call ErrorHandler('setSimpleMixing',"Invalid quantity index",idq)
   endif
!
   if ( .not.allocated(SimpleMix) ) then
      allocate( SimpleMix(NumTotalMix) )
      do id = 1,NumTotalMix
         SimpleMix(id)%Initialized = .false.
      enddo 
   endif  
!
   id = getMixID(idt,idq)
   if ( .not.SimpleMix(id)%Initialized ) then
      SimpleMix(id)%alpha = amix
   endif
!
   if ( MixingMethod /=0 ) then
      MixingMethod = 0
   endif
!
   end subroutine setSimpleMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delSimpleMixing()
!  ===================================================================
   implicit none
!
   if ( allocated(SimpleMix) ) then
!      deallocate( SimpleMix )
   endif
!
   end subroutine delSimpleMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSimpleMixing_r( list_q)
!  ===================================================================
   implicit none
!
   type ( MixListRealStruct ), target :: list_q
!
   integer(kind=IntKind) :: j, idq, vlen 
!
   real(kind=RealKind) :: alpha, beta
   real(kind=RealKind), pointer :: pvect_old(:), pvect_new(:)
!
   type ( MixListRealStruct ), pointer :: plq_tmp
!
   plq_tmp => list_q
   do idq = 1,NumTotalMix
      vlen = plq_tmp%size
      if ( .not.SimpleMix(idq)%Initialized ) then
         call setSimpleSpace(idq,vlen)
      endif
!
      SimpleMix(idq)%vector_old_r => plq_tmp%vector_old(1:vlen)
      SimpleMix(idq)%vector_new_r => plq_tmp%vector_new(1:vlen)
!
      alpha = SimpleMix(idq)%alpha
      beta = ONE-alpha
      if ( beta < zero ) beta = ZERO
!
      vlen = SimpleMix(idq)%vlen
      pvect_old => SimpleMix(idq)%vector_old_r(1:vlen)
      pvect_new => SimpleMix(idq)%vector_new_r(1:vlen)
      do j = 1,vlen
         pvect_new(j) = alpha*pvect_new(j) + beta*pvect_old(j)
      enddo
      plq_tmp => plq_tmp%next
   enddo
!
   nullify(plq_tmp)
!
   end subroutine calSimpleMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSimpleMixing_c( list_q)
!  ===================================================================
   implicit none
!
   type ( MixListCmplxStruct ), target :: list_q
!
   integer(kind=IntKind) :: j, idq, vlen 
!
   real(kind=RealKind) :: alpha, beta
   complex(kind=CmplxKind), pointer :: pvect_old(:), pvect_new(:)
!
   type ( MixListCmplxStruct ), pointer :: plq_tmp
!
   plq_tmp => list_q
   do idq = 1,NumTotalMix
      vlen = plq_tmp%size
      if ( .not.SimpleMix(idq)%Initialized ) then
         call setSimpleSpace(idq,vlen)
      endif
!
      SimpleMix(idq)%vector_old_c => plq_tmp%vector_old(1:vlen)
      SimpleMix(idq)%vector_new_c => plq_tmp%vector_new(1:vlen)
!
      alpha = SimpleMix(idq)%alpha
      beta = ONE-alpha
      if ( beta < zero ) beta = ZERO
!
      vlen = SimpleMix(idq)%vlen
      pvect_old => SimpleMix(idq)%vector_old_c(1:vlen)
      pvect_new => SimpleMix(idq)%vector_new_c(1:vlen)
      do j = 1,vlen
         pvect_new(j) = alpha*pvect_new(j) + beta*pvect_old(j)
      enddo
      plq_tmp => plq_tmp%next
   enddo
!
   nullify(plq_tmp)
!
   end subroutine calSimpleMixing_c
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setDGAMixing( idt, idq, amix )
!  ==================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: idt, idq
!  integer(kind=IntKind), optional   :: ndgaiter
!
   integer(kind=IntKind) :: id
!
   real(kind=RealKind), intent(in) :: amix
!
   if ( .not.Initialized ) then
      call ErrorHandler('setDGAMix',                                 &
          'The MixingModule needs to be initialized first')
   endif
!
   if ( idt<0 .or. idt>NumTypes ) then
      call ErrorHandler('setDGAMixing',"Invalid type index",idt)
   endif
   if ( idq<0 .or. idq>NumQuantities(idt) ) then
      call ErrorHandler('setDGAMixing',"Invalid quantity index",idq)
   endif
!
   if ( .not.allocated(DGAMix) ) then
      allocate( DGAMix(NumTotalMix) )
      do id = 1,NumTotalMix
         DGAMix(id)%Initialized = .false.
      enddo
   endif
!
   NumDGAIter = MaxNumDGAIter
!  if ( present(ndgaiter) ) then
!     if ( ndgaiter > MaxNumDGAIter ) then
!        call WarningHandler('setDGAMixing', &
!            'The number of DGA iterations exeedes the maximum limit')
!     endif 
!     NumDGAIter = ndgaiter
!  endif
!
   id = getMixID(idt,idq)
   if ( .not.DGAMix(id)%Initialized ) then
      DGAMix(id)%alpha = amix
   endif
!
   if ( MixingMethod /=2 ) then
      MixingMethod = 2
   endif
!
   end subroutine setDGAMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delDGAMixing()
!  ===================================================================
   implicit none
!
   integer(kind=IntKind) :: iq
!
   if ( allocated(DGAMix) ) then
      do iq = 1,NumTotalMix
         if ( .not.DGAMix(iq)%Initialized ) cycle
         DGAMix(iq)%Initialized = .false.
         if ( associated(DGAMix(iq)%f_old_r) ) deallocate( DGAMix(iq)%f_old_r )
         if ( associated(DGAMix(iq)%f_new_r) ) deallocate( DGAMix(iq)%f_new_r )
         if ( associated(DGAMix(iq)%f_old_c) ) deallocate( DGAMix(iq)%f_old_c )
         if ( associated(DGAMix(iq)%f_new_c) ) deallocate( DGAMix(iq)%f_new_c )
      enddo 
!      deallocate( DGAMix )
   endif
!
   end subroutine delDGAMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDGAMixing_r( list_q )
!  ===================================================================
!
   implicit none
!
   type(MixListRealStruct), target :: list_q
!
   integer(kind=IntKind) :: i, ii, jj, n, iter, isz, vlen
   integer(kind=IntKind) :: idt, idq, id
!
   real(kind=RealKind) :: alpha
   real(kind=RealKind), pointer :: pvect_old(:), pvect_new(:) , pb(:)
   real(kind=RealKind), pointer :: pf_old(:,:), pf_new(:,:)
   real(kind=RealKind), pointer :: pfi_t(:), pfo_t(:)
!
   type(MixListRealStruct), pointer :: plq_tmp
!
   iter = min( iter_count, NumDGAIter )
   n    = mod( iter_count-1, NumDGAIter ) + 1
!
   if ( .not.InitializedWkSpace ) then
      call setWorkingSpace_r(MixingMethod)
   endif
!
   plq_tmp => list_q
   do idt = 1,NumTypes
      a_r = zero
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = plq_tmp%size
         if ( .not.DGAMix(id)%Initialized ) then
            call setDGASpace_r(id,vlen)
         endif
         DGAMix(id)%vector_old_r => plq_tmp%vector_old(1:vlen)
         DGAMix(id)%vector_new_r => plq_tmp%vector_new(1:vlen)
         plq_tmp => plq_tmp%next
      enddo
!
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         pvect_old => DGAMix(id)%vector_old_r(1:vlen)
         pf_old    => DGAMix(id)%f_old_r(1:NumDGAIter,1:vlen)
         pvect_new => DGAMix(id)%vector_new_r(1:vlen)
         pf_new    => DGAMix(id)%f_new_r(1:NumDGAIter,1:vlen)
         vlen = DGAMix(id)%vlen
         do i = 1,vlen
            pf_old(n,i)  = pvect_old(i)
            pf_new(n,i) = pvect_new(i)
         enddo
!
         do jj = 1,iter
            do ii = 1,jj
!              -------------------------------------------------------
               a_r(ii,jj) = a_r(ii,jj) +                               &
                          dgarms_r(pf_old,pf_new,ii,jj,vlen)
!              -------------------------------------------------------
               if ( ii.ne.jj ) then
                  a_r(jj,ii) = a_r(jj,ii) + a_r(ii,jj)
               endif
            enddo
         enddo
!
      enddo
!     ================================================================
!     sum the a(i,j) matrix over all other atoms......................
!     ================================================================
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,a_r,NumDGAIter,NumDGAIter)
!     ----------------------------------------------------------------
!
      isz = iter+1
      pb => b_r(1:isz,1)
      do i = 1,iter
         a_r(i,isz) = one
         a_r(isz,i) = one
         pb(i) = zero
      enddo
      a_r(isz,isz) = zero
      pb(isz) = one
!
!     ================================================================
!     call the solver.................................................
!     ================================================================
!     ----------------------------------------------------------------
      call dgaleq_r(a_r,b_r,isz,NumDGAIter)
!     ----------------------------------------------------------------
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = DGAMix(id)%vlen
         alpha = DGAMix(id)%alpha
!        =============================================================
!        perform the mixing...........................................
!        =============================================================
         pf_old    => DGAMix(id)%f_old_r(1:NumDGAIter,1:vlen)
         pf_new    => DGAMix(id)%f_new_r(1:NumDGAIter,1:vlen)
         pvect_new => DGAMix(id)%vector_new_r(1:vlen)
         do i = 1,vlen
            pfi_t => pf_old(1:iter,i)
            pfo_t => pf_new(1:iter,i)
!           ----------------------------------------------------------
            pvect_new(i) = dgasad_r( NumDGAIter, pfi_t, pfo_t, &
                                   pb, iter, alpha )
!           ----------------------------------------------------------
         enddo
      enddo
   enddo
!
   end subroutine calDGAMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDGAMixing_c( list_q )
!  ===================================================================
!
   implicit none
!
   type(MixListCmplxStruct), target :: list_q
!
   integer(kind=IntKind) :: i, ii, jj, n, iter, isz, vlen
   integer(kind=IntKind) :: idt, idq, id
!
   real(kind=RealKind) :: alpha
   complex(kind=CmplxKind), pointer :: pvect_old(:), pvect_new(:) , pb(:)
   complex(kind=CmplxKind), pointer :: pf_old(:,:), pf_new(:,:)
   complex(kind=CmplxKind), pointer :: pfi_t(:), pfo_t(:)
!
   type(MixListCmplxStruct), pointer :: plq_tmp
!
   iter = min( iter_count, NumDGAIter )
   n    = mod( iter_count-1, NumDGAIter ) + 1
!
   if ( .not.InitializedWkSpace ) then
      call setWorkingSpace_c(MixingMethod)
   endif
!
   plq_tmp => list_q
   do idt = 1,NumTypes
      a_c = czero
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = plq_tmp%size
         if ( .not.DGAMix(id)%Initialized ) then
            call setDGASpace_c(id,vlen)
         endif
         DGAMix(id)%vector_old_c => plq_tmp%vector_old(1:vlen)
         DGAMix(id)%vector_new_c => plq_tmp%vector_new(1:vlen)
         plq_tmp => plq_tmp%next
      enddo
!
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         pvect_old => DGAMix(id)%vector_old_c(1:vlen)
         pf_old    => DGAMix(id)%f_old_c(1:NumDGAIter,1:vlen)
         pvect_new => DGAMix(id)%vector_new_c(1:vlen)
         pf_new    => DGAMix(id)%f_new_c(1:NumDGAIter,1:vlen)
         vlen = DGAMix(id)%vlen
         do i = 1,vlen
            pf_old(n,i) = pvect_old(i)
            pf_new(n,i) = pvect_new(i)
         enddo
!
         do jj = 1,iter
            do ii = 1,jj
!              -------------------------------------------------------
               a_c(ii,jj) = a_c(ii,jj) +                                  &
                          dgarms_c(pf_old,pf_new,ii,jj,vlen)
!              -------------------------------------------------------
               if ( ii.ne.jj ) then
                  a_c(jj,ii) = a_c(jj,ii) + a_c(ii,jj)
               endif
            enddo
         enddo
!
      enddo
!     ================================================================
!     sum the a(i,j) matrix over all other atoms......................
!     ================================================================
!     ----------------------------------------------------------------
      call GlobalSumInGroup(GroupID,a_c,NumDGAIter,NumDGAIter)
!     ----------------------------------------------------------------
!
      isz = iter+1
      pb => b_c(1:isz,1)
      do i = 1,iter
         a_c(i,isz) = cone
         a_c(isz,i) = cone
         pb(i) = czero
      enddo
      a_c(isz,isz) = zero
      pb(isz) = one
!
!     ================================================================
!     call the solver.................................................
!     ================================================================
!     ----------------------------------------------------------------
      call dgaleq_c(a_c,b_c,isz,NumDGAIter)
!     ----------------------------------------------------------------
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = DGAMix(id)%vlen
         alpha = DGAMix(id)%alpha
!        =============================================================
!        perform the mixing...........................................
!        =============================================================
         pf_old    => DGAMix(id)%f_old_c(1:NumDGAIter,1:vlen)
         pf_new    => DGAMix(id)%f_new_c(1:NumDGAIter,1:vlen)
         pvect_new => DGAMix(id)%vector_new_c(1:vlen)
         do i = 1,vlen
            pfi_t => pf_old(1:iter,i)
            pfo_t => pf_new(1:iter,i)
!           ----------------------------------------------------------
            pvect_new(i) = dgasad_c( NumDGAIter, pfi_t, pfo_t, &
                                   pb, iter, alpha )
!           ----------------------------------------------------------
         enddo
      enddo
   enddo
!
   end subroutine calDGAMixing_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setBroydenMixing( idt, idq, amix )
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: idt, idq
!  integer(kind=IntKind), optional   :: nbroyiter
!  integer(kind=IntKind), optional   :: inv_method
!
   real(kind=RealKind), intent(in) :: amix
!
   integer(kind=IntKind) :: id
!
   if ( .not.Initialized ) then
      call ErrorHandler('setBroydenMix',                              &
          'The MixingModule needs to be initialized first')
   endif
!
   if ( idt<0 .or. idt>NumTypes ) then
      call ErrorHandler('setBroydenMixing',"Invalid type index",idq)
   endif
   if ( idq<0 .or. idq>NumQuantities(idt) ) then
      call ErrorHandler('setBroydenMixing',"Invalid quantity index",idq)
   endif
!
   if ( .not.allocated(BroydenMix) ) then
      allocate( BroydenMix(NumTotalMix) )
      do id = 1, NumTotalMix
         BroydenMix(id)%Initialized = .false.
      enddo
   endif  
!
   NumBroydenIter = MaxNumBroydenIter
!  if ( present(nbroyiter) ) then
!     NumBroydenIter = nbroyiter
!     if ( nbroyiter > MaxNumBroydenIter ) then
!        call WarningHandler('setBroydenMixing', &
!            'The number of Broyden iterations exeedes the maximum limit')
!        NumBroydenIter = MaxNumBroydenIter
!     endif
!  endif
!
   id = getMixID(idt,idq)
   if ( .not.BroydenMix(id)%Initialized ) then
      BroydenMix(id)%alpha = amix
   endif
!
!  if (present(inv_method)) BroydenInvMethod = inv_method
!
   if (BroydenInvMethod<0 .or. BroydenInvMethod>1) then
      call WarningHandler('initMixing', &
           'Invalid inversion method switching to default(invert1)')
      BroydenInvMethod = 0
   endif 
!
   if ( MixingMethod /= 1 ) then
      MixingMethod = 1
   endif
!
   end subroutine setBroydenMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delBroydenMixing( )
!  ===================================================================
   implicit none
!
   integer(kind=IntKind) :: id
!
   if ( allocated(BroydenMix) ) then
      do id = 1,NumTotalMix
         if ( .not.BroydenMix(id)%Initialized ) cycle
         BroydenMix(id)%Initialized = .false.
         if ( associated(BroydenMix(id)%u_r) ) deallocate( BroydenMix(id)%u_r )
         if ( associated(BroydenMix(id)%vt_r) ) deallocate( BroydenMix(id)%vt_r )
         if ( associated(BroydenMix(id)%f_r) ) deallocate( BroydenMix(id)%f_r )
         if ( associated(BroydenMix(id)%df_r) ) deallocate( BroydenMix(id)%df_r )
         if ( associated(BroydenMix(id)%w_r) ) deallocate( BroydenMix(id)%w_r )
         if ( associated(BroydenMix(id)%vold_r) ) deallocate( BroydenMix(id)%vold_r )
         if ( associated(BroydenMix(id)%u_c) ) deallocate( BroydenMix(id)%u_c )
         if ( associated(BroydenMix(id)%vt_c) ) deallocate( BroydenMix(id)%vt_c )
         if ( associated(BroydenMix(id)%f_c) ) deallocate( BroydenMix(id)%f_c )
         if ( associated(BroydenMix(id)%df_c) ) deallocate( BroydenMix(id)%df_c )
         if ( associated(BroydenMix(id)%w_c) ) deallocate( BroydenMix(id)%w_c )
         if ( associated(BroydenMix(id)%vold_c) ) deallocate( BroydenMix(id)%vold_c )
      enddo
!      deallocate( BroydenMix )
   endif
!
   end subroutine delBroydenMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calBroydenMixing_r( list_q )
!  ===================================================================
!
   implicit none
!
   type(MixListRealStruct), target :: list_q
!
   integer(kind=IntKind) :: i, j, k, info, iter, idq, idt, id
   integer(kind=IntKind) :: lastm1, nn, vlen
!
   real(kind=RealKind) :: fac1, fac2, fnorm, dfnorm, w0, f0, df0
   real(kind=RealKind) :: aij, gmi, cmj, wtmp, bij, cij
   real(kind=RealKind) :: msgbuf(2)
   real(kind=RealKind) :: rms
   real(kind=RealKind) :: alpha
!
   real(kind=RealKind), pointer :: pf(:), pdf(:), pu(:,:), pw(:), pvt(:,:)
   real(kind=RealKind), pointer :: pvect_old(:), pvect_new(:), pvold(:)
   real(kind=RealKind), pointer :: pad(:), pbd(:)
!
   type(MixListRealStruct), pointer :: plq_tmp
!
   iter = iter_count
   plq_tmp => list_q
!
   if ( .not.InitializedWkSpace ) then
      call setWorkingSpace_r(MixingMethod)
   endif
!
   do idt = 1,NumTypes
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = plq_tmp%size
         BroydenMix(id)%rms = plq_tmp%rms
         BroydenMix(id)%weight = plq_tmp%weight
         BroydenMix(id)%vector_old_r => plq_tmp%vector_old(1:vlen)
         BroydenMix(id)%vector_new_r => plq_tmp%vector_new(1:vlen)
         if ( .not.BroydenMix(id)%Initialized ) then
            call setBroydenSpace_r(id,vlen)
         endif
         plq_tmp => plq_tmp%next
      enddo
   enddo
!
   if ( iter_count == 1) then
!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!
      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen =BroydenMix(id)%vlen
            pf => BroydenMix(id)%f_r(1:vlen)
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
!
            pvold    => BroydenMix(id)%vold_r(1:vlen)
            do k = 1,vlen
               pf(k) = pvect_new(k) - pvect_old(k)
            enddo
            do k = 1,vlen
               pvold(k) = pvect_old(k)
            enddo
!
            do k = 1,vlen          ! this is the linear mixing
               pvect_new(k) = pvect_old(k) + alpha * pf(k)
            enddo
         enddo
      enddo
!
   else
!     ===============================================================
!     iter_count > 1: this is where the non-linear mixing is done
!     ===============================================================
      lastIter = lastIter + 1      ! update pointers
      lastm1   = lastIter - 1
!
!     ===============================================================
!     set current lenght of broyden cycle
!     ===============================================================
      if ( iter > NumBroydenIter ) then 
         nn = NumBroydenIter
      else
         nn = lastIter
      endif
!
      LoopType: do idt = 1,NumTypes
!        ============================================================
!        set weighting factor for the zeroth iteration
!        ============================================================
         w0= ten2m2
!
!        ============================================================
!        find: f[i] := vector(2)[i] - vector(1)[i]
!              df   := f[i] - f[i-1]
!        ============================================================
         dfnorm = zero
         fnorm = zero
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
            pf  => BroydenMix(id)%f_r(1:vlen)
            pdf => BroydenMix(id)%df_r(1:vlen)
            if ( vlen==1 ) then
               pvect_new(1) = pvect_old(1) + alpha * pf(1)
               cycle
            endif
            do k = 1,vlen
               pdf(k) = pvect_new(k) - pvect_old(k) - pf(k)
               pf(k)  = pvect_new(k) - pvect_old(k)
            enddo
!           =========================================================
!           find: fnorm  := |f|
!                 dfnorm := |df|
!           =========================================================
            df0 = ZERO
            f0 = ZERO
            do k = 1,vlen
               df0 = df0 + pdf(k)*pdf(k)
               f0  = f0  + pf(k)*pf(k)
            enddo
            dfnorm = dfnorm + BroydenMix(id)%weight*df0
            fnorm  = fnorm  + BroydenMix(id)%weight*f0
         enddo
!
         msgbuf(1) = dfnorm
         msgbuf(2) = fnorm
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID, msgbuf,2)
!        ------------------------------------------------------------
         dfnorm = sqrt( msgbuf(1) )
         fnorm  = sqrt( msgbuf(2) )
!        ============================================================
!        set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|
!             vold      := vector(1) 
!             vector(1) := df/|df|
!        ============================================================
!
         fac2 = one/dfnorm
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen  = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            fac1 = alpha*fac2
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
            pu    => BroydenMix(id)%u_r(1:vlen,1:NumBroydenIter)
            pvt   => BroydenMix(id)%vt_r(1:vlen,1:NumBroydenIter)
            pdf   => BroydenMix(id)%df_r(1:vlen)
            pvold => BroydenMix(id)%vold_r(1:vlen)
            do k = 1,vlen
               pvect_new(k) = fac1*pdf(k) + fac2*(pvect_old(k) - pvold(k))
               pvold(k)   = pvect_old(k)
               pvect_old(k) = fac2*pdf(k)
            enddo
!           =========================================================
!           store vector(1) and vector(2) in the stacks u and vt
!           v(i) = v(i+1)
!           u(i) = u(i+1)
!           u(nn) = vector(1)
!           v(nn) = vector(2)
!           =========================================================
!           ---------------------------------------------------------
            call broy_sav_r( pu, pvt, pvect_old, pvect_new, iter-1, & 
                             NumBroydenIter, vlen )
!           ---------------------------------------------------------
         enddo
!        ============================================================
!        calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:
!           a(i,j<i) = v(j)*v(l)
!           a(i,i) = v(i)*v(i)
!           a(j,i) = a(i,j)
!           cm(i) = sum_(l=1,i) v(l)*f
!        ============================================================
         do j = 1,nn - 1          ! off diagonal elements of a(i,j)
            do i = j+1,nn
               aij = zero
               do idq = 1,NumQuantities(idt)
                  id = getMixID(idt,idq)
                  vlen = BroydenMix(id)%vlen
                  if ( vlen==1 ) then
                     cycle
                  endif
                  pvt => BroydenMix(id)%vt_r(1:vlen,1:NumBroydenIter)
                  bij = ZERO
                  do k = 1,vlen
                     bij = bij + pvt(k,j)*pvt(k,i)
                  enddo
                  aij = aij + BroydenMix(id)%weight*bij
               enddo
!              ------------------------------------------------------
!              call GlobalSumInGroup(GroupID,aij)
!              ------------------------------------------------------
               a_r(i,j) = aij
               a_r(j,i) = aij
            enddo
         enddo
!
         do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
            aij = zero
            cmj = zero
            do idq = 1, NumQuantities(idt)
               id = getMixID(idt,idq)
               vlen = BroydenMix(id)%vlen
               if ( vlen==1 ) then
                  cycle
               endif
               pvt => BroydenMix(id)%vt_r(1:vlen,1:NumBroydenIter)
               pf  => BroydenMix(id)%f_r(1:vlen)
               cij = ZERO
               bij = ZERO
               do k = 1,vlen
                  cij = cij + pvt(k,i)*pf(k)
                  bij = bij + pvt(k,i)*pvt(k,i)
               enddo
               cmj = cmj + BroydenMix(id)%weight*cij
               aij = aij + BroydenMix(id)%weight*bij
            enddo
!           msgbuf(1) = aij
!           msgbuf(2) = cmj
!           ---------------------------------------------------------
!           call GlobalSumInGroup(GroupID,msgbuf,2)
!           ---------------------------------------------------------
!           aij = msgbuf(1)
!           cmj = msgbuf(2)
            a_r(i,i) = aij
            cm_r(i) = cmj
         enddo
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,a_r,NumBroydenIter,NumBroydenIter)
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,cm_r,NumBroydenIter)
!        ------------------------------------------------------------
!        ============================================================
!        shift down weights in stack
!
!        (TCS, bug fixed 8/5/97: replace iter by iter-1 -> see broy_sav)
!        ============================================================
         do idq = 1,NumQuantities(idt)
!
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
!
            pvold => BroydenMix(id)%vold_r(1:vlen)
            pw => BroydenMix(id)%w_r(1:NumBroydenIter)
            pu => BroydenMix(id)%u_r(1:vlen,1:NumBroydenIter)
            pf => BroydenMix(id)%f_r(1:vlen)
!
            if ( iter-1 > NumBroydenIter ) then
               do i = 1,NumBroydenIter-1
                  pw(i) = pw(i+1)
               enddo
            endif
!
            rms = BroydenMix(id)%rms
            wtmp = zero
            if ( rms > ten2m9 ) wtmp = TWO*sqrt(ten2m2/rms)
            if ( wtmp < ONE ) wtmp = ONE
            if ( iter > NumBroydenIter ) then
               pw(NumBroydenIter) = wtmp
            else
               pw(lastIter) = wtmp       !w(lastm1)=wtmp
            endif
!           =========================================================
!           now calculate the b-matrix:
!           b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
!           =========================================================
            if ( BroydenInvMethod==0 .and. vlen>3*NumBroydenIter ) then
!              ======================================================
!              use invert1
!              ======================================================
               do i = 1,nn
                  do j = 1,nn
                     d_r(j,i) = a_r(j,i)*pw(j)*pw(i)
                     b_r(j,i) = zero
                  enddo
                  b_r(i,i) = one
                  d_r(i,i) = w0**2 + a_r(i,i)*pw(i)*pw(i)
               enddo
!
!              this is very unlikely
!
               if ( 3*NumBroydenIter > vlen ) then 
                  call ErrorHandler( "callBroyden", & 
                      "need larger dimension for MaxSizeBroydenVect")
               endif
               pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
               pad =>                                                 &
               BroydenMix(id)%vector_old_r(NumBroydenIter+1:2*NumBroydenIter)
               pbd =>                                                 &
               BroydenMix(id)%vector_old_r(2*NumBroydenIter+1:3*NumBroydenIter)
!              ------------------------------------------------------
               call invert1_r(d_r,b_r,nn,pvect_old,pad,pbd,NumBroydenIter)
!              ------------------------------------------------------
            else
!              ======================================================
!              use lapack
!              ======================================================
               do i = 1,nn
                  do j = 1,nn
                     b_r(j,i) = a_r(j,i)*pw(j)*pw(i)
                  enddo
                  b_r(i,i) = w0**2 + a_r(i,i)*pw(i)*pw(i)
               enddo
               if ( .not.allocated(ipiv) ) then
                  allocate( ipiv(NumBroydenIter) )
               endif
!              ------------------------------------------------------
               call dgetrf( nn, nn, b_r, NumBroydenIter,ipiv, info )
!              ------------------------------------------------------
               call dgetri( nn, b_r, NumBroydenIter, ipiv, d_r, nn, info )
!              ------------------------------------------------------
!              write(6,*) ' optimum lwork', d(1,1)
            endif
!
            do k = 1,vlen
               pvect_new(k) = pvold(k) + alpha*pf(k)
            enddo
            do i = 1,nn
               gmi = zero
               do j = 1,nn
                  gmi = gmi + cm_r(j)*b_r(j,i)*pw(j)
               enddo
               do k = 1,vlen
                  pvect_new(k) = pvect_new(k) - gmi*pu(k,i)*pw(i)
               enddo
            enddo
         enddo
      enddo LoopType
   endif
!
   end subroutine calBroydenMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calBroydenMixing_c( list_q )
!  ===================================================================
!
   implicit none
!
   type(MixListCmplxStruct), target :: list_q
!
   integer(kind=IntKind) :: i, j, k, info, iter, idq, idt, id
   integer(kind=IntKind) :: lastm1, nn, vlen
!
   real(kind=RealKind) :: fac1, fac2, fnorm, dfnorm, w0, df0, f0
   real(kind=RealKind) :: wtmp
   real(kind=RealKind) :: msgbuf(2)
   real(kind=RealKind) :: rms
   real(kind=RealKind) :: alpha
!
   complex(kind=CmplxKind) :: aij, gmi, cmj, bij, cij
   complex(kind=CmplxKind), pointer :: pf(:), pdf(:), pu(:,:), pw(:), pvt(:,:)
   complex(kind=CmplxKind), pointer :: pvect_old(:), pvect_new(:), pvold(:)
   complex(kind=CmplxKind), pointer :: pad(:), pbd(:)
!
   type(MixListCmplxStruct), pointer :: plq_tmp
!
   iter = iter_count
   plq_tmp => list_q
!
   if ( .not.InitializedWkSpace ) then
      call setWorkingSpace_c(MixingMethod)
   endif
!
   do idt = 1,NumTypes
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = plq_tmp%size
         BroydenMix(id)%rms = plq_tmp%rms
         BroydenMix(id)%weight = plq_tmp%weight
         BroydenMix(id)%vector_old_c => plq_tmp%vector_old(1:vlen)
         BroydenMix(id)%vector_new_c => plq_tmp%vector_new(1:vlen)
         if ( .not.BroydenMix(id)%Initialized ) then
            call setBroydenSpace_c(id,vlen)
         endif
         plq_tmp => plq_tmp%next
      enddo
   enddo
!
   if ( iter_count == 1) then
!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!
      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen = BroydenMix(id)%vlen
            pf => BroydenMix(id)%f_c(1:vlen)
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
!
            pvold    => BroydenMix(id)%vold_c(1:vlen)
            do k = 1,vlen
               pf(k) = pvect_new(k) - pvect_old(k)
            enddo
            do k = 1,vlen
               pvold(k) = pvect_old(k)
            enddo
!
            do k = 1,vlen          ! this is the linear mixing
               pvect_new(k) = pvect_old(k) + alpha * pf(k)
            enddo
         enddo
      enddo
!
   else
!     ===============================================================
!     iter_count > 1: this is where the non-linear mixing is done
!     ===============================================================
      lastIter = lastIter + 1      ! update pointers
      lastm1   = lastIter - 1
!
!     ===============================================================
!     set current lenght of broyden cycle
!     ===============================================================
      if ( iter > NumBroydenIter ) then
         nn = NumBroydenIter
      else
         nn = lastIter
      endif
!
      LoopType: do idt = 1,NumTypes
!        ============================================================
!        set weighting factor for the zeroth iteration
!        ============================================================
         w0= ten2m2
!
!        ============================================================
!        find: f[i] := vector(2)[i] - vector(1)[i]
!              df   := f[i] - f[i-1]
!        ============================================================
         dfnorm = zero
         fnorm = zero
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
            pf  => BroydenMix(id)%f_c(1:vlen)
            pdf => BroydenMix(id)%df_c(1:vlen)
            if ( vlen==1 ) then
               pvect_new(1) = pvect_old(1) + alpha * pf(1)
               cycle
            endif
            do k = 1,vlen
               pdf(k) = pvect_new(k) - pvect_old(k) - pf(k)
               pf(k)  = pvect_new(k) - pvect_old(k)
            enddo
!           =========================================================
!           find: fnorm  := |f|
!                 dfnorm := |df|
!           =========================================================
            df0 = ZERO
            f0 = ZERO
            do k = 1,vlen
               df0 = df0 + real(pdf(k)*conjg(pdf(k)),kind=RealKind)
               f0  = f0  + real(pf(k)*conjg(pf(k)),kind=RealKind)
            enddo
            dfnorm = dfnorm + BroydenMix(id)%weight*df0
            fnorm  = fnorm  + BroydenMix(id)%weight*f0
         enddo
!
         msgbuf(1) = dfnorm
         msgbuf(2) = fnorm
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID, msgbuf,2)
!        ------------------------------------------------------------
         dfnorm = sqrt( msgbuf(1) )
         fnorm  = sqrt( msgbuf(2) )
!        ============================================================
!        set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|
!             vold      := vector(1) 
!             vector(1) := df/|df|
!        ============================================================
!
         fac2 = one/dfnorm
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen  = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            fac1 = alpha*fac2
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
            pu    => BroydenMix(id)%u_c(1:vlen,1:NumBroydenIter)
            pvt   => BroydenMix(id)%vt_c(1:vlen,1:NumBroydenIter)
            pdf   => BroydenMix(id)%df_c(1:vlen)
            pvold => BroydenMix(id)%vold_c(1:vlen)
            do k = 1,vlen
               pvect_new(k) = fac1*pdf(k) + fac2*(pvect_old(k) - pvold(k))
               pvold(k)   = pvect_old(k)
               pvect_old(k) = fac2*pdf(k)
            enddo
!           =========================================================
!           store vector(1) and vector(2) in the stacks u and vt
!           v(i) = v(i+1)
!           u(i) = u(i+1)
!           u(nn) = vector(1)
!           v(nn) = vector(2)
!           =========================================================
!           ---------------------------------------------------------
            call broy_sav_c( pu, pvt, pvect_old, pvect_new, iter-1, & 
                             NumBroydenIter, vlen )
!           ---------------------------------------------------------
         enddo
!        ============================================================
!        calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:
!           a(i,j<i) = v(j)*v(l)
!           a(i,i) = v(i)*v(i)
!           a(j,i) = a(i,j)
!           cm(i) = sum_(l=1,i) v(l)*f
!        ============================================================
         do j = 1,nn - 1          ! off diagonal elements of a(i,j)
            do i = j+1,nn
               aij = czero
               do idq = 1,NumQuantities(idt)
                  id = getMixID(idt,idq)
                  vlen = BroydenMix(id)%vlen
                  if ( vlen==1 ) then
                     cycle
                  endif
                  pvt => BroydenMix(id)%vt_c(1:vlen,1:NumBroydenIter)
                  bij = CZERO
                  do k = 1,vlen
                     bij = bij + conjg(pvt(k,j))*pvt(k,i)
                  enddo
                  aij = aij + BroydenMix(id)%weight*bij
               enddo
!              ------------------------------------------------------
!              call GlobalSumInGroup(GroupID,aij)
!              ------------------------------------------------------
               a_c(i,j) = aij
               a_c(j,i) = conjg(aij)
            enddo
         enddo
!
         do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
            aij = czero
            cmj = czero
            do idq = 1, NumQuantities(idt)
               id = getMixID(idt,idq)
               vlen = BroydenMix(id)%vlen
               if ( vlen==1 ) then
                  cycle
               endif
               pvt => BroydenMix(id)%vt_c(1:vlen,1:NumBroydenIter)
               pf  => BroydenMix(id)%f_c(1:vlen)
               cij = CZERO
               bij = CZERO
               do k = 1,vlen
                  cij = cij + conjg(pvt(k,i))*pf(k)
                  bij = bij + conjg(pvt(k,i))*pvt(k,i)
               enddo
               cmj = cmj + BroydenMix(id)%weight*cij
               aij = aij + BroydenMix(id)%weight*bij
            enddo
!           msgbuf(1) = aij
!           msgbuf(2) = cmj
!           ---------------------------------------------------------
!           call GlobalSumInGroup(GroupID,msgbuf,2)
!           ---------------------------------------------------------
!           aij = msgbuf(1)
!           cmj = msgbuf(2)
            a_c(i,i) = aij
            cm_c(i) = cmj
         enddo
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,a_c,NumBroydenIter,NumBroydenIter)
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,cm_c,NumBroydenIter)
!        ------------------------------------------------------------
!        ============================================================
!        shift down weights in stack
!
!        (TCS, bug fixed 8/5/97: replace iter by iter-1 -> see broy_sav)
!        ============================================================
         do idq = 1,NumQuantities(idt)
!
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
!
            pvold => BroydenMix(id)%vold_c(1:vlen)
            pw => BroydenMix(id)%w_c(1:NumBroydenIter)
            pu => BroydenMix(id)%u_c(1:vlen,1:NumBroydenIter)
            pf => BroydenMix(id)%f_c(1:vlen)
!
            if ( iter-1 > NumBroydenIter ) then
               do i = 1,NumBroydenIter-1
                  pw(i) = pw(i+1)
               enddo
            endif
!
            rms = BroydenMix(id)%rms
            wtmp = zero
            if ( rms > ten2m9 ) wtmp = TWO*sqrt(ten2m2/rms)
            if ( wtmp < ONE ) wtmp = ONE
            if ( iter > NumBroydenIter ) then
               pw(NumBroydenIter) = wtmp
            else
               pw(lastIter) = wtmp       !w(lastm1)=wtmp
            endif
!           =========================================================
!           now calculate the b-matrix:
!           b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
!           =========================================================
            if ( BroydenInvMethod==0 .and. vlen>3*NumBroydenIter ) then
!              ======================================================
!              use invert1
!              ======================================================
               do i = 1,nn
                  do j = 1,nn
                     d_c(j,i) = a_c(j,i)*pw(j)*pw(i)
                     b_c(j,i) = czero
                  enddo
                  b_c(i,i) = cone
                  d_c(i,i) = w0**2 + a_c(i,i)*pw(i)*pw(i)
               enddo
!
!              this is very unlikely
!
               if ( 3*NumBroydenIter > vlen ) then 
                  call ErrorHandler( "callBroyden", & 
                      "need larger dimension for MaxSizeBroydenVect")
               endif
               pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
               pad =>                                                 &
               BroydenMix(id)%vector_old_c(NumBroydenIter+1:2*NumBroydenIter)
               pbd =>                                                 &
               BroydenMix(id)%vector_old_c(2*NumBroydenIter+1:3*NumBroydenIter)
!              ------------------------------------------------------
               call invert1_c(d_c,b_c,nn,pvect_old,pad,pbd,NumBroydenIter)
!              ------------------------------------------------------
            else
!              ======================================================
!              use lapack
!              ======================================================
               b_c = CZERO
               do i = 1,nn
                  do j = 1,nn
                     b_c(j,i) = a_c(j,i)*pw(j)*pw(i)
                  enddo
                  b_c(i,i) = w0**2 + a_c(i,i)*pw(i)*pw(i)
               enddo
               if ( .not.allocated(ipiv) ) then
                  allocate( ipiv(NumBroydenIter) )
                  ipiv = 0
               endif
!              ------------------------------------------------------
               call zgetrf( nn, nn, b_c, NumBroydenIter,ipiv, info )
!              ------------------------------------------------------
               call zgetri( nn, b_c, NumBroydenIter, ipiv, d_c, nn, info )
!              ------------------------------------------------------
!              write(6,*) ' optimum lwork', d(1,1)
            endif
!
            do k = 1,vlen
               pvect_new(k) = pvold(k) + alpha*pf(k)
            enddo
            do i = 1,nn
               gmi = czero
               do j = nn, 1, -1
                  gmi = gmi + cm_c(j)*b_c(j,i)*pw(j)
               enddo
               do k = 1,vlen
                  pvect_new(k) = pvect_new(k) - gmi*pu(k,i)*pw(i)
               enddo
            enddo
         enddo
      enddo LoopType
   endif
!
   end subroutine calBroydenMixing_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine broy_sav_r( fins, fots, vector_old, vector_new, itscf,  & 
                          istore, ivsiz )
!  ==================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: itscf
   integer(kind=IntKind), intent(in) :: ivsiz
   integer(kind=IntKind), intent(in) :: istore
!
   real(kind=RealKind), target :: vector_old(:)
   real(kind=RealKind), target :: vector_new(:)
   real(kind=RealKind), target :: fins(:,:)
   real(kind=RealKind), target :: fots(:,:)
!
   integer(kind=IntKind) :: i, j
!
!  ==================================================================
!  write(6,'('' IN BROY_SAV: istore,itscf '',2i5)') istore,itscf
!  ==================================================================
   if ( itscf <= istore ) then
!     ===============================================================
!     Load the first istore iterations in increasing iteration count
!     ===============================================================
      do i = 1,ivsiz
         fins(i,itscf) = vector_new(i)
      enddo
!
      do i = 1,ivsiz
         fots(i,itscf) = vector_old(i)
      enddo
!
   else
!     ===============================================================
!     Re-load so that the ordering is in increasing iteration count
!     ===============================================================
      do j = 1,istore - 1
!        write(6,'('' IN BROY_SAV: j,j+1 '',2i5)') j,j+1
         do i = 1,ivsiz
            fins(i,j) = fins(i,j+1)
         enddo
!
         do i = 1,ivsiz
            fots(i,j) = fots(i,j+1)
         enddo
!
      enddo
!     ===============================================================
!     Load current charge densities in the last storage location
!     ===============================================================
      do i = 1,ivsiz
         fins(i,istore) = vector_new(i)
      enddo
!
      do i = 1,ivsiz
         fots(i,istore) = vector_old(i)
      enddo
!
   endif
!
   end subroutine broy_sav_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine broy_sav_c( fins, fots, vector_old, vector_new, itscf,  & 
                          istore, ivsiz )
!  ==================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: itscf
   integer(kind=IntKind), intent(in) :: ivsiz
   integer(kind=IntKind), intent(in) :: istore
!
   complex(kind=CmplxKind), target :: vector_old(:)
   complex(kind=CmplxKind), target :: vector_new(:)
   complex(kind=CmplxKind), target :: fins(:,:)
   complex(kind=CmplxKind), target :: fots(:,:)
!
   integer(kind=IntKind) :: i, j
!
!  ==================================================================
!  write(6,'('' IN BROY_SAV: istore,itscf '',2i5)') istore,itscf
!  ==================================================================
   if ( itscf <= istore ) then
!     ===============================================================
!     Load the first istore iterations in increasing iteration count
!     ===============================================================
      do i = 1,ivsiz
         fins(i,itscf) = vector_new(i)
      enddo
!
      do i = 1,ivsiz
         fots(i,itscf) = vector_old(i)
      enddo
!
   else
!     ===============================================================
!     Re-load so that the ordering is in increasing iteration count
!     ===============================================================
      do j = 1,istore - 1
!        write(6,'('' IN BROY_SAV: j,j+1 '',2i5)') j,j+1
         do i = 1,ivsiz
            fins(i,j) = fins(i,j+1)
         enddo
!
         do i = 1,ivsiz
            fots(i,j) = fots(i,j+1)
         enddo
!
      enddo
!     ===============================================================
!     Load current charge densities in the last storage location
!     ===============================================================
      do i = 1,ivsiz
         fins(i,istore) = vector_new(i)
      enddo
!
      do i = 1,ivsiz
         fots(i,istore) = vector_old(i)
      enddo
!
   endif
!
   end subroutine broy_sav_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine invert1_r( pa, pb, m, ptd, pad, pbd, mm ) 
!  ==================================================================
!  ==================================================================
!  temporary for testting, should be lapack routine in future!
!  ==================================================================
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: m, mm
!
   real(kind=RealKind), target :: pa(:,:)
   real(kind=RealKind), target :: pb(:,:)
   real(kind=RealKind), target :: ptd(:)
   real(kind=RealKind), target :: pad(:)
   real(kind=RealKind), target :: pbd(:)
!
   integer(kind=IntKind) :: n, i, j, k
!
   real(kind=RealKind) :: atmp
!
!  ==================================================================
!  parameter (mm=5)
!  ==================================================================
!
!  ==================================================================
!  subroutine to preform gaussian elimination
!             no zeros along the diagonal
!  ==================================================================
!
   n = m
   if ( n > mm ) then
      call ErrorHandler('invert1',' invert: matrix a too large')
   endif
!
   do i = 1,n
      atmp = pa(i,i)
      if ( abs(atmp) < ten2m8 ) then
         call ErrorHandler("invert1","matrix has zero diagonal(row)",i)
      endif
   enddo
!
   if ( n /= 1 ) then
!
      do i = 1,n
!
         do j = 1,n
            ptd(j) = pa(j,i)/pa(i,i)
         enddo
!
         ptd(i) = zero
!
         do k = 1,n
            pbd(k) = pb(i,k)
            pad(k) = pa(i,k)
         enddo
!
         do k = 1,n
            do j = 1,n
               pb(j,k) = pb(j,k)-(ptd(j)*pbd(k))
               pa(j,k) = pa(j,k)-(ptd(j)*pad(k))
            enddo
         enddo
!
      enddo
!
      do i = 1,n
         do j = 1,n
            pb(j,i) = pb(j,i)/pa(j,j)
         enddo
      enddo
!
   else
      pb(1,1) = one/pa(1,1)
   endif
!
   end subroutine invert1_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine invert1_c( pa, pb, m, ptd, pad, pbd, mm ) 
!  ==================================================================
!  ==================================================================
!  temporary for testting, should be lapack routine in future!
!  ==================================================================
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: m, mm
!
   complex(kind=CmplxKind), target :: pa(:,:)
   complex(kind=CmplxKind), target :: pb(:,:)
   complex(kind=CmplxKind), target :: ptd(:)
   complex(kind=CmplxKind), target :: pad(:)
   complex(kind=CmplxKind), target :: pbd(:)
!
   integer(kind=IntKind) :: n, i, j, k
!
   complex(kind=CmplxKind) :: atmp
!
!  ==================================================================
!  parameter (mm=5)
!  ==================================================================
!
!  ==================================================================
!  subroutine to preform gaussian elimination
!             no zeros along the diagonal
!  ==================================================================
!
   n = m
   if ( n > mm ) then
      call ErrorHandler('invert1',' invert: matrix a too large')
   endif
!
   do i = 1,n
      atmp = pa(i,i)
      if ( abs(atmp) < ten2m8 ) then
         call ErrorHandler("invert1","matrix has zero diagonal(row)",i)
      endif
   enddo
!
   if ( n /= 1 ) then
!
      do i = 1,n
!
         do j = 1,n
            ptd(j) = pa(j,i)/pa(i,i)
         enddo
!
         ptd(i) = czero
!
         do k = 1,n
            pbd(k) = conjg(pb(i,k))
            pad(k) = conjg(pa(i,k))
         enddo
!
         do k = 1,n
            do j = 1,n
               pb(j,k) = pb(j,k) - ptd(j)*pbd(k)
               pa(j,k) = pa(j,k) - ptd(j)*pad(k)
            enddo
         enddo
!
      enddo
!
      do i = 1,n
         do j = 1,n
            pb(j,i) = pb(j,i)/pa(j,j)
         enddo
      enddo
!
   else
      pb(1,1) = one/pa(1,1)
   endif
!
   end subroutine invert1_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine dgaleq_r( a, y, n, ipits)
!  ==================================================================
!
   implicit   none
!
   integer(kind=IntKind), intent(in) :: n, ipits
!
   real(kind=RealKind), intent(inout) :: a(ipits+1,ipits+1)
   real(kind=RealKind), intent(inout) :: y(ipits+1)
!
   integer(kind=IntKind) :: i, j, k, ij
!
   real(kind=RealKind) :: f1,f2
!
!  ==================================================================
!  ******************************************************************
!  Important note: the following codes may will result a problem:
!         a(i-1,i-1) becomes zero!!!
!  It needs to be fixed.
!  ******************************************************************
!  ==================================================================
   do i = 2,n
      f1 = -one/a(i-1,i-1)
      do j = i,n
         f2 = f1*a(j,i-1)
         do k = 1,n
            a(j,k) = a(j,k)+f2*a(i-1,k)
         enddo
         y(j) = y(j)+f2*y(i-1)
      enddo
   enddo
   y(n) = y(n)/a(n,n)
   do ij = 1,n-1
      i = n-ij
      do j = 1,ij
         y(i) = y(i)-y(i+j)*a(i,i+j)
      enddo
      y(i) = y(i)/a(i,i)
   enddo
!
   end subroutine dgaleq_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine dgaleq_c( a, y, n, ipits)
!  ==================================================================
!
   implicit   none
!
   integer(kind=IntKind), intent(in) :: n, ipits
!
   complex(kind=CmplxKind), intent(inout) :: a(ipits+1,ipits+1)
   complex(kind=CmplxKind), intent(inout) :: y(ipits+1)
!
   integer(kind=IntKind) :: i, j, k, ij
!
   complex(kind=CmplxKind) :: f1,f2
!
!  ==================================================================
!  ******************************************************************
!  Important note: the following codes may will result a problem:
!         a(i-1,i-1) becomes zero!!!
!  It needs to be fixed.
!  ******************************************************************
!  ==================================================================
   do i = 2,n
      f1 = -cone/a(i-1,i-1)
      do j = i,n
         f2 = f1*a(j,i-1)
         do k = 1,n
            a(j,k) = a(j,k)+f2*conjg(a(i-1,k))
         enddo
         y(j) = y(j)+f2*y(i-1)
      enddo
   enddo
   y(n) = y(n)/a(n,n)
   do ij = 1,n-1
      i = n-ij
      do j = 1,ij
         y(i) = y(i)-y(i+j)*conjg(a(i,i+j))
      enddo
      y(i) = y(i)/a(i,i)
   enddo
!
   end subroutine dgaleq_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   function  dgarms_r(fins,fots,i,j,npts)         result(dga_rms)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: i, j, npts
!
   real(kind=RealKind), target :: fins(:,:), fots(:,:)
!
   integer(kind=IntKind) :: k
!
   real(kind=RealKind) :: dga_rms
!
   dga_rms = zero
   do k = 1,npts
      dga_rms = dga_rms + (fins(i,k)-fots(i,k))*(fins(j,k)-fots(j,k))
   enddo
!
   end function dgarms_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   function  dgarms_c(fins,fots,i,j,npts)         result(dga_rms)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: i, j, npts
!
   complex(kind=CmplxKind), target :: fins(:,:), fots(:,:)
!
   integer(kind=IntKind) :: k
!
   complex(kind=CmplxKind) :: dga_rms
!
   dga_rms = czero
   do k = 1,npts
      dga_rms = dga_rms + (fins(i,k)-fots(i,k))*conjg(fins(j,k)-fots(j,k))
   enddo
!
   end function dgarms_c
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function  dgasad_r(ipits,cin,cot,b,niter,amix)       result(dga_sad)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: ipits, niter
!
   real(kind=RealKind), intent(in) :: amix
   real(kind=RealKind), target :: cin(:), cot(:), b(:)
!
   integer(kind=IntKind) :: i, iter
!
   real(kind=RealKind) :: optin, optot
!
   real(kind=RealKind) :: dga_sad
!
   iter = min(niter,ipits)
   optin = zero
   optot = zero
   do i = 1,iter
      optin = optin + b(i)*cin(i)
      optot = optot + b(i)*cot(i)
   enddo
   dga_sad = optin + amix*( optot - optin )
!
   end function dgasad_r
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function  dgasad_c(ipits,cin,cot,b,niter,amix)       result(dga_sad)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: ipits, niter
!
   real(kind=RealKind), intent(in) :: amix
   complex(kind=CmplxKind), target :: cin(:), cot(:), b(:)
!
   integer(kind=IntKind) :: i, iter
!
   complex(kind=CmplxKind) :: optin, optot
!
   complex(kind=CmplxKind) :: dga_sad
!
   iter = min(niter,ipits)
   optin = czero
   optot = czero
   do i = 1,iter
      optin = optin + b(i)*cin(i)
      optot = optot + b(i)*cot(i)
   enddo
   dga_sad = optin + amix*( optot - optin )
!
   end function dgasad_c
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function  getMixID( idt, idq )                          result(id)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: idt, idq
!
   integer(kind=IntKind) :: i
!
   integer(kind=IntKind) :: id
!
   id = 0
   do i = 1,idt-1
      id = id + NumQuantities(idt)
   enddo
   id = id + idq
!
   end function getMixID
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine  setMixingAlpha( idt, idq , amix )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: idq, idt
!
   real(kind=RealKind), intent(in) :: amix 
!
   integer(kind=IntKind) :: id
!
   id = getMixID(idt,idq)
   if ( MixingMethod == 0) then
      SimpleMix(id)%alpha = amix
   else if ( MixingMethod == 1 ) then
      BroydenMix(id)%alpha = amix
   else if ( MixingMethod == 2 ) then
      DGAMix(id)%alpha = amix
   endif
!
   end subroutine setMixingAlpha
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine setSimpleSpace( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   SimpleMix(id)%vlen = size
!
   SimpleMix(id)%Initialized = .true.
!
   end subroutine setSimpleSpace
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setDGASpace_r( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( DGAMix(id)%f_old_r(NumDGAIter,size) )
   allocate( DGAMix(id)%f_new_r(NumDGAIter,size) )
   DGAMix(id)%vlen = size
!
   nullify( DGAMix(id)%f_old_c )
   nullify( DGAMix(id)%f_new_c )
!
   DGAMix(id)%Initialized = .true.
!
   end subroutine setDGASpace_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setDGASpace_c( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( DGAMix(id)%f_old_c(NumDGAIter,size) )
   allocate( DGAMix(id)%f_new_c(NumDGAIter,size) )
   DGAMix(id)%vlen = size
!
   nullify( DGAMix(id)%f_old_r )
   nullify( DGAMix(id)%f_new_r )
!
   DGAMix(id)%Initialized = .true.
!
   end subroutine setDGASpace_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setBroydenSpace_r( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( BroydenMix(id)%u_r(size,NumBroydenIter) )
   allocate( BroydenMix(id)%vt_r(size,NumBroydenIter) )
   allocate( BroydenMix(id)%f_r(size) )
   allocate( BroydenMix(id)%df_r(size) )
   allocate( BroydenMix(id)%w_r(NumBroydenIter) )
   allocate( BroydenMix(id)%vold_r(size) )
   BroydenMix(id)%u_r = ZERO
   BroydenMix(id)%vt_r = ZERO
   BroydenMix(id)%f_r = ZERO
   BroydenMix(id)%df_r = ZERO
   BroydenMix(id)%w_r = ZERO
   BroydenMix(id)%vold_r = ZERO
   BroydenMix(id)%vlen = size
!
   nullify( BroydenMix(id)%u_c )
   nullify( BroydenMix(id)%vt_c )
   nullify( BroydenMix(id)%f_c )
   nullify( BroydenMix(id)%df_c )
   nullify( BroydenMix(id)%w_c )
   nullify( BroydenMix(id)%vold_c )
!
   BroydenMix(id)%Initialized = .true.
!
   end subroutine setBroydenSpace_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setBroydenSpace_c( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( BroydenMix(id)%u_c(size,NumBroydenIter) )
   allocate( BroydenMix(id)%vt_c(size,NumBroydenIter) )
   allocate( BroydenMix(id)%f_c(size) )
   allocate( BroydenMix(id)%df_c(size) )
   allocate( BroydenMix(id)%w_c(NumBroydenIter) )
   allocate( BroydenMix(id)%vold_c(size) )
   BroydenMix(id)%u_c = CZERO
   BroydenMix(id)%vt_c = CZERO
   BroydenMix(id)%f_c = CZERO
   BroydenMix(id)%df_c = CZERO
   BroydenMix(id)%w_c = CZERO
   BroydenMix(id)%vold_c = CZERO
   BroydenMix(id)%vlen = size
!
   nullify( BroydenMix(id)%u_r )
   nullify( BroydenMix(id)%vt_r )
   nullify( BroydenMix(id)%f_r )
   nullify( BroydenMix(id)%df_r )
   nullify( BroydenMix(id)%w_r )
   nullify( BroydenMix(id)%vold_r )
!
   BroydenMix(id)%Initialized = .true.
!
   end subroutine setBroydenSpace_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setWorkingSpace_r( MixId )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: MixId 
!
   if (MixId == 1) then
      allocate( a_r( NumBroydenIter,NumBroydenIter ) )
      allocate( b_r( NumBroydenIter,NumBroydenIter ) )
      allocate( d_r( NumBroydenIter,NumBroydenIter ) )
      allocate( cm_r( NumBroydenIter ) )
      a_r = ZERO; b_r = ZERO; d_r = ZERO; cm_r = ZERO
!
      if ( BroydenInvMethod == 1 ) then
         allocate( ipiv(NumBroydenIter) )
      endif
   else if (MixId == 2) then
      allocate( a_r(NumDGAIter+1,NumDGAIter+1), b_r(NumDGAIter+1,NumDGAIter+1) )
      a_r = ZERO; b_r = ZERO
   endif
!
   InitializedWkSpace = .true.
!
   end subroutine setWorkingSpace_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setWorkingSpace_c( MixId )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: MixId
!
   if (MixId == 1) then
      allocate( a_c( NumBroydenIter,NumBroydenIter ) )
      allocate( b_c( NumBroydenIter,NumBroydenIter ) )
      allocate( d_c( NumBroydenIter,NumBroydenIter ) )
      allocate( cm_c( NumBroydenIter ) )
      a_c = CZERO; b_c = CZERO; d_c = CZERO; cm_c = CZERO
!
      if ( BroydenInvMethod == 1 ) then
         allocate( ipiv(NumBroydenIter) )
      endif
   else if (MixId == 2) then
      allocate( a_c(NumDGAIter+1,NumDGAIter+1), b_c(NumDGAIter+1,NumDGAIter+1) )
      a_c = CZERO; b_c = CZERO
   endif
!
   InitializedWkSpace = .true.
!
   end subroutine setWorkingSpace_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  delWorkingSpace()
!  ==================================================================
!
   implicit  none
!
   if ( MixingMethod == 1 ) then
      if ( allocated(a_r) ) deallocate( a_r, b_r, d_r, cm_r )
      if ( allocated(a_c) ) deallocate( a_c, b_c, d_c, cm_c )
!
      if ( BroydenInvMethod == 1 ) then
         if ( allocated(ipiv) ) deallocate( ipiv )
      endif
   else if (MixingMethod == 2) then
      if ( allocated(a_r) ) deallocate( a_r, b_r )
      if ( allocated(a_c) ) deallocate( a_c, b_c )
   endif
!
   InitializedWkSpace = .false.
!
   end subroutine delWorkingSpace
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  delMixing()
!  ==================================================================
!
   implicit  none
!
   call delWorkingSpace()
!
   if ( MixingMethod == 0 ) then
      return
   else if ( MixingMethod == 1 ) then
      call delBroydenMixing()
   else
      call delDGAMixing()
   endif
!
   end subroutine delMixing
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  resetMixing_r(pMixingList)
!  ==================================================================
!
   implicit none
!
   type (MixListRealStruct), target:: pMixingList
!
   type (MixListRealStruct), pointer:: pML
!
   pML => pMixingList
   do while ( associated(pML) )
      pML%size = 0
      pML%rms = ZERO
      nullify(pML%mesh)
      nullify(pML%vector_old)
      nullify(pML%vector_new)
      pML => pML%next
   end do
!
   call delWorkingSpace()
   call delMixing()
!
   end subroutine resetMixing_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  resetMixing_c(pMixingList)
!  ==================================================================
!
   implicit none
!
   type (MixListCmplxStruct), target:: pMixingList
!
   type (MixListCmplxStruct), pointer:: pML
!
   pML => pMixingList
   do while ( associated(pML) )
      pML%size = 0
      pML%rms = ZERO
      nullify(pML%mesh)
      nullify(pML%vector_old)
      nullify(pML%vector_new)
      pML => pML%next
   end do
!
   call delWorkingSpace()
   call delMixing()
   iter_count = 0
!
   end subroutine resetMixing_c
!  ==================================================================
end module MixingModule
