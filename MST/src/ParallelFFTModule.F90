module ParallelFFTModule
!  This module hides specific parallel FFT implementation and
!  presents a generic interface for performing parallel FFT.
! 
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
#ifdef P3DFFT
   use p3dfft
#else
   include 'fftw3.f'
   integer (kind=8), private :: plan
!  use iso_c_binding
!  include 'fftw3.f03'
!  type (C_PTR), private :: plan
#endif
!
public ::                  &
   initParallelFFT,        &
   endParallelFFT,         &
   getProcessorMesh,       &
   getNumGlobalGridPoints, & ! returns the number of overall grid points
   getNumLocalGridPoints,  & ! returns the number of the grid points allocated on my processor
   getGlobalGridIndexRange,& ! returns the global index of the grid points allocated on my processor
   getGridStep,            & ! returns the step vector of the grid points along a particular dimension
   getLocalArraySize,      & ! returns the size of the array to be allocated for storing real space data
   getGridPointCoord,      & ! returns the coordinates of a grid point
   getLocalIndexOfGridOrigin, &
   allocateFunctionSpace,  &
   aliasFunctionSpace,     &
   performTransformR2C,    & ! The returned result already includes a factor of 1/Ng
   getParaFFTCommunicator
!
   interface getGridPointCoord
      module procedure getGridPointCoord_i1, getGridPointCoord_i3
   end interface getGridPointCoord
!
   interface allocateFunctionSpace
      module procedure allocateFunctionSpace_R, allocateFunctionSpace_K
   end interface allocateFunctionSpace
!
   interface aliasFunctionSpace
      module procedure aliasFunctionSpace_R, aliasFunctionSpace_K
   end interface aliasFunctionSpace
!
   interface getGlobalGridIndexRange
      module procedure getGlobalGridIndexRange_a, getGlobalGridIndexRange_d
   end interface getGlobalGridIndexRange
!
private
   integer (kind=IntKind) :: na, nb, nc
   integer (kind=IntKind) :: istart(3), iend(3), isize(3)
   integer (kind=IntKind) :: fstart(3), fend(3), fsize(3)
   integer (kind=IntKind) :: tstart(3), tend(3), memsize(3)
   integer (kind=IntKind) :: num_total_points, numr_local, numk_local, numka
   integer (kind=IntKind) :: a_fe, a_fs
   integer (kind=IntKind) :: rgrid0(3), kgrid0(3), rgrid0_id, kgrid0_id
   integer (kind=IntKind) :: MyPE_comm, NumPEs_comm, comm
   integer (kind=IntKind) :: processor_mesh(3)
!
!  ===================================================================
!  rstep will be used to construct the uniform mesh points in the real space.
!  kstep will be used to construct the uniform mesh points in the reciprocal space.
!  ===================================================================
   real (kind=RealKind) :: rstep(3,3)
   real (kind=RealKind) :: kstep(3,3) ! kstep contain a factor of PI2
!
   logical :: Initialized = .false.
!
   character (len=6) :: MeshType  ! This depends on the parallel FFT to be used.
                                  ! It can be 'pencil', 'slab', or 'cube'
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initParallelFFT()
!  ===================================================================
   use PublicTypeDefinitionsModule, only : UniformGridStruct
!
   use VectorModule, only : getVecLength
!
   use MathParamModule, only : PI2, TEN2m8
!
   use MPPModule, only : getCommunicator, GlobalSum, getMyPE, getNumPEs
!
   use PrimeFactorsModule, only : getSubFactors
!
   use Uniform3DGridModule, only : getUniform3DGrid
!
   implicit none
!
   integer (kind=IntKind) :: dims(2), id(2)
   integer (kind=IntKind) :: idk, ig, jg, kg
   integer (kind=IntKind) :: i, j, k, iproc, jproc, m, nproc
   integer (kind=IntKind), pointer :: factors(:,:)
!
   logical :: found
!  
   real (kind=RealKind) :: vol, p1(3), p2(3)
!
   type (UniformGridStruct), pointer :: pUG
!
   comm = getCommunicator()
!  ===================================================================
!  The following lines are commented out since IBM MPI library appears
!  giving negative comm value. So, this check is not reliable.
!  ===================================================================
!  if (comm <= 0) then
!     call ErrorHandler('initParallelFFT','comm <= 0',                &
!                       'Code PotentialGenerationModule needs to be redesigned')
!  endif
!  ===================================================================
#ifdef P3DFFT
   MyPE_comm = getMyPE()
   NumPEs_comm = getNumPEs()
#else
   comm = -1
   MyPE_comm = 0
   NumPEs_comm = 1
#endif
!
   pUG => getUniform3DGrid('FFT')
   na = pUG%nga
   nb = pUG%ngb
   nc = pUG%ngc
   rstep(1:3,1) = pUG%grid_step_a(1:3)
   rstep(1:3,2) = pUG%grid_step_b(1:3)
   rstep(1:3,3) = pUG%grid_step_c(1:3)
   num_total_points = pUG%ng
!
   kstep(1,1) = (pUG%cell(2,2)*pUG%cell(3,3)-pUG%cell(3,2)*pUG%cell(2,3))
   kstep(2,1) = (pUG%cell(3,2)*pUG%cell(1,3)-pUG%cell(1,2)*pUG%cell(3,3))
   kstep(3,1) = (pUG%cell(1,2)*pUG%cell(2,3)-pUG%cell(2,2)*pUG%cell(1,3))
   kstep(1,2) = (pUG%cell(2,3)*pUG%cell(3,1)-pUG%cell(3,3)*pUG%cell(2,1))
   kstep(2,2) = (pUG%cell(3,3)*pUG%cell(1,1)-pUG%cell(1,3)*pUG%cell(3,1))
   kstep(3,2) = (pUG%cell(1,3)*pUG%cell(2,1)-pUG%cell(2,3)*pUG%cell(1,1))
   kstep(1,3) = (pUG%cell(2,1)*pUG%cell(3,2)-pUG%cell(3,1)*pUG%cell(2,2))
   kstep(2,3) = (pUG%cell(3,1)*pUG%cell(1,2)-pUG%cell(1,1)*pUG%cell(3,2))
   kstep(3,3) = (pUG%cell(1,1)*pUG%cell(2,2)-pUG%cell(2,1)*pUG%cell(1,2))
!
   vol = (pUG%cell(2,1)*pUG%cell(3,2)-pUG%cell(3,1)*pUG%cell(2,2))*pUG%cell(1,3)+    &
         (pUG%cell(3,1)*pUG%cell(1,2)-pUG%cell(1,1)*pUG%cell(3,2))*pUG%cell(2,3)+    &
         (pUG%cell(1,1)*pUG%cell(2,2)-pUG%cell(2,1)*pUG%cell(1,2))*pUG%cell(3,3)
   kstep = PI2/vol*kstep
!
#ifdef P3DFFT
   MeshType = 'pencil'
!  ===================================================================
!  nproc = the number of processors in the grid distribution mesh on
!          which the grid points are distributed.
!  ===================================================================
   if (pUG%ngb*pUG%ngc < NumPEs_comm) then
      nproc = pUG%ngb*pUG%ngc
   else
      nproc = NumPEs_comm - mod(pUG%ngb*pUG%ngc,NumPEs_comm)
   endif
!
   if (nproc > 1) then
!     ================================================================
!     nproc is divided into iproc x jproc stencles
!     This is a naive scheme, may needs optimization.
!     ================================================================
      factors => getSubFactors(nproc,2,m) ! break the total number of
                                          ! processes into two factors
!!!   k = pUG%ngb + pUG%ngc + 1
      iproc = factors(1,1)
      jproc = factors(2,1)
      found = .false.
      do i = 1, m
         if (mod(pUG%ngb,factors(1,i)) == 0 .and.                     &
             mod(pUG%ngc,factors(2,i)) == 0 .and.                     &
             mod(pUG%ngb,factors(2,i)) == 0 .and.                     &
             mod(pUG%nga/2,factors(1,i)) == 0) then
!!!         j = (pUG%nga/2+pUG%ngb)/factors(1,i) + (pUG%ngb+pUG%ngc)/factors(2,i)
!!!         if (j < k .and. factors(1,i) <= factors(2,i)) then
            if (factors(1,i) <= factors(2,i)) then
               iproc = factors(1,i)
               jproc = factors(2,i)
               found = .true.
!!!            k = j
            endif
         endif
      enddo
!
      if (.not.found) then
         call ErrorHandler('initParallelFFT',                         &
                           'A processor mesh is not established for the given grid')
      else if (MyPE_comm == 0) then
         write(6,'(/,a)')'*****************************************************'
         write(6,'(a)')  '*   Print out from subroutine createProcessorMesh   *'
         write(6,'(a,/)')'*****************************************************'
         write(6,'(a,i5)')'No. processors in the processor mesh for grid distribution (PMGD) = ',nproc
         write(6,'(a,i5)')'No. processors in the b-dimention of the PMGD = ',iproc
         write(6,'(a,i5)')'No. processors in the c-dimention of the PMGD = ',jproc
         if (NumPEs_comm-nproc > 0) then
            write(6,'(a,i5)')'No. processors not belonging to the PMGD      = ',NumPEs_comm-nproc
         endif
      endif
   else
      iproc = 1
      jproc = 1
   endif
!
   dims(1) = iproc
   dims(2) = jproc
!  -------------------------------------------------------------------
   call p3dfft_setup(dims,na,nb,nc,comm)
   call p3dfft_get_dims(istart,iend,isize,1)
   call p3dfft_get_dims(fstart,fend,fsize,2)
   call p3dfft_get_dims(tstart,tend,memsize,3)
!  -------------------------------------------------------------------
   processor_mesh(1) = 1
   processor_mesh(2) = iproc
   processor_mesh(3) = jproc
!
#else
   MeshType = 'none'
   processor_mesh = 1
   istart = 1
   iend(1) = na; iend(2) = nb; iend(3) = nc
   isize = iend - istart + 1
   memsize(1) = isize(1)+2; memsize(2) = isize(2); memsize(3) = isize(3)
   fstart = istart
   fend(1) = memsize(1)/2; fend(2) = memsize(2); fend(3) = memsize(3)
   fsize(1) = (fend(1)-fstart(1)+1); fsize(2) = memsize(2); fsize(3) = memsize(3)
#endif
!
   if (istart(1) < 1 .or. istart(2) < 1 .or. istart(3) < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'Invalid starting grid index of real space',  &
                        istart(1),istart(2),istart(3))
!     ----------------------------------------------------------------
   else if (iend(1) > na .or. iend(2) > nb .or. iend(3) > nc) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'Invalid ending grid index of real space',    &
                        iend(1),iend(2),iend(3))
!     ----------------------------------------------------------------
   else if (fstart(1) < 1 .or. fstart(2) < 1 .or. fstart(3) < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'Invalid starting grid index of reciprocal space', &
                        fstart(1),fstart(2),fstart(3))
!     ----------------------------------------------------------------
   else if (fend(1) > na .or. fend(2) > nb .or. fend(3) > nc) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'Invalid ending grid index of reciprocal space', &
                        fend(1),fend(2),fend(3))
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  The real space grid points mapped on the current processor are those with
!  global indecies (along a, b, and c dimensions) from istart(1:3) to iend(1:3),
!  with isize = iend - istart + 1. The total number of real space grid 
!  points mapped onto the current processor is numr_local.
!  The function to be Fourier transformed will be calculated at each real 
!  space grid point.
!  Note that isize(1), isize(2), and isize(3) are the number of real space
!  grid point indeces along a, b, and c dimensions, respectively.
!  ===================================================================
   numr_local = isize(1)*isize(2)*isize(3)
!
!  ===================================================================
!  The reciprocal space grid points are also distributed over the employed
!  processors. 
!  Since we intend to use in-place FFT to transform a real function from
!  the real space to the reciprocal space, the resulting complex function 
!  needs to be carefully restored on all the grid points.
!  The total number of reciprocal space points mapped on the current processor
!  is numk_local.
!  Note that, due to transpose and shuffling of the Fourier transformed results
!  among the processors employed to perform parallel FFT, numk_local is 
!  not necessarily the same as numr_local. 
!  Of course, the sum of numk_local across the processors is the same as
!  the sum of numr_local across the processors.
!  Note that numka, fsize(2), and fsize(3) are the number of reciprocal
!  space grid point indeces along a, b, and c dimensions, respectively.
!  ===================================================================
   if (fend(1) < na/2+1) then
      a_fe = fend(1)
   else
      a_fe = fend(1) - 1
   endif
   if (fstart(1) == 1) then
      a_fs = 2
   else
      a_fs = fstart(1)
   endif
   numka = a_fe-a_fs+1+fsize(1)
   numk_local = numka*fsize(2)*fsize(3)
!
!  ===================================================================
!  Determine the local index of the grid point at (0,0,0) and store the index
!  information in rgrid0(1:3) and rgrid0_id.
!  If the origin is not mapped on the local processor, set rgrid0(1:3) = 0 and rgrid_id = 0
!  ===================================================================
   rgrid0 = 0; rgrid0_id = 0
   if (istart(1) == 1 .and. istart(2) == 1 .and. istart(3) == 1) then
      rgrid0 = 1
      rgrid0_id = 1
   endif
!  ===================================================================
!  Check the consistency of the grid index data
!  ===================================================================
   p1 = getGridPointCoord_i3('R',1,1,1)
   p2 = getGridPointCoord_i1('R',1)
   if (abs(p1(1)) < TEN2m8 .and. abs(p1(2)) < TEN2m8 .and. abs(p1(3)) < TEN2m8) then
      if (abs(p2(1)) > TEN2m8 .or. abs(p2(2)) > TEN2m8 .or. abs(p2(3)) > TEN2m8) then
!        -------------------------------------------------------------
         call ErrorHandler('initParallelFFT',                         &
                           'Inconsistent index of real space grid point origin', &
                           'p1 <> p2')
!        -------------------------------------------------------------
      else if (rgrid0_id /= 1) then
!        -------------------------------------------------------------
         call ErrorHandler('initParallelFFT',                         &
                           'Inconsistent index of real space grid point origin', &
                           'p1 = 0 but rgrid0_id <> 1')
!        -------------------------------------------------------------
      endif
   endif
   kgrid0 = 0; kgrid0_id = 0
   if (fstart(1) == 1 .and. fstart(2) == 1 .and. fstart(3) == 1) then
      kgrid0 = 1
      kgrid0_id = 1
   endif
   p1 = getGridPointCoord_i3('K',1,1,1)
   p2 = getGridPointCoord_i1('K',1)
   if (abs(p1(1)) < TEN2m8 .and. abs(p1(2)) < TEN2m8 .and. abs(p1(3)) < TEN2m8) then
      if (abs(p2(1)) > TEN2m8 .or. abs(p2(2)) > TEN2m8 .or. abs(p2(3)) > TEN2m8) then
!        -------------------------------------------------------------
         call ErrorHandler('initParallelFFT',                         &
                           'Inconsistent index of reciprocal space grid point origin', &
                           'p1 <> p2')
!        -------------------------------------------------------------
      else if (kgrid0_id /= 1) then
!        -------------------------------------------------------------
         call ErrorHandler('initParallelFFT',                         &
                           'Inconsistent index of reciprocal space grid point origin', &
                           'p1 = 0 but kgrid0_id <> 1')
!        -------------------------------------------------------------
      endif
   endif
   id(1) = rgrid0_id; id(2) = kgrid0_id
#ifdef P3DFFT
!  -------------------------------------------------------------------
   call GlobalSum(id,2)
!  -------------------------------------------------------------------
#endif
   if (id(1) == 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'No (0,0,0) grid point in real space exists')
!     ----------------------------------------------------------------
   else if (id(2) == 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'No (0,0,0) grid point in reciprocal space exists')
!     ----------------------------------------------------------------
   else if (id(1) > 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'Real space (0,0,0) grid point exists on multiple processors',id(1))
!     ----------------------------------------------------------------
   else if (id(2) > 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('initParallelFFT',                            &
                        'Reciprocal space (0,0,0) grid point exists on multiple processors',id(2))
!     ----------------------------------------------------------------
   endif
!
   Initialized = .true.
!
   end subroutine initParallelFFT
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endParallelFFT()
!  ===================================================================
   implicit none
!
#ifdef P3DFFT
   call p3dfft_clean()
#endif
   Initialized = .false.
!
   end subroutine endParallelFFT
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGlobalGridPoints(dim) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: dim
   integer (kind=IntKind) :: n
!
   if (present(dim)) then
      if (dim < 1 .or. dim > 3) then
!        -------------------------------------------------------------
         call ErrorHandler('getNumGlobalGridPoints','Invalid dimension index',dim)
!        -------------------------------------------------------------
      endif
   endif
!
   if (present(dim)) then
      if (dim == 1) then
         n = na
      else if (dim == 2) then
         n = nb
      else
         n = nc
      endif
   else
      n = num_total_points
   endif
!
   end function getNumGlobalGridPoints
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getProcessorMesh() result (pm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: pm(3)
!
   pm = processor_mesh
!
   end function getProcessorMesh
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumLocalGridPoints(c,dim) result(n)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in), optional :: dim
   integer (kind=IntKind) :: n
!
   if (present(dim)) then
      if (dim < 1 .or. dim > 3) then
!        -------------------------------------------------------------
         call ErrorHandler('getNumLocalGridPoints','Invalid dimension index',dim)
!        -------------------------------------------------------------
      endif
   endif
!
   if (c == 'r' .or. c == 'R') then
      if (present(dim)) then
         n = isize(dim)
      else
         n = isize(1)*isize(2)*isize(3)
      endif
   else if (c == 'k' .or. c == 'K') then
      if (present(dim)) then
         if (dim == 1) then
            n = numka
         else
            n = fsize(dim)
         endif
      else
         n = numk_local
      endif
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getNumLocalGridPoints','Unrecognized space indicator',c)
!     ----------------------------------------------------------------
   endif
!
   end function getNumLocalGridPoints
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalGridIndexRange_a(c) result(i)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind) :: i(3,3)
   integer (kind=IntKind) :: dim
!
   if (.not.Initialized) then
      call ErrorHandler('getGlobalGridIndexRange','Module is not initialized', &
                        .true.)
   endif
!
   if (c == 'r' .or. c == 'R') then
      do dim = 1, 3
         i(dim,1) = istart(dim)
         i(dim,2) = iend(dim)
         i(dim,3) = isize(dim)
      enddo
   else if (c == 'k' .or. c == 'K') then
      do dim = 1, 3
         i(dim,1) = fstart(dim)
         i(dim,2) = fend(dim)
         if (dim == 1) then
            i(dim,3) = numka
         else
            i(dim,3) = fsize(dim)
         endif
      enddo
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getGlobalGridIndexRange','Unrecognized space indicator', &
                        c, .true.)
!     ----------------------------------------------------------------
   endif
!
   end function getGlobalGridIndexRange_a
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGlobalGridIndexRange_d(c,dim) result(i)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in) :: dim
   integer (kind=IntKind) :: i(3)
!
   if (.not.Initialized) then
      call ErrorHandler('getGlobalGridIndexRange','Need to call initParallelFFT first', &
                        .true.)
   else if (dim < 1 .or. dim > 3) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGlobalGridIndexRange','Invalid dimension index',dim,.true.)
!     ----------------------------------------------------------------
   endif
!
   if (c == 'r' .or. c == 'R') then
      i(1) = istart(dim)
      i(2) = iend(dim)
      i(3) = isize(dim)
   else if (c == 'k' .or. c == 'K') then
      i(1) = fstart(dim)
      i(2) = fend(dim)
      if (dim == 1) then
         i(3) = numka
      else
         i(3) = fsize(dim)
      endif
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getGlobalGridIndexRange','Unrecognized space indicator',c,.true.)
!     ----------------------------------------------------------------
   endif
!
   end function getGlobalGridIndexRange_d
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridStep(c,dim) result(s)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in) :: dim
!
   real (kind=RealKind) :: s(3)
!
   if (dim < 1 .or. dim > 3) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGridStep','Invalid dimension index',dim)
!     ----------------------------------------------------------------
   endif
!
   if (c == 'r' .or. c == 'R') then
      s(1:3) = rstep(1:3,dim)
   else if (c == 'k' .or. c == 'K') then
      s(1:3) = kstep(1:3,dim)
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getGridStep','Unrecognized space indicator',c)
!     ----------------------------------------------------------------
   endif
!
   end function getGridStep
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridPointCoord_i1(c,idk) result(p)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in) :: idk
   integer (kind=IntKind) :: i, j, k, n12, ig, jg, kg, ka, kb, kc
!
   real (kind=RealKind) :: p(3)
!
   if (c == 'r' .or. c == 'R') then
      if (idk < 1 .or. idk > numr_local) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index out of range',idk)
!        -------------------------------------------------------------
      endif
      n12 = isize(1)*isize(2)
      i = mod(idk-1,isize(1))+1
      j = (mod(idk-1,n12)+1-i)/isize(1)+1
      k = (idk-i-(j-1)*isize(1))/n12+1
      ig = i+istart(1)-1
      jg = j+istart(2)-1
      kg = k+istart(3)-1
      p(1:3) = (ig-1)*rstep(1:3,1)+(jg-1)*rstep(1:3,2)+(kg-1)*rstep(1:3,3)
   else if (c == 'k' .or. c == 'K') then
      if (idk < 1 .or. idk > numk_local) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index out of range',idk)
!        -------------------------------------------------------------
      endif
      n12 = numka*fsize(2)
      i = mod(idk-1,numka)+1
      j = (mod(idk-1,n12)+1-i)/numka+1
      k = (idk-i-(j-1)*numka)/n12+1
      jg = j+fstart(2)-1
      kg = k+fstart(3)-1
      if (i <= fsize(1)) then
         ig = i+fstart(1)-1
         if (ig < na/2+1 .or. na == 1) then
            ka = ig-1
         else
            ka = ig-1-na  ! negative ka
         endif
!        =============================================================
!        determine kb, kc that are mapped on the current processor
!        =============================================================
         if (jg < nb/2+1 .or. nb == 1) then
            kb = jg-1
         else
            kb = jg-1-nb  ! negative kb
         endif
         if (kg < nc/2+1 .or. nc == 1) then
            kc = kg-1
         else
            kc = kg-1-nc  ! negative kc
         endif
      else
         ig = a_fe-(i-fsize(1))+1
         ka = -(ig-1)      ! negative ka
!        =============================================================
!        determine kb, kc that are mapped on the current processor
!        =============================================================
         if (jg <= nb/2+1) then
            kb = -(jg-1)   ! negative kb
         else
            kb = nb-jg+1
         endif
         if (kg <= nc/2+1) then
            kc = -(kg-1)   ! negative kc
         else
            kc = nc-kg+1
         endif
      endif
      p(1:3) = ka*kstep(1:3,1)+kb*kstep(1:3,2)+kc*kstep(1:3,3)
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getGridPointCoord','Unrecognized space indicator',c)
!     ----------------------------------------------------------------
   endif
!
   end function getGridPointCoord_i1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridPointCoord_i3(c,i,j,k) result(p)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in) :: i, j, k
   integer (kind=IntKind) :: n12, ig, jg, kg, ka, kb, kc
!
   real (kind=RealKind) :: p(3)
!
   if (c == 'r' .or. c == 'R') then
      if (i < 1 .or. i > isize(1)) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index i out of range',i)
!        -------------------------------------------------------------
      else if (j < 1 .or. j > isize(2)) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index j out of range',j)
!        -------------------------------------------------------------
      else if (k < 1 .or. k > isize(3)) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index k out of range',k)
!        -------------------------------------------------------------
      endif
      ig = i+istart(1)-1
      jg = j+istart(2)-1
      kg = k+istart(3)-1
      p(1:3) = (ig-1)*rstep(1:3,1)+(jg-1)*rstep(1:3,2)+(kg-1)*rstep(1:3,3)
   else if (c == 'k' .or. c == 'K') then
      if (i < 1 .or. i > numka) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index i out of range',i)
!        -------------------------------------------------------------
      else if (j < 1 .or. j > fsize(2)) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index j out of range',j)
!        -------------------------------------------------------------
      else if (k < 1 .or. k > fsize(3)) then
!        -------------------------------------------------------------
         call ErrorHandler('getGridPointCoord','Local grid point index k out of range',k)
!        -------------------------------------------------------------
      endif
      jg = j+fstart(2)-1
      kg = k+fstart(3)-1
      if (i <= fsize(1)) then
         ig = i+fstart(1)-1
         if (ig < na/2+1 .or. na == 1) then
            ka = ig-1
         else
            ka = ig-1-na  ! negative ka
         endif
!        =============================================================
!        determine kb, kc that are mapped on the current processor
!        =============================================================
         if (jg < nb/2+1 .or. nb == 1) then
            kb = jg-1
         else
            kb = jg-1-nb  ! negative kb
         endif
         if (kg < nc/2+1 .or. nc == 1) then
            kc = kg-1
         else
            kc = kg-1-nc  ! negative kc
         endif
      else
         ig = a_fe-(i-fsize(1))+1
         ka = -(ig-1)      ! negative ka
!        =============================================================
!        determine kb, kc that are mapped on the current processor
!        =============================================================
         if (jg <= nb/2+1) then
            kb = -(jg-1)   ! negative kb
         else
            kb = nb-jg+1
         endif
         if (kg <= nc/2+1) then
            kc = -(kg-1)   ! negative kc
         else
            kc = nc-kg+1
         endif
      endif
      p(1:3) = ka*kstep(1:3,1)+kb*kstep(1:3,2)+kc*kstep(1:3,3)
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getGridPointCoord','Unrecognized space indicator',c)
!     ----------------------------------------------------------------
   endif
!
   end function getGridPointCoord_i3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalArraySize(c,dim) result(mem)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in), optional :: dim
   integer (kind=IntKind) :: mem
!
   if (present(dim)) then
      if (dim < 1 .or. dim > 3) then
!        -------------------------------------------------------------
         call ErrorHandler('getLocalArraySize','Invalid dimension index',dim)
!        -------------------------------------------------------------
      endif
   endif
!
   if (c == 'r' .or. c == 'R') then
      if (present(dim)) then
         mem = memsize(dim)
      else
         mem = memsize(1)*memsize(2)*memsize(3)
      endif
   else if (c == 'k' .or. c == 'K') then
      if (present(dim)) then
         if (dim == 1) then
            mem = numka
         else 
            mem = fsize(dim)
         endif
      else
         mem = numk_local
      endif
   else
!     ----------------------------------------------------------------
      call ErrorHandler('getLocalArraySize','Unrecognized space indicator',c)
!     ----------------------------------------------------------------
   endif
!
   end function getLocalArraySize
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLocalIndexOfGridOrigin(c,dim) result(i0)
!  ===================================================================
   implicit none
!
   character (len=1), intent(in) :: c
!
   integer (kind=IntKind), intent(in), optional :: dim
   integer (kind=IntKind) :: i0
!
   if (present(dim)) then
      if (dim < 1 .or. dim > 3) then
!        -------------------------------------------------------------
         call ErrorHandler('getLocalIndexOfGridOrigin','Invalid dimension index',dim)
!        -------------------------------------------------------------
      endif
   endif
!
   if (c == 'r' .or. c == 'R') then
      if (present(dim)) then
         i0 = rgrid0(dim)
      else
         i0 = rgrid0_id
      endif
   else if (c == 'k' .or. c == 'K') then
      if (present(dim)) then
         i0 = kgrid0(dim)
      else
         i0 = kgrid0_id
      endif
   endif
!
   end function getLocalIndexOfGridOrigin
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateFunctionSpace_R(rfunc,n)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind), intent(out) :: n
   real (kind=RealKind), pointer :: rfunc(:)
!
   n = memsize(1)*memsize(2)*memsize(3)
   allocate( rfunc(n) )
   rfunc = ZERO
!
   end subroutine allocateFunctionSpace_R
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateFunctionSpace_K(cfunc,n)
!  ===================================================================
   use MathParamModule, only : CZERO
   implicit none
!
   integer (kind=IntKind), intent(out) :: n
   complex (kind=CmplxKind), pointer :: cfunc(:)
!
   n = numk_local
   allocate( cfunc(n) )
   cfunc = CZERO
!
   end subroutine allocateFunctionSpace_K
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasFunctionSpace_R(rfunc,n1,n2,n3) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: n1, n2, n3
   real (kind=RealKind), target :: rfunc(:)
   real (kind=RealKind), pointer :: p(:,:,:)
!
   n1 = memsize(1)
   n2 = memsize(2)
   n3 = memsize(3)
   p => aliasArray3_r(rfunc,n1,n2,n3)
!
   end function aliasFunctionSpace_R
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasFunctionSpace_K(cfunc,n1,n2,n3) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: n1, n2, n3
   complex (kind=CmplxKind), target :: cfunc(:)
   complex (kind=CmplxKind), pointer :: p(:,:,:)
!
   n1 = numka
   n2 = fsize(2)
   n3 = fsize(3)
   p => aliasArray3_c(cfunc,n1,n2,n3)
!
   end function aliasFunctionSpace_K
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine performTransformR2C(r,c_ip)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(inout) :: r(:)
   complex (kind=CmplxKind), intent(out) :: c_ip(:) ! result of in-place FFT
   integer (kind=IntKind) :: idk, idm, ig, jg, kg, i, j, k
!
   real (kind=RealKind), pointer :: r3(:,:,:)
!
#ifndef P3DFFT
   real (kind=RealKind), allocatable ::  r1(:)
!  type (C_PTR) :: c_func_fftwr
!  real (kind=RealKind), pointer ::  r1(:)
!  complex (kind=CmplxKind), pointer ::  c1(:)
#endif
!
   if (size(r) < memsize(1)*memsize(2)*memsize(3)) then
!     ----------------------------------------------------------------
      call ErrorHandler('performTransformR2C',                        &
                        'Allocated space for the real function is insufficient', &
                        size(r),memsize(1)*memsize(2)*memsize(3))
!     ----------------------------------------------------------------
   else if (size(c_ip) < numk_local) then
!     ----------------------------------------------------------------
      call ErrorHandler('performTransformR2C',                        &
                        'Allocated space for the complex function is insufficient', &
                        size(c_ip),numk_local)
!     ----------------------------------------------------------------
   endif
!
#ifdef P3DFFT
   r3 => aliasArray3_r(r,memsize(1),memsize(2),memsize(3))
!  -------------------------------------------------------------------
   call ftran_r2c(r3,r3,'fft')
!  -------------------------------------------------------------------
   r = r/real(num_total_points,kind=RealKind)
#else
   r3 => aliasArray3_r(r,isize(1),isize(2),isize(3))
   allocate( r1(memsize(1)*memsize(2)*memsize(3)) )
!! c_func_fftwr = fftw_alloc_real((c_nx+2)*c_ny*c_nz)
!  c_func_fftwr = fftw_alloc_real(memsize(1)*memsize(2)*memsize(3))
!  call c_f_pointer(c_func_fftwr, r1, [memsize(1)*memsize(2)*memsize(3)])
!  call c_f_pointer(c_func_fftwr, c1, [memsize(1)*memsize(2)*memsize(3))/2])
   do k = 1, isize(3)
      do j = 1, isize(2)
         do i = 1, isize(1)
            idm = i + (j-1)*memsize(1) + (k-1)*isize(2)*memsize(1)
            r1(idm) = r3(i, j, k)
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call dfftw_plan_dft_r2c_3d( plan, isize(1), isize(2), isize(3), r1, r1, &
                               FFTW_ESTIMATE )
   call dfftw_execute( plan )
   call dfftw_destroy_plan( plan )
!  plan = fftw_plan_dft_r2c_3d( isize(1), isize(2), isize(3), r1, c1, &
!                               FFTW_ESTIMATE )
!  call fftw_execute_dft_r2c( plan, r1, c1 )
!  call fftw_destroy_plan( plan )
!  -------------------------------------------------------------------
   r = r1/real(num_total_points,kind=RealKind)
   deallocate(r1)
!  nullify(r1, c1)
!  call fftw_free(c_func_fftwr)
#endif
!
   idk = 0
   do kg = fstart(3), fend(3)
      k = kg-fstart(3)+1
      do jg = fstart(2), fend(2)
         j = jg-fstart(2)+1
!
!        =============================================================
!        Loop over 1, 2, ..., Na/2+1 values of ka
!        since ka's are distributed, each processor loops over the ka's 
!        mapped on it.
!        =============================================================
         do ig = fstart(1), fend(1) ! numka+fstart(1)-1
            i = ig-fstart(1)+1
            idk = idk + 1
            idm = 2*i+(j-1)*2*fsize(1)+(k-1)*fsize(2)*2*fsize(1)
!           c_ip(idk) = cmplx(r(idm-1), -r(idm), kind=CmplxKind)
            c_ip(idk) = cmplx(r(idm-1), r(idm), kind=CmplxKind)  ! 09/26/18
         enddo
!
!        =============================================================
!        Loop over Na/2+2, Na/2+3, ..., Na values of ka, which are negative
!        since ka's are distributed, each processor loops over the ka's 
!        that those corresponding -ka (>0) are mapped on it.
!        =============================================================
         do ig = a_fe, a_fs, -1
            i = ig-fstart(1)+1
            idk = idk + 1
            idm = 2*i+(j-1)*2*fsize(1)+(k-1)*fsize(2)*2*fsize(1)
!           c_ip(idk) = cmplx(r(idm-1), r(idm), kind=CmplxKind)
            c_ip(idk) = cmplx(r(idm-1), -r(idm), kind=CmplxKind)  ! 09/26/18
         enddo
      enddo
   enddo
!
   end subroutine performTransformR2C
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getParaFFTCommunicator(MyProc,NumProcs) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out), optional :: MyProc, NumProcs
   integer (kind=IntKind) :: c
!
   c = comm
!
   if (present(MyProc)) then
      MyProc = MyPE_comm
   endif
   if (present(NumProcs)) then
      NumProcs = NumPEs_comm
   endif
!
   end function getParaFFTCommunicator
!  ===================================================================
end module ParallelFFTModule
