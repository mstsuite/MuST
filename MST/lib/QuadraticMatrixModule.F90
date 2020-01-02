module QuadraticMatrixModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : CZERO, CONE, TEN2p30
!
public :: initQuadraticMatrix,    &
          solveQuadraticEquation, &  ! Solve s0+s1*delta+s2*delta^2 = 0 equation
          solveLinearEquation,    &  ! Solve    s1*delta+s2*delta^2 = 0 equation
          getEigenValue,          &
          getEigenVector,         &
          getEigenMatrix,         &
          endQuadraticMatrix
!
          interface getEigenValue
             module procedure geval1, geval2
          end interface
!
private
   logical :: isGeneral
   logical :: isLinear
!
   integer (kind=IntKind) :: ndim, ndim2
   integer (kind=IntKind), allocatable :: ipiv(:)
!
   complex (kind=CmplxKind), allocatable :: am(:,:), dm(:,:)
   complex (kind=CmplxKind), allocatable :: qm(:,:), rm(:,:), Qmatrix(:,:)
   complex (kind=CmplxKind), allocatable :: tw(:), a2m(:,:), S1inv(:,:), S2inv(:,:)
   complex (kind=CmplxKind), allocatable, target :: EigenVectorR(:), EigenValue(:)
   complex (kind=CmplxKind), allocatable, target :: EigenVectorL(:), EVLL(:)
   complex (kind=CmplxKind), allocatable, target :: EigenMatrix(:,:)
!
   integer (kind=IntKind) :: LWORK, LWMAX
   real (kind=RealKind), allocatable :: RWORK(:)
   real (kind=RealKind), allocatable :: RWORK0(:)
   complex (kind=CmplxKind), allocatable :: CWORK(:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initQuadraticMatrix(n,isGen)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: n
   logical, optional :: isGen
!
   if(present(isGen)) then
      isGeneral = isGen
   else
      isGeneral = .false.
   endif
!
   if (isGen) then
      call WarningHandler('initQuadraticMatrix','isGenral scheme is enabled!', &
                          'It needs to check the algebra for SubEigenMatrix')
   endif
!
   ndim = n
   ndim2 = ndim*2
!
   if (ndim < 1) then
      call ErrorHandler('initQuadraticMatrix','Matrix dimension < 1',n)
   endif
!
   allocate(am(ndim,ndim),dm(ndim,ndim))
   allocate(a2m(ndim,ndim),ipiv(ndim2),tw(n*16),S2inv(ndim,ndim))
   allocate(qm(ndim2,ndim2), rm(ndim2,ndim2), S1inv(ndim,ndim))
   allocate(EigenVectorR(ndim2*ndim2), EigenValue(ndim2))
   allocate(EigenVectorL(ndim2*ndim2), Qmatrix(ndim2,ndim2), EVLL(ndim))
   allocate(EigenMatrix(ndim,ndim))
!
   isLinear = .false.
!
   LWMAX = 20*ndim2
   allocate(CWORK(LWMAX))
   if(isGeneral) then
      allocate(RWORK(8*ndim2))
   else
      allocate(RWORK(2*ndim2))
   endif
   allocate(RWORK0(8*ndim))
!
   end subroutine initQuadraticMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endQuadraticMatrix()
!  ===================================================================
   implicit   none
!
   deallocate(am, a2m, dm, S1inv, S2inv, ipiv, tw, Qmatrix)
   deallocate(qm, rm, EigenVectorR, EigenValue, EigenVectorL, EVLL, EigenMatrix)
   deallocate(CWORK, RWORK, RWORK0)
!
   isLinear = .false.
!
   end subroutine endQuadraticMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine solveQuadraticEquation(s0,s1,s2,info)
!  ===================================================================
   implicit   none
!
   complex (kind=CmplxKind), intent(in) :: s0(ndim,*)
   complex (kind=CmplxKind), intent(in) :: s1(ndim,*)
   complex (kind=CmplxKind), intent(in) :: s2(ndim,*)
   complex (kind=CmplxKind) :: alpha(ndim2), beta(ndim2)
!
   integer (kind=IntKind), intent(out) :: info
   integer (kind=IntKind) :: i, j, k
!
   isLinear = .false.
!
   if(isGeneral) then
!     ================================================================
!     Set up the quadratic matrices.
!              |  -s1    -s0  |         |  s2     0  |
!         qm = |              |,  rm  = |            |
!              |   I      0   |         |   0     I  |
!     ================================================================
      do j = 1, ndim
         do i = 1, ndim
            qm(i,j)=-s1(i,j)
         enddo
         do i = 1, ndim
            qm(ndim+i,j) = CZERO
         enddo
         qm(ndim+j,j)=CONE
      enddo
      do j = 1, ndim
         do i = 1, ndim
            qm(i,ndim+j)=-s0(i,j)
         enddo
         do i = 1, ndim
            qm(ndim+i,ndim+j)=CZERO
         enddo
      enddo
!
      do j = 1, ndim
         do i = 1, ndim
            rm(i,j)=s2(i,j)
         enddo
         do i = 1, ndim
            rm(ndim+i,j) = CZERO
         enddo
      enddo
      do j = 1, ndim
         do i = 1, ndim
            rm(i,ndim+j)=CZERO
         enddo
         do i = 1, ndim
            rm(ndim+i,ndim+j)=CZERO
         enddo
         rm(ndim+j,ndim+j)=CONE
      enddo
!
!     ================================================================
!     Solve the generalized eigenvalue problem. 
!     ----------------------------------------------------------------
      LWORK = -1
      call zggev('v','v',ndim2,qm,ndim2,rm,ndim2,alpha,beta,          &
                 EigenVectorL,ndim2,EigenVectorR,ndim2,               &
                 CWORK,LWORK,RWORK,info)
!     ----------------------------------------------------------------
      LWORK = MIN( LWMAX, INT( CWORK( 1 ) ) )
      call zggev('v','v',ndim2,qm,ndim2,rm,ndim2,alpha,beta,          &
                 EigenVectorL,ndim2,EigenVectorR,ndim2,               &
                 CWORK,LWORK,RWORK,info)
!     ----------------------------------------------------------------
      if (info /= 0) then
         call ErrorHandler('solveQuadraticEquation',                  &
                          'ZGGEV failure, QM/RM is ill conditioned',info)
      endif
      do i = 1, ndim2
         write(6,'(a,4e20.12)')'alpha,beta=',alpha(i),beta(i)
         if(beta(i) == CZERO.and.alpha(i) == CZERO) then
            call ErrorHandler('solveQuadraticEquation','Undefined Eigenvalue')
         else if(beta(i) == CZERO) then
            EigenValue(i) = TEN2p30
         else
            EigenValue(i) = alpha(i)/beta(i)
         endif
      enddo
   else
!     ================================================================
!     calculate s1^(-1) and store the result in S1inv
!     ----------------------------------------------------------------
      call zcopy(ndim*ndim,s1,1,S1inv,1)
      call zgetrf( ndim, ndim, S1inv, ndim, ipiv, info )
!     ----------------------------------------------------------------
      if (info /= 0) then
         call ErrorHandler('solveQuadraticEquation',                  &
                           'ZGETRF failure, s1 is ill conditioned',info)
      endif
!     ----------------------------------------------------------------
      call zgetri( ndim, S1inv, ndim, ipiv, tw, ndim, info )
!     ----------------------------------------------------------------
      if (info /= 0) then
         call ErrorHandler('solveQuadraticEquation',                  &
                           'ZGETRI failure, s1 is ill conditioned',info)
      endif
!     ===================================================================
!     calculate s1^(-1)*s0 and store the result in am
!     ----------------------------------------------------------------
      call zgemm('n','n',ndim,ndim,ndim,CONE,S1inv,ndim,s0,ndim,CZERO,am,ndim)
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate s2^(-1) and store the result in S2inv
!     ----------------------------------------------------------------
      call zcopy(ndim*ndim,s2,1,S2inv,1)
      call zgetrf( ndim, ndim, S2inv, ndim, ipiv, info )
!     ----------------------------------------------------------------
      if (info /= 0) then
         call WarningHandler('solveQuadraticEquation',                &
                             'ZGETRF failure, s2 is ill conditioned',info)
         return
      endif
!     ----------------------------------------------------------------
      call zgetri( ndim, S2inv, ndim, ipiv, tw, ndim, info )
!     ----------------------------------------------------------------
      if (info /= 0) then
         call WarningHandler('solveQuadraticEquation',                &
                             'ZGETRI failure, s2 is ill conditioned',info)
         return
      endif
!     ================================================================
!     calculate s2^(-1)*s1 and store the result in dm
!     ----------------------------------------------------------------
      call zgemm('n','n',ndim,ndim,ndim,CONE,S2inv,ndim,s1,ndim,CZERO,dm,ndim)
!     ----------------------------------------------------------------
!
!     ================================================================
!     calculate am*am and store the result in a2m
!     ----------------------------------------------------------------
      call zgemm('n','n',ndim,ndim,ndim,CONE,am,ndim,am,ndim,CZERO,a2m,ndim)
!     ----------------------------------------------------------------
!
!     ================================================================
!     Set up the quadratic matrix.              
!              |  -am       I     |   |   D     I   |
!         qm = |                  | = |             |
!              |  -a2m   am - dm  |   | -D^2   B-D  |
!
!     where am = s1^{-1}*s0 = -D, a2m = am*am = D^2, dm = s2^{-1}*s1 = -B.
!     Note: D and B are the matrices defined in the notes.
!     ================================================================
      do j = 1, ndim
         do i = 1, ndim
            qm(i,j)=-am(i,j)
         enddo
         do i = 1, ndim
            qm(ndim+i,j) = -a2m(i,j)
         enddo
      enddo
      do j = 1, ndim
         do i = 1, ndim
            qm(i,ndim+j) = CZERO
         enddo
         qm(j,ndim+j)=CONE
         do i = 1, ndim
            qm(ndim+i,ndim+j)=am(i,j)-dm(i,j)
         enddo
      enddo
      call zcopy(ndim2*ndim2,qm,1,Qmatrix,1)
!
!     ================================================================
!     Solve the eigenvalue problem. The eigen vectors are not calculated
!     ----------------------------------------------------------------
      LWORK = -1
      call zgeev('v','v',ndim2,qm,ndim2,EigenValue,                   &
                 EigenVectorL,ndim2,EigenVectorR,ndim2,               &
                 CWORK,LWORK,RWORK,info)
!     ----------------------------------------------------------------
      LWORK = MIN( LWMAX, INT( CWORK( 1 ) ) )
      call zgeev('v','v',ndim2,qm,ndim2,EigenValue,                   &
                 EigenVectorL,ndim2,EigenVectorR,ndim2,               &
                 CWORK,LWORK,RWORK,info)
!     ----------------------------------------------------------------
      if (info /= 0) then
         call ErrorHandler('solveQuadraticEquation',                  &
                           'ZGEEV failure, QM is ill conditioned',info)
      endif
      EigenVectorL = conjg(EigenVectorL)   ! The left eigenvector given by zgeev
                                           ! is a complex conjugate with transpose
   endif
!
   end subroutine solveQuadraticEquation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine solveLinearEquation(s1,s2,info)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: s1(ndim,*)
   complex (kind=CmplxKind), intent(in) :: s2(ndim,*)
!
   integer (kind=IntKind), intent(out) :: info
   integer (kind=IntKind) :: i, j
!
   isLinear = .true.
!
!  ===================================================================
!  calculate s2^(-1) and store the result in S2inv
!  -------------------------------------------------------------------
   call zcopy(ndim*ndim,s2,1,S2inv,1)
   call zgetrf( ndim, ndim, S2inv, ndim, ipiv, info )
!  -------------------------------------------------------------------
   if (info /= 0) then
      call ErrorHandler('solveLinearEquation',                        &
                        'ZGETRF failure, s2 is ill conditioned',info)
   endif
!  -------------------------------------------------------------------
   call zgetri( ndim, S2inv, ndim, ipiv, tw, ndim, info )
!  -------------------------------------------------------------------
   if (info /= 0) then
      call ErrorHandler('solveLinearEquation',                        &
                        'ZGETRI failure, s2 is ill conditioned',info)
   endif
!  ===================================================================
!  calculate s2^(-1)*s1 and store the result in am
!  -------------------------------------------------------------------
   call zgemm('n','n',ndim,ndim,ndim,CONE,S2inv,ndim,s1,ndim,CZERO,am,ndim)
!  -------------------------------------------------------------------
   do i = 1, ndim
      do j = 1, ndim
         am(i,j) = -am(i,j)
      enddo
   enddo
!
!  ===================================================================
!  Solve the eigenvalue problem. The eigen vectors are not calculated
!            |             |   |   |
!       am = | -s2^{-1}*s1 | = | B |,  where B is the matrix defined in the notes
!            |             |   |   |
!  -------------------------------------------------------------------
   LWORK = -1
   call zgeev('v','v',ndim,am,ndim,EigenValue,                        &
              EigenVectorL,ndim,EigenVectorR,ndim,                    &
              CWORK,LWORK,RWORK0,info)
!  -------------------------------------------------------------------
   LWORK = MIN( LWMAX, INT( CWORK( 1 ) ) )
   call zgeev('v','v',ndim,am,ndim,EigenValue,                        &
              EigenVectorL,ndim,EigenVectorR,ndim,                    &
              CWORK,LWORK,RWORK0,info)
!  -------------------------------------------------------------------
   if (info /= 0) then
      call ErrorHandler('solveLinearEquation',                        &
                        'ZGEEV failure, QM is ill conditioned',info)
   endif
   EigenVectorL = conjg(EigenVectorL)
!
   end subroutine solveLinearEquation
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function geval1(nv) result(v)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(out) :: nv
!
   complex (kind=CmplxKind), pointer :: v(:)
!
   if (isLinear) then
      nv = ndim
   else
      nv = ndim2
   endif
!
   v => EigenValue(1:nv)
!
   end function geval1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function geval2(i,nv) result(v)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind), intent(out) :: nv
!
   complex (kind=CmplxKind) :: v
!
   if (isLinear) then
      nv = ndim
   else
      nv = ndim2
   endif
!
   if (i < 1 .or. i > nv) then
      call ErrorHandler('getEigenValue','Invalid eigenvalue index',i)
   endif
!
   v = EigenValue(i)
!
   end function geval2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEigenVector(a,i) result(v)
!  ===================================================================
!  Note: For eigenvector (u1,u2) of the quadratic matrix, u1 is an eigenvector 
!        of the following matrix 
!             A = lamda^2 + lamda*(S2^{inv}*S1) + S2^{inv}*S0
!        If u1 is the right-hand side eigenvector of A, it is also the
!        right-hand side eigenvector of S2*A
!        If u2 is the left-hand side eigenvector of A, u2*S2inv is then
!        the left-hand side eigenvector of S2*A, and it is therefore
!        necessary for result v to include S2inv if a = 'L'.
!  ******************************************************************
   implicit   none
!
   character (len=1), intent(in), optional :: a
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: mdim, j, k
!
   complex (kind=CmplxKind), pointer :: v(:)
   complex (kind=CmplxKind), pointer :: vp(:)
!
   if (isLinear) then
      mdim = ndim
   else
      mdim = ndim2
   endif
!
   if (i < 1 .or. i > mdim) then
      call ErrorHandler('getEigenVector','Invalid eigenvector index',i)
   endif
!
   if (present(a)) then
      if (a == 'R' .or. a == 'r') then
         v => EigenVectorR((i-1)*mdim+1:(i-1)*mdim+ndim)
      else if (a == 'L' .or. a == 'l') then
         if (isLinear) then
            vp => EigenVectorL((i-1)*mdim+1:i*mdim)
         else
            vp => EigenVectorL((i-1)*mdim+ndim+1:i*mdim)
         endif
         do k = 1, ndim
            EVLL(k) = CZERO
            do j = 1, ndim
               EVLL(k) = EVLL(k) + vp(j)*S2inv(j,k)
            enddo
         enddo
         v => EVLL
      else
         call ErrorHandler('getEigenVector','Invalid eigenvector type',a)
      endif
   else 
      v => EigenVectorR((i-1)*mdim+1:(i-1)*mdim+ndim)
   endif
!
   end function getEigenVector
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEigenMatrix(i) result(a)
!  ===================================================================
!
!  Returns: a = vr((n-1)*ndim+1:n*ndim,i) * vl(i,(m-1)*ndim+1:m*ndim)
!
!  Here: vr = EigenVectorR is the right-side eigenvector
!        vl = EigenVectorL is the left-side eigenvector
!  *******************************************************************
   implicit   none
!
   integer (kind=IntKind), intent(in) :: i
!
   integer (kind=IntKind) :: j, k, k0, m0, n0, n02
!
   complex (kind=CmplxKind), pointer :: a(:,:)
   complex (kind=CmplxKind) :: c
!
   if (i < 1 .or. i > ndim2 .or. (isLinear .and. i > ndim)) then
      call ErrorHandler('getSubEigenMatrix','Invalid eigenvector index',i)
   endif
!
   if (isLinear) then
      m0 = (i-1)*ndim
      k0 = m0
   else
      m0 = (i-1)*ndim2
      k0 = m0 + ndim
   endif
!
   do k = 1, ndim
      c = CZERO
      do j = 1, ndim
         c = c + EigenVectorL(k0+j)*S2inv(j,k)  ! Multiplied by S2^{-1}
      enddo
      do j = 1, ndim
         EigenMatrix(j,k) = EigenVectorR(m0+j)*c
      enddo
   enddo
!
   a => EigenMatrix
!
   end function getEigenMatrix
!  ===================================================================
end module QuadraticMatrixModule
