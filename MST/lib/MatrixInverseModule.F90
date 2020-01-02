module MatrixInverseModule
    use KindParamModule, only : IntKind, RealKind, CmplxKind
!
    use MathParamModule, only : one, zero, cone, czero
!
    use ErrorHandlerModule, only : ErrorHandler
!
public :: MtxInv_GE, MtxInv_LU, ABplus1Inv
!
    interface MtxInv_GE
       module procedure MtxInv_GEr0, MtxInv_GEr1, MtxInv_GEr2, MtxInv_GEr3
       module procedure MtxInv_GEc0, MtxInv_GEc1, MtxInv_GEc2, MtxInv_GEc3
    end interface
!
    interface MtxInv_LU
       module procedure MtxInv_LUr, MtxInv_LUc, MtxInv_LUc1, MtxInv_LUc2
    end interface
!
    interface ABplus1Inv
       module procedure ABplus1Inv_1, ABplus1Inv_2
    end interface
!
private
!
contains
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEc0(n,am,bmt,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n
    integer (kind=IntKind) :: i,j
!
    complex (kind=CmplxKind), intent(in) :: am(:,:)
    complex (kind=CmplxKind), intent(out) :: bmt(n,n)
    complex (kind=CmplxKind), intent(out), optional :: det
!
    complex (kind=CmplxKind), allocatable :: amt(:,:)
    complex (kind=CmplxKind), allocatable :: ad(:)
    complex (kind=CmplxKind), allocatable :: bd(:)
    complex (kind=CmplxKind), allocatable :: td(:)
    complex (kind=CmplxKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      am,       the matrix to be inversed.
!               n,        the dimension of amt.
!
!   output:     bmt = am**(-1)
!               det = Det|amt|
!
!   BLAS level 2.
!   ******************************************************************
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=czero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=cone
    enddo
!
!   ------------------------------------------------------------------
    allocate(ad(n), bd(n), td(n), amt(n,n))
!   ------------------------------------------------------------------
    do i=1,n
       amt(1:n,i)=am(1:n,i)
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=cone/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=czero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call zgeru(n,i,-cone,td,1,bd,1,bmt,n)
       call zgeru(n,n-i+1,-cone,td,1,ad(i),1,amt(1,i),n)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=cone
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, amt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEc0
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEc1(n,am,lda,bmt,ldb,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n, lda, ldb
    integer (kind=IntKind) :: i,j
!
    complex (kind=CmplxKind), intent(in) :: am(lda,n)
    complex (kind=CmplxKind), intent(out) :: bmt(ldb,n)
    complex (kind=CmplxKind), intent(out), optional :: det
!
    complex (kind=CmplxKind), allocatable :: amt(:,:)
    complex (kind=CmplxKind), allocatable :: ad(:)
    complex (kind=CmplxKind), allocatable :: bd(:)
    complex (kind=CmplxKind), allocatable :: td(:)
    complex (kind=CmplxKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      am,       the matrix to be inversed.
!               n,        the dimension of amt.
!
!   output:     bmt = am**(-1)
!
!
!   BLAS level 2.
!   ******************************************************************
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=czero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=cone
    enddo
!
!   ------------------------------------------------------------------
    allocate(ad(n), bd(n), td(n), amt(n,n))
!   ------------------------------------------------------------------
    do i=1,n
       amt(1:n,i)=am(1:n,i)
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=cone/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=czero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call zgeru(n,i,-cone,td,1,bd,1,bmt,ldb)
       call zgeru(n,n-i+1,-cone,td,1,ad(i),1,amt(1,i),n)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=cone
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, amt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEc1
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEc2(n,amt,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n
    integer (kind=IntKind) :: i,j
!
    complex (kind=CmplxKind), intent(inout) :: amt(:,:)
    complex (kind=CmplxKind), intent(out), optional :: det
!
    complex (kind=CmplxKind), allocatable :: bmt(:,:)
    complex (kind=CmplxKind), allocatable :: ad(:)
    complex (kind=CmplxKind), allocatable :: bd(:)
    complex (kind=CmplxKind), allocatable :: td(:)
    complex (kind=CmplxKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      amt,       the matrix to be inversed.
!               n,         the dimension of amt.
!
!   output:     amt = amt**(-1)
!               det = Det|amt|
!
!   WARNING:: amt will be destroyed after calling this routine.
!
!   BLAS level 2.
!   ******************************************************************
!
!   ------------------------------------------------------------------
    allocate(bmt(n,n), ad(n), bd(n), td(n))
!   ------------------------------------------------------------------
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=czero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=cone
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=cone/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=czero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call zgeru(n,i,-cone,td,1,bd,1,bmt,n)
       call zgeru(n,n-i+1,-cone,td,1,ad(i),1,amt(1,i),n)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=cone
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!
    amt=bmt
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, bmt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEc2
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEc3(n,amt,lda,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n, lda
    integer (kind=IntKind) :: i,j
!
    complex (kind=CmplxKind), intent(inout) :: amt(lda,n)
    complex (kind=CmplxKind), intent(out), optional :: det
!
    complex (kind=CmplxKind), allocatable :: bmt(:,:)
    complex (kind=CmplxKind), allocatable :: ad(:)
    complex (kind=CmplxKind), allocatable :: bd(:)
    complex (kind=CmplxKind), allocatable :: td(:)
    complex (kind=CmplxKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      amt,       the matrix to be inversed.
!               n,         the dimension of amt.
!
!   output:     amt = amt**(-1)
!
!   WARNING:: amt will be destroyed after calling this routine.
!
!   BLAS level 2.
!   ******************************************************************
!
!   ------------------------------------------------------------------
    allocate(bmt(n,n), ad(n), bd(n), td(n))
!   ------------------------------------------------------------------
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=czero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=cone
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=cone/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=czero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call zgeru(n,i,-cone,td,1,bd,1,bmt,n)
       call zgeru(n,n-i+1,-cone,td,1,ad(i),1,amt(1,i),lda)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=cone
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!
    amt=bmt
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, bmt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEc3
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEr0(n,am,bmt,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n
    integer (kind=IntKind) :: i,j
!
    real (kind=RealKind), intent(in) :: am(:,:)
    real (kind=RealKind), intent(out) :: bmt(n,n)
    real (kind=RealKind), intent(out), optional :: det
!
    real (kind=RealKind), allocatable :: amt(:,:)
    real (kind=RealKind), allocatable :: ad(:)
    real (kind=RealKind), allocatable :: bd(:)
    real (kind=RealKind), allocatable :: td(:)
    real (kind=RealKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      am,       the matrix to be inversed.
!               n,        the dimension of amt.
!
!   output:     bmt = am**(-1)
!               det = Det|amt|
!
!   BLAS level 2.
!   ******************************************************************
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=zero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=one
    enddo
!
!   ------------------------------------------------------------------
    allocate(ad(n), bd(n), td(n), amt(n,n))
!   ------------------------------------------------------------------
    do i=1,n
       amt(1:n,i)=am(1:n,i)
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=one/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=zero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call dger(n,i,-one,td,1,bd,1,bmt,n)
       call dger(n,n-i+1,-one,td,1,ad(i),1,amt(1,i),n)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=one
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, amt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEr0
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEr1(n,am,lda,bmt,ldb,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n, lda, ldb
    integer (kind=IntKind) :: i,j
!
    real (kind=RealKind), intent(in) :: am(lda,n)
    real (kind=RealKind), intent(out) :: bmt(ldb,n)
    real (kind=RealKind), intent(out), optional :: det
!
    real (kind=RealKind), allocatable :: amt(:,:)
    real (kind=RealKind), allocatable :: ad(:)
    real (kind=RealKind), allocatable :: bd(:)
    real (kind=RealKind), allocatable :: td(:)
    real (kind=RealKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      am,       the matrix to be inversed.
!               n,        the dimension of amt.
!
!   output:     bmt = am**(-1)
!
!
!   BLAS level 2.
!   ******************************************************************
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=zero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=one
    enddo
!
!   ------------------------------------------------------------------
    allocate(ad(n), bd(n), td(n), amt(n,n))
!   ------------------------------------------------------------------
    do i=1,n
       amt(1:n,i)=am(1:n,i)
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=one/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=zero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call dger(n,i,-one,td,1,bd,1,bmt,ldb)
       call dger(n,n-i+1,-one,td,1,ad(i),1,amt(1,i),n)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=one
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, amt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEr1
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEr2(n,amt,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n
    integer (kind=IntKind) :: i,j
!
    real (kind=RealKind), intent(inout) :: amt(:,:)
    real (kind=RealKind), intent(out), optional :: det
!
    real (kind=RealKind), allocatable :: bmt(:,:)
    real (kind=RealKind), allocatable :: ad(:)
    real (kind=RealKind), allocatable :: bd(:)
    real (kind=RealKind), allocatable :: td(:)
    real (kind=RealKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      amt,       the matrix to be inversed.
!               n,         the dimension of amt.
!
!   output:     amt = amt**(-1)
!               det = Det|amt|
!
!   WARNING:: amt will be destroyed after calling this routine.
!
!   BLAS level 2.
!   ******************************************************************
!
!   ------------------------------------------------------------------
    allocate(bmt(n,n), ad(n), bd(n), td(n))
!   ------------------------------------------------------------------
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=zero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=one
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=one/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=zero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call dger(n,i,-one,td,1,bd,1,bmt,n)
       call dger(n,n-i+1,-one,td,1,ad(i),1,amt(1,i),n)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=one
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!
    amt=bmt
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, bmt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEr2
!   ==================================================================
!
!   ******************************************************************
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine MtxInv_GEr3(n,amt,lda,det)
!   ==================================================================
!
    implicit   none
!
    integer (kind=IntKind), intent(in) :: n, lda
    integer (kind=IntKind) :: i,j
!
    real (kind=RealKind), intent(inout) :: amt(lda,n)
    real (kind=RealKind), intent(out), optional :: det
!
    real (kind=RealKind), allocatable :: bmt(:,:)
    real (kind=RealKind), allocatable :: ad(:)
    real (kind=RealKind), allocatable :: bd(:)
    real (kind=RealKind), allocatable :: td(:)
    real (kind=RealKind) :: amtinv
!
!   ******************************************************************
!   vectorized version for matrix inversion.........................
!
!   input:      amt,       the matrix to be inversed.
!               n,         the dimension of amt.
!
!   output:     amt = amt**(-1)
!
!   WARNING:: amt will be destroyed after calling this routine.
!
!   BLAS level 2.
!   ******************************************************************
!
!   ------------------------------------------------------------------
    allocate(bmt(n,n), ad(n), bd(n), td(n))
!   ------------------------------------------------------------------
!
!   ==================================================================
!   set up unit matrix................................................
!   ------------------------------------------------------------------
    bmt=zero
!   ------------------------------------------------------------------
    do i=1,n
       bmt(i,i)=one
    enddo
!
!   ==================================================================
!   transform amt and bmt.............................................
!   ==================================================================
    do i=1,n
       amtinv=one/amt(i,i)
       do j=1,n
          td(j)=amtinv*amt(j,i)
       enddo
       do j=1,n
          bd(j)=bmt(i,j)
       enddo
       ad(i)=amt(i,i)
       td(i)=zero
       do j=i+1,n
          ad(j)=amt(i,j)
          td(j)=amtinv*amt(j,i)
       enddo
!      ---------------------------------------------------------------
       call dger(n,i,-one,td,1,bd,1,bmt,n)
       call dger(n,n-i+1,-one,td,1,ad(i),1,amt(1,i),lda)
!      ---------------------------------------------------------------
    enddo
!
!   ==================================================================
!   fill bmt by the inverse of amt....................................
!   ==================================================================
    do j=1,n
       do i=1,n
          bmt(i,j)=bmt(i,j)/amt(i,i)
       end do
    end do
!
!   ==================================================================
!   calculate determinant ..........................................
!   ==================================================================
    if (present(det)) then
       det=one
       do i=1,n
          det=det*amt(i,i)
       enddo
    endif
!
    amt=bmt
!   ------------------------------------------------------------------
    deallocate(ad, bd, td, bmt)
!   ------------------------------------------------------------------
!
    end subroutine MtxInv_GEr3
!   ==================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MtxInv_LUr( A, m )
!  ===================================================================
   implicit none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
   integer (kind=IntKind), intent(in) :: m
!
   real (kind=RealKind) :: A(m,m)
!
!  ===================================================================
!  Local variables
!  ===================================================================
!
   integer (kind=IntKind)   :: info
   integer (kind=IntKind)   :: ipiv(m)
   integer (kind=intKind), parameter :: nb = 16 ! this number may notbe optimal
!
   real (kind=RealKind) :: work(m*nb)
!
!  ===================================================================
!   Find LU factorization of matrix
!  -------------------------------------------------------------------
   call dgetrf( m, m, A, m, ipiv, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the LU factorization worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with dgetrf: info',info)
   endif
!
!  ===================================================================
!  Find inverse of matrix
!  -------------------------------------------------------------------
   call dgetri( m, A, m, ipiv, work, m*nb, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the Inverse worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with dgetri: info',info)
   endif
!
   end subroutine MtxInv_LUr
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MtxInv_LUc( A, m )
!  ===================================================================
   implicit none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
   integer (kind=IntKind), intent(in) :: m
!
   complex (kind=CmplxKind) :: A(m,m)
!
!  ===================================================================
!  Local variables
!  ===================================================================
!
   integer (kind=IntKind)   :: info
   integer (kind=IntKind)   :: ipiv(m)
   integer (kind=intKind), parameter :: nb = 16 ! this number may notbe optimal
!
   complex (kind=CmplxKind) :: work(m*nb)
!
!  ===================================================================
!   Find LU factorization of matrix
!  -------------------------------------------------------------------
   call zgetrf( m, m, A, m, ipiv, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the LU factorization worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with zgetrf: info',info)
   endif
!
!  ===================================================================
!  Find inverse of matrix
!  -------------------------------------------------------------------
   call zgetri( m, A, m, ipiv, work, m*nb, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the Inverse worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with zgetri: info',info)
   endif
!
   end subroutine MtxInv_LUc
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MtxInv_LUc1( m, A, lda )
!  ===================================================================
   implicit none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
   integer (kind=IntKind), intent(in) :: m, lda
!
   complex (kind=CmplxKind) :: A(m*lda)
!
!  ===================================================================
!  Local variables
!  ===================================================================
!
   integer (kind=IntKind)   :: info
   integer (kind=IntKind)   :: ipiv(m)
   integer (kind=intKind), parameter :: nb = 16 ! this number may notbe optimal
!
   complex (kind=CmplxKind) :: work(m*nb)
!
!  ===================================================================
!   Find LU factorization of matrix
!  -------------------------------------------------------------------
   call zgetrf( m, m, A, lda, ipiv, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the LU factorization worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with zgetrf: info',info)
   endif
!
!  ===================================================================
!  Find inverse of matrix
!  -------------------------------------------------------------------
   call zgetri( m, A, lda, ipiv, work, m*nb, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the Inverse worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with zgetri: info',info)
   endif
!
   end subroutine MtxInv_LUc1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MtxInv_LUc2( m, A, lda )
!  ===================================================================
   implicit none
!
!  ===================================================================
!  I/O variables
!  ===================================================================
   integer (kind=IntKind), intent(in) :: m, lda
!
   complex (kind=CmplxKind) :: A(lda,m)
!
!  ===================================================================
!  Local variables
!  ===================================================================
!
   integer (kind=IntKind)   :: info
   integer (kind=IntKind)   :: ipiv(m)
   integer (kind=intKind), parameter :: nb = 16 ! this number may notbe optimal
!
   complex (kind=CmplxKind) :: work(m*nb)
!
!  ===================================================================
!   Find LU factorization of matrix
!  -------------------------------------------------------------------
   call zgetrf( m, m, A, lda, ipiv, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the LU factorization worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with zgetrf: info',info)
   endif
!
!  ===================================================================
!  Find inverse of matrix
!  -------------------------------------------------------------------
   call zgetri( m, A, lda, ipiv, work, m*nb, info )
!  -------------------------------------------------------------------
!  ===================================================================
!  Check to make sure that the Inverse worked correctly
!  ===================================================================
   if ( info /= 0 ) then
      call ErrorHandler('MtxInv_LU','problem with zgetri: info',info)
   endif
!
   end subroutine MtxInv_LUc2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ABplus1Inv_1(n,a,lda,b,ldb,c,ldc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, lda, ldb, ldc
   integer (kind=IntKind) :: i, j, k
!
   complex (kind=CmplxKind), intent(in) :: a(lda*n)
   complex (kind=CmplxKind), intent(in) :: b(ldb*n)
   complex (kind=CmplxKind), intent(out) :: c(ldc*n)
!
   if (n < 1) then
      call ErrorHandler('ABplus1Inv','Invalid n',n)
   else if (n > lda) then
      call ErrorHandler('ABplus1Inv','n > lda',n,lda)
   else if (n > ldb) then
      call ErrorHandler('ABplus1Inv','n > ldb',n,ldb)
   else if (n > ldc) then
      call ErrorHandler('ABplus1Inv','n > ldc',n,ldc)
   endif
!
   c = CZERO
!
   do i = 1, n
      c((i-1)*ldc+i) = CONE
   enddo
!
!  -------------------------------------------------------------------
   call zgemm('n', 'n', n, n, n, CONE, a, lda, b, ldb, CONE, c, ldc)
   call MtxInv_LU(n, c, ldc)
!  -------------------------------------------------------------------
!
   end subroutine ABplus1Inv_1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ABplus1Inv_2(n,a,lda,b,ldb,c,ldc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, lda, ldb, ldc
   integer (kind=IntKind) :: i, j, k
!
   complex (kind=CmplxKind), intent(in) :: a(lda,n)
   complex (kind=CmplxKind), intent(in) :: b(ldb,n)
   complex (kind=CmplxKind), intent(out) :: c(ldc,n)
!
   if (n < 1) then
      call ErrorHandler('ABplus1Inv','Invalid n',n)
   else if (n > lda) then
      call ErrorHandler('ABplus1Inv','n > lda',n,lda)
   else if (n > ldb) then
      call ErrorHandler('ABplus1Inv','n > ldb',n,ldb)
   else if (n > ldc) then
      call ErrorHandler('ABplus1Inv','n > ldc',n,ldc)
   endif
!
   c = CZERO
!
   do i = 1, n
      c(i,i) = CONE
   enddo
!
!  -------------------------------------------------------------------
   call zgemm('n', 'n', n, n, n, CONE, a, lda, b, ldb, CONE, c, ldc)
   call MtxInv_LU(n, c, ldc)
!  -------------------------------------------------------------------
!
   end subroutine ABplus1Inv_2
!  ===================================================================
end module MatrixInverseModule
