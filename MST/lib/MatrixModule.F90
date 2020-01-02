module MatrixModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, ONE, CZERO, CONE
!
   use ErrorHandlerModule, only : ErrorHandler
!
public :: setupUnitMatrix, &
          computeAStar,    & ! compute B = A^{dot}
          computeAStarInv,    & ! compute B = A^{-dot}
          computeAStarT,   & ! compute B = A^{Tdot}
          computeAStarTInv,   & ! compute B = A^{-Tdot}
          computeUAU,      & ! compute B = U1 * A * U2
          computeUAUt,     & ! compute B = U1 * A * U2^{T}
          computeUAUts,    & ! compute B = U1 * A * U2^{Tdot}
          computeUAUtc,    & ! compute B = U1 * A * conjg[U2^{T}]
          computeAprojB      ! compute C = [1 + A * B]^{-1} * A or
                             !         C = A * [1 + B * A]^{-1}
!
   interface setupUnitMatrix
      module procedure setupUMr1, setupUMc1, setupUMr12, setupUMc12,  &
                       setupUMr2, setupUMc2, setupUMr3, setupUMc3
   end interface
!
private
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMr1(n,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(out) :: a(n*n)
   real (kind=RealKind), intent(in), optional :: d
   real (kind=RealKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = ONE
   endif
!
   a = ZERO
   do i = 1, n
      a(i+(i-1)*n) = fac
   enddo
!
   end subroutine setupUMr1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMc1(n,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), intent(out) :: a(n*n)
   complex (kind=CmplxKind), intent(in), optional :: d
   complex (kind=CmplxKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = CONE
   endif
!
   a = CZERO
   do i = 1, n
      a(i+(i-1)*n) = fac
   enddo
!
   end subroutine setupUMc1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMr12(n,m,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, m
   integer (kind=IntKind) :: i, k
!
   real (kind=RealKind), intent(out) :: a(n*m)
   real (kind=RealKind), intent(in), optional :: d
   real (kind=RealKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = ONE
   endif
!
   a = ZERO
   k = min(n,m)
   do i = 1, k
      a(i+(i-1)*n) = fac
   enddo
!
   end subroutine setupUMr12
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMc12(n,m,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, m
   integer (kind=IntKind) :: i, k
!
   complex (kind=CmplxKind), intent(out) :: a(n*m)
   complex (kind=CmplxKind), intent(in), optional :: d
   complex (kind=CmplxKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = CONE
   endif
!
   a = CZERO
   k = min(n,m)
   do i = 1, k
      a(i+(i-1)*n) = fac
   enddo
!
   end subroutine setupUMc12
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMr2(n,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(out) :: a(n,n)
   real (kind=RealKind), intent(in), optional :: d
   real (kind=RealKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = ONE
   endif
!
   a = ZERO
   do i = 1, n
      a(i,i) = fac
   enddo
!
   end subroutine setupUMr2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMc2(n,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), intent(out) :: a(n,n)
   complex (kind=CmplxKind), intent(in), optional :: d
   complex (kind=CmplxKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = CONE
   endif
!
   a = CZERO
   do i = 1, n
      a(i,i) = fac
   enddo
!
   end subroutine setupUMc2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMr3(n,m,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, m
   integer (kind=IntKind) :: i, k
!
   real (kind=RealKind), intent(out) :: a(n,m)
   real (kind=RealKind), intent(in), optional :: d
   real (kind=RealKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = ONE
   endif
!
   a = ZERO
   k = min(n,m)
   do i = 1, k
      a(i,i) = fac
   enddo
!
   end subroutine setupUMr3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupUMc3(n,m,a,d)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, m
   integer (kind=IntKind) :: i, k
!
   complex (kind=CmplxKind), intent(out) :: a(n,m)
   complex (kind=CmplxKind), intent(in), optional :: d
   complex (kind=CmplxKind) :: fac
!
   if (present(d)) then
      fac = d
   else
      fac = CONE
   endif
!
   a = CZERO
   k = min(n,m)
   do i = 1, k
      a(i,i) = fac
   enddo
!
   end subroutine setupUMc3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeUAUts(Ul,nr,nc,Ur,mr,alpha,A,lda,beta,B,ldb,W)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc, mr, lda, ldb
   integer (kind=IntKind) :: n
!
   complex (kind=CmplxKind), intent(in) :: Ul(nr,nc), Ur(mr,nc), A(lda,nc)
   complex (kind=CmplxKind), intent(in) :: alpha, beta
   complex (kind=CmplxKind), intent(inout) :: B(ldb,mr)
   complex (kind=CmplxKind), intent(out), target :: W(:)
   complex (kind=CmplxKind), pointer :: Urc(:,:), UA(:,:)
!
   n = mr*nc
   if (n+nr*nc > size(W)) then
      call ErrorHandler('computeUAUts','Insufficient work space size',&
                        n+nr*nc,size(W))
   endif
   Urc => aliasArray2_c(W,mr,nc)
   UA => aliasArray2_c(W(n+1:n+nr*nc),nr,nc)
!  -------------------------------------------------------------------
   call computeAStar(Ur,mr,nc,Urc)
!  -------------------------------------------------------------------
   call zgemm('n','n',nr,nc,nc,alpha,Ul,nr,A,lda,CZERO,UA,nr)
!  -------------------------------------------------------------------
   call zgemm('n','t',nr,mr,nc,CONE,UA,nr,Urc,mr,beta,B,ldb)
!  -------------------------------------------------------------------

   end subroutine computeUAUts
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeUAUtc(Ul,nr,nc,Ur,mr,alpha,A,lda,beta,B,ldb,W)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc, mr, lda, ldb
   integer (kind=IntKind) :: i, j
!! integer (kind=IntKind) :: n
!
   complex (kind=CmplxKind), intent(in) :: Ul(nr,nc), Ur(mr,nc), A(lda,nc)
   complex (kind=CmplxKind), intent(in) :: alpha, beta
   complex (kind=CmplxKind), intent(inout) :: B(ldb,mr)
   complex (kind=CmplxKind), intent(out), target :: W(:)
!! complex (kind=CmplxKind), pointer :: Urc(:,:)
   complex (kind=CmplxKind), pointer :: UA(:,:)
!
!! Urc => aliasArray2_c(W,mr,nc)
!! n = mr*nc
!! if (n+nr*nc > size(W)) then
   if (nr*nc > size(W)) then
      call ErrorHandler('computeUAUtc','Insufficient work space size',&
                        nr*nc,size(W))
   endif
   UA => aliasArray2_c(W(1:nr*nc),nr,nc)
!  -------------------------------------------------------------------
!! Urc = conjg(Ur)
!  do j = 1, nc
!     do i = 1, mr
!        Urc(i,j) = conjg(Ur(i,j))
!     enddo
!  enddo
!  -------------------------------------------------------------------
   call zgemm('n','n',nr,nc,nc,alpha,Ul,nr,A,lda,CZERO,UA,nr)
!  -------------------------------------------------------------------
!! call zgemm('n','t',nr,mr,nc,CONE,UA,nr,Urc,mr,beta,B,nr)
   call zgemm('n','c',nr,mr,nc,CONE,UA,nr,Ur,mr,beta,B,ldb)
!  -------------------------------------------------------------------

   end subroutine computeUAUtc
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeUAU(Ul,nr,nc,Ur,mr,alpha,A,lda,beta,B,ldb,W)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc, mr, lda, ldb
   integer (kind=IntKind) :: i, j
!! integer (kind=IntKind) :: n
!
   complex (kind=CmplxKind), intent(in) :: Ul(nr,nc), Ur(mr,nc), A(lda,nc)
   complex (kind=CmplxKind), intent(in) :: alpha, beta
   complex (kind=CmplxKind), intent(inout) :: B(ldb,mr)
   complex (kind=CmplxKind), intent(out), target :: W(:)
   complex (kind=CmplxKind), pointer :: UA(:,:)
!
!! n = mr*nc
!! if (n+nr*nc > size(W)) then
   if (nr*nc > size(W)) then
      call ErrorHandler('computeUAU','Insufficient work space size',  &
                        nr*nc,size(W))
   endif
   UA => aliasArray2_c(W(1:nr*nc),nr,nc)
!  -------------------------------------------------------------------
   call zgemm('n','n',nr,nc,nc,alpha,Ul,nr,A,lda,CZERO,UA,nr)
!  -------------------------------------------------------------------
   call zgemm('n','n',nr,mr,nc,CONE,UA,nr,Ur,mr,beta,B,ldb)
!  -------------------------------------------------------------------

   end subroutine computeUAU
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeUAUt(Ul,nr,nc,Ur,mr,alpha,A,lda,beta,B,ldb,W)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc, mr, lda, ldb
   integer (kind=IntKind) :: i, j
!! integer (kind=IntKind) :: n
!
   complex (kind=CmplxKind), intent(in) :: Ul(nr,nc), Ur(mr,nc), A(lda,nc)
   complex (kind=CmplxKind), intent(in) :: alpha, beta
   complex (kind=CmplxKind), intent(inout) :: B(ldb,mr)
   complex (kind=CmplxKind), intent(out), target :: W(:)
   complex (kind=CmplxKind), pointer :: UA(:,:)
!
!! n = mr*nc
!! if (n+nr*nc > size(W)) then
   if (nr*nc > size(W)) then
      call ErrorHandler('computeUAUt','Insufficient work space size', &
                        nr*nc,size(W))
   endif
   UA => aliasArray2_c(W(1:nr*nc),nr,nc)
!  -------------------------------------------------------------------
   call zgemm('n','n',nr,nc,nc,alpha,Ul,nr,A,lda,CZERO,UA,nr)
!  -------------------------------------------------------------------
   call zgemm('n','t',nr,mr,nc,CONE,UA,nr,Ur,mr,beta,B,ldb)
!  -------------------------------------------------------------------

   end subroutine computeUAUt
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeAStar(A,nr,nc,B)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc
   integer (kind=IntKind) :: kl, l, m, klp, lp, mp
   integer (kind=IntKind) :: lmr, lmc, n, ma
!
   complex (kind=CmplxKind), intent(in) :: A(nr,nc)
   complex (kind=CmplxKind), intent(out) :: B(nr,nc)
   complex (kind=CmplxKind) :: cfac
!
   n = 0
   LOOP_lmc: do lmc = 1, nc
      if (lmc*lmc == nc) then
         n = lmc
         exit LOOP_lmc
      endif
   enddo LOOP_lmc
   if (n == 0) then
      call ErrorHandler('computeAStar','sqrt[nc] is not a integer',nc)
   endif
   lmc = n - 1
!
   n = 0
   LOOP_lmr: do lmr = 1, nr
      if (lmr*lmr == nr) then
         n = lmr
         exit LOOP_lmr
      endif
   enddo LOOP_lmr
   if (n == 0) then
      call ErrorHandler('computeAStar','sqrt[nr] is not a integer',nr)
   endif
   lmr = n - 1
!
   kl = 0
   do l = 0, lmc
      do m = -l, l
         kl = kl + 1
         klp = 0
         do lp = 0, lmr
            do mp = -lp, lp
               klp = klp + 1
               ma = abs(m+mp); n = 1-2*mod(ma,2)
               B(klp,kl) = n*A(klp-2*mp,kl-2*m)
            enddo
         enddo
      enddo
   enddo
!
   end subroutine computeAStar
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeAStarInv(A,nr,nc,B)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc
   integer (kind=IntKind) :: ipiv(2*nc), info
!
   complex (kind=CmplxKind), intent(in) :: A(nr,nc)
   complex (kind=CmplxKind), intent(out) :: B(nr,nc)
   complex (kind=CmplxKind) :: w(16*nc)
!
!  -------------------------------------------------------------------
   call computeAStar(a,nr,nc,B)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call zgetrf( nc, nc, B, nr, ipiv, info )
!  -------------------------------------------------------------------
   if (info /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('computeAStarInv','ZGETRF failure, B is ill conditioned',info)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call zgetri( nc, B, nr, ipiv, w, 16*nc, info)
!  -------------------------------------------------------------------
   if (info /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('computeAStarInv','ZGETRI failure, B is ill conditioned',info)
!     ----------------------------------------------------------------
   endif
!
   end subroutine computeAStarInv
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeAStarT(A,nr,nc,B)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc
   integer (kind=IntKind) :: kl, l, m, klp, lp, mp
   integer (kind=IntKind) :: lmr, lmc, n, ma
!
   complex (kind=CmplxKind), intent(in) :: A(nr,nc)
   complex (kind=CmplxKind), intent(out) :: B(nr,nc)
   complex (kind=CmplxKind) :: cfac
!
   n = 0
   LOOP_lmc: do lmc = 1, nc
      if (lmc*lmc == nc) then
         n = lmc
         exit LOOP_lmc
      endif
   enddo LOOP_lmc
   if (n == 0) then
      call ErrorHandler('computeAStarT','sqrt[nc] is not a integer',nc)
   endif
   lmc = n - 1
!
   n = 0
   LOOP_lmr: do lmr = 1, nr
      if (lmr*lmr == nr) then
         n = lmr
         exit LOOP_lmr
      endif
   enddo LOOP_lmr
   if (n == 0) then
      call ErrorHandler('computeAStarT','sqrt[nr] is not a integer',nr)
   endif
   lmr = n - 1
!
   kl = 0
   do l = 0, lmc
      do m = -l, l
         kl = kl + 1
         klp = 0
         do lp = 0, lmr
            do mp = -lp, lp
               klp = klp + 1
               ma = abs(m+mp); n = 1-2*mod(ma,2)
               B(klp,kl) = n*A(kl-2*m,klp-2*mp)
            enddo
         enddo
      enddo
   enddo
!
   end subroutine computeAStarT
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeAStarTInv(A,nr,nc,B)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nc
   integer (kind=IntKind) :: ipiv(2*nc), info
!
   complex (kind=CmplxKind), intent(in) :: A(nr,nc)
   complex (kind=CmplxKind), intent(out) :: B(nr,nc)
   complex (kind=CmplxKind) :: w(16*nc)
!
!  -------------------------------------------------------------------
   call computeAStarT(A,nr,nc,B)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
   call zgetrf( nc, nc, B, nr, ipiv, info )
!  -------------------------------------------------------------------
   if (info /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('computeAStarTInv','ZGETRF failure, B is ill conditioned',info)
!     ----------------------------------------------------------------
   endif
!  -------------------------------------------------------------------
   call zgetri( nc, B, nr, ipiv, w, 16*nc, info)
!  -------------------------------------------------------------------
   if (info /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('computeAStarTInv','ZGETRI failure, B is ill conditioned',info)
!     ----------------------------------------------------------------
   endif
!
   end subroutine computeAStarTInv
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeAprojB(p,n,A,B,C)
!  ===================================================================
!
!  If p = 'L' or 'l', compute C = [1 + A * B]^{-1} * A
!
!  If p = 'R' or 'r', compute C = A * [1 + B * A]^{-1}
!
!  *******************************************************************
   use MatrixInverseModule, only : MtxInv_LU
!
   implicit none
!
   character (len=1), intent(in) :: p
!
   integer (kind=IntKind), intent(in) :: n
!
   complex (kind=CmplxKind), intent(in) :: A(n,n), B(n,n)
   complex (kind=CmplxKind), intent(out) :: C(n,n)
   complex (kind=CmplxKind) :: D(n,n)
!
!  -------------------------------------------------------------------
   call setupUMc1(n,D)
!  -------------------------------------------------------------------
!
   if (p == 'L' .or. p == 'l') then
!     ----------------------------------------------------------------
      call zgemm('n','n',n,n,n,CONE,A,n,B,n,CONE,D,n)
!     ----------------------------------------------------------------
      call MtxInv_LU(D,n)
!     ----------------------------------------------------------------
      call zgemm('n','n',n,n,n,CONE,D,n,A,n,CZERO,C,n)
!     ----------------------------------------------------------------
   else if (p == 'R' .or. p == 'r') then
!     ----------------------------------------------------------------
      call zgemm('n','n',n,n,n,CONE,B,n,A,n,CONE,D,n)
!     ----------------------------------------------------------------
      call MtxInv_LU(D,n)
!     ----------------------------------------------------------------
      call zgemm('n','n',n,n,n,CONE,A,n,D,n,CZERO,C,n)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('computeAprojB','invalid projection type',p)
!     ----------------------------------------------------------------
   endif
!
   end subroutine computeAprojB
!  ===================================================================
end module MatrixModule
