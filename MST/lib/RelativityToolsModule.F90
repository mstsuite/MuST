module RelativityToolsModule

use ErrorHandlerModule, only : StopHandler
use KindParamModule, only : IntKind, RealKind, CmplxKind
use MathParamModule, only : SQRTm1, ZERO, ONE, TWO, FOUR, PI, &
                            HALF, CZERO, CONE, THREE

contains

!=====================================================================
   subroutine csbf(n,p,r,fb,fn,fh)
!=====================================================================

!  ****************************************************************
!  spherical bessel and neuman functions
!  ****************************************************************

   implicit none

   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: l
 
   real (kind=RealKind), intent(in) :: r

   complex (kind=CmplxKind), intent(out) :: fb(:)
   complex (kind=CmplxKind), intent(out) :: fn(:)
   complex (kind=CmplxKind), intent(out) :: fh(:)
   complex (kind=CmplxKind), intent(in) :: p
   complex (kind=CmplxKind) :: x1
   complex (kind=CmplxKind) :: x2

   x1=r*p
   x2=x1*x1

   fb(1)= sin(x1)/x1
   fb(2)= sin(x1)/x2 - cos(x1)/x1
   fn(1)=-cos(x1)/x1
   fn(2)=-cos(x1)/x2 - sin(x1)/x1
   fh(1)=fb(1)+sqrtm1*fn(1)
   fh(2)=fb(2)+sqrtm1*fn(2)

   do l=2,n
      fb(l+1)=(2*l-1)*fb(l)/x1-fb(l-1)
   end do
   do l=2,n
      fn(l+1)=(2*l-1)*fn(l)/x1-fn(l-1)
   end do
   do l=2,n
      fh(l+1)=fb(l+1)+sqrtm1*fn(l+1)
   end do

   return
   end subroutine csbf

!======================================================================
   subroutine gjinv(a,n,nmax,detl)
!======================================================================

   implicit none

   complex (kind=CmplxKind), intent(inout) :: a(:,:)
   complex (kind=CmplxKind), parameter :: cpi = SQRTm1*PI
   complex (kind=CmplxKind) :: q
   complex (kind=CmplxKind), intent(out) :: detl

   real (kind=RealKind) :: amaxi, atest

   integer (kind=IntKind), intent(in) :: n, nmax
   integer (kind=IntKind) :: ind(10000)
   integer (kind=IntKind) :: i, j, k

!  gauss-jordan inversion with partial pivoting

   detl = CZERO
   do 1 i=1,n
      k=i
      amaxi=cdabs(a(i,i))

      do 2 j=i+1,n
      atest=cdabs(a(i,j))
      if(atest.le.amaxi) goto 2
      k=j
      amaxi=atest
 2    continue

      ind(i)=k
      if(k.eq.i) goto 4
      detl=detl-cpi
      do 3 j=1,n
      q=a(j,i)
      a(j,i)=a(j,k)
      a(j,k)=q
 3    continue

 4    q=a(i,i)
      detl=detl+cdlog(q)

      q=CONE/q
      a(i,i)=CONE
      do 5 j=1,n
      a(j,i)=a(j,i)*q
 5    continue

      do 6 j=1,n
      if(j.eq.i) goto 6
      q=a(i,j)
      a(i,j)=CZERO
      do 7 k=1,n
      a(k,j)=a(k,j)-a(k,i)*q
 7    continue
 6    continue
 1 continue

   do 8 i=n,1,-1
      j=ind(i)
      if(i.eq.j) goto 8
      do 9 k=1,n
      q=a(i,k)
      a(i,k)=a(j,k)
      a(j,k)=q
 9    continue
 8 continue

   end subroutine gjinv

!========================================================
   subroutine repl(x,y,n,ndim)
!========================================================

!  x=y

   implicit none

   integer (kind=IntKind) :: i, j
   integer (kind=IntKind), intent(in) :: n, ndim
   complex (kind=CmplxKind), intent(out) :: x(:,:)
   complex (kind=CmplxKind), intent(in) :: y(:,:)

   do i=1,n
      do j=1,n
         x(i,j)=y(i,j)
      enddo
   enddo

   end subroutine repl

!========================================================
   subroutine replt(x,y,n,ndim)
!========================================================

!  x=y^T

   implicit none

   integer (kind=IntKind) :: i, j
   integer (kind=IntKind), intent(in) :: n, ndim
   complex (kind=CmplxKind), intent(out) :: x(:,:)
   complex (kind=CmplxKind), intent(in) :: y(:,:)

   do i=1,n
      do j=1,n
         x(i,j)=y(j,i)
      enddo
   enddo

   end subroutine replt

!=============================================================
   subroutine replrel(x,y,n,ndim)
!=============================================================
 
!  x(k,k')=(-)**(l-l')* y(k',k)*
 
   implicit none
 
   integer (kind=IntKind), intent(in) :: n, ndim
   complex (kind=CmplxKind), intent(out) :: x(:,:)
   complex (kind=CmplxKind), intent(in) :: y(:,:)
   integer (kind=IntKind) :: ldex(50), i, j, li, lj
 
   data ldex/0,0,                         &
             1,1,1,1,1,1,                 &
             2,2,2,2,2,2,2,2,2,2,         &
             3,3,3,3,3,3,3,3,3,3,3,3,3,3, &
             4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/

   do 1 i=1,n
      li=ldex(i)
      do j=1,n
         lj=ldex(j)
         x(i,j)=dconjg(y(j,i))
         if(mod(iabs(li-lj),2).eq.1) x(i,j)=-x(i,j)
      enddo
 1 continue    

   end subroutine replrel

!=================================================================
   subroutine compmat(x,y,n,ndim,tol,ic)
!=================================================================

!  check whether  x=y

   implicit none

   integer (kind=IntKind), intent(in) :: n, ndim
   integer (kind=IntKind), intent(out) :: ic
   integer (kind=IntKind) :: i, j

   real (kind=RealKind), intent(inout) :: tol
   real (kind=RealKind), parameter :: tol0 = 1.0d-12
   real (kind=RealKind) :: xnorm, dxnorm, xxre, xxim, dxre, dxim

   complex (kind=CmplxKind), intent(in) :: x(:,:), y(:,:)
   complex (kind=CmplxKind) :: xx, dx

   if(tol.le.ZERO) tol=tol0
   ic=0
   xnorm=ZERO
   dxnorm=ZERO
   do i=1,n
      do j=1,n
        xx=x(i,j)
        dx=x(i,j)-y(i,j)
        xxre=dreal(xx)
        xxim=dimag(xx)
        dxre=dreal(dx)
        dxim=dimag(dx)
        xnorm=xnorm+xxre*xxre+xxim*xxim
        dxnorm=dxnorm+dxre*dxre+dxim*dxim
      end do
   end do
   xnorm=dsqrt(xnorm)
   dxnorm=dsqrt(dxnorm)
   if(dxnorm.lt.tol*xnorm) ic=1

   end subroutine compmat

!===========================================================
   subroutine addmat(x,y,n,ndim)
!===========================================================

!  x=x+y

   implicit none

   integer (kind=IntKind), intent(in) :: n, ndim
   integer (kind=IntKind) :: i, j
   complex (kind=CmplxKind), intent(inout) :: x(:,:)
   complex (kind=CmplxKind), intent(in) :: y(:,:)

   do i=1,n
      do j=1,n
         x(i,j)=x(i,j)+y(i,j)
      enddo
   enddo

   end subroutine addmat

!============================================================
   subroutine addmat1(a,b,c,n,ndim)
!============================================================

!  c=a+b

   implicit none

   integer (kind=IntKind), intent(in) :: n, ndim
   integer (kind=IntKind) :: i, j
   complex (kind=CmplxKind), intent(in) :: a(:,:), b(:,:)
   complex (kind=CmplxKind), intent(out) :: c(:,:)

   do i=1,n
      do j=1,n
         c(i,j)=a(i,j)+b(i,j)
      enddo
   enddo

   end subroutine addmat1

!============================================================
   subroutine submat(x,y,n,ndim)
!============================================================

!  x=x-y

   implicit none

   integer (kind=IntKind), intent(in) :: n, ndim
   integer (kind=IntKind) :: i, j
   complex (kind=CmplxKind), intent(inout) :: x(:,:)
   complex (kind=CmplxKind), intent (in) :: y(:,:)

   do i=1,n
      do j=1,n
         x(i,j)=x(i,j)-y(i,j)
      enddo
   enddo

   end subroutine submat

!==============================================================
   subroutine submat1(a,b,c,n,ndim)
!==============================================================

!  c=a-b

   implicit none

   integer (kind=IntKind) :: i, j
   integer (kind=IntKind), intent(in) :: n, ndim
   complex (kind=CmplxKind), intent(in) :: a(:,:), b(:,:)
   complex (kind=CmplxKind), intent(out) :: c(:,:)

   do i=1,n
      do j=1,n
         c(i,j)=a(i,j)-b(i,j)
      enddo
   enddo

   end subroutine submat1

!==============================================================
   subroutine symmat(x,n,ndim)
!==============================================================
! Symmetrize matrix x:
! x -> 1/2*(x+xT)

   implicit none

   integer (kind=IntKind), intent(in) :: n, ndim
   integer (kind=IntKind) :: i, j
   complex (kind=CmplxKind), intent(inout) :: x(:,:)

   do i=1,n
      do j=i+1,n
         x(i,j)=HALF*(x(i,j)+x(j,i))
         x(j,i)=x(i,j)
      enddo
   enddo

   end subroutine symmat

!==============================================================
   subroutine doubmt(amat,bmat,ndim,ndimp)
!==============================================================

!  amat=amat*bmat

   implicit none

   integer (kind=IntKind), intent(in) :: ndim, ndimp
   integer (kind=IntKind) :: i, j, k
   complex (kind=CmplxKind), intent(inout) :: amat(:,:)
   complex (kind=CmplxKind), intent(in) :: bmat(:,:)
   complex (kind=CmplxKind) :: help(100), sum

   do i=1,ndim
      do k=1,ndim
         sum=CZERO
         do j=1,ndim
            sum=sum+bmat(j,k)*amat(i,j)
         end do
         help(k)=sum
      end do
      do k=1,ndim
         amat(i,k)=help(k)
      end do
   end do

   end subroutine doubmt

!================================================================
   subroutine doubmt1(amat,bmat,cmat,ndim,ndimp)
!================================================================

!  cmat=amat*bmat

   implicit none

   integer (kind=IntKind), intent(in) :: ndim, ndimp
   integer (kind=IntKind) :: i, j, k
   complex (kind=CmplxKind), intent(in) :: amat(:,:),bmat(:,:)
   complex (kind=CmplxKind), intent(out) :: cmat(:,:)
   complex (kind=CmplxKind) :: sum

   do i=1,ndim
      do k=1,ndim
         sum=CZERO
         do j=1,ndim
            sum=sum+bmat(j,k)*amat(i,j)
         end do
         cmat(i,k)=sum
      end do
   end do

   end subroutine doubmt1

!================================================================
   function trdbmt(amat,bmat,ndim,ndimp) result(tr)
!================================================================

!  trdbmt=Tr(amat*bmat)

   implicit none

   integer (kind=IntKind) :: i, j
   integer (kind=IntKind), intent(in) :: ndim, ndimp
   complex (kind=CmplxKind) :: tr
   complex (kind=CmplxKind), intent(in) :: amat(:,:), bmat(:,:)

   tr=CZERO
   do i=1,ndim
      do j=1,ndim
         tr=tr+bmat(j,i)*amat(i,j)
      end do
   end do

   end function trdbmt

!=================================================================
   subroutine tripmt(u,b,ust,ndi1,ndi2,ndim)
!=================================================================

!  vectorized routine for triple product of rectangular matrices

   implicit none

   integer (kind=IntKind), parameter :: nndim = 100
   integer (kind=IntKind), intent(in) :: ndi1, ndi2, ndim
   integer (kind=IntKind) :: i, j, k
   complex (kind=CmplxKind), intent(in) :: u(:,:), ust(:,:)
   complex (kind=CmplxKind), intent(inout) :: b(:,:)
   complex (kind=CmplxKind) :: c(nndim, nndim), x

!  left product

   do 20 i=1,ndi1
      do j=1,ndi2
         x = CZERO
         do 10 k=1,ndi1
            x = x + b(k,j)*u(i,k)
10       continue
         c(i,j)=x
      enddo
20 continue

!  right product

   do 40 i=1,ndi1
      do j=1,ndi2
         x = CZERO
         do 30 k=1,ndi2
            x = x + ust(k,j)*c(i,k)
30       continue
         b(i,j)=x
      enddo
40 continue

   end subroutine tripmt

!===================================================================
   subroutine tripmt1(u,b,ust,b1,ndi1,ndi2,ndim)
!===================================================================

! vectorized routine for triple product of rectangular matrices

   implicit none

   integer (kind=IntKind), parameter :: nndim = 100
   integer (kind=IntKind), intent(in) :: ndi1, ndi2, ndim
   integer (kind=IntKind) :: i, j, k
   complex (kind=CmplxKind), intent(in) :: u(:,:), ust(:,:) ,b(:,:)
   complex (kind=CmplxKind), intent(out) :: b1(:,:)
   complex (kind=CmplxKind) :: c(nndim,nndim), x

!  left product

   do 20 i=1,ndi1
      do j=1,ndi2
      x = CZERO
      do 10 k=1,ndi1
         x = x + b(k,j)*u(i,k)
10    continue
      c(i,j)=x
      enddo
20 continue

!  right product

   do 40 i=1,ndi1
      do j=1,ndi2
      x = CZERO
      do 30 k=1,ndi2
         x = x + ust(k,j)*c(i,k)
30    continue
      b1(i,j)=x
      enddo
40 continue

   end subroutine tripmt1

!======================================================================
   subroutine outmat(mat,n,m,ndim,nper)
!======================================================================
   implicit none
      
   integer (kind=IntKind), intent(in) :: n, m, ndim, nper
   integer (kind=IntKind) :: i, j
   real (kind=RealKind) :: r1, r2
   complex (kind=CmplxKind), intent(in) :: mat(:,:)
   complex (kind=CmplxKind) :: mat1(50,50)

   do i=1,n
      do j=1,m
         r1=dreal(mat(i,j))
         r2=dimag(mat(i,j))
         if(dabs(r1).lt.1.d-15) r1=0.d0
         if(dabs(r2).lt.1.d-15) r2=0.d0
         mat1(i,j)=dcmplx(r1,r2)
      end do
   end do

   write(nper,*) ' real part'
   do 1 i=1,n
      write(nper,10) (dreal(mat1(i,j)),j=1,m)
1  continue
   write(nper,*) ' imaginary part'
   do 2 i=1,n
      write(nper,10) (dimag(mat1(i,j)),j=1,m)
2  continue
10 format(9(1pd14.6))

   end subroutine outmat

!=====================================================================
   subroutine outmat1(mat,n,m,ndim,tol,nper)
!=====================================================================
   implicit none

   integer (kind=IntKind), intent(in) :: n, m, ndim, nper
   integer (kind=IntKind) :: i, j
   real (kind=RealKind), intent(in) :: tol
   real (kind=RealKind) :: r1, r2
   complex (kind=CmplxKind), intent(in) :: mat(ndim,ndim)

   do j=1,m
      do 1 i=1,n
         r1=dreal(mat(i,j))
         r2=dimag(mat(i,j))
         if(dabs(r1).lt.tol.and.dabs(r2).lt.tol) goto 1 
         write(nper,'(2i4,5x,1pd20.10,1pd20.10)') i,j,r1,r2
1     continue
   end do

   end subroutine outmat1

!======================================================================
   subroutine matr(tlm,lmax,dmat,dmat1)
!======================================================================

!  R=(tvec,phi)
!  dmat = D(R) and dmat1 = D(R)+ 

   implicit none

   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: kmax, kmymax, i, j, ist, icase
   integer (kind=IntKind) :: m1, m2, il, j2
   real (kind=RealKind), intent(in) :: tlm(:) ! dimension must = 4
   real (kind=RealKind) :: trd
   complex (kind=CmplxKind), intent(out) :: dmat(:,:)
   complex (kind=CmplxKind), intent(out) :: dmat1(:,:)
   complex (kind=CmplxKind) :: cmat(2*(lmax+1)*(lmax+1),2*(lmax+1)*(lmax+1))
   complex (kind=CmplxKind) :: d(2*lmax+2,2*lmax+2)

   kmax=2*lmax+1
   kmymax=2*(lmax+1)*(lmax+1)

!  Set up matrix of rotation

   do i=1,kmymax
      do j=1,kmymax
         dmat(i,j)=(0.d0,0.d0)
      end do
   end do
   ist=0
   do j2=1,2*lmax-1,2
      il=j2+1
      call rotmat(d,trd,il,tlm,1,lmax)
      do icase=1,2
         do m1=1,il
            do m2=1,il
               dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
            end do
         end do
         ist=ist+il
      end do
   end do
   j2=2*lmax+1
   il=j2+1
   call rotmat(d,trd,il,tlm,1,lmax)
   do m1=1,il
      do m2=1,il
!        fill up matrix according to ascending m index
         dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
      end do
   end do

   do i=1,kmymax
      do j=1,kmymax
         dmat1(i,j)=dconjg(dmat(j,i))
      end do
   end do

   if(0.lt.4) return
   write(6,'(/'' Matrix of rotation'')') 
   call outmat(dmat,kmymax,kmymax,2*(lmax+1)*(lmax+1),6)
   write(6,'(/'' Inverse'')') 
   call outmat(dmat1,kmymax,kmymax,2*(lmax+1)*(lmax+1),6)
   call repl(cmat,dmat,kmymax,2*(lmax+1)*(lmax+1))
   call doubmt(cmat,dmat1,kmymax,2*(lmax+1)*(lmax+1))
   write(6,'(/'' D(R) * D(R**-1) '')')
   call outmat(cmat,kmymax,kmymax,2*(lmax+1)*(lmax+1),6)

   end subroutine matr

!=================================================================
   subroutine rotmat (d,tr,ij0,lm,ig,lmax)
!=================================================================

!  rotation matrix for spherical harmonics with angular momentum
!  number j
!  input: ij0 = 2j+1 (j may be integer or half-integer)
!         lm ... lm(1)=cos(phi/2)
!                lm(2..4)=sin(phi/2)*n(x..z), where n=normalized
!                direction vector
!         ig ... >0 proper rotations, <0 improper rotations
!  output: d ... rotation matrix in condon-shortley convention
!          tr ... trace of d

   implicit none

   complex (kind=CmplxKind), intent(out) :: d(:,:)
   complex (kind=CmplxKind) :: f1,f2,f3,f4,fk,c
   real (kind=RealKind), intent(in) :: lm(4)
   real (kind=RealKind), intent(out) :: tr
   real (kind=RealKind) :: fac(2*lmax+4)
   integer (kind=IntKind) :: jdim, i, ij, im, in, k, km, kx, k0
   integer (kind=IntKind), intent(in) :: ij0, ig, lmax

   jdim=2*lmax+2

   fac(1)=ONE
   do i=1,jdim+1
      fac(i+1)=fac(i)*i
   enddo
   f1=dcmplx(lm(1),-lm(4))
   f2=dconjg(f1)
   f3=dcmplx(lm(3),-lm(2))
   f4=dconjg(f3)
   ij=iabs(ij0)
   do in=1,ij
      do im=1,ij
         km=max0(1,in-im+1)
         kx=min0(ij-im+1,in)
         d(in,im)=CZERO
         do k=km,kx
            k0=k-1
            fk=CONE
            if (ij-im-k0.eq.0) go to 20
            fk=f1**(ij-im-k0)
20          if (in-k.eq.0) go to 30
            fk=fk*f2**(in-k)
30          if (k0.eq.0) go to 40
            if (lm(2).eq.ZERO) fk=fk*lm(3)**k0
            if (lm(2).ne.ZERO) fk=fk*f3**k0
40          if (im-in+k0.eq.0) go to 50
            if (lm(2).eq.ZERO) fk=fk*(-lm(3))**(im-in+k0)
            if (lm(2).ne.ZERO) fk=fk*(-f4)**(im-in+k0)
50          d(in,im)=d(in,im)+fk/(fac(ij-im-k0+1)*fac(in-k0)*fac(k)* &
            fac(im-in+k))
         enddo
         d(in,im)=d(in,im)*dsqrt(fac(ij-in+1)*fac(in)*fac(ij-im+1)* &
         fac(im))
         if (mod(ij,4).eq.3.and.isign(1,ig).ne.1) d(in,im)=-d(in,im)
         if (ij0.lt.0.and.isign(1,ig).ne.1) d(in,im)=-d(in,im)
         if (dabs(dreal(d(in,im))).lt.1.d-14) d(in,im)=dcmplx(ZERO, &
         dimag(d(in,im)))
         if (dabs(dimag(d(in,im))).lt.1.d-14) d(in,im)=dreal(d(in,im))
      enddo
   enddo
   c=CZERO
   do i=1,ij
      c=c+d(i,i)
   enddo
   if (dabs(dimag(c)).gt.1.d-10) stop 10
   tr=c

   end subroutine rotmat

!======================================================================
   function rsimp(f,r,irn,dx) result(simp)
!======================================================================
!  radial integration via simpson

   implicit none

   integer (kind=IntKind) :: i, isw, nl, np
   integer (kind=IntKind), intent(in) :: irn
   real (kind=RealKind) :: s, simp
   real (kind=RealKind), intent(in) :: f(:), r(:), dx

   isw=0
   simp=ZERO
   if(irn.le.2) return
   if(irn/2*2.eq.irn) isw=1
   np=irn-isw
   s=f(1)*r(1)+f(np)*r(np)
   nl=np-1
   do i=2,nl,2
      s=s+FOUR*f(i)*r(i)
   enddo
   nl=nl-1
   if(nl.lt.3) goto 15
   do i=3,nl,2
      s=s+TWO*f(i)*r(i)
   enddo
15 s=s*dx/THREE
   if(isw.eq.1) goto 30
   simp=s
   return
30 simp=s+(f(irn)*r(irn)+f(irn-1)*r(irn-1))*HALF*dx

   end function rsimp

!======================================================================
   subroutine getclm(lmax,clm)
!======================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: m2
   integer (kind=IntKind) :: i

   real (kind=RealKind), intent(out) :: clm(:)
   real (kind=RealKind) :: fpi
   real (kind=RealKind) :: xfac

!  *****************************************************************
!  Coeficients for complex spherical harmonics......................
!  Calclates all the c(l,m)'s up to lmax............................
!     
!               [ (2*l+1)*(l-|m|)!]
!  c(l,m)=  sqrt[-----------------]
!               [   4*pi*(l+|m|)! ]
!
!  *****************************************************************
!
!  =================================================================
   if(lmax.lt.0 ) then
      write(6,'('' GETCLM:: bad arguments: lmax='',i5)') lmax
      stop 'getclm'
   endif
!  =================================================================
!  begin calculation of c(l,m)'s....................................
!  =================================================================
   fpi=four*PI
!  =================================================================
!  special case lmax=0..............................................
!  =================================================================
   clm(1)=sqrt(one/fpi)
   do l=1,lmax
      xfac=sqrt((2*l+1)/fpi)
      do m=0,l
         clm(l*(l+1)/2+m+1)=one
         m2=2*m
         do i=1,m2
            clm(l*(l+1)/2+m+1)=(l-m+i)*clm(l*(l+1)/2+m+1)
         enddo
!        ===========================================================
!        This version is consisten with (-1)**m being in Plm's......
!        See routine plglmax........................................
!        ===========================================================
         clm(l*(l+1)/2+m+1)=xfac*sqrt(one/clm(l*(l+1)/2+m+1))
      enddo
   enddo

   end subroutine getclm

!=====================================================================
   subroutine plglmax(lmax,x,plm)
!=====================================================================

   implicit none

   character (len=20), parameter :: sname = 'plglmax'

   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l,ll
   integer (kind=IntKind) m,mm
   integer (kind=IntKind) i

   real (kind=RealKind), intent(out) :: plm(:)
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind) :: pmm
   real (kind=RealKind) :: somx2
   real (kind=RealKind) :: fact
   real (kind=RealKind), parameter :: tol = 0.1d0**8

!  ****************************************************************
!  Associated Legendre function....................................
!  Calclates all the p(l,m)'s up to lmax...........................
!  based on the formulae given in "Numerical Recipes" pages 180-183
!  (Equations 6.6.7, 6.6.8 and 6.6.9...............................
!  W. H. Press, B. P. Flannery, S A Teukolsky and W. T. Vetterling.
!  Cambridge Univ Press 1986.......................................
!
!  N.B. The definition of p(l,m) has been modified.................
!  p(l,m) of this code = [(-1)**m]*p(l,m) of "Numerical Recipes"...
!
!  N.N.B. The definition of p(l,m) has been returned to that of NM
!  i.e the factor (-1)**m of NM is now in this routine............
!  The routine to generate Clm's for calculation of Ylm's has also
!  been modified by removing the factor of (-1)**m................
!  ****************************************************************

!  ================================================================
   if(lmax.lt.0 .or. abs(x).gt.(one+tol)) then
      write(6,'('' plgndr:: bad arguments: lmax='',i5, &
                '' x='',d14.6)') lmax,x
      call StopHandler(sname)
   endif

   if((one-abs(x)).le.tol) then
!     -------------------------------------------------------------
      call zeroout(plm,(lmax+1)*(lmax+2)/2)
!     -------------------------------------------------------------
      if(x.lt.zero) then
         do l=0,lmax
            i=(l+1)*(l+2)/2-l
            plm(i)=one-two*mod(l,2)
         enddo
      else
         do l=0,lmax
            i=(l+1)*(l+2)/2-l
            plm(i)=one
         enddo
      endif
      return
   endif

!  ================================================================
!  begin calculation of p(l,m)'s...................................
   if(lmax.eq.0) then
!     =============================================================
!     special case lmax=0..........................................
      plm(1)=one
   else
!     =============================================================
!     minus sign added to be consistant with Numerical Recipes
!     which has (-1)^m factor in plm :  July 97  by xgz.............
!     =============================================================
      somx2=-sqrt((one-x)*(one+x))
      if(lmax.eq.1) then
!        ==========================================================
!        special case lmax=1.......................................
         plm(1)=one
         plm(2)=x
         plm(3)=somx2
      else
         do m=0,lmax
!           =======================================================
!                                    m       m
!           calculate the first two P   and P
!                                    m       m+1
!           =======================================================
            if(m.eq.0) then
               plm(1)=one
               plm(2)=x
            else
               pmm=somx2
               fact=one
               do i=2,m
                  fact=fact+two
                  pmm=pmm*fact*somx2
               enddo
               mm=(m+1)*(m+2)/2
               plm(mm)=pmm
               if( mm+m+1.le.(lmax+1)*(lmax+2)/2 ) then
                  plm(mm+m+1)=x*(2*m+1)*pmm
               end if
            endif
!           =======================================================
!                               m        m
!           calculate the rest P     to P
!                               m+2      lmax
!           =======================================================
            ll=(m+2)*(m+1)/2
            fact=(two*m+one)*x
            do l=m+2,lmax
               pmm=(l+m-1)*plm(ll)
               fact=fact+two*x
               ll=ll+l-1
               plm(ll+l)=( fact*plm(ll) - pmm )/dble(l-m)
            enddo
         enddo
      endif
   endif
  
   end subroutine plglmax

!=====================================================================
   function gaunt(l1,m1,l2,m2,l3,m3,lmax,ngauss, &
                  wg,plmg,clm) result(gnt)
!=====================================================================

!            4pi  _  m1 _     m2* _     m3 _
!    gaunt = int do Y  (o) * Y   (o) * Y  (o)
!             0      l1       l2        l3   
!
!             L3
!          = C
!             L1,L2
!
!  ****************************************************************
!
!  Inputs: l1,m1,l2,m2,l3,m3   integer scalars, (l,m) indices of the gaunt
!                              number
!          lmax     integer scalar, max l values allowed, lmax >= max(l1,l2,l3)
!                   lmax also determines the array dimensions of plmg and clm
!          ngauss   integer scaler, Number of gaussian points for integration
!          wg       real*8 array of (ngauss), weights for mesh points
!                   wg is output of gauleg
!          plmg     real*8 array of ((lmax+1)*(lmax+2)/2,ngauss), where
!                   l=max(l1,l2,l3), the values of the associate Legendre
!                   function at the mesh points
!                   plmg is output of plglmax
!          clm      real*8 array of ((lmax+1)*(lmax+2)/2), prefactors for the
!                   spherical harmonics
!                   clm is output of getclm

   implicit none

   integer (kind=IntKind), intent(in) :: l1, l2, l3
   integer (kind=IntKind), intent(in) :: m1, m2, m3
   integer (kind=IntKind), intent(in) :: lmax, ngauss
   integer (kind=IntKind) :: ng
   integer (kind=IntKind) :: ifac1, ifac2, ifac3
   integer (kind=IntKind) :: jl1, jl2, jl3

   real (kind=RealKind), intent(in) :: wg(:)
   real (kind=RealKind), intent(in) :: plmg(:,:)
   real (kind=RealKind), intent(in) :: clm((lmax+1)*(lmax+2)/2)
   real (kind=RealKind) :: gnt

   real (kind=RealKind) :: fourpi = FOUR*PI

   if(l1.gt.lmax .or. l2.gt.lmax .or. l3.gt.lmax ) then
      write(6,'(''gaunt:: bad parameters: l1,l2,l3,lmax'',4i5)') &
                                          l1,l2,l3,lmax
      call StopHandler('gaunt')
   endif
   gnt=zero
   if(mod(l1+l2+l3,2).ne.0) then
      return
   else
      if(m1-m2+m3 .ne. 0) then
         return
      else
         if(l1+l2.lt.l3 .or. l2+l3.lt.l1 .or. l3+l1.lt.l2) then
            return
         else
!           -------------------------------------------------------
            call defac(l1, m1,ifac1,jl1)
            call defac(l2,-m2,ifac2,jl2)
            call defac(l3, m3,ifac3,jl3)
!           -------------------------------------------------------
            do ng=1,ngauss
               gnt=gnt+wg(ng)*plmg(jl1,ng)* &
                    plmg(jl2,ng)*           &
                    plmg(jl3,ng)
            enddo
            if(abs(gnt).lt.1.d-14) then
               gnt=0.d0
            else
               gnt=fourpi*gnt*ifac1*clm(jl1)* &
                    (1-2*mod(abs(m2),2))*ifac2*clm(jl2)*ifac3*clm(jl3)
            endif
         endif
      endif
   endif

   end function gaunt

!=====================================================================
   subroutine defac(l,m,ifac,jl)
!=====================================================================

   implicit none

   integer (kind=IntKind), intent(in) :: l, m
   integer (kind=IntKind), intent(out) :: jl
   integer (kind=IntKind), intent(out) :: ifac

   if (m.ge.0) then
      ifac= 1
      jl = (l+1)*(l+2)/2-l+m
   else
      ifac= (1-2*mod(abs(m),2))
      jl = (l+1)*(l+2)/2-l-m
   end if

   end subroutine defac

!=====================================================================
   subroutine rotr(drot,tvec,phi)
!=====================================================================
!  input:
!        tvec   normal vector of axis
!        phi    angle of rotation
!  output:
!        drot   matrix of rotation

   implicit none

   real (kind=RealKind), intent(out) :: drot(3,3)
   real (kind=RealKind), intent(in) :: tvec(3), phi
   real (kind=RealKind) :: sp, sp2, tx, ty, tz, tx2, ty2, tz2

   sp=dsin(phi)
   sp2=dsin(phi*HALF)
   tx=tvec(1)*sp
   ty=tvec(2)*sp
   tz=tvec(3)*sp
   tx2=tvec(1)*sp2
   ty2=tvec(2)*sp2
   tz2=tvec(3)*sp2

   drot(1,1)=ONE-TWO*(ty2*ty2+tz2*tz2)
   drot(2,2)=ONE-TWO*(tx2*tx2+tz2*tz2)
   drot(3,3)=ONE-TWO*(tx2*tx2+ty2*ty2)
   drot(1,2)=-tz+TWO*tx2*ty2
   drot(2,1)=tz+TWO*tx2*ty2
   drot(1,3)=ty+TWO*tx2*tz2
   drot(3,1)=-ty+TWO*tx2*tz2
   drot(2,3)=-tx+TWO*ty2*tz2
   drot(3,2)=tx+TWO*ty2*tz2

   end subroutine rotr

end module RelativityToolsModule
