      subroutine repl(x,y,n,ndim)
c====================
c
c    x=y
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=1,n
    1 x(i,j)=y(i,j)
      return
      end
      subroutine replt(x,y,n,ndim)
c====================
c
c    x=y^T
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=1,n
    1 x(i,j)=y(j,i)
      return
      end
      subroutine replrel(x,y,n,ndim)
c=======================
c
c    x(k,k')=(-)**(l-l')* y(k',k)*
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
      dimension ldex(50)
c
      data ldex/0,0,
     *          1,1,1,1,1,1,
     *          2,2,2,2,2,2,2,2,2,2,
     *          3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     *          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/
c
      do 1 i=1,n
        li=ldex(i)
      do 1 j=1,n
        lj=ldex(j)
        x(i,j)=dconjg(y(j,i))
        if(mod(iabs(li-lj),2).eq.1) x(i,j)=-x(i,j)
    1 continue
      return
      end
      subroutine compmat(x,y,n,ndim,tol,ic)
c=======================
c
c check whether  x=y
c
      implicit real*8 (a-h,o-z)
      complex*16 x(ndim,ndim),y(ndim,ndim),xx,dx
      data tol0/1.0d-12/
c
      if(tol.le.0.d0) tol=tol0
      ic=0
      xnorm=0.d0
      dxnorm=0.d0
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
      return
      end
      subroutine addmat(x,y,n,ndim)
c====================
c
c    x=x+y
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=1,n
    1 x(i,j)=x(i,j)+y(i,j)
      return
      end
      subroutine addmat1(a,b,c,n,ndim)
c====================
c
c    c=a+b
c
      complex*16 a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=1,n
    1 c(i,j)=a(i,j)+b(i,j)
      return
      end
      subroutine submat(x,y,n,ndim)
c====================
c
c    x=x-y
c
      complex*16 x(ndim,ndim),y(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=1,n
    1 x(i,j)=x(i,j)-y(i,j)
      return
      end
      subroutine submat1(a,b,c,n,ndim)
c====================
c
c    c=a-b
c
      complex*16 a(ndim,ndim),b(ndim,ndim),c(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=1,n
    1 c(i,j)=a(i,j)-b(i,j)
      return
      end
      subroutine symmat(x,n,ndim)
c======================
c Symmetrize matrix x:
c    x -> 1/2*(x+xT)
c
      complex*16 x(ndim,ndim)
c
      do 1 i=1,n
      do 1 j=i+1,n
      x(i,j)=0.5d0*(x(i,j)+x(j,i))
    1 x(j,i)=x(i,j)
      return
      end
      subroutine doubmt(amat,bmat,ndim,ndimp)
c======================
c
c  amat=amat*bmat
c
      implicit real*8 (a-h,o-z)
c
      complex*16 amat(ndimp,ndimp),bmat(ndimp,ndimp)
      complex*16 help(100),sum
c
      do i=1,ndim
        do k=1,ndim
          sum=dcmplx(0.d0,0.d0)
          do j=1,ndim
            sum=sum+bmat(j,k)*amat(i,j)
          end do
          help(k)=sum
        end do
        do k=1,ndim
          amat(i,k)=help(k)
        end do
      end do
c
      return
      end
      subroutine doubmt1(amat,bmat,cmat,ndim,ndimp)
c======================
c
c  cmat=amat*bmat
c
      implicit real*8 (a-h,o-z)
c
      complex*16 amat(ndimp,ndimp),bmat(ndimp,ndimp),cmat(ndimp,ndimp)
      complex*16 sum
c
      do i=1,ndim
        do k=1,ndim
          sum=dcmplx(0.d0,0.d0)
          do j=1,ndim
            sum=sum+bmat(j,k)*amat(i,j)
          end do
          cmat(i,k)=sum
        end do
      end do
c
      return
      end
      function trdbmt(amat,bmat,ndim,ndimp)
c===============================
c
c  trdbmt=Tr(amat*bmat)
c
      integer i,j,ndim,ndimp
      complex*16 trdbmt,amat(ndimp,ndimp),bmat(ndimp,ndimp)
c
      trdbmt=dcmplx(0.d0,0.d0)
      do i=1,ndim
      do j=1,ndim
        trdbmt=trdbmt+bmat(j,i)*amat(i,j)
      end do
      end do
c
      return
      end
      subroutine tripmt(u,b,ust,ndi1,ndi2,ndim)
c =====================
c
c vectorized routine for triple product of rectangular matrices
c
      implicit real*8 (a-h,o-z)
      parameter(nndim=100)
      complex*16 u,ust,b,c,x
      dimension u(ndim,ndim),ust(ndim,ndim),b(ndim,ndim)
      dimension c(nndim,nndim)
c
c     left product
c
      do 20 i=1,ndi1
      do 20 j=1,ndi2
c
      x= (0.0d0,0.0d0)
c
      do 10 k=1,ndi1
   10 x = x + b(k,j)*u(i,k)
c
      c(i,j)=x
   20 continue
c
c     right product
c
      do 40 i=1,ndi1
      do 40 j=1,ndi2
c
      x = (0.d0,0.0d0)
c
      do 30 k=1,ndi2
   30 x = x + ust(k,j)*c(i,k)
c
      b(i,j)=x
   40 continue
c
      return
      end
      subroutine tripmt1(u,b,ust,b1,ndi1,ndi2,ndim)
c =====================
c
c vectorized routine for triple product of rectangular matrices
c
      implicit real*8 (a-h,o-z)
      parameter(nndim=100)
      complex*16 u,ust,b,b1,c,x
      dimension u(ndim,ndim),ust(ndim,ndim)
      dimension b(ndim,ndim),b1(ndim,ndim)
      dimension c(nndim,nndim)
c
c     left product
c
      do 20 i=1,ndi1
      do 20 j=1,ndi2
c
      x= (0.0d0,0.0d0)
c
      do 10 k=1,ndi1
   10 x = x + b(k,j)*u(i,k)
c
      c(i,j)=x
   20 continue
c
c     right product
c
      do 40 i=1,ndi1
      do 40 j=1,ndi2
c
      x = (0.d0,0.0d0)
c
      do 30 k=1,ndi2
   30 x = x + ust(k,j)*c(i,k)
c
      b1(i,j)=x
   40 continue
c
      return
      end
      subroutine outmat(mat,n,m,ndim,nper)
c=====================
      implicit real*8 (a-h,o-z)
      complex*16 mat(ndim,ndim)
      complex*16 mat1(50,50)
c
      do i=1,n
       do j=1,m
        r1=dreal(mat(i,j))
        r2=dimag(mat(i,j))
        if(dabs(r1).lt.1.d-15) r1=0.d0
        if(dabs(r2).lt.1.d-15) r2=0.d0
        mat1(i,j)=dcmplx(r1,r2)
       end do
      end do
c
      write(nper,*) ' real part'
      do 1 i=1,n
   1  write(nper,10) (dreal(mat1(i,j)),j=1,m)
      write(nper,*) ' imaginary part'
      do 2 i=1,n
   2  write(nper,10) (dimag(mat1(i,j)),j=1,m)
c
  10  format(9(1pd14.6))
c
      return
      end
      subroutine outmat1(mat,n,m,ndim,tol,nper)
c=======================
      implicit real*8 (a-h,o-z)
      complex*16 mat(ndim,ndim)
c
      do j=1,m
      do 1 i=1,n
        r1=dreal(mat(i,j))
        r2=dimag(mat(i,j))
        if(dabs(r1).lt.tol.and.dabs(r2).lt.tol) goto 1 
        write(nper,'(2i4,5x,1pd20.10,1pd20.10)') i,j,r1,r2
   1  continue
      end do
c
      return
      end
