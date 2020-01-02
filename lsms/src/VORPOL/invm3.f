      subroutine invm3(xm,xminv,ierr)
c     ================================================================
c Inverts a 3x3 real*8 matrix
c
c Inputs: xm   real*8 array of (3,3), the 3x3 matrix
c Returns:xminv real*8 array of (3,3), inverse of xm
c         ierr  error code, 0 if inversion is successful, 1 if xm is singular
c
      implicit   none
c
      integer    ierr
c
      real*8     xm(3,3)
      real*8     xminv(3,3)
      real*8     det
      real*8    t1,t2,t3,t0
      real*8     tol
c
      parameter  (tol=1.d-10)
c
      t1=xm(2,2)*xm(3,3)-xm(2,3)*xm(3,2)
      t2=xm(2,1)*xm(3,3)-xm(2,3)*xm(3,1)
      t3=xm(2,1)*xm(3,2)-xm(2,2)*xm(3,1)
      det= xm(1,1)*t1
     >    -xm(1,2)*t2
     >    +xm(1,3)*t3
c
      if (abs(det) .lt. tol) then
         ierr=1
      else
         ierr=0
c
	 t0=1.d0/det
         xminv(1,1)=+t1*t0
         xminv(2,1)=-t2*t0
         xminv(3,1)=+t3*t0
         xminv(1,2)=-(xm(1,2)*xm(3,3)-xm(1,3)*xm(3,2))*t0
         xminv(2,2)=+(xm(1,1)*xm(3,3)-xm(1,3)*xm(3,1))*t0
         xminv(3,2)=-(xm(1,1)*xm(3,2)-xm(3,1)*xm(1,2))*t0
         xminv(1,3)=+(xm(1,2)*xm(2,3)-xm(1,3)*xm(2,2))*t0
         xminv(2,3)=-(xm(1,1)*xm(2,3)-xm(2,1)*xm(1,3))*t0
         xminv(3,3)=+(xm(1,1)*xm(2,2)-xm(2,1)*xm(1,2))*t0
      end if
c
      return
      end
