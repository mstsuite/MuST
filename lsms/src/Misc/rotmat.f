      subroutine rotmat (d,tr,ij0,lm,ig,lmax)
! rotation matrix for spherical harmonics with angular momentum
! number j
! input: ij0 = 2j+1 (j may be integer or half-integer)
!        lm ... lm(1)=cos(phi/2)
!               lm(2..4)=sin(phi/2)*n(x..z), where n=normalized
!               direction vector
!        ig ... >0 proper rotations, <0 improper rotations
! output: d ... rotation matrix in condon-shortley convention
!         tr ... trace of d
!
      implicit real*8 (a-h,o-z)
!     include '../param.h'
!
!     parameter (jdim=2*lmaxp+2)
!
      complex*16 d(2*lmax+2,2*lmax+2),f1,f2,f3,f4,fk,c
      real*8 lm(4),fac(2*lmax+4)
      integer jdim

      jdim=2*lmax+2
!
      fac(1)=1.d0
      do 10 i=1,jdim+1
   10    fac(i+1)=fac(i)*i
      f1=dcmplx(lm(1),-lm(4))
      f2=dconjg(f1)
      f3=dcmplx(lm(3),-lm(2))
      f4=dconjg(f3)
      ij=iabs(ij0)
      do 60 in=1,ij
         do 60 im=1,ij
         km=max0(1,in-im+1)
         kx=min0(ij-im+1,in)
         d(in,im)=(0.d0,0.d0)
         do k=km,kx
            k0=k-1
            fk=(1.d0,0.d0)
            if (ij-im-k0.eq.0) go to 20
            fk=f1**(ij-im-k0)
   20       if (in-k.eq.0) go to 30
            fk=fk*f2**(in-k)
   30       if (k0.eq.0) go to 40
            if (lm(2).eq.0.d0) fk=fk*lm(3)**k0
            if (lm(2).ne.0.d0) fk=fk*f3**k0
   40       if (im-in+k0.eq.0) go to 50
            if (lm(2).eq.0.d0) fk=fk*(-lm(3))**(im-in+k0)
            if (lm(2).ne.0.d0) fk=fk*(-f4)**(im-in+k0)
   50       d(in,im)=d(in,im)+fk/(fac(ij-im-k0+1)*fac(in-k0)*fac(k)* 
     &      fac(im-in+k))
         enddo
         d(in,im)=d(in,im)*dsqrt(fac(ij-in+1)*fac(in)*fac(ij-im+1)* 
     &   fac(im))
         if (mod(ij,4).eq.3.and.isign(1,ig).ne.1) d(in,im)=-d(in,im)
         if (ij0.lt.0.and.isign(1,ig).ne.1) d(in,im)=-d(in,im)
         if (dabs(dreal(d(in,im))).lt.1.d-14) d(in,im)=dcmplx(0.d0, 
     &   dimag(d(in,im)))
   60    if (dabs(dimag(d(in,im))).lt.1.d-14) d(in,im)=dreal(d(in,im))
      c=(0.d0,0.d0)
      do 70 i=1,ij
   70    c=c+d(i,i)
      if (dabs(dimag(c)).gt.1.d-10) then
        write(*,*) "Error in rotmat |imag(c)| = ",dabs(dimag(c))
        stop
      end if
      tr=c
      return
      end
