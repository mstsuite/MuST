      subroutine zsphbesj(lmax,x,as,ac,x1,bj,bscrat)
c_______________________________________________________________________
c calculate spherical bessel functions j_l
c input: lmax integer scalar, max l for j_l
c        x    compex*16 scalar, argument of j_l
c returns: as complex*16 scalar, sin(x)
c          ac complex*16 scalar, cos(x)
c          x1 complex*16 scalar, 1/x
c          bj complex*16 array of (0:lmax), j_l
c          bscrat complex*16 array of (0:lmax), scratch space

c  xgz ornl 1994

      implicit complex*16(a-h,o-z)
      complex*16 bj(0:lmax),bscrat(0:lmax)

      as=sin(x)
      ac=cos(x)
      x1=1.0d0/x
      bj(0)=as*x1
c   Forward recursion for small L
    1 k=ishft(int(0.75d0*(abs(dreal(x))+abs(dimag(x))))+2,-2)
      bscrat(0)=as
      if(k.ge.1) then
        if(k.gt.lmax) k=lmax
        bscrat(1)=bj(0)-ac
        bj(1)=bscrat(1)*x1
        do 710 l=2,k
          bscrat(l)=(2*l-1)*bj(l-1)-bscrat(l-2)
710       bj(l)=bscrat(l)*x1
        if(k.eq.lmax) return
      endif
c   Backward recursion from very large L down to L=k
      a=bscrat(k)
      nm=ishft(lmax,4)
      aj=2*nm+3
      x2=x*x
      do 70 l=nm,lmax+2,-1
70       aj=(2*l+1)-x2/aj
      bscrat(lmax)=(2*lmax+3)*aj*x1-x
      bj(lmax)=bscrat(lmax)*x1
      bscrat(lmax-1)=(2*lmax+1)*bj(lmax)-aj
      do 720 l=lmax-1,k+1,-1
         bj(l)=bscrat(l)*x1
720      bscrat(l-1)=(2*l+1)*bj(l)-bscrat(l+1)
c scale to give correct bj(k)
      aj=a/bscrat(k)
      do 722 l=k+1,lmax
722      bj(l)=aj*bj(l)
      return
      end
