c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bessj(nn,x,bj,bn)
c     ================================================================
c
c     computes spherical bessel function jl
c
      implicit   none
c
      integer    nn,l,k,i,nm
c
      real*8     flm,fll
c
      real*8     a,aj,x
      real*8     as,ac,x1,x2
c
      real*8     bj(nn),bn(nn)
      real*8     tol
      parameter  (tol=1.d-10)
c
      x1 = 1.d0/x
      if(abs(x1).lt.tol) then
        do l=1,nn
          bj(l)=0.d0
        enddo
        return
      endif
      as = sin(x)
      ac = cos(x)
      bj(1) = as*x1
c Forward recursion for small L
c (using bn as scratch space)
1     k=(int(abs(x))+2)/4+1
      bn(1)=as
      if(k.gt.1) then
	 if (k.gt.nn) k=nn
	 bn(2)=bj(1)-ac
	 bj(2)=bn(2)*x1
	 do l = 3, k
           flm = l+l-3
	   bn(l) = flm*bj(l-1)-bn(l-2)
	   bj(l)=bn(l)*x1
	 enddo
      endif

      if(k.lt.nn) then
c Backward recursion from very large L down to L=k
        a=bn(k)
        nm=nn*4
        aj=2*nm+1
        x2=x*x
        do l = nm, nn+2, -1
	   aj=(2*l-1)-x2/aj
        enddo
        bn(nn)=(2*nn+1)*aj*x1-x
        bj(nn)=bn(nn)*x1
        bn(nn-1)=(2*nn-1)*bj(nn)-aj
        do l = nn-1, k+1, -1
	   bj(l)=bn(l)*x1
	   bn(l-1)=(2*l-1)*bj(l)-bn(l+1)
        enddo
c scale to give correct bj(k)
        aj=a/bn(k)
        do l = k+1, nn
	   bj(l)=aj*bj(l)
        enddo
      endif

      return
      end
