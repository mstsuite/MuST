c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ricbes(nn,y,bj,bn,dj,dn)
c     ================================================================
c
c     computes spherical bessel functions j and n and its derivatives
c
      implicit   none
c
c
      integer    nn,l,k,i,nm
c
      real*8     flm,fll
c
      complex*16 a,aj,x,y
      complex*16 as,ac,xinv,x2
c
      complex*16 bj(nn),bn(nn),dj(nn),dn(nn)
c
c
      x=y
c
      if (abs(x).eq.0.0d0) then
         write(*,*) ' x equals to zero'
         call fstop('ricbes')
      endif
c
      as = sin(x)
      ac = cos(x)
      bj(1) = as
      dj(1) = ac
      dn(1) = exp((0.d0,1.d0)*x)
      bn(1) = (0.d0,-1.d0)*dn(1)
      if(nn.eq.1) return
      xinv = 1.d0/x
c Forward recursion for small L
1     k=(int(0.75*(abs(dreal(x))+abs(dimag(x))))+2)/4+1
      if(k.gt.1) then
	 if (k.gt.nn) k=nn
	 bj(2)=bj(1)*xinv-ac
	 do l = 3, k
           flm = l+l-3
	   bj(l) = flm*bj(l-1)*xinv-bj(l-2)
	 enddo
      endif

      if(k.lt.nn) then
c Backward recursion from very large L down to L=k
        a=bj(k)
        nm=nn*4
        aj=2*nm+1
        x2=x*x
        do l = nm, nn+2, -1
	   aj=(2*l-1)-x2/aj
        enddo
        bj(nn)=(2*nn+1)*aj*xinv-x
        bj(nn-1)=(2*nn-1)*bj(nn)*xinv-aj
        do l = nn-1, k+1, -1
	   bj(l-1)=(2*l-1)*bj(l)*xinv-bj(l+1)
        enddo
c scale to give correct bj(k)
        aj=a/bj(k)
        do l = k+1, nn
	   bj(l)=aj*bj(l)
        enddo
        bj(k)=a
      endif

c Find bn using forward recursion
      bn(2) = bn(1)*xinv -  dn(1)
      dj(2) = bj(1)-bj(2)*xinv
      dn(2) = bn(1)-bn(2)*xinv
c
c     ================================================================
c     recursion relations
c     ================================================================
      do l=3,nn
         flm = l+l-3
         bn(l) = flm*bn(l-1)*xinv-bn(l-2)
         fll = l-1
         dj(l) = bj(l-1)-fll*bj(l)*xinv
         dn(l) = bn(l-1)-fll*bn(l)*xinv
      enddo
c
      return
      end
