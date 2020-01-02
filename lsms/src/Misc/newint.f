      subroutine newint(nr,r,f,g,ip0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c----- this routine integrates functn f on a 1-d mesh by interpolation
c      the partial integrals (up to r(k)) are stored in b(k).  the mesh
c      is arbitrary
c      xgz  ornl 1998
c
      implicit real*8 (a-h,o-z)
c
      real*8 r(nr),cr(4)
      real*8 f(nr),g(nr)
      real*8 fr(4),rs(4)
      parameter (nj=6,njm1=nj-1)
      parameter (dnj=1.d0/nj)
      parameter (wt0=2.d0/(3.d0*nj),wt1=3.d0*wt0)
c
      ip=2*ip0+1
      pip1=1.d0/dble(ip+1)
      pip2=1.d0/dble(ip+2)
      g(1)=0.d0
      i1=1
      i2=2
      i3=3
      i4=4
      i0=1
      do i=1,4
	rs(i)=sqrt(r(i))
	fr(i)=f(i)
      enddo
      ip1=1
      ip2=4
      do i=1,nr-1
	ip1=ip1+1

	if(i.gt.2.and.ip2.le.nr) then
	  rs(i4)=sqrt(r(ip2))
          fr(i4)=f(ip2)
	endif
        call fit(rs,fr,4,i0,cr)
c First integrate a linear frunction going through f(i) and f(ip1)
c A4=fi-A1
	rip1=rs(mod(i0,4)+1)
	gj=1.d0
	fj=rip1
	do j=1,ip
	  gj=gj*rs(i0)+fj
	  fj=fj*rip1
	enddo
	fj=(gj*rs(i0)+fj)*pip2
	gj=gj*pip1
	gr=(fr(i0)-cr(3)*rs(i0))*gj+cr(3)*fj
c Using Simpson's rule. No contribution from j=0 and j=nj.
	dr=rip1-rs(i0)
        h1=dnj*dr
	h2=0.d0
	wt=wt0
	do j=1,njm1
	  h2=h2+h1
	  wt=wt1-wt
	  fj=wt
	  rj=rs(i0)+h2
	  do k=1,ip
	    fj=fj*rj
	  enddo
	  gj=h2*(dr-h2)
	  c1=(cr(3)-cr(2))*h2+cr(1)
	  c2=cr(4)*gj
	  gr=gr-fj*c1*c2/(1.d0+c2)
	enddo
	g(ip1)=g(i)+2.d0*dr*gr
	if(i.gt.1.and.ip2.lt.nr) then
	  ip2=ip2+1
	  j=i1
	  i1=i2
	  i2=i3
	  i3=i4
	  i4=j
	endif
	i0=mod(i0,4)+1
      enddo
      return
      end
