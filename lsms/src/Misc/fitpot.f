c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fitpot(r,rv,coef,nr,jmt,iend)
c     ================================================================
c
      implicit   none
c
      integer    jmt
      integer    iend
      integer    i,j,i1
      integer    l
      integer    nr
c
      real*8     r(iend)
      real*8     rv(iend)
      real*8     coef(nr,iend)
      real*8     h1,h2,r0,dr,drv
      real*8     eqn11,eqn12,eqn21,eqn22,rhs1,rhs2,drvdr,det
c
      if(nr.gt.4.or.nr.gt.jmt) then
	write(6,'(''FITPOT:: nr too large: nr='',i3)') nr
	stop'nr'
      endif
      do j=1,jmt-1
	r0=r(j)
	dr=r(j+1)-r0
	drv=rv(j+1)-rv(j)
	drvdr=drv/dr
      if(nr.eq.2) then
	  coef(1,j)=drvdr
	  coef(2,j)=rv(j)
      else if(nr.eq.3) then
	write(6,'(''fitpot: nr=3 is no longer supported'')')
	stop'fitpot'
c
      else

c fit the potential to
c f1+(f0-f1+(f1-f0)*(r-r0)/(r1-r0)+a2*(r-r0)*(r-r1))/(1+a3*(r-r0)**2)
	i1=j-1
	if(j.eq.1) i1=3
	if(j.eq.jmt-1) i1=i1-1
	h1=r(i1)-r0
	h2=r(i1)-r(j+1)
	eqn11=h2
	eqn12=-h1*(rv(i1)-rv(j+1))
	rhs1=(rv(i1)-rv(j))/h1-drvdr
	i1=i1+1
	if(i1.eq.j) i1=i1+2
	h1=r(i1)-r0
	h2=r(i1)-r(j+1)
	eqn21=h2
	eqn22=-h1*(rv(i1)-rv(j+1))
	rhs2=(rv(i1)-rv(j))/h1-drvdr
	det=eqn11*eqn22-eqn21*eqn12
	if(abs(det).gt.1.d-16) then
	  coef(2,j)=-(eqn12*rhs2-eqn22*rhs1)/det
	  coef(3,j)=(eqn11*rhs2-eqn21*rhs1)/det
c If the points are on a line or nearly a line then use linear
c interpolation
	  h1=coef(3,j)*dr*dr
          if(h1.gt.-1.d0.and.abs(h1).gt.1.d-14) then
c refit into (A1*(r-r0)+A2)/(1+A3*(r-r0)**2)+A4
	    coef(1,j)=drvdr-coef(2,j)*dr
	    coef(4,j)=rv(j+1)+coef(2,j)/coef(3,j)
	    coef(2,j)=rv(j)-coef(4,j)
	  else
            coef(1,j)=drvdr
	    coef(2,j)=rv(j)
	    coef(3,j)=0.d0
	    coef(4,j)=0.d0
	  endif
        else
          coef(1,j)=drvdr
          coef(2,j)=rv(j)
          coef(3,j)=0.d0
          coef(4,j)=0.d0
        endif
      endif
      enddo
      do j=jmt,iend
	do i=1,nr
	coef(i,j)=0.d0
	enddo
      enddo
c
      return
      end
