      subroutine initwave(srelflg,r,pot,nr,nv,ist,l,energy,y,escale,
     &     core,b,allp1)
c
      implicit none
c
c---- mode flags:
c
      logical srelflg           ! scalar relativistic if true
      logical core              ! core full relativistic if true
c
c---- input:
c
      integer l                 ! angular momentum number
      complex*16 energy         ! energy
      integer nr,ist
      real*8 pot(nr)          ! potential
      real*8 r(2)          ! radial grid
      real*8 escale
      integer nv
c
c---- output:
c
      real*8 y(nv)
c
c The fine structure constant is taken from the 1986 recommended values
c of the fundamental physical constants, Physics Today August 1992.
c c is the speed of light in Hartree.
c It is half of the speed of light in Ryd unit,
c --- xgz
      real*8 c,c2
      parameter (c=137.0359895d0,c2=c*c)
c
c---- local variables:
c
      complex*16 v0,a2,a3
      complex*16 y1,y2
      real*8 z,z0,al,alp1,allp1,b,v1,a1,rl,x
      logical rel
c
c---- the following common block is used by initwave, to comunicate with
c     dfv. It is also used by psowave of the lkkr main code for
c     the same reasons.
c
!     common /dfbits/ b,allp1
!     save /dfbits/
c
      z0=pot(2)
      if(nr.eq.4) z0=z0+pot(nr)
      z0=dble(int(abs(z0)+0.1d0))
      z=escale*z0
      rel=srelflg
      if(z.lt.0.5d0) rel=.false.
      allp1=abs(l*(l+1))
      x=r(ist)
      if(rel) then
	if(core) then
        b=c
	  allp1=l
	v0=escale*(energy-pot(1))+c2
c----- full-relativistic initial conditions
	  al=sqrt(allp1*allp1-z*z/c2)
	  rl=x**al
	  a1=al+allp1
	  a2=(-2.d0*v0*z/c2+(c2+v0)*a1/z)/(2.d0*al+1.d0)
	  a3=(-2.d0*v0*a1+c2-v0)/(b*(2.d0*al+1.d0))
          y1=rl*(1.d0+a2*x)
	  y2=rl*(a1*b/z+a3*x)
	else
        b=0.25d0/c2
c
c----- semi-relativistic initial conditions
c
c    note starting values are different from the non-relativistic case
c    the ratio of r2/r1 at the first grid point is determined by
c    kh equation 15
c
        y1=1.0d0
        y2=2.0d0*c2*(sqrt(allp1+1.0d0-z*z/c2)-1.0d0)/z
	endif
      else
        b=0.0d0
c
c----- find taylor series for v-e = -z/r - v0 -v1r
c
	v0=escale*(energy-pot(1))
	v1=2.d0*escale*pot(2)*pot(3)
	v0=v0-0.5d0*v1*(r(1)+r(2))
        al=dfloat(l)
        alp1=al+1.0d0
        a1=-z/alp1
        a2=(-z*a1-v0)/(2.0d0*al+3.0d0)
        a3=-(-a1*(z*z-(3.0d0*al+4.0d0)*v0)+(2.0d0*al+3.0d0)*
     >             v1)/(3.0d0*(al+2.0d0)*(2.0d0*al+3.0d0))
        rl=x**al
        y1=(x*rl)*(1.0d0+x*(a1+x*(a2+x*a3)))
        y2=rl*(al+x*(alp1*a1+x*((al+2.0d0)*a2+x*(al+3.0d0)*a3)))
      endif
      if(nv.eq.2) then
	y(1)=dreal(y1)
	y(2)=dreal(y2)
      else
	y(1)=dreal(y1)
	y(2)=dimag(y1)
	y(3)=dreal(y2)
	y(4)=dimag(y2)
      endif
c
      return
      end
