c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scalarr(nrelv,clight,l,
     >                  bjl,bnl,bj,bn,djl,dnl,g,f,gpl,fpl,
     >                  tandel,matom,
     >                  energy,prel,pnrel,
     >                  rv,r,h,jmt,iend,icmax,
     >                  r_sph,iprint,istop)
c     ================================================================
c
      implicit   complex*16 (a-h,o-z)
c
      character  istop*32
      character  sname*20
c
      integer    nrelv,i_out
      integer    l
      integer    jmt
      integer    iend
      integer    iprint
c
      real*8     r(iend),r_sph,rv_sph,drv_sph
      real*8     rv(iend)
      real*8     zed
      real*8     h,h24
      real*8     rvr,drvr
      real*8     v0,v1
      real*8     cinv
      real*8     clight
      real*8     c2inv
      real*8     power
c     real*8     ylag
      real*8     tole
      real*8     tzoc
      real*8     fz2oc2
      real*8     tzoc2
      real*8     one
c
      complex*16 f_sph,g_sph,df_sph,dg_sph
      complex*16 energy
      complex*16 prel
      complex*16 pnrel
      complex*16 g(iend)
      complex*16 f(iend)
      complex*16 gpl(iend)
      complex*16 fpl(iend)
      complex*16 cl
      complex*16 sl
      complex*16 dcl(4)
      complex*16 dsl(4)
      complex*16 ddcl,ddsl,dgpl,dfpl
      complex*16 bjl(2,iend)
      complex*16 bnl(2,iend)
      complex*16 bj(l+1)
      complex*16 bn(l+1)
      complex*16 djl(l+1)
      complex*16 dnl(l+1)
      complex*16 a1,a2,b1,b2
      complex*16 matom,sqrtm1
c
      parameter (tole=1.d-10)
      parameter (one=1.d0)
      parameter (sqrtm1=(0.d0,1.d0))
      parameter (sname='scalar')
c
      h24=h/24.d0
c
c     ******************************************************************
c     this subroutine sets up to solve the scalar relativistic equation.
c     ******************************************************************
c
      if(iprint.ge.1) then
         write(6,'(/'' SCALAR: nrelv,clight='',i5,d23.14)') nrelv,clight
         write(6,'('' SCALAR: energy='',2d23.14)') energy
         write(6,'('' SCALAR:   prel='',2d23.14)') prel
         write(6,'('' SCALAR: l,jmt,iend ='',4i5)')
     >   l,jmt,iend
         write(6,'('' SCALAR: r,rv at jmt ='',2d16.8)')r(jmt),rv(jmt)
      endif
c     =================================================================
c     cannot solve for the real part of the energy equal to 0.0!!!!!!!!
      if( abs(energy) .lt. tole ) then
        write(6,'('' scalar:: energy=zero'')')
        call fstop(sname)
      endif
c     =================================================================
c     set atomic number................................................
      ized=-rv(1)/2 + 0.5
      zed=ized
c     =================================================================
c     start with small r expansion of wave-fcts
      cinv=one/clight
      c2inv=cinv*cinv
      eoc2p1=one+energy*c2inv
      tzoc=2*zed*cinv
      fz2oc2=tzoc*tzoc
      tzoc2=tzoc*cinv
c     v0=(rv(1)+2*zed)/r(1)
      v1=(r(2)*(rv(1)+2*zed)-r(1)*(rv(2)+2*zed))/(r(1)*r(2)*(r(1)-r(2)))
      v0=v1*(r(1)+r(2))+(rv(1)-rv(2))/(r(1)-r(2))
      em0= 1 + (energy-v0)*c2inv
      v0me= v0-energy
      if(nrelv.le.0) then
c        ==============================================================
c        scalar-relativistic case......................................
         power=sqrt( l*(l+1) + 1 - fz2oc2 )
         a0=1.d0
         a1p=( fz2oc2 + em0*(power-1 -2*fz2oc2) )/(2*power+1)
         a1 = a0*a1p/tzoc2
         a2p=( a1p*(fz2oc2 + em0*(3*power+1 - 2*fz2oc2)) + 
     <         em0*em0*(fz2oc2 - 2*power+2) )/(4*power+4)
         a2 = a0*a2p/(tzoc2*tzoc2)
         a3p=( a2p*( fz2oc2 + em0*(5*power+5 - 2*fz2oc2) )
     <        - a1p*(4*power+1 - fz2oc2)*em0*em0
     <        + (3*power-3 - fz2oc2)*em0*em0*em0  )/(6*power+9)
         a3 = a0*a3p/(tzoc2*tzoc2*tzoc2)
      else
c        ==============================================================
c        non-relativistic case......................................
         power= l+1
c        a0 is arbitrary, but this make Rl a spherical Bessel at origin.
c        other than a factor of 1/(2l+1)!!
         a0= prel*(l+1)
         a1p= -2*zed/(2*l+2)
         a1 = a0*a1p
         a2p= (v0me - a1p*2*zed)/(4*l+6)
         a2 = a0*a2p
C        a3p= (a1p*v0me - a2p*2*zed+2*v1)/(6*l+12)
         a3p=(a1p*(4*zed*zed+(3*l+4)*v0me)/(2*l+3)+v1)/(6*l+12)
         a3 = a0*a3p
      endif
c     =================================================================
c     get first 4 points from power series expansion...................
      do j=1,4
c        ==============================================================
c        r*G = g and r*F = f, i.e. this is r*wave-fct.
         g(j)= r(j)**power*(a0 + r(j)*(a1 + r(j)*(a2 + a3*r(j))) )
         emr= eoc2p1*r(j) - rv(j)*c2inv
         f(j)= r(j)**power*( a0*(power-1) + 
     <   r(j)*(a1*power + r(j)*(a2*(power+1) + a3*(power+2)*r(j))))/emr
c        ==============================================================
c        get cl's and sl's
         cl=-bnl(2,j)*g(j) - bnl(1,j)*f(j)/prel
         sl=-bjl(2,j)*g(j) - bjl(1,j)*f(j)/prel
c        ==============================================================
c        get derivatives of cl and sl :: for predictor-corrector use 
c        constant increments on log grid therefore need dg/dx = r*dg/dr
         em= 1 + (energy - rv(j)/r(j))*c2inv
	 b1=(em-1)*f(j)*r(j)
	 b2=( l*(l+1)/(em*r(j)) + rv(j)  + 
     >           (eoc2p1-1)*energy*r(j) )*g(j)/prel
         dcl(j)= -bnl(2,j)*b1 - bnl(1,j)*b2
         dsl(j)= -bjl(2,j)*b1 - bjl(1,j)*b2
      enddo
c
c     =================================================================
c     regular solution and phase-shifts of scalar relativistic eqs.
      j1=4
      j2=3
      j3=2
      j4=1
      i_out=0
      do j=5,jmt
           if(i_out.eq.0.and.r(j).ge.r_sph)i_out=j
           clold=cl
           slold=sl
c        ==============================================================
c        evaluate predictor
         cl=clold+h24*( 55.0d+00*dcl(j1)-59.0d+00*dcl(j2)+
     >                     37.0d+00*dcl(j3)- 9.0d+00*dcl(j4) )
         sl=slold+h24*( 55.0d+00*dsl(j1)-59.0d+00*dsl(j2)+
     >                     37.0d+00*dsl(j3)- 9.0d+00*dsl(j4) )
c        ==============================================================
c        evaluate corrector
         em = eoc2p1 - c2inv*rv(j)/r(j)
         emr= em*r(j)
         do icor=1,icmax
           y1=cl*bjl(1,j)-sl*bnl(1,j)
           y2=cl*bjl(2,j)-sl*bnl(2,j)
	   b1=prel*(em-1)*y2*r(j)
	   b2=( l*(l+1)/emr + rv(j)  + (eoc2p1-1)*energy*r(j) )*y1/prel
           dcl(j4)= bnl(2,j)*b1 - bnl(1,j)*b2
           dsl(j4)= bjl(2,j)*b1 - bjl(1,j)*b2
           cl=clold+h24*( 9.0d+00*dcl(j4)+19.0d+00*dcl(j1)
     >                    - 5.0d+00*dcl(j2)+dcl(j3) )
           sl=slold+h24*( 9.0d+00*dsl(j4)+19.0d+00*dsl(j1)
     >                    - 5.0d+00*dsl(j2)+dsl(j3) )
        enddo
           y1=cl*bjl(1,j)-sl*bnl(1,j)
           y2=cl*bjl(2,j)-sl*bnl(2,j)
           g(j)=y1
           f(j)=-prel*y2
	   j0=j4
	   j4=j3
	   j3=j2
	   j2=j1
	   j1=j0
      enddo
c     =================================================================
c     get normalization etc. at sphere boundary........................
c     At R, need spherical Bessel's for all l to get normalization and 
c     physical phase-shifts to determine cl and sl from 
c     g=r*R=r*( cl*jl - sl*nl ) and f=r*R'/M for all l's. 
c     Modify expression to account for difference between spherical and 
c     Ricatti-Bessel fcts.               Recall: R=g/r and dR=f*em/r
      j=jmt
c     pr=prel*r(jmt)
      pr=prel*r_sph
c     em = eoc2p1 - c2inv*rv(j)/r(j)
      call interp(r,rv,i_out,r_sph,rv_sph,drv_sph,.false.)
      call cinterp(r,f,i_out,r_sph,f_sph,df_sph,gpl)
      call cinterp(r,g,i_out,r_sph,g_sph,dg_sph,gpl)
c     if(iprint.ge.1)write(6,'(''g_sph'',19d14.7)')r_sph,
c    >r(j-1),r(j),r(j+1),g(j-1),g(j),g(j+1),g_sph
c     em = eoc2p1 - c2inv*rv(j)/r_sph
      em = eoc2p1 - c2inv*rv_sph/r_sph
      emr= em*r_sph
c     -----------------------------------------------------------------
        call ricbes(l+1,pr,bj,bn,djl,dnl)
c     -----------------------------------------------------------------
c     =================================================================
c     sl and cl at the sphere radius
      slmt= (prel*djl(l+1)-bj(l+1)/r_sph)*g_sph - bj(l+1)*f_sph*em
      clmt= (prel*dnl(l+1)-bn(l+1)/r_sph)*g_sph - bn(l+1)*f_sph*em
c     slmt= (prel*djl(l+1)-bj(l+1)/r(j))*g(j) - bj(l+1)*f(j)*em
c     clmt= (prel*dnl(l+1)-bn(l+1)/r(j))*g(j) - bn(l+1)*f(j)*em
c     =================================================================
         gpl(j)=bn(l+1)/prel
         fpl(j)=(dnl(l+1) - bn(l+1)/pr)/eoc2p1
c     phase shifts and normalization
      cotdel=sqrtm1*(slmt-clmt)/slmt
      tandel=one/cotdel
      matom=pnrel*(sqrtm1-cotdel)
      anorm=(bj(l+1)*matom/pnrel-sqrtm1*bn(l+1))/g_sph
      if(iprint.ge.1) then
         write(6,'(/'' SCALAR::     l='',i3)') l
         write(6,'('' SCALAR::   prel='',2d23.14)') prel
         write(6,'('' SCALAR::   clmt='',2d23.14)') clmt
         write(6,'('' SCALAR::   slmt='',2d23.14)') slmt
         write(6,'('' SCALAR:: cotdel='',2d23.14)') cotdel
         write(6,'('' SCALAR:: matoml='',2d23.14)') matom
         write(6,'('' SCALAR::  anorm='',2d23.14)') anorm
         write(6,'('' SCALAR::  pnrel='',2d23.14)') pnrel
         write(6,'('' SCALAR::bj(l+1)='',2d23.14)') bj(l+1)
         write(6,'('' SCALAR::bn(l+1)='',2d23.14)') bn(l+1)
         write(6,'('' SCALAR:: g(j  )='',2d23.14)')  g(j)
         write(6,'('' SCALAR:: g_sph ='',2d23.14)')  g_sph
         write(6,'(/)')
         do j=1,jmt,20
            write(6,'('' j,f,g ='',1i5,4d15.8)')j,f(j),g(j)
         enddo
      endif
      do j=1,jmt
	g(j)=g(j)*anorm
	f(j)=f(j)*anorm
      enddo
c     =================================================================
c     get wave-fcts. beyond mt sphere
      jmtp1=jmt+1
      do j=jmtp1,iend
         pr=prel*r(j)
c        --------------------------------------------------------------
         call ricbes(l+1,pr,bj,bn,djl,dnl)
c        --------------------------------------------------------------
         g(j)=(bj(l+1)*matom/pnrel-sqrtm1*bn(l+1))
         f(j)= ( (djl(l+1)-bj(l+1)/pr)*matom/pnrel
     >         -sqrtm1*(dnl(l+1)-bn(l+1)/pr) ) /eoc2p1
         gpl(j)=bn(l+1)/prel
         fpl(j)=(dnl(l+1) - bn(l+1)/pr)/eoc2p1
      enddo
c     ===============================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scalari(nrelv,clight,l,
     >                  bjl,bnl,bj,bn,djl,dnl,g,f,gpl,fpl,
     >                  tandel,matom,
     >                  energy,prel,pnrel,
     >                  rv,r,h,jmt,iend,icmax,
     >                  r_sph,iprint,istop)
c     ================================================================
c
      implicit   complex*16 (a-h,o-z)
c
      character  istop*32
      character  sname*20
c
      integer    nrelv,i_out
      integer    l,iex
      integer    jmt
      integer    iend
      integer    iprint
c
      real*8     r(iend),r_sph,rv_sph,drv_sph
      real*8     rv(iend)
      real*8     zed
      real*8     h,h24
      real*8     rvr,drvr
      real*8     v0,v1
      real*8     cinv
      real*8     clight
      real*8     c2inv
      real*8     power
      real*8     ylag
      real*8     tole
      real*8     tzoc
      real*8     fz2oc2
      real*8     tzoc2
      real*8     one
c
      complex*16 f_sph,g_sph,df_sph,dg_sph
      complex*16 fpl_sph,gpl_sph
      complex*16 energy
      complex*16 prel
      complex*16 pnrel
      complex*16 g(iend)
      complex*16 f(iend)
      complex*16 gpl(iend)
      complex*16 fpl(iend)
      complex*16 cl
      complex*16 sl
      complex*16 dcl(4)
      complex*16 dsl(4)
      complex*16 ddcl,ddsl,dgpl,dfpl
      complex*16 bjl(2,iend)
      complex*16 bnl(2,iend)
      complex*16 bj(l+1)
      complex*16 bn(l+1)
      complex*16 djl(l+1)
      complex*16 dnl(l+1)
      complex*16 a1,a2,b1,b2
      complex*16 matom,sqrtm1
c
      parameter (tole=1.d-10)
      parameter (one=1.d0)
      parameter (sqrtm1=(0.d0,1.d0))
      parameter (sname='scalar')
c
      h24=h/24.d0
c
c     ******************************************************************
c     this subroutine sets up to solve the scalar relativistic equation.
c     ******************************************************************
c
      if(iprint.ge.1) then
         write(6,'(/'' SCALAR: nrelv,clight='',i5,d23.14)') nrelv,clight
         write(6,'('' SCALAR: energy='',2d23.14)') energy
         write(6,'('' SCALAR:   prel='',2d23.14)') prel
         write(6,'('' SCALAR: l,jmt,iend ='',4i5)')
     >   l,jmt,iend
      endif
c     =================================================================
c     cannot solve for the real part of the energy equal to 0.0!!!!!!!!
      if( abs(energy) .lt. tole ) then
        write(6,'('' scalar:: energy=zero'')')
        call fstop(sname)
      endif
c     =================================================================
c     set atomic number................................................
      ized=-rv(1)/2 + 0.5
      zed=ized
c     =================================================================
c     start with small r expansion of wave-fcts
      cinv=one/clight
      c2inv=cinv*cinv
      eoc2p1=one+energy*c2inv
      tzoc=2*zed*cinv
      fz2oc2=tzoc*tzoc
      tzoc2=tzoc*cinv
c     v0=(rv(1)+2*zed)/r(1)
      v1=(r(2)*(rv(1)+2*zed)-r(1)*(rv(2)+2*zed))/(r(1)*r(2)*(r(1)-r(2)))
      v0=v1*(r(1)+r(2))+(rv(1)-rv(2))/(r(1)-r(2))
      em0= 1 + (energy-v0)*c2inv
      v0me= v0-energy
      if(nrelv.le.0) then
c        ==============================================================
c        scalar-relativistic case......................................
         power=sqrt( l*(l+1) + 1 - fz2oc2 )
         a0=1.d0
         a1p=( fz2oc2 + em0*(power-1 -2*fz2oc2) )/(2*power+1)
         a1 = a0*a1p/tzoc2
         a2p=( a1p*(fz2oc2 + em0*(3*power+1 - 2*fz2oc2)) + 
     <         em0*em0*(fz2oc2 - 2*power+2) )/(4*power+4)
         a2 = a0*a2p/(tzoc2*tzoc2)
         a3p=( a2p*( fz2oc2 + em0*(5*power+5 - 2*fz2oc2) )
     <        - a1p*(4*power+1 - fz2oc2)*em0*em0
     <        + (3*power-3 - fz2oc2)*em0*em0*em0  )/(6*power+9)
         a3 = a0*a3p/(tzoc2*tzoc2*tzoc2)
      else
c        ==============================================================
c        non-relativistic case......................................
         power= l+1
c        a0 is arbitrary, but this make Rl a spherical Bessel at origin.
c        other than a factor of 1/(2l+1)!!
         a0= prel*(l+1)
         a1p= -2*zed/(2*l+2)
         a1 = a0*a1p
         a2p= (v0me - a1p*2*zed)/(4*l+6)
         a2 = a0*a2p
C        a3p= (a1p*v0me - a2p*2*zed+2*v1)/(6*l+12)
         a3p=(a1p*(4*zed*zed+(3*l+4)*v0me)/(2*l+3)+v1)/(6*l+12)
         a3 = a0*a3p
      endif
c     =================================================================
c     get wave-fcts. beyond mt sphere
      jmtp1=jmt+1
      do j=1,iend
c     do j=jmtp1,iend
         pr=prel*r(j)
c        --------------------------------------------------------------
         call ricbes(l+1,pr,bj,bn,djl,dnl)
c        --------------------------------------------------------------
         gpl(j)=bn(l+1)/prel
         fpl(j)=(dnl(l+1) - bn(l+1)/pr)/eoc2p1
      enddo
c
c     =================================================================
c     irregular solution of scalar relativistic eqs....................
c     start irregular solution outside m.t. sphere
c     =================================================================
c     do j=jmt,jmt+3
      i_out=0
      do j=jmt,1,-1
        if(r(j).lt.r_sph.and.i_out.eq.0)i_out=j+1
      enddo
c     if(j.eq.i_out) then
         pr=prel*r_sph
c        --------------------------------------------------------------
         call ricbes(l+1,pr,bj,bn,djl,dnl)
c        --------------------------------------------------------------
         gpl_sph=bn(l+1)/prel
         fpl_sph=(dnl(l+1) - bn(l+1)/pr)/eoc2p1
c     cl=-(-sin(pr)-cos(pr))*gpl_sph +cos(pr) *fpl_sph/prel
c     sl=-(-cos(pr)+sin(pr))*gpl_sph - sin(pr)*fpl_sph/prel
      cl=-bnl(2,i_out)*gpl(i_out) - bnl(1,i_out)*fpl(i_out)/prel
      sl=-bjl(2,i_out)*gpl(i_out) - bjl(1,i_out)*fpl(i_out)/prel
      ddcl=0.d0
      ddsl=0.d0
      dgpl=0.d0
      dfpl=0.d0
      rvr=ylag(r_sph,r,rv,0,2,i_out,iex)
      em = eoc2p1 - c2inv*rvr/r_sph
      emr= em*r_sph
      b1=(em-1)*fpl_sph*r_sph
      b2=( l*(l+1)/emr + rvr  + (eoc2p1-1)*energy*r_sph )
     &*(gpl_sph+dgpl)/prel
       dcl(1)=-(-sin(pr)-cos(pr))*b1 +cos(pr) *b2
       dsl(1)=-(-cos(pr)+sin(pr))*b1 - sin(pr)*b2
       b1=(em-1)*dfpl*r_sph
       b2=b2-( l*(l+1)/emr + (eoc2p1-1)*energy*r_sph )*gpl_sph/prel
       ddcl=-(-sin(pr)-cos(pr))*b1 +cos(pr) *b2
       ddsl=-(-cos(pr)+sin(pr))*b1 - sin(pr)*b2
c     endif
      do j=i_out,i_out+3
c           ===========================================================
c        starting with free solution ( v(r)=0 r > R), but using 
c        extrapolated potentials at and beyond R to get derivative 
c        terms of cl and sl.......................................
         pr=prel*r(j)
c           ===========================================================
c        set cl and sl beyond mt sphere radius of potential
c        extrapolate potential at m.t. radius to get starting value 
c        for irreg.
c	 if(j.eq.506) then
c        if(r(j)=r_sph)then
         dgpl=dgpl+(ddcl*bjl(1,j)-ddsl*bnl(1,j))*h
         dfpl=dfpl-prel*(ddcl*bjl(2,j)-ddsl*bnl(2,j))*h
c	 call interp(r,rv,506,r(j),rvr,drvr,.false.)
c	 call interp(r,rv,i_out,r(j),rvr,drvr,.false.)
         rvr=ylag(r(j),r,rv,0,2,i_out,iex)
         em = eoc2p1 - c2inv*rvr/r(j)
         emr= em*r(j)
c        ===========================================================
c        Note: effect of potential on fpl here ignored: it is small.
	 b1=(em-1)*(fpl(j)+dfpl)*r(j)
	 b2=( l*(l+1)/emr + rvr  + (eoc2p1-1)*energy*r(j) )
     &      *(gpl(j)+dgpl)/prel
         dcl(j-i_out+1)= -bnl(2,j)*b1 - bnl(1,j)*b2
         dsl(j-i_out+1)= -bjl(2,j)*b1 - bjl(1,j)*b2
	 b1=(em-1)*dfpl*r(j)
	 b2=b2-( l*(l+1)/emr + (eoc2p1-1)*energy*r(j) )*gpl(j)/prel
	 ddcl=-bnl(2,j)*b1-bnl(1,j)*b2
	 ddsl=-bjl(2,j)*b1-bjl(1,j)*b2
c        write(*,'(''rvr'',3i5,1f7.4,99d11.5)')j,i_out,j-i_out+1,
c    >   r(j),
c    >   rv(i_out-3),rv(i_out-2),rv(i_out-1),rv(i_out),
c    >   rvr,dsl(j-i_out+1),dcl(j-i_out+1)
c       endif
      enddo
c     ==================================================================
c     get irregular solution inside m.t. sphere
c     ==================================================================
      j1=1
      j2=2
      j3=3
      j4=4
c     jmtm1=506-1
      jmtm1=jmt-1
c     do j=jmtm1,1,-1
c     write(6,'(''cl,sl,clo,slo'',2i4,99d10.4)')
c    >0,0,0.,0.,0.,dsl(j1),dsl(j2),dsl(j3),dsl(j4)
      do j=i_out-1,1,-1
            clold=cl
            slold=sl
c       ================================================================
c        evaluate predictor
         cl=clold-h24*( 55.0d+00*dcl(j1)-59.0d+00*dcl(j2)+
     >           37.0d+00*dcl(j3)- 9.0d+00*dcl(j4) )
         sl=slold+h24*( 55.0d+00*dsl(j1)-59.0d+00*dsl(j2)+
     >           37.0d+00*dsl(j3)- 9.0d+00*dsl(j4) )
c       ================================================================
c        evaluate corrector
         em = eoc2p1 - c2inv*rv(j)/r(j)
         emr= em*r(j)
c     if(j.gt.502)write(6,'(''cl,sl,clo,slo'',2i4,99d10.4)')
c    >j,0,cl,sl,rv(j),dsl(j1),dsl(j2),dsl(j3),dsl(j4)
         do icor=1,icmax
            y1=cl*bjl(1,j)-sl*bnl(1,j)
            y2=cl*bjl(2,j)-sl*bnl(2,j)
	    b1=prel*(em-1)*y2*r(j)
	    b2=( l*(l+1)/emr + rv(j)  + (eoc2p1-1)*energy*r(j) )*y1/prel
            dcl(j4)= bnl(2,j)*b1 - bnl(1,j)*b2
            dsl(j4)= bjl(2,j)*b1 - bjl(1,j)*b2
            cl=clold-h24*( 9.0d+00*dcl(j4)+19.0d+00*dcl(j1)
     >            - 5.0d+00*dcl(j2)+dcl(j3) )
            sl=slold-h24*( 9.0d+00*dsl(j4)+19.0d+00*dsl(j1)
     >            - 5.0d+00*dsl(j2)+dsl(j3) )
c     if(j.gt.502)write(6,'(''cl,sl,clo,slo'',2i4,99d10.4)')
c    >j,icor,cl,sl,clold,slold
         enddo
         y1=cl*bjl(1,j)-sl*bnl(1,j)
         y2=cl*bjl(2,j)-sl*bnl(2,j)
         gpl(j)=y1
         fpl(j)=-prel*y2
	 j0=j4
	 j4=j3
	 j3=j2
	 j2=j1
	 j1=j0
      enddo
c
c     ===============================================================
c     check Wronskian for specified l's.............................. 
c     ===============================================================
c     if(iprint.ge.0 .and. l.eq.0) then
      if(iprint.ge.0) then
         do j=iend,iend-15,-1
            w= -( g(j)*fpl(j)  - gpl(j)*f(j) )
c           write(6,'('' j,r,w ='',1i5,7d23.14)')j,r(j),w,fpl(j),gpl(j)
c           write(6,'('' j,r,w ='',1i5,7d16.8)')j,r(j),fpl(j),gpl(j)
         enddo
c        write(6,'(/)')
c        do j=1,iend,20
c           write(6,'(''ij,g,f ='',1i5,4d23.14)')j,g(j),f(j)
c        enddo
      endif

c
c     ===============================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
