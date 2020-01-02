c
c     ===================================================================
      function alpha2(rs,dz,sp,iexch,exchg)
c     ===================================================================
c
      implicit none
c
      integer  iexch
      integer  incof
      integer  i
c
      real*8   alpha2
      real*8   rs
      real*8   dz
      real*8   sp
      real*8   exchg
      real*8   ccp
      real*8   rp
      real*8   ccf
      real*8   rf
      real*8   a(3)
      real*8   b(3)
      real*8   c(3)
      real*8   x0(3)
      real*8   cst
      real*8   aip
      real*8   fnot
      real*8   bip
      real*8   for3
      real*8   thrd
      real*8   vxx(2)
      real*8   vxcs(2)
      real*8   g(3)
      real*8   dg(3)
      real*8   tbq(3)
      real*8   tbxq(3)
      real*8   bxx(3)
      real*8   q(3)
      real*8   bb(3)
      real*8   cx(3)
      real*8   fm
      real*8   fdz
      real*8   ex
      real*8   exf
      real*8   xp
      real*8   xf
      real*8   gp
      real*8   gf
      real*8   exc
      real*8   excf
      real*8   dedz
      real*8   gpp
      real*8   gfp
      real*8   depd
      real*8   defd
      real*8   decd
      real*8   bfc
      real*8   zp1
      real*8   zm1
      real*8   xr
      real*8   pex
      real*8   xrsq
      real*8   qi
      real*8   txb
      real*8   fx
      real*8   arct
      real*8   dxs
      real*8   vcc
      real*8   facc
      real*8   ecp
      real*8   zp3
      real*8   zm3
      real*8   zp3m3
      real*8   fx1
      real*8   z4
      real*8   fz
      real*8   beta
      real*8   ec
      real*8   f3ex
c
c     data for von barth-hedin
c
      data     ccp,rp,ccf,rf/0.0450d+00,21.0d+00,
     >                       0.02250d+00,52.9166820d+00/
c
c     data for vosko-wilk-nusair
c
      data     incof/0/
      data     a/-0.0337740d+00,0.06218140d+00,0.03109070d+00/
      data     b/1.131070d+00,3.727440d+00,7.060420d+00/
      data     c/13.00450d+00,12.93520d+00,18.05780d+00/
      data     x0/-0.00475840d+00,-0.104980d+00,-0.32500d+00/
      data     cst,aip/1.923661050d+00,0.916330590d+00/
      data     fnot,bip,for3,thrd/1.709920950d+00,0.259921050d+00,
     >                            1.33333333330d+00,0.33333333330d+00/
c
c     the following are constants needed to obtain potential
c     which are in g.s. painter's paper
c     =====given here for check=====(generated below)
c
c     data q/.7123108918d+01,.6151990820d+01,.4730926910d+01/
c     data bxx/-.4140337943d-03,-.3116760868d-01,-.1446006102/
c     data tbq/.3175776232d+00,.1211783343d+01,.2984793524d+01/
c     data tbxq/.3149055315d+00,.1143525764d+01,.2710005934d+01/
c     data bb/.4526137444d+01,.1534828576d+02,.3194948257d+02/
c
c     write(6,'('' iexch,rs,dz,sp '',i5,3f10.5)') iexch,rs,dz,sp
c
c*********** for positron potential construction iexch.ge.100 *********
      if(iexch.ge.100)then
        alpha2=0.d0
        exchg=0.d0
        return
      endif
c**********************************************************************
      go to (10,20) iexch
c
c     von barth-hedin  exch-corr potential
c     j. phys. c5,1629(1972)
c
c
  10  continue
      fm=2.0d+00**(4.0d+00/3.0d+00)-2.0d+00
      fdz = ((1.0d+00+dz)**(4.0d+00/3.0d+00)
     >     +(1.0d+00-dz)**(4.0d+00/3.0d+00)-2.0d+00)/fm
      ex=-0.916330d+00/rs
      exf=ex*2.0d+00**0.333333330d+00
      xp=rs/rp
      xf=rs/rf
      gp = (1.0d+00+xp**3)*log(1.0d+00+1.0d+00/xp)
     >    -xp*xp +xp/2.0d+00 - 0.333333330d+00
      gf = (1.0d+00+xf**3)*log(1.0d+00+1.0d+00/xf)
     >    -xf*xf +xf/2.0d+00 - 0.333333330d+00
      exc = ex-ccp*gp
      excf=exf-ccf*gf
      dedz= (4.0d+00/3.0d+00)*(excf-exc)
     >     *((1.0d+00+dz)**(1.0d+00/3.0d+00)
     >     -(1.0d+00-dz)**(1.0d+00/3.0d+00))/fm
      gpp = 3.0d+00*xp*xp*log(1.0d+00+1.0d+00/xp)-1.0d+00/xp
     >     +1.50d+00-3.0d+00*xp
      gfp = 3.0d+00*xf*xf*log(1.0d+00+1.0d+00/xf)-1.0d+00/xf
     >     +1.50d+00-3.0d+00*xf
      depd=-ex/rs-ccp/rp*gpp
      defd=-exf/rs-ccf/rf*gfp
      decd=depd+(defd-depd)*fdz
c     exchange-correlation energy
      exchg= exc + (excf-exc)*fdz
c     exchange-correlation potential
      alpha2 = exc+(excf-exc)*fdz-rs*decd/3.0d+00
     >        +sp*(1.0d+00-sp*dz)*dedz
      return
c
c
  20  continue
      if(incof.ne.0) go to 30
      incof=2
c
c     vosko-wilk-nusair exch-corr potential
c     taken from g.s. painter
c     phys. rev. b24 4264,1981
c
c     generate constant coefficients for the parameterization (v-w-n)
c
      do i=1,3
         cx(i)= x0(i)**2 + b(i)*x0(i) + c(i)
         bfc= 4.0d+00*c(i) - b(i)**2.0d+00
         q(i)= sqrt(bfc)
         bxx(i)= b(i)*x0(i)/cx(i)
         tbq(i)= 2.0d+00*b(i)/q(i)
         tbxq(i)= tbq(i) + 4.0d+00*x0(i)/q(i)
         bb(i)= 4.0d+00*b(i)*( 1 - x0(i)*(b(i) + 2.0d+00*x0(i))/cx(i) )
      enddo
c
  30  continue
      zp1= 1.0d+00 + dz
      zm1= 1.0d+00 - dz
      xr=sqrt(rs)
      pex= -aip/rs
      xrsq= rs
c
c     generate g(i)=alpha,epsilon fct.s
c     and their derivatives dg(i)
c     1=alpha(spin stiffness)  2=ecp  3=ecf
c
      do i=1,3
         qi=q(i)
         txb= 2.0d+00*xr + b(i)
         fx= xrsq + xr*b(i) + c(i)
         arct= atan2(qi,txb)
         dxs= (xr-x0(i))**2/fx
         g(i)=a(i)*( log(xrsq/fx) + tbq(i)*arct-bxx(i)*(log(dxs) 
     >              +tbxq(i)*arct) )
         dg(i)=a(i)*( 2.0d+00/xr - txb/fx 
     >               -bxx(i)*(2.0d+00/(xr-x0(i))-txb/fx) 
     >               -bb(i)/(qi**2 + txb**2) )
      enddo
c
      ecp=g(2)
      zp3=zp1**thrd
      zm3=zm1**thrd
      zp3m3=zp3-zm3
c     part of last term in vx   eq(13)
      fx1=.50d+00*for3*pex*zp3m3
      z4= dz**4
      fz= cst*(zp1**for3 + zm1**for3 - 2.0d+00)
      beta= fnot*( g(3)-g(2) )/g(1) -1.0d+00
      ec= ecp + fz*g(1)*( 1.0d+00 + z4*beta )/fnot
      ex= pex*( 1.0d+00 + fz*bip )
      f3ex= for3*ex
c     echange-correlation energy
      exchg= ec + ex
c     exchange potential
      vxx(1)= f3ex + fx1*zm1
      vxx(2)= f3ex - fx1*zp1
c     correlation potential
      vcc= ec - xr*( (1.0d+00-z4*fz)*dg(2) + z4*fz*dg(3)
     >              +(1.0d+00 - z4)*fz*dg(1)/fnot )/6.0d+00
c
      facc= 4.0d+00*g(1)*( dz**3*fz*beta 
     >                    +( 1.0d+00 + beta*z4 )*zp3m3/(6.0d+00*bip) )
     >                  /fnot
c
c     exch-corr. potential for each spin as called in newpot
c
      vxcs(1)= vcc + zm1*facc + vxx(1)
      vxcs(2)= vcc - zp1*facc + vxx(2)
c
      if( sp.ge.0 ) alpha2= vxcs(1)
      if( sp.lt.0 ) alpha2= vxcs(2)
      return
      end
