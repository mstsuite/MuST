      subroutine interstitial(zj_flag,lmax,
     >                     pnrel,matom_left,matom_right,
     >                     r_sph,
     >                     ncrit,ngaussr,
     >                     grwylm,gwwylm,wylm,
     >                     cgnt,lmax_cg,
     >                     pzz,pzj,bj,bh)
c     =================================================================
c
      implicit   none
c
      integer    zj_flag
      integer    lmax,lmax_cg
      integer    ncrit
      integer    ir
      integer    lm1
      integer    lm2
      integer    l1
      integer    m1
      integer    l2
      integer    m2
      integer    l3
      integer    m3
      integer    iwylm
      integer    ng
      integer    ngaussr
c
      real*8     r_sph
      real*8     rgauss
      real*8     cgnt(lmax_cg+1,(lmax_cg+1)**2,(lmax_cg+1)**2)
      real*8     gwwylm(ngaussr,ncrit-1)
      real*8     grwylm(ngaussr,ncrit-1)
      real*8     pi,rt4pi
      parameter (pi=3.141592653589793d0)
c
      complex*16 pnrel
      complex*16 wylm((2*lmax+1)*(lmax+1),ngaussr,ncrit-1)
      complex*16 matom_left(0:lmax)
      complex*16 matom_right(0:lmax)
      complex*16 bj(0:lmax)
      complex*16 bh(0:lmax)
      complex*16 eiz
      complex*16 pzz((lmax+1)**2,(lmax+1)**2)
      complex*16 pzj
      complex*16 zlzl
      complex*16 zljl
      complex*16 zlm1
      complex*16 zlm2
      complex*16 jlm1
      complex*16 steplm
      complex*16 czero
      complex*16 sqrtm1
c
      parameter (czero=(0.0d0,0.0d0))
      parameter (sqrtm1=(0.0d0,1.0d0))
c
c     *****************************************************************
c     returns:
c                      Rs   3->     -> ^  ->       ->
c              pzz = + int d r * Z (r)*Z (r)*Sigma(r)
c                       0         L     L'
c
c                      Rs   3->     -> ^  ->       ->
c              pzj = + int d r * Z (r)*J (r)*Sigma(r) * delta
c                       0         L     L'                   LL'
c
c     *****************************************************************
c
c     =================================================================
c     For muffin-tin potentials and r=>Rmt, one has
c        ->                                     ->
c     Z (r) = -p*[j (p*r)*cot(d ) - n (p*r)]*Y (r)
c      L           l           l     l        L
c
c     ^  ->                                   * ->
c     Z (r) = -p*[j (p*r)*cot(d ) - n (p*r)]*Y (r)
c      L           l           l     l        L
c
c        ->              ->
c     J (r) = j (p*r)*Y (r)
c      L       l       L
c
c     ^  ->            * ->
c     J (r) = j (p*r)*Y (r)
c      L       l       L
c
c     calculate the interstial integrals...............................
c     loop over each region between MT-sphere and circumscibing sphere.
c     =================================================================
c     -----------------------------------------------------------------
      call zeroout(pzz,2*(lmax+1)**4)
      pzj=czero
      rt4pi=sqrt(4.d0*pi)
c     -----------------------------------------------------------------
      do ir=1,ncrit-1
c        ==============================================================
c        loop over the gaussian point in each region...................
c        ==============================================================
         do ng=1,ngaussr
            rgauss=grwylm(ng,ir)
	    if(rgauss.gt.0.d0) then
c           -----------------------------------------------------------
            call zsphbesjh(lmax,pnrel*rgauss,bj,bh,eiz,1)
c           -----------------------------------------------------------
            eiz=sqrtm1*eiz/rgauss
c           ===========================================================
c           loop over l1,m1............................................
c           ===========================================================
            lm1=0
	    do l1=0,lmax
               zlm1=matom_left(l1)*bj(l1)-eiz*bh(l1)
               jlm1=eiz*bh(l1)/matom_left(l1)
	    do m1=-l1,l1
	      lm1=lm1+1
	       zljl=0.d0
c              ========================================================
c              loop over l2,m2.........................................
	       lm2=0
	       do l2=0,lmax
                  zlm2=matom_right(l2)*bj(l2)-eiz*bh(l2)
	       do m2=-l2,l2
		 lm2=lm2+1
c                 =====================================================
c                 perform sum over l3,m3 with gaunt # and gaussian wt..
		  zlzl=0.d0
                     m3=m1-m2
                  do l3=l1+l2,max(abs(m3),abs(l1-l2)),-2
                     iwylm= (l3*(l3+1))/2 + abs(m3) +1
c theta(r) is formed with sum_lm3 steplm*Ylm
c  cgnt(lm,lmp,lmpp)=int dO Y(lm),Y(lmp),Y(lmpp)^*
                     if(m3.ge.0) then
                        steplm=wylm(iwylm,ng,ir)
                     else
                        steplm=(1-2*mod(abs(m3),2))*
     &                     conjg(wylm(iwylm,ng,ir))
                     endif
		     if(l3.eq.0) then
		     if(abs(rgauss-r_sph).lt.1.d-10) then
		       steplm=steplm-0.5d0*rt4pi
		     else if(rgauss.lt.r_sph) then
		       steplm=steplm-rt4pi
		     endif
		     endif
                     zlzl=zlzl+cgnt(l3/2+1,lm2,lm1)*steplm
                     if(lm2.eq.lm1 .and. zj_flag.eq.1) then
                        zljl=zljl+cgnt(l3/2+1,lm2,lm1)*steplm
                     endif
                  enddo
                  pzz(lm2,lm1)=pzz(lm2,lm1)+zlzl*zlm1*zlm2*gwwylm(ng,ir)
               enddo
               enddo
               pzj=pzj+zljl*zlm1*jlm1*gwwylm(ng,ir)
            enddo
            enddo
	    endif  !  rgauss.gt.zero
         enddo
      enddo
c
      return
      end
