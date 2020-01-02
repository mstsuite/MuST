c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine int_zz_zj(mtasa,zj_flag,lmax,kkrsz,
     >                     pnrel,matom_left,matom_right,
     >                     rins,r_sph,r_mesh,jmt,
     >                     zlr_left,zlr_right,jlr,
     &                     nprpts,
     >                     ngaussr,
     >                     cgnt,lmax_cg,
     >                     pzzck,pzz,pzjck,pzj,dzz,dzj,vzz,vzj,
     &                     ncrit,grwylm,gwwylm,wylm,
     >                     iprint,istop)
c     =================================================================
c
      implicit   none
c
!      include   'atom_param.h'
c
      character  sname*32
      character  istop*32
c
      integer    mtasa
      integer    zj_flag
      integer    lmax
      integer    kkrsz
      integer    jmt
      integer    l
      integer    m
      integer    lm
      integer    ir
      integer    lm1
      integer    lm2
      integer    m1
      integer    l2
      integer    m2
      integer    ngaussr
      integer    iprint
      integer    lmax_cg,nprpts
      integer ncrit
c
      real*8     rins,r_sph,r_mesh(nprpts),rtmp(0:nprpts)
      real*8     cgnt(lmax_cg+1,(lmax_cg+1)**2,(lmax_cg+1)**2)
      real*8 pi,fnpi

      real*8 grwylm(*),gwwylm(*)
c
      complex*16 pnrel
      complex*16 matom_left(0:lmax)
      complex*16 matom_right(0:lmax)
      complex*16 pzz(kkrsz,kkrsz)
      complex*16 pzj
      complex*16 pzzck(kkrsz,kkrsz)
      complex*16 pzjck
      complex*16 pzjckout(9)
      complex*16 dzz(kkrsz,kkrsz,-1:1)
      complex*16 dzj
      complex*16 vzz(kkrsz,kkrsz,-1:1)
      complex*16 vzj
      complex*16 zlr_left(nprpts,0:lmax)
      complex*16 zlr_right(nprpts,0:lmax)
      complex*16 jlr(nprpts,0:lmax)
      complex*16 fr(0:nprpts)
      complex*16 fr_int(0:nprpts)
      complex*16 zlzl
      complex*16 zljl
      complex*16 czero
      complex*16 cone
      complex*16 dummy

      complex*16 wylm(*)
c
      parameter (sname='int_zz_zj')
      parameter (czero=(0.0d0,0.0d0))
      parameter (cone=(1.0d0,0.0d0))
c
c     *****************************************************************
c     returns:
c                      Rmt     2            -> ^  ->
c            pzzck = + int dr*r * 4*pi * Z (r)*Z (r)
c                       0                 L     L'
c
c                      Rmt     2            -> ^  ->
c            pzjck = + int dr*r * 4*pi * Z (r)*J (r) * delta
c                       0                 L     L'          LL'
c
c                      Rs   3->     -> ^  ->       ->
c              pzz = + int d r * Z (r)*Z (r)*Sigma(r)
c                       0         L     L'
c
c                      Rs   3->     -> ^  ->       ->
c              pzj = + int d r * Z (r)*J (r)*Sigma(r) * delta
c                       0         L     L'                   LL'
c
c  dzz is contribution to the dipole moment
c                      Rs   3->     ->    ->       ->    *
c              dzz = + int d r * Z (r)*Z (r)*Theta(r) * Y  * r
c                       0         L     L'               1m
c
c  vzz is contribution to the gradient of the dipole potential at the origin
c                      Rs   3->     ->    ->       ->    *         2
c              vzz = + int d r * Z (r)*Z (r)*Theta(r) * Y  / (3 * r )
c                       0         L     L'               1m
c dzz and vzz are calculated only up to r_sph (the real MT radius)
c     *****************************************************************
c
c     -----------------------------------------------------------------
      pi=fnpi()
      call zeroout(pzzck,2*kkrsz*kkrsz)
      call zeroout(pzjckout,2*9)
      pzjck=czero
      rtmp(0)=0.d0
      do ir=1,jmt
	rtmp(ir)=sqrt(r_mesh(ir))
      enddo
c     -----------------------------------------------------------------
      lm=0
      do l=0,lmax
c        =============================================================
c        zlzl: integral of (rad wf)**2 over M.T.
c        =============================================================
         do ir=1,jmt
            fr(ir)=2.d0*zlr_left(ir,l)*zlr_right(ir,l)  !*r_mesh(ir)**2
         enddo
	 fr(0)=fr(1)
c        -------------------------------------------------------------
         call cnewint(jmt+1,rtmp,fr,fr_int,5)
c        -------------------------------------------------------------
c        zlzl=fr_int(jmt)
         call cinterp(rtmp,fr_int,jmt+1,sqrt(r_sph),zlzl,dummy,fr)
c        =============================================================
c        zljl: integral of (rad wf * ir rad wf) over M.T.
c        =============================================================
         if(zj_flag.eq.1) then
            do ir=1,jmt
               fr(ir)=2.d0*zlr_left(ir,l)*jlr(ir,l)  !*r_mesh(ir)**2
            enddo
            fr(0)=fr(1)
c           ----------------------------------------------------------
            call cnewint(jmt+1,rtmp,fr,fr_int,5)
c           ----------------------------------------------------------
c           zljl=fr_int(jmt)
            call cinterp(rtmp,fr_int,jmt+1,sqrt(r_sph),zljl,dummy,fr)
         else
            zljl=czero
         endif
         do m=-l,l
            lm=lm+1
ctest ctmp for l decomposed dos!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c           if(l.eq.1)then
ctest ctmp for l decomposed dos!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            pzzck(lm,lm)=zlzl
            pzjck=pzjck+zljl
c           endif
            if(lm.le.9)pzjckout(lm)=-zljl/pi
         enddo
      enddo
c     =================================================================
c     if(iprint.ge.100)write(6,'(i5,2f14.8,18f9.4,f18.9'' doss'')')
c    >iprint-100,-pzjck/pi,pzjckout,dble(pnrel**2)
c     printout if needed...............................................
      if (iprint.ge.1.and.iprint.lt.100) then
c        -------------------------------------------------------------
         write(6,'(/,'' INT_ZZ_ZJ:: pzz Muffin-tin'')')
c        -------------------------------------------------------------
         call wrtmtx(pzzck,kkrsz,istop)
c        -------------------------------------------------------------
         write(6,'(/,'' INT_ZZ_ZJ:: pzj Muffin-tin'',1p2e18.8)') pzjck
c        -------------------------------------------------------------
      endif

      call zeroout(dzz,2*kkrsz*kkrsz*3)
      call zeroout(vzz,2*kkrsz*kkrsz*3)
      dzj=czero
      vzj=czero
      do l=0,lmax
         do l2=abs(l-1),min(l+1,lmax),2
c        =============================================================
c        zlzl: integral of (rad wf)**2*r over M.T.
c        =============================================================
	 fr(0)=0.d0
         do ir=1,jmt
            fr(ir)=2.d0*zlr_left(ir,l)*zlr_right(ir,l2)*r_mesh(ir)**2
         enddo
c        -------------------------------------------------------------
	 call cnewint(jmt+1,rtmp,fr,fr_int,3)
c        -------------------------------------------------------------
c        zlzl=fr_int(jmt)
         call cinterp(rtmp,fr_int,jmt+1,sqrt(rins),zlzl,dummy,fr)
c        =============================================================
c        zljl: integral of (rad wf)**2/(3*r**2) over M.T.
c        =============================================================
	 fr(0)=0.d0
         do ir=1,jmt
            fr(ir)=2.d0*zlr_left(ir,l)*zlr_right(ir,l2)
         enddo
c        -------------------------------------------------------------
	 call cnewint(jmt+1,rtmp,fr,fr_int,1)
c        -------------------------------------------------------------
c        zljl=fr_int(jmt)/3.d0
         call cinterp(rtmp,fr_int,jmt+1,sqrt(rins),zljl,dummy,fr)
         zljl=zljl/3.d0

	 do m1=-1,1
	  lm1=m1+3
          lm=l*l
         do m=-l,l
            lm=lm+1
c dipole moment is int rho(r)*r*Ylm^*, 
c  cgnt(lm,lmp,lmpp)=int dO Y(lm),Y(lmp),Y(lmpp)^*
c   Note the wavefunction integrals should be Y(lm)Y(lm2)^*
c  so we use the complex conjugate of cgnt
c  cgnt(lm,lmp,lmpp)=int dO Y(lm)^*,Y(lmp)^*,Y(lmpp)
	    m2=m-m1
	    if(abs(m2).le.l2) then
	    lm2=l2*(l2+1)+m2+1
            dzz(lm,lm2,m1)=zlzl*cgnt(l2/2+1,lm1,lm)
            vzz(lm,lm2,m1)=zljl*cgnt(l2/2+1,lm1,lm)
	    endif
         enddo
         enddo  ! m1
         enddo  ! l2
      enddo
c
c     =================================================================
c     For ASA potentials...............................................
c     =================================================================
      if(mtasa.eq.1) then
ctest if(mtasa.ge.1) then
ctest if(mtasa.eq.2)write(6,'(''WARNING******messing with int_zz_zj'')')
c        --------------------------------------------------------------
         call mbeqa(pzzck,pzz,2*kkrsz*kkrsz)
         pzj=pzjck
c        --------------------------------------------------------------
         return
      endif

c Interstitial contributions
      call inter(zj_flag,lmax,
     >           pnrel,matom_left,matom_right,
     >           r_sph,
     >           ngaussr,
     >           cgnt,lmax_cg,
     >           pzz,pzj,
     &           ncrit,grwylm,gwwylm,wylm,
     >           iprint,istop)
c
c     =================================================================
c     add the MT and the interstial contributions......................
c     -----------------------------------------------------------------
      call zaxpy(kkrsz*kkrsz,cone,pzzck,1,pzz,1)
      pzj=pzj+pzjck
c     -----------------------------------------------------------------
c
c     =================================================================
c     printout if needed...............................................
      if (iprint.ge.1.and.iprint.lt.100) then
         write(6,'(/,'' INT_ZZ_ZJ:: pzz total cell'')')
c        -------------------------------------------------------------
         call wrtmtx(pzz,kkrsz,istop)
c        -------------------------------------------------------------
         write(6,'(/,'' INT_ZZ_ZJ:: pzj total cell'',1p2e18.9)')
     &          pzj
c        -------------------------------------------------------------
      endif
c
c     =================================================================
      if (istop.eq.sname) then
        call fstop(sname)
      else
        return
      endif
c
      end
