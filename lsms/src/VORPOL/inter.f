c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inter(zj_flag,lmax,
     >                     pnrel,matom_left,matom_right,
     >                     r_sph,
     >                     ngaussr,
     >                     cgnt,lmax_cg,
     >                     pzz,pzj,
     &                     ncrit,grwylm,gwwylm,wylm,
     >                     iprint,istop)
c     =================================================================
c
      implicit   none
c
!      include   'vorpol.h'
c
      character  sname*32
      character  istop*32
c
      integer    zj_flag
      integer    lmax,lmax_cg
      integer    ngaussr
      integer    iprint
      integer ncrit
c
      real*8     r_sph
      real*8     cgnt(lmax_cg+1,(lmax_cg+1)**2,(lmax_cg+1)**2)

      real*8 grwylm(*),gwwylm(*)
c
      complex*16 pnrel
      complex*16 matom_left(0:lmax)
      complex*16 matom_right(0:lmax)
      complex*16 bj(0:lmax)
      complex*16 bh(0:lmax)
      complex*16 pzz((lmax+1)**2,(lmax+1)**2)
      complex*16 pzj

      complex*16 wylm(*)
c
      parameter (sname='inter')
c
c     *****************************************************************
c     returns:
c
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
c     -----------------------------------------------------------------

      call interstitial(zj_flag,lmax,
     >                     pnrel,matom_left,matom_right,
     >                     r_sph,
     >                     ncrit,ngaussr,
     >                     grwylm,gwwylm,wylm,
     >                     cgnt,lmax_cg,
     >                     pzz,pzj,bj,bh)
c
c     =================================================================
c     printout if needed...............................................
      if (iprint.ge.1.and.iprint.lt.100) then
         write(6,'(/,'' INTER:: pzz interstial only'')')
c        -------------------------------------------------------------
         call wrtmtx(pzz,(lmax+1)**2,istop)
c        -------------------------------------------------------------
         write(6,'(/,'' INTER:: pzj interstial only'',1p2e18.9)')
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
