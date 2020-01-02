c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gf_local(mtasa,zj_flag,lmax,kkrsz,
     >                    rins,r_sph,r_mesh,jmt,jws,
     >                    pnrel,tau00,matom_left,matom_right,
     >                    zlr_left,zlr_right,jlr,nprpts,
     >                    ngaussr,
     >                    cgnt,lmax_cg,
     >                    dos,dosck,green,dipole,
     &                    ncrit,grwylm,gwwylm,wylm,
     >                    pi,iprint,istop)
c     ================================================================
c
c
c     ****************************************************************
c     input:
c                tau00   
c                kkrsz   (size of KKR-matrix)
c                istop   (index of subroutine prog. stops in)
c     output:
c                dos     (wigner-seitz cell density of states)
c                dosck   (muffin-tin density of states)
c                green   (Green's function)
c
c                They are all calculated in the Local frame.
c     ****************************************************************
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      include    'atom_param.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  istop*32
      character  sname*32
c
      integer    mtasa
      integer    zj_flag
      integer    lmax
      integer    kkrsz
      integer    jmt
      integer    jws
      integer    m
      integer    ngaussr
      integer    iprint,iprint_dos
      integer lmax_cg,nprpts,ncrit
c
      real*8     rins,r_sph,r_mesh(nprpts)
      real*8     cgnt(lmax_cg+1,(lmax_cg+1)**2,(lmax_cg+1)**2)
      real*8     pi
      real*8     sqr2
      real*8 grwylm(*),gwwylm(*)
c
      complex*16 pnrel
      complex*16 matom_left(lmax+1)
      complex*16 matom_right(lmax+1)
      complex*16 tau00(kkrsz*kkrsz)
      complex*16 pzz(kkrsz*kkrsz)
      complex*16 pzj
      complex*16 pzzck(kkrsz*kkrsz)
      complex*16 pzjck
      complex*16 dzz(kkrsz*kkrsz*3)
      complex*16 dzj
      complex*16 vzz(kkrsz*kkrsz*3)
      complex*16 vzj
      complex*16 zlr_left(nprpts,0:lmax)
      complex*16 zlr_right(nprpts,0:lmax)
      complex*16 jlr(nprpts,0:lmax)
      complex*16 dos
      complex*16 dosck
      complex*16 green(jws)
c dipole(m,1) is the density moment
c dipole(m,2) is the gradient of the dipole potential at the origin
      complex*16 dipole(-1:1,2)
      complex*16 wylm(*)
c
      parameter (sname='gf_local')
c
      sqr2=sqrt(2.d0)
c     =================================================================
c     call int_zz_zj to calculate:
c
c                      Rmt     2            -> ^  ->
c     pzzck = - 1/pi * int dr*r * 4*pi * Z (r)*Z (r)
c                       0                 L     L'
c
c                      Rmt     2            -> ^  ->
c     pzjck = + 1/pi * int dr*r * 4*pi * Z (r)*J (r) * delta
c                       0                 L     L'          LL'
c
c                      Rs   3->     -> ^  ->       ->
c       pzz = - 1/pi * int d r * Z (r)*Z (r)*Sigma(r)
c                       0         L     L'
c
c                      Rs   3->     -> ^  ->       ->
c       pzj = + 1/pi * int d r * Z (r)*J (r)*Sigma(r) * delta
c                       0         L     L'                   LL'
c     =================================================================
c     -----------------------------------------------------------------
      iprint_dos=iprint
      
      call int_zz_zj(mtasa,zj_flag,lmax,kkrsz,
     >               pnrel,matom_left,matom_right,
     >               rins,r_sph,r_mesh,jmt,
c    >               r_sph,r_mesh,506,
     >               zlr_left,zlr_right,jlr,nprpts,
     >               ngaussr,
     >               cgnt,lmax_cg,
     >               pzzck,pzz,pzjck,pzj,dzz,dzj,vzz,vzj,
     &               ncrit,grwylm,gwwylm,wylm,
     >               iprint_dos,istop)
c     -----------------------------------------------------------------
c
c     ================================================================
c     dos =>  ZZ(e)*tau(e) - ZJ(e)....................................
c     ----------------------------------------------------------------
      call mdosms(zj_flag,kkrsz*kkrsz,dosck,pzzck,pzjck,tau00,pi,
     >            iprint,istop)
c     ----------------------------------------------------------------
      if(iprint.ge.1 .and. zj_flag.eq.1) then
         write(6,'('' e,n(e)_mt:'',2d16.8,2f18.13)') pnrel*pnrel,dosck
      endif
c     ----------------------------------------------------------------
      call mdosms(zj_flag,kkrsz*kkrsz,dos,pzz,pzj,tau00,pi,
     >            iprint_dos,istop)
c     ----------------------------------------------------------------
      if(iprint.ge.1 .and. zj_flag.eq.1) then
         write(6,'('' e,n(e)_ws:'',2d16.8,2f18.13)') pnrel*pnrel,dos
      endif
c     ----------------------------------------------------------------
      do m=-1,1
         call mdosms(zj_flag,kkrsz*kkrsz,dipole(m,1),
     >   dzz(kkrsz*kkrsz*(m+1)+1),dzj,tau00,pi,iprint,istop)
         call mdosms(zj_flag,kkrsz*kkrsz,dipole(m,2),
     >   vzz(kkrsz*kkrsz*(m+1)+1),vzj,tau00,pi,iprint,istop)
      enddo
c change to real harmonic basis so we can take the imag part
c  dipole is Ylm^* which is
c        (x+iy)/r2
c        (  z  )
c        (-x+iy)/r2
c we want to convert it into real harmonics as
c        (y)
c        (z)
c        (x)
      do m=1,2
      dzj=(dipole(-1,m)+dipole(1,m))*dcmplx(0.d0,-1.d0/sqr2)
      dipole(1,m)=(dipole(-1,m)-dipole(1,m))/sqr2
      dipole(-1,m)=dzj
      enddo
c     ----------------------------------------------------------------
c
c     ================================================================
c     green =>  ZZ(r,e)*tau(e) - ZJ(r,e)....................................
c     ----------------------------------------------------------------
      call mgreen(zj_flag,kkrsz,lmax,jws,nprpts,
     >            zlr_left,zlr_right,jlr,green,tau00,
     >            iprint,istop)
c     if(iprint.ge.0)write(6,'(''green'',6d14.6)')green(501),green(jws)
c     ----------------------------------------------------------------
c
c     ================================================================
      if(istop.eq.sname) then
        call fstop(sname)
      else
        return
      endif
c
      end
c
