c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine green_function(mtasa,n_spin_pola,n_spin_cant,
     >                          lmax,kkrsz,wx,wy,wz,
     >                          rins,r_sph,r_mesh,jmt,jws,
     >                          pnrel,tau00_l,matom,zlr,jlr,
     &                          nprpts,nplmax,
     >                          ngaussr,
     >                          cgnt,lmax_cg,
     >                          dos,dosck,green,dipole,
     &                          ncrit,grwylm,gwwylm,wylm,
     >                          iprint,istop)
c     ================================================================
c
c
c     ****************************************************************
c     input:
c                tau00_l
c                kkrsz   (size of KKR-matrix)
c                istop   (index of subroutine prog. stops in)
c     output:
c                dos     (wigner-seitz cell density of states)
c                dosck   (muffin-tin density of states)
c                green   (Green's function)
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
      integer    n_spin_pola
      integer    n_spin_cant
      integer    lmax
      integer lmax_cg,nprpts,nplmax
      integer    kkrsz
      integer    jmt
      integer    jws
      integer    m
      integer    ir
      integer    ngaussr
      integer    iprint
      integer    zsl
      integer    zsr
      integer    zj_flag
      integer    isp
      integer ncrit
c
      real*8     rins,r_sph,r_mesh(*)
      real*8     cgnt(lmax_cg+1,(lmax_cg+1)**2,(lmax_cg+1)**2)
      real*8     fnpi
      real*8     pi
      real*8     two
      parameter (two=2.0d0)

      real*8 grwylm(*),gwwylm(*)
c
      complex*16 wx(4)
      complex*16 wy(4)
      complex*16 wz(4)
      complex*16 pnrel
      complex*16 matom(nplmax+1,n_spin_cant)
      complex*16 tau00_l(kkrsz*kkrsz*n_spin_cant*n_spin_cant)
      complex*16 zlr(nprpts*(nplmax+1),n_spin_cant)
      complex*16 jlr(nprpts*(nplmax+1),n_spin_cant)
      complex*16 dos(n_spin_cant*n_spin_cant)
      complex*16 dosck(n_spin_cant*n_spin_cant)
      complex*16 green(jws,n_spin_cant*n_spin_cant)
      complex*16 dipole(6,n_spin_cant*n_spin_cant)
      complex*16 t1
      complex*16 t2
      complex*16 t3
      complex*16 t4

      complex*16 wylm(*)
c
      parameter (sname='green_function')
c
      pi=fnpi()
c
      do isp=1,n_spin_cant*n_spin_cant
         if(isp.eq.1) then
            zsl=1
            zsr=1
            zj_flag=1
         else if(isp.eq.2) then
            zsl=2
            zsr=1
            zj_flag=0
         else if(isp.eq.3) then
            zsl=1
            zsr=2
            zj_flag=0
         else if(isp.eq.4) then
            zsl=2
            zsr=2
            zj_flag=1
         endif
!         write(*,*) "green_function: tau00_l(kkrsz*kkrsz*(isp-1)+1)"
!         write(*,*) isp,tau00_l(kkrsz*kkrsz*(isp-1)+1)
c        --------------------------------------------------------------
         call gf_local(mtasa,zj_flag,lmax,kkrsz,
     >                 rins,r_sph,r_mesh,jmt,jws,
     >                 pnrel,tau00_l(kkrsz*kkrsz*(isp-1)+1),
     >                 matom(1,zsl),matom(1,zsr),
     >                 zlr(1,zsl),zlr(1,zsr),jlr(1,zsr),
     &                 nprpts,
     >                 ngaussr,cgnt,lmax_cg,
     >                 dos(isp),dosck(isp),green(1,isp),dipole(1,isp),
     &                 ncrit,grwylm,gwwylm,wylm,
     >                 pi,iprint,istop)
c        --------------------------------------------------------------
      enddo
c
c     ================================================================
c     re-calculate dos and green so that:
c
c     dos(1) = Tr[dos],           green(r,1) = Tr[green(r)]
c                 ---                             -----
c
c     dos(2) = Tr[dos*wx],        green(r,2) = Tr[green(r)*wx]
c                 --- --                          -----    --
c
c     dos(3) = Tr[dos*wy],        green(r,3) = Tr[green(r)*wy]
c                 --- --                          -----    --
c
c     dos(4) = Tr[dos*wz],        green(r,4) = Tr[green(r)*wz]
c                 --- --                          -----    --
c
c     ================================================================
      if(n_spin_cant.eq.2) then         ! spin canting case
         t1=dos(1)+dos(4)
         t2=dos(1)*wx(1)+dos(2)*wx(3)+dos(3)*wx(2)+dos(4)*wx(4)
         t3=dos(1)*wy(1)+dos(2)*wy(3)+dos(3)*wy(2)+dos(4)*wy(4)
         t4=dos(1)*wz(1)+dos(2)*wz(3)+dos(3)*wz(2)+dos(4)*wz(4)
         dos(1)=t1
         dos(2)=t2
         dos(3)=t3
         dos(4)=t4
         t1=dosck(1)+dosck(4)
         t2=dosck(1)*wx(1)+dosck(2)*wx(3)+dosck(3)*wx(2)+dosck(4)*wx(4)
         t3=dosck(1)*wy(1)+dosck(2)*wy(3)+dosck(3)*wy(2)+dosck(4)*wy(4)
         t4=dosck(1)*wz(1)+dosck(2)*wz(3)+dosck(3)*wz(2)+dosck(4)*wz(4)
         dosck(1)=t1
         dosck(2)=t2
         dosck(3)=t3
         dosck(4)=t4
         do ir=1,jws
            t1=green(ir,1)+green(ir,4)
            t2=green(ir,1)*wx(1)+green(ir,2)*wx(3)+
     >         green(ir,3)*wx(2)+green(ir,4)*wx(4)
            t3=green(ir,1)*wy(1)+green(ir,2)*wy(3)+
     >         green(ir,3)*wy(2)+green(ir,4)*wy(4)
            t4=green(ir,1)*wz(1)+green(ir,2)*wz(3)+
     >         green(ir,3)*wz(2)+green(ir,4)*wz(4)
            green(ir,1)=t1
            green(ir,2)=t2
            green(ir,3)=t3
            green(ir,4)=t4
         enddo
	 do m=1,6
         t1=dipole(m,1)+dipole(m,4)
         t2=dipole(m,1)*wx(1)+dipole(m,2)*wx(3)+dipole(m,3)*wx(2)+
     &      dipole(m,4)*wx(4)
         t3=dipole(m,1)*wy(1)+dipole(m,2)*wy(3)+dipole(m,3)*wy(2)+
     &      dipole(m,4)*wy(4)
         t4=dipole(m,1)*wz(1)+dipole(m,2)*wz(3)+dipole(m,3)*wz(2)+
     &      dipole(m,4)*wz(4)
         dipole(m,1)=t1
         dipole(m,2)=t2
         dipole(m,3)=t3
         dipole(m,4)=t4
	 enddo
      else if(n_spin_pola.eq.1) then   ! non-spin-polarized case
         dos(1)=two*dos(1)             ! factor 2.0 is included in DOS and
         dosck(1)=two*dosck(1)         ! Green function for spin
         do ir=1,jws
            green(ir,1)=two*green(ir,1)
         enddo
	 do m=1,6
	   dipole(m,1)=two*dipole(m,1)
	 enddo
      endif
c
c     ================================================================
      if(istop.eq.sname) then
        call fstop(sname)
      else
        return
      endif
c
      end
