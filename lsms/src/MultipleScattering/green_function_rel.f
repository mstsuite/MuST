c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine green_function_rel(mtasa,
     >                          lmax,kkrsz,wx,wy,wz,
     >                          rins,r_mesh,jmt,jws,h,
     >                          pnrel,tau00,matom,
     >                          gz,fz,gj,fj,nuz,indz,
     &                          iprpts,
     >                          ngaussr,
     >                          cgnt,lmax_cg,
     >                          dos,dosck,green,dipole,
     >                          dos_orb,dosck_orb,dens_orb,
     >                          iprint,istop)
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
c                green   (Green's function) (actually complex
!                              charge & spin - magnetization densities)
!                dos_orb    )
!                dosck_orb  ) orbital dos & magn. densities
!                dens_orb   )
c     ****************************************************************
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     include    'atom_param.h'
      include    'gfill.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  istop*32
      character  sname*32

      integer nuzp
      parameter (nuzp=2)
c
      integer iprpts
      integer lmax_cg
      
      integer    mtasa
      integer    lmax
      integer    kkrsz
      integer    jmt
      integer    jws
      integer    m
      integer    ir
      integer    i
      integer    ngaussr
      integer    iprint
      integer    kmy,kmyp
      integer    nu,nup
      integer    kmy1,kmyp1
      integer    nuz(2*kkrsz)
      integer    indz(nuzp,2*kkrsz)
      integer    lambda
      integer    kmymax
c
      real*8     rsimp
      real*8     rins,r_mesh(iprpts),h
      real*8     cgnt(lmax_cg+1,(lmax_cg+1)**2,(lmax_cg+1)**2)
      real*8     rrr(iprpts+1),rir(iprpts+1)
      real*8     rri(iprpts+1),rii(iprpts+1)
      real*8     r1rr,r1ir
      real*8     r1ri,r1ii
      real*8     rlam
      real*8     fnpi
      real*8     fac
      real*8     pi
      real*8     two
      real*8     small
      parameter (two=2.0d0)
      parameter (small=1.d-15)
c
      complex*16 wx(4)
      complex*16 wy(4)
      complex*16 wz(4)
      complex*16 pnrel
      complex*16 matom(lmax+1,4)
      complex*16 tau00(kkrsz*2,kkrsz*2)
      complex*16 gz(iprpts,nuzp,kkrsz*2)
      complex*16 fz(iprpts,nuzp,kkrsz*2)
      complex*16 gj(iprpts,nuzp,kkrsz*2)
      complex*16 fj(iprpts,nuzp,kkrsz*2)
      complex*16 sum(kkrsz*2,kkrsz*2,2)
      complex*16 gf(kkrsz*2,kkrsz*2)
      complex*16 cxr,cxi
      complex*16 cox
      complex*16 zmom(kkrsz*2)
      complex*16 dos(4)
      complex*16 dosck(4)
      complex*16 green(jws,4)
      complex*16 dos_orb(3)
      complex*16 dosck_orb(3)
      complex*16 dens_orb(jws,3)
      complex*16 dipole(6,4)
c
      parameter (sname='green_function_rel')
c

! for the time beeing i set dipole to zero:
      call zeroout(dipole,6*4*2)

!     gfill & gafill are now called from zplanint_d
!     call gfill
!     call gafill

      pi=fnpi()

      lambda = 0
      fac = -1.d0/pi
      kmymax=2*(lmax+1)*(lmax+1)

!     dos(1), dos(4)
      do kmy=1,kmymax
        do kmyp=1,kmymax
          call zeroout(rrr,jmt)
          call zeroout(rir,jmt)
          call zeroout(rri,jmt)
          call zeroout(rii,jmt)
          do nu=1,nuz(kmy)
            kmy1=indz(nu,kmy)
            do nup=1,nuz(kmyp)
              kmyp1=indz(nup,kmyp)
              cox=fac*rgacoeff(kmy1,kmyp1,1)
              if(cdabs(cox).gt.small) then
                do i=1,jmt
!                 rlam=r_mesh(i)**lambda
                  rlam=1.d0
                  cxr=(gz(i,nu,kmy)*gz(i,nup,kmyp)+
     >                fz(i,nu,kmy)*fz(i,nup,kmyp))*cox*rlam
                  cxi=(gz(i,nu,kmy)*gj(i,nup,kmyp)+
     >                fz(i,nu,kmy)*fj(i,nup,kmyp))*cox*rlam
                  rrr(i)=rrr(i)+dreal(cxr)
                  rir(i)=rir(i)+dimag(cxr)
                  rri(i)=rri(i)+dreal(cxi)
                  rii(i)=rii(i)+dimag(cxi)
                end do
              end if

            end do ! nup
          end do   ! nu

          r1rr = rsimp(rrr,r_mesh,jmt,h)
          r1ir = rsimp(rir,r_mesh,jmt,h)
          r1ri = rsimp(rri,r_mesh,jmt,h)
          r1ii = rsimp(rii,r_mesh,jmt,h)

          sum(kmy,kmyp,1) = dcmplx(r1rr,r1ir)
          sum(kmy,kmyp,2) = dcmplx(r1ri,r1ii)


        end do ! kmyp
      end do   ! kmy

      call repl(gf,sum(1,1,1),kmymax,kmymax)
      call doubmt(gf,tau00,kmymax,kmymax)
      call submat(gf,sum(1,1,2),kmymax,kmymax)

!     call replms(zmom,gf,lmax,kmymax)

      dos(1) = 0.d0
      do i=1,kmymax
        dos(1) = dos(1)+gf(i,i)
      end do

      call magnet(lmax,kkrsz,
     >            rins,r_mesh,jmt,h,
     >            tau00,
     >            gz,fz,gj,fj,nuz,indz,
     >            sxcoeff,sxbcoeff,
     >            dos(2),
     &            iprpts,iplmax,
     >            iprint,istop)

      call magnet(lmax,kkrsz,
     >            rins,r_mesh,jmt,h,
     >            tau00,
     >            gz,fz,gj,fj,nuz,indz,
     >            sycoeff,sybcoeff,
     >            dos(3),
     &            iprpts,iplmax,
     >            iprint,istop)

      call magnet(lmax,kkrsz,
     >            rins,r_mesh,jmt,h,
     >            tau00,
     >            gz,fz,gj,fj,nuz,indz,
     >            szcoeff,szbcoeff,
     >            dos(4),
     &            iprpts,iplmax,
     >            iprint,istop)

      call magnet(lmax,kkrsz,
     >            rins,r_mesh,jmt,h,
     >            tau00,
     >            gz,fz,gj,fj,nuz,indz,
     >            lxcoeff,lxbcoeff,
     >            dos_orb(1),
     &            iprpts,iplmax,
     >            iprint,istop)


      call magnet(lmax,kkrsz,
     >            rins,r_mesh,jmt,h,
     >            tau00,
     >            gz,fz,gj,fj,nuz,indz,
     >            lycoeff,lybcoeff,
     >            dos_orb(2),
     &            iprpts,iplmax,
     >            iprint,istop)

      call magnet(lmax,kkrsz,
     >            rins,r_mesh,jmt,h,
     >            tau00,
     >            gz,fz,gj,fj,nuz,indz,
     >            lzcoeff,lzbcoeff,
     >            dos_orb(3),
     &            iprpts,iplmax,
     >            iprint,istop)

!     call green_4(lmax,kkrsz,
!    >          jmt,gz,fz,gj,fj,nuz,indz,tau00,
!    >          green,dens_orb,iprint,istop)

      call new_dens(lmax,kkrsz,
     >          jmt,gz,fz,gj,fj,nuz,indz,tau00,
     >          green(1,1),iprpts,iprint,istop)

      call magnetic_dens(lmax,kkrsz,
     >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
     >                sxcoeff,sxbcoeff,
     >                green(1,2),
     &                iprpts,
     >                iprint,istop)
 
      call magnetic_dens(lmax,kkrsz,
     >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
     >                sycoeff,sybcoeff,
     >                green(1,3),
     &                iprpts,
     >                iprint,istop)
 
      call magnetic_dens(lmax,kkrsz,
     >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
     >                szcoeff,szbcoeff,
     >                green(1,4),
     &                iprpts,
     >                iprint,istop)
 
! for the time beeing it seems that we don't need the orbital monent
! densities, so I just set them to zero...
      call zeroout(dens_orb,2*jws*3)

!     call magnetic_dens(lmax,kkrsz,
!    >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
!    >                lxcoeff,lxbcoeff,
!    >                dens_orb(1,1),
!    &                iprpts,
!    >                iprint,istop)
!
!     call magnetic_dens(lmax,kkrsz,
!    >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
!    >                lycoeff,lybcoeff,
!    >                dens_orb(1,2),
!    &                iprpts,
!    >                iprint,istop)
!
!     call magnetic_dens(lmax,kkrsz,
!    >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
!    >                lzcoeff,lzbcoeff,
!    >                dens_orb(1,3),
!    ?                iprpts,  
!    >                iprint,istop)


! currently we are have only implemented ASA for the relativistic case,
! therefore dos == dosck

      dosck(1) = dos(1)
      dosck(2) = dos(2)
      dosck(3) = dos(3)
      dosck(4) = dos(4)
      dosck_orb(1) = dos_orb(1)
      dosck_orb(2) = dos_orb(2)
      dosck_orb(3) = dos_orb(3)

c$$$      if(iprint.ge.0) then
c$$$        write(6,*) 'dos(1)=',dos(1)
c$$$        write(6,*) 'dos(2)=',dos(2)
c$$$        write(6,*) 'dos(3)=',dos(3)
c$$$        write(6,*) 'dos(4)=',dos(4)
c$$$        write(6,*) 'dos_orb(1)=',dos_orb(1)
c$$$        write(6,*) 'dos_orb(2)=',dos_orb(2)
c$$$        write(6,*) 'dos_orb(3)=',dos_orb(3)
c$$$      end if

      if(istop.eq.sname) then
        call fstop(sname)
      end if

      end
