c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine magnet(lmax,kkrsz,
     >                          rins,r_mesh,jmt,h,
     >                          tau00,
     >                          gz,fz,gj,fj,nuz,indz,
     >                          gcoeff,gbcoeff,
     >                          out,
     &                          iprpts,
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
c     ****************************************************************
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     include    'atom_param.h'
      include 'gfill.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  istop*32
      character  sname*32
c
      integer iprpts
!      integer kmymaxp
!      parameter (kmymaxp=2*(iplmax+1)*(iplmax+1))

      integer    lmax
      integer    kmymax
      integer    kkrsz
      integer    jmt
      integer    m
      integer    ir
      integer    i
      integer    iprint
      integer    kmy, kmyp
      integer    nu,nup
      integer    kmy1,kmyp1
      
      integer    nuzp
      parameter  (nuzp=2)

      integer    nuz(2*kkrsz)
      integer    indz(nuzp,2*kkrsz)
c
      real*8     rsimp
      real*8     rins,r_mesh(iprpts)
      real*8     h
      real*8     cgnt(lmax+1,(lmax+1)**2,(lmax+1)**2)
      real*8     fnpi
      real*8     fac
      real*8     pi
      real*8     two
      real*8     small
      real*8     rrr(iprpts),rir(iprpts)
      real*8     rri(iprpts),rii(iprpts)
      real*8     r1rr,r1ri
      real*8     r1ir,r1ii
      parameter (two=2.0d0)
      parameter (small=1.d-15)
c
      complex*16 tau00(kkrsz*2,kkrsz*2)
      complex*16 gz(iprpts,nuzp,kkrsz*2)
      complex*16 fz(iprpts,nuzp,kkrsz*2)
      complex*16 gj(iprpts,nuzp,kkrsz*2)
      complex*16 fj(iprpts,nuzp,kkrsz*2)
      complex*16 sum(kkrsz*2,kkrsz*2,2)
      complex*16 gf(kkrsz*2,kkrsz*2)
      complex*16 out
      complex*16 cxr,cxi
      complex*16 gcoeff(kmymaxp,kmymaxp)
      complex*16 gbcoeff(kmymaxp,kmymaxp)
      complex*16 gcox,gbcox
c
      parameter (sname='magnet')
c
      pi=fnpi()

      fac = -1.d0/pi
      kmymax=2*(lmax+1)*(lmax+1)

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
              gcox=fac*gcoeff(kmy1,kmyp1)
              gbcox=fac*gbcoeff(kmy1,kmyp1)
              if(abs(gcox).gt.small) then
                do i=1,jmt
                  cxr=gz(i,nu,kmy)*gz(i,nup,kmyp)*gcox
                  cxi=gz(i,nu,kmy)*gj(i,nup,kmyp)*gcox
                  rrr(i)=rrr(i)+dreal(cxr)
                  rir(i)=rir(i)+dimag(cxr)
                  rri(i)=rri(i)+dreal(cxi)
                  rii(i)=rii(i)+dimag(cxi)
                end do
              end if
              if(abs(gbcox).gt.small) then
                do i=1,jmt
                  cxr=fz(i,nu,kmy)*fz(i,nup,kmyp)*gbcox
                  cxi=fz(i,nu,kmy)*fj(i,nup,kmyp)*gbcox
                  rrr(i)=rrr(i)-dreal(cxr)
                  rir(i)=rir(i)-dimag(cxr)
                  rri(i)=rri(i)-dreal(cxi)
                  rii(i)=rii(i)-dimag(cxi)
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

      out = (0.d0,0.d0)
      do i=1,kmymax
        out = out + gf(i,i)
      end do

      if(istop.eq.sname) then
        call fstop(sname)
      end if

      end
