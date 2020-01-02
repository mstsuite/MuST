c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine magnetic_dens(lmax,kkrsz,
     >                jmt,gz,fz,gj,fj,nuz,indz,tau00,
     >                gcoeff,gbcoeff,
     >                out,
     &                iprpts,
     >                iprint,istop)
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
      include    'gfill.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  istop*32
      character  sname*32
c
      integer iprpts
!     integer kmymaxp
!     parameter (kmymaxp=2*(iplmax+1)*(iplmax+1))
      integer nuzp
      parameter (nuzp=2)

      integer    lmax
      integer    kmymax
      integer    kkrsz
      integer    jmt
      integer    m
      integer    i
      integer    kmy,kmyp
      integer    nu,nup
      integer    kmy1,kmyp1
      integer    nuz(kkrsz*2),indz(nuzp,kkrsz*2)
      integer    iprint
c
      real*8     cgnt(lmax+1,(lmax+1)**2,(lmax+1)**2)
      real*8     fnpi
      real*8     pi
      real*8     fac
      real*8     two
      real*8     small
      parameter (two=2.0d0)
      parameter (small=1.d-15)
c
      complex*16 wx(4)
      complex*16 wy(4)
      complex*16 wz(4)
      complex*16 cxr,cxi
      complex*16 pnrel
      complex*16 matom(lmax+1,4)
      complex*16 gcoeff(kmymaxp,kmymaxp)
      complex*16 gbcoeff(kmymaxp,kmymaxp)
      complex*16 gz(iprpts,nuzp,2*kkrsz),fz(iprpts,nuzp,2*kkrsz)
      complex*16 gj(iprpts,nuzp,2*kkrsz),fj(iprpts,nuzp,2*kkrsz)
      complex*16 tau00(kkrsz*2,kkrsz*2)
      complex*16 gf(kkrsz*2,kkrsz*2)
      complex*16 sum1(kkrsz*2,kkrsz*2)
      complex*16 sum2(kkrsz*2,kkrsz*2)
      complex*16 out(jmt)
      complex*16 gcox,gbcox
      complex*16 one,mone
c
      parameter (sname='magnetic_dens')
      parameter (one=(1.d0,0.d0))
      parameter (mone=(-1.d0,0.d0))
c
      pi=fnpi()

      fac = -1.d0/pi
      kmymax=2*(lmax+1)*(lmax+1)

      do i=1,jmt
        do kmy=1,kmymax
          do kmyp=1,kmymax
            cxr = (0.d0,0.d0)
            cxi = (0.d0,0.d0)
            do nu=1,nuz(kmy)
              kmy1=indz(nu,kmy)
              do nup=1,nuz(kmyp)
                kmyp1=indz(nup,kmyp)
!               gcox=fac*gcoeff(kmy1,kmyp1)
!               gbcox=fac*gbcoeff(kmy1,kmyp1)
                gcox=gcoeff(kmy1,kmyp1)
                gbcox=gbcoeff(kmy1,kmyp1)
                if(abs(gcox).gt.small) then
                  cxr=cxr+(gz(i,nu,kmy)*gz(i,nup,kmyp))*gcox
                  cxi=cxi+(gz(i,nu,kmy)*gj(i,nup,kmyp))*gcox
                end if
                if(abs(gbcox).gt.small) then
                  cxr=cxr-(fz(i,nu,kmy)*fz(i,nup,kmyp))*gbcox
                  cxi=cxi-(fz(i,nu,kmy)*fj(i,nup,kmyp))*gbcox
                end if

              end do ! nup
            end do   ! nu

            sum1(kmy,kmyp) = cxr
            sum2(kmy,kmyp) = cxi

          end do ! kmyp
        end do   ! kmy


!!      call repl(gf,sum1,kmymax,kmymax)
!!      call doubmt(gf,tau00,kmymax,kmymax)
!!      call submat(gf,sum2,kmymax,kmymax)

        call repl(gf,sum2,kmymax,kmymax)
        call zgemm('N','N',kmymax,kmymax,kmymax,one,sum1,kmymax,
     >             tau00,kmymax,mone,gf,kmymax)

        out(i) = (0.d0,0.d0)
        do kmy=1,kmymax
          out(i) = out(i) + gf(kmy,kmy)
        end do

      end do ! i

      if(istop.eq.sname) then
        do i=1,jmt
          write(6,*) i,out(i)
        end do
        call fstop(sname)
      end if

      end
