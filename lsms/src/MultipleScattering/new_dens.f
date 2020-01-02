      subroutine new_dens(lmax,kkrsz,
     >                ns,gz,fz,gj,fj,nuz,indz,tau,
     >                zrho,iprpts,iprint,istop)
c===================
c
c input:  lmax  - maximum of angular momentum index
c         ns - as in 'readpot'
c         gz,fz - big and small component of regular radial solution * r
c         gj,fj - big and small component of irregular radial solution * r
c         nuz - no. of (kap',my') components for (kap,my)
c         indz - selects (kap',my') for (kap,my)
c         tau - site-diagonal tau matrix
c output: zrho - radial density distribution
c
!     implicit real*8 (a-h,o-z)
      implicit none
c
!      include 'atom_param.h'
      include 'gfill.h'
      integer nuzp
      parameter (nuzp=2)

      integer iprpts

      character*32 istop
      character*32 sname
      parameter (sname='dens')

      real*8  pi
      real*8  fac
      real*8  small
      real*8 rlam
c
      integer iprint
      integer kkrsz
      integer lmax,lmmax,kmymax
      integer nuz(2*kkrsz),indz(nuzp,2*kkrsz)
      integer ns,i
      integer kmy,kmyp
      integer nu,nup
      integer kmy1,kmyp1
c
      complex*16 gz(iprpts,nuzp,2*kkrsz),fz(iprpts,nuzp,2*kkrsz)
      complex*16 gj(iprpts,nuzp,2*kkrsz),fj(iprpts,nuzp,2*kkrsz)
      complex*16 zrho(ns)
      complex*16 tau(2*kkrsz,2*kkrsz)
      complex*16 sum(2*kkrsz,2*kkrsz,2)
      complex*16 gf(2*kkrsz,2*kkrsz)
      complex*16 cox
c
c
      data small/1.0d-15/
c
      pi=4.d0*datan(1.d0)
      fac=-1.d0/pi
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
      call zeroout(zrho,2*ns)
c
      do i=1,ns
        call zeroout(sum,2*2*kmymax*kmymax)
c
        do kmy=1,kmymax
        do kmyp=1,kmymax

          do nu=1,nuz(kmy)
            kmy1=indz(nu,kmy)
            do nup=1,nuz(kmyp)
              kmyp1=indz(nup,kmyp)
!             cox=fac*rgacoeff(kmy1,kmyp1,1)
              cox=    rgacoeff(kmy1,kmyp1,1)
              if(cdabs(cox).gt.small) then
!                 rlam=r_mesh(i)**lambda
                  rlam=1.d0
                sum(kmy,kmyp,1)=sum(kmy,kmyp,1)+
     1               (gz(i,nu,kmy)*gz(i,nup,kmyp)+
     >                fz(i,nu,kmy)*fz(i,nup,kmyp))*cox*rlam
                sum(kmy,kmyp,2)=sum(kmy,kmyp,2)+
     1               (gz(i,nu,kmy)*gj(i,nup,kmyp)+
     >                fz(i,nu,kmy)*fj(i,nup,kmyp))*cox*rlam
              end if

            end do  ! nup
          end do    ! nu
        end do      ! kmyp
        end do      ! kmy

        call repl(gf,sum(1,1,1),kmymax,kmymax)
        call doubmt(gf,tau,kmymax,kmymax)
        call submat(gf,sum(1,1,2),kmymax,kmymax)
c
        zrho(i)=(0.d0,0.d0)
        do kmy=1,kmymax
          zrho(i)=zrho(i)+gf(kmy,kmy)
        end do
!       zrho(i)=fac*sum1
c
      end do


      end
