! originally WAFU from the SKKR code
!
      subroutine single_scatterer_rel(ce,psq,lmax,kmymax,
     >                idpot,v0,vr,br,bopr,dx,ns,rs,
     >                tminv,gz,fz,gj,fj,nuz,indz,iflag,socsc,
     >                iprpts,iprint,istop)        
c====================
c
c input:  e - any complex energy
c         lmax  - maximum of angular momentum index
c         v0 - vacuum potential
c         idpot,vr,br,bopr,dx,ns,rs - as in 'readpot'
c output: tminv - inverse of single-site scattering matrix
c         gz,fz   - big and small component of regular radial solution * r
c         nuz - no. of (kap',my') components for (kap,my)
c         indz - selects (kap',my') for (kap,my)
c
      implicit real*8 (a-h,o-z)
c
!     include '../param.h'
!      include 'atom_param.h'
      parameter (nuzp=2)
c
      integer idpot
      character*32 sname
      character*32 istop

      integer   iprint

      dimension vr(ns),br(ns),bopr(ns,2)
      complex*16 tminv2(2,2)
      complex*16 tminv(kmymax,kmymax)
      complex*16 gz2(iprpts,2,2),fz2(iprpts,2,2)
      complex*16 gj2(iprpts,2,2),fj2(iprpts,2,2)
!     complex*16 gz(ns,nuzp,kmymax),fz(ns,nuzp,kmymax)
!     complex*16 gj(ns,nuzp,kmymax),fj(ns,nuzp,kmymax)
      complex*16 gz(iprpts,nuzp,kmymax),fz(iprpts,nuzp,kmymax)
      complex*16 gj(iprpts,nuzp,kmymax),fj(iprpts,nuzp,kmymax)
      complex*16 fb(0:lmax+1),fn(0:lmax+1),fh(0:lmax+1)
      complex*16 fb1(0:lmax+1),fn1(0:lmax+1),fh1(0:lmax+1)
      complex*16 ce,psq,p,cevac,psqvac,pvac,sk,dk,xk,react,sqrtm1
      dimension nuz(kmymax),indz(nuzp,kmymax)
      data sqrtm1/(0.d0,1.d0)/

      parameter(sname='single_scatterer_rel')
c
c     write(6,'(a10)') idpot
c     write(6,'(2f10.6,i5)') dx,rs,ns
c     write(6,'(4d20.10)') (vr(j),j=1,ns)
c     write(6,'(4d20.10)') (br(j),j=1,ns)
c
!     write(6,*) 'meis'
!     write(6,*)ce,psq,lmax,idpot,v0
!     write(6,*)vr(1),vr(ns),br(1),br(ns)
!     write(6,*)dx,ns,rs

!     write(6,*) 'entering ',sname

      kmax=2*lmax+1
!     kmymax=2*(lmax+1)*(lmax+1)
      call zeroout(tminv,2*kmymax*kmymax)
c
!     idpot .eq. 0 -> Vacuum
      if(idpot.eq.0) then
c
        p=cdsqrt(psq)
        cevac=ce-v0
        pvac=cdsqrt(cevac)
        call csbf(lmax+1,p,rs,fb,fn,fh)
        call csbf(lmax+1,pvac,rs,fb1,fn1,fh1)
        do k=1,kmax
          l=k/2
          if(k.eq.2*l) then
            kap=l
            lb=l-1
            j=2*l-1
          else
            kap=-l-1
            lb=l+1
            j=2*l+1
          end if
          sk=dcmplx(dfloat(l-lb),0.d0)
          xk=sk*cevac*fb1(lb)/pvac
          xk=xk/fb1(l)
          sk=sk*ce/p
          dk=(xk*fb(l)-sk*fb(lb))/(xk*fn(l)-sk*fn(lb))
          react=-dk/p
          do my=-j,j,2
            kapmy=2*kap*kap+kap+(my+1)/2
            tminv(kapmy,kapmy)=1.d0/react+sqrtm1*p
          end do
        end do
c
      else  
c
        do l=0,lmax
          kap1=l
          kap2=-l-1
          do my=-2*l-1,2*l+1,2
c
!         write(6,*) sname,1

            call spzwafu(socsc,ce,psq,l,my,vr,br,bopr,dx,ns,rs,
     >                   tminv2,gz2,fz2,gj2,fj2,iflag,iprpts,lmax)
c
!         write(6,*) sname,2

            if(iabs(my).eq.2*l+1) then
c
              kapmy=2*kap2*kap2+kap2+(my+1)/2
              tminv(kapmy,kapmy)=tminv2(2,2)
c
              if(iflag.eq.1) then
                nuz(kapmy)=1
                indz(1,kapmy)=kapmy
                do i=1,ns
                  gz(i,1,kapmy)=gz2(i,2,2)
                  fz(i,1,kapmy)=fz2(i,2,2)
                  gj(i,1,kapmy)=gj2(i,2,2)
                  fj(i,1,kapmy)=fj2(i,2,2)
                end do
              end if
c
            else
c
              kapmy1=2*kap1*kap1+kap1+(my+1)/2
              kapmy2=2*kap2*kap2+kap2+(my+1)/2
              tminv(kapmy1,kapmy1)=tminv2(1,1)
              tminv(kapmy1,kapmy2)=tminv2(1,2)
              tminv(kapmy2,kapmy1)=tminv2(2,1)
              tminv(kapmy2,kapmy2)=tminv2(2,2)
c
              if(iflag.eq.1) then
                nuz(kapmy1)=2
                nuz(kapmy2)=2
                indz(1,kapmy1)=kapmy1
                indz(2,kapmy1)=kapmy2
                indz(1,kapmy2)=kapmy2
                indz(2,kapmy2)=kapmy1
                do i=1,ns
                  gz(i,1,kapmy1)=gz2(i,1,1)
                  fz(i,1,kapmy1)=fz2(i,1,1)
                  gz(i,2,kapmy1)=gz2(i,2,1)
                  fz(i,2,kapmy1)=fz2(i,2,1)
                  gz(i,1,kapmy2)=gz2(i,2,2)
                  fz(i,1,kapmy2)=fz2(i,2,2)
                  gz(i,2,kapmy2)=gz2(i,1,2)
                  fz(i,2,kapmy2)=fz2(i,1,2)
                  gj(i,1,kapmy1)=gj2(i,1,1)
                  fj(i,1,kapmy1)=fj2(i,1,1)
                  gj(i,2,kapmy1)=gj2(i,2,1)
                  fj(i,2,kapmy1)=fj2(i,2,1)
                  gj(i,1,kapmy2)=gj2(i,2,2)
                  fj(i,1,kapmy2)=fj2(i,2,2)
                  gj(i,2,kapmy2)=gj2(i,1,2)
                  fj(i,2,kapmy2)=fj2(i,1,2)
                end do
              end if
c
            end if
c
          end do
        end do
c
      end if
c
      if(istop.eq.sname) then
        call fstop(sname)
      end if

      end
