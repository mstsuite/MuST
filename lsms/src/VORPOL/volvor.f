c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine volvor(rcrit,ncrit,
     >                  lmax,clm,
     >                  xgq,wgq,ngaussq,ngaussr,
     >                  xp,xpt,nvplane,
     >                  cornz,dc2,indxc,ncorn,
     >                  edge,edgp,nedge,
     >                  tnode,
     >                  grwylm,gwwylm,wylm,omegint,plm,sumfi)
c     ================================================================
c
      implicit   none
c
      integer    lmax
      integer    ncrit
      integer    ngaussr
      integer    ngaussq
      integer    nvplane
      integer    ncorn
      integer    nedge
      integer    i,j,np,ng
      integer    np_index(2)
      integer    sigma
c
      real*8     rcrit(ncrit)
      real*8     clm((lmax+1)*(lmax+2)/2)
      real*8     plm((lmax+1)*(lmax+2)/2)
      real*8     xgq(ngaussq)
      real*8     wgq(ngaussq)
      real*8     xp(3,nvplane)
      real*8     xpt(3,nvplane)
      real*8     cornz(ncorn)
      real*8     dc2(ncorn)
      integer    indxc(ncorn)
      real*8     edge(3,nedge)
      real*8     edgp(3,nedge)
      real*8     rp_index(2)
      real*8     tnode(nvplane+2*nedge+2)
      real*8     rgauss
      real*8     rgausq
      real*8     omegint
      real*8     tol
      real*8     half
      real*8     zero
      real*8     one
      real*8     pi,rfpi
      parameter (pi=3.141592653589793d0)
      real*8     grwylm(ngaussr,ncrit-1)
      real*8     gwwylm(ngaussr,ncrit-1)
c
      complex*16 wylm((lmax+1)*(lmax+2)/2,ngaussr,ncrit-1)
      complex*16 sumfi(0:lmax,3)
c
      parameter (tol=1.0d-10) 
      parameter (half=0.5d0)
      parameter (zero=0.0d0)
      parameter (one=1.0d0)
c
      rfpi=sqrt(4.d0*pi)
c     ================================================================
c     generate Gaussian-Quadrature points for theta integration.......
c     ----------------------------------------------------------------
!     call gauleg(-one,one,xgq,wgq,ngaussq)
      call gauss_legendre_points(-one,one,xgq,wgq,ngaussq)
c     ----------------------------------------------------------------
c     if(iprint.ge.1) then
c        write(6,'(''           No. Gaussian pts. [cos(theta)]: '',
c    >                                               i3)') ngaussq
c        write(6,'(''           n='',i3,'' xgq='',f15.8,
c    >   '' wgq='',d16.8)') (n,xgq(n),wgq(n),n=1,ngaussq)
c     endif
c
c     ================================================================
c     generate gaussian pointsi and weights for "r" integration.......
c     use the last column of grwylm and gwwylm to store xgr and wgr
c     ----------------------------------------------------------------
!     call gauleg(-one,one,grwylm(1,ncrit-1),gwwylm(1,ncrit-1),ngaussr)
      call gauss_legendre_points
     &          (-one,one,grwylm(1,ncrit-1),gwwylm(1,ncrit-1),ngaussr)
c     ----------------------------------------------------------------
c     if(iprint.ge.1) then
c        write(6,'(''           No. Gaussian pts.          [r]: '',
c    >                                               i3)') ngaussr
c        write(6,'(''           n='',i3,'' xgr='',f15.8,
c    >   '' wgr='',d16.8)') (n,grwylm(n,ncrit-1),gwwylm(n,ncrit-1),
c    &   n=1,ngaussr)
c     endif
      i=0
      np_index(1)=0
      np_index(2)=0
      rp_index(1)=zero
      rp_index(2)=zero
      do np=1,nvplane
         if(abs(xpt(1,np)).le.tol .and. abs(xpt(2,np)).le.tol) then
            i=i+1
            np_index(i)=np
            if(sigma(zero,zero,xpt(3,np),xpt,nvplane,1).eq.1)
     &        then
               do j=ncorn,1,-1
                  if(abs(cornz(indxc(j))-xpt(3,np)).lt.tol) then
                     rp_index(i)=sqrt(dc2(indxc(j)))     ! assuming dc2 is ordered.
                  endif
               enddo
            endif
         endif
      enddo
c
c     ================================================================
c     loop over the regions between critical r-points.................
c     ================================================================
ccdir$l cncall
      do i=2,ncrit
c        =============================================================
c        loop over the gaussian mesh for current region...............
c        =============================================================
         do ng=1,ngaussr
            rgauss=half*(rcrit(i)-rcrit(i-1))*grwylm(ng,ncrit-1) +
     >             half*(rcrit(i)+rcrit(i-1))
            grwylm(ng,i-1)=rgauss
            rgausq=rgauss*rgauss
c           ==========================================================
c           analytic forms for simple cubic...........................
c           if(i.eq.1) then
c              wstefa=rfpi
c           endif
c           if(i.eq.2) then
c              wstefa=rfpi*(three/rgauss-two)
c           endif
c           if(i.eq.3) then
c              rfac1=sqrt(rgausq-two)
c              rfac2=sqrt(three*rgausq-six)
c              rfac3=sqrt(two*rgausq-two)
c              wstefa=(48.0d0/rfpi)*
c    >                ( (one/rgauss)*atan((one-rfac1)/(one+rfac1)) +
c    >                asin((rfac2-rgauss)/(two*rfac3)) )
c           endif
c
c           ==========================================================
c           calculate w(l,m) for current r-point......................
c           ----------------------------------------------------------
            call stepyll(rgauss,rcrit(1),wylm(1,ng,i-1),
     >                   lmax,clm,
     >                   ngaussq,xgq,wgq,
     >                   xp,xpt,nvplane,
     >                   dc2(indxc(ncorn)),
     >                   edge,edgp,nedge,
     >                   np_index,rp_index,tnode,plm,sumfi)
c           ----------------------------------------------------------
c           ==========================================================
c           calculate the gaussian weight.............................
            gwwylm(ng,i-1)= half*(rcrit(i)-rcrit(i-1))*
     &          gwwylm(ng,ncrit-1)*rgausq
         enddo
      enddo
c
c     ================================================================
c     calculate the cell volume by integrating wylm(l=0,m=0)..........
c     ================================================================
      omegint=zero
      do i=2,ncrit
         do ng=1,ngaussr
            omegint=omegint + rfpi*
     >                        gwwylm(ng,i-1)*dble(wylm(1,ng,i-1))
         enddo
      enddo
c
c     do i=2,ncrit
c       write(6,'(/,''VOLVOR::i,ng,gwwylm(ng),wylm(1,ng),wylm(2,ng)'',
c    &      i5)') i
c           do ng=1,ngaussr
c              write(6,'(i3,d12.5,1p4e15.6)') ng,gwwylm(ng,i-1),
c    &          wylm(1,ng,i-1),wylm(2,ng,i-1)
c           enddo
c     enddo
c     ================================================================
c
      return
      end
