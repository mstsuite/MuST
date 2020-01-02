c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine congauss(ebot,etop,eibot,egrd,dele1,ngauss,nume,
     >                    pi,ipepts,iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  sname*20
      character  istop*32
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     include    'atom_param.h'
      integer ipepts
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer    nume
      integer    ngauss
      integer    ie
      integer    iprint
c
      real*8     xgs(ipepts)
      real*8     wgs(ipepts)
      real*8     ebot
      real*8     etop
      real*8     eibot
      real*8     pi
      real*8     erad
      real*8     half
      real*8     one
c
      complex*16 egrd(ipepts)
      complex*16 dele1(ipepts)
      complex*16 zzero
      complex*16 ihalfpi
      complex*16 czero
      complex*16 cone
      complex*16 sqrtm1
c
      parameter  (half=0.5d0)
      parameter  (one=1.0d0)
      parameter  (czero=(0.0d0,0.0d0))
      parameter  (cone=(1.0d0,0.0d0))
      parameter  (sqrtm1=(0.0d0,1.0d0))
      parameter  (sname='congauss')
c
c     ****************************************************************
c     constructs semi-circle energy contour
c
c                               . ! . Im(e)
c                           .     !     .
c                       .         !         .
c                     .           !           .
c                   .             !             .
c                  .              !              .
c          -----------------------!-------------------------
c              ebot               !    Real(e)   etop
c
c
c     ebot      : bottom of the contour on the real axis..........
c     etop      : top of the contour on the real axis.............
c     ngauss    : total number of desired points and weights for...
c                 Gauss-Legendre quadrature.......................
c     nume      : = ngauss
c     iprint    : print level control (>0 for print of energies)..
c
c     requires  : ebot,etop,ngauss,iprint
c     returns   : (egrd(ie),dele1(ie),ie=1,nume)
c     ****************************************************************
c
c     ================================================================
c     set up the gaussian integration points..........................
      if(ngauss+1.gt.ipepts) then
         write(6,'(/,'' CONGAUSS:: ngauss+1.gt.ipepts'',2i5)')
     >                             ngauss,ipepts
         call fstop(sname)
      endif
c
c     ================================================================
c     set nume equal to number of points on grid i.e. ngauss..........
      nume=ngauss
c     ----------------------------------------------------------------
      call gauss_legendre_points(-one,one,xgs,wgs,ngauss)
c     ----------------------------------------------------------------
c
c     ================================================================
c     loop over energy grid for gaussian integration.................
      erad  = half*(etop-ebot)
      ihalfpi=sqrtm1*half*pi
      do ie=1,ngauss
         zzero =  erad*exp( ihalfpi*(one-xgs(ie)) )
         egrd(ie) = half*(etop+ebot) + zzero
         dele1(ie)=-ihalfpi*zzero*wgs(ie)
      enddo
c
c     ================================================================
c     add an extra energy point to be used in extrapolating E_f.......
      nume=nume+1
      egrd(nume)=cone*etop+sqrtm1*eibot
      dele1(nume)=czero
c
c     ================================================================
c     write out grid and return.......................................
      if(iprint.ge.1) then
         write(6,'(/,'' CONGAUSS:: Bottom of contour'',
     >             t40,''='',f10.5)')   ebot
         write(6,'(12x,''Top    of contour'',
     >             t40,''='',f10.5)')   etop
         write(6,'(12x,''Number of energies'',
     >             t40,''='',i5)') nume
      endif
      if(iprint.ge.2) then
         write(6,'(12x,''n='',i5,'' e='',2e12.4,'' de='',2d12.4)')
     >   (ie,egrd(ie),dele1(ie),ie=1,nume)
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
