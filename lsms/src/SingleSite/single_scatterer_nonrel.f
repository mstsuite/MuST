C> Non relativistic single scaterer
C> computes and stores the wave functions and the t matrix.
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine single_scatterer_nonrel(nrelv,clight,lmax,kkrsz,
     >                       energy,prel,pnrel,
     >                       vr,rr,h,jmt,jws,
     >                       tmat,matom,
     >                       zlr,jlr,
     >                       r_sph,iprpts,iprint,istop)
c     ================================================================
c
c     ****************************************************************
c     computes the matrix
c     computes integrals over the wigner-seitz cell.
c     revised by was for n-atoms per unit cell.
c     revised by gms to store tmat etc in matrix form
c     calls:
c            semreltz    solve radial schroedinger equation
c            bessel      get Bessel functions.
c     requires:
c            lmax        (0,1,2, or 3)
c            kkrsz       siz of kkr matrix [(lmax+1)**2]
c            istop       index of subroutine in which prgm will stop
c      returns:
c            tmat        t-matrix
c            pzz         integral of (reg. wf.)**2 over ws cell
c            pzj         integral of (reg.wf.*irreg.wf) over ws cell
c            pzzck       integral of (reg. wf.)**2 over mt sphere
c            pzjck       integral of (reg.wf.*irreg.wf) over mt sphere
c     ****************************************************************
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     include    'atom_param.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer iprpts

      character  sname*32
      character  istop*32
c
      integer    nrelv
      integer    lmax
      integer    kkrsz
      integer    jmt
      integer    jws
      integer    iprint
      integer    l
      integer    m
      integer    lm
c
      real*8     vr(iprpts)
      real*8     rr(iprpts)
      real*8     clight
      real*8     h,r_sph
c
      complex*16 matom(0:lmax)
      complex*16 tmat(kkrsz,kkrsz)
      complex*16 energy
      complex*16 prel
      complex*16 pnrel
      complex*16 zlr(iprpts,0:lmax)
      complex*16 jlr(iprpts,0:lmax)
      complex*16 cone
      complex*16 sqrtm1
c
      parameter (sname='single_site')
      parameter (cone=(1.0d0,0.0d0))
      parameter (sqrtm1=(0.0d0,1.0d0))
c
c     =================================================================
c     printout if needed...............................................
c YingWai's check
      if(iprint.ge.1) then
         write(6,'('' SINGLE_SITE:: nrelv,lmax: '',2i5)') nrelv,lmax
         write(6,'(''               energy:     '',4d12.5)') energy
         write(6,'(''               jmt,jws:    '',2i5)') jmt,jws
         write(6,'('' rr'',9d14.7)')r_sph,rr(jmt-2),rr(jmt-1),rr(jmt)
         write(6,'('' rr506'',9d14.7)')r_sph,rr(506-2),rr(506-1),rr(506)
      endif
c
c     =================================================================
c     calculate matom and zlr and jlr.................................
c     -----------------------------------------------------------------
c YingWai's check
!      write(6,*) "Inside single_scatterer_nonrel. Before semrel."
      call semrel(nrelv,clight,lmax,
     >            energy,prel,pnrel,
     >            matom,zlr,jlr,
     >            h,jmt,jws,rr,vr,
c    >            h,506,506,rr,vr,
     >            r_sph,
     >            iprpts,
     >            iprint,istop)
!      write(6,*) "Inside single_scatterer_nonrel. After semrel."
c     -----------------------------------------------------------------
c
c     =================================================================
c     calculate t-matrix and zlr and jlr...............................
c     -----------------------------------------------------------------
      call zeroout(tmat,2*kkrsz*kkrsz)
c     -----------------------------------------------------------------
      lm=0
      do l=0,lmax
         do m=-l,l
            lm=lm+1
            tmat(lm,lm)=cone/matom(l)
         enddo
c        --------------------------------------------------------------
c        call zaxpy(jws,-tmat(lm,lm),zlr(1,l),1,jlr(1,l),1)
c        --------------------------------------------------------------
      enddo
c
c     =================================================================
      if (istop.eq.sname) then
        call fstop(sname)
      else
        return
      endif
c
      end
