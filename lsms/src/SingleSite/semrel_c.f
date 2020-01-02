c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine semrel(nrelv,clight,lmax,
     >                  energy,prel,pnrel,
     >                  matom,zlr,jlr,
     >                  h,jmt,jws,r,rv,
     >                  r_sph,
     >                  iprpts,
     >                  iprint,istop)
c     ================================================================
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     include    'atom_param.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  sname*20
      character  istop*32
      character  coerr*20
c
      integer    iprpts
!
      integer    nrelv
      integer    lmax
      integer    jmt
      integer    jws
      integer    iprint
      integer    iend
      integer    i,j
      integer    l
      integer    solver
      parameter (solver=2)
c     parameter (solver=1)
      integer    icmax
      parameter (icmax=5)
      integer nr
      parameter (nr=4)
c
c     real*8     phshft(iplmax+1)
      real*8     clight
      real*8     h,r_sph
      real*8     r(iprpts)
      real*8     rv(iprpts)
      real*8     coef(nr,iprpts)
c
      complex*16 energy
      complex*16 prel
      complex*16 pnrel
      complex*16 matom(0:lmax)
      complex*16 zlr(iprpts,0:lmax)
      complex*16 jlr(iprpts,0:lmax)
      complex*16 bjl(2,iprpts)
      complex*16 bnl(2,iprpts)
c$$$      complex*16 bes(2*(iplmax+1))
c$$$      complex*16 djl(iplmax+1)
c$$$      complex*16 dnl(iplmax+1)
      complex*16 bes(2*(lmax+1))
      complex*16 djl(lmax+1)
      complex*16 dnl(lmax+1)
      complex*16 g(iprpts)
      complex*16 f(iprpts)
      complex*16 gpl(iprpts)
      complex*16 fpl(iprpts)
      complex*16 pr
      complex*16 tandel
      complex*16 sqrtm1
c
      parameter (sqrtm1=(0.0d0,1.0d0))
      parameter (sname='semrel')
c
c     ****************************************************************
c     program to solve scalar relativistic equations for complex
c     energies using Calogero's phase method     ddj/fjp  ---july 1991
c
c     Gets radial integrals, t-matrix.
c
c     In this version rv is r * v
c
c     called by: grint
c     calls:     ricbes,scalar,cqxpup,fstop
c
c      requires:
c            lmax      lmax
c            energy    energy
c            istop     switch: if istop=21 stop at end of this routine
c            h,x,r     grid info (step for log-grid, log grid, r-grid)
c
c      returns:
c            tinvl     t-1 matrix
c            zlr       rad wf
c            jlr       ir rad wf
c
c     set c for doing semi-relativistic or non-relativistic calc.
c     .....c = speed of light in proper units..........
c     Ryd units: m=1/2  hbar=1 c=2*inverse of fine structure const.
c
c     Solving a relativistic (nrelv=0) or non-relativistic (nrelv=10)
c     problem.   If non-relativistic, then set e=e-elmass equal to e
c     (where elmass = electron mass) and 1/c is set to ~zero to make
c     relativistic terms zero.
c     NOTE:  c is set to c*10**nrelv to make cod
c            fully relativistic code..............................
c
c     Using Ricatti Bessel fcts. Need only l=0 and l=1 to get wave-fcts.
c
c     To get physical phase-shifts and normalization, need get Rl in
c     terms of spherical bessel for all l's at the boundary of sphere.
c     This is properly done in SCALAR, see clmt and slmt.
c
c     With M(r)= [1 + (e-v(r))*c2inv]
c     Obtain r*(wave-fct.):   regular     g=r*R and   f=r*R'/M
c                           irregular   gpl=r*R and fpl=r*R'/M
c
c     NOTE: ====>   need to get gderiv for mom. matrix elements
c
c     ****************************************************************
c
      if(iprint.ge.1) then
         write(6,'('' semrel: lmax   ='',i5)') lmax
         write(6,'('' semrel: nrelv,clight='',i5,d12.4)') nrelv,clight
         write(6,'('' semrel: energy ='',2d12.4)') energy
         write(6,'('' semrel: prel   ='',2d12.4)') prel
         write(6,'('' semrel: pnrel  ='',2d12.4)') pnrel
         write(6,'('' semrel: jmt,jws='',2i5)') jmt,jws
         write(6,'('' rr'',9d14.7)')r_sph,r(jmt-2),r(jmt-1),r(jmt)
         write(6,'(''rr506s'',9d14.7)')r_sph,r(506-2),r(506-1),r(506)
      endif
c
c     ================================================================
c     for muffin-tins incrs and jmt should be the same
c     ================================================================
c     iend=max(506+11,jws)
      iend=max(jmt+11,jws)
c     write(6,'('' semrel: jmt,jws,iend='',3i5)') jmt,jws,iend
      if( iend .lt. jws ) then
         coerr=' iend less than jws'
         write(6,'(2(2x,a20))')coerr,sname
         write(6,'('' semrel: jmt,jws,iend='',3i5)') jmt,jws,iend
         call fstop(sname)
      else if( iend .gt. iprpts ) then
         coerr='iend > iprpts'
         write(6,'(2(2x,a10))')coerr,sname
         call fstop(sname)
      endif
      if(solver.eq.1) then
      do j=1,iend
          pr=pnrel*r(j)
          bjl(1,j)= sin(pr)
          bnl(1,j)=-cos(pr)
          bjl(2,j)= bnl(1,j) + bjl(1,j)/pr
          bnl(2,j)=-bjl(1,j) + bnl(1,j)/pr
      enddo
      else
c fit the potential to (A1*r+A2)/(1+A3*r)
      call fitpot(r,rv,coef,nr,iend-1,iend)
c     call fitpot(r,rv,coef,nr,506,iend)
c     do i=506,iend
c       do j=1,4
c         coef(j,i)=0.d0
c       enddo
c     enddo
      if(iprint.ge.2)then
      do i=jws-20,iend
      write(6,'(''coef'',i5,f10.5,9d14.6)')
     >i,r(i),rv(i),(coef(j,i),j=1,4)
      enddo
      endif
c     call fitpot(r,rv,coef,nr,506,iend)
c     call fitpot(r,rv,coef,nr,jmt,iend)
      endif
c
c     ================================================================
c     solve scalar-relativistic or non-relativistic equation.
c     ================================================================
      do l=0,lmax
        if(solver.eq.1) then
c        -------------------------------------------------------------
         call scalarr(nrelv,clight,l,
c        call scalar(nrelv,clight,l,
     >               bjl,bnl,bes,bes(lmax+2),djl,dnl,g,f,gpl,fpl,
     >               tandel,matom(l),
     >               energy,prel,pnrel,
     >               rv,r,h,jmt,iend,icmax,
c    >               rv,r,h,506,iend,icmax,
     >               r_sph,iprint,istop)
         call scalari(nrelv,clight,l,
     >               bjl,bnl,bes,bes(lmax+2),djl,dnl,g,f,gpl,fpl,
     >               tandel,matom(l),
     >               energy,prel,pnrel,
     >               rv,r,h,jmt,iend,icmax,
c    >               rv,r,h,506,iend,icmax,
     >               r_sph,iprint,istop)
c        -------------------------------------------------------------
c
c        =============================================================
c        store cotangent of phase-shift, normalization and t-matrix
c        =============================================================
c        phshft(l+1)=atan(dble(tandel))
c        =============================================================
c        get z_l(r): used in constructing green function..............
         do j=1,jws
            zlr(j,l)= g(j)/r(j)
         enddo
c        =============================================================
c        get j_l(r): used in constructing green function...........
         do j=1,jws
            jlr(j,l)= sqrtm1*prel*gpl(j)/(r(j)*matom(l))
         enddo
        else
         call rwave(nrelv.eq.0,.true.,g,gpl,f,fpl,matom(l),l,energy,
     &         coef,nr,r,(0.d0,0.d0),jmt,iend,1.d-14,1.d0,
     &         bes,2*(lmax+1),r_sph)
         do j=1,jws
            zlr(j,l)= g(j)*matom(l)/r(j)
            jlr(j,l)= sqrtm1*prel*gpl(j)/(r(j)*matom(l))
         enddo
        endif
      enddo
c
c     ================================================================
      if (istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
