c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine semrel(nrelv,clight,lmax,
     >                  energy,prel,pnrel,
     >                  matom,zlr,jlr,
     >                  h,jmt,jws,r,rv,
     >                  iprint,istop)
c     ================================================================
c
      implicit   none
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include    'atom_param.h'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      character  sname*20
      character  istop*32
      character  coerr*20
c
      integer    nrelv
      integer    lmax
      integer    jmt
      integer    jws
      integer    iprint
      integer    iend
      integer    i,j
      integer    l
      integer    solver
      parameter (solver=1)
      integer nr
      parameter (nr=4)
      integer n
      parameter (n=1)
c
c     real*8     phshft(iplmax+1)
      real*8     clight
      real*8     h
      real*8     r(iprpts)
      real*8     rv(iprpts)
      real*8     h0,dr
      real*8     r0(n*iprpts)
      real*8     rv0(n*iprpts)
      real*8     const
      real*8     coef(nr,n*iprpts)
c
      complex*16 energy
      complex*16 prel
      complex*16 pnrel
      complex*16 matom(0:lmax)
      complex*16 zlr(iprpts,0:lmax)
      complex*16 jlr(iprpts,0:lmax)
      complex*16 bjl(2,n*iprpts)
      complex*16 bnl(2,n*iprpts)
      complex*16 bes(2*(iplmax+1))
      complex*16 djl(iplmax+1)
      complex*16 dnl(iplmax+1)
      complex*16 g(n*iprpts)
      complex*16 f(n*iprpts)
      complex*16 gpl(n*iprpts)
      complex*16 fpl(n*iprpts)
      complex*16 pr
      complex*16 tandel
      complex*16 cone
      complex*16 sqrtm1
c
      parameter (cone=(1.0d0,0.0d0))
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
c            h,r     grid info (step for log-grid, log grid, r-grid)
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
      endif
c
c     ================================================================
c     for muffin-tins incrs and jmt should be the same
c     ================================================================
      iend=max(n*jmt+11,n*jws)
      if( iend .lt. n*jws ) then
         coerr=' iend less than n*jws'
         write(6,'(2(2x,a20))')coerr,sname
         write(6,'('' semrel: jmt,jws,iend='',3i5)') jmt,jws,iend
         call fstop(sname)
      else if( iend .gt. n*iprpts ) then
         coerr='iend > n*iprpts'
         write(6,'(2(2x,a10))')coerr,sname
         call fstop(sname)
      endif
      h0=h/n
c     z=dble(int(abs(rv(1))+0.1d0))
c     do j=1,jmt
crvi(j)=(z+rv(j))/r(j)
c     enddo
c fit the potential to (A1*r+A2)/(1+A3*r)
      call fitpot(r,rv,coef,nr,jmt,jmt)
      do j=1,iend
	i=(j-1)/n+1
	if(n.gt.1) then
	const=mod(j-1,n)/dble(n)
	r0(j)=r(i)*(r(i+1)/r(i))**const
	else
	r0(j)=r(j)
	endif
	if(r0(j).le.r(jmt)) then
	 dr=r0(j)-r(i)
         rv0(j)=(coef(1,i)*dr+coef(2,i))/(1.d0+dr*dr*coef(3,i))
	 if(nr.eq.4) rv0(j)=rv0(j)+coef(nr,i)
	else
	  rv0(j)=0.d0
	endif
      enddo
      if(solver.eq.1) then
      do j=1,iend
          pr=pnrel*r0(j)
          bjl(1,j)= sin(pr)
          bnl(1,j)=-cos(pr)
          bjl(2,j)= bnl(1,j) + bjl(1,j)/pr
          bnl(2,j)=-bjl(1,j) + bnl(1,j)/pr
      enddo
      else if(n.gt.1) then
      call fitpot(r0,rv0,coef,nr,n*jmt,iend)
      endif
c
c     ================================================================
c     solve scalar-relativistic or non-relativistic equation.
c     ================================================================
      do l=0,lmax
        if(solver.eq.1) then
c        -------------------------------------------------------------
         call scalar(nrelv,clight,l,
     >               bjl,bnl,bes,bes(iplmax+2),djl,dnl,g,f,gpl,fpl,
     >               tandel,matom(l),
     >               energy,prel,pnrel,
     >               rv0,r0,h0,n*jmt,iend,icmax,
     >               iprint,istop)
c        -------------------------------------------------------------
c
c        =============================================================
c        store cotangent of phase-shift, normalization and t-matrix
c        =============================================================
c        phshft(l+1)=atan(dble(tandel))
c        =============================================================
c        get z_l(r): used in constructing green function..............
         do j=1,jws
            zlr(j,l)= g(n*(j-1)+1)/r(j)
         enddo
c        =============================================================
c        get j_l(r): used in constructing green function...........
         do j=1,jws
            jlr(j,l)= sqrtm1*prel*gpl(n*(j-1)+1)/(r(j)*matom(l))
         enddo
C        write(6,*) matom(l)
C        write(6,*) zlr(jws,l),jlr(jws,l)
        else
         call rwave(nrelv.eq.0,.true.,g,gpl,f,fpl,matom(l),l,energy,
     &         coef,nr,r0,(0.d0,0.d0),n*jmt,iend,1.d-14,1.d0,
     &         bes,2*(iplmax+1))
         do j=1,jws
            zlr(j,l)= g(n*(j-1)+1)*matom(l)/r(j)
            jlr(j,l)= sqrtm1*prel*gpl(n*(j-1)+1)/(r(j)*matom(l))
         enddo
C       write(6,*) l,matom(l)
C       write(6,*) zlr(jws,l),jlr(jws,l)
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
