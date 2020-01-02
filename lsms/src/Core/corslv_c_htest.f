c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine corslv(jmt,jws,last,last2,h,r,rv,z,
     >                  numc,nc,lc,kc,ecore,
     >                  corden,semden,
     >                  ecorv,esemv,
     >                  n_spin_pola,nrelc,ndeepz,iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  sname*20
      character  istop*32
c
      include    'atom_param.h'
c
      integer    numc
      integer    nc(numc)
      integer    lc(numc)
      integer    kc(numc)
      integer    jmt
      integer    jws
      integer    n_spin_pola
      integer    nrelc
      integer    ndeepz
      integer    last
      integer    last2
      integer    nitmax
      integer    i
      integer    j
      integer    iter
      integer    iprint
      integer    ndeep
      integer    kappa
c
      real*8     ecore(numc)
      real*8     r(iprpts),rtmp(0:iprpts)
      real*8     rv(iprpts)
      real*8     corden(iprpts)
      real*8     semden(iprpts)
      real*8     f(0:iprpts)
      real*8     h
      real*8     z
      real*8     zero
      real*8     one
      real*8     two
      real*8     four
      real*8     ten
      real*8     onh
      real*8     tol
      real*8     fnstrc
      real*8     c
      real*8     esemv
      real*8     ecorv
      real*8     fac1
      real*8 gnrm,qmp(0:iprpts)
c
      parameter  (zero=0.0d+00)
      parameter  (one=1.0d+00)
      parameter  (two=2.0d+00)
      parameter  (four=4.0d+00)
      parameter  (ten=10.0d+00)
      parameter  (onh=137.0359895d+00)
      parameter  (tol=1.0d-10)
      parameter  (sname='corslv')
      integer solver
      parameter (solver=1)
      integer nr
      parameter (nr=4)
      real*8 potc(nr,iprpts)
      logical hydrogen   ! for testing the code
      parameter( hydrogen=.true.)
c
c     ****************************************************************
c     this routine drives the radial equations solution for the
c     core states. it must be called for each component(sublattice)
c     separately.
c
c     calls:  deepst(outws,inws)
c             semcor(outws,hank)
c
c     r= radial grid, rv=potential times r
c     ecore= core levels guesses: to be updated
c     g,f upper and lower components times r
c     corden= core charge density
c     nc= principal q. nos., lc= orbital q. nos., kc= kappa q.nos.
c     numc= no. of core states, z= atomic no., jmt=muffin-tin index
c     ndeepz= number of deep core electrons based on zcor
c     jws=top r index
c     ****************************************************************
c
      if(hydrogen) then
	do i=1,last
	  rv(i)=-2.d0*z
	enddo
	numc=min(numc,ndeepz)
      endif
      if (numc.le.0) then
         esemv=zero
         ecorv=zero
         return
      endif
      nitmax=50
      fnstrc=one/onh
      c=two/fnstrc
      c=c*ten**nrelc
      if(solver.eq.1) then
      else
      call fitpot(r,rv,potc,nr,jmt,last)
      endif
c
c     ================================================================
c     call deepst for each core state.................................
c     ================================================================
      esemv=zero
      ecorv=zero
      ndeep=0
      do i=1,numc
         fac1=(3-n_spin_pola)*abs(kc(i))
         ndeep=ndeep+fac1
C if(hydrogen) ecore(i)=-z*z/(nc(i)*nc(i))
	 if(solver.eq.1) then
c        -------------------------------------------------------------
         call deepst(nc(i),lc(i),kc(i),ecore(i),
     >               rv,r,f(1),h,z,c,nitmax,tol,jws,last,iter)
c        -------------------------------------------------------------
c
c        =============================================================
c        when core state energy >-10 then treat as semi-core..........
c        =============================================================
c08/04/95if( ecore(i) .ge. -ten ) then
         if(ndeep.gt.ndeepz)then
c           ----------------------------------------------------------
            call semcst(nc(i),lc(i),kc(i),ecore(i),
     >                  rv,r,f(1),h,z,c,nitmax,tol,jmt,jws,last2,iter)
c           ----------------------------------------------------------
         endif
	 else
	 if(nrelc.eq.0) then
	   kappa=kc(i)
	 else
	   kappa=lc(i)
	 endif
	 call srcore(nrelc.eq.0,ecore(i),f(1),nc(i),kappa,potc,nr,
     &        r,nitmax,tol,last,1.d0)
	 endif
c     ================================================================
c     normalize the wavefunctions
c     ================================================================
      f(0)=0.d0
      rtmp(0)=0.d0
      do j=1,last2
	rtmp(j)=sqrt(r(j))
      enddo
c     ----------------------------------------------------------------
      call newint(last2+1,rtmp,f,qmp,1)
c     ----------------------------------------------------------------
      gnrm=1.d0/(two*qmp(last2))
      do j=1,last2
         f(j)=f(j)*gnrm
      enddo    
         if(ndeep.gt.ndeepz)then
            do j=1,last
               semden(j)=semden(j) + fac1*f(j)
            enddo
            esemv=esemv+ecore(i)*fac1
	 else
            do j=1,last
               corden(j)= corden(j) + fac1*f(j)
            enddo
            ecorv=ecorv+ecore(i)*fac1
	 endif
      enddo
c
c     ================================================================
c     Major printout: core eigenvalues................................
c     ================================================================
      if(hydrogen) then
      if(iprint.ge.0) then
         write(6,'(/,'' CORSLV:: Eigenvalues:'',t25,''n'',t28,
     >               ''l'',t31,''k'',t42,''energy'')')
         do i=1,numc
	   if(nrelc.eq.0) then
	   fac1=z*fnstrc/(nc(i)-abs(kc(i))+sqrt(kc(i)*kc(i)-
     &        z*z*fnstrc*fnstrc))
	   ecorv=0.5d0*c*c*(1.d0/sqrt(1.d0+fac1*fac1)-1.d0)
	   else
	   ecorv=-z*z/(nc(i)*nc(i))
	   endif
            write(6,'(t23,i3,t26,i3,t29,i3,t35,1p2d20.12)') 
     >      nc(i),lc(i),kc(i),ecore(i),ecorv
         enddo
	 stop
      endif
      else
      if(iprint.ge.0) then
         write(6,'(/,'' CORSLV:: Eigenvalues:'',t25,''n'',t30,
     >               ''l'',t35,''k'',t48,''energy'')')
         do i=1,numc
            write(6,'(t23,i3,t28,i3,t33,i3,t41,d20.13)') 
     >      nc(i),lc(i),kc(i),ecore(i)
         enddo
      endif
      endif
      if(iprint.ge.0) then
         write(6,'(10x,''Energy     core:'',t40,''='',d16.8)')ecorv
         write(6,'(10x,''Energy semicore:'',t40,''='',d16.8)')esemv
      endif
c     ***************************************************************
      if(iprint.ge.1) then
      write(6,'(/,'' CORSOLV::'')')
         do j=1,last,50
            write(6,'('' j,corden,semden: '',i5,3d20.12)')
     >      j,corden(j),semden(j)
         enddo
      endif
c     ***************************************************************
c
c    ==================================================================
      if(istop.eq.sname) then
      call fstop(sname)
      endif
c
      return
      end
