c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine semcst(nqn,lqn,kqn,en,rv,r,rf,h,z,c,
     >                  nitmax,tol,nmt,nws,nlast,iter,iprpts,ipdeq)
c     ================================================================
c
c     ****************************************************************
c     semi-core states solver ..........................bg...june 1990
c     energy eigenvalues are found by newton-raphson method matching 
c     at the classical inversion point the f/g ratio. the energy 
c     derivative of f/g is evaluated analitically: 
c     see rose's book, chapter 4, pages 163-169, 
c     and also jp desclaux code.
c
c     calls: outws(invals) and inws
c
c     variables explanation:
c
c               nqn:      principal quantum number; 
c               lqn:      orbital quantum number;
c               kqn:      kappa quantum number; 
c               en:       energy; 
c               rv:       potential in rydbergs times r;
c               r:        radial log. grid; 
c               rg:       big component; 
c               rf:       small component;
c               h:        exp. step; 
c               z:        atomic number; 
c               nitmax:   number of iterations;
c               tol:      energy tolerance; 
c               nmt:      muffin-tin radius index;
c               nws:      bopunding sphere radius index;
c               nlast:    last tabulation point; 
c               c:        speed of light in rydbergs;
c               drg,drf:  wavefunctions derivatives times r;
c               gam:      first power for small r expansion; 
c               slp:      slope at the origin;
c               dm:       h/720;
c     ****************************************************************
      implicit   none
c 
!     include    'atom_param.h'
c
      integer    iprpts, ipdeq
      integer    nqn
      integer    lqn
      integer    kqn
      integer    nitmax
      integer    nmt
      integer    nws
      integer    nlast
      integer    iter
      integer    nmax
      integer    invp
      integer    imm
      integer    nodes
      integer    lll
      integer    j
c
      character  sname*32
c
      real*8     en
      real*8     h
      real*8     z
      real*8     c
      real*8     tol
      real*8     rg(iprpts)
      real*8     rf(iprpts)
      real*8     rv(iprpts)
      real*8     r(iprpts),rtmp(0:iprpts)
      real*8     tmp(0:iprpts)
      real*8     qmp(0:iprpts)
      real*8     drg(ipdeq*2)
      real*8     drf(ipdeq*2)
      real*8     dk
      real*8     dm
      real*8     gam
      real*8     sign
      real*8     slp
      real*8     elim
      real*8     rfm
      real*8     rgm
      real*8     rose
      real*8     gnrm
      real*8     fnrm
      real*8     de
      real*8     val
      real*8     enew
c
      parameter  (sname='semcst')
c

!     write (6,*) "semcst input en=",en
!     write (6,*) "  nmt,nws,nlast=",nmt,nws,nlast
c     ================================================================
c     initialize quantities
c     ================================================================
      iter=0
      dk=kqn
      dm=h/720.0d+00
      rtmp(0)=0.d0
      do j=1,iprpts
	rtmp(j)=sqrt(r(j))
      enddo
c
c     ================================================================
c     for semicore states you want to fix the inversion point at
c     the muffin-in, or ASA radius. set imm=1.........................
c     ================================================================
      invp=nmt
      imm=1
c
c     ================================================================
c     for MT case nmax=iprpts, but for ASA case nmax=jws.
c     this sets the charge within the proper regions.
c     ================================================================
      nmax=nlast
c
c     ================================================================
c     no. of nodes for the current states big component
c     ================================================================
      nodes=nqn-lqn
c
c     ================================================================
c     first power of small r expansion
c     ================================================================
      gam=sqrt(dk*dk-4.0d+00*z*z/(c*c))
c
c     ================================================================
c     slope of the wavefcns at the origine
c     ================================================================
      if(nodes-2*(nodes/2).ne.0) then
         sign= 1.0d+00
      else
         sign=-1.0d+00
      endif
      slp=sign*kqn/abs(kqn)
c
c     ================================================================
c     find the lowest energy limit
c     ================================================================
      lll=(lqn*(lqn+1))/2
      if (lll.ne.0) then
         elim=-z*z/(0.75d+00*nqn*nqn)
      else
         elim=(rv(1)+lll/r(1))/r(1)
         do j=2,iprpts
            elim=min((rv(j)+lll/r(j))/r(j),elim)
	 enddo
      endif
c
c     ================================================================
c     check potential
c     ================================================================
      if(elim.ge.0) then
         write(6,'('' SEMCST:: v+l*(l+1)/r**2 always positive'')')
         call fstop(sname)
      endif
c
      if(en.le.elim) en=elim*0.5d+00
c
c     ================================================================
c     routine outws performs the outward integration
c     ================================================================
c     ----------------------------------------------------------------
 2    continue
      call outws(invp,rg,rf,rv,r,en,c,drg,drf,elim,z,
     >           gam,slp,tol,imm,lll,dk,dm,nodes,nlast,iprpts,ipdeq)
c     ----------------------------------------------------------------
c
c     ================================================================
c     store the inversion point values
c     ================================================================
      rfm=rf(invp)
      rgm=rg(invp)
c
c     ================================================================
c     sets free solution directly to Riccati-Hankel for r*g and r*f.
c     ================================================================
c     ----------------------------------------------------------------
      call inwhnk(invp,iprpts,rg,rf,rv,r,en,c,drg,drf,
     >     dk,dm,nlast,imm,iprpts,ipdeq)
c     ----------------------------------------------------------------
c     ================================================================
c     routine inws performs the inward integration
c     ----------------------------------------------------------------
!      call inws (invp,nmax,rg,rf,rv,r,en,c,drg,drf,dk,dm,nws,imm,
!     >           iprpts,ipdeq)
c     ----------------------------------------------------------------
c
c     ================================================================
c     components match
c     ================================================================
      fnrm=rgm/rg(invp)
c
c     ================================================================
c     check first the sign of the big component
c     ================================================================
      if (fnrm.lt.0.0d+00) then
         write(6,'('' SEMCST:: wrong big component change'')')
         call fstop(sname)
      endif
c
      do j=invp,iprpts
         rg(j)=rg(j)*fnrm
         rf(j)=rf(j)*fnrm
      enddo   
c
c     ================================================================
c     energy derivative of the wvfcns "log. derivative"
c     ================================================================
      tmp(0)=0.d0
      do j=1,iprpts
         tmp(j)=(rg(j)**2+rf(j)**2)
      enddo   
c     ----------------------------------------------------------------
      call newint(iprpts+1,rtmp,tmp,qmp,1)
c     ----------------------------------------------------------------
      rose=2.d0*qmp(iprpts)+h*r(invp)*(rfm*rfm-rf(invp)*rf(invp))/3.0d+0
c
c     ================================================================
c     energy update
c     ================================================================
      de=rg(invp)*(rfm-rf(invp))*c/rose
      val=abs(de/en)
      if (val.le.tol) go to 8
  6   continue
      enew=en+de
      if(enew.lt.0.0d+00) go to 7
      de=de+0.5d+00
      val=val*0.5d+00
      if (val.gt.tol) go to 6
c
c     ================================================================
c     just in case the energy becomes zero
c     ================================================================
      write(6,'('' SEMCST:: zero energy'')')
      write(6,'(" n=",i2,", l=",i2,", kappa=",i2)') nqn,lqn,kqn
      call fstop(sname)
c
c     ================================================================
c     not yet convergence: try again
c     ================================================================
  7   continue
      en=enew
      iter=iter+1
      if(iter.le.nitmax) go to 2
c
c     ================================================================
c     not converged, too small tolerance or too small nitmax
c     ================================================================
      write(6,'('' SEMCST:: warning: too many energy iterations'')')
c
c     ================================================================
c     got convergence: exit
c     ================================================================
  8   continue
      do j=1,iprpts
         rf(j)=rf(j)*rf(j)+rg(j)*rg(j)
      enddo   
c
c     ================================================================
c     end of job: waiting for next state
c     ================================================================
      return
      end
