c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deepst(nqn,lqn,kqn,en,rv,r,rf,h,z,c,
     >                  nitmax,tol,nws,nlast,iter,iprpts,ipdeq)
c     ================================================================
c
c     core states solver ...............................bg...june 1990
c     energy eigenvalues are found by newton-raphson method matching 
c     at the classical inversion point the f/g ratio. the energy 
c     derivative of f/g is evaluated analitically: see rose's book, 
c     chapter 4, pages 163-169, and also jp desclaux code.............
c
c     calls: outws(invals) and inws
c
c     variables explanation:
c
c               nqn:     principal quantum number; 
c               lqn:     orbital quantum number;
c               kqn:     kappa quantum number; 
c               en:      energy; 
c               rv:      potential in rydbergs times r;
c               r:       radial log. grid; 
c               rg:      big component; 
c               rf:      small component;
c               h:       exp. step; 
c               z:       atomic number; 
c               nitmax:  number of iterations;
c               tol:     energy tolerance; 
c               nws:     bounding sphere radius index;
c               nlast:   last tabulation point; 
c               c:       speed of light in rydbergs;
c               drg,drf: wavefunctions derivatives times r;
c               gam:     first power for small r expansion; 
c               slp:     slope at the origin;
c               dm:      h/720;
c     ****************************************************************
      implicit   none
c
!     include    'atom_param.h'
c
      character  sname*32
c
      integer    ipdeq
      integer    iprpts
      integer    nqn
      integer    lqn
      integer    kqn
      integer    nitmax
      integer    nws
      integer    nlast
      integer    iter
      integer    imm
      integer    nodes
      integer    lll
      integer    invp
      integer    j
      integer    nmax
      integer    ir
c
      real*8     rg(iprpts)
      real*8     rf(iprpts)
      real*8     rv(iprpts)
      real*8     r(iprpts),rtmp(0:iprpts)
      real*8     tmp(0:iprpts)
      real*8     qmp(0:iprpts)
      real*8     drg(ipdeq*2)
      real*8     drf(ipdeq*2)
      real*8     en
      real*8     h
      real*8     z
      real*8     c
      real*8     tol
      real*8     dk
      real*8     dm
      real*8     gam
      real*8     sign
      real*8     slp
      real*8     elim
      real*8     rfm
      real*8     rgm
      real*8     rose
      real*8     fnrm
      real*8     val
      real*8     de
      real*8     gnrm
      real*8     enew
c
      parameter  (sname='deepst')
c

!     write(*,*) "Entering DEEPST"
c     ================================================================  
c     check grid and see if core routine will work on this grid
c     ================================================================  
      if(nws+ipdeq*2 .gt. iprpts) then
         write(6,'('' DEEPST::core cannot work:'',
     >   '' nws+ipdeq2.gt.iprpts'')')
         write(6,'('' nws     = '',i6)') nws
         write(6,'('' ipdeq*2 = '',i6)') ipdeq*2
         write(6,'('' iprpts  = '',i6)') iprpts
	 call fstop(sname)
      endif
c
c     ================================================================
c     initialize quantities
c     ================================================================
      imm=0
      iter=0
      dk=kqn
      dm=h/720.0d0
      rtmp(0)=0.d0
      do j=1,iprpts
	rtmp(j)=sqrt(r(j))
      enddo
c
c     ================================================================
c     no. of nodes for the current states big component
c     ================================================================
      nodes=nqn-lqn
c
c     ================================================================
c     first power of small r expansion
c     ================================================================
      gam=sqrt(dk*dk-4.0d0*z*z/(c*c))
c
c     ================================================================
c     slope of the wavefcns at the origine
c     ================================================================
      if(nodes-2*(nodes/2).ne.0) then
         sign= 1.0d0
      else
         sign=-1.0d0
      endif
      slp=sign*kqn/abs(kqn)
c
c     ================================================================
c     find the lowest energy limit
c     ================================================================
      lll=(lqn*(lqn+1))/2
      if (lll.ne.0) then
         elim=-z*z/(0.75d0*nqn*nqn)
      else
         elim=(rv(1)+lll/r(1))/r(1)
         do j=2,nlast
            elim=min((rv(j)+lll/r(j))/r(j),elim)
         enddo
      endif
c
c     ================================================================
c     check potential
c     ================================================================
      if(elim.ge.0) then
         write(6,'('' DEEPST:: v+l*(l+1)/r**2 always positive'',f5.1)')z
         call fstop(sname)
      endif
c
      if(en.le.elim) en=elim*0.5d0
c
c     ================================================================
c     routine outws performs the outward integration
c     ================================================================
c     ----------------------------------------------------------------
  2   continue
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
c     routine inws performs the inward integration
c     ================================================================
c     ----------------------------------------------------------------
      call inws(invp,nmax,rg,rf,rv,r,en,c,drg,drf,
     >          dk,dm,nlast,imm,iprpts,ipdeq)
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
      if (fnrm.lt.0.0d0) then
         write(6,'('' DEEPST:: wrong big component change'')')
         call fstop(sname)
      endif
c
      do j=invp,nmax
         rg(j)=rg(j)*fnrm
         rf(j)=rf(j)*fnrm
      enddo
c
c     ================================================================
c     energy derivative of the wvfcns "log. derivative"
c     ================================================================
      tmp(0)=0.d0
      do j=1,nmax
         tmp(j)=(rg(j)*rg(j)+rf(j)*rf(j))
      enddo   

!     write (*,*) "about to call NEWINT ",nmax
      call newint(nmax+1,rtmp,tmp,qmp,1)
      rose=2.d0*qmp(nmax)+h*r(invp)*(rfm*rfm-rf(invp)*rf(invp))/3.0d0
c
c     ================================================================
c     energy update
c     ================================================================
      de=rg(invp)*(rfm-rf(invp))*c/rose
      imm=0
      val=abs(de/en)
      if (val.le.tol) go to 7
 5    enew=en+de

      if(enew.lt.0.0d0) go to 6
      de=de+0.5d0
      val=val*0.5d0
      if (val.gt.tol) go to 5
c
c     ================================================================
c     just in case the energy becomes zero
c     ================================================================
      do ir=1,nws
        write(6,'(''ir,vr = '',i8,f20.12)') ir, rv(ir)
      end do

      write(6,'('' DEEPST:: nws = '',i4)') nws

      write(6,'('' DEEPST:: zero energy '',f5.2)')z
      write(6,'('' DEEPST:: n = '',i2)') nqn
      write(6,'('' DEEPST:: l = '',i2)') lqn
      call fstop(sname)
c
c     ================================================================
c     not yet convergence: try again
c     ================================================================
  6   en=enew
      if (val.le.0.2d0) imm=1
      iter=iter+1
      if(iter.le.nitmax) go to 2
c
c     ================================================================
c     not converged, too small tolerance or too small nitmax
c     ================================================================
      write(6,'('' DEEPST:: warning: too many energy iterations'')')
c
  7   continue
      do j=1,nmax
         rf(j)=rf(j)*rf(j)+rg(j)*rg(j)
      enddo    
      if(nmax.lt.nlast) then
         do j=nmax+1,nlast
            rf(j)=0.0d0
         enddo   
      endif
c
c     ================================================================
c     end of job: waiting for next state
c     ================================================================
      return
      end
