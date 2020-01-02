      subroutine dirmago1op(socsc,e,l,my,vr,bspr,bopr,
     >                    dx,xnot,rs,ns,g,f,gp,iprpts)
c
      implicit real*8(a-h,o-z)
!     include '../param.h'
!     include 'atom_param.h'
c
c     ***********************************************************
c     * integration of relativistic radial dirac equation       *
c     * in the presence of an internal field by adams method    *
c     * integrate outward!                                      *
c     * strange et al., j.phys.c.solid state phys.17,3355(1984) *
c     * scaling of SOC and OP-term ala Hubert Ebert included    *
c     *                                                         *
c     * my=+/-(l+1/2), kap=-l-1  case!!!                        *
c     ***********************************************************
c
      integer iprpts
      complex*16 e,p,q,pp,qp
      complex*16 pn,qn,k1,k2,k3,k4,m1,m2,m3,m4
      complex*16 psum,qsum,p0,q0,p1,q1,ppnp1,qpnp1
      complex*16 g,gp,f
c
      dimension vr(iprpts),bspr(iprpts),bopr(iprpts,2)
      dimension ba11(iprpts),ba12(iprpts)
      dimension ba22(iprpts),bb11(iprpts),bb22(iprpts)
      dimension p(iprpts),q(iprpts),pp(iprpts),qp(iprpts)
      dimension g(iprpts),f(iprpts)
c
      data nitmax/100/,test/1.0d10/,imode/0/
c
      if(iabs(my).lt.2*l+1) stop 'dirmag1: -l-1/2 < my < l+1/2'
      dkoef1=475.d0/502.d0
      dkoef2= 27.d0/502.d0
c
      xl=dfloat(l)
      xkap=-xl-1.d0
      xk=socsc*(1.d0+xkap)-1.d0
      facl=xl*(xl+1.d0)-xk*(xk+1.d0)
c
      dx2=0.5d0*dx
      c = 274.072d0
      cin = 1.0d 00/(c*c)
c
c     calculate boundary conditions
c
      if(vr(1).lt.-1.d-3)then
      hoc=-vr(1)/c
      else
c     write(6,'(99d12.4)')(vr(i),i=1,5)
      hoc=-vr(2)/c
      vr(1)=vr(2)
      endif
      u=(xk+dsqrt(xk*xk-hoc*hoc+facl))/hoc
      p(1) = 1.0d-20
      q(1) = c*u*1.0d-20
c
c     get combined effective fields
c
      call brmat(l,my,ns,bspr,bopr,ba11,ba12,ba22,bb11,bb22,iprpts)
c
      r=dexp(xnot)
      call dmag1op(pp(1),qp(1),p(1),q(1),
     >           xk,facl,cin,e,r,vr(1),ba22(1),bb22(1))
c
c---> start runge-kutta procedure (points 2, ... , 6)
c
      do n=1,5
c
        x=xnot+(n-1)*dx
c
        rn=dexp(x)
        vrn=vr(n)
        ba22n=ba22(n)
        bb22n=bb22(n)
c
        rnp1=dexp(x+dx)
        vrnp1=vr(n+1)
        ba22np1=ba22(n+1)
        bb22np1=bb22(n+1)
c
        rmid=dexp(x+dx2)
        vrmid=0.5d0*(vrn+vrnp1)
        ba22mid=0.5d0*(ba22n+ba22np1)
        bb22mid=0.5d0*(bb22n+ba22np1)
c
        pn=p(n)
        qn=q(n)
        call dmag1op(k1,m1,pn,qn,
     >             xk,facl,cin,e,rn,vrn,ba22n,bb22n)
c
        call dmag1op(k2,m2,pn+dx2*k1,qn+dx2*m1,
     >             xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)
c
        call dmag1op(k3,m3,pn+dx2*k2,qn+dx2*m2,
     >             xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)
c
        call dmag1op(k4,m4,pn+dx*k3,qn+dx*m3,
     >             xk,facl,cin,e,rnp1,vrnp1,ba22np1,bb22np1)
c
        p(n+1)=pn+dx*(k1+2.d0*k2+2.d0*k3+k4)/6.d0
        q(n+1)=qn+dx*(m1+2.d0*m2+2.d0*m3+m4)/6.d0
        call dmag1op(pp(n+1),qp(n+1),p(n+1),q(n+1),    
     >             xk,facl,cin,e,rnp1,vrnp1,ba22np1,bb22np1)
c
      end do                 
c
c---> begin adams procedure (points 7, 8, ... , ns)
c
      do n=6,ns-1
c
        x=xnot+(n-1)*dx
        r=dexp(x+dx)
        vrc=vr(n+1)
        ba22c=ba22(n+1)
        bb22c=bb22(n+1)
c
        psum=646.d0*pp(n)-264.d0*pp(n-1)+106.d0*pp(n-2)-19.d0*pp(n-3)
        qsum=646.d0*qp(n)-264.d0*qp(n-1)+106.d0*qp(n-2)-19.d0*qp(n-3)
c
c  predict for point n+1
c
        p0=p(n)+dx*(251.d0*pp(n-4)-1274.d0*pp(n-3)+2616.d0*pp(n-2)-
     >              2774.d0*pp(n-1)+1901.d0*pp(n))/720.d0
        q0=q(n)+dx*(251.d0*qp(n-4)-1274.d0*qp(n-3)+2616.d0*qp(n-2)-
     >              2774.d0*qp(n-1)+1901.d0*qp(n))/720.d0
c
c  correct
c
        if(imode.eq.1) then
c
        do nit=1,nitmax
c
          call dmag1op(ppnp1,qpnp1,p0,q0,
     >               xk,facl,cin,e,r,vrc,ba22c,bb22c)
          p1=p(n)+dx*(251.d0*ppnp1+psum)/720.d0
          q1=q(n)+dx*(251.d0*qpnp1+qsum)/720.d0
c
c  compare predictor with corrector
c
          if(test*cdabs(p1-p0).gt.cdabs(p0)) goto 15
          if(test*cdabs(q1-q0).gt.cdabs(q0)) goto 15 
          goto 10
c
   15     p0=p1   
          q0=q1   
c
        end do  
        write(6,*) n+1,r,nit,' not converged'
c
        else
c
        call dmag1op(ppnp1,qpnp1,p0,q0,
     >             xk,facl,cin,e,r,vrc,ba22c,bb22c)
        p1=p(n)+dx*(251.d0*ppnp1+psum)/720.d0
        q1=q(n)+dx*(251.d0*qpnp1+qsum)/720.d0
        q1=dkoef1*q1+dkoef2*q0
        p1=dkoef1*p1+dkoef2*p0
c
        end if
c
   10   q(n+1)=q1  
        p(n+1)=p1
        call dmag1op(pp(n+1),qp(n+1),p(n+1),q(n+1),
     >             xk,facl,cin,e,r,vrc,ba22c,bb22c)
c
      end do
c
c     store radial amplitudes times radius
c
      do n=1,ns
        g(n)=p(n)
        f(n)=q(n)/c
      end do   
      gp=(pp(ns)-p(ns))/rs
c
      return
      end
      subroutine dirmagi1op(socsc,e,l,my,vr,bspr,bopr,
     >                    dx,xnot,rs,ns,g,f,gp,iprpts)
c
      implicit real*8(a-h,o-z)
!     include '../param.h'
!     include 'atom_param.h'
c
c     ***********************************************************
c     * integration of relativistic radial dirac equation       *
c     * in the presence of an internal field by adams method    *
c     * integrate inward!                                       *
c     * strange et al., j.phys.c.solid state phys.17,3355(1984) *
c     * scaling of SOC and OP-term ala Hubert Ebert included    *
c     *                                                         *
c     * my=+/-(l+1/2), kap=-l-1  case!!!                        *
c     ***********************************************************
c
      integer iprpts
      complex*16 e,psq,pe
      complex*16 p,q,pp,qp
      complex*16 pn,qn,k1,k2,k3,k4,m1,m2,m3,m4
      complex*16 psum,qsum,p0,q0,p1,q1,ppnm1,qpnm1
      complex*16 g,gp,f
      complex*16 fb(0:5),fn(0:5),fh(0:5)
c
      dimension vr(iprpts),bspr(iprpts),bopr(iprpts,2)
      dimension xlag(iprpts),vrlag(iprpts)
      dimension bsprlag(iprpts),boprlag(iprpts,2)
      dimension ba11(iprpts),ba12(iprpts)
      dimension ba22(iprpts),bb11(iprpts),bb22(iprpts)
      dimension p(iprpts),q(iprpts),pp(iprpts),qp(iprpts)
      dimension g(iprpts),f(iprpts)
c
      data nitmax/100/,test/1.0d10/,imode/0/,nout/5/
c
      if(iabs(my).lt.2*l+1) stop 'dirmag1: -l-1/2 < my < l+1/2'
      dkoef1=475.d0/502.d0
      dkoef2= 27.d0/502.d0
c
      xl=dfloat(l)
      xkap=-xl-1.d0
      xk=socsc*(1.d0+xkap)-1.d0
      facl=xl*(xl+1.d0)-xk*(xk+1.d0)
      lb=l+1
      sk=l-lb
c
      dx2=0.5d0*dx
      c = 274.072d0
      cin = 1.0d 00/(c*c)
      ilag=3
      xs=dlog(rs)
      do i=1,ns
        xlag(i)=xs+(i-ns)*dx
        vrlag(i)=vr(i)
        bsprlag(i)=bspr(i)
        boprlag(i,1)=bopr(i,1)
        boprlag(i,2)=bopr(i,2)
      end do
      nlag=ns+1
      vrlag(nlag)=0.d0
      bsprlag(nlag)=0.d0
      boprlag(nlag,1)=0.d0
      boprlag(nlag,2)=0.d0
      xlag(nlag)=xs+nout*dx
c
c     get combined effective fields
c
      call brmat(l,my,nlag,bsprlag,boprlag,ba11,ba12,ba22,bb11,bb22,
     >           iprpts)
c
      psq=e+e*e*cin
      pe=cdsqrt(psq)
c
c---> set up the starting values outside the muffin tin
c     corresponding to the boundary condition
c
      x=xs+nout*dx
      r=dexp(x)
      call csbf(l+1,pe,r,fb,fn,fh)
      n=ns+nout
      p(n)= r*fb(l)
      q(n)= ( (xk-xkap)*fb(l) + sk*pe*r*fb(lb) ) * e/psq
      call dmag1op(pp(n),qp(n),p(n),q(n),
     >           xk,facl,cin,e,r,vrlag(nlag),ba22(nlag),bb22(nlag))
c
c---> start runge-kutta procedure (points ns+nout-1, ... , ns+nout-5)
c
      do n=ns+nout,ns+nout-4,-1
c
        x=xs+(n-ns)*dx
        rn=dexp(x)
        vrn=ylag(x,xlag,vrlag,0,ilag,nlag,iex)
        ba22n=ylag(x,xlag,ba22,0,ilag,nlag,iex)
        bb22n=ylag(x,xlag,bb22,0,ilag,nlag,iex)
c
        rmid=dexp(x-dx2)
        vrmid=ylag(x-dx2,xlag,vrlag,0,ilag,nlag,iex)
        ba22mid=ylag(x-dx2,xlag,ba22,0,ilag,nlag,iex)
        bb22mid=ylag(x-dx2,xlag,bb22,0,ilag,nlag,iex)
c
        rnm1=dexp(x-dx)
        vrnm1=ylag(x-dx,xlag,vrlag,0,ilag,nlag,iex)
        ba22nm1=ylag(x-dx,xlag,ba22,0,ilag,nlag,iex)
        bb22nm1=ylag(x-dx,xlag,bb22,0,ilag,nlag,iex)
c
        pn=p(n)
        qn=q(n)
        call dmag1op(k1,m1,pn,qn,
     >             xk,facl,cin,e,rn,vrn,ba22n,bb22n)
c
        call dmag1op(k2,m2,pn-dx2*k1,qn-dx2*m1,
     >             xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)
c
        call dmag1op(k3,m3,pn-dx2*k2,qn-dx2*m2,
     >             xk,facl,cin,e,rmid,vrmid,ba22mid,bb22mid)
c
        call dmag1op(k4,m4,pn-dx*k3,qn-dx*m3,
     >             xk,facl,cin,e,rnm1,vrnm1,ba22nm1,bb22nm1)
c
        p(n-1)=pn-dx*(k1+2.d0*k2+2.d0*k3+k4)/6.d0
        q(n-1)=qn-dx*(m1+2.d0*m2+2.d0*m3+m4)/6.d0
        call dmag1op(pp(n-1),qp(n-1),p(n-1),q(n-1),
     >             xk,facl,cin,e,rnm1,vrnm1,ba22nm1,bb22nm1)
c
      end do
c
c---> begin adams procedure (points ns+nout-6,  ... , 1)
c
      do n=ns+nout-5,2,-1
c
        x=xs+(n-ns)*dx
        r=dexp(x-dx)
        vrc=vr(n-1)
        ba22c=ba22(n-1)
        bb22c=bb22(n-1)
c
        psum=646.d0*pp(n)-264.d0*pp(n+1)+106.d0*pp(n+2)-19.d0*pp(n+3)
        qsum=646.d0*qp(n)-264.d0*qp(n+1)+106.d0*qp(n+2)-19.d0*qp(n+3)
c
c  predict for point n-1
c
        p0=p(n)-dx*(251.d0*pp(n+4)-1274.d0*pp(n+3)+2616.d0*pp(n+2)-
     >              2774.d0*pp(n+1)+1901.d0*pp(n))/720.d0
        q0=q(n)-dx*(251.d0*qp(n+4)-1274.d0*qp(n+3)+2616.d0*qp(n+2)-
     >              2774.d0*qp(n+1)+1901.d0*qp(n))/720.d0
c
c  correct
c
        if(imode.eq.1) then
c
        do nit=1,nitmax
c
          call dmag1op(ppnm1,qpnm1,p0,q0,
     >               xk,facl,cin,e,r,vrc,ba22c,bb22c)
          p1=p(n)-dx*(251.d0*ppnm1+psum)/720.d0
          q1=q(n)-dx*(251.d0*qpnm1+qsum)/720.d0
c
c  compare predictor with corrector
c
          if(test*cdabs(p1-p0).gt.cdabs(p0)) goto 15
          if(test*cdabs(q1-q0).gt.cdabs(q0)) goto 15 
          goto 10
c
   15     p0=p1   
          q0=q1   
        end do  
        write(6,*) n-1,r,nit,' not converged'
c
        else
c
        call dmag1op(ppnm1,qpnm1,p0,q0,
     >             xk,facl,cin,e,r,vrc,ba22c,bb22c)
        p1=p(n)-dx*(251.d0*ppnm1+psum)/720.d0
        q1=q(n)-dx*(251.d0*qpnm1+qsum)/720.d0
        q1=dkoef1*q1+dkoef2*q0
        p1=dkoef1*p1+dkoef2*p0

        end if
c
   10   q(n-1)=q1  
        p(n-1)=p1
        call dmag1op(pp(n-1),qp(n-1),p(n-1),q(n-1),
     >             xk,facl,cin,e,r,vrc,ba22c,bb22c)
c
      end do
c
c     store radial amplitudes times radius
c
      do n=1,ns
        g(n)=p(n)
        f(n)=q(n)/c
      end do   
      gp=(pp(ns)-p(ns))/rs
c
      return
      end
      subroutine dmag1op(pp,qp,p,q,xk,facl,cin,e,r,vr,ba22,bb22)
c
      implicit real*8 (a-h,o-z)
      complex*16 pp,qp,p,q,e,tr,sr
c
      tr=r*e-vr
      sr=cin*tr+r+cin*bb22
c
      qp= xk*q - tr*p + ba22*p + facl*p/sr
      pp=-xk*p + sr*q
c
      return
      end
