      subroutine dirmago2op(socsc,e,l,my,vr,bspr,bopr,
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
c     *                                                         *
c     * strange et al., j.phys.c.solid state phys.17,3355(1984) *
c     * scaling of SOC and OP-term ala Hubert Ebert included    *
c     ***********************************************************
c
      integer iprpts

      complex*16 e,p1,q1,p2,q2,pp1,qp1,pp2,qp2
      complex*16 p1n,q1n,k11,k12,k13,k14,m11,m12,m13,m14       
      complex*16 p2n,q2n,k21,k22,k23,k24,m21,m22,m23,m24       
      complex*16 p1sum,q1sum,p10,q10,p11,q11,pp1np1,qp1np1          
      complex*16 p2sum,q2sum,p20,q20,p21,q21,pp2np1,qp2np1          
      complex*16 g,gp,f
c
      dimension vr(iprpts),bspr(iprpts),bopr(iprpts,2)
      dimension ba11(iprpts),ba12(iprpts),ba22(iprpts)
      dimension bb11(iprpts),bb22(iprpts)
      dimension p1(iprpts),p2(iprpts),q1(iprpts),q2(iprpts)
      dimension pp1(iprpts),pp2(iprpts),qp1(iprpts),qp2(iprpts)
      dimension fk1(2),fk2(2),gk1(2),gk2(2)
      dimension g(2,2,iprpts),gp(2,2),f(2,2,iprpts)
      data nitmax/100/,test/1.0d10/,imode/0/
c
      if(iabs(my).eq.2*l+1) stop 'dirmag2: my=+/-(l+1/2)'
      dkoef1=475.d0/502.d0
      dkoef2= 27.d0/502.d0
      xl = dfloat(l)
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

c
      xl=dfloat(l)
      xkap=xl
      xk1=socsc*(1.d0+xkap)-1.d0
      facl1=xl*(xl+1.d0)-xk1*(xk1+1.d0)
      u1=(xk1+dsqrt(xk1*xk1-hoc*hoc+facl1))/hoc
      xkap=-xl-1.d0
      xk2=socsc*(1.d0+xkap)-1.d0
      facl2=xl*(xl+1.d0)-xk2*(xk2+1.d0)
      u2=(xk2+dsqrt(xk2*xk2-hoc*hoc+facl2))/hoc
c
c     get combined effective fields
c
      call brmat(l,my,ns,bspr,bopr,ba11,ba12,ba22,bb11,bb22,iprpts)
c
      do ii=1,2
      if(ii.eq.2) go to 40
c
      p1(1) = 1.0d-20
      q1(1) = c*u1*1.0d-20
      p2(1) = (0.0d 00,0.0d 00)
      q2(1) = (0.0d 00,0.0d 00)
      go to 50
c
   40 continue
      p1(1) = (0.0d 00,0.0d 00)
      q1(1) = (0.0d 00,0.0d 00)
      p2(1) = 1.0d-20
      q2(1) = c*u2*1.0d-20
   50 continue
c
      r=dexp(xnot)
      call dmag2op(pp1(1),qp1(1),pp2(1),qp2(1),
     >           p1(1),q1(1),p2(1),q2(1),
     >           xk1,xk2,facl1,facl2,cin,e,r,
     >           vr(1),ba11(1),ba12(1),ba22(1),bb11(1),bb22(1))
c
c
c---> start runge-kutta procedure (points 2, ... , 6)
c
      do n=1,5
c
        x=xnot+(n-1)*dx
        rn=dexp(x)
        vrn=vr(n)
        ba11n=ba11(n)
        ba12n=ba12(n)
        ba22n=ba22(n)
        bb11n=bb11(n)
        bb22n=bb22(n)
c
        rnp1=dexp(x+dx)
        vrnp1=vr(n+1)
        ba11np1=ba11(n+1)
        ba12np1=ba12(n+1)
        ba22np1=ba22(n+1)
        bb11np1=bb11(n+1)
        bb22np1=bb22(n+1)
c
        rmid=dexp(x+dx2)
        vrmid=0.5d0*(vrn+vrnp1)
        ba11mid=0.5d0*(ba11n+ba11np1)
        ba12mid=0.5d0*(ba12n+ba12np1)
        ba22mid=0.5d0*(ba22n+ba22np1)
        bb11mid=0.5d0*(bb11n+bb11np1)
        bb22mid=0.5d0*(bb22n+bb22np1)
c
        p1n=p1(n)
        q1n=q1(n)
        p2n=p2(n)
        q2n=q2(n)
        call dmag2op(k11,m11,k21,m21,p1n,q1n,p2n,q2n,
     >             xk1,xk2,facl1,facl2,cin,e,rn,
     >             vrn,ba11n,ba12n,ba22n,bb11n,bb22n)
c
        call dmag2op(k12,m12,k22,m22,
     >             p1n+dx2*k11,q1n+dx2*m11,p2n+dx2*k21,q2n+dx2*m21,
     >             xk1,xk2,facl1,facl2,cin,e,rmid,
     >             vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)
c
        call dmag2op(k13,m13,k23,m23,
     >             p1n+dx2*k12,q1n+dx2*m12,p2n+dx2*k22,q2n+dx2*m22,
     >             xk1,xk2,facl1,facl2,cin,e,rmid,
     >             vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)
c
        call dmag2op(k14,m14,k24,m24,
     >             p1n+dx*k13,q1n+dx*m13,p2n+dx*k23,q2n+dx*m23,
     >             xk1,xk2,facl1,facl2,cin,e,rnp1,
     >             vrnp1,ba11np1,ba12np1,ba22np1,bb11np1,bb22np1)
c
        p1(n+1)=p1n+dx*(k11+2.d0*k12+2.d0*k13+k14)/6.d0
        q1(n+1)=q1n+dx*(m11+2.d0*m12+2.d0*m13+m14)/6.d0
        p2(n+1)=p2n+dx*(k21+2.d0*k22+2.d0*k23+k24)/6.d0
        q2(n+1)=q2n+dx*(m21+2.d0*m22+2.d0*m23+m24)/6.d0
        call dmag2op(pp1(n+1),qp1(n+1),pp2(n+1),qp2(n+1),
     >             p1(n+1),q1(n+1),p2(n+1),q2(n+1),
     >             xk1,xk2,facl1,facl2,cin,e,rnp1,
     >             vrnp1,ba11np1,ba12np1,ba22np1,bb11np1,bb22np1)
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
        ba11c=ba11(n+1)
        ba12c=ba12(n+1)
        ba22c=ba22(n+1)
        bb11c=bb11(n+1)
        bb22c=bb22(n+1)
c
        p1sum=646.d0*pp1(n)-264.d0*pp1(n-1)+106.d0*pp1(n-2)-
     >        19.d0*pp1(n-3)
        q1sum=646.d0*qp1(n)-264.d0*qp1(n-1)+106.d0*qp1(n-2)-
     >        19.d0*qp1(n-3)
        p2sum=646.d0*pp2(n)-264.d0*pp2(n-1)+106.d0*pp2(n-2)-
     >        19.d0*pp2(n-3)
        q2sum=646.d0*qp2(n)-264.d0*qp2(n-1)+106.d0*qp2(n-2)-
     >        19.d0*qp2(n-3)
c
c  predict for point n+1
c
        p10=p1(n)+dx*(251.d0*pp1(n-4)-1274.d0*pp1(n-3)+
     >      2616.d0*pp1(n-2)-2774.d0*pp1(n-1)+1901.d0*pp1(n))/720.d0
        q10=q1(n)+dx*(251.d0*qp1(n-4)-1274.d0*qp1(n-3)+
     >      2616.d0*qp1(n-2)-2774.d0*qp1(n-1)+1901.d0*qp1(n))/720.d0
        p20=p2(n)+dx*(251.d0*pp2(n-4)-1274.d0*pp2(n-3)+
     >      2616.d0*pp2(n-2)-2774.d0*pp2(n-1)+1901.d0*pp2(n))/720.d0
        q20=q2(n)+dx*(251.d0*qp2(n-4)-1274.d0*qp2(n-3)+
     >      2616.d0*qp2(n-2)-2774.d0*qp2(n-1)+1901.d0*qp2(n))/720.d0
c
c  correct
c
        if(imode.eq.1) then
c
        do nit=1,nitmax
c
          call dmag2op(pp1np1,qp1np1,pp2np1,qp2np1,
     >               p10,q10,p20,q20,
     >               xk1,xk2,facl1,facl2,cin,e,r,
     >               vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
          p11=p1(n)+dx*(251.d0*pp1np1+p1sum)/720.d0
          q11=q1(n)+dx*(251.d0*qp1np1+q1sum)/720.d0
          p21=p2(n)+dx*(251.d0*pp2np1+p2sum)/720.d0
          q21=q2(n)+dx*(251.d0*qp2np1+q2sum)/720.d0
c
c  compare predictor with corrector
c
          if(test*cdabs(p11-p10).gt.cdabs(p10)) goto 15
          if(test*cdabs(q11-q10).gt.cdabs(q10)) goto 15
          if(test*cdabs(p21-p20).gt.cdabs(p20)) goto 15
          if(test*cdabs(q21-q20).gt.cdabs(q20)) goto 15
          goto 10
c
   15     p10=p11
          q10=q11
          p20=p21
          q20=q21
        end do
        write(6,*) n+1,r,nit,' not converged'
c
        else
c
        call dmag2op(pp1np1,qp1np1,pp2np1,qp2np1,
     >             p10,q10,p20,q20,
     >             xk1,xk2,facl1,facl2,cin,e,r,
     >             vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
        p11=p1(n)+dx*(251.d0*pp1np1+p1sum)/720.d0
        q11=q1(n)+dx*(251.d0*qp1np1+q1sum)/720.d0
        p21=p2(n)+dx*(251.d0*pp2np1+p2sum)/720.d0
        q21=q2(n)+dx*(251.d0*qp2np1+q2sum)/720.d0
c
        q11=dkoef1*q11+dkoef2*q10
        p11=dkoef1*p11+dkoef2*p10
        q21=dkoef1*q21+dkoef2*q20
        p21=dkoef1*p21+dkoef2*p20
c
        end if
c
   10   q1(n+1)=q11
        p1(n+1)=p11
        q2(n+1)=q21
        p2(n+1)=p21
        call dmag2op(pp1(n+1),qp1(n+1),pp2(n+1),qp2(n+1),
     >             p1(n+1),q1(n+1),p2(n+1),q2(n+1),
     >             xk1,xk2,facl1,facl2,cin,e,r,
     >             vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
c
      end do
c
c     store radial amplitudes times radius
c
      do n=1,ns
        g(1,ii,n)=p1(n)
        g(2,ii,n)=p2(n)
        f(1,ii,n)=q1(n)/c
        f(2,ii,n)=q2(n)/c
      end do
      gp(1,ii)=(pp1(ns)-p1(ns))/rs
      gp(2,ii)=(pp2(ns)-p2(ns))/rs
c
      end do  
c
      return
      end
      subroutine dirmagi2op(socsc,e,l,my,vr,bspr,bopr,
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
c     *                                                         *
c     * strange et al., j.phys.c.solid state phys.17,3355(1984) *
c     * scaling of SOC and OP-term ala Hubert Ebert included    *
c     ***********************************************************
c
      integer iprpts

      complex*16 e,psq,pe
      complex*16 p1,q1,p2,q2,pp1,qp1,pp2,qp2
      complex*16 p1n,q1n,k11,k12,k13,k14,m11,m12,m13,m14       
      complex*16 p2n,q2n,k21,k22,k23,k24,m21,m22,m23,m24       
      complex*16 p1sum,q1sum,p10,q10,p11,q11,pp1nm1,qp1nm1          
      complex*16 p2sum,q2sum,p20,q20,p21,q21,pp2nm1,qp2nm1          
      complex*16 g,gp,f
      complex*16 fb(0:5),fn(0:5),fh(0:5)
c
      dimension vr(iprpts),bspr(iprpts),bopr(iprpts,2)
      dimension xlag(iprpts),vrlag(iprpts)
      dimension bsprlag(iprpts),boprlag(iprpts,2)
      dimension ba11(iprpts),ba12(iprpts)
      dimension ba22(iprpts),bb11(iprpts),bb22(iprpts)
      dimension p1(iprpts),p2(iprpts),q1(iprpts),q2(iprpts)
      dimension pp1(iprpts),pp2(iprpts),qp1(iprpts),qp2(iprpts)
      dimension fk1(2),fk2(2),gk1(2),gk2(2)
      dimension g(2,2,iprpts),gp(2,2),f(2,2,iprpts)
      data nitmax/100/,test/1.0d10/,imode/0/,nout/5/
c
      if(iabs(my).eq.2*l+1) stop 'dirmag2: my=+/-(l+1/2)'
      dkoef1=475.d0/502.d0
      dkoef2= 27.d0/502.d0
c
      lb1=l-1
      lb2=l+1
      sk1=l-lb1
      sk2=l-lb2
c
      xl = dfloat(l)
      xkap=xl
      xk1=socsc*(1.d0+xkap)-1.d0
      facl1=xl*(xl+1.d0)-xk1*(xk1+1.d0)
      xkap=-xl-1.d0
      xk2=socsc*(1.d0+xkap)-1.d0
      facl2=xl*(xl+1.d0)-xk2*(xk2+1.d0)
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
      do ii=1,2
c
c---> set up the starting values outside the muffin tin
c     corresponding to the boundary condition
c
      x=xs+nout*dx
      r=dexp(x)
      call csbf(l+1,pe,r,fb,fn,fh)
      n=ns+nout
c
      if(ii.eq.2) go to 40
      p1(n) = r*fb(l)
      q1(n) = ( (xk1-xl)*fb(l) + sk1*pe*r*fb(lb1) ) * e/psq
      p2(n) = 0.d0             
      q2(n) = 0.d0                     
      call dmag2op(pp1(n),qp1(n),pp2(n),qp2(n),
     >           p1(n),q1(n),p2(n),q2(n),
     >           xk1,xk2,facl1,facl2,cin,e,r,
     >           vrlag(nlag),ba11(nlag),ba12(nlag),ba22(nlag),
     >                       bb11(nlag),bb22(nlag))
      go to 50
   40 continue
      p1(n) = 0.d0    
      q1(n) = 0.d0                    
      p2(n) = r*fb(l)          
      q2(n) = ( (xk2+xl+1.d0)*fb(l) + sk2*pe*r*fb(lb2) ) * e/psq 
      call dmag2op(pp1(n),qp1(n),pp2(n),qp2(n),
     >           p1(n),q1(n),p2(n),q2(n),
     >           xk1,xk2,facl1,facl2,cin,e,r,
     >           vrlag(nlag),ba11(nlag),ba12(nlag),ba22(nlag),
     >                       bb11(nlag),bb22(nlag))
   50 continue
c
c     write(6,*) n
c     write(6,'(8d15.7)') p1(n),q1(n),p2(n),q2(n)
c     write(6,'(8d15.7)') pp1(n),qp1(n),pp2(n),qp2(n)
c     write(6,*)
c
c---> start runge-kutta procedure (points ns+nout-1, ... , ns+nout-5)
c
      do n=ns+nout,ns+nout-4,-1
c
        x=xs+(n-ns)*dx
        rn=dexp(x)
        vrn=ylag(x,xlag,vrlag,0,ilag,nlag,iex)
        ba11n=ylag(x,xlag,ba11,0,ilag,nlag,iex)
        ba12n=ylag(x,xlag,ba12,0,ilag,nlag,iex)
        ba22n=ylag(x,xlag,ba22,0,ilag,nlag,iex)
        bb11n=ylag(x,xlag,bb11,0,ilag,nlag,iex)
        bb22n=ylag(x,xlag,bb22,0,ilag,nlag,iex)
c
        rmid=dexp(x-dx2)
        vrmid=ylag(x-dx2,xlag,vrlag,0,ilag,nlag,iex)
        ba11mid=ylag(x-dx2,xlag,ba11,0,ilag,nlag,iex)
        ba12mid=ylag(x-dx2,xlag,ba12,0,ilag,nlag,iex)
        ba22mid=ylag(x-dx2,xlag,ba22,0,ilag,nlag,iex)
        bb11mid=ylag(x-dx2,xlag,bb11,0,ilag,nlag,iex)
        bb22mid=ylag(x-dx2,xlag,bb22,0,ilag,nlag,iex)
c
        rnm1=dexp(x-dx)
        vrnm1=ylag(x-dx,xlag,vrlag,0,ilag,nlag,iex)
        ba11nm1=ylag(x-dx,xlag,ba11,0,ilag,nlag,iex)
        ba12nm1=ylag(x-dx,xlag,ba12,0,ilag,nlag,iex)
        ba22nm1=ylag(x-dx,xlag,ba22,0,ilag,nlag,iex)
        bb11nm1=ylag(x-dx,xlag,bb11,0,ilag,nlag,iex)
        bb22nm1=ylag(x-dx,xlag,bb22,0,ilag,nlag,iex)
c
        p1n=p1(n)
        q1n=q1(n)
        p2n=p2(n)
        q2n=q2(n)
        call dmag2op(k11,m11,k21,m21,p1n,q1n,p2n,q2n,
     >             xk1,xk2,facl1,facl2,cin,e,rn,
     >             vrn,ba11n,ba12n,ba22n,bb11n,bb22n)
c
        call dmag2op(k12,m12,k22,m22,
     >             p1n-dx2*k11,q1n-dx2*m11,p2n-dx2*k21,q2n-dx2*m21,
     >             xk1,xk2,facl1,facl2,cin,e,rmid,
     >             vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)
c
        call dmag2op(k13,m13,k23,m23,
     >             p1n-dx2*k12,q1n-dx2*m12,p2n-dx2*k22,q2n-dx2*m22,
     >             xk1,xk2,facl1,facl2,cin,e,rmid,
     >             vrmid,ba11mid,ba12mid,ba22mid,bb11mid,bb22mid)
c
        call dmag2op(k14,m14,k24,m24,
     >             p1n-dx*k13,q1n-dx*m13,p2n-dx*k23,q2n-dx*m23,
     >             xk1,xk2,facl1,facl2,cin,e,rnm1,
     >             vrnm1,ba11nm1,ba12nm1,ba22nm1,bb11nm1,bb22nm1)
c
        p1(n-1)=p1n-dx*(k11+2.d0*k12+2.d0*k13+k14)/6.d0
        q1(n-1)=q1n-dx*(m11+2.d0*m12+2.d0*m13+m14)/6.d0
        p2(n-1)=p2n-dx*(k21+2.d0*k22+2.d0*k23+k24)/6.d0
        q2(n-1)=q2n-dx*(m21+2.d0*m22+2.d0*m23+m24)/6.d0
        call dmag2op(pp1(n-1),qp1(n-1),pp2(n-1),qp2(n-1),
     >             p1(n-1),q1(n-1),p2(n-1),q2(n-1),
     >             xk1,xk2,facl1,facl2,cin,e,rnm1,
     >             vrnm1,ba11nm1,ba12nm1,ba22nm1,bb11nm1,bb22nm1)
c
      end do
c
c
c---> begin adams procedure (points ns+nout-6,  ... , 1)
c
      do n=ns+nout-5,2,-1
c
        x=xs+(n-ns)*dx
        r=dexp(x-dx)
        vrc=vr(n-1)
        ba11c=ba11(n-1)
        ba12c=ba12(n-1)
        ba22c=ba22(n-1)
        bb11c=bb11(n-1)
        bb22c=bb22(n-1)
c
        p1sum=646.d0*pp1(n)-264.d0*pp1(n+1)+106.d0*pp1(n+2)-
     >        19.d0*pp1(n+3)
        q1sum=646.d0*qp1(n)-264.d0*qp1(n+1)+106.d0*qp1(n+2)-
     >        19.d0*qp1(n+3)
        p2sum=646.d0*pp2(n)-264.d0*pp2(n+1)+106.d0*pp2(n+2)-
     >        19.d0*pp2(n+3)
        q2sum=646.d0*qp2(n)-264.d0*qp2(n+1)+106.d0*qp2(n+2)-
     >        19.d0*qp2(n+3)
c
c  predict for point n-1
c
        p10=p1(n)-dx*(251.d0*pp1(n+4)-1274.d0*pp1(n+3)+
     >      2616.d0*pp1(n+2)-2774.d0*pp1(n+1)+1901.d0*pp1(n))/720.d0
        q10=q1(n)-dx*(251.d0*qp1(n+4)-1274.d0*qp1(n+3)+
     >      2616.d0*qp1(n+2)-2774.d0*qp1(n+1)+1901.d0*qp1(n))/720.d0
        p20=p2(n)-dx*(251.d0*pp2(n+4)-1274.d0*pp2(n+3)+
     >      2616.d0*pp2(n+2)-2774.d0*pp2(n+1)+1901.d0*pp2(n))/720.d0
        q20=q2(n)-dx*(251.d0*qp2(n+4)-1274.d0*qp2(n+3)+
     >      2616.d0*qp2(n+2)-2774.d0*qp2(n+1)+1901.d0*qp2(n))/720.d0
c
c  correct
c
        if(imode.eq.1) then
c
        do nit=1,nitmax
c
          call dmag2op(pp1nm1,qp1nm1,pp2nm1,qp2nm1,
     >               p10,q10,p20,q20,
     >               xk1,xk2,facl1,facl2,cin,e,r,
     >               vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
          p11=p1(n)-dx*(251.d0*pp1nm1+p1sum)/720.d0
          q11=q1(n)-dx*(251.d0*qp1nm1+q1sum)/720.d0
          p21=p2(n)-dx*(251.d0*pp2nm1+p2sum)/720.d0
          q21=q2(n)-dx*(251.d0*qp2nm1+q2sum)/720.d0
c
c  compare predictor with corrector
c
          if(test*cdabs(p11-p10).gt.cdabs(p10)) goto 15
          if(test*cdabs(q11-q10).gt.cdabs(q10)) goto 15
          if(test*cdabs(p21-p20).gt.cdabs(p20)) goto 15
          if(test*cdabs(q21-q20).gt.cdabs(q20)) goto 15
          goto 10
c
   15     p10=p11
          q10=q11
          p20=p21
          q20=q21
c
        end do
        write(6,*) n-1,r,nit,' not converged'
c
        else
c
        call dmag2op(pp1nm1,qp1nm1,pp2nm1,qp2nm1,
     >             p10,q10,p20,q20,
     >             xk1,xk2,facl1,facl2,cin,e,r,
     >             vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
        p11=p1(n)-dx*(251.d0*pp1nm1+p1sum)/720.d0
        q11=q1(n)-dx*(251.d0*qp1nm1+q1sum)/720.d0
        p21=p2(n)-dx*(251.d0*pp2nm1+p2sum)/720.d0
        q21=q2(n)-dx*(251.d0*qp2nm1+q2sum)/720.d0
c
        q11=dkoef1*q11+dkoef2*q10
        p11=dkoef1*p11+dkoef2*p10
        q21=dkoef1*q21+dkoef2*q20
        p21=dkoef1*p21+dkoef2*p20
c
        end if
c
   10   q1(n-1)=q11
        p1(n-1)=p11
        q2(n-1)=q21
        p2(n-1)=p21
        call dmag2op(pp1(n-1),qp1(n-1),pp2(n-1),qp2(n-1),
     >             p1(n-1),q1(n-1),p2(n-1),q2(n-1),
     >             xk1,xk2,facl1,facl2,cin,e,r,
     >             vrc,ba11c,ba12c,ba22c,bb11c,bb22c)
c
      end do
c
c     store radial amplitudes times radius
c
      do n=1,ns+nout
        g(1,ii,n)=p1(n)
        g(2,ii,n)=p2(n)
        f(1,ii,n)=q1(n)/c
        f(2,ii,n)=q2(n)/c
      end do
      gp(1,ii)=(pp1(ns)-p1(ns))/rs
      gp(2,ii)=(pp2(ns)-p2(ns))/rs
c
      end do  
c
      return
      end
      subroutine dmag2op(pp1,qp1,pp2,qp2,p1,q1,p2,q2,
     >                 xk1,xk2,facl1,facl2,cin,e,r,
     >                 vr,ba11,ba12,ba22,bb11,bb22)
c
      implicit real*8 (a-h,o-z)
      complex*16 pp1,qp1,pp2,qp2,p1,q1,p2,q2,e,tr,sr1,sr2
c
      tr=r*e-vr
      sr1=cin*tr+r+cin*bb11     
      sr2=cin*tr+r+cin*bb22     
c
      qp1= xk1*q1 - tr*p1 + ba11*p1 + ba12*p2 + facl1*p1/sr1 
c
      pp1=-xk1*p1 + sr1*q1
c
      qp2= xk2*q2 - tr*p2 + ba22*p2 + ba12*p1 + facl2*p2/sr2 
c
      pp2=-xk2*p2 + sr2*q2
c
      return
      end
