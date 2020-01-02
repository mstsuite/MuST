c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine inws(invp,nmax,rg,rf,rv,r,en,c,drg,drf,dk,dm,nws,imm,
     >                iprpts,ipdeq)
c     ================================================================
c
      implicit   none
c
!     include    'atom_param.h'
c
      integer    iprpts,ipdeq
      integer    invp,nmax,nws,imm
      integer    i,j,itop,l
c
      real*8     en,c,dk,dm
      real*8     r(iprpts),rv(iprpts),rg(iprpts),rf(iprpts)
      real*8     drf(ipdeq*2),drg(ipdeq*2)
      real*8     p,pq,dpr,dqr,er,emvoc,dq,dp
      real*8     zero,one,two,three,four
      real*8     ten,twtnth,half,thrqrt
      real*8     twntsvn,fhnd,fhd,sixhnd
      real*8     coef1,coef2
c
      parameter  (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter  (ten=10.d0,twtnth=two/ten,half=one/two)
      parameter  (thrqrt=three/four)
      parameter  (twntsvn=27.d0,fhnd=475.d0,fhd=502.d0,sixhnd=600.d0)
      parameter  (coef1=.9462151394422310d+00)
      parameter  (coef2=.0537848605577689d+00)
c
c     ****************************************************************
c     adams 5 points  diff. eq. solver for dirac equations
c
c     drg,drf enrivatives of dp and dq times r;  c speed of light
c     rv potential times r; r: radial grid; en energy guess
c     coef1=475./502., coef2=27./502.
c     ****************************************************************
c
!     write(*,*) "Entering INWS"

c     ================================================================
c     at first find the last point where the wavefcns are not zero
c     ================================================================
      if( imm .ne. 1 ) then
         do j=1,nws,3
            nmax=nws+1-j
            if( (rv(nmax)-en*r(nmax))*r(nmax) .le. sixhnd ) go to 20
         enddo  
      endif
c
c     ================================================================
c     initial values for inward integration
c     ================================================================
 20   continue
      p=sqrt(-en*(one+en/(c*c)))
      pq=-p/(c+en/c)
      do i=1,ipdeq
         j=nmax+1-i
         rg(j)=exp(-p*r(j))
         drg(i)=-p*rg(j)*r(j)
         rf(j)=pq*rg(j)
         drf(i)=pq*drg(i)
      enddo
      itop=nmax-ipdeq
c
c     ================================================================
c     solve dirac equations now
c     ================================================================
      do i=invp,itop
         j=itop-i+invp
c
c        =============================================================
c        5 points predictor
c        =============================================================
         dpr= rg(j+1) - dm*(2.51d+02*drg(1)-1.274d+03*drg(2)+
     >        2.616d+03*drg(3)-2.774d+03*drg(4)+1.901d+03*drg(5) )
         dqr= rf(j+1) - dm*(2.51d+02*drf(1)-1.274d+03*drf(2)+
     >        2.616d+03*drf(3)-2.774d+03*drf(4)+1.901d+03*drf(5) )
c
c        =============================================================
c        shift derivatives
c        =============================================================
         do l=2,ipdeq
            drg(l-1)=drg(l)
            drf(l-1)=drf(l)
         enddo
c
c        =============================================================
c        dirac equations (log. mesh)
c        =============================================================
         er=en*r(j)
         emvoc=(er-rv(j))/c
         drg(ipdeq)=-dk*dpr+(c*r(j)+emvoc)*dqr
         drf(ipdeq)=dk*dqr-emvoc*dpr
c
c        =============================================================
c        5 points corrector
c        =============================================================
         dp= rg(j+1) - dm*(-1.9d+01*drg(1)+1.06d+02*drg(2)
     >              -2.64d+02*drg(3)+ 6.46d+02*drg(4)+2.51d+02*drg(5))
         dq= rf(j+1) - dm*(-1.9d+01*drf(1)+1.06d+02*drf(2)
     >              -2.64d+02*drf(3)+ 6.46d+02*drf(4)+2.51d+02*drf(5))
c
c        =============================================================
c        mixing
c        =============================================================
         dp=coef1*dp+coef2*dpr
         dq=coef1*dq+coef2*dqr
         rg(j)=dp
         rf(j)=dq
c
c        =============================================================
c        update derivative
c        =============================================================
         drg(ipdeq)=-dk*dp+(c*r(j)+emvoc)*dq
         drf(ipdeq)=dk*dq-emvoc*dp
      enddo
c
      return
      end
