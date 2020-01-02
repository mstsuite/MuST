c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outws(invp,rg,rf,rv,r,en,c,drg,drf,elim,z,
     >                 gam,slp,tol,imm,lll,dk,dm,nodes,nws,
     >                 iprpts,ipdeq)
c     ================================================================
c
      implicit   none
c
c     include    'atom_param.h'
c
      character  sname*32
c
      integer    iprpts,ipdeq
      integer    invp
      integer    imm
      integer    lll
      integer    nodes
      integer    nws
      integer    nd
      integer    k
      integer    l
      integer    j
c
      real*8     en
      real*8     c
      real*8     elim
      real*8     z
      real*8     gam
      real*8     slp
      real*8     tol
      real*8     dk
      real*8     dm
      real*8     r(iprpts)
      real*8     rv(iprpts)
      real*8     rg(iprpts)
      real*8     rf(iprpts)
      real*8     drf(ipdeq*2)
      real*8     drg(ipdeq*2)
      real*8     vor
      real*8     rpow
      real*8     emvoc
      real*8     dp
      real*8     dq
      real*8     dqr
      real*8     er
      real*8     dpr
      real*8     zero
      real*8     one
      real*8     twnsvn
      real*8     fhnd
      real*8     fnd
      real*8     coef1
      real*8     coef2
c
      parameter  (zero=0.d0)
      parameter  (one=1.d0)
      parameter  (twnsvn=27.d0)
      parameter  (fhnd=475.d0)
      parameter  (fnd=502.d0)
      parameter  (coef1=.9462151394422310d+00)
      parameter  (coef2=.0537848605577689d+00)
      parameter  (sname='outws')
c
c     ****************************************************************
c     outward solution
c     adams 5 points  diff. eq. solver for dirac equations
c     drg,drf derivatives of dp and dq times r;  c speed of light
c     rv potential times r; r: radial grid; en energy guess
c     coef1=475./502., coef2=27./502.
c     ****************************************************************

!     write(*,*) "Entering OUTWS"
c
      nd=0
      vor=rv(1)/r(1)
c
c     ================================================================
c     find classical inversion point
c     ================================================================
 10   if( imm .eq. 1 ) go to 40
c
      do j=ipdeq+2,nws,2
         invp=nws+1-j
         if( ((rv(invp)+lll/r(invp))/r(invp)-en) .le. zero ) go to 30
      enddo
c
 30   if( invp .gt. ipdeq ) go to 40
c
      en=en*0.5d+00
c
      if( en .lt. -tol .and. nd .le. nodes ) go to 10
c
c     ================================================================
      write(6,'('' OUTWS:: screwed potential'')')
      call fstop(sname)
c 
c     ================================================================
c     initial values for outward integration
c     ================================================================
c     ----------------------------------------------------------------
 40   continue
      call invals(rg,rf,r,slp,gam,vor,z,tol,drg,drf,en,c,dk,
     >            iprpts,ipdeq)
c     ----------------------------------------------------------------
c
      nd=1
      do j=1,ipdeq
         rpow=r(j)**gam
         if( j .ne. 1 ) then
            if( rg(j-1) .ne. zero ) then
               if( (rg(j)/rg(j-1)) .le. zero ) then
                  nd=nd+1
               endif
	    endif
	 endif
         rg(j)=rg(j)*rpow
         rf(j)=rf(j)*rpow
         drg(j)=drg(j)*rpow
         drf(j)=drf(j)*rpow
      enddo   
c
c     ================================================================
c     check consistence of signs of big and small component
c     ================================================================
      k=-1+2*(nodes-2*(nodes/2))
c
      if( (rg(1)*k) .gt. zero ) go to 80
c
 70   continue
c     ================================================================
      write(6,'('' OUTWS:: errors in small r expansion'')')
      write(6,*) '    rg(1), rf(1), r(1), rv(1) =',
     >           rg(1),rf(1),r(1),rv(1)
      call fstop(sname)
c
 80   if( (k*dk*rf(1)) .lt. zero ) go to 70
c
c     ================================================================
c     solve dirac eqs. now
c     ================================================================
      do j=ipdeq+1,invp
c
c        =============================================================
c        5 points predictor
c        =============================================================
         dpr= rg(j-1) + dm*(2.51d+02*drg(1)-1.274d+03*drg(2)+
     >        2.616d+03*drg(3)-2.774d+03*drg(4)+1.901d+03*drg(5) )
         dqr= rf(j-1) + dm*(2.51d+02*drf(1)-1.274d+03*drf(2)+
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
         dp= rg(j-1) + dm*(-1.9d+01*drg(1)+1.06d+02*drg(2)
     >               -2.64d+02*drg(3)+6.46d+02*drg(4)+2.51d+02*drg(5) )
         dq= rf(j-1) + dm*(-1.9d+01*drf(1)+1.06d+02*drf(2)
     >               -2.64d+02*drf(3)+6.46d+02*drf(4)+2.51d+02*drf(5) )
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
c
c        =============================================================
c        check number of nodes
c        =============================================================
         if( rg(j-1).ne.zero .and. rg(j)/rg(j-1).le.zero ) nd=nd+1
         if( nd .gt. nodes ) go to 110
      enddo   
c
c     ================================================================
c     if no. of nodes is too small increase the energy and start again
c     ================================================================
      if( nd .eq. nodes ) return
c
      en=0.8d+00*en
c
      if( en .lt. -tol ) go to 10
c
c     ================================================================
c     this energy guess has become too high
c     ================================================================
      write(6,'('' OUTWS:: not enough nodes'')')
      call fstop(sname)
c
c     ================================================================
c     if no. of nodes is too big decrease the energy and start again
c     ================================================================
 110  en=1.2d+00*en
c
      if( en .gt. elim ) go to 10
c
c     ================================================================
c     this energy guess has become too low
c     ================================================================
      write(6,'('' OUTWS:: too many nodes'')')
      call fstop(sname)
c
c     ================================================================
      return
      end
