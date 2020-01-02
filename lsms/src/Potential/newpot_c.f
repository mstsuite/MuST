c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newpot(n_spin_pola,ztotss,
     >                  rhoup,rhodn,rhot,
     >                  vrold,vrnew,vrms,vx,
     >                  vmt1,vmt,vxout,
     >                  sqr,jmt,rins,rws,
     >                  mtasa,iexch)
c     ================================================================
c
      implicit   none
c
      integer    n_spin_pola
      integer    jmt
      integer    mtasa
c
      integer    ir
      integer    iexch
c
      real*8     alpgga
      real*8     epcorr
      real*8     epcorrave
      real*8     ztotss
      real*8     rhoup(jmt)
      real*8     rhodn(jmt)
      real*8     vrold(jmt)
      real*8     vrnew(jmt)
      real*8     vrms
      real*8     vx(jmt)
      real*8     vmt1
      real*8     vmt,vxout
      real*8     sqr(0:jmt),rins,sqrmt,rws
c
      real*8     rhot(0:jmt,2),drhot,rhot_rins
      real*8     vhart
      real*8     one
      real*8     two
      real*8     width
c
      parameter  (one=1.0d0)
      parameter  (two=2.0d0)
c
c     ================================================================
c     calculate hartree and exchange-correlation potential ...........
c     then subtract muffin tin zero...................................
c     ================================================================
c
      do ir=1,jmt
c       2.d0 is the jacobian needed because we are integrating over sqrt(r)
        rhot(ir,1)=2.d0*(rhoup(ir)+(n_spin_pola-1)*rhodn(ir))/sqr(ir)**4
c YingWai's check
!       write(6,'(''ir,rhot,rhoup,rhodn = '',i8,3(1x,f20.12))')
!    &        ir, rhot(ir,1), rhoup(ir), rhodn(ir)
      enddo
      call interp(sqr(1),rhot(1,1),4,0.d0,rhot(0,1),drhot,.false.)
c     ----------------------------------------------------------------
      call newint(jmt+1,sqr,rhot(0,1),rhot(0,2),5)
      do ir=1,jmt
         vrnew(ir)=(-ztotss+rhot(ir,2))/(sqr(ir)*sqr(ir))
      enddo
      call newint(jmt+1,sqr,rhot(0,1),rhot(0,2),3)
      if(mtasa.eq.2)then
         call 
     >   interp(sqr,rhot(0,2),jmt+1,sqrt(rws),rhot_rins,drhot,.false.)
      else
         call 
     >   interp(sqr,rhot(0,2),jmt+1,sqrt(rins),rhot_rins,drhot,.false.)
      endif
c     ----------------------------------------------------------------
      do ir=1,jmt
         vhart=two*(vrnew(ir)+rhot(jmt,2)-rhot(ir,2))
!        write(6,'(''ir,vhart = '',i8,1x,f24.12)') ir, vhart
ctest    vhart=two*(vrnew(ir)+rhot_rins-rhot(ir,2))
         vrnew(ir)=(vhart+vx(ir)-vmt1-vmt-vxout)*sqr(ir)*sqr(ir)
!         if(iexch.ge.100)
!     >   vrnew(ir)=-vrnew(ir)+sqr(ir)**2*(-vxout+
!     >   epcorr(.5d0*rhot(ir,1),0.d0,sqr(ir)**2,0,alpgga,40.d0))
      enddo
!     write(6,'(''vmt1,vmt,vxout = '',3(1x,f24.12))')
!    &         vmt1,vmt,vxout
      if(mtasa.ge.2) then
      width=sqr(jmt)-sqr(jmt-1)
      sqrmt=sqrt(rins)
      do ir=1,jmt
!        if(iexch.lt.100)vrnew(ir)=vrnew(ir)+vmt*sqr(ir)*sqr(ir)
         vrnew(ir)=vrnew(ir)+vmt*sqr(ir)*sqr(ir)
     *   /(one+exp((sqrmt-sqr(ir))/width))
!        if(iexch.ge.100)vrnew(ir)=vrnew(ir)-vmt*sqr(ir)*sqr(ir)
!     *   /(one+exp((sqrmt-sqr(ir))/width))
      enddo
      endif
c use rhot as temp storage
      rhot(0,1)=0.d0
      do ir=1,jmt
         rhot(ir,1)=2.d0*(vrold(ir)-vrnew(ir))**2
c YingWai's check
!        write(6,'(''ir,rhot,vrold,vrnew = '',i8,3(1x,f20.12))')
!    &    ir, rhot(ir,1),
!    &    vrold(ir), vrnew(ir)
      enddo
c     ----------------------------------------------------------------
      call newint(jmt+1,sqr,rhot(0,1),rhot(0,2),1)
c     ----------------------------------------------------------------
      vrms=sqrt(3.d0*rhot(jmt,2)/sqr(jmt)**6)
c
      return
      end
