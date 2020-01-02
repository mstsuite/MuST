c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getqm_mt(n_spin_pola,ilast,rsmo,sqr,rho,iprpts,rhot,
     >                    i_smooth_rsmo,qtotmt,mtotmt,r_sph,iprint)
c     ================================================================
c     calculates the charge inside r_sph
c     cuts off integral smoothly at rsmo if i_smooth_rsmo.ge.2
c
      implicit   none
c
      integer    n_spin_pola
      integer    j
      integer    ilast
      integer    is
      integer    i_smooth_rsmo,iprint,iprpts
c
      real*8     rsmo,sqr(0:ilast),sqrsmo
      real*8     rho(iprpts,n_spin_pola)
      real*8     rhot(0:ilast,2),drhot
      real*8     qtotmt,qmtis,r_sph
      real*8     mtotmt
      real*8     width
      real*8     one
      parameter (one=1.d0)
c
c     ================================================================
      qtotmt=0.0d0
      qmtis=0.0d0
      mtotmt=0.0d0
      width=0.0d0
      do is=1,n_spin_pola
         if(i_smooth_rsmo.le.1) then
           do j=1,ilast
             rhot(j,1)=2.d0*rho(j,is)/(sqr(j)**4)
           enddo
         else
           width=sqr(ilast)-sqr(ilast-1)
	   sqrsmo=sqrt(rsmo)
           do j=1,ilast
             rhot(j,1)=2.d0*rho(j,is)/((one+exp((sqr(j)-sqrsmo)/width))
     &         *sqr(j)**4)
           enddo
         endif
c        call interp(sqr(1),rho(1,1),4,0.d0,rhot(0,1),drhot,.false.)
         call interp(sqr(1),rhot(1,1),4,0.d0,rhot(0,1),drhot,.false.)
c        -------------------------------------------------------------
         call newint(ilast+1,sqr,rhot(0,1),rhot(0,2),5)
c        -------------------------------------------------------------
         call interp
     >   (sqr(1),rhot(1,2),ilast,sqrt(r_sph),qmtis,drhot,.false.)
         qtotmt=qtotmt+qmtis
c        mtotmt=mtotmt+(n_spin_pola-1)*(3.0d0-2.0d0*is)*rhot(ilast,2)
         mtotmt=mtotmt+(n_spin_pola-1)*(3.0d0-2.0d0*is)*qmtis
      enddo
c
      if(iprint.ge.0) then
         write(6,'(/,'' GETQM_MT:: r_sph charge,w'',t40,''='',
     >   3f18.11,2i5)')qtotmt,width,sqr(ilast)**2,ilast,i_smooth_rsmo
         write(6,'(''            muffin-tin moment'',t40,''='',
     >                1f18.11)')mtotmt
      endif
c
c     ================================================================
      return
      end
