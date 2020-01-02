c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newexchg(n_spin_pola,sp,
     >                  rhoup,rhodn,
     >                  vx,enxc,vxout,excout,ro3,dz,
     >                  r_mesh,jmt,iexch)
c     ================================================================
c
      implicit   none
c
      integer    n_spin_pola
      integer    jmt
c
      integer    iexch
      integer    ir
c
      real*8     sp
      real*8     rhoup(jmt)
      real*8     rhodn(jmt)
      real*8     vx(jmt)
      real*8     enxc(jmt)
      real*8     vxout,excout,ro3,dz
      real*8     r_mesh(jmt)
c
      real*8     rhot
      real*8     alpha2
      real*8     dzr
      real*8     ro3r
      real*8     zero
      real*8     one
      real*8     three
      real*8     third
c
      parameter  (zero=0.0d0)
      parameter  (one=1.0d0)
      parameter  (three=3.0d0)
      parameter  (third=one/three)
c
c     ================================================================
c     calculate exchange-correlation potential and density...........
c     ================================================================
c
c        -------------------------------------------------------------
      vxout=0.d0
      if(dz.lt.-1.d0)write(*,*)'newexch, dz,ro3',dz,ro3
      if(ro3.lt.1.d9)vxout=alpha2(ro3,dz,sp,iexch,excout)
c        -------------------------------------------------------------
      do ir=1,jmt
         rhot=rhoup(ir)+(n_spin_pola-1)*rhodn(ir)
c        if(rhot.gt.zero) then
	 if(rhot.gt.1.d-9) then
           dzr=(rhoup(ir)-rhodn(ir))/rhot
           if(dzr.gt.+1.d0)dzr=+1.d0
           if(dzr.lt.-1.d0)dzr=-1.d0
           ro3r=(three*r_mesh(ir)**2/rhot)**third
           vx(ir)=alpha2(ro3r,dzr,sp,iexch,enxc(ir))
           if(.not.(enxc(ir).lt.0))then
             write(6,*)'newexch',ir,rhot,dzr,ro3r,enxc(ir)
             enxc(ir)=0.d0
           endif
	 else
	   vx(ir)=zero
	   enxc(ir)=zero
	 endif
      enddo
c
      return
      end
