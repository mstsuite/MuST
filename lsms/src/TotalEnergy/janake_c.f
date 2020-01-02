c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine janake(vrold,vrnew,
     >                  rhotot,rho,corden,
     >                  rr,rins,rmt,
     >                  jmt,jws,
     >                  komp,atcon,
     >                  ztotss,
     >                  atvol,
     >                  vx,enxc,
     >                  evalsum,ecorv,esemv,
     >                  etot,press,rspin,
     >                  iprpts,
     >                  iprint,istop)
c     >                  iprint,istop,
c     >                  iprpts)

c     ================================================================
c
      implicit   none
      integer mynod
c
      character  istop*32
      character  sname*20
      parameter (sname='janake')
c
c YingWai's trial
!      include    'atom_param.h'
c
      integer    jmt
      integer    jws
      integer    komp
      integer    iprint
      integer    iprpts
c
      integer    ir
      integer    i
      integer    ik
      integer    iformu
      integer    iex
c
      real*8     ylag
      real*8     rmt
      real*8     vrold(iprpts,komp)
      real*8     vrnew(iprpts,komp)
      real*8     rho(iprpts,komp)
      real*8     rhotot(iprpts)
      real*8     corden(iprpts,komp)
      real*8     rr(iprpts),rins,rtmp(0:iprpts)
      real*8     atcon(komp)
      real*8     ztotss(komp)
      real*8     qtmp
      real*8     atvol
      real*8     vx(iprpts,komp)
      real*8     enxc(iprpts,komp)
      real*8     evalsum(komp)
      real*8     ecorv(komp)
      real*8     esemv(komp)
      real*8     etot
      real*8     press
      real*8     rspin
      real*8     ekinetic
      real*8     erho
      real*8     ezrho
      real*8     ecoulomb
      real*8     rhoint
c
      real*8     vmvold(iprpts)
      real*8     rinv(iprpts)
      real*8     rhov(iprpts)
      real*8     derv(iprpts)
      real*8     bndint(0:iprpts),dummy
      real*8     bnd(0:iprpts)
      real*8     exchen
      real*8     exchjmt
      real*8     coren
      real*8     valen
      real*8     correc
      real*8     cor3pv
      real*8     pterm1
      real*8     ezpt
      real*8     tpzpt
      real*8     evssum
      real*8     xcmt
c
      real*8     zero
      real*8     half
      real*8     two
      real*8     three
      real*8     four
c
      parameter (iformu=2)
      parameter (zero=0.0d0)
      parameter (half=0.5d0)
      parameter (two=2.0d0)
      parameter (three=3.0d0)
      parameter (four=4.0d0)
      mynod=-iprint
c
      rtmp(0)=zero
      do ir=1,jmt
         rinv(ir)=1/rr(ir)
	 rtmp(ir)=sqrt(rr(ir))
      enddo
      do ik=1,komp
c
c        ============================================================
c        calculate the zeropoint energy...........................
c        ============================================================
c        ------------------------------------------------------------
         call zeropt(ezpt,tpzpt,atvol,ztotss(ik))
c        ------------------------------------------------------------
c        ============================================================
c        calculate the energy eigenvalue sum for both valence and
c        sem-core electrons, and store it in evssum..................
c        ============================================================
	 evssum = evalsum(ik)+esemv(ik)
c
	if(iformu.eq.2) then
c Look at some terms of interest
         do i=1,jmt
            bndint(i)=rho(i,ik)*vrold(i,ik)*rinv(i)*rinv(i)
c           bndint(i)=rho(i,ik)*rho(i,ik)*rinv(i)*rinv(i)
         enddo
c        bndint(507)=bndint(506)
c        if(iprint.ge.0)write(6,'(''v'',10f12.4)')(vrold(i,1),i=1,jmt+2)
c        if(iprint.ge.0)write(6,'(''p'',10f12.4)')(rho(i,1),i=1,jmt+2)
c        if(iprint.ge.0)write(6,'(''r'',10f12.4)')(rinv(i),i=1,jmt+2)
c        if(iprint.ge.0)write(6,'(''b'',10f12.5)')(bndint(i),i=1,jmt+2)
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
c        ------------------------------------------------------------
         call newint(jmt+1,rtmp,bndint,bnd,3)
c        if(iprint.ge.0)write(6,'(''m'',10f12.5)')
c    >   (bnd(i+1)-bnd(i),i=1,jmt+2)
c        if(iprint.ge.0)write(6,'(''i'',10f12.5)')(bnd(i),i=1,jmt+2)
c        ------------------------------------------------------------
	 call interp(rtmp,bnd,jmt+1,sqrt(rmt),ekinetic,dummy,.false.)
c        write(6,'(''ekinetic1'',d19.10)')ekinetic
         ekinetic=ylag(sqrt(rmt),rtmp,bnd,0,3,jmt+1,i)
c        write(6,'(''ekinetic2'',d19.10)')ekinetic
c        ekinetic=two*bnd(jmt)
         ekinetic=two*ekinetic
	 ekinetic=ecorv(ik)+evssum-ekinetic

         do i=1,jmt
            bndint(i)=rhotot(i)*rinv(i)*rinv(i)
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,5)
! meis: test
!         call interp(rtmp,bnd,jmt+1,sqrt(rmt),qtmp,dummy,.false.)
!         write(6,'(''janake: coulomb v:'',f22.11)') qtmp
c        ------------------------------------------------------------
         do i=1,jmt
            bndint(i)=rho(i,ik)*bnd(i)*rinv(i)**3*8.d0
         enddo
	 call interp(rtmp,bndint,jmt+1,sqrt(rmt),qtmp,dummy,.false.)
c        qtmp=bnd(jmt)-ztotss(1)
         qtmp=qtmp-ztotss(1)
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,5)
c        ------------------------------------------------------------
	 call interp(rtmp,bnd,jmt+1,sqrt(rmt),erho,dummy,.false.)
c        erho=bnd(jmt)

         do i=1,jmt
            bndint(i)=-rho(i,ik)*ztotss(ik)*four*rinv(i)*rinv(i)
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,3)
c        ------------------------------------------------------------
	 call interp(rtmp,bnd,jmt+1,sqrt(rmt),ezrho,dummy,.false.)
c        ezrho=bnd(jmt)  !-rho(1,ik)*ztotss(ik)
	 ecoulomb=erho+ezrho

         do i=1,jmt
            bndint(i)=rho(i,ik)*enxc(i,ik)*rinv(i)*rinv(i)
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,5)
	 call interp(rtmp,bnd,jmt+1,sqrt(rmt),exchen,dummy,.false.)
c        ------------------------------------------------------------
         exchjmt=two*bnd(jmt)
         exchen=two*exchen

         do i=1,jmt
            bndint(i)=rho(i,ik)*(enxc(i,ik)-vx(i,ik))*rinv(i)*rinv(i)
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,5)
	 call interp(rtmp,bnd,jmt+1,sqrt(rmt),pterm1,dummy,.false.)
c        ------------------------------------------------------------
c        pterm1=two*bnd(jmt)
         pterm1=two*pterm1

c        *************************************************************
         etot = etot+atcon(ik)*(ekinetic+ecoulomb+exchen+ezpt/rspin)
         press= press+atcon(ik)*(two*ekinetic+ecoulomb-three*pterm1+
     &          tpzpt/rspin)
	 if(iprint.ge.0) then
         write(6,'(10x,''evssum'',t30,''='',f22.11)')evssum
         write(6,'(10x,''Kinetic E'',t30,''='',f22.11)') ekinetic
         write(6,'(10x,''Coulomb E(rho)'',t30,''='',f22.11)') erho
         write(6,'(10x,''Coulomb E(Z_rho)'',t30,''='',f22.11)') ezrho
         write(6,'(10x,''Coulomb E'',t30,''='',f22.11)') ecoulomb
         write(6,'(10x,''Exch E'',t30,''='',2f22.11)') exchen,exchjmt
	 endif
ctest    ****************************** emia *************************
c     write(6,'(''emia'',i4,10f15.7)')mynod,etot,press,evssum,
c    >ekinetic,erho,ezrho,ecoulomb,exchen,qtmp
ctest    ******************************* emia ************************
	else  ! iformu.eq.2
c        ============================================================
c        note: rhov contains both valence and semi-core charge
c              density...............................................
c        ============================================================
         do ir=1,jmt
            rhov(ir)=rho(ir,ik)-corden(ir,ik)
            vmvold(ir)=vrold(ir,ik)-vrnew(ir,ik)
         enddo
c
c        ============================================================
c        calculate coren for each component where one piece of 
c        dv/dr is equal to ecorv. add this piece to derv term and 
c        calculate: coren =-int r3 dr 2pi*corden*dv/dr         iformu
c                         = int r2 dr 2pi*corden*(v-d(rv)/dr)       1
c                         = sum ec - int r2 dr 4pi*corden*d(rv)/dr  0
c        ============================================================
c        ------------------------------------------------------------
         call newder(vrold(1,ik),derv,rtmp(1),jmt)
c        ------------------------------------------------------------
         if(iformu.eq.0) then
            do i=1,jmt
               bndint(i)=-corden(i,ik)*derv(i)/(rr(i)*rtmp(i))
            enddo
c           ---------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
            call newint(jmt+1,rtmp,bndint,bnd,3)
c           ---------------------------------------------------------
            coren=ecorv(ik)+bnd(jmt)
	 else
	    do i=1,jmt
               bndint(i)=corden(i,ik)*(rinv(i)*vrold(i,ik)-
     &             derv(i)/(two*rtmp(i)))
            enddo
c           ---------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
            call newint(jmt+1,rtmp,bndint,bnd,1)
c           ---------------------------------------------------------
            coren = bnd(jmt)
         endif
c
c        ============================================================
c        calculate: pterm1 =-int r1 dr 4pi*valden d(r2v)/dr
c                          =-int r2 dr 4pi*valden (d(rv)/dr+v)
c        ============================================================
         do i=1,jmt
           bndint(i)=rhov(i)*(derv(i)/(two*rtmp(i))+vrold(i,ik)*rinv(i))
     &         *rinv(i)
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,3)
c        ------------------------------------------------------------
         pterm1=-two*bnd(jmt)
c
c        ============================================================
c        calculate: valen =-int r2 dr 4pi rhov*d(rv)/dr..............
c        ============================================================
c        ------------------------------------------------------------
         call newder(vrnew(1,ik),derv,rtmp(1),jmt)
c        ------------------------------------------------------------
         do i=1,jmt
            bndint(i)=-rhov(i)*derv(i)/(rr(i)*rr(i)*rtmp(i))
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,5)
c        ------------------------------------------------------------
         valen = bnd(jmt)
c
c        ============================================================
c        calculate: exchen = int r2 dr 4pi rho (4exc-3vxc)
c        ============================================================
         do i=1,jmt
            bndint(i)=rho(i,ik)*(four*enxc(i,ik)-three*vx(i,ik))
     &        *rinv(i)*rinv(i)
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,5)
c        ------------------------------------------------------------
         exchen=two*bnd(jmt)
c
c        ============================================================
c        calculate terms that make the energy variational: 
c                  correc = int 4pi dr r3*corden*d[vold-vnew]/dr
c                          -int 4pi dr r2*rhov*[vold-vnew]
c                         = int 4pi dr r2*corden*d[r*(vold-vnew)]/dr
c                          -int 4pi dr r2*rho*[vold-vnew]
c        calculate terms that are corrections to the pressure:
c                  cor3pv = int 8pi dr r3*corden*d[vold-vnew]/dr
c                         = int 8pi dr r2*corden*d[r*(vold-vnew)]/dr
c                          -int 8pi dr r2*corden*[vold-vnew]
c        ============================================================
c        ------------------------------------------------------------
         call newder(vmvold,derv,rtmp(1),jmt)
c        ------------------------------------------------------------
         do i=1,jmt
            bndint(i)=corden(i,ik)*derv(i)/(rr(i)*rtmp(i))
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,3)
c        ------------------------------------------------------------
         correc=bnd(jmt)
         cor3pv=two*correc
         do i=1,jmt
            bndint(i)=rho(i,ik)*vmvold(i)/(rr(i)*rr(i))
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,3)
c        ------------------------------------------------------------
         correc=correc-two*bnd(jmt)
         do i=1,jmt
            bndint(i)=corden(i,ik)*vmvold(i)/(rr(i)*rr(i))
         enddo
c        ------------------------------------------------------------
	 call interp(rtmp(1),bndint(1),4,zero,bndint(0),dummy,.false.)
         call newint(jmt+1,rtmp,bndint,bnd,3)
c        ------------------------------------------------------------
         cor3pv=cor3pv-four*bnd(jmt)
c
         xcmt=rr(jmt)*rho(jmt,ik)*( enxc(jmt,ik)-vx(jmt,ik) ) 
c
c        ============================================================
c        sum over terms to get the total energy and 3PV...........
c        ============================================================
         etot = etot+atcon(ik)*( coren+valen+exchen+correc+ezpt/rspin
     >                          +evssum-xcmt )
         press= press+atcon(ik)*( pterm1+cor3pv+tpzpt/rspin
     >                           +two*evssum-xcmt)
c
c        *************************************************************
c        Major print out for current sublattice.......................
c        *************************************************************
         if(iprint.ge.0) then
            write(6,'(/,'' JANAKE:: Component'',t30,''='',i5)') ik
c           write(6,'(10x,''ecorv'',t30,''='',f22.11)')ecorv(ik)
c           write(6,'(10x,''esemv'',t30,''='',f22.11)')esemv(ik)
c           write(6,'(10x,''evalsum'',t30,''='',f22.11)')evalsum(ik)
            write(6,'(10x,''coren'',t30,''='',f22.11)')coren
	    write(6,'(10x,''valen'',t30,''='',f22.11)')valen
            write(6,'(10x,''exchen'',t30,''='',f22.11)')exchen
            write(6,'(10x,''correc'',t30,''='',f22.11)') correc
            write(6,'(10x,''evssum'',t30,''='',f22.11)')evssum
   	    write(6,'(10x,''-xcmt'',t30,''='',f22.11)')-xcmt
	    write(6,'(10x,''c3pv'',t30,''='',f22.11)')cor3pv
         endif
ctest    ****************************** emia *************************
c     write(6,'(''emia'',i4,10f15.7)')
c    >mynod,etot,press,evalsum(1),coren,valen,exchen,correc,qtmp
ctest    ******************************* emia ************************
c        *************************************************************
	 endif  ! iformu.eq.2
	 if(iprint.ge.0) then
            write(6,'(10x,''pterm1'',t30,''='',f22.11)')pterm1
            write(6,'(10x,''ezpt/spin'',t30,''='',f22.11)')ezpt/rspin
            write(6,'(10x,''tpzpt'',t30,''='',f22.11)')tpzpt/rspin
         endif
      enddo
c
c     ===============================================================
      if( istop .eq. sname ) then
         call fstop(sname)
      else
         return
      endif
      end
