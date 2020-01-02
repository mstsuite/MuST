c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getcor(n_spin_pola,mtasa,
     >                  jmt,jws,r_mesh,h,xstart,vr,
     >                  numc,nc,lc,kc,ec,
     >                  ztotss,zsemss,zcorss,
     >                  ecorv,esemv,corden,semcor,
!ebot,etopcor,
     >                  nrelc,
     >                  qcpsc_mt,qcpsc_ws,mcpsc_mt,mcpsc_ws,
     >                  iprpts, ipcore,
     >                  iprint,istop)
c     ================================================================
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine getcor(IN n_spin_pola,IN mtasa,
!    >                  IN jmt,IN jws,IN r_mesh,IN h,IN xstart,IN vr,
!    >                  IN numc,IN nc,IN lc,IN kc,INOUT ec,
!    >                  IN ztotss,IN zsemss,IN zcorss,
!    >                  OUT ecorv,OUT esemv,OUT corden,OUT semcor,
!    >                  IN nrelc,
!    >                  OUT qcpsc_mt,OUT qcpsc_ws,OUT mcpsc_mt,OUT mcpsc_ws,
!    >                  IN iprpts, IN ipcore,
!    >                  IN iprint,IN istop)
c     ================================================================
!     use MPPModule, only : GlobalMax
c
      implicit   none
c
      character  istop*32
      character  sname*32
c
!      include    'atom_param.h'
c
      integer    iprpts,ipcore
      integer    n_spin_pola
      integer    mtasa
      integer    jmt
      integer    jws
      integer    last
      integer    last2
      integer    numc
      integer    nc(ipcore,n_spin_pola)
      integer    lc(ipcore,n_spin_pola)
      integer    kc(ipcore,n_spin_pola)
      integer    nrelc
      integer    ndeepz
      integer    iprint
      integer    is
      integer    j
c
      real*8     r_mesh(iprpts),rtmp(0:iprpts)
      real*8     h
      real*8     xstart
      real*8     vr(iprpts,n_spin_pola)
      real*8     ec(ipcore,n_spin_pola)
      real*8     ztotss
      real*8     zsemss
      real*8     zcorss
      real*8     qsemmt
      real*8     qcormt
      real*8     qsemws
      real*8     qcorws
      real*8     qcorout
      real*8     qcpsc_mt
      real*8     qcpsc_ws
      real*8     mcpsc_mt
      real*8     mcpsc_ws
      real*8     ecorv(n_spin_pola)
      real*8     esemv(n_spin_pola)
      real*8     corden(iprpts,n_spin_pola)
      real*8     semcor(iprpts,n_spin_pola)
!      real*8     ebot
      real*8     wrk1(0:iprpts),wrk2(0:iprpts)
      real*8     zero
      real*8     tolch
!     real*8     etopcor
      real*8     work
c
      parameter (zero=0.0d0)
      real*8 two
      parameter (two=2.0d0)
      parameter (tolch=.000001d0)
      parameter (sname='getcor')
c
      if(n_spin_pola.eq.2) then
         do is=1,numc
            if( nc(is,1).ne.nc(is,2) .or.
     >          lc(is,1).ne.lc(is,2) .or. kc(is,1).ne.kc(is,2) )then
               write(6,'('' i,nc(i,1),nc(i,2) ='',3i5)')
     >                      is,nc(is,1),nc(is,2)
               write(6,'('' i,lc(i,1),lc(i,2) ='',3i5)')
     >                      is,lc(is,1),lc(is,2)
               write(6,'('' i,kc(i,1),kc(i,2) ='',3i5)')
     >                      is,kc(is,1),kc(is,2)
               call fstop(sname)
            endif
         enddo
      endif
c
c     ================================================================
c     calculate the core and semi-core densities......................
c     ================================================================
      ndeepz=(zcorss/n_spin_pola+.5d0)
      last=iprpts
      if(mtasa.eq.0) then
        last2=last
      else
        last2=jws
      endif
!      etopcor=-10.0d+20
c     ----------------------------------------------------------------
      call zeroout(corden,iprpts*n_spin_pola)
      call zeroout(semcor,iprpts*n_spin_pola)
      call zeroout(ecorv,n_spin_pola)
      call zeroout(esemv,n_spin_pola)
c     ----------------------------------------------------------------
      qsemmt=zero
      qcormt=zero
      qsemws=zero
      qcorws=zero
      mcpsc_mt=zero
      mcpsc_ws=zero
      do is=1,n_spin_pola
         if(iprint.ge.0) then
            write(6,'(/,'' GETCOR:: Spin index'',t40,''='',1i5)')is
         endif
          if(iprint.ge.2) then
            write(*,*) jmt,jws,last,last2,h,r_mesh(1),vr(1,is),ztotss
            write(*,*) numc,nc(1,is),lc(1,is),kc(1,is),ec(1,is)
            write(*,*) corden(1,is),semcor(1,is)
            write(*,*) ecorv(is),esemv(is)
            write(*,*) n_spin_pola,nrelc,ndeepz,iprpts,iprint,istop
          endif
c        -------------------------------------------------------------
         call corslv(jmt,jws,last,last2,h,r_mesh,vr(1,is),ztotss,
     >               numc,nc(1,is),lc(1,is),kc(1,is),ec(1,is),
     >               corden(1,is),semcor(1,is),
     >               ecorv(is),esemv(is),
     >               n_spin_pola,nrelc,ndeepz,iprpts,iprint,istop)
c        -------------------------------------------------------------
c
c        =============================================================
c        adjust bottom of cotour: above highest semi-core state
c        =============================================================
         if(numc.gt.0) then
!            etopcor=max(ec(numc,is),etopcor)
c           ==========================================================
c           check the charges.........................................
c           ----------------------------------------------------------
c           call simpun(r_mesh,semcor(1,is),iprpts,1,wrk1)
	    rtmp(0)=0.d0
	    wrk2(0)=0.d0
	    do j=1,last2
	      wrk2(j)=semcor(j,is)/r_mesh(j)
	      rtmp(j)=sqrt(r_mesh(j))
	    enddo
            call newint(last2+1,rtmp,wrk2,wrk1,3)
c           ----------------------------------------------------------
            qsemmt=qsemmt+two*wrk1(jmt)
            qsemws=qsemws+two*wrk1(last2)
            mcpsc_mt=mcpsc_mt+two*(3-2*is)*wrk1(jmt)
            mcpsc_ws=mcpsc_ws+two*(3-2*is)*wrk1(last2)
c           ----------------------------------------------------------
c           call simpun(r_mesh,corden(1,is),iprpts,1,wrk1)
	    wrk2(0)=0.d0
	    do j=1,last2
	      wrk2(j)=corden(j,is)/r_mesh(j)
	    enddo
            call newint(last2+1,rtmp,wrk2,wrk1,3)
c           ----------------------------------------------------------
            qcormt=qcormt+two*wrk1(jmt)
            qcorws=qcorws+two*wrk1(last2)
            mcpsc_mt=mcpsc_mt+two*(3-2*is)*wrk1(jmt)
            mcpsc_ws=mcpsc_ws+two*(3-2*is)*wrk1(last2)
         endif
      enddo
      qcorout=zsemss+zcorss-qsemmt-qcormt
      qcpsc_mt=qcormt+qsemmt
      qcpsc_ws=qcorws+qsemws
      if(iprint.ge.0 .and. abs(qcorout).ge.0.001) then
         write(6,'(/,'' GETCOR:: Muffin-tin core charge'',t40,''='',
     >   f16.11)') qcormt
         write(6,'(10x,''Muffin-tin semicore charge'',t40,''='',
     >   f16.11)') qsemmt
      endif
      if(iprint.ge.0) then
         write(6,'(/,'' GETCOR:: Muffin-tin core charge'',t40,''='',
     >   f16.11)') qcormt
         write(6,'(10x,''Muffin-tin semicore charge'',t40,''='',
     >   f16.11)') qsemmt
         write(6,'(10x,''Muffin-tin core+semi moment'',t40,''='',
     >   f16.11)') mcpsc_mt
         write(6,'(10x,''Wigner-Seitz core+semi moment'',t40,''='',
     >   f16.11)') mcpsc_ws
         write(6,'(10x,''Interstitial charge core'',t40,''='',
     >   f16.11)') qcorout
      endif
c     ***************************************************************
      if(iprint.ge.1) then
         write(6,'(/,'' GETCOR::'')')
         do j=1,last2,50
            write(6,'('' j,corden[1:2]: '',i5,3d20.12)')
     >      j,corden(j,1),corden(j,2),corden(j,1)-corden(j,2)
         enddo
         do j=1,last2,50
            write(6,'('' j,semcor[1:2]: '',i5,3d20.12)')
     >      j,semcor(j,1),semcor(j,2),semcor(j,1)-semcor(j,2)
         enddo
      endif
c     ***************************************************************
c     ===============================================================
c     Find the uppermost core state across all the nodes/atoms........
c     ---------------------------------------------------------------
!      call GlobalMax(etopcor)
c     ---------------------------------------------------------------
c
!     if(iprint.eq.0) then
!        write(6,'(10x,''Upper limit of core level'',t40,''='',
!    >             f16.11)')etopcor
!     endif
c
c     ===============================================================
c     check that uppermost core state is well below ebot.............
c     ===============================================================
!     if(etopcor+0.1d0 .gt. ebot) then
!        write(6,'('' GETCOR: etopcor+0.1 .gt. ebot'',2d12.5)')
!    >         etopcor,ebot
!        call kill(0,9)
!        call fstop(sname)
!     endif
c
c     ===============================================================
c     check that semi-core and deep-core charge is not being lost....
c     ===============================================================
      if(abs(qsemws-zsemss).gt.tolch .and. iprint.ge.0) then
         write(6,'(/,10x,''Lost semi-core charge'')')
         write(6,'(10x,''z='',d17.8,'' q='',f17.8)') zsemss,qsemws
c        call fstop(sname)
      else if(abs(qcorws-zcorss).gt.tolch .and. iprint.ge.0) then
         write(6,'(/,10x,''Lost      core charge'')')
         write(6,'(10x,''z='',d17.8,'' q='',f17.8)') zcorss,qcorws
c        call fstop(sname)
      endif
c
c     ================================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
c
      end
