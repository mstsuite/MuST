

      subroutine F_readSingleAtomData_bigcell(fname_c,fname_l, &
     &     header,jmt,jws,xstart, &
     &     rmt,alat,efermi,vdif,ztotss,zcorss, &
     &     nspin,numc,xvalws, &
     &     vr,rhotot,corden,v_dim, &
     &     ec,nc,lc,kc,c_dim)
      implicit none

      byte fname_c(*)
      integer fname_l,i,j,is,jwsc
      integer jmt,jws,nspin,numc,v_dim,c_dim
      integer nc(c_dim,*),lc(c_dim,*),kc(c_dim,*)
      real*8 xvalws(2)
      real*8 xstart,rmt,alat,efermi,vdif,ztotss,zcorss
      real*8 vr(v_dim,*),rhotot(v_dim,*)
      real*8 ec(c_dim,*)
      real*8 xmt,v0
      real*8 corden(v_dim,*)
      character*80 header,jtitle
      character*127 fname
      character lst(30,2)*5
      integer funit

!     write(*,*) fname_l
!     write(*,'(127a)') (fname_c(i),i=1,fname_l)
      write(fname,'(127a)') (fname_c(i),i=1,fname_l)
!
!     write(*,*) 'reading from:',fname

      funit=55
      open(unit=funit,file=fname, &
     &     status='old',form='formatted')
!            ----------------------------------------------------------
!            write(6,'('' SWVDAT:: mynod,ib,vfname(ib):'',2i5,a60)')
!     >      mynod,ib,vfname(ib)
!
!           ==========================================================
!           read in the formatted potential data......................
!           ==========================================================

!      write(6,*) 'reading file //', trim(fname), '//'
!      flush(unit=6)

      read(funit,'(a80)') jtitle
      header=jtitle
!      header=trim(jtitle)//char(0)
      read(funit,*) nspin,vdif
!     write(6,'('' SWVDAT:: ib,n_spin1,header:'',2i5,a80,i5)')
!     >      ib,n_spin1,header
      do is=1,nspin
         read(funit,'(1a80)')jtitle
         read(funit,'(f5.0,17x,f12.5,f5.0,e20.13)') &
     &        ztotss,alat,zcorss,efermi
!        write(6,'('' SWVDAT:: mynod,ib,is,ztss,zcss,jtitle:'',
!    >            3i5,2f5.0,a40)') mynod,ib,is,ztss,zcss,jtitle
         read(funit,'(17x,2e20.13,i5)') xstart,xmt,jmt
         rmt=exp(xmt)
         if(jmt.gt.v_dim) then
            write(*,*) "!! jmt exeeds v_dim !!"
            stop
         endif
         read(funit,'(4e20.13)') (vr(j,is),j=1,jmt)
         read(funit,'(35x,e20.13)') v0
!        ====================================================
!        read in formatted total charge density..............
!        ====================================================
         read(funit,'(i5,e20.13)') jws,xvalws(is)
         if(jws.gt.v_dim) then
            write(*,*) "!! jws exeeds v_dim !!"
            stop
         endif
         read(funit,'(4e20.13)') (rhotot(j,is),j=1,jws)
!        ====================================================
!        read in formatted core state information............
!        ====================================================
         read(funit,'(2i5)') numc,jwsc
         if(numc.gt.c_dim) then
            write(*,*) "!! numc exeeds c_dim !!"
            write(*,*) numc, c_dim
            stop
         endif
         if(numc.gt.0) then
            read(funit,'(3i5,f12.5,1x,a5)') &
     &           (nc(j,is),lc(j,is), &
     &           kc(j,is),ec(j,is),lst(j,is),j=1,numc)
         endif
         if (jwsc.gt.0) then
            read(funit,'(4e20.13)') (corden(j,is),j=1,jwsc)
         endif
      enddo
      close(unit=funit)

      end
