      subroutine F_writeSingleAtomData_bigcell(fname_c,fname_l, &
     &                   header_c,header_l,jmt,jws,xstart, &
     &                   rmt,alat,efermi,vdif,ztotss,zcorss, &
     &                   n_spin_pola, numc, xvalws, &
     &                   vr,rhotot,corden, v_dim, &
     &                   atname,zsemss,zvalss, &
     &                   ec, nc,lc,kc, c_dim)
!     ================================================================

      implicit   none

      byte fname_c(*)
      byte header_c(*)
      integer fname_l
      integer header_l

      integer    n_spin_pola
      integer c_dim,v_dim

      character  header*80
      character*127  fname
      character*2 atname

      integer    jmt
      integer    jws
      integer    numc
      integer    nc(c_dim,n_spin_pola)
      integer    lc(c_dim,n_spin_pola)
      integer    kc(c_dim,n_spin_pola)
      integer    ic
      integer    ir
      integer    is
      integer    js
      integer i,j

      real*8     rmt
      real*8     rws
      real*8     h
      real*8     xstart
      real*8     vr(v_dim,n_spin_pola)
      real*8     vdif
      real*8     rhotot(v_dim,n_spin_pola)
      real*8     corden(v_dim,n_spin_pola)
      real*8     xvalws(n_spin_pola)
      real*8     ztotss
      real*8     zcorss
      real*8     zsemss
      real*8     zvalss
      real*8     ec(c_dim,n_spin_pola)
      real*8     alat
      real*8     efermi
      real*8     v0
      real*8     xmt
      real*8     zero

      parameter (zero=0.0d0)

      write(fname,'(127a)') (fname_c(i),i=1,fname_l)
      write(header,'(80a)') (header_c(i),i=1,header_l)
!        =============================================================
         xmt=log(rmt)
         v0=zero
!        -------------------------------------------------------------
         open(unit=30,file=fname,status='replace')
!        -------------------------------------------------------------
         write(30,'(a80)') header
         write(30,'(i5,2x,d20.13)') n_spin_pola,vdif
         do is=1,n_spin_pola
            write(30,'('' LSMS:'',t18,a2,t25,''z='',f4.0,t35, &
     &            ''xvalws='',f10.5)') atname,ztotss,xvalws(is)
            write(30,'(f5.0,17x,f12.5,f5.0,d20.13)') &
     &                 ztotss,alat,zcorss,efermi
            write(30,'(17x,2d20.13,i5)') xstart,xmt,jmt
            write(30,'(4d20.13)') (vr(ir,is),ir=1,jmt)
            write(30,'(35x,d20.13)') v0
            write(30,'(i5,d20.13)') jws,xvalws(is)
            write(30,'(4d20.13)') (rhotot(ir,is),ir=1,jws)
            write(30,'(2i5)') numc,jws
!           write(30,'(3i5,f12.5,2x,a5)')
!    >           (nc(j,is),lc(j,is),kc(j,is),ec(j,is),lst(j,is),
!    >                                                    j=1,numc)
            write(30,'(3i5,f12.5,7x)') &
     &           (nc(j,is),lc(j,is),kc(j,is),ec(j,is), j=1,numc)
            write(30,'(4d20.13)') (corden(ir,is),ir=1,jws)
!        -------------------------------------------------------------
         enddo

         close(unit=30)
!        -------------------------------------------------------------
      return
      end
