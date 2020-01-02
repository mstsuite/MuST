c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getstruc(mynod,vbrar,omegbra,eta,
     >                    iprslat_mad,ipknlat_mad,
     >                    rslat_x,rslat_y,rslat_z,rslatsq,nrslat,
     >                    knlat_x,knlat_y,knlat_z,knlatsq,nknlat,
     >                    iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  sname*32
      character  istop*32
c
      integer    mynod
      integer    iprslat_mad
      integer    ipknlat_mad
      integer    nrslat
      integer    nknlat
      integer    n1
      integer    nm1
      integer    nm2
      integer    nm3
      integer    iprint
c
      real*8     eta
      real*8     rslat_x(iprslat_mad)
      real*8     rslat_y(iprslat_mad)
      real*8     rslat_z(iprslat_mad)
      real*8     rslatsq(iprslat_mad)
      real*8     vbrar(3,3)
      real*8     knlat_x(ipknlat_mad)
      real*8     knlat_y(ipknlat_mad)
      real*8     knlat_z(ipknlat_mad)
      real*8     knlatsq(ipknlat_mad)
      real*8     vbrak(3,3)
      real*8     rscut
      real*8     kncut
      real*8     omegbra
      real*8     two
      real*8     fnpi
      real*8     factor
c
      parameter (two=2.0d0)
      parameter (sname='getstruc')
c
c     ****************************************************************
c     Sets up real space and receprocal space Bravais lattice vectors
c     ****************************************************************
c
c     ================================================================
c     calculate rscut, the radius of real space truncation sphere.....
c     ----------------------------------------------------------------
      call getrscut(eta,vbrar(1,1),vbrar(1,2),vbrar(1,3),
     >              rscut,nm1,nm2,nm3)
ctest
c     nm1=6
c     nm2=6
c     nm3=6
c     rscut=44.d0*8.d0/6.2831853d0
ctest
c     ----------------------------------------------------------------
c
c     ================================================================
c     generate the real space lattice vectors.........................
c     ----------------------------------------------------------------
      call lattice(vbrar,rscut,nm1,nm2,nm3,
     >             rslat_x,rslat_y,rslat_z,rslatsq,nrslat,iprslat_mad,
     >             iprint,istop)
c     ----------------------------------------------------------------
      if(iprint.ge.1) then
         write(6,'(/,'' GETSTRUC:: nm1,nm2,nm3   = '',3i5)')nm1,nm2,nm3
         write(6,'(  ''            Rs cut radius = '',1f10.5)') rscut
         write(6,'(  ''            Number of Rs  = '',i5)') nrslat
      endif
c     ================================================================
c     calculate the bravais lattice cell volume......................
      omegbra=(vbrar(2,1)*vbrar(3,2)-vbrar(3,1)*vbrar(2,2))*vbrar(1,3)+
     >        (vbrar(3,1)*vbrar(1,2)-vbrar(1,1)*vbrar(3,2))*vbrar(2,3)+
     >        (vbrar(1,1)*vbrar(2,2)-vbrar(2,1)*vbrar(1,2))*vbrar(3,3)
      omegbra=abs(omegbra)
c
c     ================================================================
c     generate basis vectors for reciprocal space....................
      factor=two*fnpi()/omegbra
      vbrak(1,1)=factor*(vbrar(2,2)*vbrar(3,3)-vbrar(3,2)*vbrar(2,3))
      vbrak(2,1)=factor*(vbrar(3,2)*vbrar(1,3)-vbrar(1,2)*vbrar(3,3))
      vbrak(3,1)=factor*(vbrar(1,2)*vbrar(2,3)-vbrar(2,2)*vbrar(1,3))
      vbrak(1,2)=factor*(vbrar(2,3)*vbrar(3,1)-vbrar(3,3)*vbrar(2,1))
      vbrak(2,2)=factor*(vbrar(3,3)*vbrar(1,1)-vbrar(1,3)*vbrar(3,1))
      vbrak(3,2)=factor*(vbrar(1,3)*vbrar(2,1)-vbrar(2,3)*vbrar(1,1))
      vbrak(1,3)=factor*(vbrar(2,1)*vbrar(3,2)-vbrar(3,1)*vbrar(2,2))
      vbrak(2,3)=factor*(vbrar(3,1)*vbrar(1,2)-vbrar(1,1)*vbrar(3,2))
      vbrak(3,3)=factor*(vbrar(1,1)*vbrar(2,2)-vbrar(2,1)*vbrar(1,2))
      if(iprint.ge.1) then
         write(6,'(/)')
         write(6,'(12x,
     >         ''    n                    rslat                  '',
     >         ''rslatsq'')')
         write(6,'(12x,56(''=''))')
         write(6,'(12x,1i5,2x,4f12.5)')
     >        (n1,rslat_x(n1),rslat_y(n1),rslat_z(n1),rslatsq(n1),
     >         n1=1,nrslat)
         write(6,'(/)')
      endif
c
c     ================================================================
c     calculate kncut, the radius of k-space truncation sphere........
c     ----------------------------------------------------------------
      call getkncut(eta,vbrak(1,1),vbrak(1,2),vbrak(1,3),
     >              kncut,nm1,nm2,nm3)
c     ----------------------------------------------------------------
c
c     ================================================================
ctest
c     nm1=6
c     nm2=6
c     nm3=6
c     kncut=.98d0
c     kncut=.98d0/8.d0*6.2831853d0
ctest
c     generate the reciprocal space lattice vectors...................
c     ----------------------------------------------------------------
      call lattice(vbrak,kncut,nm1,nm2,nm3,
     >             knlat_x,knlat_y,knlat_z,knlatsq,nknlat,ipknlat_mad,
     >             iprint,istop)
c     ----------------------------------------------------------------
      if(iprint.ge.1) then
         write(6,'(  ''            nm1,nm2,nm3   = '',3i5)')nm1,nm2,nm3
         write(6,'(  ''            Kn cut radius = '',1f10.5)') kncut
         write(6,'(  ''            Number of Kn  = '',i5)') nknlat
         write(6,'(/)')
         write(6,'(12x,
     >         ''    n                    knlat                  '',
     >         ''knlatsq'')')
         write(6,'(12x,56(''=''))')
         write(6,'(12x,1i5,2x,4f12.5)')
     >        (n1,knlat_x(n1),knlat_y(n1),knlat_z(n1),knlatsq(n1),
     >         n1=1,nknlat)
      endif
c
c     ================================================================
      if(sname.eq.istop) then
         call fstop(sname)
      else
         return
      endif
c
      end
