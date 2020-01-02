c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setup_vorpol(my_atom,num_atoms,
     >                        atom_position_1,
     >                        atom_position_2,atom_position_3,
     >                        system_bravais,
     >                        lmax,clm,ngaussq,ngaussr,
     >                        rmt,omegint,dipint,rad,
     &                        ipvp,ipnode,ipcorn,ipedge,iprcrit,
     &                        gwwylm,grwylm,
     &                        ncrit,wylm,
     >                        iprint,istop)
c     ================================================================
c
      implicit   none
c
c     ****************************************************************
!      include    'vorpol.h'
c     ****************************************************************
c
      character  sname*32
      character  istop*32
c
      integer    my_atom
      integer    num_atoms
      integer    lmax
      integer    ngaussq
      integer    ngaussr
      integer    iprint
      integer    nvplane
      integer    ncorn
      integer    nedge
      integer    i

      integer ipvp,ipnode,ipcorn,ipedge,iprcrit
      integer ncrit
c
      real*8     atom_position_1(num_atoms)
      real*8     atom_position_2(num_atoms)
      real*8     atom_position_3(num_atoms)
      real*8     system_bravais(3,3)
      real*8     vplane(3,ipvp)
      real*8     rmt
      real*8     omegint
      real*8     clm((lmax+1)*(lmax+2)/2)
c
      real*8     dp2(3*ipvp)
      real*8     xp(3*ipvp)
      real*8     corn(ipcorn,3)
      real*8     dc2(ipcorn)
      integer    indxc(ipcorn)
      real*8     edge(3*ipedge)
      real*8     edgp(3*ipedge)
      real*8     tnode(ipnode)
      real*8     xgq(ngaussq)
      real*8     wgq(ngaussq)
      real*8     rcrit(iprcrit)
      real*8     runion(ipvp+ipedge+ipcorn)
      real*8     rad(num_atoms)
c
!      real*8     plm((iplcut+1)*(iplcut+2)/2)
       real*8     plm((lmax+1)*(lmax+2)/2)
c
      real *8    gwwylm(*),grwylm(*)

!      complex*16 sumfi(0:iplcut,3)
      complex*16 sumfi(0:lmax,3)
      complex*16 dipint(-1:1,2)
c
      complex*16 wylm

      parameter (sname='setup_vorpol')
c
c     ****************************************************************
c     For each sub-lattice calculates possible boundary planes used  
c     to define the Voronoi Polyhedron................................
c     ****************************************************************
c
      if(iprint.ge.1) then
	 write(6,'(/,'' SETUP_VORPOL:: lmax   ='',i5)') lmax
	 write(6,'(  ''                ngaussq='',i5)') ngaussq
	 write(6,'(  ''                ngaussr='',i5)') ngaussr
      endif
c
c     ================================================================
c     get possible VP boundary planes {vplane}........................
c     ----------------------------------------------------------------
      call setup_boundary(my_atom,num_atoms,
     >                    atom_position_1,
     >                    atom_position_2,
     >                    atom_position_3,
     >                    system_bravais(1,1),
     >                    system_bravais(1,2),
     >                    system_bravais(1,3),
     >                    vplane,ipvp,nvplane,rad)
c     ----------------------------------------------------------------
      if(iprint.ge.0) then
         write(6,'(/,'' SETUP_VORPOL:: Boundary Planes'')')
         write(6,'(  6x,i5,3f10.5)')
     >        (i,vplane(1,i),vplane(2,i),vplane(3,i),i=1,nvplane)
      endif
c
c     ================================================================
c     get trucation function {wylm}...................................
c
c     ================================================================
c     loop over sub-lattices..........................................
c     ================================================================
c     calculate the edges and corners of VP candidate planes..........
c     ----------------------------------------------------------------
      call polyhedron(vplane,dp2,nvplane,
     >                corn,dc2,ipcorn,ncorn,indxc,
     >                edge,edgp,ipedge,nedge)
c     ----------------------------------------------------------------
c
c     ================================================================
c     set muffin-tin radius to be inscribed sphere radius.............
c     only if rmt <= 0
      if(rmt.le.0.d0) rmt=sqrt(dp2(1))
c
c     ================================================================
c     get points at which w(l,m) has discontinuous derivatives........
      if(iprint.ge.0) then
         write(6,'(/,'' SETUP_VORPOL:: nbnd,ncorn,nedge:'',3i5)')
     >   nvplane,ncorn,nedge
      endif
c     ----------------------------------------------------------------
      call rcritpts(rcrit,ncrit,
     >              rmt*rmt,dc2(indxc(ncorn)),
     >              vplane,nvplane,
     >              dc2,ncorn,
     >              edgp,nedge,
     >              runion)
      if(ncrit.gt.iprcrit) then
	write(6,'(''iprcrit too small: iprcrit='',i4)') iprcrit
	stop'rcrit'
      endif
c     ----------------------------------------------------------------
c     ================================================================
c     write list of critical points if needed.........................
      if(iprint.ge.0) then
         write(6,'('' SETUP_VORPOL:: number of critical r-points'',
     >   '' for w(r) ='',i5)') ncrit
         write(6,'(''            i='',i4,'' rcrit='',1pd16.10)')
     >   (i,rcrit(i),i=1,ncrit)
      endif
c
c     ================================================================
c     get the step function wylm on gaussin integration mesh..........
c     ----------------------------------------------------------------
      call volvor(rcrit,ncrit,
     >            lmax,clm,
     >            xgq,wgq,ngaussq,ngaussr,
     >            xp,vplane,nvplane,
     >            corn(1,3),dc2,indxc,ncorn,
     >            edge,edgp,nedge,
     >            tnode,
     >            grwylm,gwwylm,wylm,omegint,plm,sumfi)
c     ----------------------------------------------------------------
      call inter_dip(lmax,ncrit,ngaussr,omegint,
     >            grwylm,gwwylm,wylm,dipint(-1,1),dipint(-1,2))
c
c     ================================================================
c     write out the calculated volumes................................
      if(iprint.ge.0) then
         write(6,'(  ''SETUP_VORPOL:Muffin-tin radius'',t40,''='',
     >   f14.8)') rmt
c     set circumscribing sphere radius................................
         write(6,'(  ''             Bounding sphere radius'',t40,''='',
     >   f14.8)') sqrt(dc2(indxc(ncorn)))
         write(6,'(  ''             Interstial region vol.'',t40,''='',
     >   f14.8)') omegint
      endif
c     ================================================================
      if(sname.eq.istop) then
         call fstop(sname)
      else
         return
      endif
c
      end
