c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cal_madelung_matrix(mynod,num_atoms,
     >                               bravais_lattice_in,
     >                               atom_posi_x_in,
     >                               atom_posi_y_in,
     >                               atom_posi_z_in,
     >                               madmat,
     >                               iprint,istop)
c     ================================================================
c
      implicit   none
c
!     include 'atom_param.h'
      include 'madelung.h'
c
      character  istop*32
      character  sname*32
      parameter (sname='cal_madelung_matrix')
c
      integer    mynod
      integer    num_atoms
c
      integer    iprint
      integer    i
      integer    imin
      integer    nrslat
      integer    nknlat
c
      integer    j,k
c
      real*8     bravais_lattice_in(9)
      real*8     atom_posi_x_in(num_atoms)
      real*8     atom_posi_y_in(num_atoms)
      real*8     atom_posi_z_in(num_atoms)
      real*8     bravais_lattice(9)
      real*8     atom_posi_x(num_atoms)
      real*8     atom_posi_y(num_atoms)
      real*8     atom_posi_z(num_atoms)
      real*8     madmat(num_atoms)
c
      real*8     a0
      real*8     dmin
      real*8     dmin2
      real*8     dmax
      real*8     dmax2
      real*8     dx
      real*8     dy
      real*8     dz
      real*8     dis2
      real*8     dis
      real*8     a1
      real*8     a2
      real*8     a3
      real*8     omegbra
      real*8     rslat_x(iprslat_mad)
      real*8     rslat_y(iprslat_mad)
      real*8     rslat_z(iprslat_mad)
      real*8     rslatsq(iprslat_mad)
      real*8     knlat_x(ipknlat_mad)
      real*8     knlat_y(ipknlat_mad)
      real*8     knlat_z(ipknlat_mad)
      real*8     knlatsq(ipknlat_mad)
      real*8     etainv
      real*8     eta
      real*8     etainv0
      real*8     sum
      real*8     one
      real*8     zero
      parameter (one=1.0d0)
      parameter (zero=0.0d0)
c
      parameter (etainv0=0.5d0)
c
!     if(max_atoms.lt.num_atoms) then
!        write(6,'('' CAL_MADELUNG_MATRIX:: max_atoms < num_atoms'',
!    >             2i5)')max_atoms,num_atoms
!        call fstop('cal_madelung_matrix')
!     endif
c     ----------------------------------------------------------------
      call mbeqa(bravais_lattice_in,bravais_lattice,9)
      call mbeqa(atom_posi_x_in,atom_posi_x,num_atoms)
      call mbeqa(atom_posi_y_in,atom_posi_y,num_atoms)
      call mbeqa(atom_posi_z_in,atom_posi_z,num_atoms)
c     ----------------------------------------------------------------
c
      a1=sqrt(bravais_lattice(1)*bravais_lattice(1)+
     >        bravais_lattice(2)*bravais_lattice(2)+
     >        bravais_lattice(3)*bravais_lattice(3) )
      a2=sqrt(bravais_lattice(4)*bravais_lattice(4)+
     >        bravais_lattice(5)*bravais_lattice(5)+
     >        bravais_lattice(6)*bravais_lattice(6) )
      a3=sqrt(bravais_lattice(7)*bravais_lattice(7)+
     >        bravais_lattice(8)*bravais_lattice(8)+
     >        bravais_lattice(9)*bravais_lattice(9) )
      a0=min(a1,a2,a3)
ctest
      etainv=etainv0+0.01*max(a1,a2,a3)/a0
      etainv=etainv*a0
ctest
c     ================================================================
c     change units so that both bravais_lattice and atom_posi_* are in
c     in the units of a0 = 1
c     ----------------------------------------------------------------
c     call dscal(9,one/a0,bravais_lattice,1)
c     call dscal(num_atoms,one/a0,atom_posi_x,1)
c     call dscal(num_atoms,one/a0,atom_posi_y,1)
c     call dscal(num_atoms,one/a0,atom_posi_z,1)
c     ----------------------------------------------------------------
c
c     ================================================================
c     obtain the lattice vectors for the big cell.....................
c     rslat, rslatsq, knlat, and knlatsq are in the units of a0 = 1...
c     ----------------------------------------------------------------
      call getstruc(mynod,bravais_lattice,omegbra,etainv,
     >              iprslat_mad,ipknlat_mad,
     >              rslat_x,rslat_y,rslat_z,rslatsq,nrslat,
     >              knlat_x,knlat_y,knlat_z,knlatsq,nknlat,
     >              iprint,istop)
c     ----------------------------------------------------------------
c
c     ================================================================
c     set up the Madelung matrix......................................
c     ----------------------------------------------------------------
      call madewd(nrslat,nknlat,mynod+1,num_atoms,one,
     >            rslat_x,rslat_y,rslat_z,rslatsq,
     >            knlat_x,knlat_y,knlat_z,knlatsq,
     >            atom_posi_x,atom_posi_y,atom_posi_z,
     >            etainv,omegbra,madmat,
     >            iprint,istop)
c     ----------------------------------------------------------------
c     ================================================================
      if(sname.eq.istop) then
         call fstop(sname)
      else
         return
      endif
      end
