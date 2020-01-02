!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdin_old_lsms(mynod,ipfile,cpath,                       &
                               info_table,info_evec,max_adjust,          &
                               systemid,cfspath,output_to_screen,istop,  &
                               node_print,print_instr,nprint,            &
                               nrelc,nrelv,mtasa,rmt0,vshift,vshiftl,    &
                               nspin,i_vdif,system_title,num_atoms,      &
                               igrid,ebot,etop,eitop,eibot,npts,         &
                               kelvin,nument,ialgorithm,iharris,         &
                               i_potwrite,movie,ntstep,tstep,nscf,       &
                               alpdv,alpma,alpev,ctq,j_ij,               &
                               mix_quant,mix_algor,                      &
                               ngaussr,ngaussq,                          &
                               num_adj_mix,atnm_adj_mix,msgpack_mix,     &
                               num_adj_csv,atnm_adj_csv,msgpack_csv,     &
                               etol,ptol,eftol,rmstol,iexch)
!     ================================================================
      use ErrorHandlerModule, only : ErrorHandler
      use ChemElementModule, only : getZtot, getZcor, getZsem, getZval
!
      implicit   none
!
      integer    max_adjust
      integer    mynod
!
      character  cpath*50
      character  sname*20
      character  istop*32
      character  text*72
      character  ipfile*100
      character  systemid*50
      character  cfspath*70
      character  output_to_screen*1
      character  system_title*70
      character  info_table*30
      character  info_evec*30
      character  atnm_adj_mix(max_adjust)*2
      character  atnm_adj_csv(max_adjust)*2
!
      logical    new_attribute
!
      integer    node_print
      integer    print_instr
      integer    nprint
      integer    nrelc
      integer    nrelv
      integer    mtasa
      integer    nspin
      integer    i_vdif
      integer    iexch
      integer    num_atoms
      integer    igrid
      integer    npts
      integer    kelvin
      integer    nument
      integer    ialgorithm
      integer    nscf
      integer    ntstep
      integer    mix_quant
      integer    mix_algor
      integer    iharris
      integer    i_potwrite
      integer    movie
      integer    ngaussr
      integer    ngaussq
      integer    num_adj_mix
      integer    num_adj_csv
      integer    node_adj_mix
      integer    node_adj_csv
      integer    mixalgo_adj
      integer    mixquan_adj
      integer    zadj
      integer    zti
      integer    zci
      integer    zsi
      integer    zvi
      integer    j_ij
      integer    i
!
      real*8     etol
      real*8     ptol
      real*8     eftol
      real*8     rmstol
!
      real*8     ebot
      real*8     etop
      real*8     eitop
      real*8     eibot
      real*8     alpdv
      real*8     alpma
      real*8     alpev
      real*8     ctq
      real*8     alpdv_adj
      real*8     alpma_adj
      real*8     alpev_adj
      real*8     ctq_adj
      real*8     rmt0,vshift,vshiftl
      real*8     zcore_adj
      real*8     zsemi_adj
      real*8     zvale_adj
      real*8     msgpack_mix(7,max_adjust)
      real*8     msgpack_csv(4,max_adjust)
      real*8     tstep
!
      real*8     tol
      real*8     zero
!
      parameter  (sname='rdin_old_lsms')
      parameter  (tol=1.0d-8)
      parameter  (zero=0.0d0)
!
!     write(6,'('' RDIN_ATOM_FT::  cpath: '',a100)')trim(cpath)
!     write(6,'('' RDIN_ATOM_FT:: ipfile: '',a100)')ipfile
!     write(6,'('' RDIN_ATOM_FT:: t_path: '',a100)')trim(cpath)//ipfile
!     ================================================================
!     open basic input file........................................... 
!     ================================================================
!     ----------------------------------------------------------------
      open(unit=10,file=trim(cpath)//ipfile,form='formatted',status='old')
!     ----------------------------------------------------------------
!     write(6,'('' RDIN_ATOM_FT::  Finished: Open Unit-10: '')')
      do i=1,11
         read(10,'(a72)') text
!        write(6,'(a72)') text
      enddo
!     ================================================================
!     read in the structure  identity.................................
      read(10,'(a)') text
      read(10,'(a)') systemid
!     ================================================================
!     read in the standard output target switch.......................
      read(10,'(a)') text
      read(10,'(a)') output_to_screen  
!     ================================================================
!     read in subroutine stop level
      read(10,'(a)') text
      read(10,'(a)') istop
!     write(6,'(a)') istop
!     ================================================================
!     read in print level for a particular node and rest of the nodes.
      read(10,'(a)') text
      read(10,*    ) node_print,print_instr,nprint
!     write(6,'(3i5)') node_print,print_instr,nprint
!     ================================================================
!     read in the number of atoms in the system.......................
!     ================================================================
      read(10,'(a)') text
      read(10,*    ) num_atoms
!     ================================================================
!     read in indices that cotrol relativity (core & valence) and
!     muffin tin/ASA switch
      read(10,'(a)') text
      read(10,*,err=1001) nrelc,nrelv,mtasa,rmt0,vshift,vshiftl
      if(nrelc.ne.0) then
         nrelc=10
      endif
      if(nrelv.ne.0) then
         nrelv=10
      endif
      if( mtasa .ge. 3 ) then
         mtasa = 3
      else if( mtasa .lt. -3 ) then
         call ErrorHandler('rdin_old_lsms','mtasa < -3',mtasa)
      endif
!     ================================================================
!     read spin polarization index....................................
      read(10,'(a)') text
      read(10,*    ) nspin,i_vdif,iexch
      if (nspin.lt.0 .or.nspin.ge.3) then
         call ErrorHandler('rdin_old_lsms','Wrong input for nspin',nspin)
      endif
!     ================================================================
!     read in a title to identify the system .........................
      read(10,'(a)') text
      read(10,'(a)') system_title
!     ================================================================
!     Read number of Gaussian points for r and theta integrations.....
      read(10,'(a)') text
      read(10,*    ) ngaussr,ngaussq
!     ================================================================
!     Read in name of the atom........................................
!     Read in the atom position vector................................
!     Read in cut off radius for the LIZ of the atom..................
!     Read in radius steps for lmax, lmax-1, lmax-2, lmax-3 ..........
!     ================================================================
!
!     ================================================================
!     read in names of info_table & info_evec files:..................
!     ================================================================
      read(10,'(a)') text
      read(10,'(2a30)')info_table,info_evec
!
!     ================================================================
!     read in parameters that control energy inregration:.............
!     ================================================================
!     igrid   : specifies energy contour  1=slow......................
!     igrid   : specifies energy contour  2=gaussian..................
!     igrid   : specifies energy contour  3=Don's Fermi function poles
!
!     for igrid =1...[Zero Temperature Calculations Only].............
!     ebot    : bottom of contour: real axis [may be mod. by semcor]..
!     etop    : top of contour: real axis [usually reset to chempot]..
!     eitop   : top    of contour on imaginary axis...................
!     eibot   : bottom of contour on imaginary axis...................
!     npts    : number of energy points per 0.1 ry....................
!     kelvin  : not used..............................................
!
!     for igrid =2...[Zero Temperature Calculations Only].............
!     ebot    : bottom of contour: real axis [may be mod. by semcor]..
!     etop    : top of contour: real axis [usually reset to chempot]..
!     eitop   : not used..............................................
!     eibot   : not used..............................................
!     npts    : number of Gaussian distributed energy points..........
!     kelvin  : not used..............................................
!
!     for igrid =3...[Finite Temperature Calculations Only]...........
!     ebot    : bottom of contour: real axis [may be mod. by semcor]..
!                                            [then reset in congauss].
!     etop    : top of contour on real axis [usually reset to chempot]
!     eitop   : not used..............................................
!     eibot   : not used..............................................
!     npts    : not used..............................................
!     kelvin  : temperature in kelvin..................................
!     nument  : # of gaussian points on elliptical contour for Entropy
!     ================================================================
      read(10,'(a)') text
      read(10, *   ) igrid,ebot,etop,eitop,eibot,npts,kelvin,nument,ialgorithm
      if(kelvin.ne.0 .and. nument.eq.0) then
         write(6,'('' RDIN_ATOM_FT:: kelvin.ne.0 .and. nument.eq.0'')')
         write(6,'('' RDIN_ATOM_FT:: warning entropy not calculated'')')
      endif
      if(igrid.gt.3) then
         call ErrorHandler('rdin_old_lsms','Code not set up for igrid>3',igrid)
      elseif(igrid.eq.3 .and. kelvin.lt.100) then
         call ErrorHandler('rdin_old_lsms','conmatsb being used: temp.lt.100')
      endif
!
!     ================================================================
!     read in controls for performing SCF calculation:................
!     ================================================================
!     nscf        : maximum number of scf iterations requested........
!     alpdv       : mixing parameter for chg. den. or potential.......
!     alpma       : mixing parameter for moment density...............
!     alpev       : mixing parameter for moment orientation...........
!     mix_quant   : mixing charge density[potential] => 0[1]..........
!     mix_algor   : mixing simple[DGAnderson] => 0[1].................
!     iharris = 0 : do not calculate harris energy....................
!     iharris = 1 : calculate harris energy using updated chem. potl..
!     iharris >=2 : calculate harris energy at fixed chem. potl.......
!     i_potwrite  : the number of iterations between potential writes.
!     movie   = 0 : no movie data will be written.....................
!             = 1 : movie data will be written........................
!     ctq         : coefficient of torque ............................
!     ================================================================
      read(10,'(a)') text
      read(10, *   ) nscf,alpdv,alpma,alpev,mix_quant,mix_algor,      &
                     iharris,i_potwrite,movie
!     ================================================================
!     check consistencey of these parameters..........................
!     ================================================================
      if(i_potwrite.le.0 .or. i_potwrite.gt.nscf) then
         i_potwrite=nscf
      endif
      if(mix_quant.ne.0 .and. mix_quant.ne.1) then
         call ErrorHandler('rdin_old_lsms',                           &
                           'Incorrect input data for mix_quant',mix_quant)
      else if(mix_algor.gt.2) then 
         call ErrorHandler('rdin_old_lsms',                           &
                           'Incorrect input data for mix_algor',mix_algor)
      endif
!     ================================================================
!     if calculating the harris energy make sure that mix_quant switch
!     is set mix the charge density...................................
!     ================================================================
      if(iharris.ne.0) then
         mix_quant=0
      endif
!
!     ================================================================
!     read in quantities that control Spin Dynamics :.................
!     ================================================================
!     nstep        : number of time steps [for SCF only : nstep=1 ]...
!     tstep        : time step........................................
!     etol         : over-ride scf convergence tolerence : energy ....
!     ptol         : over-ride scf convergence tolerence : pressure ..
!     eftol        : over-ride scf convergence tolerence : Fermi engy.
!     rmstol       : over-ride scf convergence tolerence : rmstol ....
!     ----------------------------------------------------------------
!     etol,ptol,eftol, & rmstol can be relaxed for SD.................
!     ================================================================
      read(10,'(a)') text
      read(10,*) ntstep,tstep,etol,ptol,eftol,rmstol
!
!     ================================================================
!     read controls for calculation of torque & exchange interactions.
!     ================================================================
!     ctq        : 
!     j_ij       : 
!     ================================================================
      read(10,'(a)') text
      read(10, *   ) ctq,j_ij
!     ================================================================
!     check consistencey of these parameters..........................
!     ================================================================
      if(ctq.lt.0.0d0) then
         call ErrorHandler('rdin_old_lsms',                           &
                           'Incorrect input data for ctq < 0.0',ctq)
      else if(j_ij.ne.0 .and. j_ij.ne.1) then
         call ErrorHandler('rdin_old_lsms',                           &
                           'Incorrect input data for j_ij <> 0, and <> 1',j_ij)
      endif
!
!     ================================================================
!     start to check if any parameters need to be adjusted for some 
!     atoms...........................................................
!     ================================================================
      read(10,'(a)') text
      read(10,*) num_adj_mix
      if(num_adj_mix.gt.max_adjust) then
         call ErrorHandler('rdin_old_lsms','num_adj_mix > max_adjust', &
                           num_adj_mix,max_adjust)
      endif
!     ================================================================
      read(10,'(a)')text
      do i=1,num_adj_mix
         read(10,'(1x,a2,9x,i5,6x,3f10.5,2(3x,i3),1f10.5)')            &
              atnm_adj_mix(i),node_adj_mix,                            &
              alpdv_adj,alpma_adj,alpev_adj,mixquan_adj,mixalgo_adj,ctq_adj
         msgpack_mix(1,i)=node_adj_mix
         msgpack_mix(2,i)=alpdv_adj
         msgpack_mix(3,i)=alpma_adj 
         msgpack_mix(4,i)=alpev_adj
         msgpack_mix(5,i)=mixquan_adj
         msgpack_mix(6,i)=mixalgo_adj
         msgpack_mix(7,i)=ctq_adj
      enddo
!     ================================================================
      new_attribute=.false.
      do while(.not. new_attribute)
         read(10,'(a)')text
         if(text( :19) .eq. ' No. of atoms whose') new_attribute=.true.
      enddo
      read(10,*)num_adj_csv
      if(num_adj_mix.gt.max_adjust) then
         call ErrorHandler('rdin_old_lsms','num_adj_csv > max_adjust', &
                           num_adj_csv,max_adjust)
      endif
      read(10,'(a)')text
      do i=1,num_adj_csv
         read(10,'(1x,a2,9x,i5,8x,3(2x,f5.2))')                        &
            atnm_adj_csv(i),node_adj_csv,zcore_adj,zsemi_adj,zvale_adj
!        -------------------------------------------------------------
         zti = getZtot(atnm_adj_csv(i))
         zci = getZcor(atnm_adj_csv(i))
         zsi = getZsem(atnm_adj_csv(i))
         zvi = getZval(atnm_adj_csv(i))
!        -------------------------------------------------------------
         if(abs(zcore_adj+zsemi_adj+zvale_adj-zti).gt. tol) then
            zadj = zcore_adj+zsemi_adj+zvale_adj
            call ErrorHandler('rdin_old_lsms',                        &
                              'Incorrect Zc+Zs+Zv adjustment',zadj,zti)
         endif
         msgpack_csv(1,i)=node_adj_csv
         msgpack_csv(2,i)=zcore_adj
         msgpack_csv(3,i)=zsemi_adj
         msgpack_csv(4,i)=zvale_adj
      enddo
      read(10,'(a)')text
      if(num_adj_csv.eq.0 ) then
         read(10,'(a)')text
      end if
!
      return
!     ================================================================
      if(igrid.gt.1 .or. npts.le.300) close(unit=10)
!     ================================================================
1001  continue
      call ErrorHandler('rdin_old_lsms',                              &
                        'After mtasa, missing rmt0, vshift, or vshiftl')
!
end subroutine rdin_old_lsms
