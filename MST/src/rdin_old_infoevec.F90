!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdin_old_infoevec(mynod,info_evec,num_atoms,           &
                                   atom_position_1,atom_position_2,     &
                                   atom_position_3,                     &
                                   system_evec, system_orient_mix,      &
                                   my_evec,my_orient_mix,b_con,my_b_con)
!     ==================================================================
!
      implicit   none
!
      character  sname*20
      character  text*72
      character  info_evec*80
!
      integer    mynod
      integer    num_atoms
      integer    nodedum
      integer    i
      integer    at_num
!
      real*8     system_evec(3,num_atoms)
      real*8     system_orient_mix(num_atoms)
      real*8     b_con(3,num_atoms)
      real*8     my_evec(3)
      real*8     my_orient_mix
      real*8     my_b_con(3)
      real*8     atom_position_1(num_atoms)
      real*8     atom_position_2(num_atoms)
      real*8     atom_position_3(num_atoms)
      real*8     x
      real*8     y
      real*8     z
      real*8     tol
      real*8     zero
      real*8     one
      real*8     two
!
      parameter  (sname='rdin_info_evec')
      parameter  (zero=0.0d0)
      parameter  (one=1.0d0)
      parameter  (two=2.0d0)
      parameter  (tol=1.0d-3)
!
!     ******************************************************************
!     Reads in atom position information etc that is contained in the --
!     "info_evec"-file created by "bigcell" :...........................
!     NB:---------------------------------------------------------------
!     If the file name is read in as "default" code or if it is empty
!     does not read moment
!     orientations from "info_evec" : rather it sets system_evec=[9,0.0] 
!     and system_orient_mix=-1.0 to signal that default values are to be
!     used : then evec is taken from "v_gopen"-file and alpev is taken -
!     from the standard input file [i_b000n_"system_id"]................
!     ******************************************************************
!
!     ==================================================================
!     Open info_evec if necessary: otherwise set default values:........
!     ==================================================================
      if(info_evec.ne.'default'.and.info_evec(1:1).ne.' ') then
!        ---------------------------------------------------------------
         open(unit=15,file=info_evec,form='formatted',status='old')
!        ---------------------------------------------------------------
      else
         do i=1,num_atoms
            system_evec(1,i)=two 
            system_evec(2,i)=zero 
            system_evec(3,i)=zero 
            system_orient_mix(i)=-one
!
!           ============================================================
!           b_con(j,i) i=1,3 i=1,num_atoms: constraining magnetic field
!           defined in the right handed coordinate system:-
!           B_z-axis : along evec ......................................
!           B_y-axis : perpendicular to evec and in the plane containing 
!                      global z-axis and evec...........................
!           B_x-axis : perpendicular to z and x axes defined above......
!
!                  B_z : [sin(t)cos(p), sin(t)sin(p), cos(t)]
!                  B_y : [cos(t)cos(p), cos(t)sin(p),-sin(t)]
!                  B_x : [B_y]x[B_z]
!
!           The contraining magnetic field is initialized to [0,0,0]....
!           If B_z=0 then the constaining field is transverse to evec...
!           ============================================================
            b_con(1,i)=zero
            b_con(2,i)=zero
            b_con(3,i)=zero
         enddo
         go to 1000
      endif
!
!     ==================================================================
!     Read in head information : Not used by this subroutine:...........
!     ==================================================================
      do i=1,11
         read(15,'(a)') text
      enddo
!
!     ==================================================================
!     Read evec, and orientation mixing parameter of all atoms in system
!     ==================================================================
!     ------------------------------------------------------------------
      b_con(1:3,1:num_atoms) = zero
!     ------------------------------------------------------------------
      do i=1,num_atoms
         read(15,*) at_num,nodedum,x,y,z,                               &
                   system_evec(1,i),system_evec(2,i),system_evec(3,i),  &
                   system_orient_mix(1),b_con(1,i),b_con(2,i)
!        ===============================================================
!        Test to make sure position of atom is same as from info_table..
!        ===============================================================
         if(abs(x-atom_position_1(i)) + abs(y-atom_position_2(i))  +    &
            abs(z - atom_position_3(i) )  .gt. tol ) then
            write(6,'('' RDIN_INFO_EVEC:: Trouble with atom pos.'')')
            write(6,'('' RDIN_INFO_EVEC::   info_table atom pos.'',3f12.5)') &
                  atom_position_1(i),atom_position_2(i),atom_position_3(i)
            write(6,'('' RDIN_INFO_EVEC::    info_evec atom pos.'',3f12.5)') &
                  x,y,z
         endif
      enddo
!
!     ==================================================================
!     Close info_evec:.................................................
!     ==================================================================
!     ------------------------------------------------------------------
      close(unit=15)
!     ------------------------------------------------------------------
!
!     ==================================================================
!     Pick out evec, and orientational mixing param. for node 0.........
!     ==================================================================
1000  continue
      do i=1,3
         my_evec(i)=system_evec(i,1)
         my_b_con(i)=b_con(i,1)
      enddo
      my_orient_mix=system_orient_mix(1)
!
      end subroutine rdin_old_infoevec
