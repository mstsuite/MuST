!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rdin_old_infotable(mynod,info_table,num_atoms,iplmax,      &
                                    system_bravais,num_vacancy,num_types,   &
                                    name_type,num_atoms_in_type,head_node,  &
                                    atname,rcirclu,lmax,rsteps,my_position, &
                                    ntype,map,n_p_t,atom_position_1,        &
                                    atom_position_2,atom_position_3,rad)
!     ==================================================================
      use ErrorHandlerModule, only : ErrorHandler
      use ChemElementModule, only : getName, getZcor, getZsem, getZval
      implicit   none
!
      integer, intent(in) :: num_atoms, iplmax
!
      character  sname*20
      character  text*72
      character  info_table*30
      character  atdum*2
      character  atname*2
      character  name_type(num_atoms)*2
!
      integer    mynod
      integer    msgtyp

      integer    ntype
      integer    map(num_atoms)
      integer    n_p_t(num_atoms)
      integer    num_vacancy
      integer    num_types
      integer    num_atoms_in_type(2,num_atoms)
      integer    head_node(num_atoms)
      integer    lm0
      integer    lmax
      integer    n1
      integer    n2
      integer    iflag
      integer    node
      integer    i
      integer    at_num
      integer    ztc
      integer    zts
      integer    ztv
!
      real*8     system_bravais(3,3)
      real*8     my_position(3)
      real*8     atom_position_1(num_atoms)
      real*8     atom_position_2(num_atoms)
      real*8     atom_position_3(num_atoms)
      real*8     rcircdum
      real*8     rcirclu
      real*8     rstdum(iplmax+1)
      real*8     rsteps(iplmax+1)
      real*8     rad(num_atoms)
!
      real*8     tol
      real*8     zero
!
      parameter  (sname='rdin_old_infotable')
      parameter  (tol=1.0d-8)
      parameter  (zero=0.0d0)
!
!     ******************************************************************
!     Reads in atom position information etc that is contained in the --
!     "info_table"-file created by "bigcell" :..........................
!     ******************************************************************
!
!     ==================================================================
!     Open info_table:..................................................
!     ==================================================================
!     ------------------------------------------------------------------
      open(unit=11,file=info_table,form='formatted',status='old')
!     ------------------------------------------------------------------
!
!     ==================================================================
!     Read in Bravais lattice specifying the Big Unit Cell:.............
!     ==================================================================
      read(11,'(a)') text
      read(11,'(a)') text
      read(11,'(a)') text
      do n2=1,3
         read(11,* ) (system_bravais(n1,n2),n1=1,3)
      enddo
      read(11,'(a)')text
      read(11,'(a)')text
      read(11,'(a)')text
      read(11,'(a)')text
      read(11,'(a)')text
!
!     ==================================================================
!     Read in position (etc) of atom on node zero:......................
!     ==================================================================
      read(11,*)  at_num,node,                                          &
                  my_position(1),my_position(2),my_position(3),         &
                  lmax,rcirclu,(rsteps(i),i=1,4),rad(1)
!     ==================================================================
!     Find out the name of the atom on node zero:.......................
!     ==================================================================
!     ------------------------------------------------------------------
!     call getaninfo(atname,at_num,ztc,zts,ztv)
      atname = getName(at_num)
      ztc = getZcor(at_num)
      zts = getZsem(at_num)
      ztv = getZval(at_num)
!     ------------------------------------------------------------------
      write(6,'('' RDIN_INFO_TABLE:: at_num,atname,ztc,zts,ztv'',       &
     &       3x,i3,3x,a2,3x,4i3)') at_num,atname,ztc,zts,ztv,node
!
!     ==================================================================
!     Set array that controls stepping down of lmax LIZ shells:.......
!     ==================================================================
      do i=5,lmax+1
         rsteps(i)=rcirclu
      enddo
      if(lmax.gt.4) then
         do i=lmax,lmax-3,-1
            rsteps(i)=rsteps(i-lmax+4)
         enddo
         do i=1,lmax-4
            rsteps(i)=1.d-6
         enddo
      endif
      if(node.ne.0) then
         call ErrorHandler('rdin_old_infotable','incorrect node',node)
      else if(lmax.gt.iplmax .or. lmax.lt.0) then
         call ErrorHandler('rdin_old_infotable','incorrect lmax:',lmax)
      endif
!
!     ==================================================================
!     Initialise (for atom 1; node 0) the atom_postion array:.........
!     ==================================================================
      atom_position_1(1)=my_position(1)
      atom_position_2(1)=my_position(2)
      atom_position_3(1)=my_position(3)
!     ==================================================================
!     Initialise counting and indexing of atom types :................
!     ==================================================================
      num_types=1
      name_type(num_types)=atname
      num_atoms_in_type(1,num_types)=1
      num_atoms_in_type(2,num_types)=1
      head_node(num_types)=0
      if (atname.eq.'Va') then
         num_vacancy=1
      else
         num_vacancy=0
      endif
      ntype=1
      map(1)=ntype
      n_p_t(1)=1
!
!     ==================================================================
!     Read in position (etc) of remaining atoms :.....................
!     ==================================================================
      do n1=2,num_atoms
!        ===============================================================
!        Read in position (etc) of remaining atom:......................
!        ===============================================================
         read(11,*)  at_num,node,                                       &
         atom_position_1(n1),atom_position_2(n1),atom_position_3(n1),   &
         lm0,rcircdum,(rstdum(i),i=1,4),rad(n1)
!        ===============================================================
!        Find out the name of this atom:................................
!        ===============================================================
!        ---------------------------------------------------------------
!        call getaninfo(atdum,at_num,ztc,zts,ztv)
         atdum = getName(at_num)
         ztc = getZcor(at_num)
         zts = getZsem(at_num)
         ztv = getZval(at_num)
!        ---------------------------------------------------------------
!
!        ===============================================================
!        Set array that controls stepping down of lmax LIZ shells:....
!        ===============================================================
         do i=5,lm0+1
            rstdum(i)=rcirclu
         enddo
         if(lm0.gt.4) then
            do i=lm0,lm0-3,-1
               rstdum(i)=rstdum(i-lm0+4)
            enddo
            do i=1,lm0-4
               rstdum(i)=1.d-6
            enddo
         endif
!        ===============================================================
!        Check if this is a new type; increment type indexing:........
!        ===============================================================
         if(atdum.eq.'Va') then
            num_vacancy=num_vacancy+1
         endif
         iflag=0
         do n2=1,num_types
            if(atdum.eq.name_type(n2)) iflag=n2
         enddo
         if(iflag.eq.0) then
            num_types=num_types+1
            name_type(num_types)=atdum
            num_atoms_in_type(1,num_types)=1
            num_atoms_in_type(2,num_types)=1
            if(node.gt.0) then
               head_node(num_types)=node
            else
               call ErrorHandler('rdin_old_infotable','incorrect node',node,n1)
            endif
         else
            num_atoms_in_type(1,iflag)=num_atoms_in_type(1,iflag)+1
            if(node.gt.0) then
               num_atoms_in_type(2,iflag)=num_atoms_in_type(2,iflag)+1
            endif
         endif
         if(node.gt.0) then
            ntype=ntype+1
            n_p_t(ntype)=1
            if(node.ne.ntype-1) then
               write(6,'('' RDIN_ATOM_FT:: incorrect node:'',2i5)') node,n1
               stop 'Error'
            elseif(lm0.gt.iplmax .or. lm0.lt.0) then
               write(6,'('' RDIN_ATOM_FT:: incorrect lm0:'',1i5)')lm0
               stop 'Error'
            endif
            map(n1)=ntype
         else   ! node .gt.0
            n2=-node
            if(n2.gt.ntype) then
               call ErrorHandler('rdin_old_infotable','Trouble with nme_type', &
                                 ntype,n2)
            endif
            map(n1)=n2
            n_p_t(n2)=n_p_t(n2)+1
         endif  ! node .gt.0
      enddo
      msgtyp=msgtyp+4
!
!     ==================================================================
!     Close info_table:.................................................
!     ==================================================================
!     ------------------------------------------------------------------
      close(unit=11)
!     ------------------------------------------------------------------
!
!     **************************************************************
!     temporary code to write out the array: name_type *************
!     **************************************************************
!     write(6,'(/,'' RDIN_INFO_TABLE:: num_types,ntype,num_atoms:'',
!    >               3i5)')            num_types,ntype,num_atoms
!     do i=1,num_types
!        write(6,'('' RDIN_INFO_TABLE:: i,name_type:'',i5,3x,a2)')
!    >                                  i,name_type(i)
!      enddo
!     **************************************************************
!
      end subroutine rdin_old_infotable
