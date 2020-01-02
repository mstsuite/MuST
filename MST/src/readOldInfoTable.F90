!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readOldInfoTable(cfspath,systemid,num_atoms,nspin,   &
                                  alpdv,alpma,alpev,                  &
                                  info_table,info_evec,DataTable)
!     ================================================================
      use KindParamModule, only : IntKind, RealKind
      use ErrorHandlerModule, only : ErrorHandler
!
      use DataServiceCenterModule, only : createDataStorage,          &
                                          getDataStorage,             &
                                          RealType, RealMark,         &
                                          IntegerType, IntegerMark
!
      use ChemElementModule, only : getZtot
!
      use PublicTypeDefinitionsModule, only : InputTableStruct
!
      implicit   none
!
      type (InputTableStruct) :: DataTable
!
      integer (kind=IntKind), intent(in) :: num_atoms, nspin
      integer (kind=IntKind), parameter :: iplmax = 8
      integer (kind=IntKind) :: mynod
!
      character (len=70) ::  cfspath
      character (len=50) ::  systemid
      character (len=*)  :: info_table
      character (len=*)  :: info_evec
      character  atname*2
      character (len=2), allocatable :: name_type(:)
!
      integer    ntype
      integer, allocatable :: map(:)
      integer    n_per_type
      integer, allocatable :: n_p_t(:)
      integer    num_vacancy
      integer    num_types
      integer, allocatable :: num_atoms_in_type(:)
      integer, allocatable :: head_node(:)
      integer    lmax
      integer    n1,n2,i,n
!
      real*8     alpdv,alpma,alpev
      real*8     system_bravais(3,3)
      real*8     my_position(3)
      real*8, allocatable :: atom_position_1(:)
      real*8, allocatable :: atom_position_2(:)
      real*8, allocatable :: atom_position_3(:)
      real*8, allocatable :: system_evec(:,:)
      real*8, allocatable :: system_orient_mix(:)
      real*8, allocatable :: b_con(:,:)
      real*8     my_evec(3)
      real*8     my_orient_mix
      real*8     my_b_con(3)
      real*8     evec_tmp(3)
      real*8     o_mix_tmp
      real*8     b_con_tmp(3)
      real*8     pos_tmp
      real*8     rcirclu
      real*8     rsteps(iplmax+1)
      real*8, allocatable :: rad(:)
      real*8     radtmp
!
   allocate( name_type(1:num_atoms) )
   allocate( num_atoms_in_type(1:num_atoms*2) )
   allocate( head_node(1:num_atoms) )
   allocate( map(1:num_atoms), n_p_t(1:num_atoms) )
   allocate( atom_position_1(1:num_atoms) )
   allocate( atom_position_2(1:num_atoms) )
   allocate( atom_position_3(1:num_atoms) )
   allocate( rad(1:num_atoms) )
   allocate( system_evec(3,num_atoms) )
   allocate( system_orient_mix(num_atoms) )
   allocate( b_con(3,num_atoms) )
!
!     ================================================================
!     open info_table to read atom position information...............
!     ================================================================
!     write(6,*) 'RDIN_ATOM_FT:: call  rdin_info_table'            !new
!     ----------------------------------------------------------------
      call rdin_old_infotable(mynod,info_table,num_atoms,iplmax,      &
                              system_bravais,num_vacancy,num_types,   &
                              name_type,num_atoms_in_type,head_node,  &
                              atname,rcirclu,lmax,rsteps,my_position, &
                              ntype,map,n_p_t,atom_position_1,        &
                              atom_position_2,atom_position_3,rad)
!     ----------------------------------------------------------------
!
!     ================================================================
!     open info_evec to read moment orientation info:...[spin canting]
!     ----------------------------------------------------------------
      call rdin_old_infoevec(mynod,info_evec,num_atoms,               &
                             atom_position_1,atom_position_2,         &
                             atom_position_3,                         &
                             system_evec,system_orient_mix,           &
                             my_evec,my_orient_mix,b_con,my_b_con)
!     ----------------------------------------------------------------
!
!     ================================================================
!     move all the atoms with the real nodes to the front.............
!     ================================================================
      do n1=1,ntype
         if(map(n1).ne.n1) then
            do n2=n1+1,num_atoms
               if(map(n2).eq.n1) then
                  map(n2)=map(n1)
                  map(n1)=n1
!                 ====================================================
!                 move atom positions to front:.......................
!                 ====================================================
                  pos_tmp=atom_position_1(n1)
                  atom_position_1(n1)=atom_position_1(n2)
                  atom_position_1(n2)=pos_tmp
                  pos_tmp=atom_position_2(n1)
                  atom_position_2(n1)=atom_position_2(n2)
                  atom_position_2(n2)=pos_tmp
                  pos_tmp=atom_position_3(n1)
                  atom_position_3(n1)=atom_position_3(n2)
                  atom_position_3(n2)=pos_tmp
!                 ====================================================
!                 move moment orientations, mixing and b_con to front:
!                 ====================================================
                  do i=1,3
                     evec_tmp(i)=system_evec(i,n1)
                     system_evec(i,n1)=system_evec(i,n2)
                     system_evec(i,n2)=evec_tmp(i)
                  enddo
                  radtmp=rad(n1)
                  rad(n1)=rad(n2)
                  rad(n2)=radtmp
                  o_mix_tmp=system_orient_mix(n1)
                  system_orient_mix(n1)=system_orient_mix(n2)
                  system_orient_mix(n2)=o_mix_tmp
                  do i=1,3
                     b_con_tmp(i)=b_con(i,n1)
                     b_con(i,n1)=b_con(i,n2)
                     b_con(i,n2)=b_con_tmp(i)
                  enddo
                  goto 100
               endif
            enddo
         endif
100      continue
      enddo
!
!     ================================================================
!     send atom position and evec tables to other nodes...............
!     ================================================================
      n_per_type=n_p_t(1)
!
!     ****************************************************************
!     temporary code to write out the array: name_type ***************
!     ****************************************************************
      write(6,'(/,'' RDIN_ATOM_FT:: num_types,ntype,num_atoms'',3i5)') &
                                    num_types,ntype,num_atoms            !new
      do i=1,num_types                                                   !new
         write(6,'('' RDIN_ATOM_FT:: i,name_type:'',i5,3x,a2)')        &
                                     i,name_type(i)                      !new
      enddo                                                              !new
!     ****************************************************************
!
   n = 1
   DataTable%KeyName(n)='Bravais Vector'
   write(DataTable%KeyValue(n),'(3f19.15)')system_bravais(1:3,1)
   n = n + 1
   DataTable%KeyName(n)='Bravais Vector'
   write(DataTable%KeyValue(n),'(3f19.15)')system_bravais(1:3,2)
   n = n + 1
   DataTable%KeyName(n)='Bravais Vector'
   write(DataTable%KeyValue(n),'(3f19.15)')system_bravais(1:3,3)
   n = n + 1
   DataTable%KeyName(n)='Default Path to Potential Files'
   DataTable%KeyValue(n) = cfspath
   n = n + 1
   DataTable%KeyName(n)='Default Potential Input File Name'
   DataTable%KeyValue(n) = 'v_gopen_'//trim(systemid)
   n = n + 1
   DataTable%KeyName(n)='Default Potential Input File Form'
   DataTable%KeyValue(n) = '1'
   n = n + 1
   DataTable%KeyName(n)='Default Potential Output File Name'
   DataTable%KeyValue(n) = 'w_gopen_'//trim(systemid)
   n = n + 1
   DataTable%KeyName(n)='Default Potential Output File Form'
   DataTable%KeyValue(n) = '1'
   n = n + 1
   DataTable%KeyName(n)='Default Moment Direction'
   DataTable%KeyValue(n) = '0.00000000    0.00000000    1.00000000'
   n = n + 1
   DataTable%KeyName(n)='Default Constrain Field'
   DataTable%KeyValue(n) = '0.00000000    0.00000000    0.00000000'
   n = n + 1
   DataTable%KeyName(n)='Default Lmax-T matrix'
   write(DataTable%KeyValue(n),'(i5)')lmax
   n = n + 1
   DataTable%KeyName(n)='Default Lmax-Step Func'
   DataTable%KeyValue(n) = '0'
   n = n + 1
   DataTable%KeyName(n)='Default Lmax-Wave Func'
   DataTable%KeyValue(n) = '0'
   n = n + 1
   DataTable%KeyName(n)='Default Lmax-Potential'
   DataTable%KeyValue(n) = '0'
   n = n + 1
   DataTable%KeyName(n)='Default Lmax-Charge Den'
   DataTable%KeyValue(n) = '0'
   n = n + 1
   DataTable%KeyName(n)='Default LIZ # NN Shells'
   DataTable%KeyValue(n) = '8'
   n = n + 1
   DataTable%KeyName(n)='Default LIZ Shell Lmax'
   DataTable%KeyValue(n) = '3  3  3  3  2  2  1  1'
   n = n + 1
   DataTable%KeyName(n)='Default LIZ Cutoff Radius'
   write(DataTable%KeyValue(n),'(f12.5)')rcirclu
   n = n + 1
   DataTable%KeyName(n)='Default Rho  Mix Param.'
   write(DataTable%KeyValue(n),'(f12.5)')alpdv
   n = n + 1
   DataTable%KeyName(n)='Default Pot  Mix Param.'
   write(DataTable%KeyValue(n),'(f12.5)')alpdv
   n = n + 1
   DataTable%KeyName(n)='Default Mom  Mix Param.'
   write(DataTable%KeyValue(n),'(f12.5)')alpma
   n = n + 1
   DataTable%KeyName(n)='Default Evec Mix Param.' 
   write(DataTable%KeyValue(n),'(f12.5)')alpev
   n = n + 1
   DataTable%KeyName(n)='Default No. Rad Points ndivin'
   DataTable%KeyValue(n) = '501'
   n = n + 1
   DataTable%KeyName(n)='Default No. Rad Points ndivout'
   DataTable%KeyValue(n) = '0'
   n = n + 1
   DataTable%KeyName(n)='Default Integer Factor nmult'
   DataTable%KeyValue(n) = '1'
   n = n + 1
   DataTable%KeyName(n)='Default Screen Pot.'
   DataTable%KeyValue(n) = '4.0'
   n = n + 1
   DataTable%KeyName(n)='Default Lmax-Screen'
   DataTable%KeyValue(n) = '3'
   n = n + 1
   DataTable%KeyName(n)='Default Rcut-Screen'
   DataTable%KeyValue(n) = '4.8'
!
   DataTable%NumData = n
!
   stop 'readOldInfoTable is not fully implemented'
!
!  n = n + 1
!  do i=1,num_atoms
!     DataTable%KeyName(n)='Atom Index'
!     write(DataTable%KeyValue(n),'(i5)')i
!     n = n + 1
!  enddo
!
!  -------------------------------------------------------------------
!! call createDataStorage('Atomic Number',num_atoms,IntegerType)
!! call createDataStorage('Atomic Position',3*num_atoms,RealType)
!  -------------------------------------------------------------------
!! AtomicNumber => getDataStorage('Atomic Number',num_atoms,IntegerMark)
!! AtomPosition => getDataStorage('Atomic Position',num_atoms,3,RealMark)
!! AtomPositionX => AtomPosition(1:num_atoms,1)
!! AtomPositionY => AtomPosition(1:num_atoms,2)
!! AtomPositionZ => AtomPosition(1:num_atoms,3)
!  -------------------------------------------------------------------
!! do i=1,num_atoms
!!    AtomicNumber(i) =getZtot( name_type(i) )
!!    AtomPositionX(i) = atom_position_1(i)
!!    AtomPositionY(i) = atom_position_2(i)
!!    AtomPositionZ(i) = atom_position_3(i)
!! enddo
!! nullify( AtomicNumber )
!! nullify( AtomPosition, AtomPositionX, AtomPositionY, AtomPositionZ )
!
!! if (nspin > 2) then
!     ----------------------------------------------------------------
!!    call createDataStorage('Moment Direction',3*num_atoms,RealType)
!     ----------------------------------------------------------------
!!    Evec => getDataStorage('Moment Direction',num_atoms,3,RealMark)
!!    EvecX => Evec(1:num_atoms,1)
!!    EvecY => Evec(1:num_atoms,2)
!!    EvecZ => Evec(1:num_atoms,3)
!!    do i=1,num_atoms
!!       EvecX(i) = system_evec(1,i)
!!       EvecY(i) = system_evec(2,i)
!!       EvecZ(i) = system_evec(3,i)
!!    enddo
!!    nullify( Evec, EvecX, EvecY, EvecZ )
!
!     ----------------------------------------------------------------
!!    call createDataStorage('Constrain Field',3*num_atoms,RealType)
!     ----------------------------------------------------------------
!!    cf => getDataStorage('Constrain Field',num_atoms,3,RealMark)
!!    cfx => cf(1:num_atoms,1)
!!    cfy => cf(1:num_atoms,2)
!!    cfz => cf(1:num_atoms,3)
!!    do i=1,num_atoms
!!       cfx(i) = b_con(1,i)
!!       cfy(i) = b_con(2,i)
!!       cfz(i) = b_con(3,i)
!!    enddo
!!    nullify( cf, cfx, cfy, cfz )
!! endif
!
   deallocate( name_type )
   deallocate( num_atoms_in_type )
   deallocate( head_node )
   deallocate( map, n_p_t )
   deallocate( atom_position_1 )
   deallocate( atom_position_2 )
   deallocate( atom_position_3 )
   deallocate( rad )
   deallocate( system_evec )
   deallocate( system_orient_mix )
   deallocate( b_con )
!
end subroutine readOldInfoTable
