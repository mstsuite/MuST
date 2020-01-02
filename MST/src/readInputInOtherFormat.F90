!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readInputInOtherFormat(funit,DataTable)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : ZERO, CZERO
!
   use MPPModule, only : MyPE
!
   use PublicTypeDefinitionsModule, only : InputTableStruct
   use PublicParamDefinitionsModule, only : MaxLenFileName
!
   implicit none
!
!  type InputTableStruct
!     character (len=160) :: TableName
!     integer (kind=IntKind) :: TableIndex
!     integer (kind=IntKind) :: NumData
!     integer (kind=IntKind) :: UnitNumber
!     logical :: Open
!     character (len=50), pointer :: KeyName(:)
!     character (len=70), pointer :: KeyValue(:)
!     type (InputTableStruct), pointer :: next
!  end type InputTableStruct
!
   integer (kind=IntKind), intent(in) :: funit
!
   type (InputTableStruct), target :: DataTable
   type (InputTableStruct), pointer :: CurrentPtr
!
   integer (kind=IntKind), parameter :: max_adjust = 100
!
   character (len=50) :: cpath
   character (len=MaxLenFileName) :: ipfile
   character (len=70) :: text
   character (len=90) :: citpath
   character (len=90) :: ciepath
   character (len=70) :: cfspath
   character (len=50) :: systemid
   character (len=32) :: istop
   character (len=30) :: info_table
   character (len=30) :: info_evec
   character (len=1)  :: output_to_screen
   character (len=70) :: system_title
   character (len=2) :: atnm_adj_mix(max_adjust)
   character (len=2) :: atnm_adj_csv(max_adjust)
!
   integer (kind=IntKind) :: i, n
   integer (kind=IntKind) :: node_print
   integer (kind=IntKind) :: print_instr
   integer (kind=IntKind) :: nprint
   integer (kind=IntKind) :: nrelv
   integer (kind=IntKind) :: nrelc
   integer (kind=IntKind) :: mtasa
   integer (kind=IntKind) :: nspin
   integer (kind=IntKind) :: i_vdif
   integer (kind=IntKind) :: num_atoms
   integer (kind=IntKind) :: igrid
   integer (kind=IntKind) :: npts
   integer (kind=IntKind) :: kelvin
   integer (kind=IntKind) :: nument
   integer (kind=IntKind) :: j_ij
   integer (kind=IntKind) :: ialgorithm
   integer (kind=IntKind) :: iharris
   integer (kind=IntKind) :: i_potwrite
   integer (kind=IntKind) :: movie
   integer (kind=IntKind) :: nscf
   integer (kind=IntKind) :: ntstep
   integer (kind=IntKind) :: mix_algor
   integer (kind=IntKind) :: mix_quant
   integer (kind=IntKind) :: ngaussq
   integer (kind=IntKind) :: ngaussr
   integer (kind=IntKind) :: num_adj_mix
   integer (kind=IntKind) :: num_adj_csv
   integer (kind=IntKind) :: iexch
!
   real (kind=RealKind) :: vshift,vshiftl,rfix
   real (kind=RealKind) :: ebot
   real (kind=RealKind) :: etop
   real (kind=RealKind) :: eitop
   real (kind=RealKind) :: eibot
   real (kind=RealKind) :: tstep
   real (kind=RealKind) :: alpdv
   real (kind=RealKind) :: alpma
   real (kind=RealKind) :: alpev
   real (kind=RealKind) :: ctq
   real (kind=RealKind) :: msgpack_mix(7*max_adjust)
   real (kind=RealKind) :: msgpack_csv(4*max_adjust)
   real (kind=RealKind) :: etol
   real (kind=RealKind) :: ptol
   real (kind=RealKind) :: eftol
   real (kind=RealKind) :: rmstol
!
   do i=1,6
      read(funit,'(a)') text
   enddo
   read(funit,'(1a50)') cpath
   write(6,'('' LSMS_MAIN:: Path Standard in/out -   cpath: '',3x,a)') cpath
   read(funit,'(a)') text
   read(funit,'(1a50)') cfspath
   write(6,'('' LSMS_MAIN:: Path v(w)_gopen_xxx  - cfspath:'',3x,a50)')cfspath
   read(funit,'(a)') text
   read(funit,'(1a70)') ipfile
   write(6,'('' LSMS_MAIN:: Input file name               :'',3x,a70)')ipfile
!
!  -------------------------------------------------------------------
   call rdin_old_lsms(MyPE,ipfile,cpath,                              &
                      info_table,info_evec,max_adjust,                &
                      systemid,cfspath,output_to_screen,istop,        &
                      node_print,print_instr,nprint,                  &
                      nrelc,nrelv,mtasa,rfix,vshift,vshiftl,          &
                      nspin,i_vdif,system_title,num_atoms,            &
                      igrid,ebot,etop,eitop,eibot,npts,               &
                      kelvin,nument,ialgorithm,iharris,               &
                      i_potwrite,movie,                               &
                      ntstep,tstep,nscf,                              &
                      alpdv,alpma,alpev,ctq,j_ij,                     &
                      mix_quant,mix_algor,                            &
                      ngaussr,ngaussq,                                &
                      num_adj_mix,atnm_adj_mix,msgpack_mix,           &
                      num_adj_csv,atnm_adj_csv,msgpack_csv,           &
                      etol,ptol,eftol,rmstol,iexch)
!  -------------------------------------------------------------------
   n = 1
   DataTable%KeyName(n) = 'Current File Path'
   DataTable%KeyValue(n)=cpath
   n = n + 1
   DataTable%KeyName(n) = 'Current File Name'
   DataTable%KeyValue(n)=ipfile
   n = n + 1
   DataTable%KeyName(n) = 'Path to Info Table File'
   DataTable%KeyValue(n)=cpath
   n = n + 1
   DataTable%KeyName(n) = 'Info Table File Name'
   DataTable%KeyValue(n)=info_table
   n = n + 1
   DataTable%KeyName(n) = 'No. Iterations (> 0)'
   write(DataTable%KeyValue(n),'(i5)')nscf
   n = n + 1
   DataTable%KeyName(n) = 'Method of SCF Calculation'
   DataTable%KeyValue(n)='1'
   n = n + 1
   DataTable%KeyName(n) = 'Output to Screen (y/n)'
   DataTable%KeyValue(n)=output_to_screen
   n = n + 1
   DataTable%KeyName(n) = 'Output Level (>= -1)'
   write(DataTable%KeyValue(n),'(i5)')print_instr
   n = n + 1
   DataTable%KeyName(n) = 'Output Proc. ID (>= -1)'
   write(DataTable%KeyValue(n),'(i5)')node_print
   n = n + 1
   DataTable%KeyName(n) = 'Stop-at Routine Name'
   DataTable%KeyValue(n)=istop
   n = n + 1
   DataTable%KeyName(n) = 'No. Iter for Each Pot. Write'
   write(DataTable%KeyValue(n),'(i5)')i_potwrite
   n = n + 1
   DataTable%KeyName(n) = 'No. Iter for Each Movie'
   write(DataTable%KeyValue(n),'(i5)')movie
   n = n + 1
   DataTable%KeyName(n) = 'Calc. Harris Energy (H.E.)'
   write(DataTable%KeyValue(n),'(i5)')iharris
   n = n + 1
   DataTable%KeyName(n) = 'No. Gauss Pts. along r'
   DataTable%KeyValue(n)='80'
   n = n + 1
   DataTable%KeyName(n) = 'No. Gauss Pts. along theta'
   DataTable%KeyValue(n)='60'
   n = n + 1
   DataTable%KeyName(n) = 'Valence Band Bottom Est.'
   write(DataTable%KeyValue(n),'(f10.5)')ebot
   n = n + 1
   DataTable%KeyName(n) = 'Temperature Parameter (K)'
   write(DataTable%KeyValue(n),'(f12.4)')kelvin
   n = n + 1
   DataTable%KeyName(n) = 'Uniform Grid Parameters'
   DataTable%KeyValue(n)='10    10    10'
   n = n + 1
   DataTable%KeyName(n) = 'Energy (Ryd) Tol (> 0)'
   write(DataTable%KeyValue(n),'(f12.8)')etol
   n = n + 1
   DataTable%KeyName(n) = 'Potential Tol (> 0)'
   write(DataTable%KeyValue(n),'(f12.8)')ptol
   n = n + 1
   DataTable%KeyName(n) = 'Fermi Energy Tol (> 0)'
   write(DataTable%KeyValue(n),'(f12.8)')eftol
   n = n + 1
   DataTable%KeyName(n) = 'Other RMS Tol (> 0)'
   write(DataTable%KeyValue(n),'(f12.8)')rmstol
   n = n + 1
   DataTable%KeyName(n) = 'No. Atoms in System (> 0)'
   write(DataTable%KeyValue(n),'(i5)')num_atoms
   n = n + 1
   DataTable%KeyName(n) = 'Atom Position File Name'
   DataTable%KeyValue(n) = 'None'
   n = n + 1
   DataTable%KeyName(n) = 'Text Identification'
   DataTable%KeyValue(n)=systemid
   n = n + 1
   DataTable%KeyName(n) = 'Alloy System Description'
   DataTable%KeyValue(n)=system_title
   n = n + 1
   DataTable%KeyName(n) = 'Val. Electron Rel (>= 0)'
   if (nrelv == 0) then
       DataTable%KeyValue(n)='1'
   else
       DataTable%KeyValue(n)='0'
   endif
   n = n + 1
   DataTable%KeyName(n) = 'Core Electron Rel (>= 0)'
   if (nrelc == 0) then
       DataTable%KeyValue(n)='1'
   else
       DataTable%KeyValue(n)='0'
   endif
   n = n + 1
   DataTable%KeyName(n) = 'Potential Type (>= 0)'
   write(DataTable%KeyValue(n),'(i5)')mtasa
   n = n + 1
   DataTable%KeyName(n) = 'Exch-Corr. LDA Type (>= 0)'
   DataTable%KeyValue(n)='0'
   n = n + 1
   DataTable%KeyName(n) = 'Spin Index Param (>= 1)'
   write(DataTable%KeyValue(n),'(i5)')nspin
   n = n + 1
   DataTable%KeyName(n) = 'Canted Moment Torque Coef.'
   write(DataTable%KeyValue(n),'(f10.5)')ctq
   n = n + 1
   DataTable%KeyName(n) = 'Calculate J_ij (y/n)'
   if (j_ij == 0) then
      DataTable%KeyValue(n)='n'
   else
      DataTable%KeyValue(n)='y'
   endif
   n = n + 1
   DataTable%KeyName(n) = 'Read E-mesh from emeshs.inp'
   DataTable%KeyValue(n)='n'
   n = n + 1
   DataTable%KeyName(n) = 'Contour Type (>= 0)'
   if (igrid == 2) then
       DataTable%KeyValue(n)='0'
   else if (igrid == 3) then
       DataTable%KeyValue(n) = '3'
   else if (igrid == 1 .and. npts < 100) then
       DataTable%KeyValue(n) = '1'
   else
       DataTable%KeyValue(n) = '2'
   endif
   n = n + 1
   DataTable%KeyName(n) = 'Number of Contours (> 0)'
   DataTable%KeyValue(n)='1'
   n = n + 1
   DataTable%KeyName(n) = 'Energy Grid Type (>= 0)'
   if (igrid == 2) then
      DataTable%KeyValue(n)='1'
   else
      DataTable%KeyValue(n)='0'
   endif
   n = n + 1
   DataTable%KeyName(n) = 'No. Energy Grids'
   if (igrid == 2) then
      write(DataTable%KeyValue(n),'(i5)')npts
   else if (npts  > 300) then
      write(DataTable%KeyValue(n),'(i5)')npts-300
   else if (npts  >= 200) then
      DataTable%KeyValue(n) = '1'
   else if (npts >= 100 .and. npts < 200) then
      write(DataTable%KeyValue(n),'(i5)')max(npts-100,1)
   else
      write(DataTable%KeyValue(n),'(i5)')npts
   endif
   n = n + 1
   DataTable%KeyName(n) = 'Real Axis Bottom, erbot'
   write(DataTable%KeyValue(n),'(f10.5)')ebot
   n = n + 1
   DataTable%KeyName(n) = 'Real Axis Top, ertop'
   write(DataTable%KeyValue(n),'(f10.5)')etop
   n = n + 1
   DataTable%KeyName(n) = 'Imag Axis Bottom, eibot'
   write(DataTable%KeyValue(n),'(f10.5)')eibot
   n = n + 1
   DataTable%KeyName(n) = 'Imag Axis Top, eitop'
   write(DataTable%KeyValue(n),'(f10.5)')eitop
   n = n + 1
   DataTable%KeyName(n) = 'T-matrix inversion (>= 0)'
   write(DataTable%KeyValue(n),'(i5)')ialgorithm
   n = n + 1
   DataTable%KeyName(n) = 'M-matrix inversion (>= 0)'
   DataTable%KeyValue(n)='0'
   n = n + 1
   DataTable%KeyName(n) = 'No. Time Steps (>= 0)'
   write(DataTable%KeyValue(n),'(i5)')ntstep
   n = n + 1
   DataTable%KeyName(n) = 'Time Step'
   write(DataTable%KeyValue(n),'(f10.5)')tstep
   n = n + 1
   DataTable%KeyName(n) = 'Mixing quantity type'
   write(DataTable%KeyValue(n),'(i5)')mix_quant
   n = n + 1
   DataTable%KeyName(n) = 'Mixing algorithm'
   write(DataTable%KeyValue(n),'(i5)')mix_algor
   n = n + 1
   DataTable%KeyName(n) = 'Read K-mesh from kmeshs.inp'
   DataTable%KeyValue(n)='0'
   n = n + 1
   DataTable%KeyName(n) = 'Scheme to Generate K (>=0)'
   DataTable%KeyValue(n)='0'
   n = n + 1
   DataTable%KeyName(n) = 'No. K Meshs in IBZ (> 0)'
   DataTable%KeyValue(n)='1'
   n = n + 1
   DataTable%KeyName(n) = 'Kx, Ky, Kz Division (> 0)'
   DataTable%KeyValue(n)='16   16   16'
   n = n + 1
   DataTable%KeyName(n) = 'Symmetrize BZ Integration'
   DataTable%KeyValue(n)='1'
   DataTable%NumData = n
!
   citpath = trim(adjustl(cpath))//trim(info_table)
   ciepath = trim(adjustl(cpath))//trim(info_evec)
   allocate( DataTable%next )
   CurrentPtr => DataTable%next
   CurrentPtr%TableName = info_table
   CurrentPtr%UnitNumber = 10
   CurrentPtr%TableIndex = DataTable%TableIndex + 1
   allocate( CurrentPtr%KeyName(50), CurrentPtr%KeyValue(50) )
!  -------------------------------------------------------------------
   call readOldInfoTable(cfspath,systemid,num_atoms,nspin,            &
                         alpdv,alpma,alpev,citpath,ciepath,CurrentPtr)
!  -------------------------------------------------------------------
   CurrentPtr%Open = .false.
!
end subroutine readInputInOtherFormat
