!*********************************************************************
!* MODULE NAME    : InputModule                                      *
!*                                                                   *
!* VERSION NUMBER : 1.0                                              *
!* LAST MODIFIED  : MARCH 14, 2004                                   *
!*                                                                   *
!* VERSION NUMBER : 1.1                                              *
!* LAST MODIFIED  : JUNE 6, 2016                                     *
!* NOTES: Added function getKeyIndexValue                            *
!*        This combines calling getKeyValue and getKeyIndex into one *
!*        subroutine call.                                           *
!*        getKeyIndex should be considered obsolete                  *
!*                                                                   *
!* The module for open, readi, and close formatted input data files  *
!* The data in each input file are stored in an InputTable.          *
!* Given the key name, the data can be accessed using following      *
!* public routines:                                                  *
!*                                                                   *
!*********************************************************************
module InputModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : ZERO, CZERO
   use MPPModule, only : MyPE, bcastMessage, syncAllPEs
   use PublicTypeDefinitionsModule, only : InputTableStruct
   use PublicParamDefinitionsModule, only : StandardInputFile
!
public :: initInput,       &
          endInput,        &
          openInputFile,   &
          readInputData,   &
          closeInputFile,  &
          getTableName,    &
          getTableIndex,   &
          printKeyNames,   &
          printKeyValues,  &
          getNumData,      &
          getNumKeys,      &
          getNumKeyValues, &
          getKeyValue,     &
          getKeyIndex,     &
          getKeyIndexValue, &
          isKeyExisting
!
   interface getKeyValue
      module procedure getKeyValue_str0, getKeyValue_str1, &
                       getKeyValue_str2, getKeyValue_str3
      module procedure getKeyValue_int0, getKeyValue_int1, &
                       getKeyValue_int2, getKeyValue_int3
      module procedure getKeyValue_real0, getKeyValue_real1, &
                       getKeyValue_real2, getKeyValue_real3
      module procedure getKeyValue_cmplx0, getKeyValue_cmplx1, &
                       getKeyValue_cmplx2, getKeyValue_cmplx3
   end interface
!
   interface getKeyIndexValue
      module procedure getKeyIndexValue_str0, getKeyIndexValue_str1, &
                       getKeyIndexValue_str2, getKeyIndexValue_str3
      module procedure getKeyIndexValue_int0, getKeyIndexValue_int1, &
                       getKeyIndexValue_int2, getKeyIndexValue_int3
      module procedure getKeyIndexValue_real0, getKeyIndexValue_real1, &
                       getKeyIndexValue_real2, getKeyIndexValue_real3
   end interface
!
private
!
   logical :: Initialized = .false.
!
   integer (kind=IntKind) :: NumInputTables
!  type (InputTableStruct), target, save :: DataTable
   type (InputTableStruct), pointer :: DataTable
   type (InputTableStruct), pointer :: CurrentPtr
!
   integer (kind=IntKind), parameter :: InputPE = 0
   integer (kind=IntKind), parameter :: MessageType = 12100
!
   integer (kind=IntKind) :: pflag, fflag, status
   integer (kind=IntKind), parameter :: MaxNumData = 200
!
   character (len=80) :: file_path
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initInput()
!  ===================================================================
   use DefaultParamModule, only : initDefaultParam
   implicit none
!
   NumInputTables = 0
!  DataTable%UnitNumber = -1
!  DataTable%Open = .false.
!
   call initDefaultParam()
!
   Initialized = .true.
!
   end subroutine initInput
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endInput()
!  ===================================================================
   use DefaultParamModule, only : endDefaultParam
   implicit none
   integer (kind=IntKind) :: i
   type (InputTableStruct), pointer :: next
!
   CurrentPtr=>DataTable
   do i=1, NumInputTables
      deallocate(CurrentPtr%KeyName, CurrentPtr%KeyValue)
      next=>CurrentPtr%next
      deallocate(CurrentPtr)
      CurrentPtr=>next
   enddo
   nullify(CurrentPtr)
   nullify(next)
!
   call endDefaultParam()
!
   Initialized = .false.
   NumInputTables = 0
!
   end subroutine endInput
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultInput() result(s)
!  ===================================================================
!
   implicit none
   character (len=80) :: text
   character (len=80) :: file_name
   character (len=160) :: s
!
   integer (kind=IntKind) :: jdx
!
   if (.not.Initialized) then
      call ErrorHandler('getDefaultInput','need to call initInput first')
   endif
!
   if (MyPE == InputPE) then
      pflag=0; fflag=0
      do 
         read(5,'(a)',iostat=status)text
         if (status < 0) then
            write(6,'(a)')'File name is not found in the standard input file'
            call ErrorHandler('getDefaultInput','Invalid input data file')
         else if (text(1:17) == 'Current File Path') then
            jdx=index(text,'::')
            file_path = adjustl(text(jdx+2:))
            pflag=1
         else if (text(1:17) == 'Current File Name') then
            jdx=index(text,'::')
            file_name = adjustl(text(jdx+2:))
            fflag=1
         else if (pflag==1 .and. fflag==1) then
            s=trim(file_path)//file_name
            exit
         endif
      enddo
      close(5)
      open(unit=5,file=s)
      rewind(5)
   endif
   call bcastMessage(s,InputPE)
   call SyncAllPEs()
!
   end function getDefaultInput
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine openInputFile(funit,fname)
!  ===================================================================
   implicit none
   character (len=*),intent(in) :: fname
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: ios, i
   logical :: FileExist
!
   if (.not.Initialized) then
      call ErrorHandler('openInputFile','need to call initInput first')
   endif
!
   if (MyPE == InputPE .and. funit /= 5) then
      inquire(file=fname,exist=FileExist)
      if (FileExist) then
         open(unit=funit,file=fname,form='formatted',status='old',    &
              iostat=ios,action='read')
         if (ios > 0) then
            write(6,'(a,i5)')'File unit = ',funit
            call ErrorHandler('openInputFile','iostat > 0',fname)
         endif
      else
         call ErrorHandler('openInputFile','File does not exist',fname)
      endif
   endif
!
   if (NumInputTables == 0) then
      allocate(DataTable)
      CurrentPtr => DataTable
   else
      CurrentPtr => DataTable
      do i=2,NumInputTables
         CurrentPtr => CurrentPtr%next
      enddo
      allocate(CurrentPtr%next)
      CurrentPtr => CurrentPtr%next
   endif
!
   NumInputTables = NumInputTables+1
   CurrentPtr%TableName = fname
   CurrentPtr%UnitNumber = funit
   CurrentPtr%TableIndex = NumInputTables
   CurrentPtr%Open = .true.
   nullify(CurrentPtr%next)
!
   call SyncAllPEs()
   end subroutine openInputFile
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readInputData(funit,tindex)
!  ===================================================================
   implicit none
   logical :: opd, unit_found
   character (len=220) :: text
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(out) :: tindex
!
   integer (kind=IntKind) :: n,i,jdx
   integer (kind=IntKind) :: InputFormat
   integer (kind=IntKind) :: msgbuf(2)
   integer (kind=IntKind), parameter :: NumDataInOldFormat = 100
!
   interface
      subroutine readInputInOtherFormat(funit,DataTable)
         use PublicTypeDefinitionsModule, only : InputTableStruct
         use KindParamModule, only : IntKind
         integer (kind=IntKind), intent(in) :: funit
         type (InputTableStruct), target :: DataTable
      end subroutine readInputInOtherFormat
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('readInputData','need to call initInput first')
   endif
!
   unit_found = .false.
!
   if (NumInputTables == 0) then
      allocate(DataTable)
   endif
   CurrentPtr=>DataTable
   do i=1,NumInputTables
      if (CurrentPtr%UnitNumber==funit .and. CurrentPtr%Open) then
         unit_found = .true.
         exit
      else if (i < NumInputTables) then
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   if (.not. unit_found) then
      if (funit==5) then  ! For unit=5, openInputFile may not be called
         if (NumInputTables > 0) then
            allocate(CurrentPtr%next)
            CurrentPtr => CurrentPtr%next
         endif
         NumInputTables = NumInputTables+1
if (.false.) then
         CurrentPtr%TableName = getDefaultInput()
else
         CurrentPtr%TableName = 'stdin'
endif
         CurrentPtr%Open = .true.
         CurrentPtr%UnitNumber = funit
         CurrentPtr%TableIndex = NumInputTables
         nullify(CurrentPtr%next)
      else
         call ErrorHandler('readInputData','file unit is invalid',funit)
      endif
   else if (MyPE == InputPE .and.  funit /= 5) then
      inquire(unit=funit,opened=opd)
      if (.not.opd) then
         call ErrorHandler('readInputData','file unit is not opened',funit)
      endif
   endif
!
   tindex = CurrentPtr%TableIndex
!
   if (MyPE == InputPE) then
!     ================================================================
      read(funit,'(a)',iostat=status)text
      text=adjustl(text)
      if (text(1:1) == "*") then
          InputFormat = 0
      else
          InputFormat = 1
      endif
   endif
!  -------------------------------------------------------------------
   call bcastMessage(InputFormat,InputPE)
!  -------------------------------------------------------------------
!
   allocate(CurrentPtr%KeyName(MaxNumData), CurrentPtr%KeyValue(MaxNumData))
!
   if (MyPE == InputPE) then
      if (funit /= 5) then
         rewind(funit)
      endif
      if (InputFormat == 1) then
         n = 0
!        =============================================================
!        In case of funit=5, since one line has been read, it needs to be
!        processed to extract data from the line
!        =============================================================
         if (funit == 5 .and. text(1:1) /= '#' .and. text(1:1) /= '!') then
            jdx = index(text,'::')
            if (jdx > 1) then
               if (text(jdx+2:jdx+2) /= ':') then
                  n = n+1
                  CurrentPtr%KeyName(n)=adjustl(text(1:jdx-1))
                  CurrentPtr%KeyValue(n)=adjustl(text(jdx+2:))
               endif
            endif
         endif
!        -------------------------------------------------------------
         call readInputInStandardFormat(funit,n)
!        -------------------------------------------------------------
         CurrentPtr%NumData = n
      else
!        -------------------------------------------------------------
         call readInputInOtherFormat(funit,CurrentPtr)
!        -------------------------------------------------------------
         NumInputTables = NumInputTables+1  ! add one more table since
                                            ! info_table is also stored
      endif
   endif
!
   call bcastMessage(CurrentPtr%NumData,InputPE)
   call bcastMessage(CurrentPtr%KeyName,CurrentPtr%NumData,InputPE)
   call bcastMessage(CurrentPtr%KeyValue,CurrentPtr%NumData,InputPE)
   if (InputFormat == 0) then
      if (MyPE == InputPE) then
         msgbuf(1) = CurrentPtr%next%NumData
         msgbuf(2) = CurrentPtr%next%TableIndex
      else
         allocate(CurrentPtr%next)
         CurrentPtr%next%UnitNumber = 10
         CurrentPtr%next%Open = .false.
      endif
      CurrentPtr=>CurrentPtr%next
!     ----------------------------------------------------------------
      call bcastMessage(msgbuf,2,InputPE)
!     ----------------------------------------------------------------
      CurrentPtr%NumData = msgbuf(1); CurrentPtr%TableIndex = msgbuf(2)
!     ----------------------------------------------------------------
      call bcastMessage(CurrentPtr%TableName,InputPE)
      call bcastMessage(CurrentPtr%KeyName,CurrentPtr%NumData,InputPE)
      call bcastMessage(CurrentPtr%KeyValue,CurrentPtr%NumData,InputPE)
!     ----------------------------------------------------------------
   endif
   call SyncAllPEs()
!
   end subroutine readInputData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   recursive subroutine readInputInStandardFormat(funit,n)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind), intent(inout) :: n
   integer (kind=IntKind) :: jdx
!
   character (len=220) :: text
   character (len=80) :: stmp1, stmp2
   character (len=160) :: s
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   do 
      read(funit,'(a)',iostat=status)text
      if (status < 0) then
         exit
      endif
      text=adjustl(text)
      if (text(1:1) == '#' .or. text(1:1) == '!') then
         cycle
      else
         jdx = index(text,'::')
         if (jdx > 1) then
            if (text(jdx+2:jdx+2) /= ':') then
               stmp1 = adjustl(text(1:jdx-1))
               stmp2 = adjustl(text(jdx+2:))
               if (nocaseCompare(stmp1,'Include File')) then
                  s=trim(file_path)//stmp2
                  open(unit=funit+1,file=s,status='old',form='formatted')
!                 ----------------------------------------------------
                  call readInputInStandardFormat(funit+1,n)
!                 ----------------------------------------------------
                  close(unit=funit+1)
               else if (n >= MaxNumData) then
!                 ----------------------------------------------------
                  call ErrorHandler('readInputInStardardFormat',      &
                                    'Number of input parameters exceeds the upper limit',MaxNumData)
!                 ----------------------------------------------------
               else
                  n=n+1
                  CurrentPtr%KeyName(n)=adjustl(text(1:jdx-1))
                  CurrentPtr%KeyValue(n)=adjustl(text(jdx+2:))
               endif
            endif
         endif
      endif
   enddo
!
   end subroutine readInputInStandardFormat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine closeInputFile(funit)          
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: funit
   integer (kind=IntKind) :: i
!
   CurrentPtr=>DataTable
   do i=1,NumInputTables
      if (CurrentPtr%UnitNumber==funit .and. CurrentPtr%Open) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
   if (.not.associated(CurrentPtr)) then
      call ErrorHandler('closeInputData','file unit is invalid',funit)
   else
      CurrentPtr%Open = .false.
   endif
!
   if (MyPE == InputPE .and. funit /= 5) then
      close(unit=funit)
   endif
   call SyncAllPEs()
!
   end subroutine closeInputFile
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTableName(id) result(s)
!  ===================================================================
   implicit none
   character (len=160) :: s
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getTableName','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         s = CurrentPtr%TableName
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
   end function getTableName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getTableIndex(s) result(id)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: s
   integer (kind=IntKind) :: id
   integer (kind=IntKind) :: i
!
   id = 0
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableName == adjustl(s)) then
         id = CurrentPtr%TableIndex
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!  if (id == 0) then
!     call ErrorHandler('getTableIndex','invalid table name',s)
!  endif
!
   end function getTableIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printKeyNames(id)
!  ===================================================================
   implicit none
   character (len=50) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('printKeyName','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   key=' '
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) /= key) then
         key=CurrentPtr%KeyName(i)
         write(6,'(a15,a)')'Key Name: ',key
      endif
   enddo
!
   end subroutine printKeyNames
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printKeyValues(id)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('printKeyName','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   do i=1,CurrentPtr%NumData
      write(6,'(a15,a)')'Key Name : ',trim(CurrentPtr%KeyName(i))
      write(6,'(a15,a)')'Key Value: ',trim(CurrentPtr%KeyValue(i))
   enddo
!
   end subroutine printKeyValues
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumData(id) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: i
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getNumKeyValues','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   n=CurrentPtr%NumData
!
   end function getNumData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumKeys(id) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: i, j
   integer (kind=1), allocatable :: flag(:)
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getNumKeyValues','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   allocate(flag(1:CurrentPtr%NumData))
   flag(1:CurrentPtr%NumData)=1_1
   n=CurrentPtr%NumData
   do i=1,CurrentPtr%NumData-1
      if (flag(i) == 1_1) then
         do j=i+1,CurrentPtr%NumData
            if (flag(j) == 0_1) then
               cycle
            else if(CurrentPtr%KeyName(i) == &
                    CurrentPtr%KeyName(j)) then
               n = n-1
               flag(j)=0_1
            endif
         enddo
      endif
   enddo
   deallocate(flag)
!
   end function getNumKeys
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isKeyExisting(id,key) result(found)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i
!
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('isKeyExisting','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   LOOP_i0: do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit LOOP_i0
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo LOOP_i0
!
   found = .false.
   LOOP_i1: do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         found = .true.
         exit LOOP_i1
      endif
   enddo LOOP_i1
!
   end function isKeyExisting
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumKeyValues(id,key) result(n)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: i
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getNumKeyValues','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   n=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         n = n+1
      endif
   enddo
!
   end function getNumKeyValues
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_str0(id,key,value,default_param) result(status)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   character (len=*), intent(out) :: value
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         value=trim(adjustl(CurrentPtr%KeyValue(i)))
         found = .true.
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then
         if (getDefaultValue(key,value) == 0) then 
            status = 0
         endif
      else
         value(:) = ' '
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_str0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_str1(id,key,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: id
   character (len=*), intent(out) :: value(n)
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         value(m)=trim(adjustl(CurrentPtr%KeyValue(i)))
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:n)(:) = ' '
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_str1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_str2(id,key,k,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: k,n
   integer (kind=IntKind), intent(in) :: id
   character (len=*), intent(out) :: value(k,n)
   character (len=80) :: s
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         s=adjustl(CurrentPtr%KeyValue(i))
         read(s,'(80a)',iostat=status)value(1:k,m)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:k,1:n)(:) = ' '
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_str2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_str3(id,key,k,value,default_param) result(status)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id, k
   character (len=*), intent(out) :: value(k)
   character (len=80) :: s
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         s=adjustl(CurrentPtr%KeyValue(i))
         read(s,'(80a)',iostat=status)value(1:k)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then
         if (getDefaultValue(key,k,value) == 0) then 
            status = 0
         endif
      else
         value(1:k)(:) = ' '
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_str3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_int0(id,key,value,default_param) result(status)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: value
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)value
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then
         if (getDefaultValue(key,value) == 0) then 
            status = 0
         endif
      else
         value = 0
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_int0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_int1(id,key,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: value(n)
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(m)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:n) = 0
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_int1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_int2(id,key,k,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k,n
   integer (kind=IntKind), intent(out) :: value(k,n)
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(1:k,m)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:k,1:n) = 0
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_int2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_int3(id,key,k,value,default_param) result(status)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k
   integer (kind=IntKind), intent(out) :: value(k)
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(1:k)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then
         if (getDefaultValue(key,k,value) == 0) then 
            status = 0
         endif
      else
         value(1:k) = 0
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_int3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_real0(id,key,value,default_param) result(status)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   real (kind=RealKind), intent(out) :: value
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)value
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then
         if (getDefaultValue(key,value) == 0) then 
            status = 0
         endif
      else
         value = ZERO
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_real0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_real1(id,key,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: n
   real (kind=RealKind), intent(out) :: value(n)
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(m)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:n) = ZERO
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_real1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_real2(id,key,k,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k,n
   real (kind=RealKind), intent(out) :: value(k,n)
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(1:k,m)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:k,1:n) = ZERO
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_real2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_real3(id,key,k,value,default_param) result(status)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k
   real (kind=RealKind), intent(out) :: value(k)
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(1:k)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then
         if (getDefaultValue(key,k,value) == 0) then 
            status = 0
         endif
      else
         value(1:k) = ZERO
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_real3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_cmplx0(id,key,value,default_param) result(status)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   complex (kind=CmplxKind), intent(out) :: value
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)value
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then
         if (getDefaultValue(key,value) == 0) then 
            status = 0
         endif
      else
         value = CZERO
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_cmplx0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_cmplx1(id,key,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: n
   complex (kind=CmplxKind), intent(out) :: value(n)
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(m)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:n) = CZERO
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_cmplx1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_cmplx2(id,key,k,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k,n
   complex (kind=CmplxKind), intent(out) :: value(k,n)
   integer (kind=IntKind) :: i, m, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(1:k,m)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      value(1:k,1:n) = CZERO
      status = 1
   else
      status = 0
   endif
!
   end function getKeyValue_cmplx2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyValue_cmplx3(id,key,k,value,default_param) result(status)
!  ===================================================================
!
!  Given the data table ID: id
!        the data description: key
!        the data dimension: k
!
!  Returns the data: value(1:k)
!  ===================================================================
   use DefaultParamModule, only : getDefaultValue
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k
   complex (kind=CmplxKind), intent(out) :: value(k)
   integer (kind=IntKind) :: i, status
   logical, optional, intent(in) :: default_param
   logical :: found, dp
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   found = .false.
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)value(1:k)
         if (status < 0) then
            call ErrorHandler('getKeyValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         endif
         found = .true.
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyValue','Key not found',key)
      if (present(default_param)) then
         dp = default_param
      else
         dp = .true.
      endif
      status = 1
      if (dp) then 
         if (getDefaultValue(key,k,value) == 0) then 
            status = 0
         endif
      else
         value(1:k) = CZERO
      endif
   else
      status = 0
   endif
!
   end function getKeyValue_cmplx3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndex(id,indKey,key,k_index,n) result(status)
!  ===================================================================
!
!  Given the data table ID: id
!        the data list index description: indKey
!        the data list name description: key
!        the data list size: n
!        the data list index array: k_index(1:n)  
!
!  Returns the data: value(1:k)
!
!  Purpose: One call of this routine allows processing two input lines
!        which next to one another:
!        indKey  :: value of indKey
!        key     :: value of key
!  The "value of indKey" will be the value of an array element of index
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: indKey
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: k_index(n)
   integer (kind=IntKind) :: i, m, k, status
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndex','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index(1:n) = 0
!
   status = 0
   m = 0
   LOOP_i: do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == indKey) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)k
         if (status /= 0) then
            call ErrorHandler('getKeyIndex','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else if (k > n .or. k < 1) then
            call ErrorHandler('getKeyIndex','k is out of range',k)
         endif
      else if (CurrentPtr%KeyName(i) == key) then
         m = m + 1
         k_index(k)=m
      endif
   enddo LOOP_i
!
   end function getKeyIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_str0(id,key,k_index,value) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: k_index
   character (len=*), intent(out) :: value
   character (len=80) :: s
   integer (kind=IntKind) :: i, status, p1, p2
   logical :: found
!
   interface
      function getTokenPosition(k,s,n) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: n
         integer :: p
      end function getTokenPosition
   end interface
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(:) = ' '
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         s=trim(adjustl(CurrentPtr%KeyValue(i)))
         p1 = getTokenPosition(1,s)
         p2 = getTokenPosition(2,s)
         read(s(p1:p2-1),*,iostat=status)k_index
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else
            value = s(p2:)
            found = .true.
         endif
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_str0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_str1(id,key,k_index,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: k_index(n)
   character (len=*), intent(out) :: value(n)
   character (len=80) :: s
   integer (kind=IntKind) :: i, m, status, p1, p2, j
   logical :: found
!
   interface
      function getTokenPosition(k,s,n) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: n
         integer :: p
      end function getTokenPosition
   end interface
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!  
   k_index = 0
   found = .false.
   value(1:n)(:) = ' '
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         s=trim(adjustl(CurrentPtr%KeyValue(i)))
         p1 = getTokenPosition(1,s)
         p2 = getTokenPosition(2,s)
         read(s(p1:p2-1),*,iostat=status)j
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else if (j < 1 .or. j > n) then
            call ErrorHandler('getKeyIndexValue','Index value out of range',j)
         else
            value(m) = s(p2:)
            found = .true.
         endif
         k_index(j) = m
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_str1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_str2(id,key,k,k_index,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: k,n
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: k_index(n)
   character (len=*), intent(out) :: value(k,n)
   character (len=80) :: s
   integer (kind=IntKind) :: i, m, status, p1, p2, j
   logical :: found
!
   interface
      function getTokenPosition(k,s,n) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: n
         integer :: p
      end function getTokenPosition
   end interface
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:k,1:n)(:) = ' '
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         s=adjustl(CurrentPtr%KeyValue(i))
         p1 = getTokenPosition(1,s)
         p2 = getTokenPosition(2,s)
         read(s(p1:p2-1),*)j
         read(s(p2:),'(80a)',iostat=status)value(1:k,m)
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else if (j < 1 .or. j > n) then
            call ErrorHandler('getKeyIndexValue','Index value out of range',j)
         else
            found = .true.
         endif
         k_index(j) = m
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_str2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_str3(id,key,k,k_index,value) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id, k
   integer (kind=IntKind), intent(out) :: k_index
   character (len=*), intent(out) :: value(k)
   character (len=80) :: s
   integer (kind=IntKind) :: i, status, p1, p2
   logical :: found
!
   interface
      function getTokenPosition(k,s,n) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: n
         integer :: p
      end function getTokenPosition
   end interface
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:k)(:) = ' '
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         s=adjustl(CurrentPtr%KeyValue(i))
         p1 = getTokenPosition(1,s)
         p2 = getTokenPosition(2,s)
         read(s(p1:p2-1),*)k_index
         read(s(p2:),'(80a)',iostat=status)value(1:k)
         if (status < 0) then
            call ErrorHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else
            found = .true.
         endif
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_str3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_int0(id,key,k_index,value) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: k_index
   integer (kind=IntKind), intent(out) :: value
   integer (kind=IntKind) :: i, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value = 0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)k_index,value
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else
            found = .true.
         endif
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_int0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_int1(id,key,k_index,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: k_index(n)
   integer (kind=IntKind), intent(out) :: value(n)
   integer (kind=IntKind) :: i, m, status, j
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:n) = 0
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)j,value(m)
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else if (j < 1 .or. j > n) then
            call ErrorHandler('getKeyIndexValue','Index value out of range',j)
         else
            found = .true.
         endif
         k_index(j) = m
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_int1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_int2(id,key,k,k_index,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k,n
   integer (kind=IntKind), intent(out) :: k_index(n)
   integer (kind=IntKind), intent(out) :: value(k,n)
   integer (kind=IntKind) :: i, m, status, j
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:k,1:n) = 0
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)j,value(1:k,m)
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else if (j < 1 .or. j > n) then
            call ErrorHandler('getKeyIndexValue','Index value out of range',j)
         else
            found = .true.
         endif
         k_index(j) = m
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_int2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_int3(id,key,k,k_index,value) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k
   integer (kind=IntKind), intent(out) :: k_index
   integer (kind=IntKind), intent(out) :: value(k)
   integer (kind=IntKind) :: i, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:k) = 0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)k_index,value(1:k)
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else
            found = .true.
         endif
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_int3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_real0(id,key,k_index,value) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(out) :: k_index
   real (kind=RealKind), intent(out) :: value
   integer (kind=IntKind) :: i, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value = ZERO
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)k_index,value
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else
            found = .true.
         endif
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_real0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_real1(id,key,k_index,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: k_index(n)
   real (kind=RealKind), intent(out) :: value(n)
   integer (kind=IntKind) :: i, m, status, j
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:n) = ZERO
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)j,value(m)
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else if (j < 1 .or. j > n) then
            call ErrorHandler('getKeyIndexValue','Index value out of range',j)
         else
            found = .true.
         endif
         k_index(j) = m
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_real1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_real2(id,key,k,k_index,value,n) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k,n
   integer (kind=IntKind), intent(out) :: k_index(n)
   real (kind=RealKind), intent(out) :: value(k,n)
   integer (kind=IntKind) :: i, m, status, j
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:k,1:n) = ZERO
   m=0
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         m = m+1
         read(CurrentPtr%KeyValue(i),*,iostat=status)j,value(1:k,m)
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else if (j < 1 .or. j > n) then
            call ErrorHandler('getKeyIndexValue','Index value out of range',j)
         else
            found = .true.
         endif
         k_index(j) = m
      endif
      if (m == n) then
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_real2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getKeyIndexValue_real3(id,key,k,k_index,value) result(status)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: key
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: k
   integer (kind=IntKind), intent(out) :: k_index
   real (kind=RealKind), intent(out) :: value(k)
   integer (kind=IntKind) :: i, status
   logical :: found
!
   if (id < 1 .or. id > NumInputTables) then
      call ErrorHandler('getKeyIndexValue','invalid table index',id)
   endif
!
   CurrentPtr => DataTable
   do i=1,NumInputTables
      if (CurrentPtr%TableIndex == id) then
         exit
      else
         CurrentPtr => CurrentPtr%next
      endif
   enddo
!
   k_index = 0
   found = .false.
   value(1:k) = ZERO
   do i=1,CurrentPtr%NumData
      if (CurrentPtr%KeyName(i) == key) then
         read(CurrentPtr%KeyValue(i),*,iostat=status)k_index,value(1:k)
         if (status < 0) then
            call WarningHandler('getKeyIndexValue','Invalid value',          &
                              CurrentPtr%KeyValue(i))
         else
            found = .true.
         endif
         exit
      endif
   enddo
!
   if (.not.found) then
!     call WarningHandler('getKeyIndexValue','Key not found',key)
      status = 1
   else
      status = 0
   endif
!
   end function getKeyIndexValue_real3
!  ===================================================================
end module InputModule
