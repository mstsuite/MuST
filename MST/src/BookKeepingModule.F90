module BookKeepingModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicParamDefinitionsModule, only : MaxLenFileName
!
public :: initBookKeeping,  &
          endBookKeeping,   &
          writeHeadLine,    &
          insertColumn,     &
          isBookKeepingNew, &
          writeRow,         &
          getHeadLineValue, &
          getLastColumnValue
!
   interface writeRow
      module procedure writeRow0, writeRow1, writeRow2
   end interface
!
private
   character (len=MaxLenFileName) :: fname
!
   logical :: FileExist = .false.
!
   integer (kind=IntKind), parameter :: funit = 105
   integer (kind=IntKind), parameter :: MaxColumns = 10
   integer (kind=IntKind), parameter :: MaxHeadLines = 30
   integer (kind=IntKind) :: num_head_lines = 0
   integer (kind=IntKind) :: num_rows = 0
   integer (kind=IntKind) :: num_columns = 0
   integer (kind=IntKind) :: num_columns_old = 0
!
   character (len=12) :: ColumnTitle(MaxColumns)
   character (len=12) :: ColumnValue(MaxColumns)
   character (len=40) :: HeadLineTitle(MaxHeadLines)
   character (len=60) :: HeadLineValue(MaxHeadLines)
!
   integer (kind=IntKind) :: ColumnLength(MaxColumns)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initBookKeeping(fpath,short_label,myid)
!  ===================================================================
   use StringModule, only : initString, setString, endString,         &
                            getNumTokens, readToken
!
   implicit none
!
   character (len=*), intent(in) :: fpath
   character (len=*), intent(in) :: short_label
!
   integer (kind=IntKind), intent(in) :: myid
!
   character (len=6)  :: cnode
   character (len=80) :: text, text_sav
!
   integer (kind=IntKind), parameter :: offset = 10**(len(cnode)-1)
   integer (kind=IntKind) :: i, id, status, rflag
!
   num_head_lines=0
   num_rows = 0
   num_columns = 0
   num_columns_old = 0
   FileExist = .false.
!
   write(cnode,'(i6)')myid+offset
!
   if (myid >= offset) then
!      ---------------------------------------------------------------
       call WarningHandler('initBookKeeping','myid is too long',myid)
!      ---------------------------------------------------------------
   endif
   cnode(:1)='n'
!
   fname = trim(adjustl(fpath))//'k_'//trim(cnode)//'_'//adjustl(short_label)
!
   inquire(file=trim(fname),exist=FileExist)
   if (FileExist) then
      open(unit=funit,file=trim(fname),form='formatted',position='rewind', &
           status='old')
!
      rflag = 0
!     ----------------------------------------------------------------
      call initString(80)
!     ----------------------------------------------------------------
      LOOP_1: do
         read(funit,'(a)',iostat=status)text
         if (status < 0 .or. rflag > 2) then
            exit LOOP_1
         else if (text(1:25) == '=========================') then
            rflag = rflag + 1
         else if (text(1:25) == '-------------------------') then
            rflag = rflag + 1
         else if (rflag == 0) then
            num_head_lines = num_head_lines + 1
            id = index(text,':')
            if (id < 2 .or. id > len(text)) then
               call WarningHandler('initBookKeeping','Invalid headline',text)
               exit LOOP_1
            endif
            HeadLineTitle(num_head_lines) = adjustl(text(:id-1))
            HeadLineValue(num_head_lines) = adjustl(text(id+1:))
         else if (rflag == 1) then
!           ----------------------------------------------------------
            call setString(text)
!           ----------------------------------------------------------
            num_columns_old = getNumTokens()
            if (num_columns_old > MaxColumns) then
               call WarningHandler('initBookKeeping',                      &
                                   'Number of columns exceed upper limit',text)
               exit LOOP_1
            endif
            do i=1, num_columns_old
               call readToken(i,ColumnTitle(i),ColumnLength(i))
            enddo
         else if (rflag == 2) then
            text_sav = text
         endif
      enddo LOOP_1
!     ----------------------------------------------------------------
      call endString()
!     ----------------------------------------------------------------
!
      if (rflag /= 2) then
         close(funit)
         open(unit=funit,file=trim(fname),form='formatted',status='unknown')
         FileExist = .false.
         num_head_lines = 0
         num_columns_old = 0
         return
      else
!        -------------------------------------------------------------
         call initString(text_sav)
!        -------------------------------------------------------------
         if (getNumTokens() == num_columns_old) then
            do i=1, num_columns_old
               call readToken(i,ColumnValue(i),ColumnLength(i))
            enddo
         else if (getNumTokens() == num_columns_old+1) then
            do i=1, num_columns_old
               call readToken(i+1,ColumnValue(i),ColumnLength(i))
            enddo
         else
            call WarningHandler('initBookKeeping','Invalid line',text)
            do i=1, getNumTokens()
               call readToken(i,ColumnValue(i),ColumnLength(i))
            enddo
            do i=getNumTokens()+1,num_columns_old
               ColumnValue(i)='Invalid'
               ColumnLength(i)=ColumnLength(getNumTokens())
            enddo
         endif
!        -------------------------------------------------------------
         call endString()
!        -------------------------------------------------------------
!
         close(funit)
         open(unit=funit,file=trim(fname),form='formatted',status='old', &
              position='append')
      endif
   else
      open(unit=funit,file=trim(fname),form='formatted',status='unknown')
   endif
!
   end subroutine initBookKeeping
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endBookKeeping()
!  ===================================================================
   implicit none
!
   close(unit=funit)
!
   num_head_lines=0
   num_rows = 0
   num_columns = 0
   num_columns_old = 0
   FileExist = .false.
!
   end subroutine endBookKeeping
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isBookKeepingNew() result(bkn)
!  ===================================================================
   implicit none
!
   logical :: bkn
!
   bkn = .not.FileExist
!
   end function isBookKeepingNew
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getHeadLineValue(title) result(v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: title
!
   character (len=60) :: v
!
   integer (kind=IntKind) :: i
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t             
      end function nocaseCompare  
   end interface
!
   do i=1,num_head_lines
      if ( nocaseCompare( trim(HeadLineTitle(i)), trim(title) ) ) then
         v = HeadLineValue(i)
         return
      endif
   enddo
!
   call ErrorHandler('getHeadLineValue','Headline is not found',title)
!
   end function getHeadLineValue
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLastColumnValue(title) result(v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: title
!
   character (len=60) :: v
!
   integer (kind=IntKind) :: i
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t             
      end function nocaseCompare  
   end interface
!
   do i = 1, num_columns_old
      if ( nocaseCompare( trim(adjustl(title)),                       &
                          trim(adjustl(ColumnTitle(i))) ) ) then
         v = ColumnValue(i)
         return
      endif
   enddo
!
   call ErrorHandler('getLastColumnValue','title is not found',title)
!
   end function getLastColumnValue
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeHeadLine(title,value)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: title
   character (len=*), intent(in) :: value
!
   character (len=30) :: tt
   character (len=48) :: tv
!
   if (FileExist) then
      return
   endif
!
   tt = adjustl(title)
   tv = adjustl(value)
   write(funit,'(a,'': '',a)')tt,tv
   num_head_lines = num_head_lines + 1
!
   if (num_head_lines > MaxHeadLines) then
      call ErrorHandler('writeHeadLine',                              &
                        'Number of headlines exceed upper limit',MaxHeadLines)
   endif
   HeadLineTitle(num_head_lines) = tt
   HeadLineValue(num_head_lines) = tv
   call FlushFile(funit)
!
   end subroutine writeHeadLine
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine insertColumn(title,value)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   implicit none
!
   character (len=*), intent(in) :: title
   character (len=*), intent(in) :: value
!
   integer (kind=IntKind) :: nt, nv, n, m
!
   num_columns = num_columns + 1
!
   if (num_columns > MaxColumns) then
      call ErrorHandler('insertColumn',                               &
                        'Number of columns exceed the upper limit',MaxColumns)
   endif
!
   nt = len_trim(title)
   nv = len_trim(value)
   n = max(nt,nv)
!
   ColumnLength(num_columns) = n
!
   m = (n - nt)/2
   if (m > 0) then
      ColumnTitle(num_columns)(:m) = ' '
   endif
   ColumnTitle(num_columns)(m+1:) = title
!
   m = (n - nv)/2
   if (m > 0) then
      ColumnValue(num_columns)(:m) = ' '
   endif
   ColumnValue(num_columns)(m+1:) = value
!
   end subroutine insertColumn
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeRow0()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, n
!
   if (.not.FileExist) then
      write(funit,'(87(''=''))')
      do i = 1, num_columns
         n = ColumnLength(i)
         write(funit,'(1x,a,1x,$)')ColumnTitle(i)(1:n)
      enddo
      write(funit,'(/,87(''-''))')
      FileExist = .true.
   endif
!
   do i = 1, num_columns-1
      n = ColumnLength(i)
      write(funit,'(1x,a,1x,$)')ColumnValue(i)(1:n)
   enddo
   n = ColumnLength(num_columns)
   write(funit,'(1x,a)')ColumnValue(num_columns)(1:n)
   call FlushFile(funit)
!
   num_columns_old = num_columns
   num_columns = 0
   num_rows = num_rows + 1
!
   end subroutine writeRow0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeRow1(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: c
!
   integer (kind=IntKind) :: i, n
!
   if (.not.FileExist) then
      write(funit,'(87(''=''))')
      n = ColumnLength(1)
      write(funit,'(2x,a,1x,$)')ColumnTitle(1)(1:n)
      do i = 2, num_columns
         n = ColumnLength(i)
         write(funit,'(1x,a,1x,$)')ColumnTitle(i)(1:n)
      enddo
      write(funit,'(/,87(''-''))')
      FileExist = .true.
   endif
!
   n = ColumnLength(1)
   write(funit,'(a,1x,a,1x,$)')c(1:1),ColumnValue(1)(1:n)
   do i = 2, num_columns-1
      n = ColumnLength(i)
      write(funit,'(1x,a,1x,$)')ColumnValue(i)(1:n)
   enddo
   n = ColumnLength(num_columns)
   write(funit,'(1x,a)')ColumnValue(num_columns)(1:n)
   call FlushFile(funit)
!
   num_columns_old = num_columns
   num_columns = 0
   num_rows = num_rows + 1
!
   end subroutine writeRow1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeRow2(c,nc)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: c
   integer (kind=IntKind), intent(in) :: nc
!
   integer (kind=IntKind) :: i, n
!
   if (.not.FileExist) then
      write(funit,'(87(''=''))')
      n = ColumnLength(1)
      write(funit,'(3x,a,1x,$)')ColumnTitle(1)(1:n)
      do i = 2, num_columns
         n = ColumnLength(i)
         write(funit,'(1x,a,1x,$)')ColumnTitle(i)(1:n)
      enddo
      write(funit,'(/,87(''-''))')
      FileExist = .true.
   endif
!
   n = ColumnLength(1)
   write(funit,'(a,1x,a,1x,$)')c(1:nc),ColumnValue(1)(1:n)
   do i = 2, num_columns-1
      n = ColumnLength(i)
      write(funit,'(1x,a,1x,$)')ColumnValue(i)(1:n)
   enddo
   n = ColumnLength(num_columns)
   write(funit,'(1x,a)')ColumnValue(num_columns)(1:n)
   call FlushFile(funit)
!
   num_columns_old = num_columns
   num_columns = 0
   num_rows = num_rows + nc
!
   end subroutine writeRow2
!  ===================================================================
end module BookKeepingModule
