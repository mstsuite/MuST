module CmdLineOptionModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
!
public :: initCmdLineOption,  &
          endCmdLineOption,   &
          getCmdLineOption,   &
          getCmdLineOptionValue
!
interface getCmdLineOptionValue
   module procedure getCmdLineOptionValue_s0, getCmdLineOptionValue_s1
   module procedure getCmdLineOptionValue_i0, getCmdLineOptionValue_i1
   module procedure getCmdLineOptionValue_r0, getCmdLineOptionValue_r1
   module procedure getCmdLineOptionValue_c0, getCmdLineOptionValue_c1
end interface
!
private
   integer (kind=IntKind), parameter :: MaxOptions = 20
!
   integer (kind=intKind) :: NumOptions = 0
!
   character (len=20) :: Options(3,MaxOptions)
   character (len=50) :: OptionKeys(MaxOptions)
   character (len=50) :: OptionValues(MaxOptions)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCmdLineOption(num_args)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: num_args
!
   NumOptions = 0
!
!  -------------------------------------------------------------------
   call loadCmdLineOptions(num_args)
!  -------------------------------------------------------------------
!
   end subroutine initCmdLineOption
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCmdLineOption()
!  ===================================================================
   implicit none
!
   NumOptions = 0
!
   end subroutine endCmdLineOption
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOption(key) result(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=50) :: key_t
!
   integer (kind=IntKind) :: s, i, k
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         if (OptionValues(i) == '-1') then
            found = .false.
         else
            found = .true.
         endif
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   end function getCmdLineOption
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_s0(key,value) result(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=*), intent(inout) :: value
   character (len=50) :: key_t
!
   integer (kind=IntKind) :: s, i, k
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   value = OptionValues(k)
!
   end function getCmdLineOptionValue_s0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_s1(key,n,value) result(s)
!  ===================================================================
   use StringModule, only : initString, endString, getNumTokens, readToken
   implicit none
!
   integer (kind=intKind), intent(in) :: n
!
   character (len=*), intent(in) :: key
   character (len=*), intent(inout) :: value(n)
   character (len=50) :: key_t
!
   integer (kind=IntKind) :: s, i, k, m
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCmdLineOptionValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else if (n == 1) then
      value(1) = OptionValues(k)
   else
      call initString(OptionValues(k))
      m = getNumTokens()
      if (m < n) then
!        -------------------------------------------------------------
         call ErrorHandler('getCmdLineOptionValue','Number of string values < n',m,n)
!        -------------------------------------------------------------
      endif
      do i = 1, n
         call readToken(i,value(i))
      enddo
      call endString()
   endif
!
   end function getCmdLineOptionValue_s1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_i0(key,value) result(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=50) :: key_t
!
   integer (kind=IntKind), intent(inout) :: value
   integer (kind=IntKind) :: s, i, k, status
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   read(OptionValues(k),*,iostat=status)value
   if (status < 0) then
      call ErrorHandler('getCmdLineOptionValue','Invalid OptionValue',     &
                        OptionValues(k))
   endif
!
   end function getCmdLineOptionValue_i0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_i1(key,n,value) result(s)
!  ===================================================================
   implicit none
!
   integer (kind=intKind), intent(in) :: n
   integer (kind=IntKind), intent(inout) :: value(n)
!
   character (len=*), intent(in) :: key
   character (len=50) :: key_t
!
   integer (kind=IntKind) :: s, i, k, status
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCmdLineOptionValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else
      read(OptionValues(k),*,iostat=status)value(1:n)
      if (status < 0) then
         call ErrorHandler('getCmdLineOptionValue','Invalid OptionValue',  &
                           OptionValues(k))
      endif
   endif
!
   end function getCmdLineOptionValue_i1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_r0(key,value) result(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=50) :: key_t
!
   real (kind=RealKind), intent(inout) :: value
!
   integer (kind=IntKind) :: s, i, k, status
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   read(OptionValues(k),*,iostat=status)value
   if (status < 0) then
      call ErrorHandler('getCmdLineOptionValue','Invalid OptionValue',     &
                        OptionValues(k))
   endif
!
   end function getCmdLineOptionValue_r0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_r1(key,n,value) result(s)
!  ===================================================================
   implicit none
!
   integer (kind=intKind), intent(in) :: n
!
   real (kind=RealKind), intent(inout) :: value(n)
!
   character (len=*), intent(in) :: key
   character (len=50) :: key_t
!
   integer (kind=IntKind) :: s, i, k, status
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCmdLineOptionValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else
      read(OptionValues(k),*,iostat=status)value(1:n)
      if (status < 0) then
         call ErrorHandler('getCmdLineOptionValue','Invalid OptionValue',  &
                           OptionValues(k))
      endif
   endif
!
   end function getCmdLineOptionValue_r1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_c0(key,value) result(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=50) :: key_t
!
   complex (kind=CmplxKind), intent(inout) :: value
!
   integer (kind=IntKind) :: s, i, k, status
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   read(OptionValues(k),*,iostat=status)value
   if (status < 0) then
      call ErrorHandler('getCmdLineOptionValue','Invalid OptionValue',     &
                        OptionValues(k))
   endif
!
   end function getCmdLineOptionValue_c0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCmdLineOptionValue_c1(key,n,value) result(s)
!  ===================================================================
   implicit none
!
   integer (kind=intKind), intent(in) :: n
!
   complex (kind=CmplxKind), intent(inout) :: value(n)
!
   character (len=*), intent(in) :: key
   character (len=50) :: key_t
!
   integer (kind=IntKind) :: s, i, k, status
!
   logical :: found
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = 0
   found = .false.
!
   key_t = trim(adjustl(key))
   LOOP_i: do i = 1, NumOptions
      if (nocaseCompare(key_t,OptionKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getCmdLineOptionValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getCmdLineOptionValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else
      read(OptionValues(k),*,iostat=status)value(1:n)
      if (status < 0) then
         call ErrorHandler('getCmdLineOptionValue','Invalid OptionValue',  &
                           OptionValues(k))
      endif
   endif
!
   end function getCmdLineOptionValue_c1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine loadCmdLineOptions(num_args)
!  ===================================================================
   implicit none
!
   character (len=50) :: args
!
   integer (kind=IntKind), intent(out) :: num_args
   integer (kind=IntKind) :: i, j, n, m, ml, ieq
   integer (kind=IntKind) :: mj, mj_save
!
!  ===================================================================
!  Command line arguments list
!  -------------------------------------------------------------------
!  character (len=200) :: cmdline_arguments
!  character (len=256) :: msg
!  integer (kind=intKind) :: ios
!  integer (kind=intKind) :: NumOutputProcs, NumAtomsPerProc
!
!  namelist /cla/NumOutputProcs, NumAtomsPerProc
!  ===================================================================
!
   interface
      function nocaseCompare(s1,s2) result(t)
         implicit none
         logical :: t
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
      end function nocaseCompare
   end interface
!
   n = 0; Options = ' '; OptionKeys = ' '; OptionValues = ' '
!
   include '../src/CmdLineOptions.inc'
!
   if (n > MaxOptions) then
      call ErrorHandler('loadCmdLineOptions','n > MaxOptions',n,MaxOptions)
   endif
   NumOptions = n
!  print *,'NumOptions = ',NumOptions
!
!  ===================================================================
!  Start processing the command line arguments...........
!  -------------------------------------------------------------------
   n = command_argument_count()
!  -------------------------------------------------------------------
   num_args = 0
   mj_save = 0
   if (n > 0) then
!     cmdline_arguments = '&cla '
      do i = 1, n
!        -------------------------------------------------------------
         call get_command_argument(i,args)
!        -------------------------------------------------------------
!        print *,'Before processed: ',trim(args)
         ieq = index(args,'=')
         ml = len(trim(args))
         if (ieq == ml) then
            m = ml - 1
         else if (ieq > 0) then
            m = ieq - 1
         else
            m = ml
         endif
         mj = 0
         if (args(1:2) == '--') then
            LOOP_j2: do j = 1, NumOptions
               if (nocaseCompare(trim(args(3:m)),trim(Options(2,j)(3:))) .or. &
                   nocaseCompare(trim(args(3:m)),trim(Options(3,j)(3:)))) then
                  mj = j
                  exit LOOP_j2
               endif
            enddo LOOP_j2
         else if (args(1:1) == '-') then
            LOOP_j1: do j = 1, NumOptions
               if (nocaseCompare(trim(args(2:m)),trim(Options(1,j)(2:)))) then
                  mj = j
                  exit LOOP_j1
               endif
            enddo LOOP_j1
         endif
         if (mj > 0) then
            num_args = num_args + 1
            mj_save = mj
            if (ieq == 0 .or. ieq == ml) then
!              args = trim(OptName(j))//'='
               OptionValues(mj) = ' '
            else
!              args = trim(OptName(j))//args(ieq:ml)
               OptionValues(mj) = args(ieq+1:ml)
            endif
         else if (mj_save > 0) then
            if (ieq == 1) then
               args(1:1) = ' '
               args = adjustl(args)
            else if (ieq > 1) then
!              -------------------------------------------------------
               call ErrorHandler('loadCmdLineOptions','Error parsing the command line option',trim(args))
!              -------------------------------------------------------
            endif
            OptionValues(mj_save) = trim(OptionValues(mj_save))//args(1:ml)
         else
!           ----------------------------------------------------------
            call ErrorHandler('loadCmdLineOptions','Error parsing the command line option',trim(args))
!           ----------------------------------------------------------
         endif
!        print *,'After  processed: ',trim(args)
!        cmdline_arguments = trim(cmdline_arguments)//' '//trim(args)
      enddo
!     cmdline_arguments = trim(cmdline_arguments)//' /'
!     print *,'cmd-line args: ',trim(cmdline_arguments)
!     read(cmdline_arguments, nml=cla, iostat=ios, iomsg=msg)
!     if (ios /= 0) then
!        -------------------------------------------------------------
!        call ErrorHandler('initInput','Error parsing the command line',msg)
!        -------------------------------------------------------------
!     endif
!     print *,'NumOutputProcs  = ',NumOutputProcs
!     print *,'NumAtomsPerProc = ',NumAtomsPerProc
!  else
!     NumOutputProcs = 0
!     NumAtomsPerProc = 0
   endif
!  do i = 1, NumOptions
!     print *,trim(OptionKeys(i)),' = ',trim(OptionValues(i))
!  enddo
!  ===================================================================
!  End of processing the command line argumens
!  ===================================================================
!
   end subroutine loadCmdLineOptions
!  ===================================================================
end module CmdLineOptionModule
