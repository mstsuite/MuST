module DefaultParamModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
!
public :: initDefaultParam,  &
          endDefaultParam,   &
          getDefaultValue
!
interface getDefaultValue
   module procedure getDefaultValue_s0, getDefaultValue_s1
   module procedure getDefaultValue_i0, getDefaultValue_i1
   module procedure getDefaultValue_r0, getDefaultValue_r1
   module procedure getDefaultValue_c0, getDefaultValue_c1
end interface
!
private
   integer (kind=IntKind), parameter :: MaxParams = 150
   integer (kind=intKind) :: NumParams
!
   character (len=50), allocatable :: StoredKeys(:)
   character (len=80), allocatable :: DefaultValues(:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initDefaultParam()
!  ===================================================================
   implicit none
!
   allocate(StoredKeys(MaxParams), DefaultValues(MaxParams))
!
!  -------------------------------------------------------------------
   call loadDefaultKeyValues(StoredKeys,DefaultValues,NumParams)
!  -------------------------------------------------------------------
!
   if (NumParams > MaxParams) then
      call ErrorHandler('initDefaultParam','NumParams > MaxParams',   &
                        NumParams,MaxParams)
   endif
!
   end subroutine initDefaultParam
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endDefaultParam()
!  ===================================================================
   implicit none
!
   deallocate(StoredKeys, DefaultValues)
!
   end subroutine endDefaultParam
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_s0(key,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   value = DefaultValues(k)
!
   end function getDefaultValue_s0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_s1(key,n,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getDefaultValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else if (n == 1) then
      value(1) = DefaultValues(k)
   else
      call initString(DefaultValues(k))
      m = getNumTokens()
      if (m < n) then
!        -------------------------------------------------------------
         call ErrorHandler('getDefaultValue','Number of string values < n',m,n)
!        -------------------------------------------------------------
      endif
      do i = 1, n
         call readToken(i,value(i))
      enddo
      call endString()
   endif
!
   end function getDefaultValue_s1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_i0(key,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   read(DefaultValues(k),*,iostat=status)value
   if (status < 0) then
      call ErrorHandler('getDefaultValue','Invalid DefaultValue',     &
                        DefaultValues(k))
   endif
!
   end function getDefaultValue_i0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_i1(key,n,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getDefaultValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else
      read(DefaultValues(k),*,iostat=status)value(1:n)
      if (status < 0) then
         call ErrorHandler('getDefaultValue','Invalid DefaultValue',  &
                           DefaultValues(k))
      endif
   endif
!
   end function getDefaultValue_i1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_r0(key,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   read(DefaultValues(k),*,iostat=status)value
   if (status < 0) then
      call ErrorHandler('getDefaultValue','Invalid DefaultValue',     &
                        DefaultValues(k))
   endif
!
   end function getDefaultValue_r0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_r1(key,n,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getDefaultValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else
      read(DefaultValues(k),*,iostat=status)value(1:n)
      if (status < 0) then
         call ErrorHandler('getDefaultValue','Invalid DefaultValue',  &
                           DefaultValues(k))
      endif
   endif
!
   end function getDefaultValue_r1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_c0(key,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   read(DefaultValues(k),*,iostat=status)value
   if (status < 0) then
      call ErrorHandler('getDefaultValue','Invalid DefaultValue',     &
                        DefaultValues(k))
   endif
!
   end function getDefaultValue_c0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDefaultValue_c1(key,n,value) result(s)
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
   LOOP_i: do i = 1, NumParams
      if (nocaseCompare(key_t,StoredKeys(i))) then
         k = i
         found = .true.
         exit LOOP_i
      endif
   enddo LOOP_i
   if (.not.found) then
      s = 1
!     ----------------------------------------------------------------
!     call WarningHandler('getDefaultValue','Invalid Key',key_t)
!     ----------------------------------------------------------------
      return
   endif
!
   if (n < 1) then
!     ----------------------------------------------------------------
      call ErrorHandler('getDefaultValue','Number of values < 1',n)
!     ----------------------------------------------------------------
   else
      read(DefaultValues(k),*,iostat=status)value(1:n)
      if (status < 0) then
         call ErrorHandler('getDefaultValue','Invalid DefaultValue',  &
                           DefaultValues(k))
      endif
   endif
!
   end function getDefaultValue_c1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine loadDefaultKeyValues(Keys,Values,n)
!  ===================================================================
   implicit none
!
   character (len=50), intent(out) :: Keys(:)
   character (len=80), intent(out) :: Values(:)
!
   integer (kind=IntKind), intent(out) :: n
!
   n = 0; Keys = ' '; Values = ' '
!
   include '../src/DefaultParameters.h'
!
   if (n > size(Keys)) then
      call ErrorHandler('loadDefaultKeyValues','n > size of Keys',n,size(Keys))
   endif
!
   end subroutine loadDefaultKeyValues
!  ===================================================================
end module DefaultParamModule
