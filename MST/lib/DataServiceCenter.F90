!  *******************************************************************
!  *  DataServiceCenterModule                                        *
!  *      Author: Yang Wang, Aurelian Rusanu, Malcolm Stocks         *
!  *                                                                 *
!  *      Version: 1.0         April 5th, 2005                       *
!  *                                                                 *
!  *******************************************************************
module DataServiceCenterModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
public :: initDataServiceCenter, &
          endDataServiceCenter,  &
          isDataStorageExisting, &
          createDataStorage,     &
          getDataStorage,        &
          getDataStorageType,    &
          getDataStorageSize,    &
          deleteDataStorage,     &
          setDataStorage2Value,  &
          setDataStorageLDA,     &
          getDataStorageLDA

public :: RealType, IntegerType, ComplexType, CharacterType
   integer (kind=IntKind), parameter :: IntegerType = 1
   integer (kind=IntKind), parameter :: RealType = 2
   integer (kind=IntKind), parameter :: ComplexType = 3
   integer (kind=IntKind), parameter :: CharacterType = 4
!
public :: RealMark, IntegerMark, ComplexMark, CharacterMark
   integer (kind=IntKind), parameter :: IntegerMark = 0
   real (kind=RealKind), parameter :: RealMark = ZERO
   complex (kind=CmplxKind), parameter :: ComplexMark = CZERO
   character (len=1), parameter :: CharacterMark = '0'
!
   interface createDataStorage
      module procedure createDataStorage_v, createDataStorage_vi
      module procedure createDataStorage_a, createDataStorage_ai
   end interface
!
   interface setDataStorage2Value
      module procedure setds2v_i, setds2v_r, setds2v_c, setds2v_s
      module procedure setds2v_ii, setds2v_ri, setds2v_ci, setds2v_si
   end interface
!
   interface setDataStorageLDA
      module procedure setlda, setlda_i
   end interface
!
   interface getDataStorageLDA
      module procedure getlda, getlda_i
   end interface
!
   interface getDataStorage
      module procedure getds_i0, getds_r0, getds_c0, getds_s0
      module procedure getds_i0i, getds_r0i, getds_c0i, getds_s0i
      module procedure getds_i1, getds_r1, getds_c1, getds_s1
      module procedure getds_i1i, getds_r1i, getds_c1i, getds_s1i
      module procedure getds_i2, getds_r2, getds_c2, getds_s2
      module procedure getds_i2i, getds_r2i, getds_c2i, getds_s2i
      module procedure getds_i3, getds_r3, getds_c3, getds_s3
      module procedure getds_i3i, getds_r3i, getds_c3i, getds_s3i
      module procedure getds_i4, getds_r4, getds_c4, getds_s4
      module procedure getds_i4i, getds_r4i, getds_c4i, getds_s4i
   end interface
!
private
   integer (kind=IntKind) :: NumStorages
   integer (kind=IntKind) :: TotalIntegerSize
   integer (kind=IntKind) :: TotalRealSize
   integer (kind=IntKind) :: TotalComplexSize
   integer (kind=IntKind) :: TotalCharacterSize
!
   type DataStorageStruct
      character (len=50) :: Key
      logical :: SingleValue
      integer (kind=IntKind) :: NumIndexes
      integer (kind=IntKind) :: Index
      integer (kind=IntKind) :: DataSize
      integer (kind=IntKind) :: DataType
      integer (kind=IntKind) :: LDASize                 ! leading dimension
      integer (kind=IntKind), pointer :: value_integer
      integer (kind=IntKind), pointer :: array_integer(:)
      real (kind=RealKind), pointer :: value_real
      real (kind=RealKind), pointer :: array_real(:)
      complex (kind=CmplxKind), pointer :: value_complex
      complex (kind=CmplxKind), pointer :: array_complex(:)
      character (len=1), pointer :: value_character
      character (len=1), pointer :: array_character(:)
      type (DataStorageStruct), pointer :: next
   end type DataStorageStruct

   type (DataStorageStruct), pointer :: head
   type (DataStorageStruct), pointer :: tail
!
   character (len=1), pointer :: array_character2(:,:)
   character (len=1), pointer :: array_character3(:,:,:)
!
   integer (kind=IntKind), pointer :: array_integer2(:,:)
   integer (kind=IntKind), pointer :: array_integer3(:,:,:)
!
   real (kind=RealKind), pointer :: array_real2(:,:)
   real (kind=RealKind), pointer :: array_real3(:,:,:)
!
   complex (kind=RealKind), pointer :: array_complex2(:,:)
   complex (kind=RealKind), pointer :: array_complex3(:,:,:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initDataServiceCenter()
!  ===================================================================
   implicit none
!
   NumStorages = 0
   TotalIntegerSize = 0
   TotalRealSize = 0
   TotalComplexSize = 0
   TotalCharacterSize = 0
!
   nullify(head)
   nullify(tail)
!
   end subroutine initDataServiceCenter
!  ===================================================================
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endDataServiceCenter()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i
!
   type (DataStorageStruct), pointer :: temp
!
   temp => head
   do i=1,NumStorages
      if (temp%DataSize < 1) then
         call ErrorHandler("endDataServiceCenter","DataSize < 1",temp%DataSize)
      endif
!
      if (temp%DataType == IntegerType) then
         TotalIntegerSize = TotalIntegerSize - temp%DataSize
         if (temp%DataSize == 1 .and. associated(temp%value_integer)) then
            deallocate( temp%value_integer )
         else
            deallocate( temp%array_integer )
         endif
      else if (temp%DataType == RealType) then
         TotalRealSize = TotalRealSize - temp%DataSize
         if (temp%DataSize == 1 .and. associated(temp%value_real)) then
            deallocate( temp%value_real )
         else
            deallocate( temp%array_real )
         endif
      else if (temp%DataType == ComplexType) then
         TotalComplexSize = TotalComplexSize - temp%DataSize
         if (temp%DataSize == 1 .and. associated(temp%value_complex)) then
            deallocate( temp%value_complex )
         else
            deallocate( temp%array_complex )
         endif
      else if (temp%DataType == CharacterType) then
         TotalCharacterSize = TotalCharacterSize - temp%DataSize
         if (temp%DataSize == 1 .and. associated(temp%value_character)) then
            deallocate( temp%value_character )
         else
            deallocate( temp%array_character )
         endif
      endif
      temp => temp%next
      nullify( head%next )
      deallocate( head )
      head => temp
   enddo
!
   if (TotalIntegerSize /= 0) then
      call ErrorHandler("endDataServiceCenter","TotalIntegerSize <> 0",   &
                        TotalIntegerSize )
   else if (TotalRealSize /= 0) then
      call ErrorHandler("endDataServiceCenter","TotalRealSize <> 0",      &
                        TotalRealSize )
   else if (TotalComplexSize /= 0) then
      call ErrorHandler("endDataServiceCenter","TotalComplexSize <> 0",   &
                        TotalComplexSize )
   else if (TotalCharacterSize /= 0) then
      call ErrorHandler("endDataServiceCenter","TotalCharacterSize <> 0", &
                        TotalCharacterSize )
   endif
!
   NumStorages = 0
!
   nullify(head)
   nullify(tail)
   nullify(temp)
!
   end subroutine endDataServiceCenter
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createDataStorage_v(key,data_type)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: data_type
!
   if (NumStorages == 0) then
       allocate (head)
       tail => head
   else if (NumStorages == 1) then
       allocate (head%next)
       tail => head%next
   else
       allocate (tail%next)
       tail => tail%next
   endif
   nullify(tail%next)
!
   tail%Key = key
   tail%Index = 0
   tail%NumIndexes = 1
   tail%DataSize = 1
   tail%LDASize = -1
   tail%SingleValue = .false.
!
   if (data_type == IntegerType) then
      TotalIntegerSize = TotalIntegerSize + 1
      allocate( tail%value_integer )
      nullify( tail%array_integer )
   else if (data_type == RealType) then
      TotalRealSize = TotalRealSize + 1
      allocate( tail%value_real )
      nullify( tail%array_real )
   else if (data_type == ComplexType) then
      TotalComplexSize = TotalComplexSize + 1
      allocate( tail%value_complex )
      nullify( tail%array_complex )
   else if (data_type == CharacterType) then
      TotalCharacterSize = TotalCharacterSize + 1
      allocate( tail%value_character )
      nullify( tail%array_character )
   else
      call ErrorHandler("createDataStorage", &
                        trim(key)//": Unknown data type",data_type)
   endif
!
   tail%DataType = data_type
!
   NumStorages = NumStorages + 1
!
   end subroutine createDataStorage_v
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createDataStorage_vi(ni,key,data_type)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: ni
   integer (kind=IntKind), intent(in) :: data_type
   integer (kind=IntKind) :: i
!
   if (ni < 1) then
      call ErrorHandler("createDataStorage", &
                        trim(key)//": Invalid number of indexes",ni)
   endif
!
   do i=1,ni
      if (NumStorages == 0) then
          allocate (head)
          tail => head
      else if (NumStorages == 1) then
          allocate (head%next)
          tail => head%next
      else
          allocate (tail%next)
          tail => tail%next
      endif
      nullify(tail%next)
!
      tail%Key = key
      tail%NumIndexes = ni
      tail%Index = i
      tail%DataSize = 1
      tail%LDASize = -1
      tail%SingleValue = .false.
!
      if (data_type == IntegerType) then
         TotalIntegerSize = TotalIntegerSize + 1
         allocate( tail%value_integer )
         nullify( tail%array_integer )
      else if (data_type == RealType) then
         TotalRealSize = TotalRealSize + 1
         allocate( tail%value_real )
         nullify( tail%array_real )
      else if (data_type == ComplexType) then
         TotalComplexSize = TotalComplexSize + 1
         allocate( tail%value_complex )
         nullify( tail%array_complex )
      else if (data_type == CharacterType) then
         TotalCharacterSize = TotalCharacterSize + 1
         allocate( tail%value_character )
         nullify( tail%array_character )
      else
         call ErrorHandler("createDataStorage", &
                        trim(key)//": Unknown data type",data_type)
      endif
      tail%DataType = data_type
      NumStorages = NumStorages + 1
   enddo
!
   end subroutine createDataStorage_vi
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createDataStorage_a(key,data_size,data_type)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: data_size
   integer (kind=IntKind), intent(in) :: data_type
!
   if (data_size < 1) then
      call ErrorHandler("createDataStorage", &
                        trim(key)//": Invalid data size ",data_size)
   endif
!
   if (NumStorages == 0) then
       allocate (head)
       tail => head
   else if (NumStorages == 1) then
       allocate (head%next)
       tail => head%next
   else
       allocate (tail%next)
       tail => tail%next
   endif
   nullify(tail%next)
!
   tail%Key = key
   tail%Index = 0
   tail%NumIndexes = 1
   tail%DataSize = data_size
   tail%LDASize = -1
   tail%SingleValue = .false.
!
   if (data_type == IntegerType) then
      TotalIntegerSize = TotalIntegerSize + data_size
      allocate( tail%array_integer(data_size) )
      nullify( tail%value_integer )
   else if (data_type == RealType) then
      TotalRealSize = TotalRealSize + data_size
      allocate( tail%array_real(data_size) )
      nullify( tail%value_real )
   else if (data_type == ComplexType) then
      TotalComplexSize = TotalComplexSize + data_size
      allocate( tail%array_complex(data_size) )
      nullify( tail%value_complex )
   else if (data_type == CharacterType) then
      TotalCharacterSize = TotalCharacterSize + data_size
      allocate( tail%array_character(data_size) )
      nullify( tail%value_character )
   else
      call ErrorHandler("createDataStorage", &
                        trim(key)//": Unknown data type",data_type)
   endif
!
   tail%DataType = data_type
!
   NumStorages = NumStorages + 1
!
   end subroutine createDataStorage_a
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createDataStorage_ai(ni,key,data_size,data_type)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: ni
!  integer (kind=IntKind), intent(in) :: data_size(*)
   integer (kind=IntKind), intent(in) :: data_size(:)
   integer (kind=IntKind), intent(in) :: data_type
   integer (kind=IntKind) :: i
!
   if (ni < 1) then
      call ErrorHandler("createDataStorage", &
                        trim(key)//": Invalid number of indexes",ni)
   endif
!
   do i=1,ni
      if (NumStorages == 0) then
          allocate (head)
          tail => head
      else if (NumStorages == 1) then
          allocate (head%next)
          tail => head%next
      else
          allocate (tail%next)
          tail => tail%next
      endif
      nullify(tail%next)
!
      tail%Key = key
      tail%NumIndexes = ni
      tail%Index = i
      tail%DataSize = data_size(i)
      tail%LDASize = -1
      tail%SingleValue = .false.
!
      if (data_type == IntegerType) then
         TotalIntegerSize = TotalIntegerSize + data_size(i)
         allocate( tail%array_integer(data_size(i)) )
         nullify( tail%value_integer )
      else if (data_type == RealType) then
         TotalRealSize = TotalRealSize + data_size(i)
         allocate( tail%array_real(data_size(i)) )
         nullify( tail%value_real )
      else if (data_type == ComplexType) then
         TotalComplexSize = TotalComplexSize + data_size(i)
         allocate( tail%array_complex(data_size(i)) )
         nullify( tail%value_complex )
      else if (data_type == CharacterType) then
         TotalCharacterSize = TotalCharacterSize + data_size(i)
         allocate( tail%array_character(data_size(i)) )
         nullify( tail%value_character )
      else
         call ErrorHandler("createDataStorage", &
                        trim(key)//": Unknown data type",data_type)
      endif
      tail%DataType = data_type
      NumStorages = NumStorages + 1
   enddo
!
   end subroutine createDataStorage_ai
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i0(key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_integer)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_integer
!
   end function getds_i0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i0i(index,key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_integer)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_integer
!
   end function getds_i0i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r0(key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_real)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_real
!
   end function getds_r0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r0i(index,key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_real)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_real
!
   end function getds_r0i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c0(key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_complex)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_complex
!
   end function getds_c0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c0i(index,key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_complex)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_complex
!
   end function getds_c0i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s0(key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p
!
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_character)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_character
!
   end function getds_s0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s0i(index,key,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s /= 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is not 1",s)
   else if (.not.associated(temp%value_character)) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data storage is not created")
   endif
!
   p => temp%value_character
!
   end function getds_s0i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i1(key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p(:)
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_integer(1:n1)
!
   end function getds_i1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i1i(index,key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p(:)
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_integer(1:n1)
!
   end function getds_i1i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r1(key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_real(1:n1)
!
   end function getds_r1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r1i(index,key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_real(1:n1)
!
   end function getds_r1i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c1(key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_complex(1:n1)
!
   end function getds_c1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c1i(index,key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_complex(1:n1)
!
   end function getds_c1i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s1(key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p(:)
!
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_character(1:n1)
!
   end function getds_s1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s1i(index,key,n1,markup) result(p)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p(:)
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size < n",n1)
   else if (n1 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n",s,n1)
   endif
!
   p => temp%array_character(1:n1)
!
   end function getds_s1i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i2(key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1, n2
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p2(:,:), p1(:)
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_integer(1:n*n2)
   p2 => aliasArray2_i(p1,n,n2)
!   call relayIntegerArray2(p1,n,n2)
!   p2 => array_integer2(1:n,1:n2)
!
   end function getds_i2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i2i(index,key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p2(:,:), p1(:)
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_integer(1:n*n2)
   p2 => aliasArray2_i(p1,n,n2)
!   call relayIntegerArray2(p1,n,n2)
!   p2 => array_integer2(1:n,1:n2)
!
   end function getds_i2i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r2(key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=IntKind) :: s, n
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p2(:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_real(1:n*n2)
   p2 => aliasArray2_r(p1,n,n2)
!   call relayRealArray2(p1,n,n2)
!   p2 => array_real2(1:n,1:n2)
!
   end function getds_r2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r2i(index,key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=IntKind) :: s, n
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p2(:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_real(1:n*n2)
   p2 => aliasArray2_r(p1,n,n2)
!   call relayRealArray2(p1,n,n2)
!   p2 => array_real2(1:n,1:n2)
!
   end function getds_r2i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c2(key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=IntKind) :: s, n
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p2(:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_complex(1:n*n2)
   p2 => aliasArray2_c(p1,n,n2)
!   call relayComplexArray2(p1,n,n2)
!   p2 => array_complex2(1:n,1:n2)
!
   end function getds_c2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c2i(index,key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s, n
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p2(:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_complex(1:n*n2)
   p2 => aliasArray2_c(p1,n,n2)
!   call relayComplexArray2(p1,n,n2)
!   p2 => array_complex2(1:n,1:n2)
!
   end function getds_c2i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s2(key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p2(:,:), p1(:)
!
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_character(1:n*n2)
   p2 => aliasArray2_s(p1,n,n2)
!   call relayCharacterArray2(p1,n,n2)
!   p2 => array_character2(1:n,1:n2)
!
   end function getds_s2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s2i(index,key,n1,n2,markup) result(p2)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p2(:,:), p1(:)
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2)
   else if (n1*n2 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2",s,n1*n2)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_character(1:n*n2)
   p2 => aliasArray2_s(p1,n,n2)
!   call relayCharacterArray2(p1,n,n2)
!   p2 => array_character2(1:n,1:n2)
!
   end function getds_s2i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i3(key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1, n2,n3
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p3(:,:,:), p1(:)
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_integer(1:n*n2*n3)
   p3 => aliasArray3_i(p1,n,n2,n3)
!   call relayIntegerArray3(p1,n,n2,n3)
!   p3 => array_integer3(1:n,1:n2,1:n3)
!
   end function getds_i3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i3i(index,key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2,n3
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p3(:,:,:), p1(:)
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage","Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_integer(1:n*n2*n3)
   p3 => aliasArray3_i(p1,n,n2,n3)
!   call relayIntegerArray3(p1,n,n2,n3)
!   p3 => array_integer3(1:n,1:n2,1:n3)
!
   end function getds_i3i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r3(key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3
   integer (kind=IntKind) :: s, n
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p3(:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_real(1:n*n2*n3)
   p3 => aliasArray3_r(p1,n,n2,n3)
!   call relayRealArray3(p1,n,n2,n3)
!   p3 => array_real3(1:n,1:n2,1:n3)
!
   end function getds_r3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r3i(index,key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2,n3
   integer (kind=IntKind) :: s, n
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p3(:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_real(1:n*n2*n3)
   p3 => aliasArray3_r(p1,n,n2,n3)
!   call relayRealArray3(p1,n,n2,n3)
!   p3 => array_real3(1:n,1:n2,1:n3)
!
   end function getds_r3i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c3(key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3
   integer (kind=IntKind) :: s, n
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p3(:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_complex(1:n*n2*n3)
   p3 => aliasArray3_c(p1,n,n2,n3)
!   call relayComplexArray3(p1,n,n2,n3)
!   p3 => array_complex3(1:n,1:n2,1:n3)
!
   end function getds_c3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c3i(index,key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s, n
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p3(:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_complex(1:n*n2*n3)
   p3 => aliasArray3_c(p1,n,n2,n3)
!   call relayComplexArray3(p1,n,n2,n3)
!   p3 => array_complex3(1:n,1:n2,1:n3)
!
   end function getds_c3i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s3(key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p3(:,:,:), p1(:)
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_character(1:n*n2*n3)
   p3 => aliasArray3_s(p1,n,n2,n3)
!   call relayCharacterArray3(p1,n,n2,n3)
!   p3 => array_character3(1:n,1:n2,1:n3)
!
   end function getds_s3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s3i(index,key,n1,n2,n3,markup) result(p3)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p3(:,:,:), p1(:)
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2,n3
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3)
   else if (n1*n2*n3 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3",s,n1*n2*n3)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_character(1:n*n2*n3)
   p3 => aliasArray3_s(p1,n,n2,n3)
!   call relayCharacterArray3(p1,n,n2,n3)
!   p3 => array_character3(1:n,1:n2,1:n3)
!
   end function getds_s3i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i4(key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1, n2, n3, n4
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p4(:,:,:,:), p1(:)
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_integer(1:n*n2*n3*n4)
   p4 => aliasArray4_i(p1,n,n2,n3,n4)
!
   end function getds_i4
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_i4i(index,key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer (kind=IntKind), intent(in) :: markup
   integer (kind=IntKind), pointer :: p4(:,:,:,:), p1(:)
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage","Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_integer(1:n*n2*n3*n4)
   p4 => aliasArray4_i(p1,n,n2,n3,n4)
!
   end function getds_i4i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r4(key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer (kind=IntKind) :: s, n
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p4(:,:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_real(1:n*n2*n3*n4)
   p4 => aliasArray4_r(p1,n,n2,n3,n4)
!
   end function getds_r4
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_r4i(index,key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer (kind=IntKind) :: s, n
!
   real (kind=RealKind), intent(in) :: markup
   real (kind=RealKind), pointer :: p4(:,:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_real(1:n*n2*n3*n4)
   p4 => aliasArray4_r(p1,n,n2,n3,n4)
!
   end function getds_r4i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c4(key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer (kind=IntKind) :: s, n
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p4(:,:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_complex(1:n*n2*n3*n4)
   p4 => aliasArray4_c(p1,n,n2,n3,n4)
!
   end function getds_c4
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_c4i(index,key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s, n
!
   complex (kind=CmplxKind), intent(in) :: markup
   complex (kind=CmplxKind), pointer :: p4(:,:,:,:), p1(:)
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_complex(1:n*n2*n3*n4)
   p4 => aliasArray4_c(p1,n,n2,n3,n4)
!
   end function getds_c4i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s4(key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p4(:,:,:,:), p1(:)
!
   integer (kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_character(1:n*n2*n3*n4)
   p4 => aliasArray4_s(p1,n,n2,n3,n4)
!
   end function getds_s4
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getds_s4i(index,key,n1,n2,n3,n4,markup) result(p4)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=1), intent(in) :: markup
!
   character (len=1), pointer :: p4(:,:,:,:), p1(:)
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer (kind=IntKind) :: s, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size is 0",s)
   else if (n1 < 1 .or. n2 < 1 .or. n3 < 1 .or. n4 < 1) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Invalid dimension size",n1,n2,n3,n4)
   else if (n1*n2*n3*n4 > s) then
      call ErrorHandler("getDataStorage", &
                        trim(key)//": Data size < n1*n2*n3*n4",s,n1*n2*n3*n4)
   endif
!
   n = n1
   if ( temp%LDASize>=0 ) then
      n = temp%LDASize
   endif
!
   p1 => temp%array_character(1:n*n2*n3*n4)
   p4 => aliasArray4_s(p1,n,n2,n3,n4)
!
   end function getds_s4i
!  ===================================================================
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDataStorageType(key) result(t)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: t
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
!
   t = temp%DataType
!
   end function getDataStorageType
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDataStorageSize(key) result(s)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
!
   s = temp%DataSize
!
   end function getDataStorageSize
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteDataStorage(key)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   type (DataStorageStruct), pointer :: temp1
   type (DataStorageStruct), pointer :: temp2
!
   temp2 => searchStorage(key)
!
   if (temp2%DataSize < 1) then
      call ErrorHandler("deleteDataStorage",                           &
                        trim(key)//": Data Size < 1",temp2%DataSize)
   endif
!
   if (temp2%DataType == IntegerType) then
      TotalIntegerSize = TotalIntegerSize - temp2%DataSize
      if (temp2%DataSize == 1 .and. associated(temp2%value_integer)) then
         deallocate( temp2%value_integer )
      else
         deallocate( temp2%array_integer )
      endif
   else if (temp2%DataType == RealType) then
      TotalRealSize = TotalRealSize - temp2%DataSize
      if (temp2%DataSize == 1 .and. associated(temp2%value_real)) then
         deallocate( temp2%value_real )
      else
         deallocate( temp2%array_real )
      endif
   else if (temp2%DataType == ComplexType) then
      TotalComplexSize = TotalComplexSize - temp2%DataSize
      if (temp2%DataSize == 1 .and. associated(temp2%value_complex)) then
         deallocate( temp2%value_complex )
      else
         deallocate( temp2%array_complex )
      endif
   else if (temp2%DataType == CharacterType) then
      TotalCharacterSize = TotalCharacterSize - temp2%DataSize
      if (temp2%DataSize == 1 .and. associated(temp2%value_character)) then
         deallocate( temp2%value_character )
      else
         deallocate( temp2%array_character )
      endif
   endif
!
   if (NumStorages == 1) then
      nullify( head )
      nullify( tail )
   else
      temp1 => searchPrevStorage(key)
      if ( associated(temp1) )then
         temp1%next => temp2%next
      else
         head => temp2%next
      endif
      if (.not.associated(temp2%next) )then
         tail => temp1
      endif
   endif
   nullify( temp2%next )
   deallocate( temp2 )
!
   NumStorages = NumStorages - 1
!
   end subroutine deleteDataStorage
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isDataStorageExisting(key) result(t)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   logical :: t
!
   integer (kind=IntKind) :: i
!
   type (DataStorageStruct), pointer :: temp
!
   if (NumStorages > 0) then
      temp => head
      do i=1,NumStorages
         if (temp%Key == key) then
            t = .true.
            return
         else
            temp => temp%next
         endif
      enddo
   endif
   t = .false.
!
   end function isDataStorageExisting
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchStorage(key) result(temp)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: i
!
   type (DataStorageStruct), pointer :: temp
!
   if (NumStorages < 1) then
      call ErrorHandler("searchStorage",                              &
                         trim(key)//": Empty data storage")
   endif
!
   temp => head
   do i=1,NumStorages
      if (temp%Key == key) then
         return
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchStorage",                                 &
                     "Data storage is not found",key,force_to_print=.true.)
!
   end function searchStorage
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchStorage_i(index,key) result(temp)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: i, j, n
!
   type (DataStorageStruct), pointer :: temp
!
   if (NumStorages < 1) then
      call ErrorHandler("searchStorage_i",                            &
                        trim(key)//": Empty data storage")
   else if (index < 0) then
      call ErrorHandler("searchStorage_i",                            &
                        trim(key)//": Invalid index",index)
   endif
!
   temp => head
   do i=1,NumStorages
      if (temp%Key == key) then
         n = temp%NumIndexes
         if (index == 0 .and. n == 0) then
            return
         else if (index == 0 .and. n > 0) then
            call ErrorHandler("searchStorage_i",                     &
                              trim(key)//": Invalid index",index)
         else if (index > n) then
            call ErrorHandler("searchStorage_i",                     &
                              trim(key)//": Invalid index",index)
         else
            do j=1,n
               if (temp%Index == index) then
                  return
               else
                  temp => temp%next
               endif
            enddo
         endif
         call ErrorHandler("searchStorage_i",                         &
               trim(key)//": Data with the index is not found",index)
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchStorage_i",                                &
                    trim(key)//": Data with the key is not found",key)
!
   end function searchStorage_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchPrevStorage(key) result(temp)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: i
!
   type (DataStorageStruct), pointer :: temp
!
   if (NumStorages < 1) then
      call ErrorHandler("searchPrevStorage", &
                         trim(key)//": Empty data storage")
   else if (NumStorages == 1) then
      nullify( temp )
      return
   endif
!
   temp => head
   do i=1,NumStorages
      if (temp%next%Key == key) then
         return
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchPrevStorage", &
                    trim(key)//": Data with the key is not found",key)
!
   end function searchPrevStorage
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_i(key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: v
   integer (kind=IntKind) :: s, i, ni, ni0
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s   = temp%DataSize
   ni  = temp%NumIndexes
   ni0 = temp%Index
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("setDataStorage",                              &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage",                              &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_integer)) then
       temp%value_integer = v
       if ( ni0==0 ) return
       do i = ni0+1,ni
          temp => temp%next
          temp%value_integer = v
       enddo
   else
       temp%array_integer(1:s) = v
       if ( ni0==0 ) return
       do i = ni0+1,ni
          temp => temp%next
          s = temp%DataSize
          temp%array_integer(1:s) = v
       enddo
   endif
!
   end subroutine setds2v_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_ii(index,key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: v
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= IntegerType) then
      call ErrorHandler("setDataStorage",                              &
                        trim(key)//": Data is not of integer type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage",                              &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_integer)) then
       temp%value_integer = v
   else
       temp%array_integer(1:s) = v
   endif
!
   end subroutine setds2v_ii
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_r(key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: s, i, ni, ni0
!
   real (kind=RealKind), intent(in) :: v
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s   = temp%DataSize
   ni  = temp%NumIndexes
   ni0 = temp%Index
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("setDataStorage",                             &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage",                             &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_real)) then
       temp%value_real = v
       if ( ni0==0 ) return
       do i = ni0+1,ni
          temp => temp%next
          s = temp%DataSize
          temp%value_real = v
       enddo
   else
       temp%array_real(1:s) = v
       if ( ni0==0 ) return
       do i = ni0+1,ni
          temp => temp%next
          s = temp%DataSize
          temp%array_real(1:s) = v
       enddo
   endif
!
   end subroutine setds2v_r
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_ri(index,key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: v
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_real)) then
       temp%value_real = v
   else
       temp%array_real(1:s) = v
   endif
!
   end subroutine setds2v_ri
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_c(key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: s, i, ni, ni0
!
   complex (kind=CmplxKind), intent(in) :: v
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s   = temp%DataSize
   ni  = temp%NumIndexes
   ni0 = temp%Index
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_complex)) then
       temp%value_complex = v
       if ( ni0==0 ) return
       do i = ni0+1,ni
          temp => temp%next
          temp%value_complex = v
       enddo
   else
       temp%array_complex(1:s) = v
       if ( ni0==0 ) return
       do i = ni0+1,ni
          temp => temp%next
          s = temp%DataSize
          temp%array_complex(1:s) = v
       enddo
   endif
!
   end subroutine setds2v_c
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_ci(index,key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: v
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= ComplexType) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_complex)) then
       temp%value_complex = v
   else
       temp%array_complex(1:s) = v
   endif
!
   end subroutine setds2v_ci
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_s(key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=*), intent(in) :: v
!
   integer (kind=IntKind) :: s, i, n, j, ni, ni0
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s   = temp%DataSize
   ni  = temp%NumIndexes
   ni0 = temp%Index
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_character)) then
       temp%value_character = v(1:1)
       if ( ni0==0 ) return
       do j = ni0+1,ni
          temp => temp%next
          temp%value_character = v(1:1)
       enddo
   else
       n = min(len(v),s)
       do i = 1, n
          temp%array_character(i) = v(i:i)
       enddo
       do i = n+1, s
          temp%array_character(i) = ' '
       enddo
       if ( ni0==0 ) return
       do j = ni0+1,ni
          temp => temp%next
          s = temp%DataSize
          n = min(len(v),s)
          do i = 1, n
             temp%array_character(i) = v(i:i)
          enddo
          do i = n+1, s
             temp%array_character(i) = ' '
          enddo
       enddo
   endif
!
   end subroutine setds2v_s
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setds2v_si(index,key,v)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
   character (len=*), intent(in) :: v
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: s, i, n
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CharacterType) then
      call ErrorHandler("setDataStorage", & 
                        trim(key)//": Data is not of character type")
   else if (s < 1) then
      call ErrorHandler("setDataStorage", &
                        trim(key)//": Data size is 0",s)
   endif
!
   if (associated(temp%value_character)) then
       temp%value_character = v(1:1)
   else
       n = min(len(v),s)
       do i = 1, n
          temp%array_character(i) = v(i:i)
       enddo
       do i = n+1, s
          temp%array_character(i) = ' '
       enddo
   endif
!
   end subroutine setds2v_si
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setlda(key,lda)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: lda
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (s < 1) then
      call ErrorHandler("setDataStorageLDA", &
                        trim(key)//": Data size is 0",s)
   endif
!
   temp%LDASize = lda
!
   end subroutine setlda
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setlda_i(index,key,lda)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind), intent(in) :: lda
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (s < 1) then
      call ErrorHandler("setDataStorageLDA", &
                        trim(key)//": Data size is 0",s)
   endif
!
   temp%LDASize = lda
!
   end subroutine setlda_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getlda(key)                                    result(lda)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind) :: lda
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage(key)
   s = temp%DataSize
!
   if (s < 1) then
      call ErrorHandler("getDataStorageLDA",trim(key)//": Data size is 0",s)
   endif
!
   lda = temp%LDASize
!
   end function getlda
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getlda_i(index,key)                            result(lda)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: lda
   integer (kind=IntKind) :: s
!
   type (DataStorageStruct), pointer :: temp
!
   temp => searchStorage_i(index,key)
   s = temp%DataSize
!
   if (s < 1) then
      call ErrorHandler("getDataStorageLDA",trim(key)//": Data size is 0",s)
   endif
!
   lda = temp%LDASize
!
   end function getlda_i
!  ===================================================================
end module DataServiceCenterModule
