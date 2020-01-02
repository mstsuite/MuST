module AnaliticContinuationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   use MathParamModule, only : ZERO, ONE, TEN2m6, HALF, THREE, CZERO
!
public :: initAnaliticCont,   &
          endAnaliticCont,    &
          isGridExisting,     &
          isDataExisting,     &
          createGridData,     &
!         getGridData,        &
!         delGridData,        &
!         getDiscreteData,    &
!         getFittedData,      &
!         delDiscreteData
          createDiscreteData
!
public :: RealType, CmplxType, RealMark, CmplxMark
   integer (kind=IntKind), parameter :: RealType = 0
   integer (kind=IntKind), parameter :: CmplxType = 1
   real (kind=RealKind), parameter :: RealMark = ZERO
   complex (kind=CmplxKind), parameter :: CmplxMark = CZERO
!
    interface createGridData
       module procedure createGridData_0, createGridData_ni
    end interface
!
    interface createDiscreteData
       module procedure createDiscreteData_0, createDiscreteData_ni
    end interface
!
    interface getDiscreteData
       module procedure getDiscreteData_r, getDiscreteData_ri, getDiscreteData_r2i
       module procedure getDiscreteData_c, getDiscreteData_ci, getDiscreteData_c2i
    end interface
!
!    interface getContData
!       module procedure getContData_a, getContData_v
!    end interface
!
private
!
   logical :: Initialized = .false.
!
   integer (kind=IntKind) :: NumGrids 
   integer (kind=IntKind) :: KeyGridGen 
!   
   integer (kind=IntKind) :: NumDiscreteData
   integer (kind=IntKind) :: KeyDiscreteGen
!
   type AnaliticStruct
      integer(kind=IntKind) :: Key
      integer(kind=IntKind) :: Index
      integer(kind=IntKind) :: NumIndexes
      integer(kind=IntKind) :: DataType
      integer(kind=IntKind) :: DataSize
      real (kind=RealKind), pointer :: array_r(:)
      complex(kind=CmplxKind), pointer :: array_c(:)
      type (AnaliticStruct), pointer :: next
   end type AnaliticStruct
!
   type(AnaliticStruct), pointer :: headGridData
   type(AnaliticStruct), pointer :: tailGridData
!
   type(AnaliticStruct), pointer :: headData
   type(AnaliticStruct), pointer :: tailData
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initAnaliticCont()
!  ===================================================================
   implicit none
!
   NumGrids = 0
   KeyGridGen = 0
!   
   NumDiscreteData = 0
   KeyDiscreteGen = 0
!
   nullify(headGridData)
   nullify(headData)
   nullify(tailGridData)
   nullify(tailData)
!
   Initialized = .true.   
!
   end subroutine initAnaliticCont
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endAnaliticCont()
!  ===================================================================
   implicit none
!
  integer (kind=IntKind) :: i
!
   type (AnaliticStruct), pointer :: temp
!
   temp => headGridData
   do i = 1,NumGrids
      if (temp%DataSize < 1) then
         call ErrorHandler("endAnaliticCont","DataSize < 1",temp%DataSize)
      endif
!
      if (temp%DataType == RealType) then
         deallocate( temp%array_r )         
      else if (temp%DataType == CmplxType) then
         deallocate( temp%array_c )
      endif
      temp => temp%next
      nullify( headGridData%next )
      deallocate( headGridData )
      headGridData => temp
   enddo
!
   NumGrids = 0
!
   nullify(headGridData)
   nullify(temp)
   nullify(tailGridData)
!
   temp => headData
   do i=1,NumDiscreteData
      if (temp%DataSize < 1) then
         call ErrorHandler("endAnaliticCont","DataSize < 1",temp%DataSize)
      endif
!
      if ( temp%DataType == RealType ) then
         deallocate( temp%array_r )
      else if ( temp%DataType == CmplxType ) then
         deallocate( temp%array_c )
      endif
      temp => temp%next
      nullify( headData%next )
      deallocate( headData )
      headData => temp
   enddo
!
   NumDiscreteData = 0
!
   nullify(headData)
   nullify(temp)
   nullify(tailData)
!
   Initialized = .false.   
!
   end subroutine endAnaliticCont
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createGridData_0( key, data_size, data_type )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: key
   integer (kind=IntKind), intent(in) :: data_size
   integer (kind=IntKind), intent(in) :: data_type
!
   integer (kind=IntKind) :: i
!
   KeyGridGen = KeyGridGen+1
   key = KeyGridGen
!
   if ( NumGrids==0 ) then
      allocate( headGridData )
      tailGridData => headGridData
   else
      allocate (tailGridData%next)
      tailGridData => tailGridData%next
   endif
   nullify(tailGridData%next)
!
   tailGridData%Key = key
   tailGridData%Index = 0
   tailGridData%NumIndexes = 1
   tailGridData%DataSize = data_size
!
   if ( data_type == RealType ) then
      allocate( tailGridData%array_r(data_size) )
   else if (data_type == CmplxType) then
      allocate( tailGridData%array_c(data_size) )
   else
     call ErrorHandler("createGridData","Unknown data type",data_type)
   endif
   tailGridData%DataType = data_type
   NumGrids = NumGrids+1
!
   end subroutine createGridData_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createGridData_ni( ni, key, data_size, data_type )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ni
   integer (kind=IntKind), intent(in) :: data_size(*)
   integer (kind=IntKind), intent(in) :: data_type
   integer (kind=IntKind), intent(out) :: key
!
   integer (kind=IntKind) :: i
!
   KeyGridGen = KeyGridGen+1
   key = KeyGridGen
!
   do i = 1,ni
      if ( NumGrids==0 ) then
         allocate( headGridData )
         tailGridData => headGridData
      else
         allocate (tailGridData%next)
         tailGridData => tailGridData%next
      endif
      nullify(tailGridData%next)
!
      tailGridData%Key = key
      tailGridData%Index = i
      tailGridData%NumIndexes = ni
      tailGridData%DataSize = data_size(i)
!
      if ( data_type == RealType ) then
         allocate( tailGridData%array_r(data_size(i)) )
      else if (data_type == CmplxType) then
         allocate( tailGridData%array_c(data_size(i)) )
      else
        call ErrorHandler("createGridData","Unknown data type",data_type)
      endif
      tailGridData%DataType = data_type
      NumGrids = NumGrids+1
   enddo   
!
   end subroutine createGridData_ni
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createDiscreteData_0( key, data_size, data_type )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(out) :: key
   integer (kind=IntKind), intent(in) :: data_size
   integer (kind=IntKind), intent(in) :: data_type
!
   integer (kind=IntKind) :: i
!
   KeyDiscreteGen = KeyDiscreteGen+1
   key = KeyDiscreteGen
!
   if ( NumDiscreteData==0 ) then
      allocate( headData )
      tailData => headData
   else
      allocate (tailData%next)
      tailData => tailData%next
   endif
   nullify(tailData%next)
!
   tailData%Key = key
   tailData%Index = 0
   tailData%NumIndexes = 1
   tailData%DataSize = 1
!
   if ( data_type == RealType ) then
      allocate( tailData%array_r(data_size) )
   else if (data_type == CmplxType) then
      allocate( tailData%array_c(data_size) )
   else
     call ErrorHandler("createDiscreteData","Unknown data type",data_type)
   endif
   tailData%DataType = data_type
   NumDiscreteData = NumDiscreteData+1
!
   end subroutine createDiscreteData_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine createDiscreteData_ni( ni, key, data_size, data_type )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ni
   integer (kind=IntKind), intent(in) :: data_size(*)
   integer (kind=IntKind), intent(in) :: data_type
   integer (kind=IntKind), intent(out) :: key
!
   integer (kind=IntKind) :: i
!
   KeyDiscreteGen = KeyDiscreteGen+1
   key = KeyDiscreteGen
!
   do i = 1,ni
      if ( NumDiscreteData==0 ) then
         allocate( headData )
         tailData => headData
      else
         allocate (tailData%next)
         tailData => tailData%next
      endif
      nullify(tailData%next)
!
      tailData%Key = key
      tailData%Index = i
      tailData%NumIndexes = ni
      tailData%DataSize = data_size(i)
!
      if ( data_type == RealType ) then
         allocate( tailData%array_r(data_size(i)) )
      else if (data_type == CmplxType) then
         allocate( tailData%array_c(data_size(i)) )
      else
        call ErrorHandler("createDiscreteData", &
                          "Unknown data type",data_type)
      endif
      tailData%DataType = data_type
      NumDiscreteData = NumDiscreteData+1
   enddo   
!
   end subroutine createDiscreteData_ni
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridData_r( key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key, size 
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
!
   real (kind=RealKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchGrids(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getGridData","Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getGridData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getGridData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getGridData","Data size < n",s,size)
   endif
!
   p => temp%array_r(1:size)
!
   end function getGridData_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridData_ri( index, key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) ::  index, key, size 
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
!
   real (kind=RealKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchGrids_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getGridData","Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getGridData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getGridData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getGridData","Data size < n",s,size)
   endif
!
   p => temp%array_r(1:size)
!
   end function getGridData_ri
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridData_c( key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key, size 
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
!
   complex (kind=CmplxKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchGrids(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getGridData","Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getGridData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getGridData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getGridData","Data size < n",s,size)
   endif
!
   p => temp%array_c(1:size)
!
   end function getGridData_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGridData_ci( index, key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) ::  index, key, size 
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
!
   complex (kind=CmplxKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchGrids_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getGridData","Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getGridData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getGridData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getGridData","Data size < n",s,size)
   endif
!
   p => temp%array_c(1:size)
!
   end function getGridData_ci
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDiscreteData_r( key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key, size 
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
!
   real (kind=RealKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchDiscreteData(key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDiscreteData","Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDiscreteData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getDiscreteData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getDiscreteData","Data size < n",s,size)
   endif
!
   p => temp%array_r(1:size)
!
   end function getDiscreteData_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDiscreteData_ri( index, key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) ::  index, key, size 
   integer (kind=IntKind) :: s
!
   real (kind=RealKind), intent(in) :: markup
!
   real (kind=RealKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchDiscreteData_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDiscreteData","Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDiscreteData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getDiscreteData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getDiscreteData","Data size < n",s,size)
   endif
!
   p => temp%array_r(1:size)
!
   end function getDiscreteData_ri
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDiscreteData_r2i(index,key,size1,size2,markup)  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) ::  index, key, size1, size2 
   integer (kind=IntKind) :: s, size
!
   real (kind=RealKind), intent(in) :: markup
!
   real (kind=RealKind), pointer :: p(:,:), p1(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchDiscreteData_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= RealType) then
      call ErrorHandler("getDiscreteData","Data is not of real type")
   else if (s < 1) then
      call ErrorHandler("getDiscreteData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getDiscreteData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getDiscreteData","Data size < n",s,size)
   endif
!
   size = size1*size2
   p1 => temp%array_r(1:size)
   p  => aliasArray2_r(p1,size1,size2)
!
   end function getDiscreteData_r2i
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDiscreteData_c( key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key, size 
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
!
   complex (kind=CmplxKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchDiscreteData(key)
   s = temp%DataSize
!
   if (temp%DataType /= CmplxType) then
      call ErrorHandler("getDiscreteData","Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDiscreteData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getDiscreteData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getDiscreteData","Data size < n",s,size)
   endif
!
   p => temp%array_c(1:size)
!
   end function getDiscreteData_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDiscreteData_ci( index, key, size, markup )  result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) ::  index, key, size 
   integer (kind=IntKind) :: s
!
   complex (kind=CmplxKind), intent(in) :: markup
!
   complex (kind=CmplxKind), pointer :: p(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchDiscreteData_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CmplxType) then
      call ErrorHandler("getDiscreteData","Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDiscreteData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getDiscreteData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getDiscreteData","Data size < n",s,size)
   endif
!
   p => temp%array_c(1:size)
!
   end function getDiscreteData_ci
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDiscreteData_c2i(index,key,size1,size2,markup) result(p)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) ::  index, key, size1, size2 
   integer (kind=IntKind) :: s, size
!
   complex (kind=CmplxKind), intent(in) :: markup
!
   complex (kind=CmplxKind), pointer :: p(:,:), p1(:)
!
   type (AnaliticStruct), pointer :: temp
!
   temp => searchDiscreteData_i(index,key)
   s = temp%DataSize
!
   if (temp%DataType /= CmplxType) then
      call ErrorHandler("getDiscreteData","Data is not of complex type")
   else if (s < 1) then
      call ErrorHandler("getDiscreteData","Data size is 0",s)
   else if (size < 1) then
      call ErrorHandler("getDiscreteData","Invalid dimension size < n",size)
   else if (size > s) then
      call ErrorHandler("getDiscreteData","Data size < n",s,size)
   endif
!
   size = size1*size2
   p1 => temp%array_c(1:size)
   p => aliasArray2_c(p1,size1,size2)
!
   end function getDiscreteData_c2i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteGrid(key)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind) :: i, ni
!
   type (AnaliticStruct), pointer :: temp1
   type (AnaliticStruct), pointer :: temp2
!
   temp2 => searchGrids(key)
!
   if (temp2%DataSize < 1) then
      call ErrorHandler("deleteGrid","Data Size < 1",temp2%DataSize)
   endif
!  
   ni = temp2%NumIndexes
   do i = 1,ni
      if (temp2%DataType == RealType) then
         deallocate( temp2%array_r )
      else if (temp2%DataType == CmplxType) then
         deallocate( temp2%array_c )
      endif
!
      if (NumGrids == 1) then
         nullify( headGridData )
         nullify( tailGridData )
      else
         temp1 => searchPrevGrid(key)
         if ( associated(temp1) )then
            temp1%next => temp2%next
         else
            headGridData => temp2%next
         endif
         if (.not.associated(temp2%next) )then
            tailGridData => temp1
         endif
      endif
      nullify( temp2%next )
      deallocate( temp2 )
      temp2 => temp1%next
   enddo
!
   NumGrids = NumGrids - 1
!
   end subroutine deleteGrid
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteDiscreteData(key)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind) :: i, ni
!
   type (AnaliticStruct), pointer :: temp1
   type (AnaliticStruct), pointer :: temp2
!
   temp2 => searchDiscreteData(key)
!
   if (temp2%DataSize < 1) then
      call ErrorHandler("deleteGrid","Data Size < 1",temp2%DataSize)
   endif
!
   ni = temp2%NumIndexes
   do i = 1,ni
      if (temp2%DataType == RealType) then
         deallocate( temp2%array_r )
      else if (temp2%DataType == CmplxType) then
         deallocate( temp2%array_c )
      endif
!
      if (NumDiscreteData == 1) then
         nullify( headData )
         nullify( tailData )
      else
         temp1 => searchPrevDiscreteData(key)
         if ( associated(temp1) )then
            temp1%next => temp2%next
         else
            headData => temp2%next
         endif
         if (.not.associated(temp2%next) )then
            tailData => temp1
         endif
      endif
      nullify( temp2%next )
      deallocate( temp2 )
      temp2 => temp1%next
   enddo
!
   NumDiscreteData = NumDiscreteData - 1
!
   end subroutine deleteDiscreteData
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchGrids(key) result(temp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind) :: i
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumGrids < 1) then
      call ErrorHandler("searchGrids","Empty data grids")
   endif
!
   temp => headGridData
   do i=1,NumGrids
      if (temp%Key == key) then
         return
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchGrids","Data with the key is not found",key)
!
   end function searchGrids
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchGrids_i(index,key) result(temp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: i, j, n
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumGrids < 1) then
      call ErrorHandler("searchGrids_i","Empty data grids")
   else if (index < 0) then
      call ErrorHandler("searchGrids_i","Invalid index",index)
   endif
!
   temp => headGridData
   do i=1,NumGrids
      if (temp%Key == key) then
         n = temp%NumIndexes
         if (index == 0 .and. n == 0) then
            return
         else if (index == 0 .and. n > 0) then
            call ErrorHandler("searchGrids_i","Invalid index",index)
         else if (index > n) then
            call ErrorHandler("searchGrids_i","Invalid index",index)
         else
            do j=1,n
               if (temp%Index == index) then
                  return
               else
                  temp => temp%next
               endif
            enddo
         endif
         call ErrorHandler("searchGrids_i",                         &
                           "Data with the index is not found",index)
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchGrids_i","Data with the key is not found",key)
!
   end function searchGrids_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchDiscreteData(key) result(temp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind) :: i
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumDiscreteData < 1) then
      call ErrorHandler("searchGrids","Empty data grids")
   endif
!
   temp => headData
   do i=1,NumDiscreteData
      if (temp%Key == key) then
         return
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchGrids","Data with the key is not found",key)
!
   end function searchDiscreteData
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchDiscreteData_i(index,key) result(temp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind), intent(in) :: index
   integer (kind=IntKind) :: i, j, n
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumDiscreteData < 1) then
      call ErrorHandler("searchDiscreteData_i","Empty data grids")
   else if (index < 0) then
      call ErrorHandler("searchDiscreteData_i","Invalid index",index)
   endif
!
   temp => headData
   do i = 1,NumDiscreteData
      if (temp%Key == key) then
         n = temp%NumIndexes
         if (index == 0 .and. n == 0) then
            return
         else if (index == 0 .and. n > 0) then
            call ErrorHandler("searchDiscreteData_i","Invalid index",index)
         else if (index > n) then
            call ErrorHandler("searchDiscreteData_i","Invalid index",index)
         else
            do j=1,n
               if (temp%Index == index) then
                  return
               else
                  temp => temp%next
               endif
            enddo
         endif
         call ErrorHandler("searchDiscreteData_i",                         &
                           "Data with the index is not found",index)
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchDiscreteData_i","Data with the key is not found",key)
!
   end function searchDiscreteData_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchPrevGrid(key) result(temp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind) :: i
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumGrids < 1) then
      call ErrorHandler("searchPrevGrid","Empty data grid")
   else if (NumGrids == 1) then
      nullify( temp )
      return
   endif
!
   temp => headGridData
   do i=1,NumGrids
      if (temp%next%Key == key) then
         return
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchPrevGrid","Data with the key is not found",key)
!
   end function searchPrevGrid
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function searchPrevDiscreteData(key) result(temp)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: key
!
   integer (kind=IntKind) :: i
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumGrids < 1) then
      call ErrorHandler("searchPrevDiscreteData","Empty discrete data ")
   else if (NumGrids == 1) then
      nullify( temp )
      return
   endif
!
   temp => headData
   do i=1,NumDiscreteData
      if (temp%next%Key == key) then
         return
      else
         temp => temp%next
      endif
   enddo
!
   call ErrorHandler("searchPrevDiscreteData", & 
                     "Data with the key is not found",key)
!
   end function searchPrevDiscreteData
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isGridExisting(key) result(t)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: key
!
   logical :: t
!
   integer (kind=IntKind) :: i
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumGrids > 0) then
      temp => headGridData
      do i=1,NumGrids
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
   end function isGridExisting
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isDataExisting(key) result(t)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: key
!
   logical :: t
!
   integer (kind=IntKind) :: i
!
   type (AnaliticStruct), pointer :: temp
!
   if (NumDiscreteData > 0) then
      temp => headData
      do i=1,NumDiscreteData
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
   end function isDataExisting
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine fitData_r(GridID, DataID, DataSize, DataIndex, array,asize)
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: GridId, asize
   integer(kind=IntKind), intent(in) :: DataID, DataIndex, DataSize
   real(kind=RealKind), target :: array(asize)
!
   
!
   real(kind=RealKind), pointer :: parray(:)
!
   parray => array
!
   end subroutine fitData_r
!  ===================================================================   
end module AnaliticContinuationModule
