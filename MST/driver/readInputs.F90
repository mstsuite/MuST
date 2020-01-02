subroutine readInputs(def_id,info_id)
   use KindParamModule, only : IntKind
   use PublicParamDefinitionsModule, only : StandardInputFile
   use InputModule, only : readInputData
   use InputModule, only : openInputFile, closeInputFile
   use InputModule, only : getKeyValue, getTableIndex
   use InputModule, only : getNumKeys
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   character (len=80) :: info_table, info_path
   character (len=160) :: itname, FileName
!
   logical :: FileNamed = .false.
   logical :: StandardInputExist = .false.
!
   integer (kind=IntKind), intent(out) :: def_id
   integer (kind=IntKind), intent(out) :: info_id
!
   integer (kind=IntKind) :: fstatus
!  
!  inquire(unit=5,name=FileName,named=FileNamed)
   inquire(unit=5,name=FileName,exist=StandardInputExist)
!  if (FileNamed) then
   if (StandardInputExist) then
!     ----------------------------------------------------------------
      call readInputData(5,def_id)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call openInputFile(7,StandardInputFile)
      call readInputData(7,def_id)
!     ----------------------------------------------------------------
   endif
!  ===================================================================
!  call readInput to obtain input data................................
!      Proc 0: read input data file and broadcast to other processors
!      Other : wait and receive the data from Proc 0.
!  -------------------------------------------------------------------
!! call readInputData(5,def_id)
!  -------------------------------------------------------------------
!! fstatus =  getKeyValue(def_id,'Info Table File Name',info_table)
!! fstatus =  getKeyValue(def_id,'Current File Path',info_path)
!  -------------------------------------------------------------------
!! itname = trim(info_path)//info_table
!! info_id = getTableIndex(itname)
!! if (info_id < 1) then
!     ----------------------------------------------------------------
!!    call openInputFile(10,itname)
!!    call readInputData(10,info_id)
!!    call closeInputFile(10)
!     ----------------------------------------------------------------
!! endif
!
   if (getKeyValue(def_id,'Current File Path',info_path) /= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('main','Input parameter for Current File Path is not found!')
!     ----------------------------------------------------------------
   endif
!
   if (getKeyValue(def_id,'Info Table File Name',info_table) == 0) then
      itname = trim(info_path)//info_table
      info_id = getTableIndex(itname)
      if (info_id < 1) then
!        -------------------------------------------------------------
         call openInputFile(10,itname)
         call readInputData(10,info_id)
!        -------------------------------------------------------------
      endif
   else
      info_id = def_id
   endif
!
end subroutine readInputs
