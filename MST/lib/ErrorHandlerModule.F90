module ErrorHandlerModule
!
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use KindParamModule, only : CmplxKind
!
!  define generic procedure for ErrorHandler, StopHandler, WarningHandler
!         and MessageHandler
!
public :: ErrorHandler, &
          StopHandler,  &
          WarningHandler, &
          MessageHandler, &
          setErrorOutput
!
   interface ErrorHandler
      module procedure ErrorHandlerInt  , ErrorHandlerInt2  , &
                       ErrorHandlerInt3 , ErrorHandlerInt4  , &
                       ErrorHandlerReal , ErrorHandlerReal2 , &
                       ErrorHandlerCmplx, ErrorHandlerCmplx2, & 
                       ErrorHandlerChar , ErrorHandlerChar2 , &
                       ErrorHandlerMsg  , ErrorHandler0msg
   end interface
!
   interface StopHandler
      module procedure StopHandlerInt  , StopHandlerInt2  , &
                       StopHandlerReal , StopHandlerReal2 , &
                       StopHandlerCmplx, StopHandlerCmplx2, & 
                       StopHandlerChar , StopHandlerChar2 , &
                       StopHandlerMsg  , StopHandler0msg
   end interface
!
   interface WarningHandler
      module procedure WarningHandlerInt  , WarningHandlerInt2  , &
                       WarningHandlerReal , WarningHandlerReal2 , &
                       WarningHandlerCmplx, WarningHandlerCmplx2, & 
                       WarningHandlerChar , WarningHandlerChar2 , &
                       WarningHandlerMsg  , WarningHandler0msg
   end interface
!
   interface MessageHandler
      module procedure MessageHandlerInt  , MessageHandlerInt2  , &
                       MessageHandlerReal , MessageHandlerReal2 , &
                       MessageHandlerCmplx, MessageHandlerCmplx2, & 
                       MessageHandlerChar , MessageHandlerChar2 , &
                       MessageHandlerMsg  , MessageHandler0msg
   end interface
!
private
!
   integer (kind=IntKind) :: FileUnit = 6
   integer (kind=IntKind) :: ErrorPrintLevel = 1
   integer (kind=IntKind) :: WarningPrintLevel = 1
   integer (kind=IntKind) :: MessagePrintLevel = 1
   integer (kind=IntKind), parameter :: line_width = 60
!
   character (len=line_width), pointer :: msg_ml(:)
   integer (kind=IntKind) :: num_lines
!
   interface
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine breakLine(s,lw,nl,sml)
!     ================================================================
      use KindParamModule, only : IntKind
!
      integer (kind=IntKind), intent(in) :: lw
      integer (kind=IntKind), intent(out) :: nl
!
      character (len=*), intent(in) :: s
      character (len=*), pointer, intent(out) :: sml(:)
!
      end subroutine breakLine
!     ================================================================
   end interface
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setErrorOutput(e_print_level,w_print_level,m_print_level,funit)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in), optional :: funit
   integer (kind=IntKind), intent(in), optional :: e_print_level
   integer (kind=IntKind), intent(in), optional :: w_print_level
   integer (kind=IntKind), intent(in), optional :: m_print_level
!
   if (present(e_print_level)) then
      ErrorPrintLevel = e_print_level
   endif
!
   if (present(w_print_level)) then
      WarningPrintLevel = w_print_level
   endif
!
   if (present(m_print_level)) then
      MessagePrintLevel = m_print_level
   endif
!
   if (present(funit)) then
      FileUnit = funit
   else
      FileUnit = 6
   endif
!
   end subroutine setErrorOutput
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerInt(rname,error_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   integer (kind=IntKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,i8)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerInt
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerInt2(rname,error_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   integer (kind=IntKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,2(i8,2x))')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerInt2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerInt3(rname,error_message,value1,value2,value3,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   integer (kind=IntKind), intent(in) :: value1,value2,value3
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,3(i8,2x))')'VALUE: ',value1,value2,value3
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerInt3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerInt4(rname,error_message,value1,value2, &
                               value3,value4,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   integer (kind=IntKind), intent(in) :: value1,value2,value3,value4
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,4(i8,2x))')'VALUE: ',value1,value2,value3,value4
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerInt4
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerReal(rname,error_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   real (kind=RealKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,d18.7)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerReal
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerReal2(rname,error_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   real (kind=RealKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,d18.7,2x,d18.7)')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerReal2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerCmplx(rname,error_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   complex (kind=CmplxKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,2d18.7)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerCmplx
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerCmplx2(rname,error_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   complex (kind=CmplxKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(a,2d18.7,2x,2d18.7)')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerCmplx2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerChar(rname,error_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message,value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(2a)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerChar
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerChar2(rname,error_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message,value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(3a)')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerChar2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandlerMsg(rname,error_message,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,error_message
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'ERROR: ',error_message
      call breakLine(error_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'ERROR: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'       ',msg_ml(i)
      enddo
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandlerMsg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ErrorHandler0msg(rname,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname
   logical :: pr
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (ErrorPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(a)')'ERROR: '
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine ErrorHandler0msg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerInt(rname,stop_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message
   integer (kind=IntKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(a,i8)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerInt
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerInt2(rname,stop_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message
   integer (kind=IntKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(a,i8,2x,i8)')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerInt2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerReal(rname,stop_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message
   real (kind=RealKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(a,d18.7)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerReal
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerReal2(rname,stop_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message
   real (kind=RealKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(a,d18.7,2x,d18.7)')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerReal2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerCmplx(rname,stop_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message
   complex (kind=CmplxKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(a,2d18.7)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerCmplx
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerCmplx2(rname,stop_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message
   complex (kind=CmplxKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(a,2d18.7,2x,2d18.7)')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerCmplx2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerChar(rname,stop_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message,value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(2a)')'VALUE: ',value
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerChar
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerChar2(rname,stop_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message,value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(3a)')'VALUE: ',value1,value2
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerChar2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandlerMsg(rname,stop_message,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,stop_message
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
!     write(FileUnit,'(2a)')'MESSAGE: ',stop_message
      call breakLine(stop_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a)')'MESSAGE: ',msg_ml(1)
      do i = 2, num_lines
         write(FileUnit,'(a,a)')'         ',msg_ml(i)
      enddo
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandlerMsg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine StopHandler0msg(rname,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname
   logical :: pr
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(a)')'MESSAGE: '
      write(FileUnit,'(2a)')'STOP AT ',rname
   endif
!
   stop
   end subroutine StopHandler0msg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerInt(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   integer (kind=IntKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,i8,t79,a1)')' ! VALUE: ',value,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerInt
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerInt2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   integer (kind=IntKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,i8,2x,i8,t79,a1)')' ! VALUE: ',value1,value2,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerInt2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerReal(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   real (kind=RealKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,d20.13,t79,a1)')' ! VALUE: ',value,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerReal
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerReal2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   real (kind=RealKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,d20.13,2x,d20.13,t79,a1)')' ! VALUE: ',value1,value2,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerReal2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerCmplx(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   complex (kind=CmplxKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,2d15.8,t79,a1)')' ! VALUE: ',value,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerCmplx
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerCmplx2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   complex (kind=CmplxKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,2d15.8,2x,2d15.8,t79,a1)')' ! VALUE: ',value1,value2,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerCmplx2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerChar(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message,value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,a,t79,a1)')' ! VALUE: ',trim(value),'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerChar
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerChar2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message,value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,a,a,t79,a1)')' ! VALUE: ',trim(value1),trim(value2),'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerChar2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandlerMsg(rname,warning_message,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! WARNING: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandlerMsg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine WarningHandler0msg(rname,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname
   logical :: pr
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (WarningPrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(a,t79,a1)')' ! WARNING: ','!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine WarningHandler0msg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerInt(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   integer (kind=IntKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,i8,t79,a1)')' ! VALUE: ',value,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerInt
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerInt2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   integer (kind=IntKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,i8,2x,i8,t79,a1)')' ! VALUE: ',value1,value2,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerInt2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerReal(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   real (kind=RealKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,d20.13,t79,a1)')' ! VALUE: ',value,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerReal
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerReal2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   real (kind=RealKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,d20.13,2x,d20.13,t79,a1)')' ! VALUE: ',value1,value2,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerReal2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerCmplx(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   complex (kind=CmplxKind), intent(in) :: value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,2d15.8,t79,a1)')' ! VALUE: ',value,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerCmplx
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerCmplx2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   complex (kind=CmplxKind), intent(in) :: value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,2d15.8,2x,2d15.8,t79,a1)')' ! VALUE: ',value1,value2,'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerCmplx2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerChar(rname,warning_message,value,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message,value
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,a,t79,a1)')' ! VALUE: ',trim(value),'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerChar
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerChar2(rname,warning_message,value1,value2,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message,value1,value2
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
!     write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',warning_message,'!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,a,a,t79,a1)')' ! VALUE: ',trim(value1),trim(value2),'!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerChar2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandlerMsg(rname,warning_message,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname,warning_message
   logical :: pr
   integer (kind=IntKind) :: i
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      call breakLine(warning_message,line_width,num_lines,msg_ml)
      write(FileUnit,'(a,a,t79,a1)')' ! MESSAGE: ',msg_ml(1),'!'
      do i = 2, num_lines
         write(FileUnit,'(a,a,t79,a1)')' !          ',msg_ml(i),'!'
      enddo
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandlerMsg
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine MessageHandler0msg(rname,force_to_print)
!  ===================================================================
   implicit none
   character (len=*), intent(in) :: rname
   logical :: pr
!
   logical, optional, intent(in) :: force_to_print
!
   if (present(force_to_print)) then
      pr = force_to_print
   else if (MessagePrintLevel > 0) then
      pr = .true.
   else
      pr = .false.
   endif
!
   if (pr) then
      write(FileUnit,'(/,1x,78(''!''))')
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(a,t79,a1)')' ! MESSAGE: ','!'
      write(FileUnit,'(a,a)') ' ! AT: ',rname
      write(FileUnit,'(1x,a1,t79,a1)')'!','!'
      write(FileUnit,'(1x,78(''!''),/)')
   endif
!
   end subroutine MessageHandler0msg
!  ===================================================================
!
end module ErrorHandlerModule
