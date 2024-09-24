module XMLDataModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: initXMLData,     &
          writeStartTag,   &
          writeEndTag,     &
          writeElement,    &
          endXMLData
!
private
   character (len=120) :: blank = ' '
   character (len=50) :: root_tag
!
   integer (kind=IntKind), parameter :: xml_unit = 51
   integer (kind=IntKind) :: spacing = 0
!
   interface writeStartTag
      module procedure writeStartTag_attr, writeStartTag_noattr
   end interface writeStartTag
!
   interface writeElement
      module procedure writeElement_0, writeElement_1
   end interface writeElement
contains
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initXMLData(xml_fname,root)
!  ==================================================================
   implicit none
!
   character (len=*), intent(in) :: xml_fname
   character (len=*), intent(in) :: root
!
   open(unit=xml_unit,file=adjustl(trim(xml_fname)),status='unknown',form='formatted')
!
   write(xml_unit,'(a)')'<?xml version="1.0" encoding="UTF-8"?>'
!
   spacing = 0
   root_tag = adjustl(trim(root))
!
   call writeStartTag(root_tag)
!
   end subroutine initXMLData
!  ==================================================================
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endXMLData()
!  ==================================================================
   implicit none
!
   call writeEndTag(root_tag)
!
   if (spacing /= 0) then
      call ErrorHandler('endXMLData','A XML tag is not closed in mstrun.xml')
   endif
!
   close(unit=xml_unit)
!
   end subroutine endXMLData
!  ==================================================================
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeStartTag_noattr(tag,ename)
!  ==================================================================
   implicit none
!   
   character (len=*), intent(in) :: tag
   character (len=*), intent(in), optional :: ename
!
   if (spacing > 0) then
      write(xml_unit,'(a,$)')blank(1:spacing)
   endif
!
   if (present(ename)) then
      write(xml_unit,'(5a)')'<',trim(tag),' name="',trim(ename),'">'
   else
      write(xml_unit,'(3a)')'<',trim(tag),'>'
   endif
!
   spacing = spacing + 3
!
   end subroutine writeStartTag_noattr
!  ==================================================================
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeStartTag_attr(tag,num_attr,attr_name,attr_val)
!  ==================================================================
   implicit none
!   
   character (len=*), intent(in) :: tag
!
   integer (kind=IntKind), intent(in) :: num_attr
   integer (kind=IntKind) :: i
!
   character (len=*), intent(in) :: attr_name(:)
   character (len=*), intent(in) :: attr_val(:)
!
   if (spacing > 0) then
      write(xml_unit,'(a,$)')blank(1:spacing)
   endif
!
   write(xml_unit,'(2a,$)')'<',trim(tag)
   do i = 1, num_attr
      write(xml_unit,'(4a,$)')' ',trim(attr_name(i)),'="',trim(attr_val(i)),'"'
   enddo
   write(xml_unit,'(a)')'>'
!
   spacing = spacing + 3
!
   end subroutine writeStartTag_attr
!  ==================================================================
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeEndTag(tag)
!  ==================================================================
   implicit none
!   
   character (len=*), intent(in) :: tag
!
   spacing = spacing - 3
!
   if (spacing > 0) then
      write(xml_unit,'(4a)')blank(1:spacing),'</',trim(tag),'>'
   else
      write(xml_unit,'(3a)')'</',trim(tag),'>'
   endif
!
   end subroutine writeEndTag
!  ==================================================================
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeElement_0(tag,ename,etype,econtent)
!  ==================================================================
   implicit none
!   
   character (len=*), intent(in) :: tag
   character (len=*), intent(in), optional :: ename
   character (len=*), intent(in), optional :: etype
   character (len=*), intent(in), optional :: econtent
!
   if (spacing > 0) then
      write(xml_unit,'(a,$)')blank(1:spacing)
   endif
   write(xml_unit,'(2a,$)')'<',trim(tag)
   if (present(ename)) then
      write(xml_unit,'(3a,$)')' name="',trim(ename),'"'
   endif
   if (present(etype)) then
      write(xml_unit,'(3a,$)')' type="',trim(etype),'"'
   endif
   if (present(econtent)) then
      write(xml_unit,'(5a)')'>',trim(econtent),' </',trim(tag),'>'
   else
      write(xml_unit,'(4a)')'>',' </',trim(tag),'>'
   endif
!
   end subroutine writeElement_0
!  ==================================================================
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeElement_1(tag,num_attr,attr_name,attr_val,econtent)
!  ==================================================================
   implicit none
!   
   character (len=*), intent(in) :: tag
!
   integer (kind=IntKind), intent(in) :: num_attr
   integer (kind=IntKind) :: i
!
   character (len=*), intent(in) :: attr_name(:)
   character (len=*), intent(in) :: attr_val(:)
   character (len=*), intent(in), optional :: econtent
!
   if (spacing > 0) then
      write(xml_unit,'(a,$)')blank(1:spacing)
   endif
   write(xml_unit,'(2a,$)')'<',trim(tag)
   do i = 1, num_attr
      write(xml_unit,'(4a,$)')' ',trim(attr_name(i)),'="',trim(attr_val(i)),'"'
   enddo
   if (present(econtent)) then
      write(xml_unit,'(5a)')'>',trim(econtent),' </',trim(tag),'>'
   else
      write(xml_unit,'(4a)')'>',' </',trim(tag),'>'
   endif
!
   end subroutine writeElement_1
!  ==================================================================
end module XMLDataModule
