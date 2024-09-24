subroutine printxml_code_info(exec_date,exec_time,version)
   use KindParamModule, only : IntKind
   use XMLDataModule, only : writeStartTag, writeEndTag, writeElement
!
   implicit none
!
   character (len=*), intent(in) :: exec_date, exec_time, version
   character (len=10) :: start_date
   character (len=8)  :: start_time
!
   start_date = exec_date(5:6)//'-'//exec_date(7:8)//'-'//exec_date(1:4)
   start_time = exec_time(1:2)//':'//exec_time(3:4)//':'//exec_time(5:6)
   call writeStartTag(tag='code_info')
      call writeElement(tag='i',ename='program',etype='string',econtent='MuST')
      call writeElement(tag='i',ename='version',etype='string',econtent=version)
      call writeElement(tag='i',ename='date',etype='string',econtent=start_date)
      call writeElement(tag='i',ename='time',etype='string',econtent=start_time)
   call writeEndTag(tag='code_info')
!
end subroutine printxml_code_info
