subroutine printxml_scf_param()
   use XMLDataModule, only : writeStartTag, writeEndTag, writeElement
   use ExchCorrFunctionalModule, only : getFunctionalType
!
   implicit none
!
   character (len=10) :: xc_type
!
   xc_type = getFunctionalType()
!
   call writeStartTag(tag='scf_param')
      call writeStartTag(tag='separator',ename='electronic exchange-correlation')
         call writeElement(tag='i',ename=xc_type,etype='string')
      call writeEndTag(tag='separator')
   call writeEndTag(tag='scf_param')
!
end subroutine printxml_scf_param
