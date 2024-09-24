program testXMLData
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use XMLDataModule, only : initXMLData, endXMLData, writeStartTag, writeEndTag, writeElement
   implicit none
!
   character (len=20) :: attr_name(10), attr_val(10)
   character (len=255) :: dir
!
   integer (kind=IntKind) :: n
!
   CALL getenv("PWD", dir)
   WRITE (6,'(2a)')'PWD: ', TRIM(dir)
   call getenv("ODIR",dir)
   WRITE (6,'(2a)')'ODIR:', TRIM(dir)
!
   call initXMLData('mstrun.xml',root='mst_run')
!
   call writeStartTag(tag='code_info')
      call writeElement(tag='i',ename='program',etype='string',econtent='MuST')
      call writeElement(tag='i',ename='version',etype='string',econtent='1.9.3')
      call writeElement(tag='i',ename='date',etype='string',econtent='2024-09-16')
      call writeElement(tag='i',ename='time',etype='string',econtent='10:07:44')
   call writeEndTag(tag='code_info')
!
   call writeStartTag(tag='parameters')
      call writeStartTag(tag='separator',ename='electronic exchange-correlation')
         call writeElement(tag='i',ename='LDA',etype='string')
      call writeEndTag(tag='separator')
   call writeEndTag(tag='parameters')
!
   call writeStartTag(tag='scf_results')
!     ===========================================
      call writeStartTag(tag='varray',ename='forces')
         call writeElement(tag='v',econtent='      -0.04590874      -0.04590874      -0.00000000')
         call writeElement(tag='v',econtent='       0.02501318       0.02501318      -0.00000000')
         call writeElement(tag='v',econtent='      -0.00048468       0.01176199       0.01320202')
         call writeElement(tag='v',econtent='       0.01176199      -0.00048468      -0.01320202')
         call writeElement(tag='v',econtent='       0.01176199      -0.00048468       0.01320202')
         call writeElement(tag='v',econtent='      -0.00048468       0.01176199      -0.01320202')
         call writeElement(tag='v',econtent='      -0.00142831      -0.00142831      -0.00000000')
         call writeElement(tag='v',econtent='      -0.00023075      -0.00023075      -0.00000000')
      call writeEndTag(tag='varray')
!     ===========================================
      call writeStartTag(tag='energy')
         call writeElement(tag='i',ename='energy_band',econtent='    -29.79405929')
         call writeElement(tag='i',ename='energy_Hartree',econtent='    -29.79631394')
         call writeElement(tag='i',ename='energy_xc',econtent='    -29.79631394')
      call writeEndTag(tag='energy')
!     ===========================================
   call writeEndTag(tag='scf_results')
!
   call endXMLData()
!
end program testXMLData
