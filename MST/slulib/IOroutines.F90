!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeBinaryMatrix(matrix, dim1, dim2, id)
!  =================================================================
!
      use KindParamModule, only : IntKind, CmplxKind
!
      implicit none 
!
      integer (kind=IntKind), intent(in)   :: dim1, dim2, id
!
      complex (kind=CmplxKind), intent(in) :: matrix(dim1,dim2)
!
      character (len=20) :: outfile
      character (len=5)  :: str
!     ==============================================================
!     create output file names for each id(node)
!     ==============================================================
      write(str,'(f4.3)') id/1000.0
      outfile = 'matrix.bin'//str
!     ==============================================================
!     write matrix data
!     ==============================================================
      open(unit=30,file=outfile,status='unknown',form='unformatted')
      write(30) matrix
      close(30)
!
      return
!     ==============================================================
   end subroutine writeBinaryMatrix
!  =================================================================
!
!  *****************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine readBinaryMatrix(matrix, dim1, dim2, id)
!  =================================================================
!
      use KindParamModule, only : IntKind, CmplxKind
!
      implicit none 
!
      integer (kind=IntKind), intent(in)   :: dim1, dim2, id
!
      complex (kind=CmplxKind), intent(out) :: matrix(dim1,dim2)
!
      character (len=20) :: infile
      character (len=5)  :: str
!     ==============================================================
!     create the name of input file for id(node)
!     ==============================================================
      write(str,'(f4.3)') id/1000.0
      infile = 'matrix.bin'//str
!     ==============================================================
!     write matrix data
!     ==============================================================
      open(unit=30,file=infile,status='old',form='unformatted')
      read(30) matrix
      close(30)
!
      return
!     ==============================================================
   end subroutine readBinaryMatrix
!  =================================================================   