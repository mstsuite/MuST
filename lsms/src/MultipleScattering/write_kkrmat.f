      subroutine write_kkrmat(a,n,lda,e)
      implicit none

      integer n,lda,i,j
      complex*16 a(lda,*),e

      open(unit=42,file="kkrmat.out")

      write(42,*) n,real(e),imag(e)

      do i=1,n
      do j=1,n

      write(42,*) i,j,real(a(i,j)),imag(a(i,j))

      end do
      end do

      close(42)
      end
