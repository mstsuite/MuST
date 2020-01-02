      subroutine gen_test_mat(a,lda,na)
      implicit none
      integer lda,na
      integer i,j
      complex*16 a(lda,*)
      real*8 rr,ri


      do j = 1, na
        do i = 1, na
          rr = (11.0d0*i-3.0d0*j)/real(i+j)
          ri = (5.0d0*i-2.0d0*j)/real(i+j)
          a(i,j) = cmplx(rr,ri)
        end do
      end do

      end subroutine
