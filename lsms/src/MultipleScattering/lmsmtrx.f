      subroutine lmsmtrx(a,b,kmymax)
! transform a matrix B in kappa,my representation into
!  a matrix A in l,m,s representation
      implicit none

      include '../cgc.h'

      integer kmymax
      integer kmy1
      integer kmy2
      integer lm11,lm12
      integer lm21,lm22
!     integer ind1(50)
!     integer ind2(50)
      integer soff

      complex*16 a(kmymax,kmymax)
      complex*16 b(kmymax,kmymax)
!     real*8     u1(50)
!     real*8     u2(50)

!     common/cgc/u1,u2,ind1,ind2

      soff = kmymax/2
      call zeroout(a,2*kmymax*kmymax)
      do kmy1=1,kmymax
        lm11 = ind1(kmy1)
        lm12 = ind2(kmy1)+soff
        do kmy2=1,kmymax
          lm21 = ind1(kmy2)
          lm22 = ind2(kmy2)+soff

          a(lm11,lm21) = a(lm11,lm21)+u1(kmy1)*b(kmy1,kmy2)*u1(kmy2)
          a(lm12,lm22) = a(lm12,lm22)+u2(kmy1)*b(kmy1,kmy2)*u2(kmy2)
          a(lm11,lm22) = a(lm11,lm22)+u1(kmy1)*b(kmy1,kmy2)*u2(kmy2)
          a(lm12,lm21) = a(lm12,lm21)+u2(kmy1)*b(kmy1,kmy2)*u1(kmy2)

        end do
      end do

      end



