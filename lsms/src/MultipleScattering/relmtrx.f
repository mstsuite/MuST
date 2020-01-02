      subroutine relmtrx(a,b,kkr1,kkr2)
!     subroutine relmtrx(a,b,lmax)
c=======================
c
c *****************************************************
c * transformation from a non-relativistic matrix 'a' *
c *                to a relativistic matrix 'b'       *
c *****************************************************
c
      implicit real*8 (a-h,o-z)
!     include '../param.h'
      include '../Misc/cgc.h'
c
      complex*16 a,b
c
      dimension a(kkr1,kkr2)
!     dimension a((lmax+1)*(lmax+1),(lmax+1)*(lmax+1))
      dimension b(2*kkr1,2*kkr2)
!     dimension b(2*(lmax+1)*(lmax+1),2*(lmax+1)*(lmax+1))
!     dimension u1(50),ind1(50)
!     dimension u2(50),ind2(50)
!     common/cgc/u1,u2,ind1,ind2
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      do i=1,2*kkr1
!     do i=1,kmymax
        i1=ind1(i)
        i2=ind2(i)
        do j=1,2*kkr2
!       do j=1,kmymax
          j1=ind1(j)
          j2=ind2(j)
!         write(6,*) 'i,j',i,j
          b(i,j)=u1(i)*a(i1,j1)*u1(j)+
     >           u2(i)*a(i2,j2)*u2(j)
        end do
      end do
c
      return
      end
