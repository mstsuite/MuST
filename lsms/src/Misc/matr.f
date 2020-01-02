      subroutine matr(tlm,lmax,dmat,dmat1)
! ===================
!
! R=(tvec,phi)
! dmat = D(R) and dmat1 = D(R)+ 
!
      implicit real*8 (a-h,o-z)
!     include '../param.h'
!     parameter(jdim=2*lmaxp+2)
!
      dimension tlm(4)
      complex*16 dmat(2*(lmax+1)*(lmax+1),2*(lmax+1)*(lmax+1))
      complex*16 dmat1(2*(lmax+1)*(lmax+1),2*(lmax+1)*(lmax+1))
      complex*16 cmat(2*(lmax+1)*(lmax+1),2*(lmax+1)*(lmax+1))
      complex*16 d(2*lmax+2,2*lmax+2)
!     common/test/itest
!
      kmax=2*lmax+1
      kmymax=2*(lmax+1)*(lmax+1)
!
! Set up matrix of rotation
!
      do i=1,kmymax
      do j=1,kmymax
        dmat(i,j)=(0.d0,0.d0)
      end do
      end do
      ist=0
      do j2=1,2*lmax-1,2
        il=j2+1
        call rotmat(d,trd,il,tlm,1,lmax)
        do icase=1,2
          do m1=1,il
          do m2=1,il
            dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
          end do
          end do
          ist=ist+il
        end do
      end do
      j2=2*lmax+1
      il=j2+1
      call rotmat(d,trd,il,tlm,1,lmax)
      do m1=1,il
      do m2=1,il
! fill up matrix according to ascending m index
        dmat(ist+m1,ist+m2)=d(il+1-m1,il+1-m2)
      end do
      end do
!
      do i=1,kmymax
      do j=1,kmymax
        dmat1(i,j)=dconjg(dmat(j,i))
      end do
      end do
!
c$$$  if(0.lt.4) return
c$$$  write(6,'(/'' Matrix of rotation'')') 
c$$$  call outmat(dmat,kmymax,kmymax,2*(lmax+1)*(lmax+1),6)
c$$$  write(6,'(/'' Inverse'')') 
c$$$  call outmat(dmat1,kmymax,kmymax,2*(lmax+1)*(lmax+1),6)
c$$$  call repl(cmat,dmat,kmymax,2*(lmax+1)*(lmax+1))
c$$$  call doubmt(cmat,dmat1,kmymax,2*(lmax+1)*(lmax+1))
c$$$  write(6,'(/'' D(R) * D(R**-1) '')')
c$$$  call outmat(cmat,kmymax,kmymax,2*(lmax+1)*(lmax+1),6)
!
c$$$  return
      end
