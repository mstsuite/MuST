      subroutine replms(alms,akmy,lmax,kmymax)
c======================
c
c find diagonal elements of a matrix in (l,m,s) representation
c to be given in (kappa,my) representation
c
      implicit real*8 (a-h,o-z)
c
!     include '../param.h'
      include '../cgc.h'
c
      complex*16 alms(kmymax),akmy(kmymax,kmymax)
!     dimension u1(50),ind1(50)
!     dimension u2(50),ind2(50)
!     common/cgc/u1,u2,ind1,ind2
c
      lmmax=(lmax+1)*(lmax+1)
      kmymax=2*lmmax
c
      lm=0
      do l=0,lmax
      do m=-l,l
        lm=lm+1
        lmm=lm
        lmp=lm+lmmax
        alms(lmm)=(0.d0,0.d0)
        alms(lmp)=(0.d0,0.d0)
        do kmy=1,kmymax
        do kmyp=1,kmymax
          if(ind1(kmy).eq.lm.and.ind1(kmyp).eq.lm) 
     >    alms(lmm)=alms(lmm)+u1(kmy)*akmy(kmy,kmyp)*u1(kmyp)
          if(ind2(kmy).eq.lm.and.ind2(kmyp).eq.lm) 
     >    alms(lmp)=alms(lmp)+u2(kmy)*akmy(kmy,kmyp)*u2(kmyp)
        end do
        end do
      end do
      end do
c
      return
      end
