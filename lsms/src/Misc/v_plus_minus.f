c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine v_plus_minus(fac,vr1,vr2,vr,br,jmt,iprpts)
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit   none
c
!     include    'atom_param.h'
      integer iprpts
c
c     ****************************************************************
c     calculates either v=v_up+v_down & b=v_up-v_down (fac=1)
c     or                v_up=(v+b)/2 & v_down=(v-b)/2 (fac=1/2)
c     ****************************************************************
c
c
      integer    jmt
      integer    j
c
      real*8     fac
      real*8     vr1(iprpts)
      real*8     vr2(iprpts)
      real*8     vr(iprpts)
      real*8     br(iprpts)
c
      do j=1,jmt
         vr(j)=fac*(vr1(j)+vr2(j))
         br(j)=fac*(vr1(j)-vr2(j))
      end do
c
      return
      end

