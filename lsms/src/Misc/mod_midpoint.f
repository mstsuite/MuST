      subroutine mod_midpoint(y0,dy0,n,x0,htot,steps,y,
     &                        rn,ecomp,vr,eunit,b,allp1,nexp)
      implicit none

      integer n,steps,nexp
      real*8 y0(n),dy0(n),x0,htot,y(n)
      real*8 rn,vr(*),eunit,b,allp1
      complex*16 ecomp

      real*8 h,h2,y_i(n),y_j(n),dy(n),x,tmp
      integer i,j

      h=htot/steps
      h2=2.0*h

      x=x0
      do i=1,n
        y_j(i)=y0(i)
        y_i(i)=y0(i)+h*dy0(i)
      end do
      do j=2,steps
      x=x+h
      call dfv_m(x,y_i,dy,n,nexp,rn,ecomp,vr,eunit,b,allp1)
      do i=1,n
        tmp=y_j(i)+h2*dy(i)
        y_j(i)=y_i(i)
        y_i(i)=tmp
      end do
      end do
      x=x+h
      call dfv_m(x,y_i,dy,n,nexp,rn,ecomp,vr,eunit,b,allp1)
      do i=1,n
        y(i)=0.5d0*(y_i(i)+y_j(i)+h*dy(i))
      end do
      end



