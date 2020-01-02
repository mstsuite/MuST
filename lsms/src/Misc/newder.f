c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine newder(f,dy,r,nr)
c     ================================================================
c
      implicit real*8 (a-h,o-z)
c
      integer    i
c
      real*8 r(nr),c(4),f(nr),dy(nr)
c
      ip1=1
      do i=1,nr-1
	ip1=ip1+1

        call fit(r,f,nr,i,c)
	
      dr=r(ip1)-r(i)
      c2=c(4)*dr
      if(i.gt.1) then
c Use the previous fit and the current fit to form a five-point interpolation
c   f5(r)=0.5*(f4(r)+f4p(r))+(f4(r)-f4p(r))*(r-0.5*(rm+rp))/(rp-rm)
c So the derivative is
c   df5(r)=0.5*(df4(r)+df4p(r))+(df4(r)-df4p(r))*(r-0.5*(rm+rp))/(rp-rm)
c          +(f4(r)-f4p(r))/(rp-rm)
c   but f4(r)=f4p(r)=f(i)
	dyold=dy(i)
	dynew=c(3)-c(1)*c2
        dy(i)=0.5d0*(dy(i)+dynew)
	if(i.gt.2.and.i.lt.nr-1) then
	  dr=r(i+2)-r(i-2)
	  rj=0.5d0*(r(i+2)+r(i-2))
	  dy(i)=dy(i)+(dynew-dyold)*(r(i)-rj)/dr
	endif
      else
        dy(i)=c(3)-c(1)*c2
      endif
      dy(i+1)=c(3)+(f(ip1)-f(i)+c(1)-c(2)*dr)*c2
      enddo
      return
      end
