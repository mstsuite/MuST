!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cnewder(f,dy,r,nr)
!     ================================================================
!
      implicit real*8 (a-h,o-z)
!
      integer    i, ip1, nr
!
      real*8 r(nr), c2_r, c2_i, c_r(4), c_i(4), f_r(nr), f_i(nr)
!
      complex*16 dyold, dynew, f(nr), dy(nr)
!
      do i = 1,nr
         f_r(i) = real(f(i)) 
         f_i(i) = aimag(f(i))
      enddo
      ip1=1
      do i=1,nr-1
         ip1=ip1+1

         call fit(r,f_r,nr,i,c_r)
         call fit(r,f_i,nr,i,c_i)

         dr=r(ip1)-r(i)
         c2_r=c_r(4)*dr
         c2_i=c_i(4)*dr
         if (i.gt.1) then
! Use the previous fit and the current fit to form a five-point interpolation
!   f5(r)=0.5*(f4(r)+f4p(r))+(f4(r)-f4p(r))*(r-0.5*(rm+rp))/(rp-rm)
! So the derivative is
!   df5(r)=0.5*(df4(r)+df4p(r))+(df4(r)-df4p(r))*(r-0.5*(rm+rp))/(rp-rm)
!          +(f4(r)-f4p(r))/(rp-rm)
!   but f4(r)=f4p(r)=f(i)
            dyold=dy(i)
            x_r = c_r(3)-c_r(1)*c2_r
            x_i = c_i(3)-c_i(1)*c2_i
            dynew=cmplx(x_r,x_i)
            dy(i)=0.5d0*(dy(i)+dynew)
            if (i.gt.2.and.i.lt.nr-1) then
               dr=r(i+2)-r(i-2)
               rj=0.5d0*(r(i+2)+r(i-2))
               dy(i)=dy(i)+(dynew-dyold)*(r(i)-rj)/dr
            endif
         else
            x_r = c_r(3)-c_r(1)*c2_r
            x_i = c_i(3)-c_i(1)*c2_i
            dy(i)=cmplx(x_r,x_i)
         endif
         x_r = c_r(3)+(f_r(ip1)-f_r(i)+c_r(1)-c_r(2)*dr)*c2_r
         x_i = c_i(3)+(f_i(ip1)-f_i(i)+c_i(1)-c_i(2)*dr)*c2_i
         dy(i+1)=cmplx(x_r,x_i)
      enddo
      return
      end
