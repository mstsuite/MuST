!
!
!
      subroutine interp(r,f,nr,rs,ps,dps)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!----- this routine interpolates functn f on grid r
!      returns value ps and its derivative dps 
!
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      real*8 r(nr),f(nr),c(4)
!
      ip1=2
      do i=ip1,nr-1
      if (rs.gt.r(i)) ip1=i+1
      enddo
      i=ip1-1
      call fit(r,f,nr,i,c)
!
      h1=rs-r(i)
      h2=r(ip1)-rs
      c2=1.d0+c(4)*h1*h2
      ps=((c(3)-c(2))*h1+c(1))/c2
      dps=c(2)+(c(3)-c(2)+ps*c(4)*(h1-h2))/c2
      ps=ps+f(i)-c(1)+c(2)*h1
      
      return
      end
