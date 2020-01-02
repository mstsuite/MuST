      subroutine rzextr (iest,xest,yest,yz,dy,nv,nuse,x,d)
!      implicit none
      implicit real*8 (a-h,o-z)
      integer iest,nv,nuse,imax,nmax,ncol
      integer j,k
      real*8 xest,yest,yz,dy
      parameter (imax=11,nmax=16,ncol=7)
      dimension x(imax),yest(nv),yz(nv),dy(nv),d(nmax,ncol),fx(ncol)
! saves: x,d
!      save 
      x(iest)=xest
      if(iest.eq.1)then
        do j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
        enddo
      else
        m1=min(iest,nuse)
        s1xest=1.0d0/xest
        do k=1,m1-1
          fx(k+1)=x(iest-k)*s1xest
        enddo
        do j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do k=2,m1
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.0d0)then
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            else
              ddy=v
            endif
            if(k.ne.m1) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
          enddo
          dy(j)=ddy
          yz(j)=yy
        enddo
      endif
      return
      end
