      subroutine fit(r,f,nr,i,c)
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      real*8 r(nr),f(nr),c(4)
! does a four-point interpolation from i-1 to i+2
! returns c array 4 of real*8 which are coeficients for
! f(i)-c1+c2*(r-r(i))+((c3-c2)*(r-r(i))+c1)/(1+c4*(r-r(i))*(r(i+1)-r))
!
! It's obvious that c3=(f(i+1)-f(i))/(r(i+1)-r(i))
! c2 is adjustable and is chosen so to avoid possible singularities
! in eqns solving for c1 and c4.
!
!   c2=((f(i+2)-f(i+1))/(r(i+2)-r(i+1))-(f(i)-f(i-1))/(r(i)-r(i-1)))
! and the sign is chosen so that it is opposite of
! (f(i+1)-f(i-1))/(r(i+1)-r(i-1))
!

      if(i.lt.1) then
        write(6,'(''FIT: i<1'')')
        stop 'fit'
      endif
      ip1=i+1
      if(ip1.gt.nr) then
        if(r(1).gt.r(nr)) then   ! wrap around
          ip1=1
        else
        write(6,'(''FIT: i>n-1'')')
        stop 'fit'
        endif
      endif
      dr=r(ip1)-r(i)
      if(dr.eq.0.d0) then
        write(6,'(''FIT: degenerate grid'')')
        stop 'fit'
      endif
      drv=f(ip1)-f(i)
      c(3)=drv/dr

! two point interpolation
      if(nr.eq.2) then
        c(1)=f(i)
        c(2)=0.d0
        c(4)=0.d0
        return
      endif

! fit the function to the two outside points
      if(i.gt.1) then
        i1=i-1
      else
        i1=min(4,nr)
      endif
      if(ip1.lt.nr) then
        i2=ip1+1
      else
        i2=max(nr-3,1)
      endif
! If nr=3 then i1=i2
      if(i1.ne.i2) then
        c(2)=(f(i2)-f(ip1))/(r(i2)-r(ip1))-(f(i1)-f(i))/(r(i1)-r(i))
        if(c(2)*(f(i2)-f(i1))*(r(i2)-r(i1)).gt.0.d0) c(2)=-c(2)
! solve the 2x2 equation for c(1) and c(4)
        h1=r(i1)-r(i)
        c0=f(i1)-f(i)
        c1=c0-c(3)*h1
        c2=c0-c(2)*h1
        h1=h1*(r(ip1)-r(i1))
        eqn12=-h1*c2
        h2=r(i2)-r(i)
        c0=f(i2)-f(i)
        c3=c0-c(3)*h2
        c4=c0-c(2)*h2
        h2=h2*(r(ip1)-r(i2))
        eqn22=-h2*c4
        gj=(c4-c2)*h1*h2
        if(gj.eq.0.d0) then
          c(4)=0.d0
        else
          c(1)=(c1*eqn22-c3*eqn12)/gj
          c(4)=(c1*h2-c3*h1)/gj
        endif
! If the points are on a line or nearly a line then use linear
! interpolation
! check this by checking to see if the denominator 1+c4*(r-r(i))*(r(i+1)-r)
! can be zero or negative. The minimum is 1-c4*dr*dr/4
        gj=c(4)*dr*dr
        if(gj.gt.-4.d0.and.abs(gj).gt.1.d-14.and.gj.lt.1.d+14) then
          c(1)=c(1)/c(4)
        else
          c(1)=f(i)
          c(4)=0.d0
        endif
      else  ! i1=i2
! set A4=0.5*(f(i)+f(i+1)) and do a 3 point fit
        fj=0.5d0*(f(i)+f(ip1))
        c(1)=f(i)-fj
        c(2)=0.d0
        rj=r(i)+0.5d0*dr
        if(abs(r(i1)-rj).gt.abs(r(i2)-rj)) i1=i2
        h1=r(i1)-r(i)
        c1=f(i1)-f(i)-c(3)*h1
        h1=h1*(r(ip1)-r(i1))
        if(f(i1).ne.fj) then
          c(4)=-c1/((f(i1)-fj)*h1)
          gj=c(4)*dr*dr
          if(gj.le.-4.d0.or.abs(gj).le.1.d-14.or.gj.ge.1.d+14) then
            c(1)=f(i)
            c(4)=0.d0
          endif
        else
          c(1)=f(i)
          c(4)=0.d0
        endif
      endif  ! f(i1)=f(i2)

      return
      end
