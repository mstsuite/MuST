c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getkncut(eta,a1,a2,a3,kncut,nm1,nm2,nm3)
c     ================================================================
c
      implicit   none
c
      integer    nm1
      integer    nm2
      integer    nm3
      integer    i
      integer    j
      integer    k
c
      real*8     eta
      real*8     a1(3)
      real*8     a2(3)
      real*8     a3(3)
      real*8     r(3)
      real*8     rm2
      real*8     term
      real*8     kncut
      real*8     fac
      real*8     epsi
c
      parameter (epsi=1.0d-14)
c
      fac=eta*eta/4.0d0
c
c     ================================================================
c     calculate nm1,nm2,nm3...........................................
c     ================================================================
      r(1)=a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3)
      term=1.0d0
      nm1=0
      do while(term.gt.epsi) 
         nm1=nm1+1
         rm2=nm1*nm1*r(1)
         term=exp(-fac*rm2)/rm2
      enddo
      r(2)=a2(1)*a2(1)+a2(2)*a2(2)+a2(3)*a2(3)
      term=1.0d0
      nm2=0
      do while(term.gt.epsi) 
         nm2=nm2+1
         rm2=nm2*nm2*r(2)
         term=exp(-fac*rm2)/rm2
      enddo
      r(3)=a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3)
      term=1.0d0
      nm3=0
      do while(term.gt.epsi) 
         nm3=nm3+1
         rm2=nm3*nm3*r(3)
         term=exp(-fac*rm2)/rm2
      enddo
c
c     ================================================================
c     calculate kncut.................................................
c     ================================================================
      kncut=sqrt(r(1))*nm1
      do i=-1,1
	 r(1)=i*a1(1)*nm1
	 r(2)=i*a1(2)*nm1
	 r(3)=i*a1(3)*nm1
         do j=-1,1
	    r(1)=r(1)+j*a2(1)*nm2
	    r(2)=r(2)+j*a2(2)*nm2
	    r(3)=r(3)+j*a2(3)*nm2
            do k=-1,1
	       r(1)=r(1)+k*a3(1)*nm3
	       r(2)=r(2)+k*a3(2)*nm3
	       r(3)=r(3)+k*a3(3)*nm3
	       rm2=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
	       kncut=max(kncut,sqrt(rm2))
            enddo
         enddo
      enddo
c
      return
      end
