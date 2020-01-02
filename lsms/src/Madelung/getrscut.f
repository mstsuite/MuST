c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getrscut(eta,a1,a2,a3,rscut,nm1,nm2,nm3)
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
      real*8     rm
      real*8     term
      real*8     rscut
      real*8     epsi
      real*8     erfc
c
      parameter (epsi=1.0d-14)
c
c     ================================================================
c     calculate nm1,nm2,nm3...........................................
c     ================================================================
      r(1)=sqrt(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
      term=1.0d0
      nm1=0
      do while(term.gt.epsi) 
         nm1=nm1+1
         rm=nm1*r(1)
         term=erfc(rm/eta)/rm
      enddo
      r(2)=sqrt(a2(1)*a2(1)+a2(2)*a2(2)+a2(3)*a2(3))
      term=1.0d0
      nm2=0
      do while(term.gt.epsi) 
         nm2=nm2+1
         rm=nm2*r(2)
         term=erfc(rm/eta)/rm
      enddo
      r(3)=sqrt(a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3))
      term=1.0d0
      nm3=0
      do while(term.gt.epsi) 
         nm3=nm3+1
         rm=nm3*r(3)
         term=erfc(rm/eta)/rm
      enddo
c
c     ================================================================
c     calculate rscut.................................................
c     ================================================================
      rscut=r(1)*nm1
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
	       rm=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
	       rscut=max(rscut,rm)
            enddo
         enddo
      enddo
c
      return
      end
