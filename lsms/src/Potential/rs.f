      function rs(chd,radius,rsmax)
c
c  purpose:  convert a charge density value to the equivalent r_s 
c            where   4*pi*(r_s)**3*chd/3 = 1, 
c            i.e. r_s is the radius of a sphere containing exactly one
c                 electron in a homogeneous electron gas of density chd. 
c
c            charge density may be of the form n(r) (radius = 0.d0)
c            or 4*pi*r**2*n(r), in which case radius is the r-value.
c
c  method:   check that chd is not zero or -ve; if it is, rs is set to
c            rsmax (rsmax = 1.d3 => charge density around 2.4d-10 e/(au^3))
c
c   calls:  system function datan
c
c written:  philip sterne, sterne1@llnl.gov
c date:     April 21, 1998
c changes:  none
c comments:  
c
      implicit none
c
c  input variables
c
      real*8 chd        ! charge density - either n(r) or 4*pi*r**2*n(r)
      real*8 radius     ! radius r if input chd is 4*pi*r**2*n(r)
      real*8 rsmax      ! maximum allowable r_s value (very low cd limit)
c
c  output variables
c  
      real*8 rs         ! r_s value for this charge density
c
c  internal variables and functions
c
      real*8 pi         ! pi
      real*8 chdmin     ! minimum charge density - corresponds to rsmax
c     
c  set rsmax and chdmin values for vanishingly small charge densities
c
      rs = rsmax
      pi = 4.d0*datan(1.d0)
      if(radius.eq.0.d0) then
         chdmin = 3.d0/(4.d0*pi*rsmax**3)
         if(chd.gt.chdmin) rs = (3.d0/(4.d0*pi*chd))**(1.d0/3.d0)
      else
         chdmin = 3.d0*radius**2/(rsmax**3)
         if(chd.gt.chdmin) rs = (3.d0*radius**2/chd)**(1.d0/3.d0)
      endif
c
      return
      end

