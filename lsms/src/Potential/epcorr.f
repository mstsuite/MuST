      function epcorr(chd,gradc,radius,ivepc,alpgga,rspmax)
c
c   purpose: compute the positron correlation potential for a given
c            electron charge density
c
c   method:  
c            compute potential based on parameters and functional form.
c            if alpgga is not zero, make the gga correction based on
c            the gradient of the charge density.
c
c   calls:  fns : rs, intrinsic functions datan, exp, sqrt, log
c           
c written:  philip sterne, sterne1@llnl.gov
c date:     May 04, 1998
c changes:
c     7/13/98 : added rspmax to control maximum rs value in computing 
c               ep-corr potentials.                  sterne1@llnl.gov
c comments: 
c  ivepc    positron correlation   
c                               
c    0     Sterne-Kaiser        PRB 43, 13892 (1991)    
c    1     Boronski-Nieminen    PRB 34,  3820 (1986)
c    2     IPM                  No interaction - Independent Particle Model   
c                               
c    Sterne-Kaiser and Boronski-Nieminen are both fit to rpa, 
c    arponnen and pajanne (ann phys 121 343), and low density limits.
c                               
      implicit none
c
c     input variables
c      
      real*8 chd        ! charge density - either n(r) or 4*pi*r**2*n(r)
      real*8 gradc      ! magnitude of the gradient of the charge density
                        !  - either |del(n(r))| or 4*pi*r**2*|del(n(r))|
      real*8 radius     ! if zero, chd is just n(r), else radius r for chd
      integer ivepc     ! choice of electron-positron correlation potential
      real*8 alpgga     ! alpha for gga correction - zero means no gga
      real*8  rspmax    ! maximum rs in positron-electron correlation.
                        ! for larger rs values, epcorr(rspmax) is used
c
c     output variables
c
      real*8 epcorr     ! calculated electron-positron correlation
c
c     internal variables
c
      real*8  skc(5)    ! parameters for Sterne-Kaiser fit
      real*8 bnc(4,4)   ! parameters for Boronski-Nieminen fit
      real*8 bnrs(4)    ! rs-cutoffs for Boronski-Nieminen fit
      real*8 rs         ! function - calculate rs from charge density
      real*8 rs1        ! rs value corresponding to input chd
      real*8 rsd        ! shifted rs value - used in functional
      real*8 fwvec      ! fermi wavevector
      real*8 qtfsq      ! (Thomas-Fermi screening length)**2
      real*8 epsgga     ! contribution to gga correlation function
      real*8 expgga     ! contribution to gga correlation function
      real*8 pi         ! well-known approximation to 22/7
c
      data skc /-1.56d0, 0.7207d0, 4.092d0, 51.96d0, 0.1324d0/
c
      data bnc / -1.56d0,     0.051d0,   -0.081d0,   1.14d0,
     *           -0.92305d0, -0.05459d0,  0.d0,      0.d0,
     *          -13.15111d0,  2.5d0,      2.8655d0, -0.6298d0, 
     *       -10250.5786d0,  44.50466d0, -0.524d0,   0.d0/
c
      data bnrs/0.302d0, 0.56d0, 8.0d0, 200.d0/
c
c     convert from charge density to r_s
c
      rs1 = rs(chd,radius,rspmax)
c ---------------------------------------------------------------------------
c       electron-positron correlation potential functional
c ---------------------------------------------------------------------------
c
      if (ivepc.eq.0) then
c
c          sterne-kaiser fit
c
         if (rs1.gt.bnrs(4)) then
            epcorr=bnc(3,4)
         else 
            epcorr=skc(1)/sqrt(atan(rs1))+skc(2)+
     *             skc(5)*exp(-((rs1-skc(3))**2)/skc(4))
         endif
c     
      else if (ivepc.eq.1) then   
c
c          boronski-neimenen fit
c
         if (rs1.gt.bnrs(4)) then
            epcorr=bnc(3,4)
         else if (rs1.gt.bnrs(3)) then
            epcorr=bnc(1,4)/(rs1**6)+bnc(2,4)/(rs1**3)+bnc(3,4)
         else if (rs1.gt.bnrs(2)) then
            rsd=rs1+bnc(2,3)
            epcorr=bnc(1,3)/rsd/rsd+bnc(3,3)/rsd+bnc(4,3)
         else if (rs1.gt.bnrs(1)) then
            epcorr=bnc(1,2)+bnc(2,2)/rs1/rs1
         else 
            epcorr=bnc(1,1)/sqrt(rs1)+
     *           (bnc(2,1)*log(rs1)+bnc(3,1))*log(rs1)+bnc(4,1)
         endif
c     
      else if(ivepc.eq.2) then 
c     
c     no electron-positron correlation - IPM
c
         epcorr=0.d0
c     
      else
         write(*,*) 'error in epcorr: ivepc = ',ivepc,
     *                 '     unknown option - stopping'
         stop 'ivepc in epcorr'
      endif                                                               
c     
c   gga correction 
c
      if(alpgga.gt.0.d0) then
c
c   find local fermi wavevector fwvec and use it to calculate
c   the square of the thomas-fermi screening length, qtfsq
c
         pi = 4.d0*datan(1.d0)
         fwvec = ((9.d0*pi/4.d0)**(1.d0/3.d0))/rs1
         qtfsq = 4.d0*fwvec/pi
c
c   compute the lda gradient correction term epsgga
c   and the exponential factor expgga = exp(-alpgga*epsgga)
c
         if(rs1.lt.rspmax) then
            epsgga = (gradc/chd)**2/qtfsq
            expgga = dexp(-alpgga*epsgga)
         else
            epsgga = 0.d0
            expgga = 1.d0
         endif
c     
c   correct epcorr with gga term
c
         epcorr = epcorr*expgga**(1.d0/3.d0)
      endif
c
      return 
      end

