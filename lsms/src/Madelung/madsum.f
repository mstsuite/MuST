c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine madsum(lastkn,lastrs,ibegin,aij,
     >                  rslatmd,knlat_x,knlat_y,knlat_z,knlatsq,
     >                  eta,tau,madmat,
     >                  pi4,iprint,istop)
c     ================================================================
c
c     ****************************************************************
c     performs ewald summation for Madelung constant.
c     requires 'lastkn,lastrs,eta,knlat,aij'
c     returns 'madelung sum': madmat
c     ****************************************************************
      implicit   none
c
      character  istop*32
      character  sname*32
      parameter (sname='madsum')
c
      integer    lastkn
      integer    lastrs
      integer    ibegin
      integer    iprint
      integer    i
c
      real*8     knlat_x(lastkn)
      real*8     knlat_y(lastkn)
      real*8     knlat_z(lastkn)
      real*8     knlatsq(lastkn)
      real*8     rslatmd(lastrs)
      real*8     aij(3)
      real*8     eta
      real*8     tau
      real*8     madmat
      real*8     pi4
      real*8     fac
      real*8     zero
      real*8     four
      real*8     erfc
      real*8     erf
c
      parameter (zero=0.0d0)
      parameter (four=4.0d0)
c
c     ================================================================
c                                     2   2       -> ->
c              4*pi          1    -eta *Kq /4 - i*Kq*aij
c     term1 =  ---- * sum  ----- e
c              tau    q<>0    2
c                           Kq
c
c     note sum starts at 2 since kn=0.0 of 1/(rs-aij) is
c     canceled by kn=0.0 of the 1/rs sum.
c     ================================================================
      madmat = zero
      fac=eta*eta/four
      do i = 2,lastkn
         madmat=madmat+exp(-fac*knlatsq(i))/knlatsq(i)
     >                *cos( knlat_x(i)*aij(1)+knlat_y(i)*aij(2)
     >                     +knlat_z(i)*aij(3) )
c if(i.le.5)write(6,'(''i,a,k,f,m'',i5,10f10.4)')
c     >   i,aij(1),aij(2),aij(3),knlat_x(i),knlat_y(i)
c     >   ,knlat_z(i),exp(-fac*knlatsq(i))/knlatsq(i)*pi4/tau
c     >        ,madmat*pi4/tau*6./3.5449077d0
      enddo
      madmat = pi4/tau*madmat
c     ================================================================
c
c                           ->   ->       
c                  1 - erf(|Rn + aij|/eta) 
c     term2 = sum -------------------------
c              n         ->   ->        
c                       |Rn + aij|
c
c     note for calculation of aij=0.0 term ibegin=2.
c     ================================================================
      do i=ibegin,lastrs
c     do i=ibegin,1
	 erf=erfc(rslatmd(i)/eta)/rslatmd(i)
         madmat = madmat + erf
c     if(i.lt.20.and.iprint.ge.0)write(6,'('' i,erf'',i5,3f15.7)')i,
c    >rslatmd(i),eta,erf
      enddo
c
      return
      end
