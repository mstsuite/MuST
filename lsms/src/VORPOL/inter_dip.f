      subroutine inter_dip(lmax,ncrit,ngaussr,omegint,
     >                     grwylm,gwwylm,wylm,dipole,vdipole)
c     =================================================================
c
      implicit   none
c
      integer    lmax,ncrit
      integer    ir
      integer    m1
      integer    iwylm
      integer    ng
      integer    ngaussr
c
      real*8     rgauss,omegint
      real*8     gwwylm(ngaussr,ncrit-1)
      real*8     grwylm(ngaussr,ncrit-1)
c
      complex*16 wylm((lmax+1)*(lmax+2)/2,ngaussr,ncrit-1)
      complex*16 dipole(-1:1)
      complex*16 vdipole(-1:1)
      complex*16 steplm
      complex*16 czero
c
      parameter (czero=(0.0d0,0.0d0))
      real*8 pi,fourpi
      parameter (pi=3.14159265358979d0)
      parameter (fourpi=4.d0*pi)
c
c     =================================================================
c     -----------------------------------------------------------------
      do m1=-1,1
	dipole(m1)=czero
	vdipole(m1)=czero
      enddo
      do ir=1,ncrit-1
c        ==============================================================
c        loop over the gaussian point in each region...................
c        ==============================================================
        do ng=1,ngaussr
          rgauss=grwylm(ng,ir)
	  if(rgauss.gt.0.d0) then
	    do m1=-1,1
              iwylm= 2 + abs(m1)
              if(m1.ge.0) then
                steplm=wylm(iwylm,ng,ir)
              else
                steplm=(1-2*mod(abs(m1),2))*conjg(wylm(iwylm,ng,ir))
              endif
              dipole(m1)=dipole(m1)+gwwylm(ng,ir)*steplm*rgauss
              vdipole(m1)=vdipole(m1)+
     &              gwwylm(ng,ir)*steplm/(3.d0*rgauss*rgauss)
            enddo
	  endif  !  rgauss.gt.zero
        enddo
      enddo
      do m1=-1,1
        dipole(m1)=dipole(m1)*fourpi
        vdipole(m1)=vdipole(m1)*fourpi
      enddo
c
      return
      end
