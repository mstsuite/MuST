c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getclm(lmax,clm)
c     ================================================================
c
      implicit   none
c
      integer    lmax
      integer    l
      integer    m
      integer    m2
      integer    i
c
      real*8     clm((lmax+1)*(lmax+2)/2)
      real*8     fpi
      real*8     fnpi
      real*8     xfac
      real*8     one
      real*8     four
c
      parameter (one=1.0d0)
      parameter (four=4.0d0)
c
c     *****************************************************************
c     Coeficients for complex spherical harmonics......................
c     Calclates all the c(l,m)'s up to lmax............................
c     
c                  [ (2*l+1)*(l-|m|)!]
c     c(l,m)=  sqrt[-----------------]
c                  [   4*pi*(l+|m|)! ]
c
c     *****************************************************************
c
c     =================================================================
      if(lmax.lt.0 ) then
         write(6,'('' GETCLM:: bad arguments: lmax='',i5)') lmax
         stop 'getclm'
      endif
c     =================================================================
c     begin calculation of c(l,m)'s....................................
c     =================================================================
      fpi=four*fnpi()
c     =================================================================
c     special case lmax=0..............................................
c     =================================================================
      clm(1)=sqrt(one/fpi)
      do l=1,lmax
         xfac=sqrt((2*l+1)/fpi)
         do m=0,l
            clm(l*(l+1)/2+m+1)=one
            m2=2*m
            do i=1,m2
               clm(l*(l+1)/2+m+1)=(l-m+i)*clm(l*(l+1)/2+m+1)
            enddo
c           ===========================================================
c           This version is consisten with (-1)**m being in Plm's......
c           See plglmax.f in this directory............................
c           ===========================================================
            clm(l*(l+1)/2+m+1)=xfac*sqrt(one/clm(l*(l+1)/2+m+1))
         enddo
      enddo
c
      return
      end
