      function gaunt(l1,m1,l2,m2,l3,m3,lmax,ngauss,
     >              wg,plmg,clm)
c     ****************************************************************
c
c               4pi  _  m1 _     m2* _     m3 _
c       gaunt = int do Y  (o) * Y   (o) * Y  (o)
c                0      l1       l2        l3   
c
c                L3
c             = C
c                L1,L2
c
c     ****************************************************************
c
c  Inputs: l1,m1,l2,m2,l3,m3   integer scalars, (l,m) indices of the gaunt
c                              number
c          lmax     integer scalar, max l values allowed, lmax >= max(l1,l2,l3)
c                   lmax also determines the array dimensions of plmg and clm
c          ngauss   integer scaler, Number of gaussian points for integration
c          wg       real*8 array of (ngauss), weights for mesh points
c                   wg are the Gauss-Legendre integration weights
c          plmg     real*8 array of ((lmax+1)*(lmax+2)/2,ngauss), where
c                   l=max(l1,l2,l3), the values of the associate Legendre
c                   function at the mesh points
c                   plmg is output of plglmax
c          clm      real*8 array of ((lmax+1)*(lmax+2)/2), prefactors for the
c                   spherical harmonics
c                   clm is output of getclm
c
      implicit  none
c
      integer   l1
      integer   l2
      integer   l3
      integer   m1
      integer   m2
      integer   m3
      integer   lmax,ngauss
      integer   ng
      integer   ifac1
      integer   jl1
      integer   ifac2
      integer   jl2
      integer   ifac3
      integer   jl3
c
      real*8    wg(ngauss)
      real*8    plmg((lmax+1)*(lmax+2)/2,ngauss)
      real*8    clm((lmax+1)*(lmax+2)/2)
      real*8    gaunt,gaunt_val
      real*8    zero
      real*8    four
c
      parameter (zero=0.d0)
      parameter (four=4.d0)
      real*8    pi,fourpi
      parameter (pi=3.141592653589793d0)
      parameter (fourpi=four*pi)
c
      if(l1.gt.lmax .or. l2.gt.lmax .or. l3.gt.lmax ) then
         write(6,'(''gaunt:: bad parameters: l1,l2,l3,lmax'',4i5)')
     >                                       l1,l2,l3,lmax
         stop'gaunt'
      endif
      gaunt_val=zero
      if(mod(l1+l2+l3,2).ne.0) then
         return
      else
         if(m1-m2+m3 .ne. 0) then
            return
         else
            if(l1+l2.lt.l3 .or. l2+l3.lt.l1 .or. l3+l1.lt.l2) then
               return
            else
c              -------------------------------------------------------
               call defac(l1, m1,ifac1,jl1)
               call defac(l2,-m2,ifac2,jl2)
               call defac(l3, m3,ifac3,jl3)
c              -------------------------------------------------------
               do ng=1,ngauss
                  gaunt_val=gaunt_val+wg(ng)*plmg(jl1,ng)*
     >                 plmg(jl2,ng)*
     >                 plmg(jl3,ng)
               enddo
	       if(abs(gaunt_val).lt.1.d-14) then
		 gaunt_val=0.d0
	       else
               gaunt_val=fourpi*gaunt_val*ifac1*clm(jl1)*
     &              (1-2*mod(abs(m2),2))*ifac2*clm(jl2)*ifac3*clm(jl3)
	       endif
            endif
         endif
      endif

      gaunt=gaunt_val
      return
      end
c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defac(l,m,ifac,jl)
c     ================================================================
c
      implicit   none
c
      integer    l,m
      integer    jl
      integer    ifac
c
c     ================================================================
      if (m.ge.0) then
         ifac= 1
         jl = (l+1)*(l+2)/2-l+m
      else
         ifac= (1-2*mod(abs(m),2))
         jl = (l+1)*(l+2)/2-l-m
      end if
c
      return
      end
