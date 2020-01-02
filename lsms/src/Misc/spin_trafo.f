      subroutine spin_trafo(evec,u,ud)
c
c     ================================================================
c     set up the spin space stransformation matrix....................
c
c     u(1) = cos(theta/2)  and u(3) = sin(theta/2)*exp(-i*phi)
c     u(4) = u(1)          and u(2) = - conjg(u(3))
c
c           +
c     ud = u
c     --   -
c
c     Matrix in global frame: G;  Matrix in local frame: L
c                             -                          -
c                  +
c     L = u * G * u
c     -   -   -   -
c     ================================================================
c
      implicit   none
c
      real*8     evec(3)
      real*8     tol
      real*8     half
      real*8     one
c
      complex*16 u(4)
      complex*16 ud(4)
c
      complex*16 czero
      complex*16 cone
      complex*16 sqrtm1
c
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (tol=1.0d-08)
      parameter (czero=(0.0d0,0.0d0))
      parameter (cone=(1.0d0,0.0d0))
      parameter (sqrtm1=(0.0d0,1.0d0))
c
c     ================================================================
c
      u(1) = sqrt(half*(one+evec(3)))
      if( abs(u(1)) .lt. tol ) then
         u(1) = czero
         u(2) = cone
         u(3) = cone
      else
         u(2) =-half*(evec(1)+sqrtm1*evec(2))/u(1)
         u(3) = half*(evec(1)-sqrtm1*evec(2))/u(1)
      endif
      u(4)=u(1)
      ud(1)=conjg(u(1))
      ud(2)=conjg(u(3))
      ud(3)=conjg(u(2))
      ud(4)=conjg(u(4))
c     ================================================================
      return
      end
