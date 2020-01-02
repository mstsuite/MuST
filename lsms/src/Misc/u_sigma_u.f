c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u_sigma_u(u,ud,usu_x,usu_y,usu_z)
c     ================================================================
c
      implicit   none
c
      complex*16 u(4)
      complex*16 ud(4)
      complex*16 usu_x(4)
      complex*16 usu_y(4)
      complex*16 usu_z(4)
      complex*16 sqrtm1
c
      parameter (sqrtm1=(0.0d0,1.0d0))
c
c     ****************************************************************
c     calculate:
c                     +                   +                   +
c        u * sigma * u ,     u * sigma * u ,     u * sigma * u
c                 x                   y                   z
c     ****************************************************************
c
      usu_x(1)=u(3)*ud(1)+u(1)*ud(2)
      usu_x(2)=u(4)*ud(1)+u(2)*ud(2)
      usu_x(3)=u(3)*ud(3)+u(1)*ud(4)
      usu_x(4)=u(4)*ud(3)+u(2)*ud(4)
c
      usu_y(1)=(u(3)*ud(1)-u(1)*ud(2))*sqrtm1
      usu_y(2)=(u(4)*ud(1)-u(2)*ud(2))*sqrtm1
      usu_y(3)=(u(3)*ud(3)-u(1)*ud(4))*sqrtm1
      usu_y(4)=(u(4)*ud(3)-u(2)*ud(4))*sqrtm1
c
      usu_z(1)=u(1)*ud(1)-u(3)*ud(2)
      usu_z(2)=u(2)*ud(1)-u(4)*ud(2)
      usu_z(3)=u(1)*ud(3)-u(3)*ud(4)
      usu_z(4)=u(2)*ud(3)-u(4)*ud(4)
c
      return
      end
