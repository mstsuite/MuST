      subroutine calculate_asa_ro3(n_spin_pola,rho1,rho2,
     &           i_vdif,r_mesh,jmt,rmt,iprpts,
     &           ro3,dz)
      implicit none

      integer jmt,i_vdif,iprpts,n_spin_pola
      integer i
      real*8 rmt
      real*8 rho1(iprpts),rho2(iprpts)
      real*8 rhojmt,rhojmt1,rhojmt2,drhot
      real*8 r_mesh(iprpts),rtmp(0:iprpts)
      real*8 ro3,dz

      rtmp(0)=0.d0
      do i=1,jmt+2
        rtmp(i)=sqrt(r_mesh(i))
      enddo

      rhojmt2=0.0
      call interp(rtmp(1),rho1(1),jmt,
     >   sqrt(rmt),rhojmt1,drhot,.false.)
      if(n_spin_pola.eq.2) then
        call interp(rtmp(1),rho2(1),jmt,
     >     sqrt(rmt),rhojmt2,drhot,.false.)
      end if
      rhojmt=rhojmt1+(n_spin_pola-1)*rhojmt2

      ro3 = 1.0d10
      dz = 0.0

      if(rhojmt.gt.1.d-10)
     >   ro3=(3.0*rmt*rmt/rhojmt)**(1.0/3.0)
      if(i_vdif.ne.0 .and. n_spin_pola.eq.2)
     >   dz =(rho1(jmt)-rho2(jmt))/rhojmt

      end subroutine
