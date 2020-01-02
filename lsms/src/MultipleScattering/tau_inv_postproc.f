      subroutine tau_inv_postproc_nrel(kkrsz_ns,n_spin_cant,
     &                            wbig,delta,tmat,ipvt,tau00,
     &                            ubr,ubrd,
     &                            tau00_l)
      implicit none

      integer n_spin_cant
      integer kkrsz_ns
      integer ipvt(*)
      complex*16 wbig(*)
      complex*16 tmat(*)
      complex*16 tau00(*)
      complex*16 tau00_tmp(kkrsz_ns,kkrsz_ns)
      complex*16 tau00_l(*)
      complex*16 delta(*),ubr(*),ubrd(*)
      integer kkrsz

      integer info
      integer mtxsize

      complex*16 cmone,cone,czero
      parameter (cmone=(-1.0d0,0.0d0))
      parameter (cone=(1.0d0,0.0d0))
      parameter (czero=(0.d0,0.d0))

      kkrsz=kkrsz_ns/2
      mtxsize=kkrsz_ns*kkrsz_ns
!
!     FROM LSMS_1.9: GETTAU_CL
!
c     setup unit matrix...............................................
c     ----------------------------------------------------------------
      call cmtruni(wbig,kkrsz_ns)
c     ----------------------------------------------------------------
c     get 1-delta and put it in wbig
      call zaxpy(mtxsize,cmone,delta,1,wbig,1)
c     ================================================================
c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
c     ----------------------------------------------------------------
      call zgetrf(kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,info)
      call zcopy(kkrsz_ns*kkrsz_ns,tmat,1,tau00,1)
      call zgetrs('n',kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,tau00,
     &           kkrsz_ns,info)
!         write(*,*) delta(1), delta(kkrsz*kkrsz*2+kkrsz+1)
!
!    FROM LSMS_1.9: GETTAU
!
c     ----------------------------------------------------------------
c  Redefine tau00 to be tau00-t
c  delta is 1-t*tau00^{-1} and is calculated in gettaucl
c  and then rotated into the local frame
c     call scale_tau00(tau00_g,kkrsz,kkrsz,lofk,n_spin_cant,
c    &                 kappa_rmt)
      call zgemm('n','n',kkrsz_ns,kkrsz_ns,kkrsz_ns,cone,
     &           delta,kkrsz_ns,
     >           tau00,kkrsz_ns,czero,
     &           tau00_tmp,kkrsz_ns)
c     call inv_scale_tau00(tau00_tmp,kkrsz,kkrsz,lofk,n_spin_cant,
c    &                    kappa_rmt)
c
c     ================================================================
c     Rotate tau00 to local frame of reference
c     ================================================================
      if( n_spin_cant .eq. 2 ) then
!        Non relativistic
c        -------------------------------------------------------------
         call trgtol(kkrsz,kkrsz,ubr,ubrd,tau00_tmp,tau00_l)
      else
c        -------------------------------------------------------------
         call zcopy(kkrsz_ns*kkrsz_ns,tau00_tmp,1,tau00_l,1)
c        -------------------------------------------------------------
      endif


!      write(*,*) 'tau00_l(1,1)=',tau00_l(1)

      end subroutine

      subroutine tau_inv_postproc_rel(kkrsz_ns,
     &                            wbig,delta,tmat,ipvt,tau00,
     &                            dmatp,dmat,
     &                            tau00_l)
      implicit none

      integer kkrsz_ns
      integer ipvt(*)
      complex*16 wbig(*)
      complex*16 tmat(*)
      complex*16 tau00(*)
      complex*16 tau00_l(*)
      complex*16 delta(*)
      complex*16 dmatp,dmat

      integer info
      integer mtxsize

      complex*16 cmone,cone,czero
      parameter (cmone=(-1.0d0,0.0d0))
      parameter (cone=(1.0d0,0.0d0))
      parameter (czero=(0.d0,0.d0))

      mtxsize=kkrsz_ns*kkrsz_ns
!
!     FROM LSMS_1.9: GETTAU_CL
!
c     setup unit matrix...............................................
c     ----------------------------------------------------------------
      call cmtruni(wbig,kkrsz_ns)
c     ----------------------------------------------------------------
c     get 1-delta and put it in wbig
      call zaxpy(mtxsize,cmone,delta,1,wbig,1)
c     ================================================================
c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
c     ----------------------------------------------------------------
      call zgetrf(kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,info)
      call zcopy(kkrsz_ns*kkrsz_ns,tmat,1,tau00,1)
      call zgetrs('n',kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,tau00,
     &           kkrsz_ns,info)
c     ----------------------------------------------------------------

c     ================================================================
c     Rotate tau00 to local frame of reference

!        Relativistic
c        -------------------------------------------------------------
!      call zcopy(4*kkrsz_ns*kkrsz_ns,tau00,1,tau00_l,1)
      call zcopy(kkrsz_ns*kkrsz_ns,tau00,1,tau00_l,1)
      call tripmt(dmatp,tau00_l,dmat,kkrsz_ns,kkrsz_ns,kkrsz_ns)

      end subroutine
