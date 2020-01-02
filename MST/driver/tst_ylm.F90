program tst_ylm
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind) :: lmax, kmax
   integer (kind=IntKind) :: kl
!
   real (kind=RealKind) :: x, y, z
!
   complex (kind=CmplxKind), allocatable :: ylm(:), dylm(:,:)
!
   write(6,*)' '
   write(6,*)'          ***************************************************'
   write(6,*)'          *                                                 *'
   write(6,*)'          *    TEST CODE TO CHECK THE SPHERICAL HARMONIC    *'
   write(6,*)'          *                                                 *'
   write(6,*)'          *               FUNCTION: Y_lm(x,y,z)             *'
   write(6,*)'          *                                                 *'
   write(6,*)'          ***************************************************'
   write(6,*)' '
   write(6,'('' The maximum value of l (0, 1, 2, 3, ...) > '',$)')
   read(5,*)lmax
   write(6,'('' The (x, y, z) value, e.g., 1.0, 2.0, 3.0 > '',$)')
   read(5,*)x, y, z
   write(6,*)' '
!
   kmax=(lmax+1)*(lmax+1)
!
   allocate ( ylm(kmax), dylm(kmax,3) )
!
!  -------------------------------------------------------------------
   call calYlm(x,y,z,lmax,ylm,dylm)
!  -------------------------------------------------------------------
!
   write(6,'(/,a)')'Spherical Harmonic Function'
   write(6,'(a)')'   kl                         Ylm'
   do kl=1,kmax
      write(6,'(1i5,7x,2d20.13)')kl,ylm(kl)
   enddo
!
   write(6,'(/,a)')'Gradient of Spherical Harmonic Function'
   write(6,'(a)')'   kl               X-component of Grad Ylm'
   do kl=1,kmax
      write(6,'(1i5,7x,2d20.13)')kl,dylm(kl,1)
   enddo
   write(6,'(a)')'   kl               Y-component of Grad Ylm'
   do kl=1,kmax
      write(6,'(1i5,7x,2d20.13)')kl,dylm(kl,2)
   enddo
   write(6,'(a)')'   kl               Z-component of Grad Ylm'
   do kl=1,kmax
      write(6,'(1i5,7x,2d20.13)')kl,dylm(kl,3)
   enddo
!
   deallocate ( ylm, dylm )
!
end program tst_ylm
