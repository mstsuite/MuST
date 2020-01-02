program tst_bessel
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : sqrtm1
!
   use BesselModule, only : SphericalBessel
   use BesselModule, only : SphericalNeumann
!
   implicit none
!
   integer (kind=IntKind) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind) :: xr, xi
!
   complex (kind=CmplxKind) :: x
!
   complex (kind=CmplxKind), allocatable :: jl(:), djl(:)
   complex (kind=CmplxKind), allocatable :: nl(:), dnl(:)
!
   write(6,*)' '
   write(6,*)'          ***************************************************'
   write(6,*)'          *                                                 *'
   write(6,*)'          *    TEST CODE TO CHECK THE WRONSKIAN RELATION    *'
   write(6,*)'          *                                                 *'
   write(6,*)'          *  OF SPHERICAL BESSEL FUNCTIONS: j_l(x), n_l(x)  *'
   write(6,*)'          *                                                 *'
   write(6,*)'          ***************************************************'
   write(6,*)' '
   write(6,'('' The maximum value of l (0, 1, 2,...)?  > '',$)')
   read(5,*)lmax
   write(6,'('' The Re part of the complex argument x? > '',$)')
   read(5,*)xr
   write(6,'('' The Im part of the complex argument x? > '',$)')
   read(5,*)xi
   write(6,*)' '
   write(6,*)'The correct Wronskian for each l should be equal to 1.0!'
   write(6,*)' '
   x=cmplx(xr,xi,CmplxKind)
!
   allocate ( jl(0:lmax), djl(0:lmax), nl(0:lmax), dnl(0:lmax) )
!
!  -------------------------------------------------------------------
   call SphericalBessel(lmax,x,jl,djl)
   call SphericalNeumann(lmax,x,nl,dnl)
!  -------------------------------------------------------------------
!
   write(6,'(/,a)')'Use SphericalBessel and SphericalNeumann'
   write(6,'(a)')'l,                         Wronskian'
   do l=0,lmax
      write(6,'(1i5,2d20.13)')l,x*x*(jl(l)*dnl(l)-djl(l)*nl(l))
   enddo
!
   write(6,'(/,''  l,                      jl,                          djl'')')
   do l=0,lmax
      write(6,'(1i4,4d19.11)')l,jl(l),djl(l)
   enddo
!
   write(6,'(/,''  l,                      nl,                          dnl'')')
   do l=0,lmax
      write(6,'(1i4,4d19.11)')l,nl(l),dnl(l)
   enddo
!
!  -------------------------------------------------------------------
   call csphjy(lmax,x,l,jl,djl,nl,dnl)
!  -------------------------------------------------------------------
!
   write(6,'(/,a)')'Use csphjy'
   write(6,'(a)')'l,                         Wronskian'
   do l=0,lmax
      write(6,'(1i5,2d20.13)')l,x*x*(jl(l)*dnl(l)-djl(l)*nl(l))
   enddo
!
   write(6,'(/,''  l,                      jl,                          djl'')')
   do l=0,lmax
      write(6,'(1i4,4d19.11)')l,jl(l),djl(l)
   enddo
!
   write(6,'(/,''  l,                      nl,                          dnl'')')
   do l=0,lmax
      write(6,'(1i4,4d19.11)')l,nl(l),dnl(l)
   enddo
!
!
   deallocate ( jl, djl, nl, dnl )
!
end program tst_bessel
