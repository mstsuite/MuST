module OrthoNormModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, TEN2m8, TEN2m14, CZERO
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initOrthoNorm, &
          endOrthoNorm,   &
          orthogonalize, &
          normalize,     &
          getNorm
!
   interface orthogonalize
      module procedure orthogonalize_sp, orthogonalize_fp
   end interface orthogonalize
!
   interface normalize
      module procedure normalize_sp, normalize_fp
   end interface normalize
!
   interface getNorm
      module procedure getNorm_sp, getNorm_fp
   end interface getNorm
!
private
   logical :: isFP
!
   real (kind=RealKind), allocatable :: f2(:)
!
   complex (kind=CmplxKind), allocatable :: f2fp(:,:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initOrthoNorm(nr,lmax,FP)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, lmax
!
   logical, intent(in) :: FP
!
   isFP = FP
!
   if (isFP) then
      allocate(f2fp(nr,(2*lmax+1)*(2*lmax+1)))
   else
      allocate(f2(nr))
   endif
!
   end subroutine initOrthoNorm
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endOrthoNorm()
!  ===================================================================
   implicit none
!
   if (isFP) then
      deallocate(f2fp)
   else
      deallocate(f2)
   endif
!
   end subroutine endOrthoNorm
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine orthogonalize_sp(info, n, m, r, f)
!  ===================================================================
!
!  Orthogonalize m complex vectors of length n using Gram-Schmidt
!  Input:
!     n : length of each function (or vector)
!     m : number of functions (or vectors)
!     r : real*8 radial grid, dimension (n)
!     f : complex*16 array, dimension (n,m)
!  Output:
!     f : orthonormalized functions (or vectors), which are modified 
!         in place
!  -------------------------------------------------------------------
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(:), n, m
!
   real (kind=RealKind), intent(in) :: r(n)
!
   complex (kind=CmplxKind), intent(inout) :: f(n,m)
!
   integer :: i, j, k
!
   real (kind=RealKind) :: norm
!
   complex (kind=CmplxKind) :: dotprod
!
   do j = 1, m
      ! Subtract projections onto previous vectors
      do k = 1, j-1
         dotprod = CZERO
         do i = 1, n
            dotprod = dotprod + dconjg(f(i,k)) * f(i,j)
         end do
         do i = 1, n
            f(i,j) = f(i,j) - dotprod * f(i,k)
         end do
      end do
   end do
!
   end subroutine orthogonalize_sp
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine orthogonalize_fp(info, n, lmax, m, r, f)
!  ===================================================================
!
!  This code is an extention of the orthogonalize subroutine shown 
!  above. It intends to perform orthogonalization of the wave function
!  without spherical symmetry. The 2nd dimension has angular momentum
!  quantum number index, resulting from an expansion on spherical harmonics.
!
!  This subroutine is simply copied from the orthogonalize subroutine
!  on Jan 15, 2026, and it needs to be modified and validated.
!
!  Orthogonalize m complex vectors of length n using Gram-Schmidt
!  Input:
!    n : length of each function (or vector)
!    kmax : the angular momentum quantum number cutoff
!    m : number of functions (or vectors)
!    r : real*8 radial grid, dimension (n)
!    f : complex*16 array, dimension (n,kmax,m)
!  Output:
!    f : orthonormalized functions (or vectors), which are modified 
!        in place
!  -------------------------------------------------------------------
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(:), n, lmax, m
!
   real (kind=RealKind), intent(in) :: r(n)
!
   complex (kind=CmplxKind), intent(inout) :: f(:,:,:)
!
   integer :: i, j, k, kl, kmax
!
   real (kind=RealKind) :: norm
!
   complex (kind=CmplxKind) :: dotprod
!
   interface
      subroutine computeProdExpan(n,lf,f,lg,g,lh,h)
        use KindParamModule, only : IntKind, CmplxKind
        integer (kind=IntKind), intent(in) :: n, lf, lg, lh
        complex (kind=CmplxKind), intent(in) :: f(:,:), g(:,:)
        complex (kind=CmplxKind), intent(out) :: h(:,:)
      end subroutine computeProdExpan
   end interface

   kmax = (lmax+1)**2

   do j = 1, m
      ! Subtract projections onto previous vectors
      ! =========================================================
      ! -- the following code needs to be checked -- Jan 15, 2026
      ! ---------------------------------------------------------
      do k = 1, j-1
         dotprod = CZERO
         do kl = 1, kmax
            do i = 1, n
               dotprod = dotprod + dconjg(f(i,kl,k)) * f(i,kl,j)
            end do
         enddo
         do kl = 1, kmax
            do i = 1, n
               f(i,kl,j) = f(i,kl,j) - dotprod * f(i,kl,k)
            end do
         enddo
      end do
    end do
!
   end subroutine orthogonalize_fp
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine normalize_sp(info,nr,r,f)
!  ===================================================================
!
!  compute the modulus of a given function f
!
!  *******************************************************************
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(:), nr
!
   real (kind=RealKind), intent(in) :: r(nr)
   real (kind=RealKind) :: fnorm
!
   complex (kind=CmplxKind), intent(inout) :: f(nr)
!
   fnorm = getNorm(info,nr,r,f)
!
   if (fnorm > TEN2m8) then
      f = f/fnorm
   else
      f = CZERO
   endif
!
   end subroutine normalize_sp
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine normalize_fp(info,nr,lmax,r,f)
!  ===================================================================
!
!  compute the modulus of a given function f
!
!  *******************************************************************
   use StepFunctionModule, only : getVolumeIntegration
   use PotentialTypeModule, only : isFullPotential
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(:), nr, lmax
   integer (kind=IntKind) :: ir, site, npow, kmax
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind) :: fnorm
!
   complex (kind=CmplxKind), intent(inout) :: f(:,:)
!
   fnorm = getNorm(info,nr,lmax,r,f)
!
   if (fnorm > TEN2m8) then
      f = f/fnorm
   else
      f = CZERO
   endif
!
   end subroutine normalize_fp
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNorm_sp(info,nr,r,f) result(s)
!  ===================================================================
!
!  compute the modulus of a given function f
!
!  *******************************************************************
   use StepFunctionModule, only : getVolumeIntegration
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(:), nr
   integer (kind=IntKind) :: ir, site, npow
!
   real (kind=RealKind), intent(in) :: r(nr)
   real (kind=RealKind) :: s, fnorm_ws, fnorm_mt
!
   complex (kind=CmplxKind), intent(in) :: f(nr)
!
   site = info(1); npow = info(2)
!
   if (npow == 0) then
      do ir = 1, nr
         f2(ir) = abs(f(ir))**2
      enddo
   else
      do ir = 1, nr
         f2(ir) = abs(f(ir))**2/r(ir)**npow
      enddo
   endif
!  -------------------------------------------------------------------
   fnorm_ws = getVolumeIntegration(site,nr,r,0,f2,fnorm_mt)
!  -------------------------------------------------------------------
   s = sqrt(fnorm_mt)
!
   end function getNorm_sp
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNorm_fp(info,nr,lmax,r,f) result(s)
!  ===================================================================
!
!  compute the modulus of a given function f
!
!  *******************************************************************
   use StepFunctionModule, only : getVolumeIntegration
   use PotentialTypeModule, only : isFullPotential
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(:), nr, lmax
   integer (kind=IntKind) :: ir, site, npow, kmax, kl
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind) :: s, fnorm_ws, fnorm_mt
!
   complex (kind=CmplxKind), intent(in) :: f(:,:)
!
   interface
      subroutine computeProdExpan(n,lf,f,lg,g,lh,h)
        use KindParamModule, only : IntKind, CmplxKind
        integer (kind=IntKind), intent(in) :: n, lf, lg, lh
        complex (kind=CmplxKind), intent(in) :: f(:,:), g(:,:)
        complex (kind=CmplxKind), intent(out) :: h(:,:)
      end subroutine computeProdExpan
   end interface
!
   site = info(1); npow = info(2)
!  ===================================================================
!  This needs to be checked:
!      f2fp is an expansion of a real function, so that jmax can be used
!      when calling getVolumeIntegration
!      Alternatively, getVolumeIntegration2 could be used for the volume integration
!  -------------------------------------------------------------------
   call ErrorHandler('getNorm_fp','This code is not working yet.')
   call computeProdExpan(nr,lmax,f,lmax,f,2*lmax,f2fp)
!  -------------------------------------------------------------------
   kmax = (2*lmax+1)*(2*lmax+1)
   do kl = 1, kmax
      do ir = 1, nr
         f2fp(ir,kl) = f2fp(ir,kl)/r(ir)**npow
      enddo
   enddo
!  -------------------------------------------------------------------
   fnorm_ws = getVolumeIntegration(site,nr,r,kmax,0,f2fp,fnorm_mt)
!  -------------------------------------------------------------------
   s = sqrt(fnorm_ws)
!
   end function getNorm_fp
!  ===================================================================
end module OrthoNormModule
