module SphericalHarmonicsModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use KindParamModule, only : CmplxKind
!
   use MathParamModule, only : zero
   use MathParamModule, only : half
   use MathParamModule, only : one
   use MathParamModule, only : two
   use MathParamModule, only : ten2m8
   use MathParamModule, only : ten2m10
   use MathParamModule, only : ten2m12
   use MathParamModule, only : ten2m14
   use MathParamModule, only : czero
   use MathParamModule, only : cone
   use MathParamModule, only : sqrtm1
   use MathParamModule, only : pi
   use MathParamModule, only : pi4
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use IntegerFactorsModule, only : lofj, mofj, m1m
!
   use LegendreModule, only : legendre
!
public :: initSphericalHarmonics, &
          endSphericalHarmonics,  &
          calYlm,                 &
          calYlmConjg,            &
          getClm
!
!  ===================================================================
!  generic function for calculating spherical harmonics
!  ===================================================================
   interface calYlm
      module procedure SphericalHarmonics0, SphericalHarmonics1,      &
                       SphericalHarmonics2, SphericalHarmonics3,      &
                       SphericalHarmonics4, SphericalHarmonics5
   end interface
!
!  ===================================================================
!  generic function for calculating the complex conjugate of spherical
!  harmonics
!  ===================================================================
   interface calYlmConjg
      module procedure SphericalHarmonicsCC0, SphericalHarmonicsCC1,  &
                       SphericalHarmonicsCC2, SphericalHarmonicsCC3
   end interface
!
private
   integer (kind=IntKind) :: lmax_max = -1
   integer (kind=IntKind) :: jmax_max = 0
!
   real (kind=RealKind), allocatable, target :: clm(:)
!
   real (kind=RealKind), allocatable :: plm(:)
   complex (kind=CmplxKind), allocatable :: e_imp(:)
!
   logical :: Initialized = .false.
   real (kind=RealKind) :: tol = half*ten2m12
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSphericalHarmonics(lmax)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   if (lmax < 0) then
      call ErrorHandler('initSphericalHarmonics','lmax < 0',lmax)
   else if (Initialized) then
      call endSphericalHarmonics()
   endif
!
   lmax_max = lmax
   jmax_max = ((lmax+1)*(lmax+2))/2
!
   allocate( clm(jmax_max),plm(1:jmax_max),e_imp(-lmax_max:lmax_max) )
   clm = ZERO
   plm = ZERO
   e_imp = CZERO
!
!  -------------------------------------------------------------------
   call initIntegerFactors(lmax)
!  -------------------------------------------------------------------
   call calClm(lmax,clm)
!  -------------------------------------------------------------------
   Initialized = .true.
!
   end subroutine initSphericalHarmonics
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSphericalHarmonics()
!  ===================================================================
   implicit none
!
   deallocate( clm, plm, e_imp )
   lmax_max = -1
   jmax_max =  0
!
   end subroutine endSphericalHarmonics
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genFactors(lmax)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: jl, l, m
!
!  ===================================================================
!  set up factors lofj and mofj for l up to lmax...................
!  ===================================================================
   jl=0
   do l=0,lmax
      do m=0,l
         jl=jl+1
         lofj(jl)=l
         mofj(jl)=m
      enddo
   enddo
!
!  ===================================================================
!  load the factors (-1)**m into m1m(-lmax:lmax)...................
!  ===================================================================
   m1m(0) = 1
   do m=1,lmax
      m1m(m)= m1m(m-1)*(-1)
      m1m(-m) = m1m(m)
   enddo
!
   end subroutine genFactors
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calClm(lmax,clm_local)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: m2
   integer (kind=IntKind) :: i,j,sgn
!
   real (kind=RealKind), intent(out) :: clm_local((lmax+1)*(lmax+2)/2)
   real (kind=RealKind) :: xfac
!
!  ===================================================================
!  Coeficients for complex spherical harmonics......................
!  Calclates all the c(l,m)'s up to lmax............................
!
!              m       [ (2*l+1)*(l-|m|)!]
!  c(l,m)= (-1)  * sqrt[-----------------]
!                      [   4*pi*(l+|m|)! ]
!
!  ====================================================================
   if(lmax < 0 ) then
      call ErrorHandler('calClm','lmax < 0',lmax)
   endif
!
   clm_local(1)=sqrt(one/pi4)
   if (lmax < 50) then
      do l=1,lmax
         xfac=sqrt(real((2*l+1),kind=RealKind)/pi4)
         sgn=1
         do m=0,l
            j=(l*(l+1))/2+m+1
            clm_local(j)=one
            m2=2*m
!           The following code for calculating the factorial will overflow for large l and m values.
            do i=1,m2
               clm_local(j)=(l-m+i)*clm_local(j)
            enddo
            clm_local(j)=sgn*xfac*sqrt(one/clm_local(j))
            sgn = -sgn
         enddo
      enddo
   else
      do l=1,lmax
         xfac=sqrt(real((2*l+1),kind=RealKind)/pi4)
         sgn=1
         do m=0,l
            j=(l*(l+1))/2+m+1
            m2=2*m
            clm_local(j)=ZERO
            do i=1,m2
               clm_local(j)=log(real(l-m+i,kind=RealKind))+clm_local(j)
            enddo
            if (clm_local(j) > 256) then
               clm_local(j)=ZERO
            else
               clm_local(j)=sgn*xfac*exp(-clm_local(j)/TWO)
            endif
            sgn = -sgn
         enddo
      enddo
   endif
!
   end subroutine calClm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getClm(lmax) result(pclm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: jmax
!
   real (kind=RealKind), pointer :: pclm(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getClm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('getClm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax = ((lmax+1)*(lmax+2))/2
   pclm => clm(1:jmax)
!
   end function getClm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonics0(x,y,z,lmax,ylm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: x,y,z
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=(lmax+1)*(lmax+2)/2
   kmax=(lmax+1)*(lmax+1)
   ylm = czero
!
   q2=x*x+y*y
   r=sqrt(q2+z*z)
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
   else if (q.lt.tol .and. z.gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. z.lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)*m1m(l)
      enddo
   else
      cos_phi=x/q
      sin_phi=y/q
      cos_the=z/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp(m)=exp(m*iphi)
      enddo
      e_imp(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp(-m)
         ylm(kl+m)=cp*e_imp(m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonics0
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonics1(vec,lmax,ylm)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: vec(3)
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
   ylm = czero
!
   jmax=((lmax+1)*(lmax+2))/2
   kmax=(lmax+1)*(lmax+1)
!
   q2=vec(1)*vec(1)+vec(2)*vec(2)
   r=sqrt(q2+vec(3)*vec(3))
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
   else if (q.lt.tol .and. vec(3).gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. vec(3).lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=(l+1)*(l+2)/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)*m1m(l)
      enddo
   else
      cos_phi=vec(1)/q
      sin_phi=vec(2)/q
      cos_the=vec(3)/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp(m)=exp(m*iphi)
      enddo
      e_imp(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp(-m)
         ylm(kl+m)=cp*e_imp(m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonics1
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonics2(x,y,z,lmax,plm_t,e_imp_t,ylm)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: x,y,z
   real (kind=RealKind), intent(out) :: plm_t((lmax+1)*(lmax+2)/2)
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: e_imp_t(-lmax:lmax)
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=((lmax+1)*(lmax+2))/2
   kmax=(lmax+1)*(lmax+1)
   ylm = czero
!
   q2=x*x+y*y
   r=sqrt(q2+z*z)
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
   else if (q.lt.tol .and. z.gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=czero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=one
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. z.lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=czero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=m1m(l)
         ylm(kl)=clm(jl)*plm_t(jl)
      enddo
   else
      cos_phi=x/q
      sin_phi=y/q
      cos_the=z/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp_t(m)=exp(m*iphi)
      enddo
      e_imp_t(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm_t)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm_t(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp_t(-m)
         ylm(kl+m)=cp*e_imp_t(m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonics2
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonics3(vec,lmax,plm_t,e_imp_t,ylm)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: vec(3)
   real (kind=RealKind), intent(out) :: plm_t((lmax+1)*(lmax+2)/2)
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: e_imp_t(-lmax:lmax)
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=((lmax+1)*(lmax+2))/2
   kmax=(lmax+1)*(lmax+1)
   ylm = czero
!
   q2=vec(1)*vec(1)+vec(2)*vec(2)
   r=sqrt(q2+vec(3)*vec(3))
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
   else if (q.lt.tol .and. vec(3).gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=czero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=one
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. vec(3).lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=czero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=m1m(l)
         ylm(kl)=clm(jl)*plm_t(jl)
      enddo
   else
      cos_phi=vec(1)/q
      sin_phi=vec(2)/q
      cos_the=vec(3)/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp_t(m)=exp(m*iphi)
      enddo
      e_imp_t(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm_t)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm_t(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp_t(-m)
         ylm(kl+m)=cp*e_imp_t(m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonics3
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonicsCC0(x,y,z,lmax,ylm)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: x,y,z
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=((lmax+1)*(lmax+2))/2
   kmax=(lmax+1)*(lmax+1)
   ylm = czero
!
   q2=x*x+y*y
   r=sqrt(q2+z*z)
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
   else if (q.lt.tol .and. z.gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. z.lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)*m1m(l)
      enddo
   else
      cos_phi=x/q
      sin_phi=y/q
      cos_the=z/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp(m)=exp(m*iphi)
      enddo
      e_imp(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp(m)
         ylm(kl+m)=cp*e_imp(-m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonicsCC0
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonicsCC1(vec,lmax,ylm)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: vec(3)
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=((lmax+1)*(lmax+2))/2
   kmax=(lmax+1)*(lmax+1)
   ylm = czero
!
   q2=vec(1)*vec(1)+vec(2)*vec(2)
   r=sqrt(q2+vec(3)*vec(3))
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
   else if (q.lt.tol .and. vec(3).gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. vec(3).lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)*m1m(l)
      enddo
   else
      cos_phi=vec(1)/q
      sin_phi=vec(2)/q
      cos_the=vec(3)/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp(m)=exp(m*iphi)
      enddo
      e_imp(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp(m)
         ylm(kl+m)=cp*e_imp(-m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonicsCC1
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonicsCC2(x,y,z,lmax,plm_t,e_imp_t,ylm)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: x,y,z
   real (kind=RealKind), intent(out) :: plm_t((lmax+1)*(lmax+2)/2)
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: e_imp_t(-lmax:lmax)
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=((lmax+1)*(lmax+2))/2
   kmax=(lmax+1)*(lmax+1)
   plm_t   = zero
   e_imp_t = czero
   ylm     = czero
!
   q2=x*x+y*y
   r=sqrt(q2+z*z)
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
      plm_t(1)=one
      plm_t(2:jmax)=zero
      e_imp_t(-lmax:lmax)=cone
   else if (q.lt.tol .and. z.gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=zero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=one
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. z.lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=zero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=m1m(l)
         ylm(kl)=clm(jl)*plm_t(jl)
      enddo
   else
      cos_phi=x/q
      sin_phi=y/q
      cos_the=z/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp_t(m)=exp(m*iphi)
      enddo
      e_imp_t(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm_t)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm_t(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp_t(m)
         ylm(kl+m)=cp*e_imp_t(-m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonicsCC2
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonicsCC3(vec,lmax,plm_t,e_imp_t,ylm)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl
!
   real (kind=RealKind), intent(in) :: vec(3)
   real (kind=RealKind), intent(out) :: plm_t((lmax+1)*(lmax+2)/2)
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp
!
   complex (kind=CmplxKind), intent(out) :: e_imp_t(-lmax:lmax)
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=((lmax+1)*(lmax+2))/2
   kmax=(lmax+1)*(lmax+1)
   plm_t   = zero
   e_imp_t = czero
   ylm     = czero
!
   q2=vec(1)*vec(1)+vec(2)*vec(2)
   r=sqrt(q2+vec(3)*vec(3))
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
      plm_t(1)=one
      plm_t(2:jmax)=zero
      e_imp_t(-lmax:lmax)=cone
   else if (q.lt.tol .and. vec(3).gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=zero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=one
         ylm(kl)=clm(jl)
      enddo
   else if (q.lt.tol .and. vec(3).lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      plm_t(1:jmax)=zero
      e_imp_t(-lmax:lmax)=cone
      do l=0,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         plm_t(jl)=m1m(l)
         ylm(kl)=clm(jl)*plm_t(jl)
      enddo
   else
      cos_phi=vec(1)/q
      sin_phi=vec(2)/q
      cos_the=vec(3)/r
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      do m=-lmax,lmax
         e_imp_t(m)=exp(m*iphi)
      enddo
      e_imp_t(0)=cone
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm_t)
!     ----------------------------------------------------------------
      do jl=1,jmax
         cp=clm(jl)*plm_t(jl)
         kl=(lofj(jl)+1)*(lofj(jl)+1)-lofj(jl)
         m=mofj(jl)
!         ylm(k-m)=m1m(m)*cp*e_imp_t(m)
         ylm(kl+m) = cp*e_imp_t(-m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
      enddo
   endif
!
   end subroutine SphericalHarmonicsCC3
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonics4(x,y,z,lmax,ylm,grady)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   integer (kind=IntKind) :: jmax
   integer (kind=IntKind) :: kmax
   integer (kind=IntKind) :: l,m,jl,kl,jlp1,jlm1
!
   real (kind=RealKind), intent(in) :: x, y, z
!
   real (kind=RealKind) :: r,q,q2
   real (kind=RealKind) :: sin_phi
   real (kind=RealKind) :: cos_phi
   real (kind=RealKind) :: sin_the
   real (kind=RealKind) :: cos_the
   real (kind=RealKind) :: cp, rfac
   real (kind=RealKind) :: t_unit(3), p_unit(3)
   real (kind=RealKind) :: plmp1, plmm1, tcoef, pcoef
!
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
   complex (kind=CmplxKind), intent(out) :: grady(:,:)
!
   complex (kind=CmplxKind) :: iphi
!
   if (.not.Initialized) then
      call ErrorHandler('calYlm','Module not initialized')
   else if (lmax > lmax_max) then
      call ErrorHandler('calYlm','lmax > lmax_max',lmax,lmax_max)
   endif
!
   jmax=(lmax+1)*(lmax+2)/2
   kmax=(lmax+1)*(lmax+1)
   grady = czero
   ylm   = czero
!
   q2=x*x+y*y
   r=sqrt(q2+z*z)
   q=sqrt(q2)
!
!  ===================================================================
!  calculating Ylm => ylm
!  ===================================================================
   if (r.lt.tol) then
      ylm(1)=clm(1)
      ylm(2:kmax)=czero
   else if (q.lt.tol .and. z.gt.zero) then   ! the point on +z-axis
      ylm(1:kmax)=czero
      grady(1:kmax,1)=czero
      grady(1:kmax,2)=czero
      grady(1:kmax,3)=czero
      ylm(1)=clm(1)
      do l=1,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)
         rfac=clm(jl+1)*((l*(l+1))/2)/r
         grady(kl+1,1)=rfac
         grady(kl+1,2)=rfac*sqrtm1
         grady(kl-1,1)=-rfac
         grady(kl-1,2)=rfac*sqrtm1
      enddo
   else if (q.lt.tol .and. z.lt.zero) then   ! the point on -z-axis
      ylm(1:kmax)=czero
      grady(1:kmax,1)=czero
      grady(1:kmax,2)=czero
      grady(1:kmax,3)=czero
      ylm(1)=clm(1)
      do l=1,lmax
         jl=((l+1)*(l+2))/2-l
         kl=(l+1)*(l+1)-l
         ylm(kl)=clm(jl)*m1m(l)
         rfac=m1m(l)*clm(jl+1)*((l*(l+1))/2)/r
         grady(kl+1,1)=-rfac
         grady(kl+1,2)=rfac*sqrtm1
         grady(kl-1,1)=rfac
         grady(kl-1,2)=rfac*sqrtm1
      enddo
   else
      cos_phi=x/q
      sin_phi=y/q
      cos_the=z/r
      sin_the=sqrt(one-cos_the*cos_the)
      t_unit(1)=cos_the*cos_phi
      t_unit(2)=cos_the*sin_phi
      t_unit(3)=-sin_the
      p_unit(1)=-sin_phi/sin_the
      p_unit(2)=cos_phi/sin_the
      p_unit(3)=zero
      iphi=sqrtm1*atan2(sin_phi,cos_phi)
      e_imp(0)=cone
      do m=1,lmax
         e_imp(m)=exp(m*iphi)
      enddo
      do m=-1,-lmax,-1
         e_imp(m)=cone/e_imp(-m)
      enddo
!     ----------------------------------------------------------------
      call legendre(lmax,cos_the,plm)
!     ----------------------------------------------------------------
      do jl=1,jmax
         l=lofj(jl)
         m=mofj(jl)
         cp=clm(jl)*plm(jl)
         kl=(l+1)*(l+1)-l
!         ylm(kl-m)=m1m(m)*cp*e_imp(-m)
         ylm(kl+m) = cp*e_imp(m)
         ylm(kl-m) = m1m(m)*conjg(ylm(kl+m))
!        =============================================================
!        determine plmm1=p(l,m-1), and plmp1=p(l,m+1).................
!        =============================================================
         jlp1=jl+1
         jlm1=jl-1
         pcoef=m*plm(jl)
         if (l.eq.0) then
            plmm1=zero
            plmp1=zero
            pcoef=zero
         else if (m.eq.0) then
            pcoef=zero
            plmm1=zero
            plmp1=two*plm(jlp1)
         else if (m.eq.l) then
            plmm1=plm(jlm1)
            plmp1=zero
         else 
            plmp1=plm(jlp1)
            plmm1=plm(jlm1)
         endif
         rfac=clm(jl)/r
         tcoef=-half*(plmp1-(l+m)*(l-m+1)*plmm1)
         grady(kl+m,1)=(t_unit(1)*tcoef+sqrtm1*p_unit(1)*pcoef)*rfac*e_imp(m)
         grady(kl+m,2)=(t_unit(2)*tcoef+sqrtm1*p_unit(2)*pcoef)*rfac*e_imp(m)
         grady(kl+m,3)=(t_unit(3)*tcoef+sqrtm1*p_unit(3)*pcoef)*rfac*e_imp(m)
         if (m .ne. 0) then
            grady(kl-m,1)=m1m(m)*(t_unit(1)*tcoef-                    &
                                  sqrtm1*p_unit(1)*pcoef)*rfac*e_imp(-m)
            grady(kl-m,2)=m1m(m)*(t_unit(2)*tcoef-                    &
                                  sqrtm1*p_unit(2)*pcoef)*rfac*e_imp(-m)
            grady(kl-m,3)=m1m(m)*(t_unit(3)*tcoef-                    &
                                  sqrtm1*p_unit(3)*pcoef)*rfac*e_imp(-m)
         endif
      enddo
   endif
!
   end subroutine SphericalHarmonics4
!  ===================================================================
!
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHarmonics5(vec,lmax,ylm,grady)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   real (kind=RealKind), intent(in) :: vec(3)
!
   complex (kind=CmplxKind), intent(out) :: ylm((lmax+1)*(lmax+1))
   complex (kind=CmplxKind), intent(out) :: grady(:,:)
!
   call SphericalHarmonics4(vec(1),vec(2),vec(3),lmax,ylm,grady)
!
   end subroutine SphericalHarmonics5
!  ===================================================================
end module SphericalHarmonicsModule
