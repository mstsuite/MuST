module BesselModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
   use MathParamModule, only : ten2m12, zero, one, two, cone, sqrtm1, HALF
!
public :: SphericalBessel,  &
          SphericalNeumann, &
          SphericalHankel,  &
          IntegrateSphHankelSq
!
!  define generic procedure for bessel functions
!
   interface SphericalBessel
      module procedure SphericalBesselReal0,  SphericalBesselReal1,   &
                       SphericalBesselReal0_v1,                       &
                       SphericalBesselCmplx0, SphericalBesselCmplx1,  &
                       SphericalBesselCmplx0_v1
   end interface
!
   interface SphericalNeumann
      module procedure SphericalNeumannReal0, SphericalNeumannReal1,  &
                       SphericalNeumannCmplx0,SphericalNeumannCmplx1
   end interface
!
   interface SphericalHankel
      module procedure SphericalHankelReal0,  SphericalHankelReal1,   &
                       SphericalHankelCmplx0, SphericalHankelCmplx1
   end interface
!
   interface IntegrateSphHankelSq
      module procedure IntegrateSphHankelSq_Cs, IntegrateSphHankelSq_Ca, &
                       IntegrateSphHankelSq_Rs, IntegrateSphHankelSq_Ra
   end interface
!
private
   integer (kind=IntKind) :: lmaxr_sav = -1
   integer (kind=IntKind) :: lmaxc_sav = -1
   integer (kind=IntKind) :: lmaxjlr_sav = -1
   integer (kind=IntKind) :: lmaxdjlr_sav = -1
   integer (kind=IntKind) :: lmaxjlc_sav = -1
   integer (kind=IntKind) :: lmaxdjlc_sav = -1
!
   real (kind=RealKind) :: tol = ten2m12
   real (kind=RealKind), allocatable :: bscratr(:)
   real (kind=RealKind), allocatable :: jlr(:), djlr(:), nlr(:), dnlr(:)
!
   complex (kind=CmplxKind), allocatable :: bscratc(:)
   complex (kind=CmplxKind), allocatable :: jlc(:), djlc(:), nlc(:), dnlc(:)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalBesselCmplx0(lmax,x,bj)
!  ===================================================================
   use CmplxFunctionModule, only : ImaginaryPart, RealPart
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: k, l, nm
!
   complex (kind=CmplxKind), intent(in) :: x
   complex (kind=CmplxKind), target :: bj(0:lmax)
   complex (kind=CmplxKind) :: a, ac, as, aj, x1, x2
!
!  *******************************************************************
!  calculate spherical bessel functions j_l
!  input: lmax -- max l
!         x    -- compex argument of j_l
!  returns: bj -- j_l
!  bscrat is as a scratch space
!  *******************************************************************
!
   if ( abs(x) .lt. tol ) then
      bj(0) = ONE
      if ( lmax>0 ) then
         bj(1:lmax) = ZERO
      endif
      return
   endif
!
   as=sin(x)
   ac=cos(x)
   x1=cone/x
   bj(0)=as*x1
   if (lmax.eq.0) then
      return
   else if(lmax.eq.1) then
      bj(1) = (bj(0)-ac)*x1
      return
   endif
!
   if (.not.allocated(bscratc)) then
      allocate( bscratc(0:lmax) )
      lmaxc_sav = lmax
   else if (lmaxc_sav < lmax) then
      deallocate( bscratc )
      allocate( bscratc(0:lmax) )
      lmaxc_sav = lmax
   endif
!
!  ===================================================================
!  Forward recursion for small L
!  =================================================================== 
   k=ishft(int(0.75*(abs(RealPart(x))+abs(ImaginaryPart(x))))+2,-2)
   bscratc(0)=as
   if (k .ge. 1) then
      if (k .gt. lmax) then
         k=lmax
      endif
      bscratc(1)=bj(0)-ac
      bj(1)=bscratc(1)*x1
      do l=2,k
         bscratc(l)=(2*l-1)*bj(l-1)-bscratc(l-2)
         bj(l)=bscratc(l)*x1
      enddo
      if (k .eq. lmax) then
         return
      endif
   endif
!  ===================================================================
!  Backward recursion from very large L down to L=k
!  ===================================================================
   a=bscratc(k)
   nm=ishft(lmax,4)
   aj=2*nm+3
   x2=x*x
   do l=nm,lmax+2,-1
      aj=(2*l+1)-x2/aj
   enddo
   bscratc(lmax)=(2*lmax+3)*aj*x1-x
   bj(lmax)=bscratc(lmax)*x1
   bscratc(lmax-1)=(2*lmax+1)*bj(lmax)-aj
   do l=lmax-1,k+1,-1
      bj(l)=bscratc(l)*x1
      bscratc(l-1)=(2*l+1)*bj(l)-bscratc(l+1)
   enddo
!  ===================================================================
!  scale to give correct bj(k)
!  ===================================================================
   aj=a/bscratc(k)
   do l=k+1,lmax
      bj(l)=aj*bj(l)
   enddo
!
   end subroutine SphericalBesselCmplx0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalBesselCmplx0_v1(lmax,nx,x,bj)
!  ===================================================================
   use CmplxFunctionModule, only : ImaginaryPart, RealPart
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax, nx
   integer (kind=IntKind) :: k, l, nm, n
!
   complex (kind=CmplxKind), target :: x(nx)
   complex (kind=CmplxKind), target :: bj(1:nx,0:lmax)
   complex (kind=CmplxKind) :: a, ac, as, aj, x1, x2
!
!  *******************************************************************
!  calculate spherical bessel functions j_l
!  input: lmax -- max l
!         x    -- compex argument of j_l
!  returns: bj -- j_l
!  bscrat is as a scratch space
!  *******************************************************************
!
   if (.not.allocated(bscratc)) then
      allocate( bscratc(0:lmax) )
      lmaxc_sav = lmax
   else if (lmaxc_sav < lmax) then
      deallocate( bscratc )
      allocate( bscratc(0:lmax) )
      lmaxc_sav = lmax
   endif
!
   LOOP_N : do n = 1,nx
      if ( abs(x(n)) .lt. tol ) then
         bj(n,0) = ONE
         if ( lmax>0 ) then
            do l = 1,lmax
               bj(n,l) = ZERO
            enddo
         endif
         cycle LOOP_N
      endif
!
      as=sin(x(n))
      ac=cos(x(n))
      x1=cone/x(n)
      bj(n,0)=as*x1
      if (lmax.eq.0) then
         cycle LOOP_N
      else if(lmax.eq.1) then
         bj(n,1) = (bj(n,0)-ac)*x1
         cycle LOOP_N
      endif
!
!     ================================================================
!     Forward recursion for small L
!     ================================================================ 
      k=ishft(int(0.75*(abs(RealPart(x(n)))+abs(ImaginaryPart(x(n)))))+2,-2)
      bscratc(0)=as
      if (k .ge. 1) then
         if (k .gt. lmax) then
            k=lmax
         endif
         bscratc(1)=bj(n,0)-ac
         bj(n,1)=bscratc(1)*x1
         do l=2,k
            bscratc(l)=(2*l-1)*bj(n,l-1)-bscratc(l-2)
            bj(n,l)=bscratc(l)*x1
         enddo
         if (k .eq. lmax) then
            cycle LOOP_N
         endif
      endif
!     ================================================================
!     Backward recursion from very large L down to L=k
!     ================================================================
      a=bscratc(k)
      nm=ishft(lmax,4)
      aj=2*nm+3
      x2=x(n)*x(n)
      do l=nm,lmax+2,-1
         aj=(2*l+1)-x2/aj
      enddo
      bscratc(lmax)=(2*lmax+3)*aj*x1-x(n)
      bj(n,lmax)=bscratc(lmax)*x1
      bscratc(lmax-1)=(2*lmax+1)*bj(n,lmax)-aj
      do l=lmax-1,k+1,-1
         bj(n,l)=bscratc(l)*x1
         bscratc(l-1)=(2*l+1)*bj(n,l)-bscratc(l+1)
      enddo
!     ================================================================
!     scale to give correct bj(k)
!     ================================================================
      aj=a/bscratc(k)
      do l=k+1,lmax
         bj(n,l)=aj*bj(n,l)
      enddo
   enddo LOOP_N
!
   end subroutine SphericalBesselCmplx0_v1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalBesselCmplx1(lmax,x,bj,dj)
!  ===================================================================
   use CmplxFunctionModule, only : ImaginaryPart, RealPart
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: k, l, nm
!
   complex (kind=CmplxKind), intent(in) :: x
   complex (kind=CmplxKind), target :: bj(0:lmax)
   complex (kind=CmplxKind), target :: dj(0:lmax)
   complex (kind=CmplxKind) :: a, ac, as, aj, x1, x2
!
!  *******************************************************************
!  calculate spherical bessel functions j_l
!  input: lmax -- max l
!         x    -- compex argument of j_l
!  returns: bj -- j_l
!  bscrat is as a scratch space
!  *******************************************************************
!
   if ( abs(x) .lt. tol ) then
      bj(0) = ONE
      dj(0) = ONE
      if ( lmax>0 ) then
         bj(1:lmax) = ZERO
         dj(1:lmax) = ZERO
      endif
      return
   endif
!
   as=sin(x)
   ac=cos(x)
   x1=cone/x
   bj(0)=as*x1
   dj(0)=(ac-bj(0))*x1
   if (lmax.eq.0) then
      return
   else if(lmax.eq.1) then
      bj(1) = -dj(0)
      dj(1) = bj(0)-two*bj(1)*x1
      return
   endif
!
   if (.not.allocated(bscratc)) then
      allocate( bscratc(0:lmax) )
      lmaxc_sav = lmax
   else if (lmaxc_sav < lmax) then
      deallocate( bscratc )
      allocate( bscratc(0:lmax) )
      lmaxc_sav = lmax
   endif
!
!  ===================================================================
!  Forward recursion for small L
!  =================================================================== 
   k=ishft(int(0.75*(abs(RealPart(x))+abs(ImaginaryPart(x))))+2,-2)
   bscratc(0)=as
   if (k .ge. 1) then
      if (k .gt. lmax) then
         k=lmax
      endif
      bscratc(1)=bj(0)-ac
      bj(1)=bscratc(1)*x1
      dj(1)=bj(0)-two*bj(1)*x1
      do l=2,k
         bscratc(l)=(2*l-1)*bj(l-1)-bscratc(l-2)
         bj(l)=bscratc(l)*x1
         dj(l)=bj(l-1)-(l+1)*bj(l)*x1
      enddo
      if (k .eq. lmax) then
         return
      endif
   endif
!  ===================================================================
!  Backward recursion from very large L down to L=k
!  ===================================================================
   a=bscratc(k)
   nm=ishft(lmax,4)
   aj=2*nm+3
   x2=x*x
   do l=nm,lmax+2,-1
      aj=(2*l+1)-x2/aj
   enddo
   bscratc(lmax)=(2*lmax+3)*aj*x1-x
   bj(lmax)=bscratc(lmax)*x1
   bscratc(lmax-1)=(2*lmax+1)*bj(lmax)-aj
   do l=lmax-1,k+1,-1
      bj(l)=bscratc(l)*x1
      bscratc(l-1)=(2*l+1)*bj(l)-bscratc(l+1)
   enddo
!  ===================================================================
!  scale to give correct bj(k)
!  ===================================================================
   aj=a/bscratc(k)
   do l=k+1,lmax
      bj(l)=aj*bj(l)
      dj(l)=bj(l-1)-(l+1)*bj(l)*x1
   enddo
!
   end subroutine SphericalBesselCmplx1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalBesselReal0(lmax,x,bj)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: k, l, nm
!
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind), target :: bj(0:lmax)
   real (kind=RealKind) :: a, ac, as, aj, x1, x2
!
!  *******************************************************************
!  calculate spherical bessel functions j_l
!  input: lmax -- max l
!         x    -- compex argument of j_l
!  returns: bj -- j_l
!  bscrat is as a scratch space
!  *******************************************************************
!
   if ( abs(x) .lt. tol ) then
      bj(0) = ONE
      if ( lmax>0 ) then
         bj(1:lmax) = ZERO
      endif
      return
   endif
!
   as=sin(x)
   ac=cos(x)
   x1=one/x
   bj(0)=as*x1
   if (lmax.eq.0) then
      return
   else if(lmax.eq.1) then
      bj(1) = (bj(0)-ac)*x1
      return
   endif
!
   if (.not.allocated(bscratr)) then
      allocate( bscratr(0:lmax) )
      lmaxr_sav = lmax
   else if (lmaxr_sav < lmax) then
      deallocate( bscratr )
      allocate( bscratr(0:lmax) )
      lmaxr_sav = lmax
   endif
!
!  ===================================================================
!  Forward recursion for small L
!  =================================================================== 
   k=ishft(int(0.75*abs(x))+2,-2)
   bscratr(0)=as
   if (k .ge. 1) then
      if (k .gt. lmax) then
         k=lmax
      endif
      bscratr(1)=bj(0)-ac
      bj(1)=bscratr(1)*x1
      do l=2,k
         bscratr(l)=(2*l-1)*bj(l-1)-bscratr(l-2)
         bj(l)=bscratr(l)*x1
      enddo
      if (k .eq. lmax) then
         return
      endif
   endif
!  ===================================================================
!  Backward recursion from very large L down to L=k
!  ===================================================================
   a=bscratr(k)
   nm=ishft(lmax,4)
   aj=2*nm+3
   x2=x*x
   do l=nm,lmax+2,-1
      aj=(2*l+1)-x2/aj
   enddo
   bscratr(lmax)=(2*lmax+3)*aj*x1-x
   bj(lmax)=bscratr(lmax)*x1
   bscratr(lmax-1)=(2*lmax+1)*bj(lmax)-aj
   do l=lmax-1,k+1,-1
      bj(l)=bscratr(l)*x1
      bscratr(l-1)=(2*l+1)*bj(l)-bscratr(l+1)
   enddo
!  ===================================================================
!  scale to give correct bj(k)
!  ===================================================================
   aj=a/bscratr(k)
   do l=k+1,lmax
      bj(l)=aj*bj(l)
   enddo
!
   end subroutine SphericalBesselReal0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalBesselReal0_v1(lmax,nx,x,bj)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax, nx
   integer (kind=IntKind) :: k, l, nm, n
!
   real (kind=RealKind), target :: x(nx)
   real (kind=RealKind), target :: bj(1:nx,0:lmax)
   real (kind=RealKind) :: a, ac, as, aj, x1, x2
!
!  *******************************************************************
!  calculate spherical bessel functions j_l
!  input: lmax -- max l
!         x    -- compex argument of j_l
!  returns: bj -- j_l
!  bscrat is as a scratch space
!  *******************************************************************
!
   if (.not.allocated(bscratr)) then
      allocate( bscratr(0:lmax) )
      lmaxr_sav = lmax
   else if (lmaxr_sav < lmax) then
      deallocate( bscratr )
      allocate( bscratr(0:lmax) )
      lmaxr_sav = lmax
   endif
!
   LOOP_N : do n = 1,nx
      if ( abs(x(n)) .lt. tol ) then
         bj(n,0) = ONE
         if ( lmax>0 ) then
            do l = 1,lmax
               bj(n,l) = ZERO
            enddo
         endif
         cycle LOOP_N
      endif
!
      as=sin(x(n))
      ac=cos(x(n))
      x1=one/x(n)
      bj(n,0)=as*x1
      if (lmax.eq.0) then
         cycle LOOP_N
      else if(lmax.eq.1) then
         bj(n,1) = (bj(n,0)-ac)*x1
         cycle LOOP_N
      endif
!
!     ================================================================
!     Forward recursion for small L
!     ================================================================ 
      k = ishft(int(0.75*abs(x(n)))+2,-2)
      bscratr(0)=as
      if (k .ge. 1) then
         if (k .gt. lmax) then
            k=lmax
         endif
         bscratr(1)=bj(n,0)-ac
         bj(n,1)=bscratr(1)*x1
         do l=2,k
            bscratr(l)=(2*l-1)*bj(n,l-1)-bscratr(l-2)
            bj(n,l)=bscratr(l)*x1
         enddo
         if (k .eq. lmax) then
            cycle LOOP_N
         endif
      endif
!     ================================================================
!     Backward recursion from very large L down to L=k
!     ================================================================
      a=bscratr(k)
      nm=ishft(lmax,4)
      aj=2*nm+3
      x2=x(n)*x(n)
      do l=nm,lmax+2,-1
         aj=(2*l+1)-x2/aj
      enddo
      bscratr(lmax)=(2*lmax+3)*aj*x1-x(n)
      bj(n,lmax)=bscratr(lmax)*x1
      bscratr(lmax-1)=(2*lmax+1)*bj(n,lmax)-aj
      do l=lmax-1,k+1,-1
         bj(n,l)=bscratr(l)*x1
         bscratr(l-1)=(2*l+1)*bj(n,l)-bscratr(l+1)
      enddo
!     ================================================================
!     scale to give correct bj(k)
!     ================================================================
      aj=a/bscratr(k)
      do l=k+1,lmax
         bj(n,l)=aj*bj(n,l)
      enddo
   enddo LOOP_N
!
   end subroutine SphericalBesselReal0_v1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalBesselReal1(lmax,x,bj,dj)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: k, l, nm
!
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind), target :: bj(0:lmax)
   real (kind=RealKind), target :: dj(0:lmax)
   real (kind=RealKind) :: a, ac, as, aj, x1, x2
!
!  *******************************************************************
!  calculate spherical bessel functions j_l
!  input: lmax -- max l
!         x    -- compex argument of j_l
!  returns: bj -- j_l
!  bscrat is as a scratch space
!  *******************************************************************
!
   if (abs(x) .lt. tol) then
      bj(0) = ONE
      dj(0) = ONE
      if ( lmax>0 ) then
         bj(1:lmax) = ZERO
         dj(1:lmax) = ZERO
      endif
      return
   endif
!
   as=sin(x)
   ac=cos(x)
   x1=one/x
   bj(0)=as*x1
   dj(0)=(ac-bj(0))*x1
   if (lmax.eq.0) then
      return
   else if(lmax.eq.1) then
      bj(1) = (bj(0)-ac)*x1
      dj(1) = bj(0)-two*bj(1)*x1
      return
   endif
!
   if (.not.allocated(bscratr)) then
      allocate( bscratr(0:lmax) )
      lmaxr_sav = lmax
   else if (lmaxr_sav < lmax) then
      deallocate( bscratr )
      allocate( bscratr(0:lmax) )
      lmaxr_sav = lmax
   endif
!
!  ===================================================================
!  Forward recursion for small L
!  =================================================================== 
   k=ishft(int(0.75*abs(x))+2,-2)
   bscratr(0)=as
   if (k .ge. 1) then
      if (k .gt. lmax) then
         k=lmax
      endif
      bscratr(1)=bj(0)-ac
      bj(1)=bscratr(1)*x1
      dj(1)=bj(0)-two*bj(1)*x1
      do l=2,k
         bscratr(l)=(2*l-1)*bj(l-1)-bscratr(l-2)
         bj(l)=bscratr(l)*x1
         dj(l)=bj(l-1)-(l+1)*bj(l)*x1
      enddo
      if (k .eq. lmax) then
         return
      endif
   endif
!  ===================================================================
!  Backward recursion from very large L down to L=k
!  ===================================================================
   a=bscratr(k)
   nm=ishft(lmax,4)
   aj=2*nm+3
   x2=x*x
   do l=nm,lmax+2,-1
      aj=(2*l+1)-x2/aj
   enddo
   bscratr(lmax)=(2*lmax+3)*aj*x1-x
   bj(lmax)=bscratr(lmax)*x1
   bscratr(lmax-1)=(2*lmax+1)*bj(lmax)-aj
   do l=lmax-1,k+1,-1
      bj(l)=bscratr(l)*x1
      bscratr(l-1)=(2*l+1)*bj(l)-bscratr(l+1)
   enddo
!  ===================================================================
!  scale to give correct bj(k)
!  ===================================================================
   aj=a/bscratr(k)
   do l=k+1,lmax
      bj(l)=aj*bj(l)
      dj(l)=bj(l-1)-(l+1)*bj(l)*x1
   enddo
!
   end subroutine SphericalBesselReal1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalNeumannCmplx0(lmax,x,bn)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   complex (kind=CmplxKind), intent(in) :: x
   complex (kind=CmplxKind), target :: bn(0:lmax)
   complex (kind=CmplxKind) :: x1
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalNeumann', 'x is near zero ',x)
   endif
!
   x1=cone/x
   bn(0) = -cos(x)*x1
   if (lmax .eq. 0) then
      return
   endif
   bn(1) = (bn(0) - sin(x))*x1
!  ===================================================================
!  recursion relations
!  ===================================================================
   do l=2,lmax
      bn(l) = (2*l-1)*bn(l-1)*x1-bn(l-2)
   enddo
!
   end subroutine SphericalNeumannCmplx0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalNeumannCmplx1(lmax,x,bn,dn)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   complex (kind=CmplxKind), intent(in) :: x
   complex (kind=CmplxKind), target :: bn(0:lmax)
   complex (kind=CmplxKind), target :: dn(0:lmax)
   complex (kind=CmplxKind) :: x1
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalNeumann', 'x is near zero ',x)
   endif
!
   x1=cone/x
   bn(0) = -cos(x)*x1
   dn(0) = (sin(x)-bn(0))*x1
   if (lmax .eq. 0) then
      return
   endif
   bn(1) = -dn(0)
   dn(1) = bn(0)-two*bn(1)*x1
!  ===================================================================
!  recursion relations
!  ===================================================================
   do l=2,lmax
      bn(l) = (2*l-1)*bn(l-1)*x1-bn(l-2)
      dn(l) = bn(l-1)-(l+1)*bn(l)*x1
   enddo
!
   end subroutine SphericalNeumannCmplx1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalNeumannReal0(lmax,x,bn)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind), target :: bn(0:lmax)
   real (kind=RealKind) :: x1
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalNeumann', 'x is near zero ',x)
   endif
!
   x1=one/x
   bn(0) = -cos(x)*x1
   if (lmax .eq. 0) then
      return
   endif
   bn(1) = (bn(0) - sin(x))*x1
!  ===================================================================
!  recursion relations
!  ===================================================================
   do l=2,lmax
      bn(l) = (2*l-1)*bn(l-1)*x1-bn(l-2)
   enddo
!
   end subroutine SphericalNeumannReal0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalNeumannReal1(lmax,x,bn,dn)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: x
   real (kind=RealKind), target :: bn(0:lmax)
   real (kind=RealKind), target :: dn(0:lmax)
   real (kind=RealKind) :: x1
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalNeumann', 'x is near zero ',x)
   endif
!
   x1=one/x
   bn(0) = -cos(x)*x1
   dn(0) = (sin(x)-bn(0))*x1
   if (lmax .eq. 0) then
      return
   endif
   bn(1) = -dn(0)
   dn(1) = bn(0)-two*bn(1)*x1
!  ===================================================================
!  recursion relations
!  ===================================================================
   do l=2,lmax
      bn(l) = (2*l-1)*bn(l-1)*x1-bn(l-2)
      dn(l) = bn(l-1)-(l+1)*bn(l)*x1
   enddo
!
   end subroutine SphericalNeumannReal1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHankelReal0(lmax,x,hl)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   real (kind=RealKind), intent(in) :: x
!
   complex (kind=CmplxKind), target :: hl(0:lmax)
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalHankel', 'x is near zero ',x)
   endif
!
   hl(0) = (sin(x)-sqrtm1*cos(x))/x
   if (lmax .eq. 0) then
      return
   else if (lmax.eq.1) then
      hl(1)=(-sqrtm1+cone/x)*hl(0)
      return
   endif
!
   if (.not.allocated(jlr)) then
      allocate( jlr(0:lmax), nlr(0:lmax) )
      lmaxjlr_sav = lmax
   else if (lmaxjlr_sav < lmax) then
      deallocate( jlr, nlr )
      allocate( jlr(0:lmax), nlr(0:lmax) )
      lmaxjlr_sav = lmax
   endif
!
!  -------------------------------------------------------------------
   call SphericalBessel(lmax,x,jlr)
   call SphericalNeumann(lmax,x,nlr)
!  -------------------------------------------------------------------
   hl(1:lmax)=jlr(1:lmax)+sqrtm1*nlr(1:lmax)
!
   end subroutine SphericalHankelReal0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHankelReal1(lmax,x,hl,dhl)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   real (kind=RealKind), intent(in) :: x
!
   complex (kind=CmplxKind), target :: hl(0:lmax), dhl(0:lmax)
   complex (kind=CmplxKind) :: x1
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalHankel', 'x is near zero ',x)
   endif
!
   x1=cone/x
   hl(0) = (sin(x)-sqrtm1*cos(x))*x1
   dhl(0)= (sqrtm1-x1)*hl(0)
   if (lmax .eq. 0) then
      return
   else if (lmax.eq.1) then
      hl(1)=-dhl(0)
      dhl(1)=hl(0)-two*hl(1)*x1
      return
   endif
!
   if (.not.allocated(jlr)) then
      allocate( jlr(0:lmax), nlr(0:lmax) )
      lmaxjlr_sav = lmax
   else if (lmaxjlr_sav < lmax) then
      deallocate( jlr, nlr )
      allocate( jlr(0:lmax), nlr(0:lmax) )
      lmaxjlr_sav = lmax
   endif
!
   if (.not.allocated(djlr)) then
      allocate( djlr(0:lmax), dnlr(0:lmax) )
      lmaxdjlr_sav = lmax
   else if (lmaxdjlr_sav < lmax) then
      deallocate( djlr, dnlr )
      allocate( djlr(0:lmax), dnlr(0:lmax) )
      lmaxdjlr_sav = lmax
   endif
!
!  -------------------------------------------------------------------
   call SphericalBessel(lmax,x,jlr,djlr)
   call SphericalNeumann(lmax,x,nlr,dnlr)
!  -------------------------------------------------------------------
   hl(1:lmax)=jlr(1:lmax)+sqrtm1*nlr(1:lmax)
   dhl(1:lmax)=djlr(1:lmax)+sqrtm1*dnlr(1:lmax)
!
   end subroutine SphericalHankelReal1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHankelCmplx0(lmax,x,hl)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   complex (kind=CmplxKind), intent(in) :: x
   complex (kind=CmplxKind) :: x1
!
   complex (kind=CmplxKind), target :: hl(0:lmax)
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalHankel', 'x is near zero ',x)
   endif
!
   x1=cone/x
   hl(0) = (sin(x)-sqrtm1*cos(x))*x1
   if (lmax .eq. 0) then
      return
   else if (lmax.eq.1) then
      hl(1)=(-sqrtm1+x1)*hl(0)
      return
   endif
!
   if (.not.allocated(jlc)) then
      allocate( jlc(0:lmax), nlc(0:lmax) )
      lmaxjlc_sav = lmax
   else if (lmaxjlc_sav < lmax) then
      deallocate( jlc, nlc )
      allocate( jlc(0:lmax), nlc(0:lmax) )
      lmaxjlc_sav = lmax
   endif
!
!  -------------------------------------------------------------------
   call SphericalBessel(lmax,x,jlc)
   call SphericalNeumann(lmax,x,nlc)
!  -------------------------------------------------------------------
   hl(1:lmax)=jlc(1:lmax)+sqrtm1*nlc(1:lmax)
!
   end subroutine SphericalHankelCmplx0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine SphericalHankelCmplx1(lmax,x,hl,dhl)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
!
   complex (kind=CmplxKind), intent(in) :: x
   complex (kind=CmplxKind) :: x1
!
   complex (kind=CmplxKind), target :: hl(0:lmax), dhl(0:lmax)
!
!  *******************************************************************
!  computes spherical Neumann function n_l............................
!  *******************************************************************
   if (abs(x) .lt. tol) then
      call ErrorHandler('SphericalHankel', 'x is near zero ',x)
   endif
!
   x1=cone/x
   hl(0) = (sin(x)-sqrtm1*cos(x))/x
   dhl(0)= (sqrtm1-x1)*hl(0)
   if (lmax .eq. 0) then
      return
   else if (lmax.eq.1) then
      hl(1)=-dhl(0)
      dhl(1)=hl(0)-two*hl(1)*x1
      return
   endif
!
   if (.not.allocated(jlc)) then
      allocate( jlc(0:lmax), nlc(0:lmax) )
      lmaxjlc_sav = lmax
   else if (lmaxjlc_sav < lmax) then
      deallocate( jlc, nlc )
      allocate( jlc(0:lmax), nlc(0:lmax) )
      lmaxjlc_sav = lmax
   endif
!
   if (.not.allocated(djlc)) then
      allocate( djlc(0:lmax), dnlc(0:lmax) )
      lmaxdjlc_sav = lmax
   else if (lmaxdjlc_sav < lmax) then
      deallocate( djlc, dnlc )
      allocate( djlc(0:lmax), dnlc(0:lmax) )
      lmaxdjlc_sav = lmax
   endif
!
!  -------------------------------------------------------------------
   call SphericalBessel(lmax,x,jlc,djlc)
   call SphericalNeumann(lmax,x,nlc,dnlc)
!  -------------------------------------------------------------------
   hl(1:lmax)=jlc(1:lmax)+sqrtm1*nlc(1:lmax)
   dhl(1:lmax)=djlc(1:lmax)+sqrtm1*dnlc(1:lmax)
!
   end subroutine SphericalHankelCmplx1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine IntegrateSphHankelSq_Ca(lmax,rc,energy,fint)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: rc
!
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind), intent(out) :: fint(0:lmax)
!  complex (kind=CmplxKind) :: bessjl(0:lmax+1), bessnl(0:lmax+1), x, besshl, besshlp1, besshlm1
   complex (kind=CmplxKind) :: besshl(0:lmax+1), x
   complex (kind=CmplxKind) :: Ul, Ulp1
!
!  x = rc*sqrt(energy)+SQRTm1*0.00001d0
   x = rc*sqrt(energy)
!  -------------------------------------------------------------------
!  call SphericalBessel(lmax+1,x,bessjl)
!  call SphericalNeumann(lmax+1,x,bessnl)
   call SphericalHankel(lmax+1,x,besshl)
!  -------------------------------------------------------------------
!  besshl = bessjl(0)+SQRTm1*bessnl(0)
!  besshlp1 = bessjl(1)+SQRTm1*bessnl(1)
!  fint(0) = (besshl*besshlp1/x-besshl**2-besshlp1**2)*HALF
   fint(0) = (besshl(0)*besshl(1)/x-besshl(0)**2-besshl(1)**2)*HALF
   do l = 1, lmax
!     besshlm1 = besshl
!     besshl = besshlp1
!     besshlp1 = bessjl(l+1)+SQRTm1*bessnl(l+1)
!!    fint(l) = ((2*l+1)*besshl*besshlp1/x-besshl**2-besshlp1**2)*HALF
!     fint(l) = (besshlm1*besshlp1-besshl**2)*HALF
      fint(l) = (besshl(l-1)*besshl(l+1)-besshl(l)**2)*HALF
   enddo
!
!  ===================================================================
!  Another way of calculating fint is to use recursive relation.
!  Both ways have shown to give the same results.
!  ===================================================================
!! Ul = -SQRTm1*exp(SQRTm1*TWO*x)/(TWO*x**3)
!! write(6,'(a,2d16.8,2x,2d16.8)')'fint(0), Ul = ',fint(0), Ul
!! fint(0) = Ul
!! do l = 1, lmax
!!    besshlm1 = bessjl(l-1)+SQRTm1*bessnl(l-1)
!!    besshl = bessjl(l)+SQRTm1*bessnl(l)
!!    Ulp1 = (besshl**2+besshlm1**2+(2*l+1)*fint(l-1))/(2*l-1.0d0)
!!    write(6,'(a,2d16.8,2x,2d16.8)')'fint(l), Ulp1 = ',fint(l), Ulp1
!!    fint(l) = Ulp1
!! enddo
!  ===================================================================
!
   do l = 0, lmax
      fint(l) = fint(l)*rc**3
   enddo
!
   end subroutine IntegrateSphHankelSq_Ca
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine IntegrateSphHankelSq_Cs(lmax,rc,energy,fint)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: rc
!
   complex (kind=CmplxKind), intent(in) :: energy
   complex (kind=CmplxKind), intent(out) :: fint
   complex (kind=CmplxKind) :: besshl(0:lmax+1), x
   complex (kind=CmplxKind) :: Ul, Ulp1
!
!  x = rc*sqrt(energy)+SQRTm1*0.00001d0
   x = rc*sqrt(energy)
!  -------------------------------------------------------------------
   call SphericalHankel(lmax+1,x,besshl)
!  -------------------------------------------------------------------
   if (lmax == 0) then
      fint = (besshl(0)*besshl(1)/x-besshl(0)**2-besshl(1)**2)*HALF
   else
      fint = (besshl(lmax-1)*besshl(lmax+1)-besshl(lmax)**2)*HALF
   endif
!
!  ===================================================================
!  Another way of calculating fint is to use recursive relation.
!  Both ways have shown to give the same results.
!  ===================================================================
!! Ul = -SQRTm1*exp(SQRTm1*TWO*x)/(TWO*x**3)
!! write(6,'(a,2d16.8,2x,2d16.8)')'fint(0), Ul = ',fint(0), Ul
!! fint(0) = Ul
!! do l = 1, lmax
!!    besshlm1 = bessjl(l-1)+SQRTm1*bessnl(l-1)
!!    besshl = bessjl(l)+SQRTm1*bessnl(l)
!!    Ulp1 = (besshl**2+besshlm1**2+(2*l+1)*fint(l-1))/(2*l-1.0d0)
!!    write(6,'(a,2d16.8,2x,2d16.8)')'fint(l), Ulp1 = ',fint(l), Ulp1
!!    fint(l) = Ulp1
!! enddo
!  ===================================================================
!
   fint = fint*rc**3
!
   end subroutine IntegrateSphHankelSq_Cs
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine IntegrateSphHankelSq_Ra(lmax,rc,energy,fint)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: rc
!
   real (kind=RealKind), intent(in) :: energy
   real (kind=RealKind), intent(out) :: fint(0:lmax)
!
!  complex (kind=CmplxKind) :: bessjl(0:lmax+1), bessnl(0:lmax+1), x, besshl, besshlp1, besshlm1
   complex (kind=CmplxKind) :: besshl(0:lmax+1), x, ec
   complex (kind=CmplxKind) :: Ul, Ulp1
!
!  x = rc*sqrt(energy)+SQRTm1*0.00001d0
   ec = energy
   x = rc*sqrt(ec)
!  -------------------------------------------------------------------
   call SphericalHankel(lmax+1,x,besshl)
!  -------------------------------------------------------------------
   fint(0) = abs(besshl(0)*besshl(1)/x-besshl(0)**2-besshl(1)**2)*HALF
   do l = 1, lmax
!     fint(l) = ((2*l+1)*besshl*besshlp1/x-besshl**2-besshlp1**2)*HALF
      fint(l) = abs(besshl(l-1)*besshl(l+1)-besshl(l)**2)*HALF
   enddo
!
!  ===================================================================
!  Another way of calculating fint is to use recursive relation.
!  Both ways have shown to give the same results.
!  ===================================================================
!! Ul = -SQRTm1*exp(SQRTm1*TWO*x)/(TWO*x**3)
!! write(6,'(a,2d16.8,2x,2d16.8)')'fint(0), Ul = ',fint(0), Ul
!! fint(0) = Ul
!! do l = 1, lmax
!!    besshlm1 = bessjl(l-1)+SQRTm1*bessnl(l-1)
!!    besshl = bessjl(l)+SQRTm1*bessnl(l)
!!    Ulp1 = (besshl**2+besshlm1**2+(2*l+1)*fint(l-1))/(2*l-1.0d0)
!!    write(6,'(a,2d16.8,2x,2d16.8)')'fint(l), Ulp1 = ',fint(l), Ulp1
!!    fint(l) = Ulp1
!! enddo
!  ===================================================================
!
   do l = 0, lmax
      fint(l) = fint(l)*rc**3
   enddo
!
   end subroutine IntegrateSphHankelSq_Ra
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine IntegrateSphHankelSq_Rs(lmax,rc,energy,fint)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind), intent(in) :: rc
!
   real (kind=RealKind), intent(in) :: energy
   real (kind=RealKind), intent(out) :: fint
!
!  complex (kind=CmplxKind) :: bessjl(0:lmax+1), bessnl(0:lmax+1), x, besshl, besshlp1, besshlm1
   complex (kind=CmplxKind) :: besshl(0:lmax+1), x, ec
   complex (kind=CmplxKind) :: Ul, Ulp1
!
!  x = rc*sqrt(energy)+SQRTm1*0.00001d0
   ec = energy
   x = rc*sqrt(ec)
!  -------------------------------------------------------------------
   call SphericalHankel(lmax+1,x,besshl)
!  -------------------------------------------------------------------
   if (lmax == 0) then
      fint = abs(besshl(0)*besshl(1)/x-besshl(0)**2-besshl(1)**2)*HALF
   else
      fint = abs(besshl(lmax-1)*besshl(lmax+1)-besshl(lmax)**2)*HALF
   endif
!
!  ===================================================================
!  Another way of calculating fint is to use recursive relation.
!  Both ways have shown to give the same results.
!  ===================================================================
!! Ul = -SQRTm1*exp(SQRTm1*TWO*x)/(TWO*x**3)
!! write(6,'(a,2d16.8,2x,2d16.8)')'fint(0), Ul = ',fint(0), Ul
!! fint(0) = Ul
!! do l = 1, lmax
!!    besshlm1 = bessjl(l-1)+SQRTm1*bessnl(l-1)
!!    besshl = bessjl(l)+SQRTm1*bessnl(l)
!!    Ulp1 = (besshl**2+besshlm1**2+(2*l+1)*fint(l-1))/(2*l-1.0d0)
!!    write(6,'(a,2d16.8,2x,2d16.8)')'fint(l), Ulp1 = ',fint(l), Ulp1
!!    fint(l) = Ulp1
!! enddo
!  ===================================================================
!
   fint = fint*rc**3
!
   end subroutine IntegrateSphHankelSq_Rs
!  ===================================================================
end module BesselModule
