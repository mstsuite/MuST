module SurfElementsModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initSurfElements,        &
          endSurfElements,         &
          genGaussPoints,          &
          getNumGaussPoints,       &
          getGaussPosi,            &
          getGaussWght,            &
          getGaussR,               &
          getGaussNvec,            &
          getGaussRvecDotNvec,     &
          genGaussSphHarmonics,    &
          getGaussYlm,             &
          getGaussGradYlmDotNvec,  &
          getNumDiffSurfR,         &
          getDiffSurfR,            & 
          getGaussR2SurfR,         & 
          getWronskianIntegration, &
          testWronskian,           &
          printSurfaceElements,    &
          getSurfAverage
!
   interface getWronskianIntegration
      module procedure getWronskianIntegration1, getWronskianIntegration2
   end interface
!
private
   type PlaneStruct
      integer (kind=IntKind) :: ng
      real (kind=RealKind) :: surf_norm(3)
      real (kind=RealKind), pointer :: position(:,:)
      real (kind=RealKind), pointer :: weight(:)
      real (kind=RealKind), pointer :: r(:)
      real (kind=RealKind), pointer :: ndotr(:)
      real (kind=RealKind) :: area
   end type PlaneStruct
!
   type PolyhedronSurfStruct
      integer (kind=IntKind) :: NumPlanes
      integer (kind=IntKind) :: NumGaussPoints
      integer (kind=IntKind) :: NumDiffR
      integer (kind=IntKind), pointer :: PlaneIndex(:)
      integer (kind=IntKind), pointer :: id_r2diff_r(:)
      real (kind=RealKind) :: area
      real (kind=RealKind) :: epsilon
      real (kind=RealKind), pointer :: position(:,:)
      real (kind=RealKind), pointer :: weight(:)
      real (kind=RealKind), pointer :: r(:)
      real (kind=RealKind), pointer :: ndotr(:)
      real (kind=RealKind), pointer :: surf_norm(:,:)
      real (kind=RealKind), pointer :: diff_r(:)
      complex (kind=CmplxKind), pointer :: Ylm(:,:)
      complex (kind=CmplxKind), pointer :: NdotGradYlm(:,:)
   end type PolyhedronSurfStruct
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: lmax = -1
   integer (kind=IntKind) :: jmax = 0
   integer (kind=IntKind) :: kmax = 0
!
   integer (kind=IntKind) :: NumPolyhedra
   integer (kind=IntKind) :: print_level
!
   type (PolyhedronSurfStruct), allocatable, target :: PolyhedronSurface(:)
   type (PlaneStruct), allocatable :: Plane(:)
!
   logical :: Initialized = .false.
   logical :: YlmGenerated = .false.
   logical :: GPGenerated = .false.
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSurfElements(istop,iprint)
!  ===================================================================
   use PolyhedraModule, only : getNumPolyhedra
   use PolyhedraModule, only : getNumPlanes
   use PolyhedraModule, only : getSurfaceArea
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: ip
!
   real (kind=RealKind) :: delta = 0.0004d0
!
   if ( Initialized .eqv. .true. ) return
!
   NumPolyhedra = getNumPolyhedra()
!
   if (NumPolyhedra < 1) then
      call ErrorHandler('initSurfElements','NumPolyhedra < 1',NumPolyhedra)
   endif
!
   allocate( PolyhedronSurface(NumPolyhedra) )
   stop_routine = istop
   print_level = iprint
!
   do ip=1,NumPolyhedra
      PolyhedronSurface(ip)%NumGaussPoints = 0
      PolyhedronSurface(ip)%NumPlanes = getNumPlanes(ip)
      PolyhedronSurface(ip)%area = getSurfaceArea(ip)
      PolyhedronSurface(ip)%epsilon = PolyhedronSurface(ip)%area*delta
      nullify( PolyhedronSurface(ip)%PlaneIndex )
      nullify( PolyhedronSurface(ip)%position )
      nullify( PolyhedronSurface(ip)%weight )
      nullify( PolyhedronSurface(ip)%r )
      nullify( PolyhedronSurface(ip)%ndotr )
      nullify( PolyhedronSurface(ip)%surf_norm )
      nullify( PolyhedronSurface(ip)%Ylm )
      nullify( PolyhedronSurface(ip)%NdotGradYlm )
   enddo
!
   Initialized = .true.
   YlmGenerated = .false.
   GPGenerated = .false.
!
   end subroutine initSurfElements
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSurfElements()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i
!
   if ( .not.Initialized ) return
!
   do i=1,NumPolyhedra
      if (GPGenerated) then
         deallocate( PolyhedronSurface(i)%position )
         deallocate( PolyhedronSurface(i)%weight )
         deallocate( PolyhedronSurface(i)%r )
         deallocate( PolyhedronSurface(i)%ndotr )
         deallocate( PolyhedronSurface(i)%PlaneIndex )
         deallocate( PolyhedronSurface(i)%surf_norm )
         deallocate( PolyhedronSurface(i)%id_r2diff_r )
         deallocate( PolyhedronSurface(i)%diff_r )
      endif
      if (YlmGenerated) then
         deallocate( PolyhedronSurface(i)%Ylm )
         deallocate( PolyhedronSurface(i)%NdotGradYlm )
      endif
   enddo
   deallocate( PolyhedronSurface )
!
   Initialized = .false.
   YlmGenerated = .false.
   GPGenerated = .false.
!
   end subroutine endSurfElements
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genGaussPoints()
!  ===================================================================
   use MathParamModule, only : ZERO, TEN2m8, ONE
!
   use SortModule, only : HeapSort
!
   use PolyhedraModule, only : getNumPlaneCorners
   use PolyhedraModule, only : getPlaneCorner
   use PolyhedraModule, only : getPlaneArea
   use PolyhedraModule, only : getPlaneNormal
   use PolyhedraModule, only : getPlane
!
   use PolygonModule, only : initPolygon, endPolygon
   use PolygonModule, only : genPolygonGaussPoints => genGaussPoints
   use PolygonModule, only : getNumPolygonGaussPoints => getNumGaussPoints
   use PolygonModule, only : getPolygonGaussPosi => getGaussPosi
   use PolygonModule, only : getPolygonGaussWght => getGaussWght
   use PolygonModule, only : getArea
   use PolygonModule, only : getNormal
!
   implicit none
!
   integer (kind=IntKind) :: method = 5
!
   integer (kind=IntKind) :: i, j, k, ng, np, numc, ntot, id, poly, nr
   integer (kind=IntKind), allocatable :: indx(:)
!
   real (kind=RealKind) :: r, efac
   real (kind=RealKind) :: c(3), dif(3), posi(3)
   real (kind=RealKind), allocatable :: cnx(:), cny(:), cnz(:)
   real (kind=RealKind), allocatable :: rtmp(:)
   real (kind=RealKind), pointer :: pp(:,:)
!
   if (.not.Initialized) then
       call ErrorHandler('genGaussPoints',                           &
                         'Need to initialize SurfElementsModule')
   else if (GPGenerated) then
       return
   endif
!
   do i=1,NumPolyhedra
      pp => getPlane(i)
      np = PolyhedronSurface(i)%NumPlanes
      allocate( Plane(np) )
      ntot = 0
      do j=1,np
         Plane(j)%surf_norm(1:3) = getPlaneNormal(i,j)
         Plane(j)%area = getPlaneArea(i,j)
         numc = getNumPlaneCorners(i,j)
         allocate( cnx(numc), cny(numc), cnz(numc) )
         do k=1,numc
            c(1:3) = getPlaneCorner(i,j,k)
            cnx(k) = c(1)
            cny(k) = c(2)
            cnz(k) = c(3)
         enddo
!        -------------------------------------------------------------
         call initPolygon(method,numc,cnx,cny,cnz,stop_routine,print_level)
!        -------------------------------------------------------------
         if (abs(getArea()-Plane(j)%area) > TEN2m8) then
            call ErrorHandler('genGaussPoints','Invalid PlaneArea',   &
                              getArea(),Plane(j)%area)
         else 
            dif(1:3) = getNormal()-Plane(j)%surf_norm(1:3)
            if (abs(dif(1))+abs(dif(2))+abs(dif(3)) > TEN2m8) then
               print *,'Normal from PolygonModule: ',getNormal()
               print *,'Normal from PolyhedraModule: ', Plane(j)%surf_norm(1:3)
               call ErrorHandler('genGaussPoints','Invalid PlaneNormal')
            endif
         endif
!
         r = sqrt(pp(1,j)*pp(1,j)+pp(2,j)*pp(2,j)+pp(3,j)*pp(3,j))
         efac = ( pp(1,j)*Plane(j)%surf_norm(1) +                     &
                  pp(2,j)*Plane(j)%surf_norm(2) +                     &
                  pp(3,j)*Plane(j)%surf_norm(3) )/r
         if (efac <= ZERO) then
            call ErrorHandler('genGaussPoints','efac <= 0',efac)
         else if (efac > ONE+TEN2m8) then
            call ErrorHandler('genGaussPoints','efac > 1',efac)
         endif
!
!        -------------------------------------------------------------
         call genPolygonGaussPoints(PolyhedronSurface(i)%epsilon*efac)
!        -------------------------------------------------------------
         ng = getNumPolygonGaussPoints()
         ntot = ntot + ng
         Plane(j)%ng = ng
         allocate( Plane(j)%position(1:3,ng) )
         allocate( Plane(j)%weight(ng) )
         allocate( Plane(j)%r(ng) )
         allocate( Plane(j)%ndotr(ng) )
         do k=1,ng
            posi(1:3) = getPolygonGaussPosi(k)
            Plane(j)%position(1:3,k) = posi(1:3)
            Plane(j)%weight(k) = getPolygonGaussWght(k)
            Plane(j)%r(k)=sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
            Plane(j)%ndotr(k) = (posi(1)*Plane(j)%surf_norm(1)+ &
                                 posi(2)*Plane(j)%surf_norm(2)+ &
                                 posi(3)*Plane(j)%surf_norm(3))/Plane(j)%r(k)
         enddo
!        -------------------------------------------------------------
         call endPolygon()
!        -------------------------------------------------------------
         deallocate( cnx, cny, cnz )
      enddo
      PolyhedronSurface(i)%NumGaussPoints = ntot
!
      allocate( PolyhedronSurface(i)%position(3,ntot) )
      allocate( PolyhedronSurface(i)%weight(ntot) )
      allocate( PolyhedronSurface(i)%r(ntot) )
      allocate( PolyhedronSurface(i)%ndotr(ntot) )
      allocate( PolyhedronSurface(i)%PlaneIndex(ntot) )
      allocate( PolyhedronSurface(i)%surf_norm(3,np) )
!
      id = 0
      do j=1,np
         PolyhedronSurface(i)%surf_norm(1:3,j)=Plane(j)%surf_norm(1:3)
         do k=1,Plane(j)%ng
            id = id + 1
            PolyhedronSurface(i)%position(1:3,id) = Plane(j)%position(1:3,k)
            PolyhedronSurface(i)%weight(id) = Plane(j)%weight(k)
            PolyhedronSurface(i)%r(id) = Plane(j)%r(k)
            PolyhedronSurface(i)%ndotr(id) = Plane(j)%ndotr(k)
            PolyhedronSurface(i)%PlaneIndex(id) = j
         enddo
         deallocate( Plane(j)%position )
         deallocate( Plane(j)%weight )
         deallocate( Plane(j)%r )
         deallocate( Plane(j)%ndotr )
      enddo
      deallocate( Plane )
   enddo
!
   do poly=1,NumPolyhedra
      ng=PolyhedronSurface(poly)%NumGaussPoints
      allocate(rtmp(1:ng), indx(1:ng))
!
      rtmp(1:ng)=PolyhedronSurface(poly)%r(1:ng)
!     ----------------------------------------------------------------
      call HeapSort(ng,rtmp,indx)
!     ----------------------------------------------------------------
      nr=1
      r=rtmp(1)
      do i=2,ng
         if (rtmp(i)-r > TEN2m8) then
             nr=nr+1
             r=rtmp(i)
         endif
      enddo
      PolyhedronSurface(poly)%NumDiffR=nr
!
      allocate( PolyhedronSurface(poly)%diff_r(1:nr) )
      allocate( PolyhedronSurface(poly)%id_r2diff_r(1:ng) )
!
      PolyhedronSurface(poly)%diff_r(1)=rtmp(1)
      k=indx(1)
      PolyhedronSurface(poly)%id_r2diff_r(k)=1
      j=1
      do i=2, ng
         if (rtmp(i)-PolyhedronSurface(poly)%diff_r(j) > TEN2m8) then
             j=j+1
             PolyhedronSurface(poly)%diff_r(j)=rtmp(i)
         endif
         k=indx(i)
         PolyhedronSurface(poly)%id_r2diff_r(k)=j
      enddo
      deallocate(rtmp, indx)
   enddo
!
   GPGenerated = .true.
!
   end subroutine genGaussPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGaussPoints(poly) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: n
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getNumGaussPoints','Invalid polyhedron index',poly)
   endif
!
   n = PolyhedronSurface(poly)%NumGaussPoints
!
   end function getNumGaussPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussPosi(poly) result(pos)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
!
   real (kind=RealKind), pointer :: pos(:,:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussPosi','Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   pos => PolyhedronSurface(poly)%position(1:3,1:ng)
!
   end function getGaussPosi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussWght(poly) result(w)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
!
   real (kind=RealKind), pointer :: w(:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussWght','Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   w => PolyhedronSurface(poly)%weight(1:ng)
!
   end function getGaussWght
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussR(poly) result(r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
!
   real (kind=RealKind), pointer :: r(:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussR','Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   r => PolyhedronSurface(poly)%r(1:ng)
!
   end function getGaussR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussNvec(poly) result(nvec)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
!
   real (kind=RealKind), pointer :: nvec(:,:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussNvec','Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   nvec => PolyhedronSurface(poly)%surf_norm(1:3,1:ng)
!
   end function getGaussNvec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussRvecDotNvec(poly) result(r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
!
   real (kind=RealKind), pointer :: r(:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussRvecDotNvec','Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   r => PolyhedronSurface(poly)%ndotr(1:ng)
!
   end function getGaussRvecDotNvec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genGaussSphHarmonics(l)
!  ===================================================================
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: l
   integer (kind=IntKind) :: poly, ig, ip, kl, ng
!
   real (kind=RealKind) :: x, y, z
!
   complex (kind=CmplxKind), allocatable :: ylm(:)
   complex (kind=CmplxKind), allocatable :: grady(:,:)
!
   if (.not.GPGenerated) then
      call ErrorHandler('genGaussSphHarmonics','call genGaussPoints first')
   else if (l < 0) then
      call ErrorHandler('genGaussSphHarmonics','l < 0',l)
   endif
!
   if ( lmax >= l ) then
      return
   else
      lmax = l
      jmax = (l+1)*(l+2)/2
      kmax = (l+1)*(l+1)
      allocate( ylm(kmax), grady(kmax,3) )
   endif
!
   do poly=1,NumPolyhedra
      ng = PolyhedronSurface(poly)%NumGaussPoints
      if ( associated(PolyhedronSurface(poly)%Ylm) ) then
         deallocate( PolyhedronSurface(poly)%Ylm )
      endif
      if ( associated(PolyhedronSurface(poly)%NdotGradYlm) ) then
         deallocate( PolyhedronSurface(poly)%NdotGradYlm )
      endif
      allocate( PolyhedronSurface(poly)%Ylm(kmax,ng) )
      allocate( PolyhedronSurface(poly)%NdotGradYlm(kmax,ng) )
      do ig=1,ng
         x=PolyhedronSurface(poly)%position(1,ig)
         y=PolyhedronSurface(poly)%position(2,ig)
         z=PolyhedronSurface(poly)%position(3,ig)
!        -------------------------------------------------------------
         call calYlm(x,y,z,l,ylm,grady)
!        -------------------------------------------------------------
         ip = PolyhedronSurface(poly)%PlaneIndex(ig)
         do kl=1,kmax
            PolyhedronSurface(poly)%Ylm(kl,ig) = ylm(kl)
            PolyhedronSurface(poly)%NdotGradYlm(kl,ig) =              &
                 PolyhedronSurface(poly)%surf_norm(1,ip)*grady(kl,1)+ &
                 PolyhedronSurface(poly)%surf_norm(2,ip)*grady(kl,2)+ &
                 PolyhedronSurface(poly)%surf_norm(3,ip)*grady(kl,3)
         enddo
      enddo
   enddo
!
   deallocate( ylm, grady )
   YlmGenerated = .true.
!
   end subroutine genGaussSphHarmonics
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussYlm(poly) result(ylm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
!
   complex (kind=CmplxKind), pointer :: ylm(:,:)
!
   if(.not.YlmGenerated) then
      call ErrorHandler('getGaussYlm','call genGaussSphHarmonics first')
   else if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussYlm','Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   ylm => PolyhedronSurface(poly)%ylm(1:kmax,1:ng)
!
   end function getGaussYlm
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussGradYlmDotNvec(poly) result(ndylm)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
!
   complex (kind=CmplxKind), pointer :: ndylm(:,:)
!
   if(.not.YlmGenerated) then
      call ErrorHandler('getGaussGradYlmDotNvec',                     &
                        'call genGaussSphHarmonics first')
   else if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussGradYlmDotNvec',                     &
                        'Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   ndylm => PolyhedronSurface(poly)%NdotGradYlm(1:kmax,1:ng)
!
   end function getGaussGradYlmDotNvec
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumDiffSurfR(poly) result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: n
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getNumDiffSurfR','Invalid polyhedron index',poly)
   endif
!
   n = PolyhedronSurface(poly)%NumDiffR
!
   end function getNumDiffSurfR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDiffSurfR(poly) result(r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: n
!
   real (kind=RealKind), pointer :: r(:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getDiffSurfR','Invalid polyhedron index',poly)
   endif
!
   n = PolyhedronSurface(poly)%NumDiffR
   r => PolyhedronSurface(poly)%diff_r(1:n)
!
   end function getDiffSurfR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussR2SurfR(poly) result(i2j)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind) :: ng
   integer (kind=IntKind), pointer :: i2j(:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getGaussR2SurfR','Invalid polyhedron index',poly)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
   i2j => PolyhedronSurface(poly)%id_r2diff_r(1:ng)
!
   end function getGaussR2SurfR
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSurfAverage(poly,nr,lmax_in,r_mesh,f) result(av)
!  ===================================================================
   use MathParamModule, only : ZERO, TWO
   use InterpolationModule, only : PolyInterp, FitInterp
!
   implicit none
!
   integer, intent(in) :: poly, nr, lmax_in
   integer :: l, m, jl, kl, ng, ig, jh
   integer (kind=IntKind), parameter :: n_inter=5
!
   real (kind=RealKind), target, intent(in) :: r_mesh(1:nr)
   real (kind=RealKind) :: av, err_est
   real (kind=RealKind), pointer :: rg(:)
   real (kind=RealKind), pointer :: wg(:)
   real (kind=RealKind), pointer :: xp(:)
!
   complex (kind=CmplxKind), target, intent(in) :: f(:,:)
   complex (kind=CmplxKind), pointer :: ylm(:,:)
   complex (kind=CmplxKind), pointer :: yp(:)
   complex (kind=CmplxKind) :: fp
!
   interface
      subroutine hunt(n,xx,x,jlo)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(inout) :: jlo
         real (kind=RealKind), intent(in) :: xx(n)
         real (kind=RealKind), intent(in) :: x
      end subroutine hunt
   end interface
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('getSurfAverage','Invalid polyhedron index',poly)
   else if (lmax_in > lmax .or. lmax_in < 0) then
      call ErrorHandler('getSurfAverage','Input lmax exceeds the limit',lmax_in)
   else if (.not.YlmGenerated) then
      call ErrorHandler('getSurfAverage','Need to call genGaussSphHarmonics first')
   else if (.not.GPGenerated) then
      call ErrorHandler('getSurfAverage','Need to call genGaussSphHarmonics first')
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
!
   ylm => PolyhedronSurface(poly)%Ylm(1:kmax,1:ng)
   rg => PolyhedronSurface(poly)%r(1:ng)
   wg => PolyhedronSurface(poly)%weight(1:ng)
!
   jh = nr/2
   av = ZERO
   do ig = 1, ng
!     ----------------------------------------------------------------
      call hunt(nr, r_mesh, rg(ig), jh)
!     ----------------------------------------------------------------
      if (jh > nr-(n_inter-1)/2) then
         jh=nr-n_inter+1
      else if (2*jh+1 > n_inter) then
         jh=jh-(n_inter-1)/2
      else
         jh=1
      endif
      xp => r_mesh(jh:jh+n_inter-1)
      jl = 0
      do l = 0, lmax_in
         kl = (l+1)*l
         do m = 0, l
            jl = jl + 1
            kl = kl + 1
            yp => f(jh:jh+n_inter-1,jl)
!           ----------------------------------------------------------
            call PolyInterp(n_inter,xp,yp,rg(ig),fp,err_est)
!           ----------------------------------------------------------
            if (m == 0) then
               av = av + wg(ig)*real(fp*ylm(kl,ig),kind=RealKind)
            else
               av = av + wg(ig)*TWO*real(fp*ylm(kl,ig),kind=RealKind)
            endif
         enddo
      enddo
   enddo
   av = av/PolyhedronSurface(poly)%area
!
   nullify( ylm, rg, wg, xp, yp )
!
   end function getSurfAverage
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine testWronskian(poly,ltest,kappa)
!  ===================================================================
   use MathParamModule, only : CZERO, TEN2m6, SQRTm1, CONE, TEN2m8, TEN2m16
   use MathParamModule, only : HALF
   use BesselModule, only : SphericalBessel, SphericalNeumann
   use PolyhedraModule, only : getInscrSphRadius
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: poly
   integer (kind=IntKind), intent(in) :: ltest
   integer (kind=IntKind) :: ig, ng, kl, klp, l, lp, ktest
!
   real (kind=RealKind), pointer :: n_dot_r(:)
   real (kind=RealKind), pointer :: r(:)
   real (kind=RealKind), pointer :: wi(:)
   real (kind=RealKind) :: rmt
!
   complex (kind=CmplxKind), intent(in) :: kappa
!
   complex (kind=CmplxKind) :: ydy, jdn, ndn, jdj, h1dn, h2dn
   complex (kind=CmplxKind), pointer :: nDotGradY(:,:), ylm(:,:)
   complex (kind=CmplxKind), allocatable :: sjl(:), sdjl(:), snl(:), sdnl(:)
   complex (kind=CmplxKind), allocatable :: jl(:,:), djl(:,:), nl(:,:), dnl(:,:)
   complex (kind=CmplxKind), allocatable :: h1l(:,:),dh1l(:,:)
   complex (kind=CmplxKind), allocatable :: h2l(:,:),dh2l(:,:)
   complex (kind=CmplxKind) :: wronski_jn, wronski_nn, wronski_jj
   complex (kind=CmplxKind) :: wronski_h1n, wronski_h2n, h1ph2, h1mh2
   complex (kind=CmplxKind) :: a1, a2, b1, b2, a, b
   complex (kind=CmplxKind) :: nl1, nlp1, nfac
   complex (kind=CmplxKind), allocatable :: nlrmt(:), dnlrmt(:)
!
   if (poly < 1 .or. poly > NumPolyhedra) then
      call ErrorHandler('testWronskian','Invalid polyhedron index',poly)
   else if (ltest < 0) then
      call ErrorHandler('testWronskian','l < 0',ltest)
   endif
!
   ktest = (ltest+1)*(ltest+1)
!
   if (lmax < ltest) then
      call genGaussSphHarmonics(ltest)
   endif
!
   ng = PolyhedronSurface(poly)%NumGaussPoints
!
   allocate( sjl(0:ltest), sdjl(0:ltest), jl(1:ng,0:ltest), djl(1:ng,0:ltest) )
   allocate( snl(0:ltest), sdnl(0:ltest), nl(1:ng,0:ltest), dnl(1:ng,0:ltest) )
   allocate( h1l(1:ng,0:ltest), dh1l(1:ng,0:ltest) )
   allocate( h2l(1:ng,0:ltest), dh2l(1:ng,0:ltest) )
!
   rmt = getInscrSphRadius(poly)
   allocate( nlrmt(0:ltest), dnlrmt(0:ltest) )
!  -------------------------------------------------------------------
   call SphericalNeumann(ltest,kappa*rmt,nlrmt,dnlrmt)
!  -------------------------------------------------------------------
!
   r => PolyhedronSurface(poly)%r(1:ng)
   do ig=1,ng
!     ----------------------------------------------------------------
      call SphericalBessel(ltest,kappa*r(ig),sjl,sdjl)
!     ----------------------------------------------------------------
      call SphericalNeumann(ltest,kappa*r(ig),snl,sdnl)
!     ----------------------------------------------------------------
      do l=0,ltest
         jl(ig,l) = sjl(l)
         djl(ig,l) = sdjl(l)
         nl(ig,l) = snl(l)
         dnl(ig,l) = sdnl(l)
         wronski_jn=r(ig)*r(ig)*kappa*kappa*(jl(ig,l)*dnl(ig,l)-      &
                                             djl(ig,l)*nl(ig,l))
         h1l(ig,l) = sjl(l)+SQRTm1*snl(l)
         dh1l(ig,l) = sdjl(l)+SQRTm1*sdnl(l)
         h2l(ig,l) = sjl(l)-SQRTm1*snl(l)
         dh2l(ig,l) = sdjl(l)-SQRTm1*sdnl(l)
      enddo
   enddo
   deallocate( sjl, sdjl, snl, sdnl )
!
   ylm => PolyhedronSurface(poly)%Ylm(1:kmax,1:ng)
   nDotGradY => PolyhedronSurface(poly)%NdotGradYlm(1:kmax,1:ng)
   n_dot_r => PolyhedronSurface(poly)%ndotr(1:ng)
   wi => PolyhedronSurface(poly)%weight(1:ng)
!
   do klp=1,ktest
      lp = sqrt(klp-1.0d0)
      nlp1 = dnlrmt(lp)
      do kl=1,ktest
         l = sqrt(kl-1.0d0)
         nl1 = dnlrmt(l)
         nfac = nlp1*nl1
         wronski_jn=czero
         wronski_jj=czero
         wronski_nn=czero
         wronski_h1n=czero
         wronski_h2n=czero
         a = czero; b = czero
         do ig=1,ng
            jdn=kappa*wi(ig)*n_dot_r(ig)*( jl(ig,lp)*dnl(ig,l)-       &
                                           djl(ig,lp)*nl(ig,l) )
            jdj=kappa*wi(ig)*n_dot_r(ig)*( jl(ig,lp)*djl(ig,l)-       &
                                           djl(ig,lp)*jl(ig,l) )
            ndn=kappa*wi(ig)*n_dot_r(ig)*( nl(ig,lp)*dnl(ig,l)/nfac-  &
                                           dnl(ig,lp)*nl(ig,l)/nfac )
            h1dn=kappa*wi(ig)*n_dot_r(ig)*( h1l(ig,lp)*dnl(ig,l)-     &
                                            dh1l(ig,lp)*nl(ig,l) )
            h2dn=kappa*wi(ig)*n_dot_r(ig)*( h2l(ig,lp)*dnl(ig,l)-     &
                                            dh2l(ig,lp)*nl(ig,l) )
            ydy=wi(ig)*( conjg(ylm(klp,ig))*nDotGradY(kl,ig)-         &
                         conjg(nDotGradY(klp,ig))*ylm(kl,ig) )
            wronski_jn=wronski_jn+(jl(ig,lp)*nl(ig,l)*ydy+            &
                                   conjg(ylm(klp,ig))*ylm(kl,ig)*jdn)/nl1
            wronski_nn=wronski_nn+(nl(ig,lp)*nl(ig,l)/nfac*ydy+       &
                                   conjg(ylm(klp,ig))*ylm(kl,ig)*ndn)
            wronski_jj=wronski_jj+jl(ig,lp)*jl(ig,l)*ydy+             &
                                  conjg(ylm(klp,ig))*ylm(kl,ig)*jdj
            wronski_h1n=wronski_h1n+(h1l(ig,lp)*nl(ig,l)*ydy+         &
                                     conjg(ylm(klp,ig))*ylm(kl,ig)*h1dn)/nl1
            wronski_h2n=wronski_h2n+(h2l(ig,lp)*nl(ig,l)*ydy+         &
                                     conjg(ylm(klp,ig))*ylm(kl,ig)*h2dn)/nl1
            a1 = kappa*wi(ig)*n_dot_r(ig)*jl(ig,lp)*dnl(ig,l)
            a2 = wi(ig)*conjg(ylm(klp,ig))*nDotGradY(kl,ig)
            b1 = kappa*wi(ig)*n_dot_r(ig)*djl(ig,lp)*nl(ig,l)
            b2 = wi(ig)*conjg(nDotGradY(klp,ig))*ylm(kl,ig)
            a=a+jl(ig,lp)*nl(ig,l)*a2/nl1+conjg(ylm(klp,ig))*ylm(kl,ig)*a1/nl1
            b=b+jl(ig,lp)*nl(ig,l)*b2/nl1+conjg(ylm(klp,ig))*ylm(kl,ig)*b1/nl1
         enddo
         wronski_jn=kappa*wronski_jn*nl1
         wronski_jj=kappa*wronski_jj
         wronski_nn=kappa*wronski_nn*nfac
         h1ph2 = kappa*(wronski_h1n+wronski_h2n)*nl1*HALF
         h1mh2 = SQRTm1*kappa*(wronski_h2n-wronski_h1n)*nl1*HALF
         if (abs(wronski_jn) > TEN2m6) then
            print *,'kl = ',kl,', klp = ',klp,', wronski_jn = ',wronski_jn
         endif
         if (abs(wronski_jj) > TEN2m6) then
            print *,'kl = ',kl,', klp = ',klp,', wronski_jj = ',wronski_jj
         endif
         if (abs(wronski_nn) > TEN2m6) then
            print *,'kl = ',kl,', klp = ',klp,', wronski_nn = ',wronski_nn
         endif
         if (abs(h1ph2) > TEN2m6) then
            print *,'kl = ',kl,', klp = ',klp,', wronski_(h1n+h2n)/2 = ', h1ph2
         endif
         if (abs(h1mh2) > TEN2m6) then
            print *,'kl = ',kl,', klp = ',klp,', wronski_(h1n-h2n)/2 = ', h1mh2
         endif
         if (abs(kappa*(a-b)*nl1) > TEN2m6) then
            print *,'kl = ',kl,', klp = ',klp,', a-b = ',kappa*(a-b)*nl1
         endif
      enddo
   enddo
!
   deallocate( nlrmt, dnlrmt )
   deallocate( jl, djl, nl, dnl, h1l, dh1l, h2l, dh2l )
   nullify( n_dot_r, r, wi, nDotGradY, ylm )
!
   end subroutine testWronskian
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWronskianIntegration1(nr,kmax,fr,rdf,gr,rdg) result(w)
!  ===================================================================
!
!  calculate the surface integration of the wronskian:
!         [f^*, g] = f^{*) grad[g] - g grad[f^{*}]
!
!  the input function fr and gr are functions f and g multiplied by r, 
!  and are expanded in Y_{l,m} as follows:
!         r*f(vec{r}) = r * sum_{l,m} f(r,{l,m})*Y_{l,m}
!         r*g(vec{r}) = r * sum_{l,m} g(r,{l,m})*Y_{l,m}
!  where {l,m} is cut off at kmax.
!
!  the input function rdf and rdg are the radial derivative of f and g
!  multiplied by r.
!  *******************************************************************
   use MathParamModule, only : CZERO
   implicit none
!
   integer, intent(in) :: nr
   integer, intent(in) :: kmax
!
   complex (kind=CmplxKind), intent(in) :: fr(nr,kmax)
   complex (kind=CmplxKind), intent(in) :: rdf(nr,kmax)
   complex (kind=CmplxKind), intent(in) :: gr(nr,kmax)
   complex (kind=CmplxKind), intent(in) :: rdg(nr,kmax)
   complex (kind=CmplxKind) :: w
!
   w = CZERO
!
   end function getWronskianIntegration1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWronskianIntegration2(nr,kl,fr,rdf,kmax,gr,rdg) result(w)
!  ===================================================================
!
!  calculate the surface integration of the wronskian:
!         [f^*, g] = f^{*) grad[g] - g grad[f^{*}]
!
!  the input function fr and gr are functions f and g multiplied by r, 
!  and are expanded in Y_{l,m} as follows:
!         r*f(vec{r}) = r * f(r)*Y_{kl}
!         r*g(vec{r}) = r * sum_{l,m} g(r,{l,m})*Y_{l,m}
!  where {l,m} is cut off at kmax.
!
!  the input function rdf and rdg are the radial derivative of f and g
!  multiplied by r.
!  *******************************************************************
   use MathParamModule, only : CZERO
   implicit none
!
   integer, intent(in) :: nr
   integer, intent(in) :: kl
   integer, intent(in) :: kmax
!
   complex (kind=CmplxKind), intent(in) :: fr(nr)
   complex (kind=CmplxKind), intent(in) :: rdf(nr)
   complex (kind=CmplxKind), intent(in) :: gr(nr,kmax)
   complex (kind=CmplxKind), intent(in) :: rdg(nr,kmax)
   complex (kind=CmplxKind) :: w
!
   w = CZERO
!
   end function getWronskianIntegration2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printSurfaceElements( id )
!  ===================================================================
   use MPPModule, only: MyPE
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
   character(len=30) :: file_wfl
   integer (kind=IntKind) ::  ir, NumRs, plane
   integer (kind=IntKind) ::  funit, nc
   integer (kind=IntKind) :: offset = 100000
   integer (kind=IntKind) :: offset_at = 1000
!
   real (kind=RealKind), pointer :: r_mesh(:), posi(:,:), weight(:)
!
   type (PolyhedronSurfStruct), pointer :: PSS
!
   file_wfl(1:30) = " "
   nc = 5
   file_wfl(1:3)  = "SPT"
   file_wfl(4:nc)  = "S_"
   write(file_wfl(nc+1:nc+6),'(i6)') offset+MyPE
   file_wfl(nc+1:nc+1) = 'n'
   file_wfl(nc+7:nc+7) = '_'
   write(file_wfl(nc+8:nc+11),'(i4)') offset_at+id
   file_wfl(nc+8:nc+8) = 'a'
   funit = 55+id
   open(unit=funit,file=trim(file_wfl),status='unknown')
!
   PSS => PolyhedronSurface(id)
   write(funit,'(a)') "# Ind       r, posi, weight"
!
   NumRs = PSS%NumGaussPoints
   r_mesh => PSS%r
   posi => PSS%position(:,:)
   weight => PSS%weight(:)
   do ir = 1, NumRs
      plane = PSS%PlaneIndex(ir)
      write(funit,'(i4,5(1x,d16.8),1x,i5)') id, r_mesh(ir), posi(1:3,ir), weight(ir), plane
   enddo
   write(funit,'(a)') "  "
!
   close(funit)
!
   end subroutine printSurfaceElements
!  ===================================================================
end module SurfElementsModule
