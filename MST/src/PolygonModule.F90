module PolygonModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler
!
public :: initPolygon,       &
          endPolygon,        &
          genGaussPoints,    &
          getNumGaussPoints, &
          getGaussPosi,      &
          getGaussWght,      &
          getNumTriangles,   &
          getArea,           &
          getNumCorners,     &
          getCorner,         &
          getNormal
!
private
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level 
   integer (kind=IntKind) :: SetupMethod 
   integer (kind=IntKind) :: NumCorners 
   integer (kind=IntKind) :: NumGaussPoints 
   integer (kind=IntKind) :: NumTriangles   ! no. of finest triangles 
   integer (kind=IntKind) :: NumRootTris    ! no. of root triangles
!
   real (kind=RealKind), allocatable :: Corners(:,:)
   real (kind=RealKind) :: SurfNorm(3)
   real (kind=RealKind) :: SurfArea
!
   Type TriGaussStruct
      integer (kind=IntKind) :: ng 
      real (kind=RealKind), pointer :: xi(:,:)
      real (kind=RealKind), pointer :: wi(:)
   end type TriGaussStruct
!
   type (TriGaussStruct), allocatable :: TriGauss(:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPolygon(method,numc,cnx,cny,cnz,istop,iprint)
!  ===================================================================
   use MathParamModule, only : ZERO, HALF, TEN2m4, TEN2m6, TEN2m8, TEN2m10
   implicit none
!
   character (len=*) :: istop
!
   integer (kind=IntKind), intent(in) :: method, numc, iprint
   integer (kind=IntKind) :: ic
!
   real (kind=RealKind), intent(in) :: cnx(numc), cny(numc), cnz(numc)
   real (kind=RealKind) :: vpm, v1(3), v2(3), v(3), px, py, pz
   real (kind=RealKind) :: tol1 = Half*TEN2m8
   real (kind=RealKind) :: tol2 = Half*TEN2m10
!
   if (numc < 3) then
      call ErrorHandler('initPolygon','numc < 3',numc)
   else if (numc == 3) then
      NumRootTris = 1
   else
      NumRootTris = numc
   endif
!
   stop_routine = istop
   print_level = iprint
!
   allocate(Corners(3,numc), TriGauss(NumRootTris))
   do ic = 1,numc
      Corners(1,ic)=cnx(ic)
      Corners(2,ic)=cny(ic)
      Corners(3,ic)=cnz(ic)
   enddo
   NumCorners = numc
   SetupMethod=method
!
!  ===================================================================
!  Calculate the surface normal vector and the surface area
!  Also check if corners are on the same plane
!  ===================================================================
   px=ZERO
   py=ZERO
   pz=ZERO
   do ic=2,numc-1
!     ================================================================
!     Calculate:
!
!           ->    ->       ->    ->
!         ( c  -  c  ) x ( c   - c  )
!            i     1        i+1   1
!
!     and store the result in v(1:3) ...........................
!     ================================================================
      v1(1)=cnx(ic)-cnx(1)
      v1(2)=cny(ic)-cny(1)
      v1(3)=cnz(ic)-cnz(1)
      v2(1)=cnx(ic+1)-cnx(1)
      v2(2)=cny(ic+1)-cny(1)
      v2(3)=cnz(ic+1)-cnz(1)
!     ----------------------------------------------------------------
      call vcross(v1,v2,v)
!     ----------------------------------------------------------------
!
      if (cnx(1)*v(1)+cny(1)*v(2)+cnz(1)*v(3) < ZERO) then
         v(1)=-v(1)
         v(2)=-v(2)
         v(3)=-v(3)
      endif
      vpm = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      if (vpm < tol1) then
!        -------------------------------------------------------------
         call ErrorHandler('initPolygon','At least 3 corners are on a line')
!        -------------------------------------------------------------
      else 
         px=px+v(1)
         py=py+v(2)
         pz=pz+v(3)
         v(1:3) = v(1:3)/vpm
         if (ic == 2) then
            SurfNorm(1:3) = v(1:3)
         else if (abs(v(1)-SurfNorm(1))+abs(v(2)-Surfnorm(2))+        &
                  abs(v(3)-SurfNorm(3)) > tol2) then
!           ----------------------------------------------------------
            call ErrorHandler('initPolygon','Bad corners')
!           ----------------------------------------------------------
         endif
      endif
   enddo
!  ===================================================================
!  Calculate the plane area using formular:
!
!                1     ->   ->     ->   ->           ->   ->
!        Area = --- | (c1 x c2) + (c2 x c3) + ... + (cn x c1) |
!                2
!
!  ===================================================================
   vpm = sqrt(px*px+py*py+pz*pz)
   SurfArea = HALF*vpm
!
   NumGaussPoints = 0
   NumTriangles = 0
!
   end subroutine initPolygon
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endPolygon
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: it
!
   deallocate(Corners)
!
   do it = 1, NumRootTris
      deallocate(TriGauss(it)%xi, TriGauss(it)%wi)
   enddo
   deallocate(TriGauss)
   NumGaussPoints = 0
   NumTriangles = 0
   NumRootTris = 0
   NumCorners = 0
!
   end subroutine endPolygon
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGaussPoints() result(ng)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: ng
!
   ng = NumGaussPoints
   end function getNumGaussPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussPosi(i) result(x)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: ic, jc
!
   real (kind=RealKind) :: x(3)
!
   if (i < 1 .or. i > NumGaussPoints) then
      call ErrorHandler('getGaussPosi','Invalid Gauss Point index',i)
   endif
!
   jc = 0
   do ic=1,NumRootTris
      if (i <= TriGauss(ic)%ng+jc) then
         x(1:3) = TriGauss(ic)%xi(1:3,i-jc)
         return
      else
         jc = jc+TriGauss(ic)%ng
      endif
   enddo
   call ErrorHandler('getGaussPosi','Unknown error')
   end function getGaussPosi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussWght(i) result(w)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: ic, jc
!
   real (kind=RealKind) :: w
!
   if (i < 1 .or. i > NumGaussPoints) then
      call ErrorHandler('getGaussWght','Invalid Gauss Point index',i)
   endif
!
   jc = 0
   do ic=1,NumRootTris
      if (i <= TriGauss(ic)%ng+jc) then
         w = TriGauss(ic)%wi(i-jc)
      else
         jc = jc+TriGauss(ic)%ng
      endif
   enddo
   end function getGaussWght
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumTriangles() result(nt)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: nt
!
   nt = NumTriangles
   end function getNumTriangles
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getArea() result(a)
!  ===================================================================
   implicit none
   real (kind=RealKind) :: a
!
   a = SurfArea
   end function getArea
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCorners() result(nc)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: nc
!
   nc = NumCorners
   end function getNumCorners
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCorner(i) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
!
   real (kind=RealKind) :: c(3)
!
   if (i < 1 .or. i > NumCorners) then
      call ErrorHandler('getCorners','Invalid corner index',i)
   endif
!
   c(1:3) = Corners(1:3,i)
   end function getCorner
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNormal() result(s)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: s(3)
!
   s(1:3)=SurfNorm(1:3)
!
   end function getNormal
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genGaussPoints(epsilon)
!  ===================================================================
   use MathParamModule, only : ZERO, HALF, ONE, THIRD, TEN2m6
   use TriangleModule, only : initTriangle, endTriangle
   use TriangleModule, only : genTriGaussPoints => genGaussPoints
   use TriangleModule, only : getTriArea => getArea
   use TriangleModule, only : getNumTriGaussPoints => getNumGaussPoints
   use TriangleModule, only : getNumFinestTris
   use TriangleModule, only : getTriGaussPosi => getGaussPosi
   use TriangleModule, only : getTriGaussWght => getGaussWght
!
   implicit   none
!
   character (len=14), parameter :: sname='genGaussPoints'
!
   integer (kind=IntKind) :: ic, jc, kc, ng, ig
!
   real (kind=Realkind), intent(in) :: epsilon
   real (kind=Realkind) :: TriCorner1(3), TriCorner2(3), TriCorner3(3)
   real (kind=Realkind) :: p(3)
   real (kind=Realkind) :: PolyArea
!
!  ===================================================================
!  Find the center of the polygon
!  ===================================================================
   if (NumRootTris > 1) then  ! Note: NumRootTris is either 1 or numc (> 3)
      p(1:3) = ZERO
      do ic=1,NumCorners
         p(1:3) = p(1:3)+Corners(1:3,ic)
      enddo
      TriCorner1(1:3)=p(1:3)/dble(NumCorners)
   else
      TriCorner1(1:3)=Corners(1:3,1)
   endif
!
!  ===================================================================
!  loop over the corners forming triangles
!  ===================================================================
   NumGaussPoints = 0
   NumTriangles = 0
   PolyArea = ZERO
   do ic=1,NumRootTris         ! no. of corner per boundary plane
      if (NumRootTris == 1) then
         jc = 2
         kc = 3
      else if (ic < NumCorners) then
         jc = ic
         kc = ic+1
      else
         jc = ic
         kc = 1
      endif
      TriCorner2(1:3)=Corners(1:3,jc)
      TriCorner3(1:3)=Corners(1:3,kc)
!     ================================================================
!     call initTriangle to setup a root triangle......................
!     ----------------------------------------------------------------
      call initTriangle(TriCorner1,TriCorner2,TriCorner3,SetupMethod)
!     ----------------------------------------------------------------
!
!     ================================================================
!     call genTriGaussPoints to divide the root triangle and generate
!     Gaussian points.................................................
!     epsilon tells when to stop division: area/TriArea <= epsilon
!     ----------------------------------------------------------------
      call genTriGaussPoints(epsilon)
!     ----------------------------------------------------------------
      ng = getNumTriGaussPoints()
      NumTriangles = NumTriangles + getNumFinestTris()
!     ----------------------------------------------------------------
      TriGauss(ic)%ng = ng
      if (ng < 1) then
         call ErrorHandler(sname,'ng < 1',ng)
      else
         allocate(TriGauss(ic)%xi(3,ng), TriGauss(ic)%wi(ng))
      endif
      do ig=1,ng
         TriGauss(ic)%xi(1:3,ig)=getTriGaussPosi(ig)
         TriGauss(ic)%wi(ig)=getTriGaussWght(ig)
         PolyArea = PolyArea + TriGauss(ic)%wi(ig)
      enddo
      NumGaussPoints = NumGaussPoints + ng
!
!     ================================================================
!     call endTriangle to clear up the space.........................
!     ----------------------------------------------------------------
      call endTriangle()
!     ----------------------------------------------------------------
   enddo 
!
   if (abs(PolyArea-SurfArea) > TEN2m6) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'PolyArea <> SurfArea',PolyArea,SurfArea)
!     ----------------------------------------------------------------
   endif
!
   if (print_level >= 1) then
      write(6,*)' '
      write(6,'('' Number of Gaussian Points   = '',i5)')NumGaussPoints
      write(6,'('' Number of Finest Triangles  = '',i5)')NumTriangles
   endif
!
!  ===================================================================
   if(stop_routine.eq.'genGaussPoints') then
      call StopHandler(sname)
   endif
!
   end subroutine genGaussPoints
!
end module PolygonModule
