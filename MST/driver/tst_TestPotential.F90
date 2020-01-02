program tst_TestPotential
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : THIRD, ONE, TWO, PI, PI4, CZERO, ZERO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius
   use PolyhedraModule, only : getOutscrSphRadius
   use PolyhedraModule, only : printPolyhedron
!
   use RadialGridModule, only : initRadialGrid, endRadialGrid
   use RadialGridModule, only : genRadialGrid, printRadialGrid
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType
!
   use TestPotentialModule, only : initTestPotential,     &
                                   endTestPotential,      &
                                   readTestPotential,     &
                                   testPotential
!
   implicit none
!
   integer (kind=IntKind) :: NumAtoms, ptype
   integer (kind=IntKind) :: i
   integer (kind=IntKind), allocatable :: lmax(:)
!
   integer (kind=IntKind), parameter :: ndivin = 101
   integer (kind=IntKind), parameter :: ndivout = 50
   integer (kind=IntKind), parameter :: nmult = 4
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: AtomRadius(:)
   real (kind=RealKind), allocatable :: x(:), y(:), z(:), v0(:), a(:,:), u(:,:)
!
   read(5,*)bravais(1:3,1)
   read(5,*)bravais(1:3,2)
   read(5,*)bravais(1:3,3)
   read(5,*)NumAtoms, ptype
!
   if (NumAtoms < 1) then
      call ErrorHandler('main','invalid NumAtoms',NumAtoms)
   endif
!
!  ===================================================================
!  initialize radial grid
!  -------------------------------------------------------------------
   call initRadialGrid(NumAtoms, 'main', 0)
!  -------------------------------------------------------------------
   call initPolyhedra(NumAtoms,bravais,'main',0)
!  -------------------------------------------------------------------
   call initPotentialType(ptype)
!  -------------------------------------------------------------------
   call initTestPotential(NumAtoms,1)
!  -------------------------------------------------------------------
!
   allocate(AtomPosition(3,NumAtoms), AtomRadius(NumAtoms))
   allocate(lmax(NumAtoms),x(NumAtoms),y(NumAtoms),z(NumAtoms))
!
   do i=1,NumAtoms
      read(5,*)AtomPosition(1:3,i),AtomRadius(i)
      read(5,*)lmax(i)
      read(5,*)x(i),y(i),z(i)
   enddo
!
   call readTestPotential(lmax)
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call genPolyhedron(i,i,NumAtoms,AtomPosition,AtomRadius)
!     ----------------------------------------------------------------
      call printPolyhedron(i)
!     ----------------------------------------------------------------
      call genRadialGrid(i,getInscrSphRadius(i),getOutscrSphRadius(i), &
                         ndivin,ndivout,nmult)
!     ----------------------------------------------------------------
      call printRadialGrid(i)
!     ----------------------------------------------------------------
      call testPotential(i,1,x(i),y(i),z(i))
!     ----------------------------------------------------------------
   enddo
!
!  -------------------------------------------------------------------
   call endTestPotential()
   call endPotentialType()
   call endPolyhedra()
   call endRadialGrid()
!  -------------------------------------------------------------------
!
   deallocate(AtomPosition, AtomRadius)
   deallocate(lmax, v0, u, a, x, y, z)
!
   stop 'Ok'
end program tst_TestPotential
