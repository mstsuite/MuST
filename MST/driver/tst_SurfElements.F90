program tst_SurfElements
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : printPolyhedron
!
   use SurfElementsModule, only : initSurfElements
   use SurfElementsModule, only : endSurfElements
   use SurfElementsModule, only : genGaussPoints
   use SurfElementsModule, only : testWronskian
!
   use TimerModule, only : initTimer, getTime
!
   implicit none
!
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: i, lmax
!
   real (kind=RealKind) :: t0
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: AtomRadius(:)
!
   complex (kind=CmplxKind) :: kappa
!
   call initTimer()
   t0 = getTime()
!
   read(5,*)bravais(1:3,1)
   read(5,*)bravais(1:3,2)
   read(5,*)bravais(1:3,3)
   read(5,*)NumAtoms
   if (NumAtoms < 1) then
      call ErrorHandler('main','invalid NumAtoms',NumAtoms)
   endif
!  -------------------------------------------------------------------
   call initPolyhedra(NumAtoms,bravais,'main',0)
!  -------------------------------------------------------------------
!
   allocate(AtomPosition(3,NumAtoms), AtomRadius(NumAtoms))
   do i=1,NumAtoms
      read(5,*)AtomPosition(1:3,i),AtomRadius(i)
   enddo
   read(5,*)lmax
   read(5,*)kappa
!
!  -------------------------------------------------------------------
   call initGauntFactors(lmax,'main',0)
!  -------------------------------------------------------------------
!
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call genPolyhedron(i,i,NumAtoms,AtomPosition,AtomRadius)
!     ----------------------------------------------------------------
   enddo
   deallocate(AtomPosition, AtomRadius)
!
!  -------------------------------------------------------------------
   call initSurfElements('main',0)
   call genGaussPoints()
!  -------------------------------------------------------------------
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call testWronskian(i,lmax,kappa)
!     ----------------------------------------------------------------
   enddo
!
   call endSurfElements()
   call endPolyhedra()
   call endGauntFactors()
!
   write(6,'(a,f10.5)')'Elapsed CPU time = ',getTime()-t0
!
   stop 'Ok'
end program tst_SurfElements
