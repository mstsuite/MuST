program tst_polyhedra
   use KindParamModule, only : IntKind, RealKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : printPolyhedron
!
   implicit none
!
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: i
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: AtomRadius(:)
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
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call genPolyhedron(i,i,NumAtoms,AtomPosition,AtomRadius)
!     ----------------------------------------------------------------
      call printPolyhedron(i)
!     ----------------------------------------------------------------
   enddo
   deallocate(AtomPosition, AtomRadius)
   call endPolyhedra()
   stop 'Ok'
end program tst_polyhedra
