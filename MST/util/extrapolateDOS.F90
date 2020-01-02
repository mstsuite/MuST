program extrapolateDOS
   use KindParamModule, only : IntKind, RealKind
   use InterpolationModule, only : PolyInterp
!
   implicit none
!
   character (len=150) :: title
!
   integer (kind=IntKind) :: nr, ne, ns, is, ie, i, fu
!
   real (kind=RealKind) :: rlarge, fd
   real (kind=RealKind), allocatable :: r(:), e(:)
   real (kind=RealKind), allocatable, target :: dos(:,:,:)
   real (kind=RealKind) :: exdos(2)
   real (kind=RealKind), pointer :: pdos(:)
!
   fu = 10
!
Data: Averaged over      4 atoms
--------------------------------
Fermi energy:  0.68638901
LIZ Radius: :q
Num Energy Points:   500
--------------------------------
        Energy           DOS

   nr = 
   ne =
   ns = 
!
   allocate(r(nr), e(ne), dos(nr,ns,ne))
!
   do i = 1, 6
      read(fu,'(a)')title
      write(6,'(a)')trim(title)
   enddo
   close(fu)
!
   do i = 1, nr
      open(unit=fu,file=,status='old',form='formatted')
      do i = 1, 6
         read(fu,'(a)')title
         write(6,'(a)')trim(title)
      enddo
      r = 
      if (ns == 1) then
         do ie = 1, ne
            read(fu,*) e(ie),dos(i,1,ie)
         enddo
      else
         do ie = 1, ne
            read(fu,*) e(ie),dos(i,1,ie),dos(i,2,ie)
         enddo
      endif
      close(fu)
   enddo
!
!
   rlarge = 1.0d2
!
   do ie = 1, ne
      do is = 1, ns
         pdos => dos(:,is,ie)
         call PolyInterp(r,pdos,nr,rlarge,exdos(is))
      enddo
      if (ns == 1) then
         write(fu, '(f14.6,e18.6)')e(ie),exdos(1)
      else
         write(fu, '(f10.6,2e14.6)')e(ie),exdos(1),exdos(2)
      endif
   enddo
!
   deallocate(r, e, dos)
!
end program extrapolateDOS
