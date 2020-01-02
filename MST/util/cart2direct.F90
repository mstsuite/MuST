!  convert the cartition coordinates to direct cooridnates
program cart2direct
   use KindParamModule, only : IntKind, RealKind
   use Matrix3dModule, only : invm3
   implicit none
!
   character (len=80) :: text, tmp
   character (len=2) :: atom_name
!
   real (kind=RealKind) :: b0(3), b1(3), fac
   real (kind=RealKind) :: bra(3,3), brainv(3,3)
!
   integer (kind=IntKind) :: i, j, k, status
!
   interface
      function getTokenPosition(k,s) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer :: p
      end function getTokenPosition
   end interface
!
   j = 0
   do 
      read(5,'(a)',iostat=status) text
      if (status /= 0) then
         exit
      else
         tmp = adjustl(text)
         if (tmp(1:1) == '#' .or. tmp(1:1) == ' ') then
            write(6,'(a)')text
            cycle
         endif
      endif
      if (j <= 3) then
         if (j == 0) then
            read(text,*,iostat=status)fac
            if (status /= 0) then
               write(6,'(a)')text
               stop 'Error'
            endif
!           write(6,'(/,1d12.6)')1.88973d0
            write(6,'(/,a)')'# Insert the scaling factor in the following line:'
            write(6,'(/,a)')'# Bravais Lattice..........'
         else
            read(text,*)bra(1:3,j)
            bra(1:3,j) = bra(1:3,j)*fac
            write(6,'(6x,3d24.16)')bra(1:3,j)
            if (j == 3) then
               write(6,'(/,a)')'DirectCoordinates'
               call invm3(bra,brainv,status)
            endif
         endif
         j = j + 1
      else
         i = getTokenPosition(1,text)
         k = getTokenPosition(2,text)
         atom_name = trim(text(i:k-1))
         read(text(k:),*) b0(1:3)
         b0(1:3) = b0(1:3)*fac
         b1(1:3) = brainv(1:3,1)*b0(1)+brainv(1:3,2)*b0(2)+brainv(1:3,3)*b0(3)
         write(6,'(2x,a,2x,3d24.16)')atom_name,b1(1:3)
         b0(1:3) = bra(1:3,1)*b1(1)+bra(1:3,2)*b1(2)+bra(1:3,3)*b1(3)
!        write(6,'(a,2x,3d24.16)')'check: ',b0(1:3)/fac
      endif
   enddo
!
   stop 'Ok'
end program cart2direct
