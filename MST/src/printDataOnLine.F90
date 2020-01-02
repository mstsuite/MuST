!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDataOnLine( data_name, id, getData)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
   use MathParamModule, only : ZERO
!
   use GroupCommModule, only : getGroupID, getMyPEinGroup
!
   use ScfDataModule, only : TableID, n_spin_pola
!
   use InputModule, only : getKeyValue, getNumKeyValues
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SystemModule, only : getAtomName, getNumAtoms
!
   use RadialGridModule, only : getGridRadius
!
   implicit none
!
   character(len=*), intent(in) :: data_name
!
   integer (kind=IntKind), intent(in) :: id
!
   character (len=80) :: fname, svalue
   character (len=10) :: key_name
!
   integer (kind=IntKind) :: ig, i, j, k, n
   integer (kind=IntKind) :: funit, num_points, num_vectors
   integer (kind=IntKind), parameter :: MaxVecs = 50
!
   real (kind=RealKind) :: r, rg, posi(3), val1, val2
   real (kind=RealKind) :: x0, y0, z0, x1, y1, z1, dx, dy, dz
   real (kind=RealKind) :: vec(3,MaxVecs)
!
   logical :: print_pe = .false.
!
   interface
      function getData(title,id,ia,is,posi) result(val)
         use KindParamModule, only : IntKInd, RealKind
!
         character (len=*), intent(in) :: title
         integer (kind=IntKind), intent(in) :: id, ia, is
         real (kind=RealKind), intent(in) :: posi(3)
         real (kind=RealKind) :: val
      end function getData
   end interface
!
   if ( getKeyValue(TableID,'Origin Grid Vector', svalue) == 0 ) then
      read(svalue,*)x0, y0, z0
   else
      call ErrorHandler('printDataOnLine','Origin Grid Vector','Not exist')
   endif
   num_vectors = getNumKeyValues(TableID,'Grid Vector')
   if (num_vectors < 1) then
      call ErrorHandler('printDataOnLine','Num vectors < 1',num_vectors)
   else if (num_vectors > MaxVecs) then
      call ErrorHandler('printDataOnLine','Num vectors > MaxVecs',num_vectors,MaxVecs)
   endif
   if ( getKeyValue(TableID,'Grid Points', svalue) == 0 ) then
      read(svalue,*)num_points
   else
      call ErrorHandler('printDataOnLine','Grid Points','Not exist')
   endif
!
   if ( getMyPEinGroup(getGroupID('K-Mesh')) + getMyPEinGroup(getGroupID('Energy Mesh')) == 0) then
      print_pe = .true.
   endif
!
   if ( print_pe ) then
      ig = getGlobalIndex(id)
      rg = getGridRadius(id)
      if ( getKeyValue(TableID,'Grid Vector', 3, vec, num_vectors) /= 0 ) then
         call ErrorHandler('printDataOnLine','Grid Vector','Invalid data')
      endif
      do k = 1, num_vectors
         r = sqrt(vec(1,k)**2+vec(2,k)**2+vec(3,k)**2)
         x1 = vec(1,k)*rg/r; y1 = vec(2,k)*rg/r; z1 = vec(3,k)*rg/r;
         dx = (x1-x0)/real(num_points)
         dy = (y1-y0)/real(num_points)
         dz = (z1-z0)/real(num_points)
!
         r = getNumAtoms(); i = floor(log10(r))+1
         r = num_vectors; j = floor(log10(r))+1
         n = 10**(i+j+1) + ig*10**(j+1) + 10**j + k
         write(key_name,'(i10)')n; key_name = adjustl(key_name)
         key_name(1:1) = 'a'; key_name(i+2:i+2) = 'v'
         fname = trim(adjustl(data_name))//'_'//key_name; funit = 10+k
         open(unit=funit, file=trim(fname), status='unknown')
         write(funit,'(a,i5,3a,3f12.5)') '# atom index: ',ig,', atom name: ',getAtomName(ig), &
                                         ', ending point: ', vec(1:3,k)
         write(funit,'(a)') '# x, y, z, r,         data value'
!
         do i = 1, num_points
            posi(1) = dx*(i-1) + x0
            posi(2) = dy*(i-1) + y0
            posi(3) = dz*(i-1) + z0
            r = sqrt(posi(1)**2+posi(2)**2+posi(3)**2)
            val1 = getData(data_name,id,1,1,posi)
            if (val1 < -99999.99999999d0) then
               val1 = -99999.99999999d0
            else if (val1 > 99999.99999999d0) then
               val1 = 99999.99999999d0
            endif
            if (n_spin_pola == 1) then
               write(funit,'(4f15.8,2x,f15.8)')posi(1:3), r, val1
            else
               val2 = getData(data_name,id,1,2,posi)
               if (val2 < -99999.99999999d0) then
                  val2 = -99999.99999999d0
               else if (val2 > 99999.99999999d0) then
                  val2 = 99999.99999999d0
               endif
               write(funit,'(4f15.8,2x,f15.8,2x,f15.8)')posi(1:3), r, val1, val2
            endif
         enddo
!
         call FlushFile(funit)
         close(funit)
      enddo
   endif
!
!  ==================================================================
   end subroutine printDataOnLine
