module WriteFunctionModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind   
   use ErrorHandlerModule, only : ErrorHandler
!
public :: writeFunction
!
   interface writeFunction
      module procedure writeFunction_r1, writeFunction_r2,            &
                       writeFunction_c1, writeFunction_c2
      module procedure writeFunction_2r1, writeFunction_2r2,          &
                       writeFunction_2c1, writeFunction_2c2
      module procedure writeFunction_3r1, writeFunction_3r2
   end interface writeFunction
!
private
   integer (kind=IntKind), parameter :: funit = 103
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_r1(file_name,nrs,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f(:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (present(rpow)) then
      if (rpow >= 0 .and. rpow < 10) then
         write(funit,'(a,i1)') ' ir        r              f(r)*r^',rpow
      else
         write(funit,'(a,i2)') ' ir        r             f(r)*r^',rpow
      endif
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir)*r(ir)**rpow
      enddo
   else
      write(funit,'(a)') ' ir        r               f(r)'
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir)
      enddo
   endif
!
   close(funit)
!
   end subroutine writeFunction_r1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_r2(file_name,nrs,ns,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs, ns
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f(:,:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (ns < 1 .or. ns > 2) then
      call ErrorHandler('writeFunction_r2','Not implemented for a real function with ns other than 1 and 2',ns)
   endif
!
   if (present(rpow)) then
      if (ns == 1) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1)')' ir        r              f(r)*r^',rpow
         else
            write(funit,'(a,i2)')' ir        r             f(r)*r^',rpow
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir,1)*r(ir)**rpow
         enddo
      else if (ns == 2) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1,a)')' ir        r            f(r)*r^',rpow,'-up             f(r)*r^',rpow,'-down'
         else
            write(funit,'(a,i2,a,i2,a)')' ir        r           f(r)*r^',rpow,'-up            f(r)*r^',rpow,'-down'
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8,6x,d15.8)')ir,r(ir),f(ir,1)*r(ir)**rpow,f(ir,2)*r(ir)**rpow
         enddo
      endif
   else
      if (ns == 1) then
         write(funit,'(a)')' ir        r               f(r)'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8)')ir,r(ir),f(ir,1)
         enddo
      else if (ns == 2) then
         write(funit,'(a)')' ir        r             f(r)-up              f(r)-down'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,5x,d15.8,6x,d15.8)')ir,r(ir),f(ir,1),f(ir,2)
         enddo
      endif
   endif
!
   close(funit)
!
   end subroutine writeFunction_r2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_c1(file_name,nrs,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   complex (kind=CmplxKind), intent(in) :: f(:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (present(rpow)) then
      if (rpow >= 0 .and. rpow < 10) then
         write(funit,'(a,i1,a,i1)') ' ir        r           Re[f(r)]*r^',rpow, &
                                                '           Im[f(r)]r^]',rpow
      else
         write(funit,'(a,i2,a,i2)') ' ir        r          Re[f(r)]*r^',rpow,  &
                                                '          Im[f(r)]*r^',rpow
      endif
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir)*r(ir)**rpow
      enddo
   else
      write(funit,'(a)') ' ir        r            Re[f(r)]            Im[f(r)]'
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir)
      enddo
   endif
!
   close(funit)
!
   end subroutine writeFunction_c1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_c2(file_name,nrs,ns,r,f,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs, ns
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   complex (kind=CmplxKind), intent(in) :: f(:,:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (ns < 1 .or. ns > 2) then
      call ErrorHandler('writeFunction_c2','Not implemented for a real function with ns other than 1 and 2',ns)
   endif
!
   if (present(rpow)) then
      if (ns == 1) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1)')' ir        r             Re[f(r)]*r^',rpow, &
                                                  '             Im[f(r)]*r^',rpow
         else
            write(funit,'(a,i2,a,i2)')' ir        r             Re[f(r)]*r^',rpow, &
                                                  '             Im[f(r)]*r^',rpow
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir,1)*r(ir)**rpow
         enddo
      else if (ns == 2) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1,a,i1,a,i1,a)')' ir        r            Re[f(r)]*r^',rpow,    &
                                                              '-up            Im[f(r)]*r^',rpow, &
                                                              '-up            Re[f(r)]*r^',rpow,  &
                                                              '-down          Im[f(r)]*r^',rpow,'-down'
         else
            write(funit,'(a,i2,a,i2,a)')' ir        r            Re[f(r)]*r^',rpow,  &
                                                    '-up          Re[f(r)]*r^',rpow, &
                                                    '-up          Im[f(r)]*r^',rpow, &
                                                    '-down        Im[f(r)]*r^',rpow,'-down'
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8,6x,d15.8))')ir,r(ir),f(ir,1)*r(ir)**rpow,f(ir,2)*r(ir)**rpow
         enddo
      endif
   else
      if (ns == 1) then
         write(funit,'(a)')' ir        r            Re[f(r)]            Im[f(r)]'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f(ir,1)
         enddo
      else if (ns == 2) then
         write(funit,'(a)')' ir        r           Re[f(r)]-up         Im[f(r)]-up      Re[f(r)]-down      Im[f(r)]-down'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8,6x,d15.8))')ir,r(ir),f(ir,1),f(ir,2)
         enddo
      endif
   endif
!
   close(funit)
!
   end subroutine writeFunction_c2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_2r1(file_name,nrs,r,f1,f2,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f1(:)
   real (kind=RealKind), intent(in) :: f2(:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (present(rpow)) then
      if (rpow >= 0 .and. rpow < 10) then
         write(funit,'(a,i1,a,i1)') ' ir        r             f1(r)*r^',rpow,'             f2(r)*r^',rpow
      else
         write(funit,'(a,i2,a,i2)') ' ir        r            f1(r)*r^',rpow,'            f2(r)*r^',rpow
      endif
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f1(ir)*r(ir)**rpow,f2(ir)*r(ir)**rpow
      enddo
   else
      write(funit,'(a)') ' ir        r              f1(r)              f2(r)'
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f1(ir),f2(ir)
      enddo
   endif
!
   close(funit)
!
   end subroutine writeFunction_2r1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_2r2(file_name,nrs,ns,r,f1,f2,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs, ns
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f1(:,:)
   real (kind=RealKind), intent(in) :: f2(:,:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (ns < 1 .or. ns > 2) then
      call ErrorHandler('writeFunction_2r2','Not implemented for a real function with ns other than 1 and 2',ns)
   endif
!
   if (present(rpow)) then
      if (ns == 1) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1)')' ir        r              f1(r)*r^',rpow,'              f2(r)*r^',rpow
         else
            write(funit,'(a,i2,a,i2)')' ir        r             f1(r)*r^',rpow,'              f2(r)*r^',rpow
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f1(ir,1)*r(ir)**rpow,f2(ir,1)*r(ir)**rpow
         enddo
      else if (ns == 2) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(2(a,i1,a,i1,a))')' ir        r            f1(r)*r^',rpow,'-up             f1(r)*r^',rpow,'-down', &
                                                       '            f2(r)*r^',rpow,'-up             f2(r)*r^',rpow,'-down'
         else
            write(funit,'(2(a,i2,a,i2,a))')' ir        r           f1(r)*r^',rpow,'-up            f1(r)*r^',rpow,'-down', &
                                                       '           f2(r)*r^',rpow,'-up            f2(r)*r^',rpow,'-down'
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8,6x,d15.8))')ir,r(ir),f1(ir,1)*r(ir)**rpow,f1(ir,2)*r(ir)**rpow, &
                                                                      f2(ir,1)*r(ir)**rpow,f2(ir,2)*r(ir)**rpow
         enddo
      endif
   else
      if (ns == 1) then
         write(funit,'(a)')' ir        r               f1(r)               f2(r)'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f1(ir,1),f2(ir,1)
         enddo
      else if (ns == 2) then
         write(funit,'(a)')' ir        r            f1(r)-up             f1(r)-down            f2(r)-up             f2(r)-down'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8,6x,d15.8))')ir,r(ir),f1(ir,1),f1(ir,2),f2(ir,1),f2(ir,2)
         enddo
      endif
   endif
!
   close(funit)
!
   end subroutine writeFunction_2r2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_2c1(file_name,nrs,r,f1,f2,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   complex (kind=CmplxKind), intent(in) :: f1(:)
   complex (kind=CmplxKind), intent(in) :: f2(:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (present(rpow)) then
      if (rpow >= 0 .and. rpow < 10) then
         write(funit,'(2(a,i1,a,i1))') ' ir        r           Re[f1(r)]*r^',rpow, &
                                                   '           Im[f1(r)]*r^',rpow, &
                                                   '           Re[f2(r)]*r^',rpow, &
                                                   '           Im[f2(r)]*r^',rpow
      else
         write(funit,'(2(a,i2,a,i2))') ' ir        r          Re[f1(r)]*r^',rpow,  &
                                                   '          Im[f1(r)]*r^',rpow,  &
                                                   '          Re[f2(r)]*r^',rpow,  &
                                                   '          Im[f2(r)]*r^',rpow
      endif
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,4(5x,d15.8))')ir,r(ir),f1(ir)*r(ir)**rpow,f2(ir)*r(ir)**rpow
      enddo
   else
      write(funit,'(a,a)') ' ir        r            Re[f1(r)]            Im[f1(r)]', &
                                       '            Re[f2(r)]            Im[f2(r)]'
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,4(5x,d15.8))')ir,r(ir),f1(ir),f2(ir)
      enddo
   endif
!
   close(funit)
!
   end subroutine writeFunction_2c1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_2c2(file_name,nrs,ns,r,f1,f2,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs, ns
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   complex (kind=CmplxKind), intent(in) :: f1(:,:)
   complex (kind=CmplxKind), intent(in) :: f2(:,:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (ns < 1 .or. ns > 2) then
      call ErrorHandler('writeFunction_2c2','Not implemented for a real function with ns other than 1 and 2',ns)
   endif
!
   if (present(rpow)) then
      if (ns == 1) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(2(a,i1,a,i1))')' ir        r             Re[f1(r)]*r^',rpow, &
                                                     '             Im[f1(r)]*r^',rpow, &
                                                     '             Re[f2(r)]*r^',rpow, &
                                                     '             Im[f2(r)]*r^',rpow
         else
            write(funit,'(2(a,i2,a,i2))')' ir        r             Re[f1(r)]*r^',rpow, &
                                                     '             Im[f1(r)]*r^',rpow, &
                                                     '             Re[f2(r)]*r^',rpow, &
                                                     '             Im[f2(r)]*r^',rpow
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,4(5x,d15.8))')ir,r(ir),f1(ir,1)*r(ir)**rpow,f2(ir,1)*r(ir)**rpow
         enddo
      else if (ns == 2) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(2(a,i1,a,i1,a,i1,a,i1,a))')' ir        r            Re[f1(r)]*r^',rpow,    &
                                                                 '-up            Im[f1(r)]*r^',rpow, &
                                                                 '-up            Re[f1(r)]*r^',rpow, &
                                                                 '-down          Im[f1(r)]*r^',rpow,'-down', &
                                                                 '               Re[f2(r)]*r^',rpow, &
                                                                 '-up            Im[f2(r)]*r^',rpow, &
                                                                 '-up            Re[f2(r)]*r^',rpow, &
                                                                 '-down          Im[f2(r)]*r^',rpow,'-down'
         else
            write(funit,'(2(a,i2,a,i2,a))')' ir        r            Re[f(r)]*r^',rpow,  &
                                                      '-up          Re[f(r)]*r^',rpow, &
                                                      '-up          Im[f(r)]*r^',rpow, &
                                                      '-down        Im[f(r)]*r^',rpow,'-down', &
                                                      '             Re[f(r)]*r^',rpow, &
                                                      '-up          Im[f(r)]*r^',rpow, &
                                                      '-up          Im[f(r)]*r^',rpow, &
                                                      '-down        Im[f(r)]*r^',rpow,'-down'
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,4(5x,d15.8,6x,d15.8))')ir,r(ir),f1(ir,1)*r(ir)**rpow,f1(ir,2)*r(ir)**rpow, &
                                                                      f2(ir,1)*r(ir)**rpow,f2(ir,2)*r(ir)**rpow
         enddo
      endif
   else
      if (ns == 1) then
         write(funit,'(a)')' ir        r            Re[f1(r)]            Im[f1(r)]            Re[f2(r)]            Im[f2(r)]'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,2(5x,d15.8))')ir,r(ir),f1(ir,1),f2(ir,1)
         enddo
      else if (ns == 2) then
         write(funit,'(a,a)')' ir        r           Re[f1(r)]-up      Im[f1(r)]-up      Re[f1(r)]-down      Im[f1(r)]-down', &
                                              '      Re[f2(r)]-up      Im[f2(r)]-up      Re[f2(r)]-down      Im[f2(r)]-down'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,4(5x,d15.8,6x,d15.8))')ir,r(ir),f1(ir,1),f1(ir,2),f2(ir,1),f2(ir,2)
         enddo
      endif
   endif
!
   close(funit)
!
   end subroutine writeFunction_2c2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_3r1(file_name,nrs,r,f1,f2,f3,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f1(:),f2(:),f3(:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (present(rpow)) then
      if (rpow >= 0 .and. rpow < 10) then
         write(funit,'(3(a,i1))') ' ir        r             f1(r)*r^',rpow, &
                                              '             f2(r)*r^',rpow, &
                                              '             f3(r)*r^',rpow
      else
         write(funit,'(3(a,i2))') ' ir        r            f1(r)*r^',rpow,  &
                                              '            f2(r)*r^',rpow,  &
                                              '            f3(r)*r^',rpow
      endif
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,3(5x,d15.8))')ir,r(ir),f1(ir)*r(ir)**rpow, &
                                                          f2(ir)*r(ir)**rpow, &
                                                          f3(ir)*r(ir)**rpow
      enddo
   else
      write(funit,'(a)') ' ir        r              f1(r)              f2(r)              f3(r)'
      do ir = 1, nrs
         write(funit,'(i5,2x,f10.8,3(5x,d15.8))')ir,r(ir),f1(ir),f2(ir),f3(ir)
      enddo
   endif
!
   close(funit)
!
   end subroutine writeFunction_3r1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine writeFunction_3r2(file_name,nrs,ns,r,f1,f2,f3,rpow)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: file_name
!
   integer (kind=IntKind), intent(in) :: nrs, ns
   integer (kind=IntKind), intent(in), optional :: rpow
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind), intent(in) :: r(:)
   real (kind=RealKind), intent(in) :: f1(:,:)
   real (kind=RealKind), intent(in) :: f2(:,:)
   real (kind=RealKind), intent(in) :: f3(:,:)
!
   open(unit=funit,file=trim(file_name),status='unknown')
!
   if (ns < 1 .or. ns > 2) then
      call ErrorHandler('writeFunction_3r2','Not implemented for a real function with ns other than 1 and 2',ns)
   endif
!
   if (present(rpow)) then
      if (ns == 1) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(a,i1,a,i1,a,i1)')' ir        r              f1(r)*r^',rpow, &
                                                       '              f2(r)*r^',rpow, &
                                                       '              f3(r)*r^',rpow
         else
            write(funit,'(a,i2,a,i2,a,i2)')' ir        r             f1(r)*r^',rpow,  &
                                                       '             f2(r)*r^',rpow,  &
                                                       '             f3(r)*r^',rpow
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,3(5x,d15.8))')ir,r(ir),f1(ir,1)*r(ir)**rpow,    &
                                                             f2(ir,1)*r(ir)**rpow,    &
                                                             f3(ir,1)*r(ir)**rpow
         enddo
      else if (ns == 2) then
         if (rpow >= 0 .and. rpow < 10) then
            write(funit,'(3(a,i1,a,i1,a))')' ir        r            f1(r)*r^',rpow,'-up             f1(r)*r^',rpow,'-down', &
                                                       '            f2(r)*r^',rpow,'-up             f2(r)*r^',rpow,'-down', &
                                                       '            f3(r)*r^',rpow,'-up             f3(r)*r^',rpow,'-down'
         else
            write(funit,'(3(a,i2,a,i2,a))')' ir        r           f1(r)*r^',rpow,'-up            f1(r)*r^',rpow,'-down', &
                                                       '           f2(r)*r^',rpow,'-up            f2(r)*r^',rpow,'-down', &
                                                       '           f3(r)*r^',rpow,'-up            f3(r)*r^',rpow,'-down'
         endif
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,3(5x,d15.8,6x,d15.8))')ir,r(ir),f1(ir,1)*r(ir)**rpow,f1(ir,2)*r(ir)**rpow, &
                                                                      f2(ir,1)*r(ir)**rpow,f2(ir,2)*r(ir)**rpow, &
                                                                      f3(ir,1)*r(ir)**rpow,f3(ir,2)*r(ir)**rpow
         enddo
      endif
   else
      if (ns == 1) then
         write(funit,'(a)')' ir        r               f1(r)               f2(r)               f3(r)'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,3(5x,d15.8))')ir,r(ir),f1(ir,1),f2(ir,1),f3(ir,1)
         enddo
      else if (ns == 2) then
         write(funit,'(a,a,a)')' ir        r            f1(r)-up             f1(r)-down', &
                                           '            f2(r)-up             f2(r)-down', &
                                           '            f3(r)-up             f3(r)-down'
         do ir = 1, nrs
            write(funit,'(i5,2x,f10.8,3(5x,d15.8,6x,d15.8))')ir,r(ir),f1(ir,1),f1(ir,2),f2(ir,1),f2(ir,2),f3(ir,1),f3(ir,2)
         enddo
      endif
   endif
!
   close(funit)
!
   end subroutine writeFunction_3r2
!  ===================================================================
end module WriteFunctionModule
