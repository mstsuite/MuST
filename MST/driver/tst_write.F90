program tst_write
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : PI, FOUR
   use MPPModule, only : initMPP, endMPP, getNumPEs, getMyPE, syncAllPEs
!
   implicit none
!
   character (len=80) :: file_name
!
   integer (kind=IntKind) :: NumNodes, MyNode
   integer (kind=IntKind) :: funit, i, j, fp_pos, MaxRSize, offset
   integer (kind=IntKind), parameter :: Repeats = 36
   integer (kind=IntKind), parameter :: IntSize = 10
   integer (kind=IntKind), parameter :: MaxSize = 3000
   integer (kind=IntKind) :: iarray(10)
   integer (kind=IntKind) :: rsize, integer4_size, real8_size
!
   real (kind=RealKind), parameter :: h = PI/FOUR
   real (kind=RealKind) :: a0
   real (kind=RealKind) :: rarray(MaxSize)
!
   call initMPP()
!
   NumNodes = getNumPEs()
   MyNode = getMyPE()
!
   call c_dtsize(integer4_size,real8_size)
!
   file_name = 'data_w'
   call c_gopen(funit,trim(file_name),len_trim(file_name),           &
                'unknown',len('unknown'))
!
!  ====================================================================
   if (MyNode == 0) then
      fp_pos = 1
      iarray(1)=NumNodes
      iarray(2)=Repeats
      iarray(3)=MaxSize
      iarray(4)=IntSize
      call c_fseek(funit,fp_pos,0)
      call c_write_integer(funit,iarray,4*integer4_size)
!     =================================================================
      fp_pos = 4*integer4_size+1
      call c_fseek(funit,fp_pos,0)
      call c_write_double(funit,h,real8_size)
   endif
!  ====================================================================
   offset = 4*integer4_size+real8_size
!
   MaxRSize = 1000+10*Repeats+2
   iarray(:)=-1
   LOOP_i: do i = 1, Repeats
      rsize = 1000+10*i+2
      if (rsize > MaxSize) then
         write(6,'(a,2i8)')'rsize > MaxSize',rsize,MaxSize
         exit LOOP_i
      endif
!     =================================================================
      iarray(1)=NumNodes
      iarray(2)=MyNode
      iarray(3)=Repeats
      iarray(4)=i
      iarray(5)=rsize
!
      fp_pos = offset+(MyNode*Repeats+i-1)*IntSize*integer4_size+1
      call c_fseek(funit,fp_pos,0)
      call c_write_integer(funit,iarray,IntSize*integer4_size)
!     =================================================================
      a0=(MyNode*Repeats+i)/real(NumNodes,kind=RealKind)
      rarray(1)=h
      rarray(2)=a0
      do j=3,rsize
         rarray(j)=a0+h*(j-1)
      enddo
!
      fp_pos = offset+NumNodes*Repeats*IntSize*integer4_size+          &
               (MyNode*Repeats+i-1)*real8_size*MaxRSize+1
      call c_fseek(funit,fp_pos,0)
      call c_write_double(funit,rarray,rsize*real8_size)
!     =================================================================
   enddo LOOP_i
!
   call syncAllPEs()
!
   call c_close(funit)
!
   call endMPP()
!
   if (MyNode == 0) then
      stop 'Ok'
   else
      stop
   endif
end program tst_write
