program tst_read
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : PI, FOUR, TEN2m8
!
   implicit none
!
   character (len=80) :: file_name
!
   integer (kind=IntKind) :: NumNodes, MyNode
   integer (kind=IntKind) :: funit, i, j, fp_pos, MaxRSize, offset
   integer (kind=IntKind) :: Repeats
   integer (kind=IntKind) :: IntSize
   integer (kind=IntKind) :: MaxSize
   integer (kind=IntKind) :: rsize, integer4_size, real8_size
   integer (kind=IntKind) :: itemp(4)
   integer (kind=IntKind), allocatable :: iarray(:)
!
   real (kind=RealKind) :: h
   real (kind=RealKind) :: a0
   real (kind=RealKind), allocatable :: rarray(:)
!
   call c_dtsize(integer4_size,real8_size)
!
   file_name = 'data_v'
   call c_gopen(funit,trim(file_name),len_trim(file_name),'old',len('old'))
!
!  ====================================================================
   fp_pos = 1
   call c_fseek(funit,fp_pos,0)
   call c_read_integer(funit,itemp,4*integer4_size)
   NumNodes = itemp(1)
   Repeats = itemp(2)
   MaxSize = itemp(3)
   IntSize = itemp(4)
   write(6,'(a,4i8)')'NumNodes, Repeats, MaxSize, IntSize = ',         &
                      NumNodes, Repeats, MaxSize, IntSize
   allocate( iarray(IntSize), rarray(MaxSize) )
!  ====================================================================
   fp_pos = 4*integer4_size+1
   call c_fseek(funit,fp_pos,0)
   call c_read_double(funit,h,real8_size)
   write(6,'(a,d15.8)')'h = ',h
!  ====================================================================
   offset = 4*integer4_size+real8_size
!
   MaxRSize = 1000+10*Repeats+2
   iarray(:)=-1
   LOOP_k: do MyNode = 0, NumNodes-1
   LOOP_i: do i = 1, Repeats
      rsize = 1000+10*i+2
      if (rsize > MaxSize) then
         write(6,'(a,2i8)')'rsize > MaxSize',rsize,MaxSize
         exit LOOP_i
      endif
!     =================================================================
      fp_pos = offset+(MyNode*Repeats+i-1)*IntSize*integer4_size+1
      call c_fseek(funit,fp_pos,0)
      call c_read_integer(funit,iarray,IntSize*integer4_size)
!
      if (NumNodes /= iarray(1)) then
         write(6,'(a,2i8)')'iarray(1),NumNodes = ',iarray(1),NumNodes
         exit LOOP_i
      else if (MyNode /= iarray(2)) then
         write(6,'(a,2i8)')'iarray(2),MyNode = ',iarray(2),MyNode
         exit LOOP_i
      else if (Repeats /= iarray(3)) then
         write(6,'(a,2i8)')'iarray(3),Repeats = ',iarray(3),Repeats
         exit LOOP_i
      else if (i /= iarray(4)) then
         write(6,'(a,2i8)')'iarray(4),i = ',iarray(4),i
         exit LOOP_i
      else if (rsize /= iarray(5)) then
         write(6,'(a,2i8)')'iarray(5),rsize = ',iarray(5),rsize
         exit LOOP_i
      endif
!     =================================================================
      fp_pos = offset+NumNodes*Repeats*IntSize*integer4_size+          &
               (MyNode*Repeats+i-1)*real8_size*MaxRSize+1
      call c_fseek(funit,fp_pos,0)
      call c_read_double(funit,rarray,rsize*real8_size)
!
      a0=(MyNode*Repeats+i)/real(NumNodes,kind=RealKind)
      if (abs(rarray(1)-h) >= TEN2m8) then
         write(6,'(a,2d15.8)')'rarray(1),h = ',rarray(1),h
         exit LOOP_i
      else if (abs(rarray(2)-a0) >= TEN2m8) then
         write(6,'(a,2d15.8)')'rarray(2),a0 = ',rarray(2),a0
         exit LOOP_i
      else 
         do j=3,rsize
            if (abs(rarray(j)-(a0+h*(j-1))) >= TEN2m8) then
               write(6,'(a,2d15.8)')'rarray(j),v = ',rarray(j),a0+h*(j-1)
               exit LOOP_i
            endif
         enddo
      endif
!     =================================================================
   enddo LOOP_i
   enddo LOOP_k
!
   call c_close(funit)
!
   deallocate( iarray, rarray )
!
   stop 'Ok'
end program tst_read
