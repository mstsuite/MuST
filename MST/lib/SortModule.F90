module SortModule
   use KindParamModule, only : IntKind, RealKind
!
!  define generic procedure for sort
!
public :: QuickSort, &
          HeapSort
!
   interface QuickSort
      module procedure qsort_int, qsort_real
   end interface
!
   interface HeapSort
      module procedure hpsort_int, hpsort_real
   end interface
!
private
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qsort_int(n,arr)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), parameter :: M=7
   integer (kind=IntKind), parameter :: NSTACK=50
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: ir
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jstack
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: istack(NSTACK)
   integer (kind=IntKind) :: flag
!
   integer (kind=intKind), intent(inout) :: arr(n)
   integer (kind=intKind) :: a
   integer (kind=intKind) :: temp
!
   if (n < 2) then
      return
   endif
!
   jstack=0
   l=1
   ir=n
   do
      if(ir-l.lt.M)then
         do j=l+1,ir
            a=arr(j)
            flag=0
            do i=j-1,1,-1
               if (arr(i).le.a) then
                  flag=1
                  exit
               endif
               arr(i+1)=arr(i)
            enddo
            if (flag.eq.0) then
               i=0
            endif
            arr(i+1)=a
         enddo
         if (jstack.eq.0) return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
         endif
         i=l+1
         j=ir
         a=arr(l)
         do 
            i=i+1
            if (arr(i).lt.a) cycle
            flag=0
            do while(flag.eq.0)
               j=j-1
               if (arr(j).le.a) flag=1
            enddo
            if (j.lt.i) then
               exit
            endif
            temp=arr(i)
            arr(i)=arr(j)
            arr(j)=temp
         enddo
         arr(l)=arr(j)
         arr(j)=a
         jstack=jstack+2
         if(jstack.gt.NSTACK)stop 'NSTACK too small in sort'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
   enddo
   end subroutine qsort_int
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qsort_real(n,arr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), parameter :: M=7
   integer (kind=IntKind), parameter :: NSTACK=50
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: ir
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jstack
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: l
   integer (kind=IntKind) :: istack(NSTACK)
   integer (kind=IntKind) :: flag
!
   real (kind=RealKind), intent(inout) :: arr(n)
   real (kind=RealKind) :: a
   real (kind=RealKind) :: temp
!
   if (n < 2) then
      return
   endif
!
   jstack=0
   l=1
   ir=n
   do
      if(ir-l.lt.M)then
         do j=l+1,ir
            a=arr(j)
            flag=0
            do i=j-1,1,-1
               if (arr(i).le.a) then
                  flag=1
                  exit
               endif
               arr(i+1)=arr(i)
            enddo
            if (flag.eq.0) then
               i=0
            endif
            arr(i+1)=a
         enddo
         if (jstack.eq.0) return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
         endif
         i=l+1
         j=ir
         a=arr(l)
         do 
            i=i+1
            if (arr(i).lt.a) cycle
            flag=0
            do while(flag.eq.0)
               j=j-1
               if (arr(j).le.a) flag=1
            enddo
            if (j.lt.i) then
               exit
            endif
            temp=arr(i)
            arr(i)=arr(j)
            arr(j)=temp
         enddo
         arr(l)=arr(j)
         arr(j)=a
         jstack=jstack+2
         if(jstack.gt.NSTACK)stop 'NSTACK too small in sort'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
   enddo
   end subroutine qsort_real
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hpsort_int(n,ra,idx)
!  ===================================================================
!
!  *******************************************************************
!  Heapsorting routine.
!
!  Sorts an array ra of length n into ascending numerical order using
!  the Heapsort algorithm, while make the corresponding rearrangement
!  of the array idx.
!  *******************************************************************
!     
   implicit   none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: idx(n)
   integer (kind=IntKind) :: i,j,l,ir,irb
!
   integer (kind=IntKind), intent(inout) :: ra(n)
   integer (kind=IntKind) :: rra
!
!  ===================================================================
!
   do i=1,n
      idx(i)=i
   end do
!
   if (n < 2) then
      return
   endif
!
   l=n/2+1
   ir=n
   do
      if (l.gt.1)then
         l=l-1
         rra=ra(l)
         irb=idx(l)
      else
         rra=ra(ir)
         irb=idx(ir)
         ra(ir)=ra(1)
         idx(ir)=idx(1)
         ir=ir-1
         if (ir.eq.1)then
            ra(1)=rra
            idx(1)=irb
            return
         end if
      end if
      i=l
      j=l+l
      do while(j.le.ir)
         if (j.lt.ir)then
            if (ra(j).lt.ra(j+1))j=j+1
         end if
         if (rra.lt.ra(j))then
            ra(i)=ra(j)
            idx(i)=idx(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
      end do
      ra(i)=rra
      idx(i)=irb
   enddo
   end subroutine hpsort_int
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine hpsort_real(n,ra,idx)
!  ===================================================================
!
!  *******************************************************************
!     Heapsorting routine.
!
!     Sorts an array ra of length n into ascending numerical order using
!     the Heapsort algorithm, while make the corresponding rearrangement
!     of the array idx.
!  *******************************************************************
!     
   implicit   none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), intent(out) :: idx(n)
   integer (kind=IntKind) :: i,j,l,ir,irb
!
   real (kind=RealKind), intent(inout) :: ra(n)
   real (kind=RealKind) :: rra
!
!  ===================================================================
!
   do i=1,n
      idx(i)=i
   end do
!
   if (n < 2) then
      return
   endif
!
   l=n/2+1
   ir=n
   do
      if (l.gt.1)then
         l=l-1
         rra=ra(l)
         irb=idx(l)
      else
         rra=ra(ir)
         irb=idx(ir)
         ra(ir)=ra(1)
         idx(ir)=idx(1)
         ir=ir-1
         if (ir.eq.1)then
            ra(1)=rra
            idx(1)=irb
            return
         end if
      end if
      i=l
      j=l+l
      do while(j.le.ir)
         if (j.lt.ir)then
            if (ra(j).lt.ra(j+1))j=j+1
         end if
         if (rra.lt.ra(j))then
            ra(i)=ra(j)
            idx(i)=idx(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
      end do
      ra(i)=rra
      idx(i)=irb
   enddo
   end subroutine hpsort_real
!  ===================================================================
end module SortModule
