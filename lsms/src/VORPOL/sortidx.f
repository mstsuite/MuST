      subroutine sortidx(n,ra,idx)
c     ================================================================
c
c     ****************************************************************
c     Heapsorting routine.
c
c     Sorts an array ra of length n into ascending numerical order using
c     the Heapsort algorithm, while make the corresponding rearrangement
c     of the array idx.
c     ****************************************************************
c  input: n  integer scalar, number of elements in ra
c         ra real*8 array of (n), array to be sorted (unchanged on return)
c  returns: idx integer array of (n), ra(idx(i)) is the ith element of
c               the sorted array
c     
      implicit   none
c
      integer    n
      integer    idx(n)
      integer    i,j,l,ir,irb
c
      real*8     ra(n)
      real*8     rra
c
c     ================================================================
c
      do i=1,n
         idx(i)=i
      end do
c
      if(n.le.1) return
      l=n/2+1
      ir=n
10    continue
         if (l.gt.1)then
            l=l-1
            irb=idx(l)
            rra=ra(irb)
         else
            irb=idx(ir)
            rra=ra(irb)
            idx(ir)=idx(1)
            ir=ir-1
            if (ir.eq.1)then
               idx(1)=irb
               return
            end if
         end if
         i=l
         j=l+l
20       if (j.le.ir)then
            if (j.lt.ir)then
               if (ra(idx(j)).lt.ra(idx(j+1)))j=j+1
            end if
            if (rra.lt.ra(idx(j)))then
               idx(i)=idx(j)
               i=j
               j=j+j
            else
               j=ir+1
            end if
            go to 20
         end if
         idx(i)=irb
         go to 10
      end
