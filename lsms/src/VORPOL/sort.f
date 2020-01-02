c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sort(n,ra)
c     ================================================================
c
      implicit none
c
      integer l
      integer n
      integer ir
      integer i
      integer j
c
      real*8 ra(n)
      real*8 rra
c
c     ****************************************************************
c     Sorting routine.
c     Sorts an array ra of length n into ascending numerical order
c     See: Numerical Recipes. page 231.
c     ****************************************************************
c
c     ================================================================
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
      go to 10
      end
