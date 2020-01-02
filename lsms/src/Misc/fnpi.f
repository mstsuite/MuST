c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fnpi()
c     ==========================================================
c
      implicit none
c
      real*8   fnpi
      real*8   one
      real*8   two
c
      parameter (one=1.0,two=2.0)
c
      fnpi=two*asin(one)
c
      return
      end
