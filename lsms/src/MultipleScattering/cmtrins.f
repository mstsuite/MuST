c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmtrins(ma,mb,kr,kc,rst,cst,fac,a,lda,b,ldb)
c     ================================================================
c
      implicit   none
c
      integer    ma
      integer    mb
      integer    kr
      integer    kc
      integer    rst
      integer    cst
      integer    lda
      integer    ldb
      integer    i,j
c
      complex*16 fac
      complex*16 a(lda,ma)
      complex*16 b(ldb,mb)
c
c     ****************************************************************
c     a :     a matrix with dimension (lda x ma)
c     =
c
c     b :     a matrix assigned within dimension (kr x kc)
c     =
c     construct a new matrix b with dimension (ldb x mb):
c                            =
c
c             | b       |
c             | =       |
c     b   =   |         |
c     =       |   fac*a |
c             |       = |     
c     ****************************************************************
c
      do j=1,kc
c        -------------------------------------------------------------
	 call zcopy(kr,a(1,j),1,b(rst+1,cst+j),1)
	 call zscal(kr,fac,b(rst+1,cst+j),1)
c        -------------------------------------------------------------
      enddo
c
      return
      end
