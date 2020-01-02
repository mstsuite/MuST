c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trltog(nr,nc,u,ud,ml1,ml2,mg)
c     ================================================================
c
c     ****************************************************************
c     input::
c          nr:  number of rows of matrix in L-space
c          nc:  number of columns of matrix in L-space
c          u:   transformation matrix for spin space
c          ud:  complex conjugate and transpose of u matrix
c          ml1: (1,1) element of matrix in the local spin space
c          ml1: (2,2) element of matrix in the local spin space
c
c     mg = ud * ml * u
c
c     output::
c          mg:  matrix in the global spin space
c     ****************************************************************
c
      implicit   none
c
      integer    nr
      integer    nc
      integer    i
c
      complex*16 u(4)
      complex*16 ud(4)
      complex*16 ml1(nr,nc)
      complex*16 ml2(nr,nc)
      complex*16 mg(2*nr,2*nc)
      complex*16 a
      complex*16 b
c
c     ================================================================
c     zero out mg.....................................................
c     ----------------------------------------------------------------
      call zeroout(mg,8*nc*nr)
c     ----------------------------------------------------------------
c
      a=ud(1)*u(1)
      b=ud(3)*u(2)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,ml1(1,i),1,mg(1,i),1)
	 call zaxpy(nr,b,ml2(1,i),1,mg(1,i),1)
c        -------------------------------------------------------------
      enddo
      a=ud(2)*u(1)
      b=ud(4)*u(2)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,ml1(1,i),1,mg(nr+1,i),1)
	 call zaxpy(nr,b,ml2(1,i),1,mg(nr+1,i),1)
c        -------------------------------------------------------------
      enddo
      a=ud(1)*u(3)
      b=ud(3)*u(4)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,ml1(1,i),1,mg(1,nc+i),1)
	 call zaxpy(nr,b,ml2(1,i),1,mg(1,nc+i),1)
c        -------------------------------------------------------------
      enddo
      a=ud(2)*u(3)
      b=ud(4)*u(4)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,ml1(1,i),1,mg(nr+1,nc+i),1)
	 call zaxpy(nr,b,ml2(1,i),1,mg(nr+1,nc+i),1)
c        -------------------------------------------------------------
      enddo
c
c     ================================================================
      return
      end
