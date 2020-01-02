c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trgtol(nr,nc,u,ud,mg,ml)
c     ================================================================
c
c     ****************************************************************
c     input::
c          nr:  number of rows of matrix in L-space
c          nc:  number of columns of matrix in L-space
c          u:   transformation matrix for spin space
c          ud:  complex conjugate and transpose of u matrix
c          mg:  matrix in the global spin space
c
c     ml = u * mg * ud
c
c     output::
c          ml1: (1,1) element of matrix in the local spin space
c          ml1: (2,2) element of matrix in the local spin space
c     ****************************************************************
c
      implicit   none
c
      integer    nr
      integer    nc
      integer    i
c
      complex*16 mg(2*nr,2*nc)
      complex*16 ml(nr,nc,4)
      complex*16 u(4)
      complex*16 ud(4)
      complex*16 a
      complex*16 b
c
c     ================================================================
c     zero out ml.....................................................
c     ----------------------------------------------------------------
      call zeroout(ml,8*nc*nr)
c     ----------------------------------------------------------------
c
      a=u(1)*ud(1)
      b=u(3)*ud(1)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,i),   1,ml(1,i,1),1)
	 call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,1),1)
c        -------------------------------------------------------------
      enddo
      a=u(1)*ud(2)
      b=u(3)*ud(2)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,1),1)
	 call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,1),1)
c        -------------------------------------------------------------
      enddo
c
      a=u(2)*ud(1)
      b=u(4)*ud(1)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,i),   1,ml(1,i,2),1)
	 call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,2),1)
c        -------------------------------------------------------------
      enddo
      a=u(2)*ud(2)
      b=u(4)*ud(2)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,2),1)
	 call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,2),1)
c        -------------------------------------------------------------
      enddo
c
      a=u(1)*ud(3)
      b=u(3)*ud(3)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,i),   1,ml(1,i,3),1)
	 call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,3),1)
c        -------------------------------------------------------------
      enddo
      a=u(1)*ud(4)
      b=u(3)*ud(4)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,3),1)
	 call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,3),1)
c        -------------------------------------------------------------
      enddo
c
      a=u(2)*ud(3)
      b=u(4)*ud(3)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,i),   1,ml(1,i,4),1)
	 call zaxpy(nr,b,mg(nr+1,i),1,ml(1,i,4),1)
c        -------------------------------------------------------------
      enddo
      a=u(2)*ud(4)
      b=u(4)*ud(4)
      do i=1,nc
c        -------------------------------------------------------------
	 call zaxpy(nr,a,mg(1,nc+i),   1,ml(1,i,4),1)
	 call zaxpy(nr,b,mg(nr+1,nc+i),1,ml(1,i,4),1)
c        -------------------------------------------------------------
      enddo
c
c     ================================================================
      return
      end
