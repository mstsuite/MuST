C**********************************************************************
C
      SUBROUTINE WASINV(NLEN,KKRSZ,AWAS,ndim,rhs,ldr,VECS,NLIM,
     >  work,r0,rho,pointer,TOL,ALG,IQMR)
C
C     Purpose:
C
C     Parameters:
C     NDIM  = the declared size of A (input).
C     NLEN  = the size of the matrix (input).
C     KKRSZ =
C     A     = the matrix, NLEN by NLEN (input).
C     rhs   = the right hand sides, nlen by kkrsz
C     VECS  =
C     NLIM  = the maximum number of iterations (input).
c     work  = work space, nlen by 10
C     TOL   = the relative convergence tolerance (input).
C     ALG   = the choice of algorithm (input).
C     IQMR  = indicates an error and will switch to an LU based algorithm
C
C     External routines used:
C
C     Noel M. Nachtigal
C     March 4, 1994
C
C**********************************************************************
C

      integer   ALG
      integer   KKRSZ
      integer   NDIM,ldr
      integer   NLEN
      integer   NLIM
      integer   pointer(kkrsz)
c
      complex*16 AWAS(NDIM,nlen)
      complex*16 rhs(ldr,kkrsz)
      complex*16 VECS(nlen,kkrsz,6)
      complex*16 work(nlen,6)
c
      real*8     TOL
C
C     Miscellaneous parameters.
C
      complex*16 ZONE, ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0), ZZERO = (0.0D0,0.0D0))
      real*8     DZERO
      PARAMETER (DZERO = 0.0D0)
C
C     Variables used by reverse communication.
C
      integer  CB
      integer  CX,cx0
      integer  IERR
      integer  REVCOM
      integer  IQMR
      integer  n
      complex*16 rho(kkrsz)
      real*8   r0(kkrsz)
C
C     Local variables.
C
      integer  I
      integer  J
      integer  INFO(4)
      integer  cols, cols_old
C
c     ==================================================================
      if(kkrsz.le.0) then
	write(6,'(''WASINV: kkrsz incorrect: kkrsz='',i5)') kkrsz
	stop'wasinv'
      endif
c     if alg.le.3 then there is no block-code.
      if(alg.le.3) then
	 call WASINV_p(NLEN,KKRSZ,AWAS,ndim,rhs,ldr,VECS,NLIM,
     &                 vecs(1,1,2),TOL,ALG,IQMR)
	 return
      endif
c     ==================================================================
C     Compute the modified right hand side.
C
c     ==================================================================
      do j=1,kkrsz
	call zcopy(nlen,rhs(1,j),1,vecs(1,j,2),1)
      enddo
      CALL ZGEMM('n','N',nlen,kkrsz,NLEN,-ZONE,AWAS,NDIM,rhs,ldr,
     >               zone,VECS(1,1,2),nlen)
c precondition: note that diagonal is 1
      call ztrtrs('l','n','u',nlen,kkrsz,awas,ndim,vecs(1,1,2),nlen,
     &     info)
c     ==================================================================
C
      cols = kkrsz
      cx = 2
      do j=1,kkrsz
	pointer(j)=j
      enddo

      do n=1,nlim
C
      cols_old = cols
      cols = 0
      cx0 = cx

C     Loop over all the columns.
      DO J = 1, cols_old
C
C     Compute Jth column of the inverse.
C
C
         if(n.eq.1) then
C     Set up call to linear systems solver.
C     Compute true residual norms, generate second starting vector.
C
         INFO(1) = 0
         INFO(2) = 0
         call zcopy(nlen,vecs(1,j,2),1,work(1,2),1)
	 else
         INFO(1) = 8
         INFO(2) = 1
	 do i=1,5
         call zcopy(nlen,vecs(1,j,i),1,work(1,i),1)
	 enddo
c no preconditioner so we fake it
         call zcopy(nlen,work(1,cx0),1,work(1,6),1)
	 endif
C
C     Call the solver.
C
        CALL zmar1 (nlen,NLEN,work(1,1),TOL,INFO,
     &              n,rho(j),r0(j))
C
C     Perform matrix-vector multiplications when needed.
C
        IERR   = INFO(1)
        REVCOM = INFO(2)

        IF (REVCOM.EQ.1) THEN
c     ==================================================================
        CX     = INFO(3)
        CB     = INFO(4)
          cols = cols + 1
	  if(j.gt.cols) then
	     i = pointer(cols)
	     pointer(cols) = pointer(j)
	     rho(cols) = rho(j)
	     r0(cols) = r0(j)
	     pointer(j) = i
             call zcopy(nlen,vecs(1,cols,1),1,vecs(1,j,1),1)
	  endif
	  do i=1,5
             call zcopy(nlen,work(1,i),1,vecs(1,cols,i),1)
	  enddo
        ELSE
	   IF (REVCOM.GT.1) THEN
            WRITE (6,'(''REVCOM is greater than 1.'')')
            IQMR=1
            return
           END IF
C
C     Check why the solver stopped (this could be more compact).
C
         IF (IERR.EQ.0) THEN
c           WRITE (6,'(A32)') 'The residual norm has converged.'
	    GO TO 20
         ELSE IF (IERR.EQ.1) THEN
            WRITE (6,'(A35)') 'Invalid reverse communication call.'
         ELSE IF (IERR.EQ.2) THEN
            WRITE (6,'(A27)') 'Invalid inputs encountered.'
         ELSE IF (IERR.EQ.4) THEN
            WRITE (6,'(A31)') 'The algorithm did not converge.'
         ELSE IF (IERR.EQ.8) THEN
            WRITE (6,'(A25)') 'The algorithm broke down.'
         END IF
         WRITE (6,'(A19,I5)') 'Unknown INFO code: ', IERR
            IQMR=IERR
            return
c     save the converged solution
 20      CALL Zcopy (NLEN,work(1,1),1,vecs(1,j,1),1)
        END IF
C
C     Do the next column.
C
      enddo

      if(cols.eq.0) goto 30
C
C     Multiply work(1,CX) with the preconditioned matrix.
C
C     if(cols.gt.1) then
C     CALL ZGEMm('N','n',nlen,cols,NLEN,ZONE,AWAS,NDIM,
C    >           vecs(1,1,cx),nlen,ZZERO,vecs(1,1,CB),nlen)
C     else
C     CALL ZGEMV('N',nlen,NLEN,ZONE,AWAS,NDIM,
C    >           vecs(1,1,cx),1,ZZERO,vecs(1,1,CB),1)
C     endif
      call zcopy(nlen*cols,vecs(1,1,cx),1,vecs(1,1,6),1)
      call ztrtrs('u','n','u',nlen,cols,awas,ndim,vecs(1,1,6),nlen,
     &     info)
      do i=1,cols
	do j=1,nlen
	  vecs(j,i,cb)=vecs(j,i,cx)-vecs(j,i,6)
	enddo
      enddo
      call ztrtrs('l','n','u',nlen,cols,awas,ndim,vecs(1,1,cb),nlen,
     &     info)
      do i=1,cols
	do j=1,nlen
	  vecs(j,i,cb)=vecs(j,i,cb)+vecs(j,i,6)
	enddo
      enddo
c     ==================================================================
      enddo
      iqmr=1
      return
C
 30   continue
C
C     Compute the unpreconditioned solution.
C
      call ztrtrs('u','n','u',nlen,kkrsz,awas,ndim,vecs,nlen,info)
      do j=1,kkrsz
	call zaxpy(nlen,zone,vecs(1,j,1),1,rhs(1,pointer(j)),1)
      enddo
      RETURN
      END
C
C**********************************************************************
