C**********************************************************************
C
      SUBROUTINE WASINV_p(NLEN,KKRSZ,AWAS,ndim,rhs,ldr,VECS,NLIM,
     >work,TOL,ALG,IQMR)
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
      INTRINSIC CDABS
c
      EXTERNAL ZUCPX, ZUQMX, ZUTFX
C
      integer   ALG
      integer   KKRSZ
      integer   NDIM,ldr
      integer   NLEN
      integer   NLIM
c
      complex*16 AWAS(NDIM,nlen)
      complex*16 rhs(ldr,kkrsz)
      complex*16 VECS(nlen,kkrsz)
      complex*16 work(nlen,9)
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
      integer  CX
      integer  IERR
      integer  REVCOM
      integer  IQMR
C
C     Local variables.
C
      integer  I
      integer  INLIM
      integer  J
      integer  INFO(4)
C
C     Compute the modified right hand side.
C
c     =====================================================================
      do j=1,kkrsz
	call zcopy(nlen,rhs(1,j),1,vecs(1,j),1)
      enddo
      CALL ZGEMM('n','N',nlen,kkrsz,NLEN,-ZONE,AWAS,NDIM,rhs,ldr,
     >               zone,VECS,nlen)
c     =====================================================================
C
C     Loop over all the columns.
C
      DO 30 J = 1, KKRSZ
C
C     Compute Jth column of the inverse.
C
         call zcopy(nlen,vecs(1,j),1,work(1,2),1)
C
C     Set up call to linear systems solver.
C     Compute true residual norms, generate second starting vector.
C
         INFO(2) = 0
c        INFO(1) = 010006
c        INFO(1) = 000006
         INFO(1) = 000000
         INLIM   = NLIM
C
C     Call the solver.
C
 10     IF (ALG.EQ.1) THEN
           CALL ZUCPX (nlen,NLEN,INLIM,work(1,1),TOL,INFO)
        ELSE IF (ALG.EQ.2) THEN
           CALL ZUQMX (nlen,NLEN,INLIM,work(1,1),TOL,INFO)
        ELSE IF (ALG.EQ.3) THEN
           CALL ZUTFX (nlen,NLEN,INLIM,work(1,1),TOL,INFO)
        END IF
C
C     Perform matrix-vector multiplications when needed.
C
        IERR   = INFO(1)
        REVCOM = INFO(2)
        CX     = INFO(3)
        CB     = INFO(4)
C
C     Multiply work(1,CX) with the preconditioned matrix.
C
        IF (REVCOM.EQ.1) THEN
c     ==================================================================
         CALL ZGEMV('N',nlen,NLEN,ZONE,AWAS,NDIM,work(1,CX),1,
     >               ZZERO,work(1,CB),1)
c     ==================================================================
           GO TO 10
C
C     Multiply work(1,CX) with the preconditioned transpose.
C
        ELSE IF (REVCOM.EQ.2) THEN
         CALL ZGEMV('t',nlen,NLEN,ZONE,AWAS,NDIM,work(1,CX),1,
     >               ZZERO,work(1,CB),1)
           GO TO 10
        END IF
C
C     Check why the solver stopped (this could be more compact).
C
         IF (IERR.EQ.0) THEN
c           WRITE (6,'(A32)') 'The residual norm has converged.'
            GO TO 20
         ELSE IF (IERR.EQ.1) THEN
            WRITE (6,'(A35)') 'Invalid reverse communication call.'
            IQMR=IERR
            return
         ELSE IF (IERR.EQ.2) THEN
            WRITE (6,'(A27)') 'Invalid inputs encountered.'
            IQMR=IERR
            return
         ELSE IF (IERR.EQ.4) THEN
            WRITE (6,'(A31)') 'The algorithm did not converge.'
            IQMR=IERR
            return
         ELSE IF (IERR.EQ.8) THEN
            WRITE (6,'(A25)') 'The algorithm broke down.'
            IQMR=IERR
            return
         END IF
         IF ((ALG.EQ.1).OR.(ALG.EQ.2)) THEN
            IF (IERR.EQ.16) THEN
               WRITE (6,'(A39)') 'An A-invariant subspace has been found
     $.'
               GO TO 20 
            ELSE IF (IERR.EQ.32) THEN
               WRITE (6,'(A41)') 'An A^T-invariant subspace has been fou
     $nd.'
               GO TO 20
            ELSE IF (IERR.EQ.48) THEN
               WRITE (6,'(A41)') 'Both invariant subspaces have been fou
     $nd.'
               GO TO 20
            END IF
         END IF
         WRITE (6,'(A19,I5)') 'Unknown INFO code: ', IERR
C
C     Compute the unpreconditioned solution.
C
c20      CALL ZAXPY (NLEN,zone,work(1,1),1,VECS(1,J),1)
 20      CALL Zcopy (NLEN,work(1,1),1,VECS(1,J),1)
C
C     Do the next column.
C
 30   CONTINUE
C
      do j=1,kkrsz
	call zaxpy(nlen,zone,vecs(1,j),1,rhs(1,j),1)
      enddo
      RETURN
      END
C
C**********************************************************************
