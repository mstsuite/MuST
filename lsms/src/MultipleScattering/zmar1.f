      SUBROUTINE Zmar1 (NDIM,NLEN,VECS,TOL,INFO,
     &           N,RHO,R0)
C
C     Purpose:
C     This subroutine uses the MAR1 algorithm to solve linear systems.
C     It runs the  algorithm to convergence  or until  a user-specified
C     limit on the number of iterations is reached.
C
C     The  code is  set up  to solve  the system  A x = b  with initial
C     guess x_0 = 0.  Here  A x = b  denotes the preconditioned system,
C     and  it  is  connected with the  original system as follows.  Let
C     B y = c be the original unpreconditioned system to be solved, and
C     let y_0 be an arbitrary initial guess for its solution.  Then:
C          A x = b, where  A = M_1^{-1} B M_2^{-1},
C                          x = M_2 (y - y_0), b = M_1^{-1} (c - B y_0).
C     Here M = M_1 M_2 is the preconditioner.
C
C     To recover the final iterate y_n for  the original system B y = c
C     from the final iterate x_n for the preconditioned system A x = b,
C     set
C               y_n = y_0 + M_2^{-1} x_n.
C
C     Parameters:
C     For a description of  the parameters, see the file `zutfx.doc' in
C     the current directory.
C
C     External routines used:
C     double precision dlamch(ch)
C        LAPACK routine, computes machine-related constants.
C     double precision dznrm2(n,x,incx)
C        BLAS-1 routine, computes the 2-norm of x.
C     subroutine zaxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     double precision zdotu(n,x,incx,y,incy)
C        BLAS-1 routine, computes y' * x.
C     subroutine zrandn(n,x,seed)
C        Library routine, fills x with random numbers.
C
C     Noel M. Nachtigal
C     April 13, 1993
C
C**********************************************************************
C
c
      implicit none
      real*8        DLAMCH
      real*8        DZNRM2
      real*8        TOL
C
      integer       INFO(4)
      integer       NDIM
      integer       NLEN
c
      complex*16    ZDOTc
      complex*16    VECS(NDIM,9)
C
C     Miscellaneous parameters.
C
      complex*16    ZONE
      complex*16    ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
c
      real*8        DONE
      real*8        DZERO
      PARAMETER (DONE = 1.0D0,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER  N
c
      complex*16     ALPHA, BETA, RHO
      real*8         R0, RESN
C
C     Local variables, transient.
C
      complex*16     ZTMP(3)
      real*8         det
C
C
      INFO(2) = 0
         IF (N.EQ.1) GO TO 10
         IF (N.EQ.2) GO TO 20
            GO TO 30
C
C     Check whether the inputs are valid.
C
 10   info(1) = 0
      IF (NDIM.LT.1)    info(1) = 2
      IF (NLEN.LT.1)    info(1) = 2
      IF (NLEN.GT.NDIM) info(1) = 2
      IF (info(1).NE.0) return
C
C     Check the convergence tolerance.
C
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E'))
C
C     Set x_0 and compute the norm of the initial residual.
C
      R0 = DZNRM2(NLEN,VECS(1,2),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) return
C
C     Initialize the variables.
C
      RHO  = ZONE
      info(1) = 8

      INFO(2) = 1
      INFO(3) = 2
      INFO(4) = 3
      RETURN
C
C     This is one step of the MAR algorithm.
C     Compute \rho_{n-1}.
C
 20   continue
      call zgemv('c',nlen,2,zone,vecs(1,2),ndim,vecs(1,3),1,
     &     zzero,ztmp,1)
      rho = ztmp(2)
      BETA = conjg(ZTMP(1)) / RHO
C
      call zcopy(nlen,vecs(1,6),1,vecs(1,1),1)
      CALL Zscal (NLEN,beta,vecs(1,1),1)
      call zcopy(nlen,vecs(1,2),1,vecs(1,4),1)
      CALL ZAXPY (NLEN,-beta,vecs(1,3),1,vecs(1,4),1)
      call zcopy(nlen,vecs(1,6),1,vecs(1,5),1)
C
 25   INFO(2) = 1
      INFO(3) = 4
      INFO(4) = 2
      RETURN
 30   continue
      call zgemv('c',nlen,3,zone,vecs(1,2),ndim,vecs(1,2),1,
     &     zzero,ztmp,1)
      det = ztmp(1)*rho-ztmp(2)*conjg(ztmp(2))
      ztmp(3)=conjg(ztmp(3))
      alpha = -ztmp(2)*ztmp(3)/det
      beta = rho*ztmp(3)/det
C
      call zscal(nlen,alpha,vecs(1,5),1)
      CALL ZAXPY (NLEN,beta,vecs(1,6),1,vecs(1,5),1)
      CALL ZAXPY (NLEN,zone,vecs(1,5),1,vecs(1,1),1)
      call zscal(nlen,alpha,vecs(1,3),1)
      CALL ZAXPY (NLEN,beta,vecs(1,2),1,vecs(1,3),1)
      CALL ZAXPY (NLEN,-zone,vecs(1,3),1,vecs(1,4),1)
      rho=ztmp(1)
      resn = DZNRM2(NLEN,VECS(1,4),1)
C
C     Check for convergence or termination.  Stop if:
C         1. algorithm converged;
C         2. algorithm exceeded the iterations limit.
C
      IF (RESN.LE.TOL*r0) THEN
         info(1) = 0
         return
      END IF
C
      GO TO 25
C
      END
C
C**********************************************************************
