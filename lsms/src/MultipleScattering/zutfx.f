C**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
C     This file contains the routine for the TFQMR algorithm.
C
C**********************************************************************
C
      SUBROUTINE ZUTFX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
C
C     Purpose:
C     This subroutine uses the TFQMR algorithm to solve linear systems.
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
C     The algorithm was first described  in the RIACS Technical Report
C     91.18, "A Transpose-Free  Quasi-Minimal  Residual Algorithm  for
C     Non-Hermitian Linear Systems", by Roland Freund, September 1991,
C     which subsequently appeared in  SIAM J. Sci. Comput., 14 (1993),
C     pp. 470--482.
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
      INTRINSIC CDABS, DBLE, DSQRT, MAX0
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN
c
      real*8        DLAMCH
      real*8        DZNRM2
      real*8        TOL
C
      integer       INFO(4)
      integer       NDIM
      integer       NLEN
      integer       NLIM
c
      complex*16    ZDOTU
      complex*16    VECS(NDIM,9)
C
C     Miscellaneous parameters.
C
      complex*16    ZONE
      complex*16    ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
c
      real*8        DHUN
      real*8        DONE
      real*8        DTEN
      real*8        DZERO
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, N, RETLBL, TF, TRES, VF
      SAVE    IERR, N, RETLBL, TF, TRES, VF
c
      complex*16     ALPHA, BETA, ETA, RHO
      SAVE           ALPHA, BETA, ETA, RHO
      real*8         COS, VAR, R0, RESN, TAU, UCHK, UNRM
      SAVE           COS, VAR, R0, RESN, TAU, UCHK, UNRM
C
C     Local variables, transient.
C
      integer        INIT
      integer        REVCOM
      complex*16     ZTMP
      real*8         DTMP
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C        REVCOM   RETLBL      Comment
C           0        0    first call, go to label 10
C           1       30    returning from AXB, go to label 30
C           1       40    returning from AXB, go to label 40
C           1       60    returning from AXB, go to label 60
C           1       70    returning from AXB, go to label 70
C
      REVCOM  = INFO(2)
      INFO(2) = 0
      IF (REVCOM.EQ.0) THEN
         N = 0
         IF (RETLBL.EQ.0) GO TO 10
      ELSE IF (REVCOM.EQ.1) THEN
         IF (RETLBL.EQ.30) THEN
            GO TO 30
         ELSE IF (RETLBL.EQ.40) THEN
            GO TO 40
         ELSE IF (RETLBL.EQ.60) THEN
            GO TO 60
         ELSE IF (RETLBL.EQ.70) THEN
            GO TO 70
         END IF
      END IF
      IERR = 1
      GO TO 90
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)    IERR = 2
      IF (NLEN.LT.1)    IERR = 2
      IF (NLIM.LT.1)    IERR = 2
      IF (NLEN.GT.NDIM) IERR = 2
      IF (IERR.NE.0) GO TO 90
C
C     Extract from INFO the output units TF and VF, the true residual
C     flag TRES, and the left starting vector flag INIT.
C
      VF   = MAX0(INFO(1),0)
      INIT = VF / 100000
      VF   = VF - INIT * 100000
      TRES = VF / 10000
      VF   = VF - TRES * 10000
      TF   = VF / 100
      VF   = VF - TF * 100
C
C     Check the convergence tolerance.
C
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E'))
C
C     Start the trace messages and convergence history.
C
      IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') 0, 0, DONE, DONE
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') 0, 0, DONE, DONE
C
C     Set x_0 = 0 and compute the norm of the initial residual.
C
c     CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,2),ZZERO,VECS(1,5))
c     CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1))
      CALL Zcopy (NLEN,VECS(1,2),1,VECS(1,5),1)
      CALL Zscal (NLEN,ZZERO,VECS(1,1),1)
      R0 = DZNRM2(NLEN,VECS(1,5),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 90
C
C     Check whether the auxiliary vector must be supplied.
C
      IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,3),1)
C
C     Initialize the variables.
C
      N    = 1
      RESN = DONE
      RHO  = ZONE
      VAR  = DZERO
      ETA  = ZZERO
      TAU  = R0 * R0
      IERR = 8
c     CALL ZAXPBY (NLEN,VECS(1,8),ZZERO,VECS(1,8),ZZERO,VECS(1,8))
c     CALL ZAXPBY (NLEN,VECS(1,4),ZZERO,VECS(1,4),ZZERO,VECS(1,4))
c     CALL ZAXPBY (NLEN,VECS(1,6),ZZERO,VECS(1,6),ZZERO,VECS(1,6))
      CALL Zscal (NLEN,ZZERO,VECS(1,4),1)
      CALL Zscal (NLEN,ZZERO,VECS(1,6),1)
      CALL Zscal (NLEN,ZZERO,VECS(1,8),1)
C
C     This is one step of the TFQMR algorithm.
C     Compute \beta_{n-1} and \rho_{n-1}.
C
 20   ZTMP = ZDOTU(NLEN,VECS(1,3),1,VECS(1,5),1)
      BETA = ZTMP / RHO
      RHO  = ZTMP
C
C     Compute y_{2n-1}, v_{n-1}, and A y_{2n-1}.
C
      CALL ZAXPBY (NLEN,VECS(1,4),BETA,VECS(1,4),ZONE,VECS(1,8)) 
      CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,5),BETA,VECS(1,6))
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,6),VECS(1,9))
C
      INFO(2) = 1
      INFO(3) = 6
      INFO(4) = 9
      RETLBL  = 30
      RETURN
 30   CALL ZAXPBY (NLEN,VECS(1,4),BETA,VECS(1,4),ZONE,VECS(1,9))
C
C     Compute \sigma{n-1} and check for breakdowns.
C
      ZTMP = ZDOTU(NLEN,VECS(1,3),1,VECS(1,4),1)
      IF ((CDABS(ZTMP).EQ.DZERO).OR.(CDABS(RHO).EQ.DZERO)) THEN
         IERR = 8
         GO TO 90
      END IF
C
C     Compute \alpha_{n-1}, d_{2n-1} and w_{2n}.
C
      ALPHA = RHO / ZTMP
      ZTMP  = VAR * ETA / ALPHA
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,6),ZTMP,VECS(1,7))
c     CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,5),-ALPHA,VECS(1,9))
      CALL ZAXPY (NLEN,-ALPHA,VECS(1,9),1,VECS(1,5),1)
C
C     Compute \varepsilon_{2n-1}^2, \eta_{2n-1}^2, c_{2n-1}^2, and
C     \tau_{2n-1}^2.
C
      DTMP = DZNRM2(NLEN,VECS(1,5),1)
      DTMP = DTMP * DTMP
      VAR  = DTMP / TAU
      COS  = DONE / ( DONE + VAR )
      TAU  = DTMP * COS
      ETA  = ALPHA * COS
C
C     Compute x_{2n-1} and the upper bound for its residual norm.
C
c     CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ETA,VECS(1,7))
      CALL ZAXPY (NLEN,ETA,VECS(1,7),1,VECS(1,1),1)
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      UNRM = DSQRT(DBLE(2*N) * TAU) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN)) GO TO 50
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,1),VECS(1,9))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 9
      RETLBL  = 40
      RETURN
 40   CALL ZAXPBY (NLEN,VECS(1,9),ZONE,VECS(1,2),-ZONE,VECS(1,9))
      RESN = DZNRM2(NLEN,VECS(1,9),1) / R0
      UCHK = RESN
C
C     Output the trace messages and convergence history.
C
 50   IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') N, 2*N-1, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') N, 2*N-1, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C         1. algorithm converged;
C         2. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 90
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 90
      END IF
C
C     Compute y_{2n}, A y_{2n}, d_{2n}, and w_{2n+1}.
C
c     CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,6),-ALPHA,VECS(1,4))
      CALL ZAXPY (NLEN,-ALPHA,VECS(1,4),1,VECS(1,6),1)
      ZTMP = VAR * COS
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,6),ZTMP,VECS(1,7))
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,6),VECS(1,8))
C
      INFO(2) = 1
      INFO(3) = 6
      INFO(4) = 8
      RETLBL  = 60
      RETURN
c60   CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,5),-ALPHA,VECS(1,8))
 60   CALL ZAXPY (NLEN,-ALPHA,VECS(1,8),1,VECS(1,5),1)
C
C     Compute \varepsilon_{2n}^2, \eta_{2n}^2, c_{2n}^2, and
C     \tau_{2n}^2.
C
      DTMP = DZNRM2(NLEN,VECS(1,5),1)
      DTMP = DTMP * DTMP
      VAR  = DTMP / TAU
      COS  = DONE / ( DONE + VAR )
      TAU  = DTMP * COS
      ETA  = ALPHA * COS
C
C     Compute x_{2n}.
C
c     CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ETA,VECS(1,7))
      CALL ZAXPY (NLEN,ETA,VECS(1,7),1,VECS(1,1),1)
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      UNRM = DSQRT(DBLE(2*N+1) * TAU) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 80
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,1),VECS(1,9))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 9
      RETLBL  = 70
      RETURN
 70   CALL ZAXPBY (NLEN,VECS(1,9),ZONE,VECS(1,2),-ZONE,VECS(1,9))
      RESN = DZNRM2(NLEN,VECS(1,9),1) / R0
      UCHK = UNRM
C
C     Output the trace messages and convergence history.
C
 80   IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') N, 2*N, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') N, 2*N, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C         1. algorithm converged;
C         2. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100;
C         3. algorithm exceeded the iterations limit.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 90
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 90
      ELSE IF (N.GE.NLIM) THEN
         IERR = 4
         GO TO 90
      END IF
C
C     Update the running counter.
C
      N = N + 1
      GO TO 20
C
C     That's all.
C
 90   NLIM = N
      RETLBL  = 0
      INFO(1) = IERR
C
      RETURN
      END
C
C**********************************************************************
