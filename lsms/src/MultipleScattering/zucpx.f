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
C     This  file  contains  the  routines  for  the  QMR algorithm  for
C     unsymmetric  matrices,  using  the  coupled  two-term  recurrence
C     variant of the Lanczos algorithm without look-ahead.
C
C**********************************************************************
C
      SUBROUTINE ZUCPX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
C
C     Purpose:
C     This subroutine uses the QMR algorithm based  on the coupled two-
C     term variant of the Lanczos  process without look-ahead  to solve
C     linear systems.  It runs the algorithm to convergence  or until a
C     user-specified limit on the number of  iterations is reached.
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
C     The algorithm  was first described  in the RIACS Technical Report
C     92.15, "An Implementation of the QMR Method Based on Coupled Two-
C     Term Recurrences",  June 1992.  This implementation does not have
C     look-ahead, so it is less robust than the full version.
C
C     Parameters:
C     For a description of  the parameters, see the file `zucpx.doc' in
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
C     subroutine zrotg(a,b,cos,sin)
C        BLAS-1 routine, computes the Givens rotation which rotates the
C        vector [a; b] into [ sqrt(a**2 + b**2); 0 ].
C     double precision zucpxo(n)
C        User-supplied routine, specifies the QMR scaling factors.
C
C     Noel M. Nachtigal
C     March 30, 1993
C
C**********************************************************************
C
      INTRINSIC CDABS, DABS, DBLE, DCMPLX, DCONJG, DMAX1, DSQRT, MAX0
c
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN, ZROTG, ZUCPXO
c
      integer        INFO(4)
      integer        NDIM
      integer        NLEN
      integer        NLIM
c
      real*8         DLAMCH
      real*8         DZNRM2
      real*8         ZUCPXO
      real*8         TOL
c
      complex*16     ZDOTU
      complex*16     VECS(NDIM,8)
C
C     Miscellaneous parameters.
C
      complex*16     ZONE
      complex*16     ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
c
      real*8         DHUN
      real*8         DONE
      real*8         DTEN
      real*8         DZERO
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, N, RETLBL, TF, TRES, VF
      SAVE    IERR, N, RETLBL, TF, TRES, VF
c
      complex*16     DNN, ENN, SCS, SINN, RHSN
      SAVE           DNN, ENN, SCS, SINN, RHSN
c
      real*8           COSN, GAMN, LNP1N, MAXOMG, OMG, R0, SCPN, SCQN
      SAVE             COSN, GAMN, LNP1N, MAXOMG, OMG, R0, SCPN, SCQN
c
      real*8           SCV, RESN, TMAX, TMIN, TNRM, UCHK, UNRM
      SAVE             SCV, RESN, TMAX, TMIN, TNRM, UCHK, UNRM
C
C     Local variables, transient.
C
      integer        INIT
      integer        REVCOM
c
      complex*16     LNN
      complex*16     RHN
      complex*16     RHNM1
      complex*16     RHNP1
      complex*16     RHSNP1
      complex*16     UNM1N
      complex*16     ZTMP
c
      real*8         GAMNM1
      real*8         SCW
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C        REVCOM   RETLBL      Comment
C           0        0    first call, go to label 10
C           1       30    returning from AXB, go to label 30
C           1       50    returning from AXB, go to label 50
C           2       40    returning from ATXB, go to label 40
C
      REVCOM  = INFO(2)
      INFO(2) = 0
      IF (REVCOM.EQ.0) THEN
         N = 0
         IF (RETLBL.EQ.0) GO TO 10
      ELSE IF (REVCOM.EQ.1) THEN
         IF (RETLBL.EQ.30) THEN
            GO TO 30
         ELSE IF (RETLBL.EQ.50) THEN
            GO TO 50
         END IF
      ELSE IF (REVCOM.EQ.2) THEN
         IF (RETLBL.EQ.40) GO TO 40
      END IF
      IERR = 1
      GO TO 70
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)        IERR = 2
      IF (NLEN.LT.1)        IERR = 2
      IF (NLIM.LT.1)        IERR = 2
      IF (NLEN.GT.NDIM)     IERR = 2
      IF (IERR.NE.0) GO TO 70
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
C     Extract and check the various tolerances.
C
      TNRM = DLAMCH('E') * DTEN
      TMIN = DSQRT(DSQRT(DLAMCH('S')))
      TMAX = DONE / TMIN
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E'))
C
C     Start the trace messages and convergence history.
C
      IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') 0, DONE, DONE
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') 0, DONE, DONE
C
C     Set x_0 = 0 and compute the norm of the initial residual.
C
c     CALL ZAXPBY (NLEN,VECS(1,3),ZONE,VECS(1,2),ZZERO,VECS(1,3)) 
c     CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1))
      CALL Zcopy (NLEN,VECS(1,2),1,VECS(1,3),1)
      CALL Zscal (NLEN,ZZERO,VECS(1,1),1)
      R0 = DZNRM2(NLEN,VECS(1,3),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 70
C
C     Check whether the auxiliary vector must be supplied.
C
      IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,7),1)
C
C     Initialize the variables.
C
      N      = 1
      SCV    = R0
      ENN    = ZONE
      COSN   = DONE
      GAMN   = DONE
      RESN   = DONE
      SCPN   = DONE
      SCQN   = DONE
      SCS    = ZZERO
      SINN   = ZZERO
      LNP1N  = DZERO
      OMG    = ZUCPXO(N)
      RHSN   = OMG * R0
      MAXOMG = DONE / OMG
C
C     This is one step of the coupled two-term Lanczos algorithm.
C     Check whether E_n is nonsingular.
C
 20   IF (ENN.EQ.DZERO) THEN
         IERR = 8
         GO TO 70
      END IF
C
C     Compute scale factor for the vector w_{n}.
C     Check for invariant subspaces, and scale the vectors if needed.
C
      IERR = 0
      SCW  = DZNRM2(NLEN,VECS(1,7),1)
      IF (SCPN*SCV.LT.TNRM) IERR = IERR + 16
      IF (SCQN*SCW.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 70
      GAMNM1 = GAMN
      GAMN   = GAMN * SCPN / SCQN * SCV / SCW
      DNN    = ZDOTU(NLEN,VECS(1,3),1,VECS(1,7),1) / ( SCV * SCW )
      IF ((SCV.GE.TMAX).OR.(SCV.LE.TMIN)) THEN
         ZTMP = DCMPLX(DONE / SCV,DZERO)
c        CALL ZAXPBY (NLEN,VECS(1,3),ZTMP,VECS(1,3),ZZERO,VECS(1,3))
         CALL Zscal (NLEN,ZTMP,VECS(1,3),1)
         SCV = DONE
      END IF
      IF ((SCW.GE.TMAX).OR.(SCW.LE.TMIN)) THEN
         ZTMP = DCMPLX(DONE / SCW,DZERO)
c        CALL ZAXPBY (NLEN,VECS(1,7),ZTMP,VECS(1,7),ZZERO,VECS(1,7))
         CALL Zscal (NLEN,ZTMP,VECS(1,7),1)
         SCW = DONE
      END IF
      SCV = DONE / SCV
      SCW = DONE / SCW
C
C     Build the vectors p_n and q_n.
C
      UNM1N = DNN * LNP1N * GAMNM1 / ( GAMN * ENN )
      ZTMP  = UNM1N * SCPN / SCV
      CALL ZAXPBY (NLEN,VECS(1,4),ZONE,VECS(1,3),-ZTMP,VECS(1,4))
      ZTMP  = UNM1N * SCQN / SCW * GAMN / GAMNM1
      CALL ZAXPBY (NLEN,VECS(1,8),ZONE,VECS(1,7),-ZTMP,VECS(1,8))
      SCPN  = SCV
      SCQN  = SCW
C
C     Check whether D_n is nonsingular.
C
      IF (CDABS(DNN).EQ.DZERO) THEN
         IERR = 8
         GO TO 70
      END IF
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,4),VECS(1,6))
C
      INFO(2) = 1
      INFO(3) = 4
      INFO(4) = 6
      RETLBL  = 30
      RETURN
C
C     Compute q_n^T A p_n.
C
 30   ENN = SCPN * SCQN * ZDOTU(NLEN,VECS(1,8),1,VECS(1,6),1)
C
C     Build the vector v_{n+1}.
C
      LNN = ENN / DNN
      CALL ZAXPBY (NLEN,VECS(1,3),ZONE,VECS(1,6),-LNN,VECS(1,3))
C
C     Have the caller carry out ATXB, then return here.
C        CALL ATXB (VECS(1,8),VECS(1,6))
C
      INFO(2) = 2
      INFO(3) = 8
      INFO(4) = 6
      RETLBL  = 40
      RETURN
C
C     Build the vector w_{n+1}.
C
 40   LNN = ENN / DNN
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,6),-LNN,VECS(1,7))
C
C     Compute scale factor for the vector v_{n+1}.
C
      SCV    = DZNRM2(NLEN,VECS(1,3),1)
      LNP1N  = SCPN * SCV
C
C     The QMR code starts here.
C     Multiply the new column by the previous omega's.
C     Get the next scaling factor omega(i) and update MAXOMG.
C
      RHN    = OMG * LNN
      OMG    = ZUCPXO(N+1)
      RHNP1  = OMG * LNP1N
      MAXOMG = DMAX1(MAXOMG,DONE/OMG)
C
C     Apply the previous rotation.
C
      RHNM1 = SINN * RHN
      RHN   = COSN * RHN
C
C     Compute the rotation for the last element (this also applies it).
C
      CALL ZROTG (RHN,RHNP1,COSN,SINN)
C
C     Apply the new rotation to the right-hand side vector.
C
      RHSNP1 = -DCONJG(SINN) * RHSN
      RHSN   =  COSN * RHSN
C
C     Compute the next search direction s_i.
C
      ZTMP = RHNM1 * SCS / SCPN
      CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,4),-ZTMP,VECS(1,5))
C
C     Compute the new QMR iterate, then scale the search direction.
C
      SCS  = SCPN / RHN
      ZTMP = SCS * RHSN
c     CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ZTMP,VECS(1,5))
      CALL ZAXPY (NLEN,ZTMP,VECS(1,5),1,VECS(1,1),1)
      IF ((CDABS(SCS).GE.TMAX).OR.(CDABS(SCS).LE.TMIN)) THEN
c        CALL ZAXPBY (NLEN,VECS(1,5),SCS,VECS(1,5),ZZERO,VECS(1,5))
         CALL Zscal (NLEN,SCS,VECS(1,5),1)
         SCS = ZONE
      END IF
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      RHSN = RHSNP1
      UNRM = DSQRT(DBLE(N+1)) * MAXOMG * CDABS(RHSNP1) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 60
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,1),VECS(1,6))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 6
      RETLBL  = 50
      RETURN
 50   CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,2),-ZONE,VECS(1,6))
      RESN = DZNRM2(NLEN,VECS(1,6),1) / R0
      UCHK = RESN
C
C     Output the convergence history.
C
 60   IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') N, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') N, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C         1. algorithm converged;
C         2. there is an error condition;
C         3. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100;
C         4. algorithm exceeded the iterations limit.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 70
      ELSE IF (IERR.NE.0) THEN
         GO TO 70
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 70
      ELSE IF (N.GE.NLIM) THEN
         IERR = 4
         GO TO 70
      END IF
C
C     Update the running counter.
C
      N = N + 1
      GO TO 20
C
C     That's all.
C
 70   NLIM = N
      RETLBL  = 0
      INFO(1) = IERR
C
      RETURN
      END
C
C**********************************************************************
C
      FUNCTION ZUCPXO (I)
C
C     Purpose:
C     Returns the scaling parameter OMEGA(I).
C
C     Parameters:
C     I = the index of the parameter OMEGA (input).
C
C     Noel M. Nachtigal
C     March 30, 1993
C
C**********************************************************************
C
      integer I
c
      real*8 ZUCPXO
C
      ZUCPXO = 1.0D0
C
      RETURN
      END
C
C**********************************************************************
