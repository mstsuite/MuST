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
C     This  file  contains  the  routine  for  the  QMR  algorithm  for
C     unsymmetric matrices, using the  three-term recurrence variant of
C     the Lanczos algorithm without look-ahead.
C
C**********************************************************************
C
      SUBROUTINE ZUQMX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
C
C     Purpose:
C     This subroutine uses  the QMR algorithm to solve  linear systems.
C     It runs  the algorithm  to convergence or until  a user-specified
C     limit on the number of  iterations is reached.
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
C     The implementation does not have look-ahead, so it is less robust
C     than the full version.
C
C     Parameters:
C     For a description of  the parameters, see the file `zuqmx.doc' in
C     the current directory.
C
C     External routines used:
C     double precision dlamch(ch)
C        LAPACK routine, computes machine-related constants.
C     double precision dznrm2(n,x,incx)
C        BLAS-1 routine, computes the 2-norm of x.
C     subroutine zaxpy(n,a,x,incx,y,incy)
C        BLAS-1 routine, computes y = a * x + y.
C     subroutine zaxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     double precision zdotu(n,x,incx,y,incy)
C        BLAS-1 routine, computes y' * x.
C     double precision zqmxom(n)
C        User-supplied routine, specifies the QMR scaling factors.
C     subroutine zrandn(n,x,seed)
C        Library routine, fills x with random numbers.
C     subroutine zrotg(a,b,cos,sin)
C        BLAS-1 routine, computes the Givens rotation which rotates the
C        vector [a; b] into [ sqrt(a**2 + b**2); 0 ].
C
C     Noel M. Nachtigal
C     May 25, 1993
C
C**********************************************************************
C
      INTRINSIC CDABS, DABS, DBLE, DCMPLX, DCONJG, DMAX1, DSQRT, MAX0
      INTRINSIC MOD
c
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN, ZROTG, ZUQMXO
C
      integer    INFO(4)
      integer    NDIM
      integer    NLEN
      integer    NLIM
c
      real*8     DLAMCH
      real*8     DZNRM2
      real*8     ZUQMXO
      real*8     TOL
c
      complex*16 ZDOTU
      complex*16 VECS(NDIM,9)
C
C     Miscellaneous parameters.
C
      complex*16  ZONE
      complex*16  ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
c
      real*8      DHUN
      real*8      DONE
      real*8      DTEN
      real*8      DZERO
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, ISN, ISNM1, ISNM2, IVN, IVNM1, IVNP1, IWN, IWNM1
      SAVE    IERR, ISN, ISNM1, ISNM2, IVN, IVNM1, IVNP1, IWN, IWNM1
c
      INTEGER IWNP1, N, RETLBL, TF, TRES, VF
      SAVE    IWNP1, N, RETLBL, TF, TRES, VF

      complex*16     DN, DNM1, DNP1, RHSN, SCSN, SCSNM1, SINN, SINNM1
      SAVE           DN, DNM1, DNP1, RHSN, SCSN, SCSNM1, SINN, SINNM1

      real*8           COSN, COSNM1, CSIN, CSINP1, MAXOMG, OMGN
      SAVE             COSN, COSNM1, CSIN, CSINP1, MAXOMG, OMGN

      real*8           OMGNP1, R0, RESN, RHON, RHONP1, SCVN, SCVNP1
      SAVE             OMGNP1, R0, RESN, RHON, RHONP1, SCVN, SCVNP1

      real*8           SCWN, SCWNP1, TMAX, TMIN, TNRM, UCHK, UNRM
      SAVE             SCWN, SCWNP1, TMAX, TMIN, TNRM, UCHK, UNRM
C
C     Local variables, transient.
C
      integer       INIT
      integer       REVCOM
c
      complex*16    RHN
      complex*16    RHNM1
      complex*16    RHNM2
      complex*16    RHNP1
      complex*16    RHSNP1
      complex*16    ZTMP
c
      real*8        DTMP
      real*8        SCVNM1
      real*8        SCWNM1
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C        REVCOM   RETLBL      Comment
C           0        0    first call, go to label 10
C           1       30    returning from AXB, go to label 30
C           1       60    returning from AXB, go to label 60
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
         ELSE IF (RETLBL.EQ.60) THEN
            GO TO 60
         END IF
      ELSE IF (REVCOM.EQ.2) THEN
         IF (RETLBL.EQ.40) GO TO 40
      END IF
      IERR = 1
      GO TO 80
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)        IERR = 2
      IF (NLEN.LT.1)        IERR = 2
      IF (NLIM.LT.1)        IERR = 2
      IF (NLEN.GT.NDIM)     IERR = 2
      IF (IERR.NE.0) GO TO 80
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
C     Initialize the wrapped indices.
C
      ISNM1 = 5
      ISN   = ISNM1
      IVN   = 3
      IVNP1 = IVN
      IWN   = 8
      IWNP1 = IWN
C
C     Set x_0 = 0 and compute the norm of the initial residual.
C
c     CALL ZAXPBY (NLEN,VECS(1,IVN),ZONE,VECS(1,2),ZZERO,VECS(1,IVN))
c     CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1))
      CALL Zcopy (NLEN,VECS(1,2),1,VECS(1,IVN),1)
      CALL Zscal (NLEN,zzero,VECS(1,1),1)
      R0 = DZNRM2(NLEN,VECS(1,IVN),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 80
C
C     Check whether the auxiliary vector must be supplied.
C
      IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,IWN),1)
C
C     Compute scale factors and check for invariant subspaces.
C
      SCVNP1 = R0
      SCWNP1 = DZNRM2(NLEN,VECS(1,IWN),1)
      IF (SCVNP1.LT.TNRM) IERR = IERR + 16
      IF (SCWNP1.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 80
      DNP1 = ZDOTU(NLEN,VECS(1,IWN),1,VECS(1,IVN),1) / ( SCVNP1 * 
     $SCWNP1 )
      IF ((SCVNP1.GE.TMAX).OR.(SCVNP1.LE.TMIN)) THEN
         ZTMP = DCMPLX(DONE / SCVNP1,DZERO)
c        CALL ZAXPBY (NLEN,VECS(1,IVN),ZTMP,VECS(1,IVN),ZZERO,VECS(1,
c    $IVN))
         CALL Zscal (NLEN,ZTMP,VECS(1,IVN),1)
         SCVNP1 = DONE
      END IF
      IF ((SCWNP1.GE.TMAX).OR.(SCWNP1.LE.TMIN)) THEN
         ZTMP = DCMPLX(DONE / SCWNP1,DZERO)
c        CALL ZAXPBY (NLEN,VECS(1,IWN),ZTMP,VECS(1,IWN),ZZERO,VECS(1,
c    $IWN))
         CALL Zscal (NLEN,ZTMP,VECS(1,IWN),1)
         SCWNP1 = DONE
      END IF
      RHONP1 = SCVNP1
      CSINP1 = SCWNP1
      SCVNP1 = DONE / SCVNP1
      SCWNP1 = DONE / SCWNP1
C
C     Initialize the variables.
C
      N      = 1
      DN     = ZONE
      COSN   = DONE
      RESN   = DONE
      COSNM1 = DONE
      OMGN   = DZERO
      SCSN   = DZERO
      SCVN   = DZERO
      SCWN   = DZERO
      SCSNM1 = DZERO
      SINN   = ZZERO
      SINNM1 = ZZERO
      OMGNP1 = ZUQMXO(N)
      RHSN   = OMGNP1 * R0
      MAXOMG = DONE / OMGNP1
C
C     This is one step of the classical Lanczos algorithm.
C
 20   IVNM1 = IVN
      IVN   = IVNP1
      IVNP1 = MOD(N,2) + 3
      IWNM1 = IWN
      IWN   = IWNP1
      IWNP1 = MOD(N,2) + 8
C
C     Check whether D_n is nonsingular.
C
      DNM1 = DN
      DN   = DNP1
      IF (CDABS(DN).EQ.DZERO) THEN
         IERR = 8
         GO TO 80
      END IF
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,IVN),VECS(1,7))
C
      INFO(2) = 1
      INFO(3) = IVN
      INFO(4) = 7
      RETLBL  = 30
      RETURN
 30   RETLBL = 0
C
C     Compute H_{n-1,n} and build part of the vector v_{n+1}.
C
      SCVNM1 = SCVN
      CSIN   = CSINP1
      SCVN   = SCVNP1
      ZTMP   = CSIN * DN / DNM1 * SCVNM1 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZONE,VECS(1,7),-ZTMP,VECS(1,IVNM1)
     $)
C
C     Have the caller carry out ATXB, then return here.
C        CALL ATXB (VECS(1,IWN),VECS(1,7))
C
      INFO(2) = 2
      INFO(3) = IWN
      INFO(4) = 7
      RETLBL  = 40
      RETURN
 40   RETLBL = 0
C
C     Build part of the vector w_{n+1}.
C
      SCWNM1 = SCWN
      RHON   = RHONP1
      SCWN   = SCWNP1
      ZTMP   = RHON * DN / DNM1 * SCWNM1 / SCWN
      CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZONE,VECS(1,7),-ZTMP,VECS(1,IWNM1)
     $)
C
C     Compute H_{nn} and finish the new vectors.
C
      RHN = SCVN * SCWN * ZDOTU(NLEN,VECS(1,IWN),1,VECS(1,IVNP1),1) / DN
c     CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZONE,VECS(1,IVNP1),-RHN,VECS(1,
c    $IVN))
c     CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZONE,VECS(1,IWNP1),-RHN,VECS(1,
c    $IWN))
      CALL ZAXPY (NLEN,-RHN,VECS(1,IVN),1,VECS(1,IVNP1),1)
      CALL ZAXPY (NLEN,-RHN,VECS(1,IWN),1,VECS(1,IWNP1),1)
C
C     Compute scale factors and check for invariant subspaces.
C
      IERR   = 0
      SCVNP1 = DZNRM2(NLEN,VECS(1,IVNP1),1)
      SCWNP1 = DZNRM2(NLEN,VECS(1,IWNP1),1)
      RHONP1 = SCVN * SCVNP1
      CSINP1 = SCWN * SCWNP1
      RHNP1  = RHONP1
      IF (SCVNP1.LT.TNRM) IERR = IERR + 16
      IF (SCWNP1.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 50
      DNP1   = ZDOTU(NLEN,VECS(1,IWNP1),1,VECS(1,IVNP1),1) / ( SCVNP1 * 
     $SCWNP1 )
      IF ((SCVNP1.GE.TMAX).OR.(SCVNP1.LE.TMIN)) THEN
         ZTMP = DCMPLX(DONE / SCVNP1,DZERO)
c        CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZTMP,VECS(1,IVNP1),ZZERO,
c    $VECS(1,IVNP1))
         CALL Zscal (NLEN,ZTMP,VECS(1,IVNP1),1)
         SCVNP1 = DONE
      END IF
      IF ((SCWNP1.GE.TMAX).OR.(SCWNP1.LE.TMIN)) THEN
         ZTMP = DCMPLX(DONE / SCWNP1,DZERO)
c        CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZTMP,VECS(1,IWNP1),ZZERO,
c    $VECS(1,IWNP1))
         CALL Zscal (NLEN,ZTMP,VECS(1,IWNP1),1)
         SCWNP1 = DONE
      END IF
      SCVNP1 = DONE / SCVNP1
      SCWNP1 = DONE / SCWNP1
C
C     The QMR code starts here.
C     Multiply the new column by the previous omega's.
C     Get the next scaling factor omega(i) and update MAXOMG.
C
 50   RHNM1  = CSIN * DN * OMGN / DNM1
      OMGN   = OMGNP1
      RHN    = OMGN * RHN
      OMGNP1 = ZUQMXO(N+1)
      RHNP1  = OMGNP1 * RHNP1
      MAXOMG = DMAX1(MAXOMG,DONE/OMGN)
C
C     Apply the previous rotations.
C
      RHNM2  = SINNM1 * RHNM1
      RHNM1  = COSNM1 * RHNM1
      COSNM1 = COSN
      SINNM1 = SINN
      ZTMP   = RHNM1
      RHNM1  =  COSNM1 * ZTMP + SINNM1 * RHN
      RHN    = -DCONJG(SINNM1) * ZTMP + COSNM1 * RHN
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
C     Compute the next search direction s_n.
C
      ISNM2  = ISNM1
      ISNM1  = ISN
      ISN    = MOD(N-1,2) + 5
      ZTMP   = SCSNM1 * RHNM2 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,ISN),ZONE,VECS(1,IVN),-ZTMP,VECS(1,ISNM2)
     $)
      SCSNM1 = SCSN
      ZTMP   = SCSNM1 * RHNM1 / SCVN
c     CALL ZAXPBY (NLEN,VECS(1,ISN),ZONE,VECS(1,ISN),-ZTMP,VECS(1,ISNM1)
c    $)
      CALL ZAXPY (NLEN,-ZTMP,VECS(1,ISNM1),1,VECS(1,ISN),1)
      SCSN   = SCVN / RHN
C
C     Compute the new QMR iterate, then scale the search direction.
C
      ZTMP = SCSN * RHSN
c     CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ZTMP,VECS(1,ISN))
      CALL ZAXPY (NLEN,ZTMP,VECS(1,ISN),1,VECS(1,1),1)
      DTMP = CDABS(SCSN)
      IF ((DTMP.GE.TMAX).OR.(DTMP.LE.TMIN)) THEN
c        CALL ZAXPBY (NLEN,VECS(1,ISN),SCSN,VECS(1,ISN),ZZERO,VECS(1,
c    $ISN))
         CALL Zscal (NLEN,SCSN,VECS(1,ISN),1)
         SCSN = ZONE
      END IF
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      RHSN = RHSNP1
      UNRM = DSQRT(DBLE(N+1)) * MAXOMG * CDABS(RHSNP1) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 70
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,1),VECS(1,7))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 7
      RETLBL  = 60
      RETURN
 60   RETLBL = 0
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,2),-ZONE,VECS(1,7))
      RESN = DZNRM2(NLEN,VECS(1,7),1) / R0
      UCHK = RESN
C
C     Output the trace messages and convergence history.
C
 70   IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') N, UNRM, RESN
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
         GO TO 80
      ELSE IF (IERR.NE.0) THEN
         GO TO 80
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 80
      ELSE IF (N.GE.NLIM) THEN
         IERR = 4
         GO TO 80
      END IF
C
C     Update the running counter.
C
      N = N + 1
      GO TO 20
C
C     That's all.
C
 80   ZLIM = N
      RETLBL  = 0
      INFO(1) = IERR
C
      RETURN
      END
C
C**********************************************************************
C
      FUNCTION ZUQMXO (I)
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
      integer  I
c
      real*8   ZUQMXO
C
      ZUQMXO = 1.0D0
C
      RETURN
      END
C
C**********************************************************************
