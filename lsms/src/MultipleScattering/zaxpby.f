**********************************************************************
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
      SUBROUTINE ZAXPBY (N,ZZ,ZA,ZX,ZB,ZY)
C
C     Purpose:
C     This subroutine computes ZZ = ZA * ZX + ZB * ZY.  Several special
C     cases are handled separately:
C        ZA =  0.0, ZB =  0.0 => ZZ = 0.0
C        ZA =  0.0, ZB =  1.0 => ZZ = ZY  (this is COPY)
C        ZA =  0.0, ZB = -1.0 => ZZ = -ZY
C        ZA =  0.0, ZB =   ZB => ZZ = ZB * ZY  (this is SCAL)
C        ZA =  1.0, ZB =  0.0 => ZZ = ZX  (this is COPY)
C        ZA =  1.0, ZB =  1.0 => ZZ = ZX + ZY
C        ZA =  1.0, ZB = -1.0 => ZZ = ZX - ZY
C        ZA =  1.0, ZB =   ZB => ZZ = ZX + ZB * ZY (this is AXPY)
C        ZA = -1.0, ZB =  0.0 => ZZ = -ZX
C        ZA = -1.0, ZB =  1.0 => ZZ = -ZX + ZY
C        ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY
C        ZA = -1.0, ZB =   ZB => ZZ = -ZX + ZB * ZY
C        ZA =   ZA, ZB =  0.0 => ZZ = ZA * ZX  (this is SCAL)
C        ZA =   ZA, ZB =  1.0 => ZZ = ZA * ZX + ZY  (this is AXPY)
C        ZA =   ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY
C        ZA =   ZA, ZB =   ZB => ZZ = ZA * ZX + ZB * ZY
C     ZZ may be the same as ZX or ZY.
C
C     Parameters:
C     N  = the dimension of the vectors (input).
C     ZZ = the vector result (output).
C     ZA = scalar multiplier for ZX (input).
C     ZX = one of the vectors (input).
C     ZB = scalar multiplier for ZY (input).
C     ZY = the other vector (input).
C
C     Noel M. Nachtigal
C     March 23, 1993
C
C**********************************************************************
C
      INTRINSIC DIMAG, DREAL
C
      integer N
c
      complex*16   ZA
      complex*16   ZB
      complex*16   ZX(N)
      complex*16   ZY(N)
      complex*16   ZZ(N)
C
C     Local variables.
C
      integer      I
c
      real*8       DAI
      real*8       DAR
      real*8       DBI
      real*8       DBR
C
      IF (N.LE.0) RETURN
C
      DAI = DIMAG(ZA)
      DAR = DREAL(ZA)
      DBI = DIMAG(ZB)
      DBR = DREAL(ZB)
      IF ((DAR.EQ.0.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 0.0, ZB = 0.0 => ZZ = 0.0.
            DO 10 I = 1, N
               ZZ(I) = (0.0D0,0.0D0)
 10         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 0.0, ZB = 1.0 => ZZ = ZY (this is COPY).
            DO 20 I = 1, N
               ZZ(I) = ZY(I)
 20         CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 0.0, ZB = -1.0 => ZZ = -ZY.
            DO 30 I = 1, N
               ZZ(I) = -ZY(I)
 30         CONTINUE
         ELSE
C           ZA = 0.0, ZB = ZB => ZZ = ZB * ZY (this is SCAL).
            DO 40 I = 1, N
               ZZ(I) = ZB * ZY(I)
 40         CONTINUE
         END IF
      ELSE IF ((DAR.EQ.1.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 1.0, ZB = 0.0 => ZZ = ZX (this is COPY).
            DO 50 I = 1, N
               ZZ(I) = ZX(I)
 50         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 1.0, ZB = 1.0 => ZZ = ZX + ZY.
            DO 60 I = 1, N
               ZZ(I) = ZX(I) + ZY(I)
 60         CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = 1.0, ZB = -1.0 => ZZ = ZX - ZY.
            DO 70 I = 1, N
               ZZ(I) = ZX(I) - ZY(I)
 70         CONTINUE
         ELSE
C           ZA = 1.0, ZB = ZB => ZZ = ZX + ZB * ZY (this is AXPY).
            DO 80 I = 1, N
               ZZ(I) = ZX(I) + ZB * ZY(I)
 80         CONTINUE
         END IF
      ELSE IF ((DAR.EQ.-1.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = -1.0, ZB = 0.0 => ZZ = -ZX
            DO 90 I = 1, N
               ZZ(I) = -ZX(I)
 90         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = -1.0, ZB = 1.0 => ZZ = -ZX + ZY
            DO 100 I = 1, N
               ZZ(I) = -ZX(I) + ZY(I)
 100        CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY.
            DO 110 I = 1, N
               ZZ(I) = -ZX(I) - ZY(I)
 110        CONTINUE
         ELSE
C           ZA = -1.0, ZB = ZB => ZZ = -ZX + ZB * ZY
            DO 120 I = 1, N
               ZZ(I) = -ZX(I) + ZB * ZY(I)
 120        CONTINUE
         END IF
      ELSE
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = ZA, ZB = 0.0 => ZZ = ZA * ZX (this is SCAL).
            DO 130 I = 1, N
               ZZ(I) = ZA * ZX(I)
 130        CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = ZA, ZB = 1.0 => ZZ = ZA * ZX + ZY (this is AXPY)
            DO 140 I = 1, N
               ZZ(I) = ZA * ZX(I) + ZY(I)
 140        CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C           ZA = ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY.
            DO 150 I = 1, N
               ZZ(I) = ZA * ZX(I) - ZY(I)
 150        CONTINUE
         ELSE
C           ZA = ZA, ZB = ZB => ZZ = ZA * ZX + ZB * ZY.
            DO 160 I = 1, N
               ZZ(I) = ZA * ZX(I) + ZB * ZY(I)
 160        CONTINUE
         END IF
      END IF
C
      RETURN
      END
C
C**********************************************************************
