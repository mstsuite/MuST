C**********************************************************************
C
      SUBROUTINE ZRANDN (N,ZX,SEED)
C
C     Purpose:
C     Fills the vector ZX with random numbers  between 0 and 1.  If the
C     SEED is given, it should be odd and positive.  The generator is a
C     fairly unsophisticated one, from Pearson's  "Numerical methods in
C     engineering and science" book.
C
C     Parameters:
C     N    = the dimension of the vector (input).
C     ZX   = the vector to fill with random numbers (output).
C     SEED = the seed for the generator (input).
C
C     Noel M. Nachtigal
C     April 23, 1993
C
C**********************************************************************
C
      INTRINSIC DBLE, DCMPLX, IABS, MOD
C
      integer      N
      integer      SEED
      integer      I
      integer      J
C
C     Local variables.
c
      real*8       IMAGX
      real*8       REALX
c
      complex*16   ZX(N)
C
C     Local variables that are saved from one call to the next.
C
      real*8       DMAX
      integer      IM
      integer      IMAX
      integer      IS
      SAVE         DMAX, IM, IMAX, IS
c
      DATA IM/0/
C
C     Initialize the generator data.
C
      IF (IM.EQ.0) THEN
         J  = 0
         IM = 1
         DO 10 I = 1, 31
            J = J + 1
            IF (IM*2.LE.IM) GO TO 20
            IM = IM * 2
 10      CONTINUE
 20      IMAX = (IM-1) * 2 + 1
         DMAX = DBLE(IMAX)
         DO 30 I = 1, MOD(J,3)
            J = J - 1
            IM = IM / 2
 30      CONTINUE
         IM = IM + 5
         IS = IABS(MOD(IM*30107,IMAX))
      END IF
C
C     Check whether there is a new seed.
C
      IF (SEED.GT.0) IS = (SEED / 2) * 2 + 1
C
C     Here goes the rest.
C
      DO 40 I = 1, N
         REALX = DBLE(IS) / DMAX
         IS    = IABS(MOD(IM*IS,IMAX))
         IMAGX = DBLE(IS) / DMAX
         IS    = IABS(MOD(IM*IS,IMAX))
         ZX(I) = DCMPLX(REALX,IMAGX)
 40   CONTINUE
C
      RETURN
      END
C
C**********************************************************************
