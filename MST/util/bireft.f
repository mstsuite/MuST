      PROGRAM EXBIRFT
C     Fits energy versus volume data to the form
C
C     E(V) = Sum [ a(n) v^(-2n/3) , {n,0,N}]
C
C     and then calculates E(Vo), Vo, Ko, and Ko'
      PARAMETER (MAXN=100,MAXM=10)
      IMPLICIT DOUBLE PRECISION (A-H,K,O-Z)
      LOGICAL CONA,DISB,ENGR,MINFOUND
      DIMENSION V(MAXN),E(MAXN),A(0:MAXM),P(0:MAXM,MAXN),CONV(3),
     1          EARRAY(0:3)
      CHARACTER*50 DATAFL,OUTFL
      CHARACTER*1 BLANK,ANS,YESU,YESL
      DATA CONV/2.5D-1,1D0,0.50d0/
C! Conversion factor:  V=CONV*a*a*a
      DATA BLANK/' '/,YESU/'Y'/,YESL/'y'/
      DATA AUOAN/5.2917706D-1/
C! 1 au in Angstroms
      DATA RYOEV/1.3605804D1/,THIRD/3.33333333333333333D-1/,
     1     TWO3/6.666666666666666667D-1/
C! 1 Ry in eV
      AU3=AUOAN*AUOAN*AUOAN
100   PRINT 105
105   FORMAT(' Enter name of data file:  ')
      READ(5,115,ERR=1000,END=1000) DATAFL
115   FORMAT(A50)
      IF(DATAFL.EQ.BLANK) STOP
      OPEN(UNIT=11,FILE=DATAFL,STATUS='OLD',ERR=100)
      PRINT 117
117   FORMAT(' Enter 1 for volume, 2 for lattice constant:  ')
      READ(5,*) IWHICH
      CONA=IWHICH.EQ.2
      IF(.NOT.CONA) IWHICH=1
C! Default
      PRINT 119
119   FORMAT(' Enter 1 for FCC, 2 for SC, 3 for BCC:  ')
      READ(5,*) LAT
      PRINT 121
121   FORMAT(' Distances in (1) Angstroms, or (2) Bohrs?  ')
      READ(5,*) IAUNIT
      DISB=IAUNIT.EQ.2
      IF(.NOT.DISB) IAUNIT=1
C! Default
      IF(DISB) THEN
         AUNIT=1D0
      ELSE
         AUNIT=AUOAN
      END IF
      VUNIT=AUNIT*AUNIT*AUNIT
      PRINT 123
123   FORMAT(' Energy units in (1) eV, or (2) Rydbergs?  ')
      READ(5,*) IEUNIT
      ENGR=IEUNIT.EQ.2
      IF(.NOT.ENGR) IEUNIT=1
C! Default
      IF(ENGR) THEN
         EUNIT=1D0
      ELSE
         EUNIT=RYOEV
      END IF
C
C     To find the minimum energy, we'll need an estimated starting value.
C       Since we've got the energies on file, we might as well use the
C       volume of the lowest energy as our estimated Vo.
      EMIN=1D10
      VMIN=0D0
C     Read data in the form X,E (x=volume or lattice constant, E=energy)
C     where Angstroms and eV or atomic units (Bohrs and Rydbergs) are
C     used depending upont the setting of IUNIT
      N=0
199   CONTINUE
200   READ(11,*,ERR=199,END=300) X,EE
      IF(CONA) THEN
         VA=CONV(LAT)*X*X*X
      ELSE
         VA=X
      END IF
      N=N+1
      V(N)=VA/VUNIT
C! Volume in au**3
      E(N)=EE/EUNIT
C! Energy in Rydbergs
      WRITE(6,255) N,VA,EE,V(N),E(N)
255   FORMAT(1X,I5,4F15.5)
      IF(E(N).LT.EMIN) THEN
        EMIN=E(N)
        VMIN=V(N)
      END IF
      GO TO 200
300   CLOSE(11)
C     XSTART is the starting value for the minimum search:
      XSTART=VMIN**(-TWO3)
      WRITE(6,305) EMIN,VMIN,XSTART
305   FORMAT(/' Minimum energy in file = ',F15.5,' at V = ',
     1        2F15.7/)
C
C     Set up the fit.  Note that the functions to be fitted are
C       v^(-2n/3),n=0,1,2,...MAXM
310   PRINT 315, MAXM
315   FORMAT(' Input order of v^(-2/3) polynomial (<=',I3,') :  ')
      READ(5,*,ERR=310,END=100) M
      IF(M.LT.0) GO TO 100
      IF(M.GT.MAXM) GO TO 310
      M1=M+1
      DO 400 I=1,N
C       Establish the basis functions:
        P(0,I)=1D0
        X=V(I)**(-TWO3)
        DO J=1,M
          P(J,I)=X**J
        enddo
400   continue
      CALL LSTSQR(N,M1,E,A,P)
      WRITE(6,415) (I,A(I),I=0,M)
415   FORMAT(/' Fitting coefficients:'/(1X,I5,1PE16.8))
C
C     Now for the error checking:
C
      ERMS=0D0
      EMAX=0D0
      DO 600 I=1,N
        XO=V(I)**(-TWO3)
        CALL PLYEVL(M,0D0,A,XO,0,EARRAY)
        ECK=EARRAY(0)
        ERR=ECK-E(I)
        IF(ABS(ERR).GT.ABS(EMAX)) EMAX=ERR
        ERMS=ERMS+ERR*ERR
        WRITE(6,605) I,V(I),E(I),ECK,ERR
600   continue
605     FORMAT(1X,I5,F12.5,3F15.5)
      ERMS=SQRT(ERMS/N)
C     Convert the errors to eV
      ERMSEV=ERMS*RYOEV
      EMAXEV=EMAX*RYOEV
C     Now we must find the equilibrium volume, VO, if we can.  We'll
C       use Newton's method.  Note that we can write the energy
C       as E(v)=f(v^(-2/3)).  Then every extremum of E is also an
C       extremum of f.
C     Pick a trial value for XO.  Eventually, VO=XO**(-3/2).
C       The logical starting value is VMIN**(-2/3), which we
C       calculated above.  After the first minimization, we should
C       have a good estimate of VO, so we'll use that.
      XO=XSTART
C     Since we can calculate exact second derivatives of the
C       fit for M>1, we can use Newton's method to find
C       the root.  Try no more than 100 times
      DO 700 IT=1,100
        CALL PLYEVL(M,0D0,A,XO,2,EARRAY)
C       EARRAY contains f and its derivatives
        DIFF=EARRAY(1)/EARRAY(2)
        WRITE(6,695) XO,EARRAY(1),EARRAY(2),DIFF
695     FORMAT(1X,1P4E16.7)
        XO=XO-DIFF
        IF(ABS(DIFF).LT.1D-5) GO TO 800
700     CONTINUE
C     No roots found.  Ask if fit should be printed out anyway:
      PRINT 705
705   FORMAT(/' No minimum found after 100 iterations.'
     1       /' Print out fitting parameters?  ',$)
      READ(5,875) ANS
      MINFOUND=.FALSE.
      GO TO 900
C
C     Use this new value of XO as the new starting value
800   XSTART=XO
      MINFOUND=.TRUE.
C     XOI is the 2/3rd root of the volume:
      XOI=1D0/XO
      VO=XOI**1.5D0
C     Now for the other equilibrium constants:
      CALL PLYEVL(M,0D0,A,XO,3,EARRAY)
      WRITE(6,711) XO,(EARRAY(I),I=0,3)
711   FORMAT(1X,1P5E15.6)
      EO=EARRAY(0)
C     Check that we've really found an extremum:
      PO=-TWO3*EARRAY(1)*(XO**2.5D0)
      WRITE(6,715) PO
715   FORMAT(' "Equilibrium" pressure = ',1PE15.6)
      KO=(1D1*EARRAY(1)+4D0*XO*EARRAY(2))*XO/(9D0*VO)
      KOP=(2.5D1*EARRAY(1)+4D0*XO*(6D0*EARRAY(2)+XO*EARRAY(3)))/
     1    (3D0*(5D0*EARRAY(1)+2D0*XO*EARRAY(2)))
      ALAT=(VO/CONV(LAT))**3.33333333333333333D-1
      WRITE(6,835) VO,ALAT,EO,KO,KOP,ERMS,EMAX
835   FORMAT(///' Equilibrium parameters for the Birch-Murnaghan',
     1          ' equation:'
     2//'    Vo = ',F15.5,' Bohr**3  a =',F15.5,' Bohrs',
     3//'    Eo = ',F15.5,' Rydbergs'
     4//'    Ko = ',F15.5,' Rydbergs/Bohr**3'
     5//'    Ko''= ',F15.5
     6//'    RMS error in energy fit = ',F15.5,' Rydbergs'
     7//'    Maximum error in energy fit = ',F15.5,' Rydbergs'/)
C     Now convert to standard units:
      EO=EO*RYOEV
C! In eV
      KO=1.4710756D2*KO
C! In MBar
      VO=VO*AU3
C! In Angstroms **3
      ALAT=ALAT*AUOAN
      WRITE(6,855) VO,ALAT,EO,KO,KOP,ERMSEV,EMAXEV
855   FORMAT(///' Equilibrium parameters for the Birch-Murnaghan',
     1          ' equation:'
     2//'    Vo = ',F15.5,' Angstroms**3   a =',F15.5,' Angstroms',
     3//'    Eo = ',F15.5,' eV'
     4//'    Ko = ',F15.5,' Mbar'
     5//'    Ko''= ',F15.5
     6//'    RMS error in energy fit = ',F15.5,' eV'
     7//'    Maximum error in energy fit = ',F15.5,' eV'//)
      PRINT 865
865   FORMAT(' Print results in output file?  ')
      READ(5,875) ANS
875   FORMAT(A1)
900   IF((ANS.EQ.YESL).OR.(ANS.EQ.YESU)) THEN
        PRINT 885
885     FORMAT(' Name of output file?  ')
        READ(5,115) OUTFL
        OPEN(UNIT=21,FILE=OUTFL,STATUS='NEW')
        WRITE(21,895) M,(A(I),I=0,M)
895     FORMAT(1X,I5/(1P4E20.12))
C       Add the equilibrium functions as a reference:
        IF(MINFOUND) WRITE(21,897) EO,VO,ALAT,KO,KOP
897     FORMAT(' Eo = ',F15.5,' eV'/' Vo = ',F15.5,' A^3'/
     1    ' ao = ',F15.5,' Angstroms'/' Ko = ',F15.5,' MBar'/
     2    ' Ko''= ',F15.5)
C       Also print the errors in the calculation
        WRITE(21,925) ERMSEV,EMAXEV
925     FORMAT(' RMS error = ',1PE11.3,' eV'/' Max error = ',
     1    E11.3,' eV')
        CLOSE(21)
      END IF
      GO TO 310
1000  STOP
      END
C_______________________________________________________________
C
      SUBROUTINE LSTSQR(NN,MM,F,A,P)
C     FITS THE FUNCTION F(X), DEFINED AT THE N POINTS X(I)
C     ( F(I)=F(X(I)) ), TO THE M FUNCTIONS P(I,X)
C     ( P(I,J)=P(I,X(J)) ), USING A LINEARIZED LEAST SQUARES
C     FITTING ROUTINE.  THE COEFFICIENT OF P(I,X) WILL BE
C     RETURNED AS A(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=100,MAXM=11)
      DIMENSION X(MAXN),F(MAXN),A(MAXM),P(MAXM,MAXN)
      DIMENSION G(MAXM,MAXN),AVG2(MAXM),B(MAXM),
     1  C(MAXM,MAXM),SGP(MAXM),D(MAXM,MAXM)
      M=MM
      N=NN
      RN=1D0/DBLE(N)
C     THE G'S LINEAR COMBINATIONS OF THE P'S, CHOSEN BY THE
C     GRAM-SCHMIDT ORTHOGONALIZATION PROCEDURE SO THAT
C     <G(M)*G(N)>=0 IF M.NE.N, WHERE
C     <H>=[H(X(1))+H(X(2))+H(X(3))+...+H(X(N))]/N
C     FOR ANY FUNCTION H(X)
C     CALCULATE THE ITH FUNCTION
      DO 40 I=1,M
         IM=I-1
         IP=I+1
         SG2=0D0
         DO 10 K=IP,M
           SGP(K)=0.0D0
10       continue
C        AT THE JTH POINT
         DO 30 J=1,N
           SUM=0D0
           DO 20 K=1,IM
             SUM=SUM+C(I,K)*G(K,J)
20         continue
           G(I,J)=P(I,J)+SUM
           SG2=SG2+G(I,J)*G(I,J)
           DO K=IP,M
             SGP(K)=SGP(K)+P(K,J)*G(I,J)
           enddo
30       continue
C        AVG2(I)=<G(I)*G(I)>
         AVG2(I)=RN*SG2
C        C(K,I)= -<P(K)*G(I)>/<G(I)*G(I)>
         DO K=IP,M
           C(K,I)=-SGP(K)/SG2
         enddo
40    continue
C     SINCE <G(I)*G(J)>=0 FOR I.NE.J, IT'S TRIVIAL TO FIND
C     THE COEFFICIENTS FOR A LEAST SQUARES FIT OF F TO THE G'S
      DO 60 I=1,M
         SUM=0D0
         DO 50 J=1,N
           SUM=SUM+G(I,J)*F(J)
50       continue
         B(I)=RN*SUM/AVG2(I)
60    continue
C     TO CONVERT THE B'S INTO A'S, WE FIRST NEED TO FIND THE
C     COEFFICIENTS FOR EXPANDING THE G'S IN TERMS OF THE P'S
      DO 80 I=1,M
         D(I,I)=1D0
         IM=I-1
         DO K=1,IM
           SUM=0D0
           DO 70 L=K,IM
             SUM=SUM+C(I,L)*D(L,K)
70         continue
           D(I,K)=SUM
         enddo
80    continue
C     FINALLY, WE CAN CHANGE THE B'S INTO A'S
      DO 100 I=1,M
         SUM=0D0
         DO 90 J=I,M
           SUM=SUM+B(J)*D(J,I)
90       continue
         A(I)=SUM
100   continue
      RETURN
      END
C_______________________________________________________________________
C
      SUBROUTINE  PLYEVL(M,X0,A,X,N,P)
C     Evaluates the polynomial given by
C
C       p(x) = Sum [ a(i) (x-x0)^i ,{i,0,m}]
C
C       and its first N derivatives
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Note the dummy indexing on A and P:
      DIMENSION A(0:1),P(0:1)
      Y=X-X0
C     Zeroth order term (in Y)
      IPROD=1
      DO 100 J=0,N
        P(J)=IPROD*A(M)
        IPROD=IPROD*(M-J)
100   continue
      DO 200 I=M-1,0,-1
        IPROD=1
        DO J=0,N
          IF(IPROD.GT.0D0) P(J)=P(J)*Y+IPROD*A(I)
          IPROD=IPROD*(I-J)
        enddo
200   continue
      RETURN
      END
