
      PROGRAM MURN
C
C     modified from murn2.f : R.S 19.04.91
C     modified from murn1.f : gps 6.12.90
C     least-squares fit of E vs V by Murnaghan equation of state
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EDATA(50),VDATA(50),X(40),AUX(40)
      EXTERNAL FUN
C CONVERSION FACTOR FROM RYDBERG TO EV:
      data CONVYY /13.605826d0/
      COMMON /C/ EDATA,VDATA, VOLMIN,VOLMAX,B0MIN,B0MAX,B0PMIN,B0PMAX,U
     &LA,NNN
C
      IN=5
      IOUT=6
      write(iout,*) 'least-squares fit of etot(v) by Murnaghan eq.'
      write(iout,*)
c      WRITE(IOUT,1002)
c 1002 FORMAT(' LEAST-SQUARES FIT OF ETOT(V) BY MURNAGHAN EQUATION'/1X)
C
C READING THE INPUT DATA:
C
C UNIT OF LENGTH, IN ANGSTROEM, WHICH WILL BE USED THROUGHOUT
C THIS PROGRAM:
      ULA=0.52917715d0
      WRITE(IOUT,1010)ULA
 1010 FORMAT(12X,'UNIT OF LENGTH:   ULA =',F15.6,2X,'ANGSTROEM')
C
C UNITS OF ENERGY USED TO INPUT THE DATA: 1 MEANS RYDBERG, 2 IS EV., 3 
C  = HARTREE
C UNITS OF
      READ(IN,*)IUNIT
      WRITE(IOUT,1060)IUNIT
 1060 FORMAT(12X,'UNIT OF ENERGY: IUNIT =',I4,':')
      IF(IUNIT .EQ. 1) then 
        CONV_inp = 13.605826
        WRITE(IOUT,1070)
      else IF(IUNIT .EQ. 2) then 
        CONV_inp = 1.0
        WRITE(IOUT,1080)
      else IF(IUNIT .EQ. 3) then
        CONV_inp = 2 * 13.605826
        WRITE(IOUT,1081)
      else
        WRITE(IOUT,1090)
        STOP
      ENDIF
 1070 FORMAT(40X,' (ENERGIES INPUT IN RYDBERGS)'/1X)
 1080 FORMAT(40X,' (ENERGIES INPUT IN ELECTRONVOLTS)'/1X)
 1081 FORMAT(40X,' (ENERGIES INPUT IN HARTREES)'/1X)
 1090 FORMAT(' ONLY IUNIT = 1 AND 2 ARE RECOGNIZED; STOP.')
C
C CONVERSION FACTOR FROM A**3 TO THE VOLUME OF UNIT CELL
C (E.G. 0.25, IF THE STRUCTURE IS SUCH THAT V = A**3/4 ):
      READ(IN,*)CONVAV
      WRITE(IOUT,1020)CONVAV
 1020 FORMAT(' IN THE STRUCTURE CONSIDERED, THE UNIT CELL VOLUME = A**3
     & *',F15.6/1X)
c read in alat boundaries and number of points for that fitted Etot sho
c  uld be computed
c read in 
      read(in,*) alat_min, alat_max, nr_alat
      write(iout,
     & '("min and max alat(bohr) to compute Etot for at",i4," points ="
     &,f20.7, " and", f20.7)') 
     & nr_alat, alat_min, alat_max
      if (nr_alat .lt. 2 .or. alat_max .le. alat_min) stop 'something w
     &rong with this range'
C NUMBER OF DATA POINTS A, E(A):
      READ(IN,*) NNN
      WRITE(IOUT,1030)NNN
 1030 FORMAT(9X,'NUMBER OF DATA POINTS: N =',I5)
      IF(NNN .GT. 50) THEN
        WRITE(IOUT,*)' N TOO LARGE, NMAX=50'
        STOP
      ENDIF
C
      WRITE(IOUT,1120)
      WRITE(IOUT,1040)
 1040 FORMAT(1X/5X,'I',5X,'ALATT',9X,'VOLUME',8X,
     1  'ENERGY (EV)',12X,'ENERGY (RY)'/1X)
C
      DO 10 I=1,NNN
      READ(IN,*)ALATT,ENINP
      ALATT=ALATT
      VDATA(I)=ALATT**3*CONVAV
      EDATA(I)=ENINP*CONV_inp
      WRITE(IOUT,1110)I,ALATT,VDATA(I),EDATA(I),ENINP
 1110 FORMAT(1X,I5,2F15.6,2F19.10)
   10 CONTINUE
      WRITE(IOUT,1120)
 1120 FORMAT(1X,79(1H-)/1X)
C
C  STARTING VALUES FOR ITERATION:
C
      WRITE(IOUT,1130)
 1130 FORMAT(6X,'ALATT0',3X, 'VOL0 (ULA**3)',2X,'B0 (MBAR)',5X,'B0PRIME
     &',6X,'E0 (EV)',6X,'E0 (RY)'/1X)
C
C     MIMIMUM IN EDATA(I):
      E0MIN=1.D6
      VOLMIN=1.D6
      DO 100 I=1,NNN
        IF(EDATA(I) .GE. E0MIN) GO TO 100
        E0MIN=EDATA(I)
        VOLMIN=VDATA(I)
  100 CONTINUE
C
      VOL0=VOLMIN
      E0EV=E0MIN
      E0RY=E0*CONVYY
      ALATT0=(VOL0/CONVAV)**(1.D0/3.D0)
C
C     FOR OTHER VARIABLES, WE CHOOSE:
c
      B0=1.D0
      B0PRIM=4.D0
C
      WRITE(IOUT,1140)ALATT0,VOL0,B0,B0PRIM,E0EV,E0RY
 1140 FORMAT(5X,F9.4,3F13.5,2F13.5)
      WRITE(IOUT,1150)
 1150 FORMAT(64X,'(=INITIAL GUESS)'/1X)
C
      ALMIN=0.5D0*ALATT0
      ALMAX=1.5D0*ALATT0
      B0MIN=0.01D0
      B0MAX=10.D0
      B0PMIN=0.01D0
      B0PMAX=10.D0
C
      VOLMIN=ALMIN**3.D0*CONVAV
      VOLMAX=ALMAX**3.D0*CONVAV
      WRITE(IOUT,1210)ALMIN,VOLMIN,B0MIN,B0PMIN,ALMAX,VOLMAX,B0MAX,B0PM
     &AX
 1210 FORMAT(' MIN:',F9.4,3F13.5/' MAX:',F9.4,3F13.5)
C
C     A UNIFORM SHIFT OF ALL THE ENERGIES
C     (WILL NOT INFLUENCE THE B0, B0PRIME, VOL0, ONLY SHIFTS E0):
C
      SHIFT=-E0EV-1.D0
      DO 20 I=1,NNN
   20 EDATA(I)=EDATA(I)+SHIFT
C
C     MURNAGHAN LEAST-SQUARES FIT:
      WRITE(IOUT,1120)
      WRITE(IOUT,1280)
 1280 FORMAT(' FIT BY MURNAGHAN EQUATION EQUATION OF STATE')
      WRITE(IOUT,1120)
      WRITE(IOUT,1190)
 1190 FORMAT(1X,'ITER',1X,'ALATT0',3X, 'VOL0 (ULA**3)',2X,
     1   'B0 (MBAR)',5X,'B0PRIME',6X,'E0 (EV)',4X,'SUM OF SQUARES'/
     2   75X,'IERR'/1X)
C
      X(1)=E0EV+SHIFT
      X(2)=B0
      X(3)=B0PRIM
      X(4)=VOL0
      NVAR=4
      LIM=20
C
      DO 70 I = 1, 75
      CALL DFMND(FUN,X,FFF,NVAR,LIM,AUX,IERR)
C
C     THE RESULT OF 20 ITERATION :
      WRITE(IOUT,1250) 20*I,(X(4)/CONVAV)**(1.D0/3.D0),X(4),X(2),X(3),
     1                   X(1)-SHIFT,FFF,IERR
 1250 FORMAT(1X,I4,F9.4,5F13.5/64X,I15)
C
      IF(IERR .EQ. 0) GO TO 80
C
   70 CONTINUE
   80 CONTINUE
C     THE RESULT OF THE LAST ITERATION:
      WRITE(IOUT,1250)20*I,(X(4)/CONVAV)**(1.D0/3.D0),X(4),X(2),X(3),
     1                   X(1)-SHIFT,FFF,IERR
C
      WRITE(IOUT,1120)
      WRITE(IOUT,1130)
      VOL0=X(4)
      ALATT0=(VOL0/CONVAV)**(1.D0/3.D0)
      B0=X(2)
      B0PRIM=X(3)
      E0EV=X(1)-SHIFT
      E0RY=E0EV/CONVYY
      write(iout,'("alat=",f10.5, " b0=", f10.4, " b0p=", f10.3," E0=",
     &f11.6)') 
     &  alatt0, b0, b0prim, E0EV / 27.212
      WRITE(IOUT,1140)ALATT0,VOL0,B0,B0PRIM,E0EV,E0RY
      alatt0=ula*alatt0
      write(iout,1140)alatt0
      WRITE(IOUT,1120)
c
c     print out fitted energy values for nr_alat lattice
c     constants around the equilibrium value given by alat_min and alat
c  _max
c     cons
c
      write(iout,500)
  500 format(8x,'alat(A)',10x,'energy(Ryd)',12x,'alat(bohr)', 10x,'ener
     &gy(H)')
      write(iout,1120)
      do i=1,nr_alat
c	vol=0.6d0*vol0 + (i-1)*0.05*vol0
c	alatt=(vol/convav)**(1./3.)
c	alatt=alatt*ula
	alatt=alat_min + (i-1)*(alat_max-alat_min)/(nr_alat-1)
        vol = alatt**3 * convav
	alatt=alatt*ula
	call murng1(ula,vol,vol0,b0,b0prim,e0ev,etot)
	etotry=etot/convyy
	write(iout,502) alatt, etotry, alatt/ula, etotry/2
      enddo
  502 format(2(f15.5,f20.7))
c
      STOP
      END
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DOUBLE PRECISION FUNCTION FUN(X)
C FUNCTION TO BE MINIMIZED IN THE LEAST-SQUARES FIT (BY SUBROUTINE DFMN
C  D)
C FUNCTION
C CALCULATES THE SUM OF THE  A B S O L U T E  DEVIATIONS
C            (E(THEOR)-E(EXP))**2
C DIVIDED BY THE NUMBER OF EXP. POINTS,
C ASSUMING FOR EQUATION OF STATE THE MURNAGHAN EXPRESSION.
C
C MEANING OF VARIABLES:
C      X(1) .... E0
C      X(2) .... B0
C      X(3) .... B0PRIM
C      X(4) .... VOL0
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EDATA(50),VDATA(50),X(40)
      COMMON /C/ EDATA,VDATA, VOLMIN,VOLMAX,B0MIN,B0MAX,
     1                                       B0PMIN,B0PMAX,ULA,NNN
C
C         *         *         *         *         *        *
      E0=X(1)
      B0=X(2)
      B0PRIM=X(3)
      VOL0=X(4)
C
      SUM=0.D0
C     THE SUM OF SQUARES:
      DO 10 I=1,NNN
      VOLACT=VDATA(I)
      CALL MURNG1(ULA,VOLACT,VOL0,B0,B0PRIM,E0,ETOT)
      SUM=SUM+(ETOT-EDATA(I))**2
   10 CONTINUE
      FUN=SUM/DFLOAT(NNN)
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE MURNG1(ULA,VOL,VOL0,B0,B0PRIM,E0,ETOT)
C EVALUATION OF THE MURNAGHAN EXPRESSION FOR ENERGY AS A FUNCTION
C OF VOLUME.
C
C INPUT DATA:
C
C      ULA ..... UNIT OF LENGTH, IN ANGSTROEMS, USED HERE.
C      VOL ..... VOLUME, IN THE ABOVE UNITS OF LENGTH CUBED.
C      VOL0 .... VOLUME AT THE ENERGY MINIMUM, IN THE SAME UNITS.
C      B0 ...... BULK MODULUS, IN UNITS MEGABAR.
C      B0PRIM .. PRESSURE DERIVATIVE OF THE BULK MODULUS.
C                SHOULD BE CLOSE NEITHER TO 0 NOR TO 1.
C      E0 ...... AN ADDITIVE CONSTANT (IN ELECTRONVOLTS), ADDED
C                TO THE ENERGY EXPRESSION.
C                (SEE,  PR B28, p 5484: Fu and Ho)
C
C OUTPUT DATA:
C
C      ETOT .... THE ENERGY, INCLUDING THE E0 CONSTANT, IN ELECTRONVOLTS
C
C IF B0 DIFFERS FROM 0. OR FROM 1. BY LESS THAN 1.d-6, THEN
C ETOT IS SET AT +111111.111111 ELECTRONVOLTS.
C
C      *      *      *      *      *      *      *      *
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     CONVERSION FACTOR FROM ERG TO ELECTRONVOLT:
      PARAMETER( CONVXX = 1.60209D0 )
C     1 ELECTRONVOLT = 1.60209 D-12 ERG
C
C
      IF(DABS(B0PRIM)-1.D0 .LT. 1.d-6   .OR.
     1                             DABS(B0PRIM) .LT. 1.d-6) THEN
      ETOT=+111111.111111D0
      RETURN
      END IF
C
      IF(VOL .LT. 0.D0   .OR.   VOL0 .LT. 0.D0) THEN
      ETOT=+111111.111111D0
      RETURN
      END IF
C
      ETOT = E0 - B0*VOL0/B0PRIM *
     1 (((VOL/VOL0)**(1.D0-B0PRIM)-1.D0)/(1.D0-B0PRIM)-VOL/VOL0+1.D0)
     2 *ULA**3/CONVXX
C
      RETURN
      END
C ---------------------------------------------------------------------
C  --
C --------
      SUBROUTINE DFMND(F,X,Y,N,LIM,AUX,IER)
C
C     ******************************************************************
C     *   MINIMIZATION OF A FUNCTION OF SEVERAL VARIABLES              *
C     *   USING POWELL'S ALGORITHM WITHOUT DERIVATIVES                 *
C     ******************************************************************
C
      DOUBLE PRECISION F,X,Y,AUX,EPS,ETA,TOL,
     1     DA,DB,DC,DM,DQ,FA,FB,FC,FL,FM,FS,HD,HQ,HX
      DIMENSION X(1),AUX(1)
C
C SUBROUTINES REQUIRED: THE EXTERNAL FUNCTION F.
C
C INPUT DATA:
C     F .... THE FUNCTION OF N VARIABLES X(1)...X(N) TO BE MINIMIZED
C     X .... X(1) ... X(N) = STARTING GUESS FOR THE VARIABLES;
C            BEWARE: X HAS TO BE DIMENSIONED IN THE MAIN PROGRAM
C            TO CONSIDERABLY MORE THAN N.
C     N .... NUMBER OF VARIABLES; THE DIMENSION OF X AND AUX IN THE
C            CALLING PROGRAM MUST BE, HOWEVER, MUCH HIGHER:
C            - PERHAPS 10 TIMES?
C     LIM .. MAXIMUM NUMBER OF ITERATIONS
C     AUX .. AUXILIARY ARRAY OF THE SAME DIMENSION AS X.
C OUTPUT DATA:
C     X .... X(1) ... X(N) THE RESULTING MINIMUM
C     Y .... VALUE OF THE FUNCTION AT THE MINIMUM
C     IER .. SOME ERROR INDICATION 
C            IERR=0 MEANS 'CONVERGENCE ACHIEVED'.
C
C      *      *      *      *      *      *      *      *
C
      ISW  =IER
      IER  =0
      IF (N) 1,1,2
    1 IER  =1000
      GOTO 109
    2 IF (LIM) 3,3,4
    3 IER  =2000
      GOTO 109
C
C     SET INITIAL VALUES AND SAVE INITIAL ARGUMENT
C
    4 N1   =N+1
      N2   =N+N
      N3   =N+N2
      N4   =N+N3
      N5   =N*N+N4
      EPS  =1.D-15
      ETA  =N*EPS
      DO 5 K=1,N
         AUX(K)=X(K)
         J    =N3+K
    5    AUX(J)=1.D0
      FS   =F(AUX)
      FB   =FS
      I    =1
      IC   =1
      IT   =0
      M    =N4
      MF   =0
      IS   =1
C
C     START ITERATION CYCLE
C
    6 IT   =IT+1
      FL   =FS
      DO 7 K=N1,N2
    7    AUX(K)=0.D0
      ID   =I
      I    =IC
      IW   =0
C
C     START MINIMIZATION ALONG NEXT DIRECTION
C
    8 NS   =0
      IP   =0
      DB   =0.D0
      IF (IW) 10,9,10
    9 HX   =AUX(I)
   10 IF (IT-1) 11,11,14
   11 IF (IS-1) 14,12,14
   12 DM   =.1D0
      IF (DABS(HX)-1.D0) 38,38,13
   13 DM   =-DM*HX
      GOTO 38
   14 IF (IS-2) 18,15,18
   15 IF (IT-1) 17,16,17
   16 DM   =HQ
      GOTO 38
   17 DM   =DQ
      GOTO 38
C
C     INTERPOLATION USING ESTIMATE OF SECOND DERIVATIVE
C
   18 IF (IW-1) 20,19,20
   19 J    =N2+I
      GOTO 21
   20 J    =N3+I
   21 HD   =AUX(J)
      DC   =1.D-2
      IF (IT-2) 23,23,22
   22 DC   =HQ
   23 DM   =DC
      MK   =1
      GOTO 51
   24 DM   =DC*HD
      IF (DM) 26,25,26
   25 DM   =1.D0
   26 DM   =.5D0*DC-(FM-FB)/DM
      MK   =2
      IF (FM-FB) 27,29,29
   27 FC   =FB
      FB   =FM
      DB   =DC
      IF (DM-DB) 28,67,28
   28 DC   =0.D0
      GOTO 51
   29 IF (DM-DB) 31,30,31
   30 DA   =DC
      FA   =FM
      GOTO 37
   31 FC   =FM
      GOTO 51
C
C     ANALYSE INTERPOLATED FUNCTION VALUE
C
   32 IF (FM-FB) 34,33,33
   33 DA   =DM
      FA   =FM
      GOTO 35
   34 DA   =DB
      FA   =FB
      DB   =DM
      FB   =FM
   35 IF ((DC-DA)/(DB-DA)) 36,36,50
   36 IF (DB) 67,37,67
   37 NS   =1
      DM   =-DC
C
C     LINEAR SEARCH FOR SMALLER FUNCTION VALUES
C     ALONG CURRENT DIRECTION
C
   38 IF (NS-15) 43,43,39
   39 IF (FS-FM) 41,40,41
   40 MF   =N+2
      DB   =0.D0
      GOTO 67
   41 IF (DABS(DM)-1.D6) 43,43,42
   42 IER  =100
      GOTO 67
   43 NS   =NS+1
      MK   =3
      GOTO 51
   44 IF (FM-FB) 45,46,47
   45 DA   =DB
      FA   =FB
      DB   =DM
      FB   =FM
      DM   =DM+DM
      GOTO 38
   46 IF (FS-FB) 47,45,47
   47 IF (NS-1) 48,48,49
   48 DA   =DM
      FA   =FM
      DM   =-DM
      GOTO 38
   49 DC   =DM
      FC   =FM
C
C     REFINE MINIMUM USING QUADRATIC INTERPOLATION
C
   50 HD   =(FC-FB)/(DC-DB)+(FA-FB)/(DB-DA)
      DM   =.5D0*(DA+DC)+(FA-FC)/(HD+HD)
      IP   =IP+1
      MK   =4
C
C     STEP ARGUMENT VECTOR AND CALCULATE FUNCTION VALUE
C
   51 IF (IW-1) 54,52,54
   52 DO 53 K=1,N
         L    =M+K
   53    AUX(K)=X(K)+DM*AUX(L)
      GOTO 55
   54 AUX(I)=HX+DM
   55 FM   =F(AUX)
      GOTO (24,32,44,56),MK
C
C     ANALYSE INTERPOLATED FUNCTION VALUE
C
   56 IF (FM-FB) 61,61,57
   57 IF (IP-3) 58,62,62
   58 IF ((DC-DB)/(DM-DB)) 60,60,59
   59 DC   =DM
      FC   =FM
      GOTO 50
   60 DA   =DM
      FA   =FM
      GOTO 50
   61 DB   =DM
      FB   =FM
C
C     CALCULATE NEW ESTIMATE OF SECOND DERIVATIVE
C     ALONG THE CURRENT DIRECTION
C
   62 HD=(HD+HD)/(DC-DA)
      IF (IW-1) 64,63,64
   63 J    =N2+I
      GOTO 65
   64 J    =N3+I
   65 AUX(J)=HD
      IF (FB-FS) 67,66,67
   66 DB   =0.D0
C
C     SAVE ARGUMENT VECTOR WITH SMALLEST FUNCTION VALUE FOUND
C
   67 IF (IW-1) 70,68,70
   68 DO 69 K=1,N
         L    =M+K
         J    =N+K
         HD   =DB*AUX(L)
         AUX(J)=AUX(J)+HD
         HD   =X(K)+HD
         AUX(K)=HD
   69    X(K) =HD
      GOTO 71
   70 J    =N+I
      AUX(J)=AUX(J)+DB
      HD   =HX+DB
      AUX(I)=HD
      X(I) =HD
   71 IF (IER-100) 72,108,72
C
C     DETERMINE DIRECTION FOR NEXT LINEAR SEARCH
C
   72 FS   =FB
      IF (I-N) 74,73,73
   73 I    =0
   74 I    =I+1
      IF (IS) 75,75,80
   75 IF (DB) 77,76,77
   76 IF (I-IC) 8,77,8
   77 IC   =I
      IS   =1
      IF (IT-N) 79,79,78
   78 IW   =1
   79 I    =ID
      GOTO 8
   80 M    =M+N
      IF (M-N5) 82,81,81
   81 M    =N4
   82 IF (IS-1) 83,83,94
   83 IF (I-1) 84,84,85
   84 IW   =1
   85 IF (I-ID) 8,86,8
   86 HQ=0.D0
      DO 87 K=N1,N2
   87    HQ=HQ+AUX(K)*AUX(K)
      IF (HQ) 90,88,90
   88 IF (MF-N1) 108,108,89
   89 IER  =200
      GOTO 108
   90 DQ   =DSQRT(HQ)
      HQ   =DQ
      IF (HQ-1.D0) 92,92,91
   91 HQ   =1.D0
   92 DO 93 K=N1,N2
         L    =M+K-N
   93    AUX(L)=AUX(K)/DQ
      IS   =2
      GOTO 8
C
C     END OF ITERATION CYCLE
C     TEST FOR TERMINATION OF MINIMIZATION
C
   94 IS   =0
      TOL  =EPS
      IF (DABS(FS)-1.D0) 96,96,95
   95 TOL  =EPS*DABS(FS)
   96 IF (FL-FS-TOL) 100,100,97
   97 IF (IT-LIM) 99,99,98
   98 IER  =10
      GOTO 108
   99 MF   =0
      GOTO 6
  100 IF (MF-N1) 102,102,101
  101 IER  =200
      GOTO 108
  102 MF   =MF+1
      DQ   =0.D0
      DO 105 K=1,N
         J    =N+K
         IF (DABS(AUX(K))-1.D0) 103,103,104
  103    DQ   =DQ+DABS(AUX(J))
         GOTO 105
  104    DQ   =DQ+DABS(AUX(J)/AUX(K))
  105    CONTINUE
      IF (DQ-ETA) 108,108,106
  106 IF (MF-N) 6,6,107
  107 IER  =1
  108 Y    =FB
      IF (IER) 111,111,109
  109 IF (ISW+12345) 110,111,110
C 110 CALL WIER(IER,20212)
  110 CONTINUE
  111 RETURN
      END
