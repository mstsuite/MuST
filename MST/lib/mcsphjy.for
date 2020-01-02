        PROGRAM MCSPHJY
C
C       ================================================================
C       Purpose: This program computes the spherical Bessel functions 
C                jn(z), yn(z), and their derivatives for a complex
C                argument using subroutine CSPHJY
C       Input :  z --- Complex argument
C                n --- Order of jn(z) & yn(z) ( 0 ó n ó 250 )
C       Output:  CSJ(n) --- jn(z)
C                CDJ(n) --- jn'(z)
C                CSY(n) --- yn(z)
C                CDY(n) --- yn'(z)
C       Example: z = 4.0+i 2.0
C
C     n     Re[jn(z)]       Im[jn(z)]       Re[jn'(z)]      Im[jn'(z)]
C   --------------------------------------------------------------------
C     0  -.80651523D+00  -.18941093D+00  -.37101203D-01   .75210758D+00
C     1   .37101203D-01  -.75210758D+00  -.67093420D+00   .11885235D+00
C     2   .60314368D+00  -.27298399D+00  -.24288981D+00  -.40737409D+00
C     3   .42955048D+00   .17755176D+00   .18848259D+00  -.24320520D+00
C     4   .12251323D+00   .22087111D+00   .19660170D+00   .17937264D-01
C     5  -.10242676D-01   .10975433D+00   .68951842D-01   .83020305D-01
C
C     n     Re[yn(z)]       Im[yn(z)]       Re[yn'(z)]      Im[yn'(z)]
C   --------------------------------------------------------------------
C     0   .21734534D+00  -.79487692D+00  -.77049661D+00  -.87010064D-02
C     1   .77049661D+00   .87010064D-02  -.92593503D-01  -.64425800D+00
C     2   .24756293D+00   .56894854D+00   .45127429D+00  -.25839924D+00
C     3  -.23845941D+00   .43646607D+00   .26374403D+00   .12439192D+00
C     4  -.27587985D+00   .20902555D+00  -.67092335D-01   .89500599D-01
C     5  -.70001327D-01   .18807178D+00  -.30472133D+00  -.58661384D-01
C       ================================================================
C
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION X,Y
        DIMENSION CSJ(0:250),CDJ(0:250),CSY(0:250),CDY(0:250)
        WRITE(*,*)'Please enter n,x,y (z=x+iy) '
        READ(*,*)N,X,Y
        WRITE(*,30)N,X,Y
30      FORMAT(3X,6HNmaz =,I3,',     ','z = ',F8.1,'+ i',F8.1)
        Z=DCMPLX(X,Y)
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'Please enter order step Ns '
           READ(*,*)NS
        ENDIF
        CALL CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY)
        WRITE(*,*)
        do k=0,n
           write(6,'(1i5,2d20.13)')k,z*z*(csj(k)*cdy(k)-cdj(k)*csy(k))
        enddo
        WRITE(*,*)
        WRITE(*,*)'  n      Re[jn(z)]        Im[jn(z)]',
     &  '        Re[jn''(z)]       Im[jn''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO K=0,NM,NS
           WRITE(*,'(x,i3,4d19.11)')K,CSJ(K),CDJ(K)
        enddo
        WRITE(*,*)
        WRITE(*,*)'  n      Re[yn(z)]        Im[yn(z)]',
     &  '        Re[yn''(z)]       Im[yn''(z)]'
        WRITE(*,*)'--------------------------------------------',
     &  '----------------------------'
        DO K=0,NM,NS
           WRITE(*,'(x,i3,4d19.11)')K,CSY(K),CDY(K)
        enddo
        END


        SUBROUTINE CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY)
C
C       ==========================================================
C       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
C                and their derivatives for a complex argument
C       Input :  z --- Complex argument
C                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
C       Output:  CSJ(n) --- jn(z)
C                CDJ(n) --- jn'(z)
C                CSY(n) --- yn(z)
C                CDY(n) --- yn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ==========================================================
C
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION A0
        DIMENSION CSJ(0:N),CDJ(0:N),CSY(0:N),CDY(0:N)
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-60) THEN
           DO 10 K=0,N
              CSJ(K)=0.0D0
              CDJ(K)=0.0D0
              CSY(K)=-1.0D+300
10            CDY(K)=1.0D+300
           CSJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(.333333333333333D0,0.0D0)
           RETURN
        ENDIF
        CSJ(0)=CDSIN(Z)/Z
        CSJ(1)=(CSJ(0)-CDCOS(Z))/Z
        IF (N.GE.2) THEN
           CSA=CSJ(0)
           CSB=CSJ(1)
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF0=0.0D0
           CF1=1.0D0-100
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF1/Z-CF0
              IF (K.LE.NM) CSJ(K)=CF
              CF0=CF1
15            CF1=CF
           IF (CDABS(CSA).GT.CDABS(CSB)) CS=CSA/CF
           IF (CDABS(CSA).LE.CDABS(CSB)) CS=CSB/CF0
           DO 20 K=0,NM
20            CSJ(K)=CS*CSJ(K)
        ENDIF
        CDJ(0)=(CDCOS(Z)-CDSIN(Z)/Z)/Z
        DO 25 K=1,NM
25         CDJ(K)=CSJ(K-1)-(K+1.0D0)*CSJ(K)/Z
        CSY(0)=-CDCOS(Z)/Z
        CSY(1)=(CSY(0)-CDSIN(Z))/Z
        CDY(0)=(CDSIN(Z)+CDCOS(Z)/Z)/Z
        CDY(1)=(2.0D0*CDY(0)-CDCOS(Z))/Z
        DO 30 K=2,NM
           IF (CDABS(CSJ(K-1)).GT.CDABS(CSJ(K-2))) THEN
              CSY(K)=(CSJ(K)*CSY(K-1)-1.0D0/(Z*Z))/CSJ(K-1)
           ELSE
              CSY(K)=(CSJ(K)*CSY(K-2)-(2.0D0*K-1.0D0)/Z**3)/CSJ(K-2)
           ENDIF
30      CONTINUE
        DO 35 K=2,NM
35         CDY(K)=CSY(K-1)-(K+1.0D0)*CSY(K)/Z
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
