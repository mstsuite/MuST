SUBROUTINE CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY)
!
!  ===================================================================
!       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!                and their derivatives for a complex argument
!       Input :  z --- Complex argument
!                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!       Output:  CSJ(n) --- jn(z)
!                CDJ(n) --- jn'(z)
!                CSY(n) --- yn(z)
!                CDY(n) --- yn'(z)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: N
   integer (kind=IntKind), intent(out) :: NM
   integer (kind=IntKind) :: K, M, MSTA1, MSTA2
   real (kind=RealKind) :: A0
   complex (kind=CmplxKind), intent(in) :: Z
   complex (kind=CmplxKind), intent(out) :: CSJ(0:N),CDJ(0:N),CSY(0:N),CDY(0:N)
   complex (kind=CmplxKind) :: CF, CF0, CF1, CS, CSA, CSB
!
   A0=CDABS(Z)
   NM=N
   IF (A0.LT.1.0D-60) THEN
      DO K=0,N
         CSJ(K)=0.0D0
         CDJ(K)=0.0D0
         CSY(K)=-1.0D+300
         CDY(K)=1.0D+300
      enddo
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
      DO K=M,0,-1
         CF=(2.0D0*K+3.0D0)*CF1/Z-CF0
         IF (K.LE.NM) CSJ(K)=CF
         CF0=CF1
         CF1=CF
      enddo
      IF (CDABS(CSA).GT.CDABS(CSB)) CS=CSA/CF
      IF (CDABS(CSA).LE.CDABS(CSB)) CS=CSB/CF0
      DO K=0,NM
         CSJ(K)=CS*CSJ(K)
      enddo
   ENDIF
   CDJ(0)=(CDCOS(Z)-CDSIN(Z)/Z)/Z
   DO K=1,NM
      CDJ(K)=CSJ(K-1)-(K+1.0D0)*CSJ(K)/Z
   enddo
   CSY(0)=-CDCOS(Z)/Z
   CSY(1)=(CSY(0)-CDSIN(Z))/Z
   CDY(0)=(CDSIN(Z)+CDCOS(Z)/Z)/Z
   CDY(1)=(2.0D0*CDY(0)-CDCOS(Z))/Z
   DO K=2,NM
      IF (CDABS(CSJ(K-1)).GT.CDABS(CSJ(K-2))) THEN
         CSY(K)=(CSJ(K)*CSY(K-1)-1.0D0/(Z*Z))/CSJ(K-1)
      ELSE
         CSY(K)=(CSJ(K)*CSY(K-2)-(2.0D0*K-1.0D0)/Z**3)/CSJ(K-2)
      ENDIF
   enddo
   DO K=2,NM
      CDY(K)=CSY(K-1)-(K+1.0D0)*CSY(K)/Z
   enddo
END SUBROUTINE CSPHJY


FUNCTION MSTA1(X,MP) result(M)
!
!  ===================================================================
!       Purpose: Determine the starting point for backward  
!                recurrence such that the magnitude of    
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point   
!  ===================================================================
!
   use KindParamModule, only : IntKind, RealKind
   implicit none
   integer (kind=IntKind), intent(in) :: MP
   integer (kind=IntKind) :: M
   integer (kind=IntKind) :: N0, N1, NN, IT
   real (kind=RealKind), intent(in) :: X
   real (kind=RealKind) :: A0, F, F0, F1, ENVJ
!
   A0=DABS(X)
   N0=INT(1.1*A0)+1
   F0=ENVJ(N0,A0)-MP
   N1=N0+5
   F1=ENVJ(N1,A0)-MP
   DO IT=1,20             
      NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
      F=ENVJ(NN,A0)-MP
      IF(ABS(NN-N1).LT.1) exit
      N0=N1
      F0=F1
      N1=NN
      F1=F
   enddo
   M=NN
END function MSTA1


FUNCTION MSTA2(X,N,MP) result(M)
!
!  ===================================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!  ===================================================================
!
   use KindParamModule, only : IntKind, RealKind
   implicit none
   integer (kind=IntKind), intent(in) :: N
   integer (kind=IntKind), intent(in) :: MP
   integer (kind=IntKind) :: M
   integer (kind=IntKind) :: N0, N1, NN, IT
   real (kind=RealKind), intent(in) :: X
   real (kind=RealKind) :: A0, EJN, HMP, OBJ, F, F0, F1, ENVJ
!
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
   DO IT=1,20
      NN=N1-(N1-N0)/(1.0D0-F0/F1)
      F=ENVJ(NN,A0)-OBJ
      IF (ABS(NN-N1).LT.1) exit
      N0=N1
      F0=F1
      N1=NN
      F1=F
   enddo
   M=NN+10
END FUNCTION MSTA2


FUNCTION ENVJ(N,X) result(e)
   use KindParamModule, only : IntKind, RealKind
   integer (kind=IntKind), intent(in) :: N
   real (kind=RealKind), intent(in) :: X
   real (kind=RealKind) :: e
!
   E=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
END function ENVJ
