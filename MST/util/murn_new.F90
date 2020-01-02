!  ====================================================================
!
!  modified from murn.f: Yang Wang June 1st, 2002
!  modified from murn2.f : R.S 19.04.91
!  modified from murn1.f : gps 6.12.90
!  least-squares fit of E vs V by Murnaghan equation of state
!
!  ====================================================================
program murn_new
   implicit none
   integer, parameter :: IntKind = kind(1)
   integer, parameter :: RealKind = kind(1.0d0)
!
   integer (kind=IntKind), parameter :: NVAR=4
   integer (kind=IntKind), parameter :: LIM=20
   integer (kind=IntKind), parameter :: MAXITER=75
!
   integer (kind=IntKind) :: in,iout,iunit,aorv,num_av,nnn,i,ierr
!
   real (kind=RealKind), parameter :: CONVYY=13.605826d0 ! conversion factor
                                                         ! from Ryd. to ev
   real (kind=RealKind), parameter :: ULA=0.52917715d0   ! unit of length
   real (kind=RealKind), parameter :: third=1.0d0/3.0d0
!
   real (kind=RealKind) :: convav, conv_inp
   real (kind=RealKind) :: avinp, eninp, av_min, av_max, av
   real (kind=RealKind) :: alatt0, vol0, e0min, b0, b0prim, a0min, der
   real (kind=RealKind) :: alatt, vol, etot, etotry, pr, bm
   real (kind=RealKind) :: e0ev, e0ry, almin, almax, shift
   real (kind=RealKind) :: VOLMIN,VOLMAX,B0MIN,B0MAX,B0PMIN,B0PMAX,FFF
   real (kind=RealKind), allocatable :: EDATA(:),VDATA(:), ADATA(:)
   real (kind=RealKind), allocatable :: X(:),AUX(:)
!
   IN=5
   IOUT=6
   write(iout,*) 'Least-squares fit of etot(v) by Murnaghan eq.'
   write(iout,*)
!
!  ====================================================================
!  READING THE INPUT DATA:
!
!  UNIT OF LENGTH, IN ANGSTROEM, WHICH WILL BE USED THROUGHOUT
!  THIS PROGRAM:
!  ====================================================================
   write(IOUT,'(12X,''UNIT OF LENGTH: ULA ='',F10.6,X,''ANGSTROEM'')')ULA
!
!  ====================================================================
!  UNITS OF ENERGY USED TO INPUT THE DATA:
!      1 means Rydberg, 2 is electron-volts, and 3 is Hartree units
!  ====================================================================
   READ(IN,*)IUNIT
   WRITE(IOUT,'(12X,''UNIT OF ENERGY INPUT: '',$)')
   IF(IUNIT .EQ. 1) then 
      CONV_inp = 13.605826
      WRITE(IOUT,'('' IN RYDBERGS''/1X)')
   else IF(IUNIT .EQ. 2) then 
      CONV_inp = 1.0
      WRITE(IOUT,'('' IN ELECTRONVOLTS''/1X)')
   else IF(IUNIT .EQ. 3) then
      CONV_inp = 2 * 13.605826
      WRITE(IOUT,'('' IN HARTREES''/1X)')
   else
      WRITE(IOUT,'('' ONLY IUNIT = 1 AND 2 ARE RECOGNIZED; STOP.'')')
      STOP
   ENDIF
!
!  ====================================================================
!  Data attribute
!      1: Lattice constant versus Energy
!      2: Volume versus Energy
!  ====================================================================
   read(in,*)aorv
   if (aorv /= 1 .and. aorv /= 2) then
      WRITE(IOUT,'('' ONLY AORV = 1 AND 2 ARE RECOGNIZED; STOP.'')')
      STOP
   endif
!
!  ====================================================================
!  CONVERSION FACTOR FROM A**3 TO THE VOLUME OF UNIT CELL
!  (E.G. 0.25, IF THE STRUCTURE IS SUCH THAT V = A**3/4 ):
!  ====================================================================
   READ(IN,*)CONVAV
   WRITE(IOUT,'('' IN THE STRUCTURE CONSIDERED, THE UNIT CELL VOLUME'',&
   &            '' = A**3 *'',f15.6/x)')CONVAV
!
!  ====================================================================
!  read in alat boundaries and number of points for that fitted Etot
!  should be computed
!  ====================================================================
   read(in,*) av_min, av_max, num_av
   if (num_av < 2 .or. av_max <= av_min) then
      stop 'something wrong with this range'
   else if (aorv == 1) then
      write(iout,'('' Min and Max alat to compute Etot for at'',i4, &
   &               '' points ='',f10.5,'' and'',f10.5)')num_av,av_min,av_max
   else
      write(iout,'('' Min and Max vol to compute Etot for at'',i4, &
   &               '' points ='',f10.5,'' and'',f10.5)')num_av,av_min,av_max
   endif
!
!  ====================================================================
!  NUMBER OF DATA POINTS A, E(A):
!  ====================================================================
   READ(IN,*) nnn
   WRITE(IOUT,'(9X,''NUMBER OF DATA POINTS: N ='',I5)')NNN
   IF (NNN < 3) THEN
      STOP ' N TOO SMALL'
   ENDIF
   allocate(EDATA(nnn),VDATA(nnn),ADATA(nnn))
   allocate(X(10*NVAR),AUX(10*NVAR))
!
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   WRITE(IOUT,'(5X,''I'',5X,''ALATT'',9X,''VOLUME'',8X,''ENERGY(EV)'', &
   &            10X,''ENERGY(RY)'')')
!
!  ====================================================================
!  Read in alat vs e or vol vs e data
!  ====================================================================
   DO I=1,NNN
      READ(IN,*)ADATA(I),VDATA(I),eninp
!     READ(IN,*)avinp,eninp
!     if (aorv == 1) then
!        VDATA(I)=avinp**3*CONVAV
!     else
!        VDATA(I)=avinp
!        avinp = (VDATA(I)/CONVAV)**third
!     endif
      EDATA(I)=eninp*CONV_inp
      WRITE(IOUT,'(1X,I5,F10.5,6x,F10.5,4x,F12.5,8x,f12.5)') &
            I,ADATA(I),VDATA(I),EDATA(I),eninp
   enddo
   WRITE(IOUT,'(1X,79(1H-)/1X)')
!
!  ====================================================================
!  STARTING VALUES FOR ITERATION:
!  ====================================================================
   WRITE(IOUT,'(8X,''ALATT0'',3X, ''VOL0 (ULA**3)'',2X,''B0 (MBAR)'',4X, &
                   ''B0PRIME'',6X,''E0 (EV)'',6X,''E0 (RY)'')')
!
!  ====================================================================
!  MIMIMUM IN EDATA(I):
!  ====================================================================
   E0MIN=1.D6
   VOLMIN=1.D6
   A0MIN=1.D6
   LOOP_I: DO I=1,NNN
      IF(EDATA(I) >= E0MIN) then
         cycle LOOP_I
      endif
      E0MIN=EDATA(I)
      VOLMIN=VDATA(I)
      A0MIN=ADATA(I)
   enddo LOOP_I
!
   VOL0=VOLMIN
   E0EV=E0MIN
   E0RY=E0EV/CONVYY
!  ALATT0=(VOL0/CONVAV)**third
   ALATT0=A0MIN
!
!  ====================================================================
!  FOR OTHER VARIABLES, WE CHOOSE:
!  ====================================================================
   B0=1.D0
   B0PRIM=4.D0
!
   WRITE(IOUT,'(5X,F9.4,3F13.5,2F13.5)')ALATT0,VOL0,B0,B0PRIM,E0EV,E0RY
   WRITE(IOUT,'(64X,''(=INITIAL GUESS)''/1X)')
!
   ALMIN=0.5D0*ALATT0
   ALMAX=1.5D0*ALATT0
!  VOLMIN=ALMIN**3.D0*CONVAV
   call interp(ADATA,VDATA,NNN,ALMIN,VOLMIN,der)
!  VOLMAX=ALMAX**3.D0*CONVAV
   call interp(ADATA,VDATA,NNN,ALMAX,VOLMAX,der)
   B0MIN=0.01D0
   B0MAX=10.D0
   B0PMIN=0.01D0
   B0PMAX=10.D0
   WRITE(IOUT,'('' MIN:'',F9.4,3F13.5/'' MAX:'',F9.4,3F13.5)') &
   &     ALMIN,VOLMIN,B0MIN,B0PMIN,ALMAX,VOLMAX,B0MAX,B0PMAX
!
!  ====================================================================
!  A UNIFORM SHIFT OF ALL THE ENERGIES
!  (WILL NOT INFLUENCE THE B0, B0PRIME, VOL0, ONLY SHIFTS E0):
!  ====================================================================
   SHIFT=-E0EV-1.D0
   DO I=1,NNN
      EDATA(I)=EDATA(I)+SHIFT
   enddo
!
!  ====================================================================
!  MURNAGHAN LEAST-SQUARES FIT:
!  ====================================================================
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   WRITE(IOUT,'('' FIT BY MURNAGHAN EQUATION EQUATION OF STATE'')')
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   WRITE(IOUT,'(1X,''ITER'',2X,''ALATT0'',3X, ''VOL0(ULA**3)'',2X,  &
   &  ''B0(MBAR)'',2X,''B0PRIME'',4X,''E0(EV)'',2X,''SUM OF SQUARES'',&
   &  2x,''IERR''/1X)')
!
   X(1)=E0EV+SHIFT
   X(2)=B0
   X(3)=B0PRIM
   X(4)=VOL0
   IERR = 0
!
   LOOP_i2: DO I = 1, MAXITER
!     CALL DFMND(FUN,X,FFF,NVAR,LIM,AUX,IERR)
      CALL DFMND(X,FFF,NVAR,LIM,AUX,IERR)
      call interp(VDATA,ADATA,NNN,X(4),ALATT0,der)
!     ================================================================
!     THE RESULT OF 20 ITERATION :
!     ================================================================
      WRITE(IOUT,'(1X,I4,F8.4,2x,2F11.5,F10.5,F13.5,1x,F11.5,2x,I5)')  &
   &     20*I,ALATT0,X(4),X(2),X(3),X(1)-SHIFT,FFF,IERR
      IF(IERR .EQ. 0) then
         exit LOOP_i2
      endif
   enddo LOOP_i2
!
!  ====================================================================
!  THE RESULT OF THE LAST ITERATION:
!  ====================================================================
!  call interp(VDATA,ADATA,NNN,X(4),ALATT0,der)
!  WRITE(IOUT,'(1X,I4,F8.4,2x,2F11.5,F10.5,F13.5,1x,F11.5,2x,I5)')  &
!  &  20*I,ALATT0,X(4),X(2),X(3),X(1)-SHIFT,FFF,IERR
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   WRITE(IOUT,'(6X,''ALATT0'',3X, ''VOL0 (ULA**3)'',3X,''B0(MBAR)'',6X, &
   &               ''B0PRIME'',7X,''E0(EV)'',7X,''E0(RY)''/1X)')
   VOL0=X(4)
!  ALATT0=(VOL0/CONVAV)**third
   call interp(VDATA,ADATA,NNN,VOL0,ALATT0,der)
   B0=X(2)
   B0PRIM=X(3)
   E0EV=X(1)-SHIFT
   E0RY=E0EV/CONVYY
   WRITE(IOUT,'(3X,F9.4,X,3F13.5,2F13.5)')ALATT0,VOL0,B0,B0PRIM,E0EV,E0RY
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   write(iout,'(''Alat = '',f12.5,''(Bohr)   = '',f12.5,''(A)'')') &
                alatt0,ula*alatt0
   write(iout,'('' Vol = '',f12.5,''(Bohr^3) = '',f12.5,''(A^3)'')') &
                VOL0,ula**3*VOL0
   write(iout,'('' B0  = '',f12.5,''(MBar)   = '',f12.5,''(GPa)'')') &
                b0,b0*100.0d0
   write(iout,'('' B0p = '',f12.5)')b0prim
   write(iout,'('' E0  = '',f12.5,''(EV)     = '',f12.5,''(Ryd) = '',f12.5,&
   &            ''(Hartree)'')') E0EV,E0RY,E0RY*0.5d0
   WRITE(IOUT,'(1X,79(1H-)/1X)')
!
!  ===================================================================
!  print out fitted energy values for num_av lattice constants or volume 
!  around the equilibrium value given by av_min and av_max cons
!  ===================================================================
   if (aorv == 1) then
      write(iout,'(2x,''alat(A)'',5x,''alat(bohr)'',5x,''energy(ev)'',   &
      &            5x,''energy(Ryd)'',4x,''P(MBar)'',5x,''B(MBar)'')')
   else
      write(iout,'(2x,''vol(A^3)'',5x,''vol(bohr^3)'',4x,''energy(ev)'', &
      &            4x,''energy(Ryd)'',4x,''P(MBar)'',5x,''B(MBar)'')')
   endif
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   do i=1,num_av
      av=av_min + (i-1)*(av_max-av_min)/(num_av-1)
      if (aorv == 1) then
         alatt=av
!        vol = alatt**3 * convav
         call interp(ADATA,VDATA,NNN,alatt,vol,der)
      else
         vol = av
!        alatt = (vol/convav)**third
         call interp(VDATA,ADATA,NNN,vol,alatt,der)
      endif
      call murng1(ula,vol,vol0,b0,b0prim,e0ev,etot,pr,bm)
      etotry=etot/convyy
      if (aorv == 1) then
         write(iout,'(f10.5,4x,f10.5,2x,2(2x,f12.5),1x,2(2x,f10.5))') &
               alatt*ula,alatt,etot,etotry,pr,bm
      else
         write(iout,'(f10.5,4x,f10.5,2x,2(2x,f12.5),1x,2(2x,f10.5))') &
               vol*ula**3,vol,etot,etotry,pr,bm
      endif
   enddo
   deallocate(EDATA,VDATA,ADATA)
   deallocate(X,AUX)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   FUNCTION FUN(X_loc) result(f)
!  ===================================================================
!  FUNCTION TO BE MINIMIZED IN THE LEAST-SQUARES FIT (BY SUBROUTINE DFMND)
!  FUNCTION CALCULATES THE SUM OF THE  A B S O L U T E  DEVIATIONS
!            (E(THEOR)-E(EXP))**2
!  DIVIDED BY THE NUMBER OF EXP. POINTS,
!  ASSUMING FOR EQUATION OF STATE THE MURNAGHAN EXPRESSION.
!
!  MEANING OF VARIABLES:
!      X_loc(1) .... E0
!      X_loc(2) .... B0
!      X_loc(3) .... B0PRIM
!      X_loc(4) .... VOL0
!
!
!  ===================================================================
   IMPLICIT none
   real (kind=RealKind), intent(in) :: X_loc(4)
   real (kind=RealKind) :: f
!
   integer (kind=IntKind) :: i_loc
   real (kind=RealKind) :: E0, B0, B0PRIM, VOL0, SUM, VOLACT
   real (kind=RealKind) :: p,b
!
   E0=X_loc(1)
   B0=X_loc(2)
   B0PRIM=X_loc(3)
   VOL0=X_loc(4)
!
   SUM=0.D0
   DO i_loc=1,NNN
      VOLACT=VDATA(i_loc)
      CALL MURNG1(ULA,VOLACT,VOL0,B0,B0PRIM,E0,ETOT,p,b)
      SUM=SUM+(ETOT-EDATA(i_loc))**2
   enddo
   f=SUM/real(NNN,kind=RealKind)
   end FUNCTION FUN
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   SUBROUTINE MURNG1(ULA,VOL,VOL0,B0,B0PRIM,E0,ETOT,pressure,bmv)
!  ===================================================================
!  EVALUATION OF THE MURNAGHAN EXPRESSION FOR ENERGY AS A FUNCTION
!  OF VOLUME.
!
!  INPUT DATA:
!
!      ULA ..... UNIT OF LENGTH, IN ANGSTROEMS, USED HERE.
!      VOL ..... VOLUME, IN THE ABOVE UNITS OF LENGTH CUBED.
!      VOL0 .... VOLUME AT THE ENERGY MINIMUM, IN THE SAME UNITS.
!      B0 ...... BULK MODULUS, IN UNITS MEGABAR.
!      B0PRIM .. PRESSURE DERIVATIVE OF THE BULK MODULUS.
!                SHOULD BE CLOSE NEITHER TO 0 NOR TO 1.
!      E0 ...... AN ADDITIVE CONSTANT (IN ELECTRONVOLTS), ADDED
!                TO THE ENERGY EXPRESSION.
!                (SEE,  PR B28, p 5484: Fu and Ho)
!
!  OUTPUT DATA:
!
!      ETOT .... THE ENERGY, INCLUDING THE E0 CONSTANT, IN ELECTRONVOLTS
!      pressure: The pressure. It is the 1st derivative of the energy.
!                Unit: MPa.
!      bmv:      The bulk modulus at volume VOL. It is the 2nd derivative 
!                of the energy times VOL. Unit: MPa.
!
!
!  Note:
!      1(ev)/1(bohr^3) = 1.60209*10^5/(0.529177^3) MPa 
!                      = 1.60209/(0.529177^3) MBar
!  IF B0 DIFFERS FROM 0. OR FROM 1. BY LESS THAN 1.d-6, THEN
!  ETOT IS SET AT +111111.111111 ELECTRONVOLTS.
!  ===================================================================
   IMPLICIT none
!
!  ===================================================================
!  CONVERSION FACTOR FROM ERG TO ELECTRONVOLT:
!     1 ELECTRONVOLT = 1.60209 D-12 ERG
!  ===================================================================
   integer, parameter :: IntKind = kind(1)
   integer, parameter :: RealKind = kind(1.0d0)
!
   real (kind=RealKind), PARAMETER :: CONVXX = 1.60209D0
!
   real (kind=RealKind), intent(in) :: ULA,VOL,VOL0,B0,B0PRIM,E0
   real (kind=RealKind), intent(out) :: ETOT
   real (kind=RealKind), intent(out) :: pressure, bmv
!
   IF(ABS(B0PRIM)-1.D0 < 1.d-6 .OR. ABS(B0PRIM) < 1.d-6) THEN
      ETOT=+111111.111111D0
      RETURN
   END IF
!
   IF(VOL < 0.D0 .OR. VOL0 < 0.D0 .OR. VOL0 > 1.0d+10) THEN
      ETOT=+111111.111111D0
      RETURN
   END IF
!
   ETOT = E0  &
         -B0*VOL0/B0PRIM  &
            *(((VOL/VOL0)**(1.D0-B0PRIM)-1.D0)/(1.D0-B0PRIM)-VOL/VOL0+1.D0) &
            *ULA**3/CONVXX
!
   pressure = B0/B0PRIM*((VOL0/VOL)**B0PRIM-1.0d0)
   bmv = B0*(VOL0/VOL)**B0PRIM
!
   END SUBROUTINE MURNG1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  SUBROUTINE DFMND(F,X,Y,N,LIM,AUX,IER)
   SUBROUTINE DFMND(X,Y,N,LIM,AUX,IER)
!  ===================================================================
!
!  *******************************************************************
!  *    MINIMIZATION OF A FUNCTION OF SEVERAL VARIABLES              *
!  *    USING POWELL'S ALGORITHM WITHOUT DERIVATIVES                 *
!  *******************************************************************
!
   implicit none
!
   integer, parameter :: IntKind = kind(1)
   integer, parameter :: RealKind = kind(1.0d0)
!
   integer (kind=IntKind), intent(in) :: N, LIM
   integer (kind=IntKind), intent(inout) :: IER
   integer (kind=IntKind) :: ISW, N1, N2, N3, N4, N5, I, J, K, L, IC, IT
   integer (kind=IntKind) :: M, MF, IS, ID, IW, NS, IP, MK
!
   real (kind=RealKind), intent(inout) :: X(1:)
   real (kind=RealKind), intent(out) :: AUX(1:)
   real (kind=RealKind), intent(out) :: Y
!
!  real (kind=RealKind) :: F,EPS,ETA,TOL
   real (kind=RealKind) :: EPS,ETA,TOL
   real (kind=RealKind) :: DA,DB,DC,DM,DQ,FA,FB,FC,FL,FM,FS,HD,HQ,HX
!
!  interface
!     function F(x_loc) result(f0)
!        integer, parameter :: RealKind = kind(1.0d0)
!        real (kind=RealKind), intent(in) :: x_loc(4)
!        real (kind=RealKind) :: f0
!     end function F
!  end interface
!
!  ===================================================================
!  SUBROUTINES REQUIRED: THE EXTERNAL FUNCTION F.
!
!  INPUT DATA:
!     F .... THE FUNCTION OF N VARIABLES X(1)...X(N) TO BE MINIMIZED
!     X .... X(1) ... X(N) = STARTING GUESS FOR THE VARIABLES;
!            BEWARE: X HAS TO BE DIMENSIONED IN THE MAIN PROGRAM
!            TO CONSIDERABLY MORE THAN N.
!     N .... NUMBER OF VARIABLES; THE DIMENSION OF X AND AUX IN THE
!            CALLING PROGRAM MUST BE, HOWEVER, MUCH HIGHER:
!            - PERHAPS 10 TIMES?
!     LIM .. MAXIMUM NUMBER OF ITERATIONS
!     AUX .. AUXILIARY ARRAY OF THE SAME DIMENSION AS X.
!  OUTPUT DATA:
!     X .... X(1) ... X(N) THE RESULTING MINIMUM
!     Y .... VALUE OF THE FUNCTION AT THE MINIMUM
!     IER .. SOME ERROR INDICATION 
!            IERR=0 MEANS 'CONVERGENCE ACHIEVED'.
!  ===================================================================
!
   ISW  =IER
   IER  =0
!     IF (N) 1,1,2
      IF (N <= 0) THEN
         GOTO 1
      ELSE
         GOTO 2
      ENDIF
    1 IER  =1000
      GOTO 109
!   2 IF (LIM) 3,3,4
    2 IF (LIM <= 0) THEN
         GOTO 3
      ELSE
         GOTO 4
      ENDIF
    3 IER  =2000
      GOTO 109
!
!     SET INITIAL VALUES AND SAVE INITIAL ARGUMENT
!
    4 N1   =N+1
      N2   =N+N
      N3   =N+N2
      N4   =N+N3
      N5   =N*N+N4
      EPS  =1.D-15
      ETA  =N*EPS
      DO K=1,N
         AUX(K)=X(K)
         J    =N3+K
         AUX(J)=1.D0
      enddo
!     FS   =F(AUX)
      FS   =FUN(AUX)
      FB   =FS
      I    =1
      IC   =1
      IT   =0
      M    =N4
      MF   =0
      IS   =1
!
!     START ITERATION CYCLE
!
    6 IT   =IT+1
      FL   =FS
      DO K=N1,N2
         AUX(K)=0.D0
      enddo
      ID   =I
      I    =IC
      IW   =0
!
!     START MINIMIZATION ALONG NEXT DIRECTION
!
    8 NS   =0
      IP   =0
      DB   =0.D0
!     IF (IW) 10,9,10
      IF (IW == 0) THEN
         GOTO 9
      ELSE
         GOTO 10
      ENDIF
    9 HX   =AUX(I)
!  10 IF (IT-1) 11,11,14
   10 IF (IT-1 <= 0) THEN
         GOTO 11
      ELSE
         GOTO 14
      ENDIF
!  11 IF (IS-1) 14,12,14
   11 IF (IS-1 == 0) THEN
         GOTO 12
      ELSE
         GOTO 14
      ENDIF
   12 DM   =.1D0
!     IF (ABS(HX)-1.D0) 38,38,13
      IF (ABS(HX)-1.D0 > 0.0D0) THEN
         GOTO 13
      ELSE
         GOTO 38
      ENDIF
   13 DM   =-DM*HX
      GOTO 38
!  14 IF (IS-2) 18,15,18
   14 IF (IS-2 == 0) THEN
         GOTO 15
      ELSE
         GOTO 18
      ENDIF
!  15 IF (IT-1) 17,16,17
   15 IF (IT-1 == 0) THEN
         GOTO 16
      ELSE
         GOTO 17
      ENDIF
   16 DM   =HQ
      GOTO 38
   17 DM   =DQ
      GOTO 38
!
!     INTERPOLATION USING ESTIMATE OF SECOND DERIVATIVE
!
!  18 IF (IW-1) 20,19,20
   18 IF (IW-1 == 0) THEN
         GOTO 19
      ELSE
         GOTO 20
      ENDIF
   19 J    =N2+I
      GOTO 21
   20 J    =N3+I
   21 HD   =AUX(J)
      DC   =1.D-2
!     IF (IT-2) 23,23,22
      IF (IT-2 <= 0) THEN
         GOTO 23
      ELSE
         GOTO 22
      ENDIF
   22 DC   =HQ
   23 DM   =DC
      MK   =1
      GOTO 51
   24 DM   =DC*HD
!     IF (DM) 26,25,26
      IF (DM == 0.0D0) THEN
         GOTO 25
      ELSE
         GOTO 26
      ENDIF
   25 DM   =1.D0
   26 DM   =.5D0*DC-(FM-FB)/DM
      MK   =2
!     IF (FM-FB) 27,29,29
      IF (FM-FB < 0.0D0) THEN
         GOTO 27
      ELSE
         GOTO 29
      ENDIF
   27 FC   =FB
      FB   =FM
      DB   =DC
!     IF (DM-DB) 28,67,28
      IF (DM-DB == 0.0D0) THEN
         GOTO 67
      ELSE
         GOTO 28
      ENDIF
   28 DC   =0.D0
      GOTO 51
!  29 IF (DM-DB) 31,30,31
   29 IF (DM-DB == 0.0D0) THEN
         GOTO 30
      ELSE
         GOTO 31
      ENDIF
   30 DA   =DC
      FA   =FM
      GOTO 37
   31 FC   =FM
      GOTO 51
!
!     ANALYSE INTERPOLATED FUNCTION VALUE
!
!  32 IF (FM-FB) 34,33,33
   32 IF (FM-FB < 0.0D0) THEN
         GOTO 34
      ELSE
         GOTO 33
      ENDIF
   33 DA   =DM
      FA   =FM
      GOTO 35
   34 DA   =DB
      FA   =FB
      DB   =DM
      FB   =FM
!  35 IF ((DC-DA)/(DB-DA)) 36,36,50
   35 IF ((DC-DA)/(DB-DA) <= 0.0D0) THEN
         GOTO 36
      ELSE
         GOTO 50
      ENDIF
!  36 IF (DB) 67,37,67
   36 IF (DB == 0.0D0) THEN
         GOTO 37
      ELSE
         GOTO 67
      ENDIF
   37 NS   =1
      DM   =-DC
!
!     LINEAR SEARCH FOR SMALLER FUNCTION VALUES
!     ALONG CURRENT DIRECTION
!
!  38 IF (NS-15) 43,43,39
   38 IF (NS-15 <= 0) THEN
         GOTO 43
      ELSE
         GOTO 39
      ENDIF
!  39 IF (FS-FM) 41,40,41
   39 IF (FS-FM == 0.0D0) THEN
         GOTO 40
      ELSE
         GOTO 41
      ENDIF
   40 MF   =N+2
      DB   =0.D0
      GOTO 67
!  41 IF (ABS(DM)-1.D6) 43,43,42
   41 IF (ABS(DM)-1.D6 > 0.0D0) THEN
         GOTO 42
      ELSE
         GOTO 43
      ENDIF
   42 IER  =100
      GOTO 67
   43 NS   =NS+1
      MK   =3
      GOTO 51
!  44 IF (FM-FB) 45,46,47
   44 IF (FM-FB < 0.0D0) THEN
         GOTO 45
      ELSE IF (FM-FB == 0.0D0) THEN
         GOTO 46
      ELSE
         GOTO 47
      ENDIF
   45 DA   =DB
      FA   =FB
      DB   =DM
      FB   =FM
      DM   =DM+DM
      GOTO 38
!  46 IF (FS-FB) 47,45,47
   46 IF (FS-FB == 0.0D0) THEN
         GOTO 45
      ELSE
         GOTO 47
      ENDIF
!  47 IF (NS-1) 48,48,49
   47 IF (NS-1 <= 0) THEN
         GOTO 48
      ELSE
         GOTO 49
      ENDIF
   48 DA   =DM
      FA   =FM
      DM   =-DM
      GOTO 38
   49 DC   =DM
      FC   =FM
!
!     REFINE MINIMUM USING QUADRATIC INTERPOLATION
!
   50 HD   =(FC-FB)/(DC-DB)+(FA-FB)/(DB-DA)
      DM   =.5D0*(DA+DC)+(FA-FC)/(HD+HD)
      IP   =IP+1
      MK   =4
!
!     STEP ARGUMENT VECTOR AND CALCULATE FUNCTION VALUE
!
!  51 IF (IW-1) 54,52,54
   51 IF (IW-1 == 0) THEN
         GOTO 52
      ELSE
         GOTO 54
      ENDIF
   52 DO 53 K=1,N
         L    =M+K
         AUX(K)=X(K)+DM*AUX(L)
   53 CONTINUE
      GOTO 55
   54 AUX(I)=HX+DM
!  55 FM   =F(AUX)
   55 FM   =FUN(AUX)
      GOTO (24,32,44,56),MK
!
!     ANALYSE INTERPOLATED FUNCTION VALUE
!
!  56 IF (FM-FB) 61,61,57
   56 IF (FM-FB > 0.0D0) THEN
         GOTO 57
      ELSE
         GOTO 61
      ENDIF
!  57 IF (IP-3) 58,62,62
   57 IF (IP-3 < 0) THEN
         GOTO 58
      ELSE
         GOTO 62
      ENDIF
!  58 IF ((DC-DB)/(DM-DB)) 60,60,59
   58 IF ((DC-DB)/(DM-DB) > 0.0D0) THEN
         GOTO 59
      ELSE
         GOTO 60
      ENDIF
   59 DC   =DM
      FC   =FM
      GOTO 50
   60 DA   =DM
      FA   =FM
      GOTO 50
   61 DB   =DM
      FB   =FM
!
!     CALCULATE NEW ESTIMATE OF SECOND DERIVATIVE
!     ALONG THE CURRENT DIRECTION
!
   62 HD   =(HD+HD)/(DC-DA)
!     IF (IW-1) 64,63,64
      IF (IW-1 == 0) THEN
         GOTO 63
      ELSE
         GOTO 64
      ENDIF
   63 J    =N2+I
      GOTO 65
   64 J    =N3+I
   65 AUX(J)=HD
!     IF (FB-FS) 67,66,67
      IF (FB-FS == 0.0D0) THEN
         GOTO 66
      ELSE
         GOTO 67
      ENDIF
   66 DB   =0.D0
!
!     SAVE ARGUMENT VECTOR WITH SMALLEST FUNCTION VALUE FOUND
!
!  67 IF (IW-1) 70,68,70
   67 IF (IW-1 == 0) THEN
         GOTO 68
      ELSE
         GOTO 70
      ENDIF
   68 DO K=1,N
         L    =M+K
         J    =N+K
         HD   =DB*AUX(L)
         AUX(J)=AUX(J)+HD
         HD   =X(K)+HD
         AUX(K)=HD
         X(K) =HD
      enddo
      GOTO 71
   70 J    =N+I
      AUX(J)=AUX(J)+DB
      HD   =HX+DB
      AUX(I)=HD
      X(I) =HD
!  71 IF (IER-100) 72,108,72
   71 IF (IER-100 == 0) THEN
         GOTO 108
      ELSE
         GOTO 72
      ENDIF
!
!     DETERMINE DIRECTION FOR NEXT LINEAR SEARCH
!
   72 FS   =FB
!     IF (I-N) 74,73,73
      IF (I-N < 0) THEN
         GOTO 74
      ELSE
         GOTO 73
      ENDIF
   73 I    =0
   74 I    =I+1
!     IF (IS) 75,75,80
      IF (IS <= 0) THEN
         GOTO 75
      ELSE
         GOTO 80
      ENDIF
!  75 IF (DB) 77,76,77
   75 IF (DB == 0.0D0) THEN
         GOTO 76
      ELSE
         GOTO 77
      ENDIF
!  76 IF (I-IC) 8,77,8
   76 IF (I-IC == 0) THEN
         GOTO 77
      ELSE
         GOTO 8
      ENDIF
   77 IC   =I
      IS   =1
!     IF (IT-N) 79,79,78
      IF (IT-N <= 0) THEN
         GOTO 79
      ELSE
         GOTO 78
      ENDIF
   78 IW   =1
   79 I    =ID
      GOTO 8
   80 M    =M+N
!     IF (M-N5) 82,81,81
      IF (M-N5 < 0) THEN
         GOTO 82
      ELSE
         GOTO 81
      ENDIF
   81 M    =N4
!  82 IF (IS-1) 83,83,94
   82 IF (IS-1 <= 0) THEN
         GOTO 83
      ELSE
         GOTO 94
      ENDIF
!  83 IF (I-1) 84,84,85
   83 IF (I-1 <= 0) THEN
         GOTO 84
      ELSE
         GOTO 85
      ENDIF
   84 IW   =1
!  85 IF (I-ID) 8,86,8
   85 IF (I-ID == 0) THEN
         GOTO 86
      ELSE
         GOTO 8
      ENDIF
   86 HQ   =0.D0
      DO 87 K=N1,N2
         HQ   =HQ+AUX(K)*AUX(K)
   87 CONTINUE
!     IF (HQ) 90,88,90
      IF (HQ == 0) THEN
         GOTO 88
      ELSE
         GOTO 90
      ENDIF
!  88 IF (MF-N1) 108,108,89
   88 IF (MF-N1 <= 0) THEN
         GOTO 108
      ELSE
         GOTO 89
      ENDIF
   89 IER  =200
      GOTO 108
   90 DQ   =DSQRT(HQ)
      HQ   =DQ
!     IF (HQ-1.D0) 92,92,91
      IF (HQ-1.D0 <= 0.0D0) then
         goto 92
      else
         goto 91
      endif
   91 HQ   =1.D0
   92 DO 93 K=N1,N2
         L    =M+K-N
         AUX(L)=AUX(K)/DQ
   93 CONTINUE
      IS   =2
      GOTO 8
!
!     END OF ITERATION CYCLE
!     TEST FOR TERMINATION OF MINIMIZATION
!
   94 IS   =0
      TOL  =EPS
!     IF (ABS(FS)-1.D0) 96,96,95
      IF (ABS(FS)-1.D0 <= 0.0D0) then
         goto 96
      ELSE
         goto 95
      ENDIF
   95 TOL  =EPS*ABS(FS)
!  96 IF (FL-FS-TOL) 100,100,97
   96 IF (FL-FS-TOL <= 0.0D0) THEN
         goto 100
      ELSE
         goto 97
      ENDIF
!  97 IF (IT-LIM) 99,99,98
   97 IF (IT-LIM <= 0) then
         goto 99
      ELSE
         goto 98
      ENDIF
   98 IER  =10
      GOTO 108
   99 MF   =0
      GOTO 6
! 100 IF (MF-N1) 102,102,101
  100 IF (MF-N1 <= 0) then
         goto 102
      ELSE
         goto 101
      ENDIF
  101 IER  =200
      GOTO 108
  102 MF   =MF+1
      DQ   =0.D0
      DO 105 K=1,N
         J    =N+K
!        IF (ABS(AUX(K))-1.D0) 103,103,104
         IF (ABS(AUX(K))-1.D0 <= 0.0D0) then
            goto 103
         ELSE
            goto 104
         ENDIF
  103    DQ   =DQ+ABS(AUX(J))
         GOTO 105
  104    DQ   =DQ+ABS(AUX(J)/AUX(K))
  105    CONTINUE
!     IF (DQ-ETA) 108,108,106
      IF (DQ-ETA <= 0.0D0) then
         goto 108
      else
         goto 106
      endif
! 106 IF (MF-N) 6,6,107
  106 IF (MF-N <= 0) then
         goto 6
      else
         goto 107
      endif
  107 IER  =1
  108 Y    =FB
!     IF (IER) 111,111,109
      IF (IER <= 0) then
         goto 111
      else
         goto 109
      endif
! 109 IF (ISW+12345) 110,111,110
  109 IF (ISW+12345 == 0) then
         goto 111
      ELSE
         goto 110
      ENDIF
! 110 CALL WIER(IER,20212)
  110 CONTINUE
  111 RETURN
  END SUBROUTINE DFMND
!
end program murn_new
