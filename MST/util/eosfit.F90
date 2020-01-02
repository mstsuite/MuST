!  ====================================================================
!  This utility code is intended to fit the following Equation of State
!       1. Murnaghan equation
!       2. Birch-Murnaghan equation
!       3. Vinet equation
!       4. Poirier-Tarantola equation
!  given a set of total energy versus the lattice constant data
!
!  -- By Yang Wang @ PSC, Oct 29, 2012
!
!  This code is based on a code originally called murn2.f
!      modified from murn.f: Yang Wang June 1st, 2002
!      modified from murn2.f : R.S 19.04.91
!      modified from murn1.f : gps 6.12.90
!      least-squares fit of E vs V by Murnaghan equation of state
!
!  ====================================================================
program eosfit
   implicit none
   integer, parameter :: IntKind = kind(1)
   integer, parameter :: RealKind = kind(1.0d0)
!
   character (len=50) :: text
!
   integer (kind=IntKind), parameter :: LIM=20
   integer (kind=IntKind), parameter :: MAXITER=75
!
   integer (kind=IntKind) :: in, iout, ierr, status, NVAR, I0
   integer (kind=IntKind) :: lunit, eunit, latparam, num_av, nnn, i, na, eos
!
   real (kind=RealKind), parameter :: CONVYY=13.605826d0 ! conversion factor
                                                         ! from Ryd. to ev
   real (kind=RealKind), parameter :: ULA=0.52917715d0   ! unit of length
   real (kind=RealKind), parameter :: third=1.0d0/3.0d0
!
   real (kind=RealKind) :: convav, ECONV, LCONV
   real (kind=RealKind) :: eninp, av_min, av_max, av
   real (kind=RealKind) :: ALATT0, VOL0, b0, b0prim
   real (kind=RealKind) :: alatt, etot, etotry, pr, bm
   real (kind=RealKind) :: e0ev, e0ry, shift
   real (kind=RealKind) :: FFF,VOL
   real (kind=RealKind), allocatable :: EDATA(:),VDATA(:), pinp(:,:)
   real (kind=RealKind) :: X(100), AUX(100)
!
   IN=5
   IOUT=6
   write(iout,'(/,8x,a)')'A Utility for Fitting the Equation of State'
   write(iout,'(  8x,a)')'==========================================='
   write(iout,'(/)')
!
!  ====================================================================
!  READING THE INPUT DATA:
!
!  UNIT OF LENGTH, IN ANGSTROEM, WHICH WILL BE USED THROUGHOUT
!  THIS PROGRAM.
!  ====================================================================
!  write(IOUT,'(8X,''UNIT OF LENGTH: ULA ='',F10.6,X,''ANGSTROEM'')')ULA
   lunit = 0
   do while (lunit < 1 .or. lunit > 2)
      write(IOUT,'( 8X,a)')'Options for the units of length:'
      write(IOUT,'(12X,a)')'1. Bohr radius'
      write(IOUT,'(12X,a)')'2. Angstrom'
      write(IOUT,'( 8X,a,$)')'Enter the input >>> '
      READ(IN,*)lunit
      IF(lunit .EQ. 1) then 
         LCONV = ULA
         WRITE(IOUT,'('' IN Bohr Radius'')')
      else IF(lunit .EQ. 2) then 
         LCONV = 1.0d0
         WRITE(IOUT,'('' IN ANGSTROMS'')')
      else
         WRITE(IOUT,'('' ONLY 1 or 2 ARE RECOGNIZED.'')')
      ENDIF
   enddo
   write(IOUT,'(/)')
!
!  ====================================================================
!  UNITS OF ENERGY USED TO INPUT THE DATA:
!      1 means Rydberg, 2 is electron-volts, and 3 is Hartree units
!  ====================================================================
   eunit = 0
   do while (eunit < 1 .or. eunit > 3)
      write(IOUT,'( 8X,a)')'Options for the units of energy:'
      write(IOUT,'(12X,a)')'1. Rydberg'
      write(IOUT,'(12X,a)')'2. eV'
      write(IOUT,'(12X,a)')'3. Hartree'
      write(IOUT,'( 8X,a,$)')'Enter the input >>> '
      READ(IN,*)eunit
      IF(eunit .EQ. 1) then 
         ECONV = CONVYY
         WRITE(IOUT,'('' IN RYDBERGS'')')
      else IF(eunit .EQ. 2) then 
         ECONV = 1.0d0
         WRITE(IOUT,'('' IN ELECTRONVOLTS'')')
      else IF(eunit .EQ. 3) then
         ECONV = 2.0d0 * CONVYY
         WRITE(IOUT,'('' IN HARTREES'')')
      else
         WRITE(IOUT,'('' ONLY 1, 2, or 3 ARE RECOGNIZED.'')')
      endif
   enddo
   write(IOUT,'(/)')
!
!  ====================================================================
!  Data attribute
!      1-9: Lattice constant versus Energy
!      10: Volume versus Energy
!  ====================================================================
   latparam = 0
   do while (latparam < 1 .or. latparam > 10)
      write(IOUT,'( 8X,a)')'Options for the input data form:'
      write(IOUT,'(12X,a)')'1. a0 versus Energy/Atom for CUBIC system'
      write(IOUT,'(12X,a)')'2. (a0, c0) versus Energy/Atom for HCP structure'
      write(IOUT,'(12X,a)')'3. (a0, c0/a0) versus Energy/Atom for HCP structure'
      write(IOUT,'(12X,a)')'4. (a0, c0) versus Energy/Atom for TETRAGONAL system'
      write(IOUT,'(12X,a)')'5. (a0, c0/a0) versus Energy/Atom for TETRAGONAL system'
      write(IOUT,'(12X,a)')'6. (a0, b0, c0) versus Energy/Atom for ORTHORHOMBIC system'
      write(IOUT,'(12X,a)')'7. (a0, b0/a0, c0/a0) versus Energy/Atom for ORTHORHOMBIC system'
      write(IOUT,'(12X,a)')'8. (a0, b0, c0, degree, degree, degree) versus Energy/Atom'
      write(IOUT,'(12X,a)')'9. (a0, b0/a0, c0/a0, degree, degree, degree) versus Energy/Atom'
      write(IOUT,'(12X,a)')'10. Unit Cell Volume versus Energy/Atom for CUBIC system'
      write(IOUT,'( 8X,a,$)')'Enter the input >>> '
      read(in,*)latparam
      if (latparam < 1 .and. latparam > 10) then
         WRITE(IOUT,'('' ONLY 1 - 10 ARE RECOGNIZED.'')')
      else
         WRITE(IOUT,'('' The option entered is'',i2)')latparam
      endif
   enddo
   write(IOUT,'(/)')
!
!  ====================================================================
!  CONVERSION FACTOR FROM A**3 TO THE ATOMIC VOLUME
!  (E.G. 0.25, IF THE STRUCTURE IS SUCH THAT V = A**3/4 ):
!  ====================================================================
   na = 0
   do while (na < 1) 
      write(IOUT,'( 8X,a,$)')'Enter the number of atoms per unit cell >>> '
      READ(IN,*)na
      if (na < 1) then
         WRITE(IOUT,'('' ONLY VALUE = 1, 2, ...,etc IS RECOGNIZED.'')')
      else
         WRITE(IOUT,'('' The number of atoms per unit cell ='',i5)')na
      endif
   enddo
   write(IOUT,'(/)')
   CONVAV = 1.0d0/real(na,kind=RealKind)
   WRITE(IOUT,'(8x,''IN THE STRUCTURE CONSIDERED, THE ATOMIC VOLUME'', &
   &            '' = Unit Cell Volume *'',f15.6,/)')CONVAV
!
!  ====================================================================
!  Equation of State for fitting the data.             
!       1. Murnaghan equation
!       2. Birch-Murnaghan equation
!       3. Vinet equation
!       4. Poirer-Tarantola equation
!  ====================================================================
   eos = 0
   do while (eos < 1 .or. eos > 4)
      write(IOUT,'( 8X,a)')'Options for the Equation of State:'
      write(IOUT,'(12X,a)')'1. Murnaghan equation'
      write(IOUT,'(12X,a)')'2. Birch-Murnaghan equation'
      write(IOUT,'(12X,a)')'3. Vinet equation'
      write(IOUT,'(12X,a)')'4. Poirer-Tarantola equation'
      write(IOUT,'( 8X,a,$)')'Enter the input >>> '
      read(in,*)eos
      if (eos == 1) then
         WRITE(IOUT,'(8x,a)')'Murnaghan equation is used to fit the data'
      else if (eos == 2) then
         WRITE(IOUT,'(8x,a)')'Birch-Murnaghan equation is used to fit the data'
      else if (eos == 3) then
         WRITE(IOUT,'(8x,a)')'Vinet equation is used to fit the data'
      else if (eos == 4) then
         WRITE(IOUT,'(8x,a)')'Poirer-Tarantola equation is used to fit the data'
      else
         WRITE(IOUT,'('' ONLY 1 - 4 ARE RECOGNIZED.'')')
      endif
   enddo
   write(IOUT,'(/)')
!
!  ====================================================================
!  read in alat boundaries and number of points for that fitted Etot
!  should be computed
!  ====================================================================
   num_av = 0; av_min = 0.0d0; av_max = -1.0d0
   do while (num_av < 1 .or. av_min < 1.0d-6 .or. av_max < av_min)
      if (latparam == 1) then
         write(iout,'( 8X,a,$)')'Enter lattice constant range: a0_min, a0_max, and num_a0 for plotting >>> '
      else
         write(iout,'( 8X,a,$)')'Enter unit cell volume range: vol_min, vol_max, and num_vol for plotting >>> '
      endif
      read(in,*) av_min, av_max, num_av
      if (num_av < 1 .or. av_max < av_min .or. av_min < 1.0d-6) then
         WRITE(IOUT,'('' INVALID LATTICE CONSTANT/VOLUME RANGE.'')')
      endif
   enddo
   if (latparam == 1) then
      write(iout,'('' Min and Max alat to compute Etot for at'',i4, &
   &               '' points ='',f10.5,'' and'',f10.5)')num_av,av_min,av_max
   else
      write(iout,'('' Min and Max vol to compute Etot for at'',i4, &
   &               '' points ='',f10.5,'' and'',f10.5)')num_av,av_min,av_max
   endif
   write(IOUT,'(/)')
!
!  ====================================================================
!  NUMBER OF DATA POINTS A, E(A):
!  ====================================================================
   nnn = 0
   do while (nnn < 3)
      write(iout,'( 8X,a,$)')'Enter the number of (a vs E) or (vol vs E) data >>> '
      READ(IN,*) nnn
      if (nnn < 3) then
         WRITE(IOUT,'( 8X,''NUMBER OF DATA POINTS SHOULD BE GREATER THAN 2'')')
      endif
   enddo
   WRITE(IOUT,'(9X,''NUMBER OF DATA POINTS: N ='',I5)')NNN
   write(IOUT,'(/)')
!
   allocate(EDATA(nnn),VDATA(nnn),pinp(6,nnn))
!
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   WRITE(IOUT,'(5X,''I'',7X,''VOLUME(Angstrom^3)'',8X,''ENERGY(EV)'')')
!
   if (latparam == 1) then
      av_min = av_min*LCONV
      av_max = av_max*LCONV
   else
      av_min = av_min*LCONV**3*CONVAV
      av_max = av_max*LCONV**3*CONVAV
   endif
!
!  ====================================================================
!  Read in alat vs e or vol vs e data
!  ====================================================================
   pinp = 0.0d0
   DO I=1,NNN
      read(IN,'(a)')text
      if (latparam == 1 .or. latparam == 10) then
         READ(text,*,iostat=status)pinp(1,I),eninp
         pinp(1,I) = pinp(1,I)*LCONV
      else if (latparam == 2) then
         READ(text,*,iostat=status)pinp(1:2,I),eninp
         pinp(1:2,I) = pinp(1:2,I)*LCONV
      else if (latparam == 3) then
         READ(text,*,iostat=status)pinp(1:2,I),eninp
         pinp(1,I) = pinp(1,I)*LCONV
         pinp(2,I) = pinp(1,I)*pinp(2,I)
      else if (latparam == 4) then
         READ(text,*,iostat=status)pinp(1:2,I),eninp
         pinp(1:2,I) = pinp(1:2,I)*LCONV
      else if (latparam == 5) then
         READ(text,*,iostat=status)pinp(1:2,I),eninp
         pinp(1,I) = pinp(1,I)*LCONV
         pinp(2,I) = pinp(1,I)*pinp(2,I)
      else if (latparam == 6) then
         READ(text,*,iostat=status)pinp(1:3,I),eninp
         pinp(1:3,I) = pinp(1:3,I)*LCONV
      else if (latparam == 7) then
         READ(text,*,iostat=status)pinp(1:3,I),eninp
         pinp(1,I) = pinp(1,I)*LCONV
         pinp(2,I) = pinp(1,I)*pinp(2,I)
         pinp(3,I) = pinp(1,I)*pinp(3,I)
      else if (latparam == 8) then
         READ(text,*,iostat=status)pinp(1:6,I),eninp
         pinp(1:3,I) = pinp(1:3,I)*LCONV
      else if (latparam == 9) then
         READ(text,*,iostat=status)pinp(1:6,I),eninp
         pinp(1,I) = pinp(1,I)*LCONV
         pinp(2,I) = pinp(1,I)*pinp(2,I)
         pinp(3,I) = pinp(1,I)*pinp(3,I)
      endif
      if (status /= 0) then
         write(IOUT,'(/,a,a,a,i2)')'This data format: ',trim(text),  &
   &        ' is inconsistent with data format option = ',latparam
         stop 'Input Error Encountered.'
      endif
!     ----------------------------------------------------------------
      VDATA(I) = getCellVolume(pinp(1:6,I))*CONVAV ! This is volume/atom Angstrom units.
!     ----------------------------------------------------------------
      EDATA(I)=eninp*ECONV
      WRITE(IOUT,'(1X,I5,8x,F12.5,12x,f12.5)') I,VDATA(I),EDATA(I)
   enddo
   WRITE(IOUT,'(1X,79(1H-)/1X)')
!
!  ====================================================================
!  Determine the STARTING VALUES FOR ITERATION.
!
!  MIMIMUM IN EDATA(I):
!  ====================================================================
   E0EV=1.D6
   VOL0=1.D6
   I0 = 1
   DO I=1,NNN
      IF(EDATA(I) < E0EV) then
         E0EV=EDATA(I)
         VOL0=VDATA(I)
         I0 = I
      endif
   enddo
   E0RY=E0EV/CONVYY
!
!  ====================================================================
!  FOR OTHER VARIABLES, WE CHOOSE:
!  ====================================================================
   B0=1.D0
   B0PRIM=4.D0
!
   if (latparam == 1 .or. latparam == 10) then
      ALATT0=(VOL0*na)**third
      WRITE(IOUT,'(8X,''A0(Ang)'',3X, ''VOL0 (Ang**3)'',4X,''B0 (MBAR)'',4X, &
                      ''B0PRIME'',6X,''E0 (EV)'',6X,''E0 (RY)'')')
      WRITE(IOUT,'(5X,F9.4,2x,3F13.5,2x,2F13.5)')ALATT0,VOL0,B0,B0PRIM,E0EV,E0RY
   else
      WRITE(IOUT,'(8X,''VOL0 (Ang**3)'',2X,''B0 (MBAR)'',4X, &
                      ''B0PRIME'',6X,''E0 (EV)'',6X,''E0 (RY)'')')
      WRITE(IOUT,'(5X,3F13.5,2F13.5)')VOL0,B0,B0PRIM,E0EV,E0RY
   endif
   WRITE(IOUT,'(64X,''(=INITIAL GUESS)''/1X)')
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
   if (latparam == 1 .or. latparam == 10) then
      NVAR = 4
   else if (latparam >= 2 .and. latparam <= 5) then
      NVAR = 5
   else if (latparam == 6 .or. latparam == 7) then
      NVAR = 6
   else if (latparam == 8 .or. latparam == 9) then
      NVAR = 9
   endif 
!
!  *******************************************************************
!  Here is the previous statement:
!       allocate(X(10*NVAR),AUX(10*NVAR))
!  Now, the dimension of these two arrays is set to be 100, which should
!  be sufficient.
!  *******************************************************************
   X(:) = 0.0d0; AUX(:) = 0.0d0
!
   X(1)=E0EV+SHIFT
   X(2)=B0
   X(3)=B0PRIM
   do I = 4, NVAR
      X(I) = pinp(I-3,I0)
   enddo
   IERR = 0
!
!  ====================================================================
!  Eqation of State LEAST-SQUARES FIT:
!  ====================================================================
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   if (eos == 1) then
      WRITE(IOUT,'('' FIT BY MURNAGHAN EQUATION EQUATION OF STATE'')')
   else if (eos == 2) then
      WRITE(IOUT,'('' FIT BY BIRCH-MURNAGHAN EQUATION EQUATION OF STATE'')')
   else if (eos == 3) then
      WRITE(IOUT,'('' FIT BY VINET EQUATION EQUATION OF STATE'')')
   else if (eos == 4) then
      WRITE(IOUT,'('' FIT BY Poirier-Tarantola EQUATION EQUATION OF STATE'')')
   endif
   WRITE(IOUT,'(1X,79(1H-)/1X)')
   if (latparam == 1 .or. latparam == 10) then
      WRITE(IOUT,'(1X,''ITER'',2X,''A(Ang)'',3X, ''VOL0(Ang**3)'',3X,    &
   &     ''B0(MBAR)'',2X,''B0PRIME'',4X,''E0(EV)'',4X,''SUM OF SQUARES'',&
   &     2x,''IERR'')')
   else if (latparam >= 2 .and. latparam <= 5) then
      WRITE(IOUT,'(1X,''ITER'',2X,''A(Ang)'',4X,''C(Ang)'',3X,           &
   &                  ''VOL0(Ang**3)'',3X,                               &
   &     ''B0(MBAR)'',2X,''B0PRIME'',4X,''E0(EV)'',4X,''SUM OF SQUARES'',&
   &     2x,''IERR'')')
   else if (latparam == 6 .and. latparam == 7) then
      WRITE(IOUT,'(1X,''ITER'',2X,''A(Ang)'',4X,''B(Ang)'',4X,''C(Ang)'',&
   &               3X, ''VOL0(Ang**3)'',3X,                              &
   &     ''B0(MBAR)'',2X,''B0PRIME'',4X,''E0(EV)'',4X,''SUM OF SQUARES'',&
   &     2x,''IERR'')')
   else if (latparam == 8 .and. latparam == 9) then
      WRITE(IOUT,'(1X,''ITER'',2X,''A(Ang)'',4X,''B(Ang)'',4X,''C(Ang)'',&
   &                  ''T1(Deg)'',5X,''T2(Deg)'',5X,''T3(Deg)'',3X,      &
   &                  ''VOL0(Ang**3)'',3X,                               &
   &     ''B0(MBAR)'',2X,''B0PRIME'',4X,''E0(EV)'',4X,''SUM OF SQUARES'',&
   &     2x,''IERR'')')
   endif
!
   LOOP_i: DO I = 1, MAXITER
!     ----------------------------------------------------------------
      CALL DFMND(X,FFF,NVAR,LIM,AUX,IERR)
!     ----------------------------------------------------------------
      VOL = getCellVolume(X(4:9))*CONVAV
      if (latparam == 1) then
         WRITE(IOUT,'(1X,I4,F8.4,2x,2F11.5,F10.5,F13.5,1x,F11.5,2x,I5)') &
   &        20*I,X(4),VOL,X(2),X(3),X(1)-SHIFT,FFF,IERR
      else if (latparam >= 2 .and. latparam <= 5) then
         WRITE(IOUT,'(1X,I4,2(F8.4,2X),X,2F11.5,F10.5,F13.5,1x,F11.5,3x,I5)') &
   &        20*I,X(4),X(5),VOL,X(2),X(3),X(1)-SHIFT,FFF,IERR
      else if (latparam == 6 .and. latparam == 7) then
         WRITE(IOUT,'(1X,I4,3(F8.4,2X),X,2F11.5,F10.5,F13.5,1x,F11.5,3x,I5)') &
   &        20*I,X(4:6),VOL,X(2),X(3),X(1)-SHIFT,FFF,IERR
      else if (latparam == 8 .and. latparam == 9) then
         WRITE(IOUT,'(1X,I4,6(F8.4,2X),X,2F11.5,F10.5,F13.5,1x,F11.5,3x,I5)') &
   &        20*I,X(4:9),VOL,X(2),X(3),X(1)-SHIFT,FFF,IERR
      else if (latparam == 10) then
         WRITE(IOUT,'(1X,I4,F8.4,2x,2F11.5,F10.5,F13.5,1x,F11.5,3x,I5)') &
   &        20*I,(VOL*na)**third,VOL,X(2),X(3),X(1)-SHIFT,FFF,IERR
      endif
!
      IF(IERR .EQ. 0) then
         exit LOOP_i
      endif
   enddo LOOP_i
   WRITE(IOUT,'(1X,79(1H-)/1X)')
!
!  ====================================================================
!  THE RESULT OF THE LAST ITERATION:
!  ====================================================================
   VOL0 = getCellVolume(X(4:9))*CONVAV
   B0=X(2)
   B0PRIM=X(3)
   E0EV=X(1)-SHIFT
   E0RY=E0EV/CONVYY
   if (latparam == 1) then
      write(iout,'(''Alat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(4)/ULA,X(4)
   else if (latparam >= 2 .or. latparam <= 5) then
      write(iout,'(''Alat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(4)/ULA,X(4)
      write(iout,'(''Clat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(5)/ULA,X(5)
   else if (latparam == 6 .or. latparam == 7) then
      write(iout,'(''Alat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(4)/ULA,X(4)
      write(iout,'(''Blat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(5)/ULA,X(5)
      write(iout,'(''Clat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(6)/ULA,X(6)
   else if (latparam == 8 .or. latparam == 9) then
      write(iout,'(''Alat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(4)/ULA,X(4)
      write(iout,'(''Blat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(5)/ULA,X(5)
      write(iout,'(''Clat = '',f12.5,''( Bohr ) = '',f12.5,''( A )'')') X(6)/ULA,X(6)
      write(iout,'(''Alpha = '',f12.5,''(Degree) '')') X(7)
      write(iout,'(''Beta  = '',f12.5,''(Degree) '')') X(8)
      write(iout,'(''Gamma = '',f12.5,''(Degree) '')') X(9)
   endif
   write(iout,'('' Vol = '',f12.5,''(Bohr^3) = '',f12.5,''(A^3)'')') &
                VOL0/ULA**3,VOL0
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
   if (latparam == 1) then
      write(iout,'(3x,''alat(A)'',5x,''alat(bohr)'',4x,''energy(ev)'',                   &
      &            5x,''energy(Ryd)'',5x,''P(MBar)'',5x,''B(MBar)'')')
!  else if (latparam >= 2 .or. latparam <= 5) then
!     write(iout,'(3x,''alat(A)'',5x,''alat(bohr)'',5x,''clat(A)'',5x,''clat(bohr)'',5x, &
!     write(iout,'(3x,''alat(A)'',5x,''alat(bohr)'',5x,''clat(A)'',5x,''clat(bohr)'',5x, &
!     &               ''energy(ev)'',5x,''energy(Ryd)'',4x,''P(MBar)'',5x,''B(MBar)'')')
!  else if (latparam == 6 .or. latparam == 7) then
!     write(iout,'(3x,''alat(A)'',5x,''alat(bohr)'',5x,''blat(A)'',5x,''blat(bohr)'',5x, &
!     &               ''clat(A)'',5x,''clat(bohr)'',5x,                                  &
!     &               ''energy(ev)'',5x,''energy(Ryd)'',5x,''P(MBar)'',5x,''B(MBar)'')')
!  else if (latparam == 8 .or. latparam == 9) then
!     write(iout,'(3x,''alat(A)'',5x,''alat(bohr)'',5x,''blat(A)'',5x,''blat(bohr)'',5x, &
!     &               ''clat(A)'',5x,''clat(bohr)'',5x,                                  &
!     &               ''t1(deg)'',5x,''t2(deg)'',5x,''t3(deg)'',5x,                      &
!     &               ''energy(ev)'',5x,''energy(Ryd)'',5x,''P(MBar)'',5x,''B(MBar)'')')
!  else if (latparam == 10) then
   else
      write(iout,'(2x,''vol(A^3)'',4x,''vol(bohr^3)'',4x,''energy(ev)'', &
      &            4x,''energy(Ryd)'',6x,''P(MBar)'',5x,''B(MBar)'')')
   endif
   WRITE(IOUT,'(1X,79(1H-))')
   do i=1,num_av
      av=av_min + (i-1)*(av_max-av_min)/(num_av-1)
      if (latparam == 1) then
         VOL = av**3 * convav
      else
         VOL = av
      endif
!     ----------------------------------------------------------------
      call eosfun(VOL,VOL0,b0,b0prim,e0ev,etot,pr,bm)
!     ----------------------------------------------------------------
      etotry=etot/CONVYY
      if (latparam == 1) then
         write(iout,'(f10.5,4x,f10.5,2x,2(2x,f12.5),1x,2(2x,f10.5))') &
               av,av/ULA,etot,etotry,pr,bm
      else
         write(iout,'(f10.5,4x,f10.5,2x,2(2x,f12.5),1x,2(2x,f10.5))') &
               VOL,VOL/ULA**3,etot,etotry,pr,bm
      endif
   enddo
   deallocate(EDATA,VDATA,pinp)
!  deallocate(X,AUX)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   FUNCTION getCellVolume(param) result(v)
!  ===================================================================
   IMPLICIT none
!
   real (kind=RealKind), intent(in) :: param(6)
   real (kind=RealKind) :: v, r1, r2, r3
   real (kind=RealKind), parameter :: d2r=3.14159265358979d0/180.0d0
!
   if (latparam == 1) then
      v=param(1)**3
   else if (latparam == 2 .or. latparam == 3) then
      v=0.5d0*sqrt(3.0d0)*param(1)**2*param(2)
   else if (latparam == 4 .or. latparam == 5) then
      v=param(1)**2*param(2)
   else if (latparam == 6 .or. latparam == 7) then
      v=param(1)*param(2)*param(3)
   else if (latparam == 8 .or. latparam == 9) then
      r1 = param(4)*d2r
      r2 = param(5)*d2r
      r3 = param(6)*d2r
      v=param(1)*param(2)*param(3)                                    &
                *sqrt( 1.0d0-cos(r1)*cos(r1)                          &
                            -cos(r2)*cos(r2)                          &
                            -cos(r3)*cos(r3)                          &
                            +2.0d0*cos(r1)*cos(r2)*cos(r3) )
   else if (latparam == 10) then
      v=param(1)
   endif
!
   end FUNCTION getCellVolume
!  ===================================================================
!
!  *******************************************************************
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
!      X_loc(4:9) .... Lattice parameters
!
!
!  ===================================================================
   IMPLICIT none
!
   real (kind=RealKind), intent(in) :: X_loc(9)
   real (kind=RealKind) :: f
!
   integer (kind=IntKind) :: i_loc
   real (kind=RealKind) :: E0, B0, B0PRIM, VOL0, SUM, VOLACT
   real (kind=RealKind) :: p,b
!
   E0=X_loc(1)
   B0=X_loc(2)
   B0PRIM=X_loc(3)
!  VOL0=X_loc(4)
   VOL0=getCellVolume(X_loc(4:9))/real(na,kind=RealKind)
!
   SUM=0.D0
   DO i_loc=1,NNN
      VOLACT=VDATA(i_loc)
      CALL eosfun(VOLACT,VOL0,B0,B0PRIM,E0,ETOT,p,b)
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
   SUBROUTINE eosfun(VOL,VOL0,B0,B0PRIM,E0,ETOT,pressure,bmv)
!  ===================================================================
!  EVALUATION OF THE MURNAGHAN EXPRESSION FOR ENERGY AS A FUNCTION
!  OF VOLUME.
!
!  INPUT DATA:
!
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
   real (kind=RealKind), PARAMETER :: CONVXX = 1.60209D0
!
   real (kind=RealKind), intent(in) :: VOL,VOL0,B0,B0PRIM,E0
   real (kind=RealKind), intent(out) :: ETOT
   real (kind=RealKind), intent(out) :: pressure, bmv
   real (kind=RealKind) :: fac, fac2, r
!
   IF(ABS(B0PRIM-1.D0) < 1.0d-6 .AND. (eos == 1 .or. eos == 3)) THEN
      ETOT=+111111.111111D0
      RETURN
   else IF(ABS(B0PRIM) < 1.0d-6) THEN
      ETOT=+111111.111111D0
      RETURN
   END IF
!
   IF(VOL < 1.0D-6 .OR. VOL0 < 1.0d-6 .OR. VOL0 > 1.0d+8) THEN
      ETOT=+111111.111111D0
      RETURN
   END IF
!
   if (eos == 1) then
      ETOT = E0-B0*VOL0/B0PRIM*( ((VOL/VOL0)**(1.D0-B0PRIM)-1.D0)/(1.D0-B0PRIM) &
                                -VOL/VOL0+1.D0 )/CONVXX
!                               -VOL/VOL0+1.D0 )*ULA**3/CONVXX
      pressure = B0/B0PRIM*( (VOL0/VOL)**B0PRIM-1.0d0 )
   else if (eos == 2) then
      fac = (VOL0/VOL)**(1.0d0/3.0d0)
      fac2 = (VOL0/VOL)**(2.0d0/3.0d0)
      ETOT = E0 + 0.5625d0*VOL0*B0*(fac2-1.0d0)**2*((fac2-1.0d0)*B0PRIM+6.0d0-4.0d0*fac2)/CONVXX
!     ETOT = E0 + 0.5625d0*VOL0*B0*(fac2-1.0d0)**2*((fac2-1.0d0)*B0PRIM+6.0d0-4.0d0*fac2)*ULA**3/CONVXX
      pressure = 1.5d0*B0*fac**5*(fac2-1.0d0)*(1.0d0+0.75d0*(B0PRIM-4.0d0)*(fac2-1.0d0))
   else if (eos == 3) then
      fac = (VOL0/VOL)**(1.0d0/3.0d0)
      fac2 = (VOL0/VOL)**(2.0d0/3.0d0)
      ETOT = E0 + 2.0d0*VOL0*B0/(B0PRIM-1.0d0)**2                    &
                       *( 2.0d0-(5.0d0+3.0d0*B0PRIM*(fac-1.0d0)-3.0d0*fac)*exp(-1.5d0*(B0PRIM-1.0d0)*(fac-1.0d0)) ) &
                       /CONVXX
!                      *ULA**3/CONVXX
      pressure = -3.0d0*B0*(fac-1.0d0)/fac2*exp(-1.5d0*(B0PRIM-1.0d0)*(fac-1.0d0))
   else if (eos == 4) then
      fac = (VOL0/VOL)**(1.0d0/3.0d0)
      r = log(fac)
      ETOT = E0 + 4.5d0*VOL0*B0*r**2*(1.0d0-r*(B0PRIM-2.0d0))/CONVXX
!     ETOT = E0 + 4.5d0*VOL0*B0*r**2*(1.0d0-r*(B0PRIM-2.0d0))*ULA**3/CONVXX
      pressure = B0*r*(1.0d0+(B0PRIM-2.0d0)/2.0d0*r)
   else if (eos == 5) then
      write(6,'(a,i5)')'Unrecognized equation of state option:',eos
      stop 'Error'
   endif
   bmv = B0*(VOL0/VOL)**B0PRIM
!
   END SUBROUTINE eosfun
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!     function FUN(x_loc) result(f0)
!        integer, parameter :: IntKind = kind(1)
!        integer, parameter :: RealKind = kind(1.0d0)
!        real (kind=RealKind), intent(in) :: x_loc(4)
!        real (kind=RealKind) :: f0
!     end function FUN
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
      IF (N > 0) THEN
         GOTO 2
      ELSE
         GOTO 1
      ENDIF
    1 IER  =1000
      GOTO 109
!   2 IF (LIM) 3,3,4
    2 IF (LIM > 0) THEN
         GOTO 4
      ELSE
         GOTO 3
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
   10 IF (IT-1 > 0) THEN
         GOTO 14
      ELSE
         GOTO 11
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
      IF (IT-2 > 0) THEN
         GOTO 22
      ELSE
         GOTO 23
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
   35 IF ((DC-DA)/(DB-DA) > 0.0D0) THEN
         GOTO 50
      ELSE
         GOTO 36
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
   38 IF (NS-15 > 0) THEN
         GOTO 39
      ELSE
         GOTO 43
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
   47 IF (NS-1 > 0) THEN
         GOTO 49
      ELSE
         GOTO 48
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
      IF (IS > 0) THEN
         GOTO 80
      ELSE
         GOTO 75
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
      IF (IT-N > 0) THEN
         GOTO 78
      ELSE
         GOTO 79
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
   82 IF (IS-1 > 0) THEN
         GOTO 94
      ELSE
         GOTO 83
      ENDIF
!  83 IF (I-1) 84,84,85
   83 IF (I-1 > 0) THEN
         GOTO 85
      ELSE
         GOTO 84
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
      IF (HQ == 0.0D0) THEN
         GOTO 88
      ELSE
         GOTO 90
      ENDIF
!  88 IF (MF-N1) 108,108,89
   88 IF (MF-N1 > 0) THEN
         GOTO 89
      ELSE
         GOTO 108
      ENDIF
   89 IER  =200
      GOTO 108
   90 DQ   =DSQRT(HQ)
      HQ   =DQ
!     IF (HQ-1.D0) 92,92,91
      IF (HQ-1.D0 > 0.0D0) THEN
         GOTO 91
      ELSE
         GOTO 92
      ENDIF
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
      IF (ABS(FS)-1.D0 > 0.0D0) THEN
         GOTO 95
      ELSE
         GOTO 96
      ENDIF
   95 TOL  =EPS*ABS(FS)
!  96 IF (FL-FS-TOL) 100,100,97
   96 IF (FL-FS-TOL > 0.0D0) THEN
         GOTO 97
      ELSE
         GOTO 100
      ENDIF
!  97 IF (IT-LIM) 99,99,98
   97 IF (IT-LIM > 0) THEN
         GOTO 98
      ELSE
         GOTO 99
      ENDIF
   98 IER  =10
      GOTO 108
   99 MF   =0
      GOTO 6
! 100 IF (MF-N1) 102,102,101
  100 IF (MF-N1 > 0) THEN
         GOTO 101
      ELSE
         GOTO 102
      ENDIF
  101 IER  =200
      GOTO 108
  102 MF   =MF+1
      DQ   =0.D0
      DO 105 K=1,N
         J    =N+K
!        IF (ABS(AUX(K))-1.D0) 103,103,104
         IF (ABS(AUX(K))-1.D0 > 0.0D0) THEN
            GOTO 104
         ELSE
            GOTO 103
         ENDIF
  103    DQ   =DQ+ABS(AUX(J))
         GOTO 105
  104    DQ   =DQ+ABS(AUX(J)/AUX(K))
  105    CONTINUE
!     IF (DQ-ETA) 108,108,106
      IF (DQ-ETA > 0.0D0) THEN
         GOTO 106
      ELSE
         GOTO 108
      ENDIF
! 106 IF (MF-N) 6,6,107
  106 IF (MF-N > 0) THEN
         GOTO 107
      ELSE
         GOTO 6
      ENDIF
  107 IER  =1
  108 Y    =FB
!     IF (IER) 111,111,109
      IF (IER > 0) THEN
         GOTO 109
      ELSE
         GOTO 111
      ENDIF
! 109 IF (ISW+12345) 110,111,110
  109 IF (ISW+12345 == 0) THEN
         GOTO 111
      ELSE
         GOTO 110
      ENDIF
! 110 CALL WIER(IER,20212)
  110 CONTINUE
  111 RETURN
  END SUBROUTINE DFMND
!
end program eosfit
