module ContourModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicParamDefinitionsModule, only : HalfCircle, RectBox, HorizLine, VertLine, &
                                            ButterFly, MatsubaraPoles
   use PublicParamDefinitionsModule, only : EqualInterval, GaussianPoints, &
                                            LogInterval, NicholsonPoints
!
public :: initContour,      &
          endContour,       &
          setupContour,     &
          setupOffsetE,     &
          getNumEs,         &
          getEPoint,        &
          getEWeight,       &
          getEPointOffset,  &
          getEWeightOffset, &
          getNumContours,   &
          getContourKPoint, &
          isHorizontalContour, &
          isNotSearchingEf, &
          isMatsubaraContour, &
          printContour, &
          eGridType_Lloyd_Ef_search
!
   interface initContour
      module procedure initContour0, initContour1
   end interface
!
   interface getNumEs
      module procedure getNes, getNesatik
   end interface
!
   interface getEPoint
      module procedure getEall, getEatie, getEatikie
   end interface
!
   interface getEWeight
      module procedure getEWall, getEWatie
   end interface
!
   interface getEPointOffset
      module procedure getEOffsetall, getEOffsetatie
   end interface
!
   integer (kind=IntKind), parameter :: eGridType_Lloyd_Ef_search=3
!
private
   character (len=50) :: stop_routine
!
   logical :: NoLastE = .false.
!
   integer (kind=IntKind), parameter :: DefaultNumEsOfSmallContour = 10
   integer (kind=IntKind) :: NumEs_sc1 = 0, NumEs_sc2 = 0
   integer (kind=IntKind) :: NumEs
   integer (kind=IntKind) :: NumContours = 0
   integer (kind=IntKind) :: NumKGroups = 0
   integer (kind=IntKind) :: ContourType
   integer (kind=IntKind) :: eGridType
   integer (kind=IntKind) :: print_level
!
   real (kind=RealKind) :: Temperature
   real (kind=RealKind) :: ErBottom
   real (kind=RealKind) :: ErTop
   real (kind=RealKind) :: EiBottom
   real (kind=RealKind) :: EiTop
! 
!  ========================
!  Contour type parameters:
!  ========================
   integer (kind=IntKind), parameter :: ReadEmesh=-1
!
   complex (kind=CmplxKind), allocatable, target :: EnergyMesh(:)
   complex (kind=CmplxKind), allocatable, target :: EnergyWght(:)
!
!  ===================================================================
!  EnergyKMesh is only used for band structure calculations
!  ===================================================================
   type EvsKStruct
      integer (kind=IntKind) :: nume
      real (kind=RealKind) :: kmesh(3)
      complex (kind=CmplxKind), pointer :: emesh(:)
   end type EvsKStruct
!
   type (EvsKStruct), allocatable :: EnergyKMesh(:)
!
!  ==========================================
!  E-offset fitting energy points parameters:
!  ==========================================
   integer (kind=IntKind) :: NumOffsetE = 0
   integer (kind=IntKind) :: Offset = 0
!
   complex (kind=CmplxKind), allocatable, target :: EnergyOffset(:)
!
   logical :: Initialized = .false.
   logical :: SearchingFermiEnergy = .true.
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initContour0(filename, istop, iprint)
!  ===================================================================
   use MathParamModule, only : TEN2m8, ZERO
   use MPPModule, only : MyPE, bcastMessage
   use ScfDataModule, only : OffsetE, NumExtraEs
   implicit none
!
   character (len=*), intent(in) :: istop
   character (len=*), intent(in) :: filename
!
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind) :: ie, ik, n, ios, EfSearch, msg(2)
!
   integer (kind=IntKind), parameter :: ReadPE = 0
!
   real (kind=RealKind), parameter :: etol = TEN2m8
!
   ContourType = ReadEmesh
   stop_routine = istop
   print_level = iprint
   Temperature = ZERO
!
   if (MyPE == ReadPE) then
      open(unit=12,file=filename,form='formatted',status='old',iostat=ios)
      if (ios /= 0) then
         call ErrorHandler('initContour0','invalid file name',filename)
      endif
!
      read(12,*)msg(1:2)
      if (msg(1) <= 0) then
         msg(1) = 0
      endif
      if (msg(2) <= 0) then
         msg(2) = 0
      endif
!     ----------------------------------------------------------------
      call bcastMessage(msg(1:2),2,ReadPE)
!     ----------------------------------------------------------------
      NumKGroups = msg(1)
      EfSearch = msg(2)
!
      if (NumKGroups == 0) then
         NumContours = 1
         read(12,*)NumEs
!        -------------------------------------------------------------
         call bcastMessage(NumEs,ReadPE)
!        -------------------------------------------------------------
         allocate( EnergyMesh(NumEs), EnergyWght(NumEs) )
         do ie=1,NumEs
            read(12,*)EnergyMesh(ie),EnergyWght(ie)
            if (abs(EnergyMesh(ie)) .lt. etol) then
               EnergyMesh(ie)=EnergyMesh(ie)+0.001
            endif
         enddo
!        -------------------------------------------------------------
         call bcastMessage(EnergyMesh,NumEs,ReadPE)
         call bcastMessage(EnergyWght,NumEs,ReadPE)
!        -------------------------------------------------------------
      else
         NumContours = NumKGroups
         allocate( EnergyKMesh(NumContours) )
         NumEs=0
         do ik=1,NumContours
            read(12,*)EnergyKMesh(ik)%kmesh(1),EnergyKMesh(ik)%kmesh(2), &
                      EnergyKMesh(ik)%kmesh(3)
            read(12,*)n
!           ----------------------------------------------------------
            call bcastMessage(EnergyKMesh(ik)%kmesh(1:3),3,ReadPE)
            call bcastMessage(n,ReadPE)
!           ----------------------------------------------------------
            NumEs=NumEs+n
            EnergyKMesh(ik)%nume = n
            allocate( EnergyKMesh(ik)%emesh(n) )
            do ie=1,n
!              -------------------------------------------------------
               read(12,*)EnergyKMesh(ik)%emesh(ie)
!              -------------------------------------------------------
               if (abs(EnergyKMesh(ik)%emesh(ie)) .lt. etol) then
                  EnergyKMesh(ik)%emesh(ie)=EnergyKMesh(ik)%emesh(ie)+0.001
               endif
            enddo
!           ----------------------------------------------------------
            call bcastMessage(EnergyKMesh(ik)%emesh,n,ReadPE)
!           ----------------------------------------------------------
         enddo
      endif
      close(12)
   else
!     ----------------------------------------------------------------
      call bcastMessage(msg(1:2),2,ReadPE)
!     ----------------------------------------------------------------
      NumKGroups = msg(1)
      EfSearch = msg(2)
      if (NumKGroups == 0) then
         NumContours = 1
!        -------------------------------------------------------------
         call bcastMessage(NumEs,ReadPE)
!        -------------------------------------------------------------
         allocate( EnergyMesh(NumEs), EnergyWght(NumEs) )
!        -------------------------------------------------------------
         call bcastMessage(EnergyMesh,NumEs,ReadPE)
         call bcastMessage(EnergyWght,NumEs,ReadPE)
!        -------------------------------------------------------------
      else
         NumContours = NumKGroups
         allocate( EnergyKMesh(NumContours) )
         NumEs=0
         do ik=1,NumContours
!           ----------------------------------------------------------
            call bcastMessage(EnergyKMesh(ik)%kmesh(1:3),3,ReadPE)
            call bcastMessage(n,ReadPE)
!           ----------------------------------------------------------
            NumEs=NumEs+n
            EnergyKMesh(ik)%nume = n
            allocate( EnergyKMesh(ik)%emesh(n) )
!           ----------------------------------------------------------
            call bcastMessage(EnergyKMesh(ik)%emesh,n,ReadPE)
!           ----------------------------------------------------------
         enddo
      endif
   endif
!
   if (NumEs == 1 .or. EfSearch == 0) then
      SearchingFermiEnergy = .false.
   endif
!
   if (OffsetE > NumEs) then
      call ErrorHandler('initContour1','OffsetE > NumEs',OffsetE, NumEs)
   endif
!
   NoLastE = .true.
!  Offset = OffsetE
!  NumOffsetE = NumExtraEs
!  allocate( EnergyOffset(NumOffsetE) )
!
   Initialized = .true.
!
   end subroutine initContour0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initContour1(ctype,eg,ne,Temp,istop,iprint,nle)
!  ===================================================================
   use MathParamModule, only : ZERO, CZERO, ONE, TEN2m6
   use MPPModule, only : MyPE
   use ScfDataModule, only : OffsetE, NumExtraEs
!   use ScfDataModule, only : isLloyd,getLloydMode
   use MatsubaraModule, only : initMatsubara, getNumMatsubaraPoles
!
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: ctype
   integer (kind=IntKind), intent(in) :: iprint
   integer (kind=IntKind), intent(in) :: eg
   integer (kind=IntKind), intent(in) :: ne
!
   real (kind=RealKind), intent(in) :: Temp
!
   logical, intent(in), optional :: nle
!
   if (ne < 1) then
      call ErrorHandler('initContour1','invalid no. of energies',ne)
   else if (Temp < ZERO) then
      call ErrorHandler('initContour1','Temperature < 0',Temp)
   endif
!
   if (ctype == ReadEmesh) then
      call ErrorHandler('initContour1','Invalid contour type')
   endif
!
   Temperature = Temp
!
   if (Temperature < TEN2m6) then
      ContourType = ctype
   else
      ContourType = MatsubaraPoles
   endif
!
   if (ContourType == MatsubaraPoles) then
!     ================================================================
!     If temperature is too low, e.g. less than 300.0 Kevin, while contour
!     type is determined to be Matsubara poles, the XG Zhang's Gaussian
!     points scheme is used so to avoid generating too many poles in
!     using Nicholson's scheme.
!     ================================================================
      if (Temperature > 299.0d0 .and. eg == NicholsonPoints) then
         eGridType = NicholsonPoints
!        -------------------------------------------------------------
         call initMatsubara(eGridType,Temperature,eWidth=1.5d0,iprint=iprint)
!        -------------------------------------------------------------
      else
         eGridType = GaussianPoints
!        -------------------------------------------------------------
         call initMatsubara(eGridType,Temperature,NumPs=ne,iprint=iprint)
!        -------------------------------------------------------------
      endif
      NumEs = getNumMatsubaraPoles()
   else if (eg == NicholsonPoints) then
      if (MyPE == 0) then
!        -------------------------------------------------------------
         call WarningHandler('initContour1',                          &
                             'Energy grid type is set to Gaussian')
!        -------------------------------------------------------------
      endif
      eGridType = GaussianPoints
      NumEs = ne
   else
      eGridType = eg
      NumEs = ne
   endif
!
   NumContours = 0
   stop_routine = istop
   print_level = iprint
!
   if (present(nle)) then
      NoLastE = nle ! If NoLastE = true, there is no extra energy point added to the contour
                    ! Due to the legacy code, by default, it is set to be false.
   endif
!
   if (eGridType == GaussianPoints .or. eGridType == NicholsonPoints) then
      if (.not.NoLastE) then
         NumEs = NumEs + 1
      endif
!      if (isLloyd().and.(getLloydMode().eq.1)) then
!         NumEs = NumEs + 1
!      end if
   endif
!
   if (OffsetE > NumEs) then
      call ErrorHandler('initContour1','OffsetE > NumEs',OffsetE, NumEs)
   endif
!
   if (NumEs == 1) then
      call WarningHandler('initContour','NumEs = 1')
   endif
!
   Offset = OffsetE
   NumOffsetE = NumExtraEs
!  allocate( EnergyOffset(NumOffsetE) )
!
   if (ContourType == ButterFly) then
      NumEs_sc1 = max(DefaultNumEsOfSmallContour,NumExtraEs)
      NumEs_sc2 = NumEs
      NumEs = NumEs_sc1 + NumEs_sc2
   else
      NumEs_sc1 = 0
      NumEs_sc2 = NumEs
   endif
!
   allocate( EnergyMesh(NumEs), EnergyWght(NumEs) )
   EnergyMesh = CZERO
   EnergyWght = CZERO
!
   Initialized = .true.
!
   end subroutine initContour1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupContour(erb,ert,eib,eit,core_top)
!  ===================================================================
   use MathParamModule, only : ZERO, CONE, CZERO
   use ScfDataModule, only : NumExtraEs
   use MatsuBaraModule, only : calMatsubaraPoles
   implicit none
!
   real (kind=RealKind), intent(in) :: erb
   real (kind=RealKind), intent(in) :: ert
   real (kind=RealKind), intent(in) :: eib
   real (kind=RealKind), intent(in) :: eit
   real (kind=RealKind), intent(in), optional :: core_top
   real (kind=RealKind) :: etopcor
!
   if (.not.Initialized) then
      call ErrorHandler('setupContour','ContourModule is not initialized')
   endif
!
   ErBottom = erb
   ErTop = ert
   EiBottom = eib
   EiTop = eit
!
   if (present(core_top)) then
      etopcor = core_top
   else
      etopcor = -2.5d0
   endif
!
   if (ContourType == ButterFly) then
      NumEs_sc1 = max(DefaultNumEsOfSmallContour,NumExtraEs)
      NumEs_sc2 = NumEs - NumEs_sc1  ! number of E's in the 2nd semi-circle contour
   else
      NumEs_sc1 = 0
      NumEs_sc2 = NumEs
   endif
!
   NumContours = 1
!
   if (ContourType == ReadEmesh) then
      call WarningHandler('setupContour','Energy contour is read in')
      return
   else if (ErTop < ErBottom) then
      call ErrorHandler('initContour1','ErTop < ErBottom',ErTop,ErBottom)
   else if (EiTop < EiBottom) then
      call ErrorHandler('initContour1','EiTop < EiBottom',EiTop,EiBottom)
   else if (EiBottom < ZERO) then
      call ErrorHandler('initContour1','EiBottom < 0',EiBottom)
   else if (NumEs == 1) then
      EnergyMesh(1) = cmplx(ErBottom,EiBottom,CmplxKind)
      EnergyWght(1) = CONE
      return
   endif
!
   if (ContourType == MatsubaraPoles) then
!     ----------------------------------------------------------------
      call calMatsubaraPoles(ErBottom,ErTop,EnergyMesh,EnergyWght,    &
                             core_top=etopcor)
!     ----------------------------------------------------------------
      if (.not.NoLastE) then
         EnergyMesh(NumEs) = cmplx(ErTop,EiBottom,CmplxKind)
         EnergyWght(NumEs) = CZERO
      endif
   else
      if (eGridType == EqualInterval) then
!        -------------------------------------------------------------
         call setupEqlGrid()
!        -------------------------------------------------------------
      else if (eGridType == GaussianPoints) then
!        -------------------------------------------------------------
         call setupGauGrid()
!        -------------------------------------------------------------
      else if (eGridType == LogInterval) then
!        -------------------------------------------------------------
         call setupLogGrid()
!        -------------------------------------------------------------
      else if (eGridType == eGridType_Lloyd_Ef_search) then
!        -------------------------------------------------------------
         call setupLloyd_Grid_Ef_search()
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call ErrorHandler('setupContour','Invalid eGridType',eGridType)
!        -------------------------------------------------------------
      endif
   endif
!
   end subroutine setupContour
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupEqlGrid()
!  ===================================================================
   use MathParamModule, only : HALF, PI, CZERO, CONE, SQRTm1, ten2m8
!   use ScfDataModule, only : isLloyd,getLloydMode
   implicit none
!
   integer (kind=IntKind) :: ie, n1, n2, n3
!
   real (kind=RealKind) :: de, el, dp, er, phi
   real (kind=RealKind), parameter :: etol = TEN2m8
!
   complex (kind=CmplxKind) :: es
!
   integer (kind=IntKind) :: iLloyd
!
   iLloyd=0
!   if (isLloyd().and.(getLloydMode().eq.1)) then
!      iLloyd=1
!   end if
!
   if (ContourType == HalfCircle) then
      dp = PI/real(NumEs-1-iLloyd,kind=RealKind)
      er = (ErTop-ErBottom)*HALF
      es = cmplx(HALF*(ErBottom+ErTop),EiBottom,CmplxKind)
      do ie=1,NumEs-iLloyd
         phi=PI-dp*(ie-1)
         EnergyMesh(ie)=er*exp(SQRTm1*phi) + es
      enddo
      EiTop = EiBottom + er
   else if (ContourType == RectBox) then
      if (NumEs < 4) then
!        -------------------------------------------------------------
         call ErrorHandler('setupEGrid','NumEs < 4',NumEs)
!        -------------------------------------------------------------
      endif
      el = 2*(EiTop-EiBottom)+ErTop-ErBottom
      n1 = int((EiTop-EiBottom)/el+HALF)*NumEs-1
      if (n1 <= 1) then
         n1 = 1
         EnergyMesh(n1)=cmplx(ErBottom,EiBottom,CmplxKind)
      else
         de = (EiTop-EiBottom)/real(n1,kind=RealKind)
         do ie=1,n1
            EnergyMesh(ie)=cmplx(ErBottom,EiBottom+de*(ie-1),CmplxKind)
         enddo
      endif
      n2 = int((ErTop-ErBottom)/el+HALF)*NumEs-1
      if (n2 <= 1) then
         n2 = 1
         EnergyMesh(n1+n2)=cmplx(ErBottom,EiTop,CmplxKind)
      else
         de = (ErTop-ErBottom)/real(n2,kind=RealKind)
         do ie=1,n2
            EnergyMesh(n1+ie)=cmplx(ErBottom+de*(ie-1),EiTop,CmplxKind)
         enddo
      endif
      n3 = NumEs - n1 - n2
      if (n3 <= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('setupEGrid','n3 <= 0',n3)
!        -------------------------------------------------------------
      else if (n3 == 1) then
         EnergyMesh(NumEs)=cmplx(ErTop,EiTop,CmplxKind)
         call WarningHandler('setupEGrid','(ErTop,EiBottom) is not included')
      else
         de = -(EiTop-EiBottom)/real(n3-1,kind=RealKind)
         do ie=1,n3
            EnergyMesh(n1+n2+ie)=cmplx(ErTop,EiTop+de*(ie-1),CmplxKind)
         enddo
      endif
   else if (ContourType == VertLine) then
      ErTop = ErBottom
      de = -(EiTop-EiBottom)/real(NumEs-1,kind=RealKind)
      do ie=1,NumEs
         EnergyMesh(ie)=cmplx(ErBottom,EiTop+de*(ie-1),CmplxKind)
      enddo
   else if (ContourType == HorizLine) then
      EiTop = EiBottom
      de =  (ErTop-ErBottom)/real(NumEs-1,kind=RealKind)
      do ie=1,NumEs
         EnergyMesh(ie)=cmplx(ErBottom+de*(ie-1),EiBottom,CmplxKind)
         if (abs(EnergyMesh(ie)) < etol) then
            EnergyMesh(ie) = EnergyMesh(ie) + 0.001d0
         endif
      enddo
   else if (ContourType == ButterFly) then
      dp = PI/real(NumEs_sc1-1,kind=RealKind)
      er = abs(ErBottom+EiBottom)*HALF
      es = HALF*(ErBottom-EiBottom) ! center location of the 1st circle
      do ie = 1, NumEs_sc1
         phi=PI-dp*(ie-1)
         EnergyMesh(ie)=er*exp(SQRTm1*phi) + es
      enddo
!
      dp = PI/real(NumEs_sc2-1-iLloyd,kind=RealKind)
      er = abs(ErTop+EiBottom)*HALF
      es = HALF*(ErTop-EiBottom) ! center location of the 2nd circle
      do ie=1,NumEs_sc2-iLloyd
         phi=PI-dp*(ie-1)
         EnergyMesh(NumEs_sc1+ie)=er*exp(SQRTm1*phi) + es
      enddo
      EiTop = er
   else
      call ErrorHandler('setupEqlGrid','Invalid Contour type',ContourType)
   endif ! countour type
   if (ContourType == ButterFly) then
      EnergyWght(1)=( EnergyMesh(2)-EnergyMesh(1) )*HALF
      do ie=2,NumEs_sc1-1
         EnergyWght(ie)=( EnergyMesh(ie+1)-EnergyMesh(ie-1) )*HALF
      enddo
      EnergyWght(NumEs_sc1)=CZERO
      EnergyWght(NumEs_sc1+1)=( EnergyMesh(NumEs_sc1+2)-EnergyMesh(NumEs_sc1+1) )*HALF
      do ie=NumEs_sc1+2,NumEs-1-iLloyd
         EnergyWght(ie)=( EnergyMesh(ie+1)-EnergyMesh(ie-1) )*HALF
      enddo
      do ie=NumEs-iLloyd,NumEs
         EnergyWght(ie)=CZERO
      enddo
   else
      EnergyWght(1)=( EnergyMesh(2)+EnergyMesh(1) )*HALF-ErBottom
      EnergyWght(NumEs)=CZERO
      do ie=2,NumEs-1-iLloyd
         EnergyWght(ie)=( EnergyMesh(ie+1)-EnergyMesh(ie-1) )*HALF
      enddo
   endif
!
   if (iLloyd.eq.1) then
      EnergyMesh(NumEs)=EnergyMesh(NumEs-1)
      EnergyWght(NumEs)=EnergyWght(NumEs-1)
      EnergyMesh(NumEs-1)=cmplx(ErTop,ten2m8,CmplxKind)
      EnergyWght(NumEs-1)=CZERO
   endif
!
   end subroutine setupEqlGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupLogGrid()
!  ===================================================================
   use MathParamModule, only : HALF, PI, CZERO, CONE, SQRTm1, TEN2m8
   implicit none
!
   integer (kind=IntKind) :: ie, n1, n2, n3
!
   real (kind=RealKind) :: de, el, dp, er, phi
   real (kind=RealKind), parameter :: etol = TEN2m8
!
   complex (kind=CmplxKind) :: es
!
   if (ContourType == HalfCircle) then
      dp = log(PI)/real(NumEs-1, kind=RealKind)
      er = (ErTop-ErBottom)*HALF
      es = cmplx(HALF*(ErBottom+ErTop),EiBottom,CmplxKind)
      do ie=1,NumEs
         phi=PI-exp(dp*(ie-1))
         EnergyMesh(ie)=er*exp(SQRTm1*phi) + es
      enddo
      EiTop = EiBottom + er
   else if (ContourType == RectBox) then
      if (NumEs < 4) then
!        -------------------------------------------------------------
         call ErrorHandler('setupEGrid','NumEs < 4',NumEs)
!        -------------------------------------------------------------
      endif
      el = 2*(EiTop-EiBottom)+ErTop-ErBottom
      n1 = int((EiTop-EiBottom)/el+HALF)*NumEs-1
      if (n1 <= 1) then
         n1 = 1
         EnergyMesh(n1)=cmplx(ErBottom,EiBottom,CmplxKind)
      else
         de = (EiTop-EiBottom)/real(n1,kind=RealKind)
         do ie=1,n1
            EnergyMesh(ie)=cmplx(ErBottom,EiBottom+de*(ie-1),CmplxKind)
         enddo
      endif
      n2 = int((ErTop-ErBottom)/el+HALF)*NumEs-1
      if (n2 <= 1) then
         n2 = 1
         EnergyMesh(n1+n2)=cmplx(ErBottom,EiTop,CmplxKind)
      else
         de = (ErTop-ErBottom)/real(n2,kind=RealKind)
         do ie=1,n2
            EnergyMesh(n1+ie)=cmplx(ErBottom+de*(ie-1),EiTop,CmplxKind)
         enddo
      endif
      n3 = NumEs - n1 - n2
      if (n3 <= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('setupEGrid','n3 <= 0',n3)
!        -------------------------------------------------------------
      else if (n3 == 1) then
         EnergyMesh(NumEs)=cmplx(ErTop,EiTop,CmplxKind)
         call WarningHandler('setupEGrid','(ErTop,EiBottom) is not included')
      else
         de = -log(EiTop/EiBottom)/real(n3-1,kind=RealKind)
         do ie=1,n3
            EnergyMesh(n1+n2+ie)=cmplx(ErTop,EiTop*exp(de*(ie-1)),CmplxKind)
         enddo
      endif
   else if (ContourType == VertLine) then
      ErTop = ErBottom
      de = -log(EiTop/EiBottom)/real(NumEs-1,kind=RealKind)
      do ie=1,NumEs
         EnergyMesh(ie)=cmplx(ErBottom,EiTop*exp(de*(ie-1)),CmplxKind)
      enddo
   else if (ContourType == HorizLine) then
      EiTop = EiBottom
      de =  log(ErTop/ErBottom)/real(NumEs-1,kind=RealKind)
      do ie=1,NumEs
         EnergyMesh(ie)=cmplx(ErBottom*exp(de*(ie-1)),EiBottom,CmplxKind)
         if (abs(EnergyMesh(ie)) < etol) then
            EnergyMesh(ie) = EnergyMesh(ie) + 0.001d0
         endif
      enddo
   else if (ContourType == ButterFly) then
      dp = log(PI)/real(NumEs_sc1-1, kind=RealKind)
      er = abs(ErBottom+EiBottom)*HALF
      es = HALF*(ErBottom-EiBottom)
      do ie=1,NumEs_sc1
         phi=PI-exp(dp*(ie-1))
         EnergyMesh(ie)=er*exp(SQRTm1*phi) + es
      enddo
!
      dp = log(PI)/real(NumEs_sc2-1, kind=RealKind)
      er = abs(ErTop+EiBottom)*HALF
      es = HALF*(ErTop-EiBottom)
      do ie=1,NumEs_sc2
         phi=PI-exp(dp*(ie-1))
         EnergyMesh(NumEs_sc1+ie)=er*exp(SQRTm1*phi) + es
      enddo
      EiTop = er
   else
      call ErrorHandler('setupEqlGrid','Invalid Contour type',ContourType)
   endif
   if (ContourType == ButterFly) then
      EnergyWght(1)=( EnergyMesh(2)-EnergyMesh(1) )*HALF
      do ie=2,NumEs_sc1-1
         EnergyWght(ie)=( EnergyMesh(ie+1)-EnergyMesh(ie-1) )*HALF
      enddo
      EnergyWght(NumEs_sc1)=CZERO
      EnergyWght(NumEs_sc1+1)=( EnergyMesh(NumEs_sc1+2)-EnergyMesh(NumEs_sc1+1) )*HALF
      do ie=NumEs_sc1+2,NumEs-1
         EnergyWght(ie)=( EnergyMesh(ie+1)-EnergyMesh(ie-1) )*HALF
      enddo
      EnergyWght(NumEs)=CZERO
   else
      EnergyWght(1)=( EnergyMesh(2)+EnergyMesh(1) )*HALF-ErBottom
      EnergyWght(NumEs)=CZERO
      do ie=2,NumEs-1
         EnergyWght(ie)=( EnergyMesh(ie+1)-EnergyMesh(ie-1) )*HALF
      enddo
   endif
!
   end subroutine setupLogGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupGauGrid()
!  ===================================================================
   use MathParamModule, only : ZERO, HALF, ONE, PI, CZERO, CONE, SQRTm1
   use MathParamModule, only : ten2m8, TEN2m6
!   use ScfDataModule, only : isLloyd,getLloydMode
   implicit none
!
   integer (kind=IntKind) :: ie, ne, n1, n2, n3
!
   real (kind=RealKind) :: el, dp, er
   real (kind=RealKind), allocatable :: xg(:), wg(:)
   real (kind=RealKind), parameter :: etol = TEN2m8
!
   complex (kind=CmplxKind) :: ec, es, check_wght
!
   integer (kind=IntKind) :: iLloyd
!
   iLloyd=0
!   if (isLloyd().and.(getLloydMode().eq.1)) then
!      iLloyd=1
!   end if
!
   if (NoLastE) then
      ne = NumEs - iLloyd
   else
      ne = NumEs - 1 - iLloyd
   endif
!
   allocate(xg(ne), wg(ne))
!
   if (ContourType == HalfCircle) then
!     ----------------------------------------------------------------
      call gauleg(-ONE, ONE, xg, wg, ne)
!     ----------------------------------------------------------------
      er = (ErTop-ErBottom)*HALF
      ec = PI*HALF*SQRTm1
      check_wght = CZERO
      do ie=1,ne
         es = er*exp(ec*(ONE-xg(ie)))
         EnergyMesh(ie)=HALF*(ErTop+ErBottom) + es
         EnergyWght(ie)=-ec*es*wg(ie)
         check_wght = check_wght + EnergyWght(ie)
      enddo
      EiTop = EiBottom + er
      if (abs(check_wght-(ErTop-ErBottom)) > TEN2m6) then
         call ErrorHandler('setupGauGrid','Weight check is failed',check_wght)
      endif
!      if ( Offset>0 ) then
!         call setupOffsetE()
!      endif
   else if (ContourType == RectBox) then
      if (ne < 3) then
!        -------------------------------------------------------------
         call ErrorHandler('setupGauGrid','NumEs <= 3',NumEs)
!        -------------------------------------------------------------
      endif
      dp = 2*(EiTop-EiBottom)+ErTop-ErBottom
      n1 = int((EiTop-EiBottom)/dp+HALF)*ne
!     ----------------------------------------------------------------
      call gauleg(-ONE, ONE, xg, wg, n1)
!     ----------------------------------------------------------------
      if (n1 <= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('setupGauGrid','n1 <= 0',n1)
!        -------------------------------------------------------------
      else if (n1 == 1) then
         EnergyMesh(n1)=cmplx(ErBottom,HALF*(EiBottom+EiTop),CmplxKind)
         EnergyWght(n1)=EiTop-EiBottom
      else
!        -------------------------------------------------------------
         call gauleg(-ONE, ONE, xg, wg, n1)
!        -------------------------------------------------------------
         el = HALF*(EiBottom+EiTop)
         er = HALF*(EiTop-EiBottom)
         do ie=1,n1
            EnergyMesh(ie)=cmplx(ErBottom,el+er*xg(ie),CmplxKind)
            EnergyWght(ie)=wg(ie)*er
         enddo
      endif
      n2 = int((ErTop-ErBottom)/dp+HALF)*ne
      if (n2 <= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('setupGauGrid','n2 <= 0',n2)
!        -------------------------------------------------------------
      else if (n2 == 1) then
         EnergyMesh(n1+n2)=cmplx(HALF*(ErBottom+ErTop),EiTop,CmplxKind)
         EnergyWght(n1+n2)=ErTop-ErBottom
      else
!        -------------------------------------------------------------
         call gauleg(-ONE, ONE, xg, wg, n2)
!        -------------------------------------------------------------
         el = HALF*(ErBottom+ErTop)
         er = HALF*(ErTop-ErBottom)
         do ie=1,n2
            EnergyMesh(n1+ie)=cmplx(el+er*xg(ie),EiTop,CmplxKind)
            EnergyWght(n1+ie)=wg(ie)*er
         enddo
      endif
      n3 = ne - n1 - n2
      if (n3 <= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('setupGauGrid','n3 <= 0',n3)
!        -------------------------------------------------------------
      else if (n3 == 1) then
         EnergyMesh(ne)=cmplx(ErTop,HALF*(EiTop+EiBottom),CmplxKind)
         EnergyWght(ne)=EiBottom-EiTop
      else
!        -------------------------------------------------------------
         call gauleg(-ONE, ONE, xg, wg, n3)
!        -------------------------------------------------------------
         el = HALF*(EiBottom+EiTop)
         er = HALF*(EiBottom-EiTop)
         do ie=1,n3
            EnergyMesh(n1+n2+ie)=cmplx(ErTop,el+er*xg(ie),CmplxKind)
            EnergyWght(n1+n2+ie)=wg(ie)*er
         enddo
      endif
   else if (ContourType == VertLine) then
!     ----------------------------------------------------------------
      call gauleg(-ONE, ONE, xg, wg, ne)
!     ----------------------------------------------------------------
      ErTop = ErBottom
      el = HALF*(EiBottom+EiTop)
      er = HALF*(EiBottom-EiTop)
      do ie=1,ne
         EnergyMesh(ie)=cmplx(ErBottom,el+er*xg(ie),CmplxKind)
         EnergyWght(ie)=wg(ie)*er
      enddo
   else if (ContourType == HorizLine) then
!     ----------------------------------------------------------------
      call gauleg(-ONE, ONE, xg, wg, ne)
!     ----------------------------------------------------------------
      EiTop = EiBottom
      el = HALF*(ErBottom+ErTop)
      er = HALF*(ErTop-ErBottom)
      do ie=1,ne
         EnergyMesh(ie)=cmplx(el+er*xg(ie),EiBottom,CmplxKind)
         EnergyWght(ie)=wg(ie)*er
         if (abs(EnergyMesh(ie)) < etol) then
            EnergyMesh(ie) = EnergyMesh(ie) + 0.001d0
         endif
      enddo
   else if (ContourType == ButterFly) then
      ne = NumEs_sc1
!     ----------------------------------------------------------------
      call gauleg(-ONE, ONE, xg, wg, ne)
!     ----------------------------------------------------------------
      er = abs(ErBottom)*HALF
      ec = PI*HALF*SQRTm1
      check_wght = CZERO
      do ie=1,ne
         es = er*exp(ec*(ONE-xg(ie)))
         EnergyMesh(ie)=HALF*ErBottom + es
         EnergyWght(ie)=-ec*es*wg(ie)
         check_wght = check_wght + EnergyWght(ie)
      enddo
      if (abs(check_wght+ErBottom) > TEN2m6) then
         call ErrorHandler('setupGauGrid','Weight check is failed',check_wght)
      endif
!
      if (NoLastE) then
         ne = NumEs_sc2 - iLloyd
      else
         ne = NumEs_sc2 - 1 - iLloyd
      endif
!     ----------------------------------------------------------------
      call gauleg(-ONE, ONE, xg, wg, ne)
!     ----------------------------------------------------------------
      er = ErTop*HALF
      ec = PI*HALF*SQRTm1
      check_wght = CZERO
      do ie=1,ne
         es = er*exp(ec*(ONE-xg(ie)))
         EnergyMesh(NumEs_sc1+ie)=HALF*ErTop + es
         EnergyWght(NumEs_sc1+ie)=-ec*es*wg(ie)
         check_wght = check_wght + EnergyWght(NumEs_sc1+ie)
      enddo
      EiTop = er
      if (abs(check_wght-ErTop) > TEN2m6) then
         call ErrorHandler('setupGauGrid','Weight check is failed',check_wght)
      endif
   else
      call ErrorHandler('setupGauGrid','Invalid Contour type',ContourType)
   endif
!
   if (.not.NoLastE) then
!     ================================================================
!     The last point is reserved for extrapolating E_f
!     ================================================================
      EnergyMesh(NumEs)=cmplx(ErTop,EiBottom,CmplxKind)
      EnergyWght(NumEs)=CZERO
   end if
!  ===================================================================
!  The last or the 2nd last point is reserved for Lloyd Value close to ErTop
!  ===================================================================
   if (iLloyd.eq.1) then
      EnergyMesh(ne+1)=cmplx(ErTop,ten2m8,CmplxKind)
      EnergyWght(ne+1)=CZERO
   end if
!
   deallocate(xg, wg)
!
   end subroutine setupGauGrid
!  ===================================================================
!
!  *******************************************************************
!
!  generates emesh where the distance between point follows a power law
!  mesh is symmetric: [ xmid-dx : xmid+dx ] , end points are included
!  power=1.0 corresponds to a linear mesh
!  ei is the imaginary part (constant) for all energy points
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupLloyd_Grid_Ef_search()
     implicit none

     integer(kind=IntKind) :: np,i
     real(kind=RealKind)   :: xmid,dx,power,ei 
     real(kind=RealKind)   :: x1,t1
    
     np=NumEs
     xmid=ErTop
     dx=0.20d0      !!! to be changed later as input parameter
     ei=0.000010d0 !!! to be changed later as input parameter
     power=2.30d0
     
     do i = 1, np
        x1 = -1.0d0 + 2.0d0*dble(i-1)/dble(np-1)
        t1 = Sign(1.0d0,x1)*Abs(x1)**power
        EnergyMesh(i)= cmplx(xmid + t1*dx,ei,RealKind)
        EnergyWght(i)=1.00d0
     end do ! i

   end subroutine setupLloyd_Grid_Ef_search
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endContour()
!  ===================================================================
   use MatsubaraModule, only : endMatsubara
!
   implicit none
!
   integer (kind=IntKind) :: ik
!
   if (NumKGroups > 0) then
      do ik=1,NumContours
         deallocate(EnergyKMesh(ik)%emesh)
      enddo
      deallocate(EnergyKMesh)
   else
      deallocate(EnergyMesh, EnergyWght)
   endif
   NumContours = 0
   NumKGroups = 0
   NumEs = 0
   if ( allocated(EnergyOffset) ) then
      deallocate(EnergyOffset)
      Offset = -1
   endif
!
   if (ContourType == MatsubaraPoles) then
!     ----------------------------------------------------------------
      call endMatsubara()
!     ----------------------------------------------------------------
   endif
!
   end subroutine endContour
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumContours() result(nc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: nc
!
   if (.not.Initialized) then
      call ErrorHandler('getNumContours','ContourModule is not initialized')
   endif
!
   nc = NumContours
   end function getNumContours
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getContourKPoint(ik) result(kp)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   character (len=16), parameter :: sname='getContourKPoint'
   integer (kind=IntKind), intent(in) :: ik
   real (kind=RealKind) :: kp(3)
!
   if (.not.Initialized) then
      call ErrorHandler('getContourKPoint','ContourModule is not initialized')
   endif
!
   if (NumKGroups < 1) then
      kp(1:3)=ZERO
      call WarningHandler(sname,'Cannot return valid k-point')
   else if (ik>NumContours .or. ik<1) then
      call ErrorHandler(sname,'ik out of range',ik)
   else
      kp(1:3)=EnergyKMesh(ik)%kmesh(1:3)
   endif
   end function getContourKPoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNes(LowContour,HighContour) result(ne)
!  ===================================================================
   implicit none
!
   logical, optional, intent(in) :: LowContour, HighContour
   integer (kind=IntKind) :: ne, n1, n2
!
   if (.not.Initialized) then
      call ErrorHandler('getNumEs','ContourModule is not initialized')
   else if (.not.present(LowContour) .and. .not.present(HighContour)) then
      ne = NumEs
      return
   endif
!
   n1 = NumEs_sc1
   n2 = NumEs_sc2
!
   if (present(LowContour)) then
      if (.not.LowContour) then
         n1 = 0
      endif
   endif
!
   if (present(HighContour)) then
      if (.not.HighContour) then
         n2 = 0
      endif
   endif
!
   ne = n1 + n2
!
   end function getNes
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNesatik(ik) result(ne)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ik
   integer (kind=IntKind) :: ne
!
   if (.not.Initialized) then
      call ErrorHandler('getNumEs','ContourModule is not initialized')
   else if (NumKGroups < 1) then
      call ErrorHandler('getNumEs','Invalid function call')
   endif
!
   if (ik>NumContours .or. ik<1) then
      call ErrorHandler('getNumEs','ik out of range',ik)
   endif
   ne = EnergyKMesh(ik)%nume
   end function getNesatik
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEWall() result(wp)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), pointer :: wp(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getEWeight','ContourModule is not initialized')
   else if (NumContours < 1 .or. NumEs < 1) then
      call ErrorHandler('getEWeight','NumContours or NumEs < 1')
   else if (NumKGroups > 0) then
      call ErrorHandler('getEWeight','Energy weight is not defined')
   endif
   wp => EnergyWght(1:NumEs)
!
   end function getEWall
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEWatie(ie) result(w)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ie
!
   complex (kind=CmplxKind) :: w
!
   if (.not.Initialized) then
      call ErrorHandler('getEWeight','ContourModule is not initialized')
   else if (NumContours < 1) then
      call ErrorHandler('getEWeight','NumContours < 1', NumContours)
   else if (NumKGroups > 0) then
      call ErrorHandler('getEWeight','Energy weight is not defined')
   else if (ie < 1 .or. ie > NumEs) then
      call ErrorHandler('getEWeight','Energy index out of bound',ie)
   endif
   w = EnergyWght(ie)
!
   end function getEWatie
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEall() result(ep)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), pointer :: ep(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getEPoint','ContourModule is not initialized')
   endif
!
   if (NumKGroups > 0) then
      call ErrorHandler('getEPoint','Invalid function argument(s)')
   else if (NumContours < 1 .or. NumEs < 1) then
      call ErrorHandler('getEPoint','NumContours or NumEs < 1')
   endif
   ep => EnergyMesh(1:NumEs)
!
   end function getEall
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEatie(ie) result(e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ie
   complex (kind=CmplxKind) :: e
!
   if (.not.Initialized) then
      call ErrorHandler('getEPoint','ContourModule is not initialized')
   endif
!
   if (NumKGroups > 0) then
      call ErrorHandler('getEPoint','Invalid interface arguments')
   else if (NumContours < 1) then
      call ErrorHandler('getEPoint','NumContours < 1', NumContours)
   else
      if (ie>NumEs .or. ie<1) then
         call ErrorHandler('getEPoint','ie out of range',ie)
      else
         e = EnergyMesh(ie)
      endif
   endif
!
   end function getEatie
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEatikie(ik,ie) result(e)
!  ===================================================================
   implicit none
! 
   integer (kind=IntKind), intent(in) :: ik
   integer (kind=IntKind), intent(in) :: ie
   integer (kind=IntKind) :: n
   complex (kind=CmplxKind) :: e
!
   if (.not.Initialized) then
      call ErrorHandler('getEPoint','ContourModule is not initialized')
   endif
!
   if (NumKGroups > 0) then
      if (ik>NumContours .or. ik<1) then
         call ErrorHandler('getEPoint','ik out of range',ik)
      else
         n = EnergyKMesh(ik)%nume
         if (ie>n .or. ie<1) then
            call ErrorHandler('getEPoint','ie out of range',ie)
         else
            e = EnergyKMesh(ik)%emesh(ie)
         endif
      endif
   else
      call ErrorHandler('getEPoint','Invalid interface arguments')
   endif
   end function getEatikie
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupOffsetE()
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind) :: i
   real (kind=RealKind) :: d_rImE
!
   if ( mod(NumOffsetE,2) == 0 ) then
      d_rImE = 0.8d0*aimag(EnergyMesh(Offset))/ &
               real((NumOffsetE/2),kind=RealKind)
      do i = -NumOffsetE/2,NumOffsetE/2
!           EnergyMesh(indOffset+ie) = EnergyMesh(indOffset)+    &
!              real(ie,kind=RealKind)*cmplx(dImEo,ZERO,kind=CmplxKind)
         EnergyOffset(NumOffsetE/2+1+i) = & 
               cmplx(ErTop,aimag(EnergyMesh(Offset)),kind=CmplxKind)+ &
               real(i,kind=RealKind)*cmplx(d_rImE,ZERO,kind=CmplxKind)
      enddo
   else
      d_rImE = 0.8d0*aimag(EnergyMesh(Offset))/ &
               real((NumOffsetE+1)/2,kind=RealKind)
      do i = -(NumOffsetE-1)/2,(NumOffsetE-1)/2
         EnergyOffset((NumOffsetE-1)/2+1+i) = & 
               cmplx(ErTop,aimag(EnergyMesh(Offset)),kind=CmplxKind)+ &
               real(i,kind=RealKind)*cmplx(d_rImE,ZERO,kind=CmplxKind)
      enddo
      i= (NumOffsetE+1)/2
      EnergyOffset(NumOffsetE) = & 
            cmplx(ErTop,aimag(EnergyMesh(Offset)),kind=CmplxKind)+ &
            real(i,kind=RealKind)*cmplx(d_rImE,ZERO,kind=CmplxKind)
   endif
!
   end subroutine setupOffsetE
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEOffsetall() result(ep)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), pointer :: ep(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getEPointOffset','ContourModule is not initialized')
   endif
!
   if (NumKGroups > 0) then
      call ErrorHandler('getEPoint','Invalid function argument(s)')
   else if (NumContours < 1 .or. NumOffsetE < 1) then
      call ErrorHandler('getEPointOffset','NumContours or NumOffsetE < 1')
   endif
   ep => EnergyOffset(1:NumOffsetE)
!
   end function getEOffsetall
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEOffsetatie(ie) result(e)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ie
   complex (kind=CmplxKind) :: e
!
   if (.not.Initialized) then
      call ErrorHandler('getEPointOffset','ContourModule is not initialized')
   endif
!
   if (NumKGroups > 0) then
      call ErrorHandler('getEPointOffset','Invalid interface arguments')
   else if (NumContours < 1) then
      call ErrorHandler('getEPointOffset','NumContours < 1', NumContours)
   else
      if (ie>NumOffsetE .or. ie<1) then
         call ErrorHandler('getEPointOffset','ie out of range',ie)
      else
         e = EnergyOffset(ie)
      endif
   endif
!
   end function getEOffsetatie
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getEWeightOffset() result(w)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind) :: w
!
   if (.not.Initialized) then
      call ErrorHandler('getEWeightOffset','ContourModule is not initialized')
   else if (NumContours < 1) then
      call ErrorHandler('getEWeightOffset','NumContours < 1', NumContours)
   else if (NumKGroups > 0) then
      call ErrorHandler('getEWeightOffset','Energy weight is not defined')
   else if (Offset<1 .or. Offset > NumEs) then
      call ErrorHandler('getEWeightOffset','Energy weight is not defined')
   endif
   w = EnergyWght(Offset)
!
   end function getEWeightOffset
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isHorizontalContour() result(hc)
!  ===================================================================
   implicit none
!
   logical hc
!
   if (ContourType == HorizLine) then
      hc = .true.
   else
      hc = .false.
   endif
!
   end function isHorizontalContour
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printContour()
!  ===================================================================
   implicit none
! 
   integer (kind=IntKind) :: ik
   integer (kind=IntKind) :: ie
   integer (kind=IntKind) :: n
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')'********************************'
   write(6,'( 24x,a )')'*   Output from printContour   *'
   write(6,'(24x,a,/)')'********************************'
!
   if (NumKGroups == 0) then
      write(6,'(''Number of Eneregy Mesh = '',i5)')NumEs
      write(6,'(/,80(''=''))')
      write(6,'(a)')   &
      &   ' Index                 Energy                               Weight'
      write(6,'(80(''-''))')
      do ie=1,NumEs
         write(6,'(i5,5x,2d16.8,5x,2d16.8)')ie,EnergyMesh(ie),EnergyWght(ie)
      enddo
      write(6,'(80(''=''))')
   else
      write(6,'(''Number of Contours = '',i5)')NumContours
      do ik=1,NumContours
         n=EnergyKMesh(ik)%nume
         write(6,'(''K-vector = '',3f10.5)')EnergyKMesh(ik)%kmesh(1:3)
         write(6,'(''Number of Energy Mesh = '',i5)')n
         write(6,'(/,80(''=''))')
         write(6,'(a)')  ' Index                  Energy'
         write(6,'(80(''-''))')
         do ie=1,n
            write(6,'(i5,6x,2d16.8)')ie,EnergyKMesh(ik)%emesh(ie)
         enddo
         write(6,'(80(''=''))')
      enddo
   endif
   end subroutine printContour
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isNotSearchingEf() result(ef)
!  ===================================================================
   implicit none
!
   logical :: ef
!
   ef = .not.SearchingFermiEnergy
!
   end function isNotSearchingEf
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isMatsubaraContour() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if ( ContourType == MatsubaraPoles ) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isMatsubaraContour
!  ===================================================================
end module ContourModule
