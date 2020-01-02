program testPolyFermi
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use SystemModule, only : getNumAtoms, getAtomPosition
!
   use ScfDataModule, only : n_spin_pola, n_spin_cant, Temperature, NumEs
   use ScfDataModule, only : istop, EvBottom, isNonRelativisticCore
   use ScfDataModule, only : eGridType, Contourtype
!
   use PublicParamDefinitionsModule, only : NicholsonPoints
!
   use AtomModule, only : getPhiLmax, getStepFuncLmax, getPotLmax
!
   use MathParamModule, only : ZERO, TEN2m6, ONE, TWO, TEN2m8, PI4, PI
   use MathParamModule, only : CZERO, SQRTm1, CONE
!
   use PhysParamModule, only : Boltzmann
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use SphericalHarmonicsModule, only : endSphericalHarmonics
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use MadelungModule, only : initMadelung, endMadelung
!
   use SystemSymmetryModule, only : initSystemSymmetry, endSystemSymmetry,   &
                                    calSymmetryFlags
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : getVolume, getInscrSphVolume
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use OutputModule, only : getStandardOutputLevel
!
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
!
   use PotentialModule, only : initPotential, endPotential
   use PotentialModule, only : readPotential
!
   use DataServiceCenterModule, only : isDataStorageExisting, &
                                       getDataStorage, RealMark
!
   use ContourModule, only : setupContour, isMatsubaraContour
   use ContourModule, only : getEPoint, getEWeight, getNumEs
   use ContourModule, only : initContour, endContour
!
   implicit   none
!
   integer (kind=IntKind) :: iprint = 0
   integer (kind=IntKind) :: NumAtoms, LocalNumAtoms
   integer (kind=IntKind) :: lmax_max
   integer (kind=IntKind) :: i, id, ig, ne, ie, NumEs_don
   integer (kind=IntKind), parameter :: MaxNumEs_don = 100000
!
   integer (kind=IntKind), allocatable :: atom_print_level(:)
   integer (kind=IntKind), allocatable :: lmax_pot(:), lmax_step(:)
   integer (kind=IntKind), allocatable :: GlobalIndex(:)
!
   real (kind=RealKind) :: Temp, Chempot, w, fac, E_don, N_don, E_xg, N_xg
   real (kind=RealKind) :: N_band, E_band
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
!
   complex (kind=CmplxKind), allocatable :: xg(:), wg(:)
   complex (kind=CmplxKind) :: Epole_don(MaxNumEs_don), Weight_don(MaxNumEs_don)
   complex (kind=CmplxKind) :: c, d, gf
   complex (kind=CmplxKind), pointer :: ep(:), ew(:)
!
   interface
      subroutine calNicholsonPoles(etopcor,ebot,chempot,temperature,  &
                                   nume,epole,weight,iprint,istop)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         character (len=*), intent(in) :: istop
         real (kind=RealKind), intent(in) :: etopcor,ebot,chempot,temperature
         integer (kind=IntKind), intent(in) :: iprint
         integer (kind=IntKind), intent(out) :: nume
         complex (kind=CmplxKind), intent(out) :: epole(:), weight(:)
      end subroutine calNicholsonPoles
   end interface
!
   interface
      subroutine polyfermi(Temp,xg,wg,ne,npoles,mu)
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         implicit none
         integer (kind=IntKind), intent(in) :: ne
         integer (kind=IntKind), intent(in), optional :: npoles
         real (kind=RealKind), intent(in) :: Temp
         real (kind=RealKind), intent(in), optional :: mu
         complex (kind=CmplxKind), intent(out) :: xg(ne),wg(ne)
      end subroutine polyfermi
   end interface
!
!  -------------------------------------------------------------------
   call startProcess()
!  NumAtoms = getNumAtoms()
!  LocalNumAtoms=getLocalNumAtoms()
!  -------------------------------------------------------------------
!
!  allocate(atom_print_level(1:LocalNumAtoms))
!  allocate(lmax_step(LocalNumAtoms),lmax_pot(LocalNumAtoms))
!
!  lmax_max = 0
!  do id = 1, LocalNumAtoms
!     lmax_pot(id) = getPotLmax(id)
!     lmax_step(id)  = getStepFuncLmax(id)
!     lmax_max = max(lmax_max,2*getPhiLmax(id),lmax_step(id))
!     atom_print_level(id) = getStandardOutputLevel(id)
!  enddo
!
!  allocate(AtomPosition(1:3,1:NumAtoms))
!  do ig=1,NumAtoms
!     AtomPosition(1:3,ig)=getAtomPosition(ig)
!  enddo
!  allocate(GlobalIndex(LocalNumAtoms))
!  do id=1,LocalNumAtoms
!     GlobalIndex(id)=getGlobalIndex(id)
!  enddo
!
!  -------------------------------------------------------------------
!  call initSphericalHarmonics(2*lmax_max)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initilize GauntFactors:
!  ===================================================================
!  call initGauntFactors(lmax_max,istop,iprint)
!  -------------------------------------------------------------------
!
!  -------------------------------------------------------------------
!  call setupRadGridAndCell(LocalNumAtoms,lmax_max)
!  -------------------------------------------------------------------
!
!  if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
!     bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
!  else
!     ----------------------------------------------------------------
!     call ErrorHandler('testSineMatrixPole','Bravais vector data does not exist')
!     ----------------------------------------------------------------
!  endif
!
!  -------------------------------------------------------------------
!  call initMadelung(LocalNumAtoms,NumAtoms,GlobalIndex,              &
!                    maxval(lmax_pot),maxval(lmax_pot),bravais,AtomPosition,0)
!  -------------------------------------------------------------------
!  call initSystemSymmetry(NumAtoms,LocalNumAtoms,lmax_pot,lmax_step, &
!                          atom_print_level)
!  -------------------------------------------------------------------
!  call calSymmetryFlags()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize and read potential data
!  -------------------------------------------------------------------
!  call initPotential(LocalNumAtoms,lmax_pot,lmax_step,               &
!                     n_spin_pola,n_spin_cant,istop,atom_print_level)
!  -------------------------------------------------------------------
!  call readPotential()
!  -------------------------------------------------------------------
!
   Temp = Temperature
   Chempot = 0.65d0
!
!  -------------------------------------------------------------------
   call calNicholsonPoles(-3.0d0,Chempot-1.5d0,Chempot,Temp,          &
                          NumEs_don,Epole_don,Weight_don,0,'none')
!  -------------------------------------------------------------------
!
   ne = NumEs
   write(6,'(/,a,f15.5)')'Temperature = ',Temp
   write(6,'(a,i5)')'Number of Gaussian points on the contour = ',ne
   allocate(xg(ne),wg(ne))
!  -------------------------------------------------------------------
   call polyfermi(Temp,xg,wg,ne,mu=Chempot)
!  -------------------------------------------------------------------
   do ie = 1, ne
      write(6,'(a,i5,2x,2d15.7,2x,2d15.7)')'ie, xg, wg = ',ie,xg(ie),wg(ie)
   enddo
!
!  ===================================================================
!  Testing the band energy of an 1-D tight-binding model, whose Green 
!  function is given by
!     G(z) = -i * 20/w * 1/sqrt(1-4*(z-mu)^2/w^2)
!  ===================================================================
   w = 1.0d+5*Boltzmann
   fac = 20.0d0/w
!
   c = CZERO
   d = CZERO
   do ie = 1, NumEs_don
      gf = -SQRTm1*fac/sqrt(CONE-(TWO*(Epole_don(ie)-Chempot)/w)**2)
      c = c + Epole_don(ie)*Weight_don(ie)*gf
      d = d + Weight_don(ie)*gf
   enddo
   N_don = -aimag(d)/PI
   E_don = -aimag(c)/PI
!
   c = CZERO
   d = CZERO
   do ie = 1, ne
!     ================================================================
!     Note: If calling polyfermi is made as follows:
!        call polyfermi(Temp,xg,wg,ne)
!     xg is relative to mu=0, the Green function is given by
!        G(z) = -i * 20/w * 1/sqrt(1-4*z^2/w^2)
!     and the band energy should be cauculated using
!        E_xg = - Im {sum_i (z_i + Chempot)*w_i*G(z_i)}/PI
!     ================================================================
      gf = -SQRTm1*fac/sqrt(CONE-(TWO*(xg(ie)-Chempot)/w)**2)
      c = c + xg(ie)*wg(ie)*gf
      d = d + wg(ie)*gf
   enddo
   N_xg = -aimag(d)/PI
   E_xg = -aimag(c)/PI
!
   write(6,'(/,a,d15.8)')"Nicholson''s method: Band energy = ",E_don
   write(6,'(  a,d15.8)')"Nicholson''s method: No. of elec = ",N_don
   write(6,'(/,a,d15.8)')"Zhang''s method:     Band energy = ",E_xg
   write(6,'(  a,d15.8)')"Zhang''s method:     No. of elec = ",N_xg
!
!  deallocate(atom_print_level,lmax_step,lmax_pot,GlobalIndex,AtomPosition)
!
!  ===================================================================
!  Test the new ContourModule, which uses MatsubaraModule
!  ===================================================================
   if (isMatsubaraContour()) then
do i = 1, 10
      call endContour()
      call initContour( ContourType, eGridType, NumEs, Temperature, 'none', &
                        maxval(atom_print_level))
!     ----------------------------------------------------------------
      call setupContour( Chempot-1.5d0, Chempot, 0.0d0, 0.0d0 )
!     ----------------------------------------------------------------
      ep => getEPoint()
      ew => getEWeight()
      c = CZERO
      d = CZERO
      w = 1.0d+5*Boltzmann
      fac = 20.0d0/w
      do ie = 1, getNumEs()
         if (eGridType == NicholsonPoints .and. ie <= NumEs_don) then
            if (abs(ep(ie)-Epole_don(ie)) > TEN2m8) then
               write(6,'(a,i5,2x,2d15.8,2x,2d15.8)')'ie, weight = ',  &
                     ie,ew(ie),Weight_don(ie)
               call ErrorHandler('testPolyFermi','ep <> Epole_don',   &
                                 ep(ie),Epole_don(ie))
            else if (abs(ew(ie)-Weight_don(ie)) > TEN2m8) then
               write(6,'(a,i5,2x,2d15.8,2x,2d15.8)')'ie, epoint = ',  &
                     ie,ep(ie),Epole_don(ie)
               call ErrorHandler('testPolyFermi','ew <> Weight_don',  &
                                 ew(ie),Weight_don(ie))
            endif
         else if (ie <= ne) then
            if (abs(ep(ie)-xg(ie)) > TEN2m8) then
               write(6,'(a,i5,2x,2d15.8,2x,2d15.8)')'ie, weight = ',  &
                     ie,ew(ie),wg(ie)
               call ErrorHandler('testPolyFermi','ep <> xg',ep(ie),xg(ie))
            else if (abs(ew(ie)-wg(ie)) > TEN2m8) then
               write(6,'(a,i5,2x,2d15.8,2x,2d15.8)')'ie, epoint = ',  &
                     ie,ep(ie),wg(ie)
               call ErrorHandler('testPolyFermi','ew <> wg',ew(ie),wg(ie))
            endif
         endif
         gf = -SQRTm1*fac/sqrt(CONE-(TWO*(ep(ie)-Chempot)/w)**2)
         c = c + ep(ie)*ew(ie)*gf
         d = d + ew(ie)*gf
      enddo
      N_band = -aimag(d)/PI
      E_band = -aimag(c)/PI
!
      write(6,'(/,a,d15.8)')"Band energy = ",E_band
      write(6,'(  a,d15.8)')"No. of elec = ",N_band
!
enddo
      nullify(ep,ew)
   endif
!
   deallocate(xg,wg)
!
!  -------------------------------------------------------------------
!  call endPotential()
!  call endSystemSymmetry()
!  call endMadelung()
!  call endGauntFactors()
!  call endSphericalHarmonics()
   call finishProcess()
!  -------------------------------------------------------------------
   stop 'Ok'
!
end program testPolyFermi
