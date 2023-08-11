module ChemElementModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
public :: getAtomInfo,      &
          getName,          &
          getZtot,          &
          getZcor,          &
          getZsem,          &
          getZval,          &
          setConfiguration, &
          setNumCoreStates, &
          getNumCoreStates, &
          getCoreStateN,    &
          getCoreStateL,    &
          getCoreStateKappa,&
          getCoreStateIndex,&
          getCoreStateSymbol,&
          getAtomicRadius,  &
          getAtomicMass,    &
          getImplicitMuffinTinRadius,  &
          getImplicitCoreRadius,  &
          getDebyeTemperature, &
          isAtomName ! check if the given string a chemical element name.
!
   interface setNumCoreStates
      module procedure setNumCoreStates_a, setNumCoreStates_n
   end interface
!
   interface setConfiguration
      module procedure setConfig_a, setConfig_n
   end interface
!
   interface getZcor
      module procedure getZcor_a, getZcor_n
   end interface
!
   interface getZsem
      module procedure getZsem_a, getZsem_n
   end interface
!
   interface getZval
      module procedure getZval_a, getZval_n
   end interface
!
   interface getNumCoreStates
      module procedure getNumCoreStates_a, getNumCoreStates_n
   end interface
!
   interface getCoreStateN
      module procedure getCoreStateN_a, getCoreStateN_n
   end interface
!
   interface getCoreStateL
      module procedure getCoreStateL_a, getCoreStateL_n
   end interface
!
   interface getCoreStateKappa
      module procedure getCoreStateKappa_a, getCoreStateKappa_n
   end interface
!
   interface getCoreStateSymbol
      module procedure getCoreStateSymbol_a, getCoreStateSymbol_n
   end interface
!
   interface getDebyeTemperature
      module procedure getDebyeTemperature_a, getDebyeTemperature_n
   end interface
!
   interface getAtomicRadius
      module procedure getAtomicRadius_a, getAtomicRadius_n
   end interface
!
   interface getAtomicMass
      module procedure getAtomicMass_a, getAtomicMass_n
   end interface
!
   interface getImplicitMuffinTinRadius
      module procedure getImplicitMuffinTinRadius_a, getImplicitMuffinTinRadius_n
   end interface
!
   interface getImplicitCoreRadius
      module procedure getImplicitCoreRadius_a, getImplicitCoreRadius_n
   end interface
!
   integer (kind=IntKind), parameter, public :: MaxLenOfAtomName = 16
!
   integer (kind=IntKind), parameter, public :: MaxNumc = 26
!
private
   type ElementStruc
      character (len=MaxLenOfAtomName) :: AtomName
      integer (kind=IntKind) :: AtomicNumber
      integer (kind=IntKind) :: NumDeepCoreElectrons
      integer (kind=IntKind) :: NumSemiCoreElectrons
      integer (kind=IntKind) :: NumValenceElectrons
      integer (kind=IntKind) :: NumCoreStates
      integer (kind=IntKind) :: NumVariations
      real (kind=RealKind) :: DebyeT
      real (kind=RealKind) :: AtomicRadius
      real (kind=RealKind) :: AtomicMass
      real (kind=RealKind) :: ImplicitMuffinTinRadius
      real (kind=RealKind) :: ImplicitCoreRadius
      type (ElementStruc), pointer :: Variation(:)
!     integer (kind=IntKind), pointer :: QuantumNumber_n(:)
!     integer (kind=IntKind), pointer :: QuantumNumber_l(:)
!     integer (kind=IntKind), pointer :: QuantumNumber_kappa(:)
   end type ElementStruc
!
   integer (kind=IntKind), parameter :: NumElements = 103
   integer (kind=IntKind), parameter :: MinZtot = -1
!
   type (ElementStruc), save :: Element(MinZtot:NumElements)
!
   logical :: Initialized = .false.
!
   type AtomicLevel
      character (len=2) :: Symbol
      integer (kind=IntKind) :: n
      integer (kind=IntKind) :: l
      integer (kind=IntKind) :: kappa
   end type AtomicLevel
!
   type (AtomicLevel) :: CoreState(1:MaxNumc)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initChemElement()
!  ===================================================================
   use MathParamModule, only : ZERO
   use PhysParamModule, only : Angstrom2Bohr
!
   implicit none
!
!  ===================================================================
!
   CoreState( 1)%Symbol = '1s'
   CoreState( 1)%n = 1
   CoreState( 1)%l = 0
   CoreState( 1)%kappa = -1
!
   CoreState( 2)%Symbol = '2s'
   CoreState( 2)%n = 2
   CoreState( 2)%l = 0
   CoreState( 2)%kappa = -1
!
   CoreState( 3)%Symbol = '2p'
   CoreState( 3)%n = 2
   CoreState( 3)%l = 1
   CoreState( 3)%kappa = 1
!
   CoreState( 4)%Symbol = '2p'
   CoreState( 4)%n = 2
   CoreState( 4)%l = 1
   CoreState( 4)%kappa = -2
!
   CoreState( 5)%Symbol = '3s'
   CoreState( 5)%n = 3
   CoreState( 5)%l = 0
   CoreState( 5)%kappa = -1
!
   CoreState( 6)%Symbol = '3p'
   CoreState( 6)%n = 3
   CoreState( 6)%l = 1
   CoreState( 6)%kappa = 1
!
   CoreState( 7)%Symbol = '3p'
   CoreState( 7)%n = 3
   CoreState( 7)%l = 1
   CoreState( 7)%kappa = -2
!
   CoreState( 8)%Symbol = '3d'
   CoreState( 8)%n = 3
   CoreState( 8)%l = 2
   CoreState( 8)%kappa = 2
!
   CoreState(9)%Symbol = '3d'
   CoreState(9)%n = 3
   CoreState(9)%l = 2
   CoreState(9)%kappa = -3
!
   CoreState(10)%Symbol = '4s'
   CoreState(10)%n = 4
   CoreState(10)%l = 0
   CoreState(10)%kappa = -1
!
   CoreState(11)%Symbol = '4p'
   CoreState(11)%n = 4
   CoreState(11)%l = 1
   CoreState(11)%kappa = 1
!
   CoreState(12)%Symbol = '4p'
   CoreState(12)%n = 4
   CoreState(12)%l = 1
   CoreState(12)%kappa = -2
!
   CoreState(13)%Symbol = '4d'
   CoreState(13)%n = 4
   CoreState(13)%l = 2
   CoreState(13)%kappa = 2
!
   CoreState(14)%Symbol = '4d'
   CoreState(14)%n = 4
   CoreState(14)%l = 2
   CoreState(14)%kappa = -3
!
   CoreState(15)%Symbol = '4f'
   CoreState(15)%n = 4
   CoreState(15)%l = 3
   CoreState(15)%kappa = 3
!
   CoreState(16)%Symbol = '4f'
   CoreState(16)%n = 4
   CoreState(16)%l = 3
   CoreState(16)%kappa = -4
!
   CoreState(17)%Symbol = '5s'
   CoreState(17)%n = 5
   CoreState(17)%l = 0
   CoreState(17)%kappa = -1
!
   CoreState(18)%Symbol = '5p'
   CoreState(18)%n = 5
   CoreState(18)%l = 1
   CoreState(18)%kappa = 1
!
   CoreState(19)%Symbol = '5p'
   CoreState(19)%n = 5
   CoreState(19)%l = 1
   CoreState(19)%kappa = -2
!
   CoreState(20)%Symbol = '5d'
   CoreState(20)%n = 5
   CoreState(20)%l = 2
   CoreState(20)%kappa = 2
!
   CoreState(21)%Symbol = '5d'
   CoreState(21)%n = 5
   CoreState(21)%l = 2
   CoreState(21)%kappa = -3
!
   CoreState(22)%Symbol = '5f'
   CoreState(22)%n = 5
   CoreState(22)%l = 3
   CoreState(22)%kappa = 3
!
   CoreState(23)%Symbol = '5f'
   CoreState(23)%n = 5
   CoreState(23)%l = 3
   CoreState(23)%kappa = -4
!
   CoreState(24)%Symbol = '6s'
   CoreState(24)%n = 6
   CoreState(24)%l = 0
   CoreState(24)%kappa = -1
!
   CoreState(25)%Symbol = '6p'
   CoreState(25)%n = 6
   CoreState(25)%l = 1
   CoreState(25)%kappa = 1
!
   CoreState(26)%Symbol = '6p'
   CoreState(26)%n = 6
   CoreState(26)%l = 1
   CoreState(26)%kappa = -2
!
!  ===================================================================
!  Note: Atomic radius values are mostly taken from the "empirical" 
!    column of the table in Atomic radii of the elements (data page) page
!    (except for few that no data available in the "empirical" column)
!    https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
!  where the empirical data in the page are taken from 
!    J.C. Slater (1964). "Atomic Radii in Crystals". J. Chem. Phys. 41: 3199.
!    doi:10.1063/1.1725697
!  ===================================================================
!
   Element( MinZtot)%AtomName='CPA'
   Element( MinZtot)%AtomicNumber=MinZtot
   Element( MinZtot)%NumDeepCoreElectrons=0
   Element( MinZtot)%NumSemiCoreElectrons=0
   Element( MinZtot)%NumValenceElectrons=0
   Element( MinZtot)%NumCoreStates=0
   Element( MinZtot)%NumVariations=0
   Element( MinZtot)%DebyeT=ZERO
   Element( MinZtot)%AtomicRadius=ZERO
!
   Element(  0)%AtomName='Va'
   Element(  0)%AtomicNumber=0
   Element(  0)%NumDeepCoreElectrons=0
   Element(  0)%NumSemiCoreElectrons=0
   Element(  0)%NumValenceElectrons=0
   Element(  0)%NumCoreStates=0
   Element(  0)%NumVariations=0
   Element(  0)%DebyeT=ZERO
   Element(  0)%AtomicRadius=ZERO
   Element(  0)%AtomicMass=ZERO
!
   Element(  1)%AtomName='H '
   Element(  1)%AtomicNumber=1
   Element(  1)%NumDeepCoreElectrons=0
   Element(  1)%NumSemiCoreElectrons=0
   Element(  1)%NumValenceElectrons=1
   Element(  1)%NumCoreStates=0
   Element(  1)%NumVariations=0
   Element(  1)%DebyeT=ZERO
   Element(  1)%AtomicRadius=0.25d0*Angstrom2Bohr
   Element(  1)%AtomicMass=1.00797d0
   Element(  1)%ImplicitMuffinTinRadius=0.75d0*Element(  1)%AtomicRadius
   Element(  1)%ImplicitCoreRadius=0.75d0*Element(  1)%AtomicRadius
!
   Element(  2)%AtomName='He'
   Element(  2)%AtomicNumber=2
   Element(  2)%NumDeepCoreElectrons=2
   Element(  2)%NumSemiCoreElectrons=0
   Element(  2)%NumValenceElectrons=0
   Element(  2)%NumCoreStates=1
   Element(  2)%NumVariations=0
   Element(  2)%DebyeT=ZERO
   Element(  2)%AtomicRadius=1.20d0*Angstrom2Bohr
   Element(  2)%AtomicMass=4.00260d0
   Element(  2)%ImplicitMuffinTinRadius=0.75d0*Element(  2)%AtomicRadius
   Element(  2)%ImplicitCoreRadius=0.75d0*Element(  2)%AtomicRadius
!
   Element(  3)%AtomName='Li'
   Element(  3)%AtomicNumber=3
   Element(  3)%NumDeepCoreElectrons=2
   Element(  3)%NumSemiCoreElectrons=0
   Element(  3)%NumValenceElectrons=1
   Element(  3)%NumCoreStates=1
   Element(  3)%NumVariations=0
   Element(  3)%DebyeT=344.7d0
   Element(  3)%AtomicRadius=1.45d0*Angstrom2Bohr
   Element(  3)%AtomicMass=6.941d0
   Element(  3)%ImplicitMuffinTinRadius=0.75d0*Element(  3)%AtomicRadius
   Element(  3)%ImplicitCoreRadius=0.75d0*Element(  3)%AtomicRadius
!
   Element(  4)%AtomName='Be'
   Element(  4)%AtomicNumber=4
   Element(  4)%NumDeepCoreElectrons=2
   Element(  4)%NumSemiCoreElectrons=0
   Element(  4)%NumValenceElectrons=2
   Element(  4)%NumCoreStates=1
   Element(  4)%NumVariations=0
   Element(  4)%DebyeT=ZERO
   Element(  4)%AtomicRadius=1.05d0*Angstrom2Bohr
   Element(  4)%AtomicMass=9.01218d0
   Element(  4)%ImplicitMuffinTinRadius=0.75d0*Element(  4)%AtomicRadius
   Element(  4)%ImplicitCoreRadius=0.75d0*Element(  4)%AtomicRadius
!
   Element(  5)%AtomName='B '
   Element(  5)%AtomicNumber=5
   Element(  5)%NumDeepCoreElectrons=2
   Element(  5)%NumSemiCoreElectrons=0
   Element(  5)%NumValenceElectrons=3
   Element(  5)%NumCoreStates=1
   Element(  5)%NumVariations=0
   Element(  5)%DebyeT=ZERO
   Element(  5)%AtomicRadius=0.85d0*Angstrom2Bohr
   Element(  5)%AtomicMass=10.81d0
   Element(  5)%ImplicitMuffinTinRadius=0.75d0*Element(  5)%AtomicRadius
   Element(  5)%ImplicitCoreRadius=0.75d0*Element(  5)%AtomicRadius
!
   Element(  6)%AtomName='C '
   Element(  6)%AtomicNumber=6
   Element(  6)%NumDeepCoreElectrons=2
   Element(  6)%NumSemiCoreElectrons=0
   Element(  6)%NumValenceElectrons=4
   Element(  6)%NumCoreStates=1
   Element(  6)%NumVariations=0
   Element(  6)%DebyeT=ZERO
   Element(  6)%AtomicRadius=0.70d0*Angstrom2Bohr
   Element(  6)%AtomicMass=12.011d0
   Element(  6)%ImplicitMuffinTinRadius=0.75d0*Element(  6)%AtomicRadius
   Element(  6)%ImplicitCoreRadius=0.75d0*Element(  6)%AtomicRadius
!
   Element(  7)%AtomName='N '
   Element(  7)%AtomicNumber=7
   Element(  7)%NumDeepCoreElectrons=2
   Element(  7)%NumSemiCoreElectrons=0
   Element(  7)%NumValenceElectrons=5
   Element(  7)%NumCoreStates=1
   Element(  7)%NumVariations=0
   Element(  7)%DebyeT=ZERO
   Element(  7)%AtomicRadius=0.65d0*Angstrom2Bohr
   Element(  7)%AtomicMass=14.0067d0
   Element(  7)%ImplicitMuffinTinRadius=0.75d0*Element(  7)%AtomicRadius
   Element(  7)%ImplicitCoreRadius=0.75d0*Element(  7)%AtomicRadius
!
   Element(  8)%AtomName='O '
   Element(  8)%AtomicNumber=8
   Element(  8)%NumDeepCoreElectrons=2
   Element(  8)%NumSemiCoreElectrons=0
   Element(  8)%NumValenceElectrons=6
   Element(  8)%NumCoreStates=1
   Element(  8)%NumVariations=0
   Element(  8)%DebyeT=ZERO
   Element(  8)%AtomicRadius=0.60d0*Angstrom2Bohr
   Element(  8)%AtomicMass=15.9994d0
   Element(  8)%ImplicitMuffinTinRadius=0.75d0*Element(  8)%AtomicRadius
   Element(  8)%ImplicitCoreRadius=0.75d0*Element(  8)%AtomicRadius
!
   Element(  9)%AtomName='F '
   Element(  9)%AtomicNumber=9
   Element(  9)%NumDeepCoreElectrons=2
   Element(  9)%NumSemiCoreElectrons=0
   Element(  9)%NumValenceElectrons=7
   Element(  9)%NumCoreStates=1
   Element(  9)%NumVariations=0
   Element(  9)%DebyeT=ZERO
   Element(  9)%AtomicRadius=0.50d0*Angstrom2Bohr
   Element(  9)%AtomicMass=18.998403d0
   Element(  9)%ImplicitMuffinTinRadius=0.75d0*Element(  9)%AtomicRadius
   Element(  9)%ImplicitCoreRadius=0.75d0*Element(  9)%AtomicRadius
!
   Element( 10)%AtomName='Ne'
   Element( 10)%AtomicNumber=10
   Element( 10)%NumDeepCoreElectrons=2
   Element( 10)%NumSemiCoreElectrons=8
   Element( 10)%NumValenceElectrons=0
   Element( 10)%NumCoreStates=4
   Element( 10)%NumVariations=0
   Element( 10)%DebyeT=ZERO
   Element( 10)%AtomicRadius=1.60d0*Angstrom2Bohr
   Element( 10)%AtomicMass=20.179d0
   Element( 10)%ImplicitMuffinTinRadius=0.75d0*Element( 10)%AtomicRadius
   Element( 10)%ImplicitCoreRadius=0.75d0*Element( 10)%AtomicRadius
!
   Element( 11)%AtomName='Na'
   Element( 11)%AtomicNumber=11
   Element( 11)%NumDeepCoreElectrons=2
   Element( 11)%NumSemiCoreElectrons=8
   Element( 11)%NumValenceElectrons=1
   Element( 11)%NumCoreStates=4
   Element( 11)%NumVariations=0
   Element( 11)%DebyeT=162.6d0
   Element( 11)%AtomicRadius=1.80d0*Angstrom2Bohr
   Element( 11)%AtomicMass=22.98977d0
   Element( 11)%ImplicitMuffinTinRadius=0.75d0*Element( 11)%AtomicRadius
   Element( 11)%ImplicitCoreRadius=0.75d0*Element( 11)%AtomicRadius
!
   Element( 12)%AtomName='Mg'
   Element( 12)%AtomicNumber=12
   Element( 12)%NumDeepCoreElectrons=2
   Element( 12)%NumSemiCoreElectrons=8
   Element( 12)%NumValenceElectrons=2
   Element( 12)%NumCoreStates=4
   Element( 12)%NumVariations=0
   Element( 12)%DebyeT=ZERO
   Element( 12)%AtomicRadius=1.50d0*Angstrom2Bohr
   Element( 12)%AtomicMass=24.305d0
   Element( 12)%ImplicitMuffinTinRadius=0.75d0*Element( 12)%AtomicRadius
   Element( 12)%ImplicitCoreRadius=0.75d0*Element( 12)%AtomicRadius
!
   Element( 13)%AtomName='Al'
   Element( 13)%AtomicNumber=13
   Element( 13)%NumDeepCoreElectrons=2
   Element( 13)%NumSemiCoreElectrons=8
   Element( 13)%NumValenceElectrons=3
   Element( 13)%NumCoreStates=4
   Element( 13)%NumVariations=0
   Element( 13)%DebyeT=393.6d0
   Element( 13)%AtomicRadius=1.25d0*Angstrom2Bohr
   Element( 13)%AtomicMass=26.98154d0
   Element( 13)%ImplicitMuffinTinRadius=0.75d0*Element( 13)%AtomicRadius
   Element( 13)%ImplicitCoreRadius=0.75d0*Element( 13)%AtomicRadius
!
   Element( 14)%AtomName='Si'
   Element( 14)%AtomicNumber=14
   Element( 14)%NumDeepCoreElectrons=2
   Element( 14)%NumSemiCoreElectrons=8
   Element( 14)%NumValenceElectrons=4
   Element( 14)%NumCoreStates=4
   Element( 14)%NumVariations=0
   Element( 14)%DebyeT=ZERO
   Element( 14)%AtomicRadius=1.10d0*Angstrom2Bohr
   Element( 14)%AtomicMass=28.0855d0
   Element( 14)%ImplicitMuffinTinRadius=0.75d0*Element( 14)%AtomicRadius
   Element( 14)%ImplicitCoreRadius=0.75d0*Element( 14)%AtomicRadius
!
   Element( 15)%AtomName='P '
   Element( 15)%AtomicNumber=15
   Element( 15)%NumDeepCoreElectrons=2
   Element( 15)%NumSemiCoreElectrons=8
   Element( 15)%NumValenceElectrons=5
   Element( 15)%NumCoreStates=4
   Element( 15)%NumVariations=0
   Element( 15)%DebyeT=ZERO
   Element( 15)%AtomicRadius=1.00d0*Angstrom2Bohr
   Element( 15)%AtomicMass=30.97376d0
   Element( 15)%ImplicitMuffinTinRadius=0.75d0*Element( 15)%AtomicRadius
   Element( 15)%ImplicitCoreRadius=0.75d0*Element( 15)%AtomicRadius
!
   Element( 16)%AtomName='S '
   Element( 16)%AtomicNumber=16
   Element( 16)%NumDeepCoreElectrons=2
   Element( 16)%NumSemiCoreElectrons=8
   Element( 16)%NumValenceElectrons=6
   Element( 16)%NumCoreStates=4
   Element( 16)%NumVariations=0
   Element( 16)%DebyeT=ZERO
   Element( 16)%AtomicRadius=1.00d0*Angstrom2Bohr
   Element( 16)%AtomicMass=32.06d0
   Element( 16)%ImplicitMuffinTinRadius=0.75d0*Element( 16)%AtomicRadius
   Element( 16)%ImplicitCoreRadius=0.75d0*Element( 16)%AtomicRadius
!
   Element( 17)%AtomName='Cl'
   Element( 17)%AtomicNumber=17
   Element( 17)%NumDeepCoreElectrons=2
   Element( 17)%NumSemiCoreElectrons=8
   Element( 17)%NumValenceElectrons=7
   Element( 17)%NumCoreStates=4
   Element( 17)%NumVariations=0
   Element( 17)%DebyeT=ZERO
   Element( 17)%AtomicRadius=1.00d0*Angstrom2Bohr
   Element( 17)%AtomicMass=35.453d0
   Element( 17)%ImplicitMuffinTinRadius=0.75d0*Element( 17)%AtomicRadius
   Element( 17)%ImplicitCoreRadius=0.75d0*Element( 17)%AtomicRadius
!
   Element( 18)%AtomName='Ar'
   Element( 18)%AtomicNumber=18
   Element( 18)%NumDeepCoreElectrons=10
   Element( 18)%NumSemiCoreElectrons=8
   Element( 18)%NumValenceElectrons=0
   Element( 18)%NumCoreStates=7
   Element( 18)%NumVariations=0
   Element( 18)%DebyeT=ZERO
   Element( 18)%AtomicRadius=0.71d0*Angstrom2Bohr
   Element( 18)%AtomicMass=39.948d0
   Element( 18)%ImplicitMuffinTinRadius=0.75d0*Element( 18)%AtomicRadius
   Element( 18)%ImplicitCoreRadius=0.75d0*Element( 18)%AtomicRadius
!
   Element( 19)%AtomName='K '
   Element( 19)%AtomicNumber=19
   Element( 19)%NumDeepCoreElectrons=10
   Element( 19)%NumSemiCoreElectrons=8
   Element( 19)%NumValenceElectrons=1
   Element( 19)%NumCoreStates=7
   Element( 19)%NumVariations=0
   Element( 19)%DebyeT=98.8d0
   Element( 19)%AtomicRadius=2.20d0*Angstrom2Bohr
   Element( 19)%AtomicMass=39.0983d0
   Element( 19)%ImplicitMuffinTinRadius=0.75d0*Element( 19)%AtomicRadius
   Element( 19)%ImplicitCoreRadius=0.75d0*Element( 19)%AtomicRadius
!
   Element( 20)%AtomName='Ca'
   Element( 20)%AtomicNumber=20
   Element( 20)%NumDeepCoreElectrons=10
   Element( 20)%NumSemiCoreElectrons=8
   Element( 20)%NumValenceElectrons=2
   Element( 20)%NumCoreStates=7
   Element( 20)%NumVariations=0
   Element( 20)%DebyeT=190.8d0
   Element( 20)%AtomicRadius=1.80d0*Angstrom2Bohr
   Element( 20)%AtomicMass=40.078d0
   Element( 20)%ImplicitMuffinTinRadius=0.75d0*Element( 20)%AtomicRadius
   Element( 20)%ImplicitCoreRadius=0.75d0*Element( 20)%AtomicRadius
!
   Element( 21)%AtomName='Sc'
   Element( 21)%AtomicNumber=21
   Element( 21)%NumDeepCoreElectrons=10
   Element( 21)%NumSemiCoreElectrons=8
   Element( 21)%NumValenceElectrons=3
   Element( 21)%NumCoreStates=7
   Element( 21)%NumVariations=0
   Element( 21)%DebyeT=273.7d0
   Element( 21)%AtomicRadius=1.60d0*Angstrom2Bohr
   Element( 21)%AtomicMass=44.956d0
   Element( 21)%ImplicitMuffinTinRadius=0.75d0*Element( 21)%AtomicRadius
   Element( 21)%ImplicitCoreRadius=0.75d0*Element( 21)%AtomicRadius
!
   Element( 22)%AtomName='Ti'
   Element( 22)%AtomicNumber=22
   Element( 22)%NumDeepCoreElectrons=10
   Element( 22)%NumSemiCoreElectrons=8
   Element( 22)%NumValenceElectrons=4
   Element( 22)%NumCoreStates=7
   Element( 22)%NumVariations=0
   Element( 22)%DebyeT=347.9d0
   Element( 22)%AtomicRadius=1.40d0*Angstrom2Bohr
   Element( 22)%AtomicMass=47.867d0
   Element( 22)%ImplicitMuffinTinRadius=0.75d0*Element( 22)%AtomicRadius
   Element( 22)%ImplicitCoreRadius=0.75d0*Element( 22)%AtomicRadius
!
   Element( 23)%AtomName='V '
   Element( 23)%AtomicNumber=23
   Element( 23)%NumDeepCoreElectrons=10
   Element( 23)%NumSemiCoreElectrons=8
   Element( 23)%NumValenceElectrons=5
   Element( 23)%NumCoreStates=7
   Element( 23)%NumVariations=0
   Element( 23)%DebyeT=432.7d0
   Element( 23)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 23)%AtomicMass=50.942d0
   Element( 23)%ImplicitMuffinTinRadius=0.75d0*Element( 23)%AtomicRadius
   Element( 23)%ImplicitCoreRadius=0.75d0*Element( 23)%AtomicRadius
!
   Element( 24)%AtomName='Cr'
   Element( 24)%AtomicNumber=24
   Element( 24)%NumDeepCoreElectrons=10
   Element( 24)%NumSemiCoreElectrons=8
   Element( 24)%NumValenceElectrons=6
   Element( 24)%NumCoreStates=7
   Element( 24)%NumVariations=0
   Element( 24)%DebyeT=500.7d0
   Element( 24)%AtomicRadius=1.40d0*Angstrom2Bohr
   Element( 24)%AtomicMass=51.996d0
   Element( 24)%ImplicitMuffinTinRadius=0.75d0*Element( 24)%AtomicRadius
   Element( 24)%ImplicitCoreRadius=0.75d0*Element( 24)%AtomicRadius
!
   Element( 25)%AtomName='Mn'
   Element( 25)%AtomicNumber=25
   Element( 25)%NumDeepCoreElectrons=10
   Element( 25)%NumSemiCoreElectrons=8
   Element( 25)%NumValenceElectrons=7
   Element( 25)%NumCoreStates=7
   Element( 25)%NumVariations=0
   Element( 25)%DebyeT=ZERO
   Element( 25)%AtomicRadius=1.40d0*Angstrom2Bohr
   Element( 25)%AtomicMass=54.938d0
   Element( 25)%ImplicitMuffinTinRadius=0.75d0*Element( 25)%AtomicRadius
   Element( 25)%ImplicitCoreRadius=0.75d0*Element( 25)%AtomicRadius
!
   Element( 26)%AtomName='Fe'
   Element( 26)%AtomicNumber=26
   Element( 26)%NumDeepCoreElectrons=10
   Element( 26)%NumSemiCoreElectrons=8
   Element( 26)%NumValenceElectrons=8
   Element( 26)%NumCoreStates=7
   Element( 26)%NumVariations=0
   Element( 26)%DebyeT=442.2d0
   Element( 26)%AtomicRadius=1.40d0*Angstrom2Bohr
   Element( 26)%AtomicMass=55.845d0
   Element( 26)%ImplicitMuffinTinRadius=0.75d0*Element( 26)%AtomicRadius
   Element( 26)%ImplicitCoreRadius=0.75d0*Element( 26)%AtomicRadius
!
   Element( 27)%AtomName='Co'
   Element( 27)%AtomicNumber=27
   Element( 27)%NumDeepCoreElectrons=10
   Element( 27)%NumSemiCoreElectrons=8
   Element( 27)%NumValenceElectrons=9
   Element( 27)%NumCoreStates=7
   Element( 27)%NumVariations=0
   Element( 27)%DebyeT=422.1d0
   Element( 27)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 27)%AtomicMass=58.933d0
   Element( 27)%ImplicitMuffinTinRadius=0.75d0*Element( 27)%AtomicRadius
   Element( 27)%ImplicitCoreRadius=0.75d0*Element( 27)%AtomicRadius
!
   Element( 28)%AtomName='Ni'
   Element( 28)%AtomicNumber=28
   Element( 28)%NumDeepCoreElectrons=10
   Element( 28)%NumSemiCoreElectrons=8
   Element( 28)%NumValenceElectrons=10
   Element( 28)%NumCoreStates=7
   Element( 28)%NumVariations=2
   Element( 28)%DebyeT=422.9d0
   Element( 28)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 28)%AtomicMass=58.693d0
   Element( 28)%ImplicitMuffinTinRadius=0.75d0*Element( 28)%AtomicRadius
   Element( 28)%ImplicitCoreRadius=0.75d0*Element( 28)%AtomicRadius
   allocate( Element( 28)%Variation(1:Element( 28)%NumVariations) )
!
   Element( 28)%Variation(1)%AtomName='Ni&3p'
   Element( 28)%Variation(1)%AtomicNumber=28
   Element( 28)%Variation(1)%NumDeepCoreElectrons=10
   Element( 28)%Variation(1)%NumSemiCoreElectrons=2
   Element( 28)%Variation(1)%NumValenceElectrons=16
   Element( 28)%Variation(1)%NumCoreStates=5
   Element( 28)%Variation(1)%NumVariations=0
!
   Element( 28)%Variation(2)%AtomName='Ni&3p3s'
   Element( 28)%Variation(2)%AtomicNumber=28
   Element( 28)%Variation(2)%NumDeepCoreElectrons=10
   Element( 28)%Variation(2)%NumSemiCoreElectrons=0
   Element( 28)%Variation(2)%NumValenceElectrons=18
   Element( 28)%Variation(2)%NumCoreStates=4
   Element( 28)%Variation(2)%NumVariations=0
!
   Element( 29)%AtomName='Cu'
   Element( 29)%AtomicNumber=29
   Element( 29)%NumDeepCoreElectrons=10
   Element( 29)%NumSemiCoreElectrons=8
   Element( 29)%NumValenceElectrons=11
   Element( 29)%NumCoreStates=7
   Element( 29)%NumVariations=0
   Element( 29)%DebyeT=369.0d0
   Element( 29)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 29)%AtomicMass=63.546d0
   Element( 29)%ImplicitMuffinTinRadius=0.75d0*Element( 29)%AtomicRadius
   Element( 29)%ImplicitCoreRadius=0.75d0*Element( 29)%AtomicRadius
!
   Element( 30)%AtomName='Zn'
   Element( 30)%AtomicNumber=30
   Element( 30)%NumDeepCoreElectrons=10
   Element( 30)%NumSemiCoreElectrons=8
   Element( 30)%NumValenceElectrons=12
   Element( 30)%NumCoreStates=7
   Element( 30)%NumVariations=0
   Element( 30)%DebyeT=ZERO
   Element( 30)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 30)%AtomicMass=65.38d0
   Element( 30)%ImplicitMuffinTinRadius=0.75d0*Element( 30)%AtomicRadius
   Element( 30)%ImplicitCoreRadius=0.75d0*Element( 30)%AtomicRadius
!
   Element( 31)%AtomName='Ga'
   Element( 31)%AtomicNumber=31
   Element( 31)%NumDeepCoreElectrons=10
   Element( 31)%NumSemiCoreElectrons=8
   Element( 31)%NumValenceElectrons=13
   Element( 31)%NumCoreStates=7
   Element( 31)%NumVariations=0
   Element( 31)%DebyeT=ZERO
   Element( 31)%AtomicRadius=1.30d0*Angstrom2Bohr
   Element( 31)%AtomicMass=69.723d0
   Element( 31)%ImplicitMuffinTinRadius=0.75d0*Element( 31)%AtomicRadius
   Element( 31)%ImplicitCoreRadius=0.75d0*Element( 31)%AtomicRadius
!
   Element( 32)%AtomName='Ge'
   Element( 32)%AtomicNumber=32
   Element( 32)%NumDeepCoreElectrons=10
   Element( 32)%NumSemiCoreElectrons=8
   Element( 32)%NumValenceElectrons=14
   Element( 32)%NumCoreStates=7
   Element( 32)%NumVariations=0
   Element( 32)%DebyeT=ZERO
   Element( 32)%AtomicRadius=1.25d0*Angstrom2Bohr
   Element( 32)%AtomicMass=72.630d0
   Element( 32)%ImplicitMuffinTinRadius=0.75d0*Element( 32)%AtomicRadius
   Element( 32)%ImplicitCoreRadius=0.75d0*Element( 32)%AtomicRadius
!
   Element( 33)%AtomName='As'
   Element( 33)%AtomicNumber=33
   Element( 33)%NumDeepCoreElectrons=10
   Element( 33)%NumSemiCoreElectrons=8
   Element( 33)%NumValenceElectrons=15
   Element( 33)%NumCoreStates=7
   Element( 33)%NumVariations=0
   Element( 33)%DebyeT=ZERO
   Element( 33)%AtomicRadius=1.15d0*Angstrom2Bohr
   Element( 33)%AtomicMass=74.922d0
   Element( 33)%ImplicitMuffinTinRadius=0.75d0*Element( 33)%AtomicRadius
   Element( 33)%ImplicitCoreRadius=0.75d0*Element( 33)%AtomicRadius
!
   Element( 34)%AtomName='Se'
   Element( 34)%AtomicNumber=34
   Element( 34)%NumDeepCoreElectrons=18
   Element( 34)%NumSemiCoreElectrons=10
   Element( 34)%NumValenceElectrons=6
   Element( 34)%NumCoreStates=9
   Element( 34)%NumVariations=0
   Element( 34)%DebyeT=ZERO
   Element( 34)%AtomicRadius=1.15d0*Angstrom2Bohr
   Element( 34)%AtomicMass=78.971d0
   Element( 34)%ImplicitMuffinTinRadius=0.75d0*Element( 34)%AtomicRadius
   Element( 34)%ImplicitCoreRadius=0.75d0*Element( 34)%AtomicRadius
!
   Element( 35)%AtomName='Br'
   Element( 35)%AtomicNumber=35
   Element( 35)%NumDeepCoreElectrons=10
   Element( 35)%NumSemiCoreElectrons=8
   Element( 35)%NumValenceElectrons=17
   Element( 35)%NumCoreStates=7
   Element( 35)%NumVariations=0
   Element( 35)%DebyeT=ZERO
   Element( 35)%AtomicRadius=1.15d0*Angstrom2Bohr
   Element( 35)%AtomicMass=79.904d0
   Element( 35)%ImplicitMuffinTinRadius=0.75d0*Element( 35)%AtomicRadius
   Element( 35)%ImplicitCoreRadius=0.75d0*Element( 35)%AtomicRadius
!
   Element( 36)%AtomName='Kr'
   Element( 36)%AtomicNumber=36
   Element( 36)%NumDeepCoreElectrons=28
   Element( 36)%NumSemiCoreElectrons=8
   Element( 36)%NumValenceElectrons=0
   Element( 36)%NumCoreStates=12
   Element( 36)%NumVariations=0
   Element( 36)%DebyeT=ZERO
   Element( 36)%AtomicRadius=1.10d0*Angstrom2Bohr
   Element( 36)%AtomicMass=83.798d0
   Element( 36)%ImplicitMuffinTinRadius=0.75d0*Element( 36)%AtomicRadius
   Element( 36)%ImplicitCoreRadius=0.75d0*Element( 36)%AtomicRadius
!
   Element( 37)%AtomName='Rb'
   Element( 37)%AtomicNumber=37
   Element( 37)%NumDeepCoreElectrons=28
   Element( 37)%NumSemiCoreElectrons=8
   Element( 37)%NumValenceElectrons=1
   Element( 37)%NumCoreStates=12
   Element( 37)%NumVariations=0
   Element( 37)%DebyeT=60.1d0
   Element( 37)%AtomicRadius=2.35d0*Angstrom2Bohr
   Element( 37)%AtomicMass=85.468d0
   Element( 37)%ImplicitMuffinTinRadius=0.75d0*Element( 37)%AtomicRadius
   Element( 37)%ImplicitCoreRadius=0.75d0*Element( 37)%AtomicRadius
!
   Element( 38)%AtomName='Sr'
   Element( 38)%AtomicNumber=38
   Element( 38)%NumDeepCoreElectrons=28
   Element( 38)%NumSemiCoreElectrons=8
   Element( 38)%NumValenceElectrons=2
   Element( 38)%NumCoreStates=12
   Element( 38)%NumVariations=0
   Element( 38)%DebyeT=93.6d0
   Element( 38)%AtomicRadius=2.00d0*Angstrom2Bohr
   Element( 38)%AtomicMass=87.62d0
   Element( 38)%ImplicitMuffinTinRadius=0.75d0*Element( 38)%AtomicRadius
   Element( 38)%ImplicitCoreRadius=0.75d0*Element( 38)%AtomicRadius
!
   Element( 39)%AtomName='Y '
   Element( 39)%AtomicNumber=39
   Element( 39)%NumDeepCoreElectrons=28
   Element( 39)%NumSemiCoreElectrons=8
   Element( 39)%NumValenceElectrons=3
   Element( 39)%NumCoreStates=12
   Element( 39)%NumVariations=0
   Element( 39)%DebyeT=ZERO
   Element( 39)%AtomicRadius=1.80d0*Angstrom2Bohr
   Element( 39)%AtomicMass=88.906d0
   Element( 39)%ImplicitMuffinTinRadius=0.75d0*Element( 39)%AtomicRadius
   Element( 39)%ImplicitCoreRadius=0.75d0*Element( 39)%AtomicRadius
!
   Element( 40)%AtomName='Zr'
   Element( 40)%AtomicNumber=40
   Element( 40)%NumDeepCoreElectrons=28
   Element( 40)%NumSemiCoreElectrons=8
   Element( 40)%NumValenceElectrons=4
   Element( 40)%NumCoreStates=12
   Element( 40)%NumVariations=0
   Element( 40)%DebyeT=225.2d0
   Element( 40)%AtomicRadius=1.55d0*Angstrom2Bohr
   Element( 40)%AtomicMass=91.224d0
   Element( 40)%ImplicitMuffinTinRadius=0.75d0*Element( 40)%AtomicRadius
   Element( 40)%ImplicitCoreRadius=0.75d0*Element( 40)%AtomicRadius
!
   Element( 41)%AtomName='Nb'
   Element( 41)%AtomicNumber=41
   Element( 41)%NumDeepCoreElectrons=28
   Element( 41)%NumSemiCoreElectrons=8
   Element( 41)%NumValenceElectrons=5
   Element( 41)%NumCoreStates=12
   Element( 41)%NumVariations=0
   Element( 41)%DebyeT=319.4d0
   Element( 41)%AtomicRadius=1.45d0*Angstrom2Bohr
   Element( 41)%AtomicMass=92.906d0
   Element( 41)%ImplicitMuffinTinRadius=0.75d0*Element( 41)%AtomicRadius
   Element( 41)%ImplicitCoreRadius=0.75d0*Element( 41)%AtomicRadius
!
   Element( 42)%AtomName='Mo'
   Element( 42)%AtomicNumber=42
   Element( 42)%NumDeepCoreElectrons=28
   Element( 42)%NumSemiCoreElectrons=8
   Element( 42)%NumValenceElectrons=6
   Element( 42)%NumCoreStates=12
   Element( 42)%NumVariations=0
   Element( 42)%DebyeT=373.3d0
   Element( 42)%AtomicRadius=1.45d0*Angstrom2Bohr
   Element( 42)%AtomicMass=95.95d0
   Element( 42)%ImplicitMuffinTinRadius=0.75d0*Element( 42)%AtomicRadius
   Element( 42)%ImplicitCoreRadius=0.75d0*Element( 42)%AtomicRadius
!
   Element( 43)%AtomName='Tc'
   Element( 43)%AtomicNumber=43
   Element( 43)%NumDeepCoreElectrons=28
   Element( 43)%NumSemiCoreElectrons=8
   Element( 43)%NumValenceElectrons=7
   Element( 43)%NumCoreStates=12
   Element( 43)%NumVariations=0
   Element( 43)%DebyeT=ZERO
   Element( 43)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 43)%AtomicMass=98.0d0
   Element( 43)%ImplicitMuffinTinRadius=0.75d0*Element( 43)%AtomicRadius
   Element( 43)%ImplicitCoreRadius=0.75d0*Element( 43)%AtomicRadius
!
   Element( 44)%AtomName='Ru'
   Element( 44)%AtomicNumber=44
   Element( 44)%NumDeepCoreElectrons=28
   Element( 44)%NumSemiCoreElectrons=8
   Element( 44)%NumValenceElectrons=8
   Element( 44)%NumCoreStates=12
   Element( 44)%NumVariations=0
   Element( 44)%DebyeT=ZERO
   Element( 44)%AtomicRadius=1.30d0*Angstrom2Bohr
   Element( 44)%AtomicMass=101.07d0
   Element( 44)%ImplicitMuffinTinRadius=0.75d0*Element( 44)%AtomicRadius
   Element( 44)%ImplicitCoreRadius=0.75d0*Element( 44)%AtomicRadius
!
   Element( 45)%AtomName='Rh'
   Element( 45)%AtomicNumber=45
   Element( 45)%NumDeepCoreElectrons=28
   Element( 45)%NumSemiCoreElectrons=8
   Element( 45)%NumValenceElectrons=9
   Element( 45)%NumCoreStates=12
   Element( 45)%NumVariations=0
   Element( 45)%DebyeT=361.8d0
   Element( 45)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 45)%AtomicMass=102.91d0
   Element( 45)%ImplicitMuffinTinRadius=0.75d0*Element( 45)%AtomicRadius
   Element( 45)%ImplicitCoreRadius=0.75d0*Element( 45)%AtomicRadius
!
   Element( 46)%AtomName='Pd'
   Element( 46)%AtomicNumber=46
   Element( 46)%NumDeepCoreElectrons=28
   Element( 46)%NumSemiCoreElectrons=8
   Element( 46)%NumValenceElectrons=10
   Element( 46)%NumCoreStates=12
   Element( 46)%NumVariations=0
   Element( 46)%DebyeT=365.8d0
   Element( 46)%AtomicRadius=1.40d0*Angstrom2Bohr
   Element( 46)%AtomicMass=106.42d0
   Element( 46)%ImplicitMuffinTinRadius=0.75d0*Element( 46)%AtomicRadius
   Element( 46)%ImplicitCoreRadius=0.75d0*Element( 46)%AtomicRadius
!
   Element( 47)%AtomName='Ag'
   Element( 47)%AtomicNumber=47
   Element( 47)%NumDeepCoreElectrons=28
   Element( 47)%NumSemiCoreElectrons=8
   Element( 47)%NumValenceElectrons=11
   Element( 47)%NumCoreStates=12
   Element( 47)%NumVariations=0
   Element( 47)%DebyeT=240.5d0
   Element( 47)%AtomicRadius=1.60d0*Angstrom2Bohr
   Element( 47)%AtomicMass=107.87d0
   Element( 47)%ImplicitMuffinTinRadius=0.75d0*Element( 47)%AtomicRadius
   Element( 47)%ImplicitCoreRadius=0.75d0*Element( 47)%AtomicRadius
!
   Element( 48)%AtomName='Cd'
   Element( 48)%AtomicNumber=48
   Element( 48)%NumDeepCoreElectrons=28
   Element( 48)%NumSemiCoreElectrons=8
   Element( 48)%NumValenceElectrons=12
   Element( 48)%NumCoreStates=12
   Element( 48)%NumVariations=0
   Element( 48)%DebyeT=ZERO
   Element( 48)%AtomicRadius=1.55d0*Angstrom2Bohr
   Element( 48)%AtomicMass=112.41d0
   Element( 48)%ImplicitMuffinTinRadius=0.75d0*Element( 48)%AtomicRadius
   Element( 48)%ImplicitCoreRadius=0.75d0*Element( 48)%AtomicRadius
!
   Element( 49)%AtomName='In'
   Element( 49)%AtomicNumber=49
   Element( 49)%NumDeepCoreElectrons=28
   Element( 49)%NumSemiCoreElectrons=18
   Element( 49)%NumValenceElectrons=3
   Element( 49)%NumCoreStates=14
   Element( 49)%NumVariations=0
   Element( 49)%DebyeT=ZERO
   Element( 49)%AtomicRadius=1.55d0*Angstrom2Bohr
   Element( 49)%AtomicMass=114.82d0
   Element( 49)%ImplicitMuffinTinRadius=0.75d0*Element( 49)%AtomicRadius
   Element( 49)%ImplicitCoreRadius=0.75d0*Element( 49)%AtomicRadius
!
   Element( 50)%AtomName='Sn'
   Element( 50)%AtomicNumber=50
   Element( 50)%NumDeepCoreElectrons=28
   Element( 50)%NumSemiCoreElectrons=18
   Element( 50)%NumValenceElectrons=4
   Element( 50)%NumCoreStates=14
   Element( 50)%NumVariations=0
   Element( 50)%DebyeT=180.5d0
   Element( 50)%AtomicRadius=1.45d0*Angstrom2Bohr
   Element( 50)%AtomicMass=118.71d0
   Element( 50)%ImplicitMuffinTinRadius=0.75d0*Element( 50)%AtomicRadius
   Element( 50)%ImplicitCoreRadius=0.75d0*Element( 50)%AtomicRadius
!
   Element( 51)%AtomName='Sb'
   Element( 51)%AtomicNumber=51
   Element( 51)%NumDeepCoreElectrons=28
   Element( 51)%NumSemiCoreElectrons=18
   Element( 51)%NumValenceElectrons=5
   Element( 51)%NumCoreStates=14
   Element( 51)%NumVariations=0
   Element( 51)%DebyeT=ZERO
   Element( 51)%AtomicRadius=1.45d0*Angstrom2Bohr
   Element( 51)%AtomicMass=121.76d0
   Element( 51)%ImplicitMuffinTinRadius=0.75d0*Element( 51)%AtomicRadius
   Element( 51)%ImplicitCoreRadius=0.75d0*Element( 51)%AtomicRadius
!
   Element( 52)%AtomName='Te'
   Element( 52)%AtomicNumber=52
   Element( 52)%NumDeepCoreElectrons=28
   Element( 52)%NumSemiCoreElectrons=18
   Element( 52)%NumValenceElectrons=6
   Element( 52)%NumCoreStates=14
   Element( 52)%NumVariations=0
   Element( 52)%DebyeT=ZERO
   Element( 52)%AtomicRadius=1.40d0*Angstrom2Bohr
   Element( 52)%AtomicMass=127.6d0
   Element( 52)%ImplicitMuffinTinRadius=0.75d0*Element( 52)%AtomicRadius
   Element( 52)%ImplicitCoreRadius=0.75d0*Element( 52)%AtomicRadius
!
   Element( 53)%AtomName='I '
   Element( 53)%AtomicNumber=53
   Element( 53)%NumDeepCoreElectrons=28
   Element( 53)%NumSemiCoreElectrons=18
   Element( 53)%NumValenceElectrons=7
   Element( 53)%NumCoreStates=14
   Element( 53)%NumVariations=0
   Element( 53)%DebyeT=ZERO
   Element( 53)%AtomicRadius=1.40d0*Angstrom2Bohr
   Element( 53)%AtomicMass=126.9d0
   Element( 53)%ImplicitMuffinTinRadius=0.75d0*Element( 53)%AtomicRadius
   Element( 53)%ImplicitCoreRadius=0.75d0*Element( 53)%AtomicRadius
!
   Element( 54)%AtomName='Xe'
   Element( 54)%AtomicNumber=54
   Element( 54)%NumDeepCoreElectrons=46
   Element( 54)%NumSemiCoreElectrons=8
   Element( 54)%NumValenceElectrons=0
   Element( 54)%NumCoreStates=17
   Element( 54)%NumVariations=0
   Element( 54)%DebyeT=ZERO
   Element( 54)%AtomicRadius=1.30d0*Angstrom2Bohr
   Element( 54)%AtomicMass=131.29d0
   Element( 54)%ImplicitMuffinTinRadius=0.75d0*Element( 54)%AtomicRadius
   Element( 54)%ImplicitCoreRadius=0.75d0*Element( 54)%AtomicRadius
!
   Element( 55)%AtomName='Cs'
   Element( 55)%AtomicNumber=55
   Element( 55)%NumDeepCoreElectrons=46
   Element( 55)%NumSemiCoreElectrons=8
   Element( 55)%NumValenceElectrons=1
   Element( 55)%NumCoreStates=17
   Element( 55)%NumVariations=0
   Element( 55)%DebyeT=ZERO
   Element( 55)%AtomicRadius=2.60d0*Angstrom2Bohr
   Element( 55)%AtomicMass=132.91d0
   Element( 55)%ImplicitMuffinTinRadius=0.75d0*Element( 55)%AtomicRadius
   Element( 55)%ImplicitCoreRadius=0.75d0*Element( 55)%AtomicRadius
!
   Element( 56)%AtomName='Ba'
   Element( 56)%AtomicNumber=56
   Element( 56)%NumDeepCoreElectrons=46
   Element( 56)%NumSemiCoreElectrons=8
   Element( 56)%NumValenceElectrons=2
   Element( 56)%NumCoreStates=17
   Element( 56)%NumVariations=0
   Element( 56)%DebyeT=60.5d0
   Element( 56)%AtomicRadius=2.15d0*Angstrom2Bohr
   Element( 56)%AtomicMass=137.33d0
   Element( 56)%ImplicitMuffinTinRadius=0.75d0*Element( 56)%AtomicRadius
   Element( 56)%ImplicitCoreRadius=0.75d0*Element( 56)%AtomicRadius
!
   Element( 57)%AtomName='La'
   Element( 57)%AtomicNumber=57
   Element( 57)%NumDeepCoreElectrons=36
   Element( 57)%NumSemiCoreElectrons=18
   Element( 57)%NumValenceElectrons=3
   Element( 57)%NumCoreStates=17
   Element( 57)%NumVariations=3
   Element( 57)%DebyeT=ZERO
   Element( 57)%AtomicRadius=1.95d0*Angstrom2Bohr
   Element( 57)%AtomicMass=138.91d0
   Element( 57)%ImplicitMuffinTinRadius=0.75d0*Element( 57)%AtomicRadius
   Element( 57)%ImplicitCoreRadius=0.75d0*Element( 57)%AtomicRadius
!
   allocate( Element( 57)%Variation(1:Element( 57)%NumVariations) )
!
   Element( 57)%Variation(1)%AtomName='La0'
   Element( 57)%Variation(1)%AtomicNumber=57
   Element( 57)%Variation(1)%NumDeepCoreElectrons=36
   Element( 57)%Variation(1)%NumSemiCoreElectrons=12
   Element( 57)%Variation(1)%NumValenceElectrons=9
   Element( 57)%Variation(1)%NumCoreStates=15
   Element( 57)%Variation(1)%NumVariations=0
!
   Element( 57)%Variation(2)%AtomName='La1'
   Element( 57)%Variation(2)%AtomicNumber=57
   Element( 57)%Variation(2)%NumDeepCoreElectrons=36
   Element( 57)%Variation(2)%NumSemiCoreElectrons=18
   Element( 57)%Variation(2)%NumValenceElectrons=3
   Element( 57)%Variation(2)%NumCoreStates=17
   Element( 57)%Variation(2)%NumVariations=0
!
   Element( 57)%Variation(3)%AtomName='La2'
   Element( 57)%Variation(3)%AtomicNumber=57
   Element( 57)%Variation(3)%NumDeepCoreElectrons=36
   Element( 57)%Variation(3)%NumSemiCoreElectrons=10
   Element( 57)%Variation(3)%NumValenceElectrons=11
   Element( 57)%Variation(3)%NumCoreStates=14
   Element( 57)%Variation(3)%NumVariations=0
!
   Element( 58)%AtomName='Ce'
   Element( 58)%AtomicNumber=58
   Element( 58)%NumDeepCoreElectrons=36
   Element( 58)%NumSemiCoreElectrons=18
   Element( 58)%NumValenceElectrons=4
   Element( 58)%NumCoreStates=17
   Element( 58)%NumVariations=3
   Element( 58)%DebyeT=ZERO
   Element( 58)%AtomicRadius=1.85d0*Angstrom2Bohr
   Element( 58)%AtomicMass=140.12d0
   Element( 58)%ImplicitMuffinTinRadius=0.75d0*Element( 58)%AtomicRadius
   Element( 58)%ImplicitCoreRadius=0.75d0*Element( 58)%AtomicRadius
!
   allocate( Element( 58)%Variation(1:Element( 58)%NumVariations) )
!
   Element( 58)%Variation(1)%AtomName='Ce0'
   Element( 58)%Variation(1)%AtomicNumber=58
   Element( 58)%Variation(1)%NumDeepCoreElectrons=36
   Element( 58)%Variation(1)%NumSemiCoreElectrons=12
   Element( 58)%Variation(1)%NumValenceElectrons=10
   Element( 58)%Variation(1)%NumCoreStates=15
   Element( 58)%Variation(1)%NumVariations=0
!
   Element( 58)%Variation(2)%AtomName='Ce1'
   Element( 58)%Variation(2)%AtomicNumber=58
   Element( 58)%Variation(2)%NumDeepCoreElectrons=36
   Element( 58)%Variation(2)%NumSemiCoreElectrons=18
   Element( 58)%Variation(2)%NumValenceElectrons=4
   Element( 58)%Variation(2)%NumCoreStates=17
   Element( 58)%Variation(2)%NumVariations=0
!
   Element( 58)%Variation(3)%AtomName='Ce2'
   Element( 58)%Variation(3)%AtomicNumber=58
   Element( 58)%Variation(3)%NumDeepCoreElectrons=36
   Element( 58)%Variation(3)%NumSemiCoreElectrons=10
   Element( 58)%Variation(3)%NumValenceElectrons=12
   Element( 58)%Variation(3)%NumCoreStates=14
   Element( 58)%Variation(3)%NumVariations=0
!
   Element( 59)%AtomName='Pr'
   Element( 59)%AtomicNumber=59
   Element( 59)%NumDeepCoreElectrons=46
   Element( 59)%NumSemiCoreElectrons=11
   Element( 59)%NumValenceElectrons=2
   Element( 59)%NumCoreStates=19
   Element( 59)%NumVariations=0
   Element( 59)%DebyeT=ZERO
   Element( 59)%AtomicRadius=1.85d0*Angstrom2Bohr
   Element( 59)%AtomicMass=140.91d0
   Element( 59)%ImplicitMuffinTinRadius=0.75d0*Element( 59)%AtomicRadius
   Element( 59)%ImplicitCoreRadius=0.75d0*Element( 59)%AtomicRadius
!
   Element( 60)%AtomName='Nd'
   Element( 60)%AtomicNumber=60
   Element( 60)%NumDeepCoreElectrons=46
   Element( 60)%NumSemiCoreElectrons=12
   Element( 60)%NumValenceElectrons=2
   Element( 60)%NumCoreStates=19
   Element( 60)%NumVariations=0
   Element( 60)%DebyeT=ZERO
   Element( 60)%AtomicRadius=1.85d0*Angstrom2Bohr
   Element( 60)%AtomicMass=144.24d0
   Element( 60)%ImplicitMuffinTinRadius=0.75d0*Element( 60)%AtomicRadius
   Element( 60)%ImplicitCoreRadius=0.75d0*Element( 60)%AtomicRadius
!
   Element( 61)%AtomName='Pm'
   Element( 61)%AtomicNumber=61
   Element( 61)%NumDeepCoreElectrons=46
   Element( 61)%NumSemiCoreElectrons=13
   Element( 61)%NumValenceElectrons=2
   Element( 61)%NumCoreStates=19
   Element( 61)%NumVariations=0
   Element( 61)%DebyeT=ZERO
   Element( 61)%AtomicRadius=1.85d0*Angstrom2Bohr
   Element( 61)%AtomicMass=145d0
   Element( 61)%ImplicitMuffinTinRadius=0.75d0*Element( 61)%AtomicRadius
   Element( 61)%ImplicitCoreRadius=0.75d0*Element( 61)%AtomicRadius
!
   Element( 62)%AtomName='Sm'
   Element( 62)%AtomicNumber=62
   Element( 62)%NumDeepCoreElectrons=36  ! 1s2s2p3s3p4s4p3d
   Element( 62)%NumSemiCoreElectrons=12  ! 4d5s
   Element( 62)%NumValenceElectrons=14
   Element( 62)%NumCoreStates=15
   Element( 62)%NumVariations=2
   Element( 62)%DebyeT=ZERO
   Element( 62)%AtomicRadius=1.85d0*Angstrom2Bohr
   Element( 62)%AtomicMass=150.36d0
   Element( 62)%ImplicitMuffinTinRadius=0.75d0*Element( 62)%AtomicRadius
   Element( 62)%ImplicitCoreRadius=0.75d0*Element( 62)%AtomicRadius
!
   allocate( Element( 62)%Variation(1:Element( 62)%NumVariations) )
!
   Element( 62)%Variation(1)%AtomName='Sm0'
   Element( 62)%Variation(1)%AtomicNumber=62
   Element( 62)%Variation(1)%NumDeepCoreElectrons=36  ! 1s2s2p3s3p4s4p3d
   Element( 62)%Variation(1)%NumSemiCoreElectrons=24  ! 4d5s5p4f(6)
   Element( 62)%Variation(1)%NumValenceElectrons=2
   Element( 62)%Variation(1)%NumCoreStates=19
   Element( 62)%Variation(1)%NumVariations=0
!
   Element( 62)%Variation(2)%AtomName='Sm1'
   Element( 62)%Variation(2)%AtomicNumber=62
   Element( 62)%Variation(2)%NumDeepCoreElectrons=46  ! 1s2s2p3s3p4s4p3d4d
   Element( 62)%Variation(2)%NumSemiCoreElectrons=14  ! 5s5p4f(6)
   Element( 62)%Variation(2)%NumValenceElectrons=2
   Element( 62)%Variation(2)%NumCoreStates=19
   Element( 62)%Variation(2)%NumVariations=0
!
   Element( 63)%AtomName='Eu'
   Element( 63)%AtomicNumber=63
   Element( 63)%NumDeepCoreElectrons=46
   Element( 63)%NumSemiCoreElectrons=15
   Element( 63)%NumValenceElectrons=2
   Element( 63)%NumCoreStates=19
   Element( 63)%NumVariations=0
   Element( 63)%DebyeT=ZERO
   Element( 63)%AtomicRadius=1.85d0*Angstrom2Bohr
   Element( 63)%AtomicMass=151.96d0
   Element( 63)%ImplicitMuffinTinRadius=0.75d0*Element( 63)%AtomicRadius
   Element( 63)%ImplicitCoreRadius=0.75d0*Element( 63)%AtomicRadius
!
   Element( 64)%AtomName='Gd'
   Element( 64)%AtomicNumber=64
   Element( 64)%NumDeepCoreElectrons=36
   Element( 64)%NumSemiCoreElectrons=10
   Element( 64)%NumValenceElectrons=18
   Element( 64)%NumCoreStates=14
   Element( 64)%NumVariations=0
   Element( 64)%DebyeT=ZERO
   Element( 64)%AtomicRadius=1.80d0*Angstrom2Bohr
   Element( 64)%AtomicMass=157.25d0
   Element( 64)%ImplicitMuffinTinRadius=0.75d0*Element( 64)%AtomicRadius
   Element( 64)%ImplicitCoreRadius=0.75d0*Element( 64)%AtomicRadius
!
   Element( 65)%AtomName='Tb'
   Element( 65)%AtomicNumber=65
   Element( 65)%NumDeepCoreElectrons=46
   Element( 65)%NumSemiCoreElectrons=17
   Element( 65)%NumValenceElectrons=2
   Element( 65)%NumCoreStates=19
   Element( 65)%NumVariations=0
   Element( 65)%DebyeT=ZERO
   Element( 65)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 65)%AtomicMass=158.93d0
   Element( 65)%ImplicitMuffinTinRadius=0.75d0*Element( 65)%AtomicRadius
   Element( 65)%ImplicitCoreRadius=0.75d0*Element( 65)%AtomicRadius
!
   Element( 66)%AtomName='Dy'
   Element( 66)%AtomicNumber=66
   Element( 66)%NumDeepCoreElectrons=46
   Element( 66)%NumSemiCoreElectrons=18
   Element( 66)%NumValenceElectrons=2
   Element( 66)%NumCoreStates=19
   Element( 66)%NumVariations=0
   Element( 66)%DebyeT=ZERO
   Element( 66)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 66)%AtomicMass=162.5d0
   Element( 66)%ImplicitMuffinTinRadius=0.75d0*Element( 66)%AtomicRadius
   Element( 66)%ImplicitCoreRadius=0.75d0*Element( 66)%AtomicRadius
!
   Element( 67)%AtomName='Ho'
   Element( 67)%AtomicNumber=67
   Element( 67)%NumDeepCoreElectrons=46
   Element( 67)%NumSemiCoreElectrons=19
   Element( 67)%NumValenceElectrons=2
   Element( 67)%NumCoreStates=19
   Element( 67)%NumVariations=0
   Element( 67)%DebyeT=ZERO
   Element( 67)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 67)%AtomicMass=164.93d0
   Element( 67)%ImplicitMuffinTinRadius=0.75d0*Element( 67)%AtomicRadius
   Element( 67)%ImplicitCoreRadius=0.75d0*Element( 67)%AtomicRadius
!
   Element( 68)%AtomName='Er'
   Element( 68)%AtomicNumber=68
   Element( 68)%NumDeepCoreElectrons=46
   Element( 68)%NumSemiCoreElectrons=20
   Element( 68)%NumValenceElectrons=2
   Element( 68)%NumCoreStates=19
   Element( 68)%NumVariations=0
   Element( 68)%DebyeT=ZERO
   Element( 68)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 68)%AtomicMass=167.26d0
   Element( 68)%ImplicitMuffinTinRadius=0.75d0*Element( 68)%AtomicRadius
   Element( 68)%ImplicitCoreRadius=0.75d0*Element( 68)%AtomicRadius
!
   Element( 69)%AtomName='Tm'
   Element( 69)%AtomicNumber=69
   Element( 69)%NumDeepCoreElectrons=46
   Element( 69)%NumSemiCoreElectrons=21
   Element( 69)%NumValenceElectrons=2
   Element( 69)%NumCoreStates=19
   Element( 69)%NumVariations=0
   Element( 69)%DebyeT=ZERO
   Element( 69)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 69)%AtomicMass=168.93d0
   Element( 69)%ImplicitMuffinTinRadius=0.75d0*Element( 69)%AtomicRadius
   Element( 69)%ImplicitCoreRadius=0.75d0*Element( 69)%AtomicRadius
!
   Element( 70)%AtomName='Yb'
   Element( 70)%AtomicNumber=70
   Element( 70)%NumDeepCoreElectrons=46
   Element( 70)%NumSemiCoreElectrons=22
   Element( 70)%NumValenceElectrons=2
   Element( 70)%NumCoreStates=19
   Element( 70)%NumVariations=0
   Element( 70)%DebyeT=ZERO
   Element( 70)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 70)%AtomicMass=173.05d0
   Element( 70)%ImplicitMuffinTinRadius=0.75d0*Element( 70)%AtomicRadius
   Element( 70)%ImplicitCoreRadius=0.75d0*Element( 70)%AtomicRadius
!
   Element( 71)%AtomName='Lu'
   Element( 71)%AtomicNumber=71
   Element( 71)%NumDeepCoreElectrons=46
   Element( 71)%NumSemiCoreElectrons=22
   Element( 71)%NumValenceElectrons=3
   Element( 71)%NumCoreStates=19
   Element( 71)%NumVariations=0
   Element( 71)%DebyeT=ZERO
   Element( 71)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 71)%AtomicMass=174.97d0
   Element( 71)%ImplicitMuffinTinRadius=0.75d0*Element( 71)%AtomicRadius
   Element( 71)%ImplicitCoreRadius=0.75d0*Element( 71)%AtomicRadius
!
   Element( 72)%AtomName='Hf'
   Element( 72)%AtomicNumber=72
   Element( 72)%NumDeepCoreElectrons=46
   Element( 72)%NumSemiCoreElectrons=22
   Element( 72)%NumValenceElectrons=4
   Element( 72)%NumCoreStates=19
   Element( 72)%NumVariations=1
   Element( 72)%DebyeT=167.1d0
   Element( 72)%AtomicRadius=1.55d0*Angstrom2Bohr
   Element( 72)%AtomicMass=178.49d0
   Element( 72)%ImplicitMuffinTinRadius=0.75d0*Element( 72)%AtomicRadius
   Element( 72)%ImplicitCoreRadius=0.75d0*Element( 72)%AtomicRadius
!
   allocate( Element( 72)%Variation(1:Element( 72)%NumVariations) )
   Element( 72)%Variation(1)%AtomName='Hf-4f'
   Element( 72)%Variation(1)%AtomicNumber=72
   Element( 72)%Variation(1)%NumDeepCoreElectrons=46
   Element( 72)%Variation(1)%NumSemiCoreElectrons=8
   Element( 72)%Variation(1)%NumValenceElectrons=18
   Element( 72)%Variation(1)%NumCoreStates=17
   Element( 72)%Variation(1)%NumVariations=0
!
   Element( 73)%AtomName='Ta'
   Element( 73)%AtomicNumber=73
   Element( 73)%NumDeepCoreElectrons=46
   Element( 73)%NumSemiCoreElectrons=22
   Element( 73)%NumValenceElectrons=5
   Element( 73)%NumCoreStates=19
   Element( 73)%NumVariations=0
   Element( 73)%DebyeT=233.2d0
   Element( 73)%AtomicRadius=1.45d0*Angstrom2Bohr
   Element( 73)%AtomicMass=180.95d0
   Element( 73)%ImplicitMuffinTinRadius=0.75d0*Element( 73)%AtomicRadius
   Element( 73)%ImplicitCoreRadius=0.75d0*Element( 73)%AtomicRadius
!
   Element( 74)%AtomName='W '
   Element( 74)%AtomicNumber=74
   Element( 74)%NumDeepCoreElectrons=46
   Element( 74)%NumSemiCoreElectrons=22
   Element( 74)%NumValenceElectrons=6
   Element( 74)%NumCoreStates=19
   Element( 74)%NumVariations=0
   Element( 74)%DebyeT=281.1d0
   Element( 74)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 74)%AtomicMass=183.84d0
   Element( 74)%ImplicitMuffinTinRadius=0.75d0*Element( 74)%AtomicRadius
   Element( 74)%ImplicitCoreRadius=0.75d0*Element( 74)%AtomicRadius
!
   Element( 75)%AtomName='Re'
   Element( 75)%AtomicNumber=75
   Element( 75)%NumDeepCoreElectrons=46
   Element( 75)%NumSemiCoreElectrons=22
   Element( 75)%NumValenceElectrons=7
   Element( 75)%NumCoreStates=19
   Element( 75)%NumVariations=0
   Element( 75)%DebyeT=ZERO
   Element( 75)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 75)%AtomicMass=186.21d0
   Element( 75)%ImplicitMuffinTinRadius=0.75d0*Element( 75)%AtomicRadius
   Element( 75)%ImplicitCoreRadius=0.75d0*Element( 75)%AtomicRadius
!
   Element( 76)%AtomName='Os'
   Element( 76)%AtomicNumber=76
   Element( 76)%NumDeepCoreElectrons=46
   Element( 76)%NumSemiCoreElectrons=22
   Element( 76)%NumValenceElectrons=8
   Element( 76)%NumCoreStates=19
   Element( 76)%NumVariations=0
   Element( 76)%DebyeT=ZERO
   Element( 76)%AtomicRadius=1.30d0*Angstrom2Bohr
   Element( 76)%AtomicMass=190.23d0
   Element( 76)%ImplicitMuffinTinRadius=0.75d0*Element( 76)%AtomicRadius
   Element( 76)%ImplicitCoreRadius=0.75d0*Element( 76)%AtomicRadius
!
   Element( 77)%AtomName='Ir'
   Element( 77)%AtomicNumber=77
   Element( 77)%NumDeepCoreElectrons=46
   Element( 77)%NumSemiCoreElectrons=22
   Element( 77)%NumValenceElectrons=9
   Element( 77)%NumCoreStates=19
   Element( 77)%NumVariations=0
   Element( 77)%DebyeT=277.1d0
   Element( 77)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 77)%AtomicMass=192.22d0
   Element( 77)%ImplicitMuffinTinRadius=0.75d0*Element( 77)%AtomicRadius
   Element( 77)%ImplicitCoreRadius=0.75d0*Element( 77)%AtomicRadius
!
   Element( 78)%AtomName='Pt'
   Element( 78)%AtomicNumber=78
   Element( 78)%NumDeepCoreElectrons=46
   Element( 78)%NumSemiCoreElectrons=22
   Element( 78)%NumValenceElectrons=10
   Element( 78)%NumCoreStates=19
   Element( 78)%NumVariations=0
   Element( 78)%DebyeT=244.5d0
   Element( 78)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 78)%AtomicMass=195.08d0
   Element( 78)%ImplicitMuffinTinRadius=0.75d0*Element( 78)%AtomicRadius
   Element( 78)%ImplicitCoreRadius=0.75d0*Element( 78)%AtomicRadius
!
   Element( 79)%AtomName='Au'
   Element( 79)%AtomicNumber=79
   Element( 79)%NumDeepCoreElectrons=46
   Element( 79)%NumSemiCoreElectrons=22
   Element( 79)%NumValenceElectrons=11
   Element( 79)%NumCoreStates=19
   Element( 79)%NumVariations=0
   Element( 79)%DebyeT=212.5d0
   Element( 79)%AtomicRadius=1.35d0*Angstrom2Bohr
   Element( 79)%AtomicMass=196.97d0
   Element( 79)%ImplicitMuffinTinRadius=0.75d0*Element( 79)%AtomicRadius
   Element( 79)%ImplicitCoreRadius=0.75d0*Element( 79)%AtomicRadius
!
   Element( 80)%AtomName='Hg'
   Element( 80)%AtomicNumber=80
   Element( 80)%NumDeepCoreElectrons=46
   Element( 80)%NumSemiCoreElectrons=22
   Element( 80)%NumValenceElectrons=12
   Element( 80)%NumCoreStates=19
   Element( 80)%NumVariations=0
   Element( 80)%DebyeT=ZERO
   Element( 80)%AtomicRadius=1.50d0*Angstrom2Bohr
   Element( 80)%AtomicMass=200.59d0
   Element( 80)%ImplicitMuffinTinRadius=0.75d0*Element( 80)%AtomicRadius
   Element( 80)%ImplicitCoreRadius=0.75d0*Element( 80)%AtomicRadius
!
   Element( 81)%AtomName='Tl'
   Element( 81)%AtomicNumber=81
   Element( 81)%NumDeepCoreElectrons=46
   Element( 81)%NumSemiCoreElectrons=32
   Element( 81)%NumValenceElectrons=3
   Element( 81)%NumCoreStates=19
   Element( 81)%NumVariations=0
   Element( 81)%DebyeT=ZERO
   Element( 81)%AtomicRadius=1.90d0*Angstrom2Bohr
   Element( 81)%AtomicMass=204.38d0
   Element( 81)%ImplicitMuffinTinRadius=0.75d0*Element( 81)%AtomicRadius
   Element( 81)%ImplicitCoreRadius=0.75d0*Element( 81)%AtomicRadius
!
   Element( 82)%AtomName='Pb'
   Element( 82)%AtomicNumber=82
   Element( 82)%NumDeepCoreElectrons=60
   Element( 82)%NumSemiCoreElectrons=18
   Element( 82)%NumValenceElectrons=4
   Element( 82)%NumCoreStates=21
   Element( 82)%NumVariations=0
   Element( 82)%DebyeT=143.0d0
   Element( 82)%AtomicRadius=1.80d0*Angstrom2Bohr
   Element( 82)%AtomicMass=207.2d0
   Element( 82)%ImplicitMuffinTinRadius=0.75d0*Element( 82)%AtomicRadius
   Element( 82)%ImplicitCoreRadius=0.75d0*Element( 82)%AtomicRadius
!
   Element( 83)%AtomName='Bi'
   Element( 83)%AtomicNumber=83
   Element( 83)%NumDeepCoreElectrons=60
   Element( 83)%NumSemiCoreElectrons=18
   Element( 83)%NumValenceElectrons=5
   Element( 83)%NumCoreStates=21
   Element( 83)%NumVariations=0
   Element( 83)%DebyeT=ZERO
   Element( 83)%AtomicRadius=1.60d0*Angstrom2Bohr
   Element( 83)%AtomicMass=208.98d0
   Element( 83)%ImplicitMuffinTinRadius=0.75d0*Element( 83)%AtomicRadius
   Element( 83)%ImplicitCoreRadius=0.75d0*Element( 83)%AtomicRadius
!
   Element( 84)%AtomName='Po'
   Element( 84)%AtomicNumber=84
   Element( 84)%NumDeepCoreElectrons=60
   Element( 84)%NumSemiCoreElectrons=18
   Element( 84)%NumValenceElectrons=6
   Element( 84)%NumCoreStates=21
   Element( 84)%NumVariations=0
   Element( 84)%DebyeT=ZERO
   Element( 84)%AtomicRadius=1.90d0*Angstrom2Bohr
   Element( 84)%AtomicMass=209d0
   Element( 84)%ImplicitMuffinTinRadius=0.75d0*Element( 84)%AtomicRadius
   Element( 84)%ImplicitCoreRadius=0.75d0*Element( 84)%AtomicRadius
!
   Element( 85)%AtomName='At'
   Element( 85)%AtomicNumber=85
   Element( 85)%NumDeepCoreElectrons=60
   Element( 85)%NumSemiCoreElectrons=18
   Element( 85)%NumValenceElectrons=7
   Element( 85)%NumCoreStates=21
   Element( 85)%NumVariations=0
   Element( 85)%DebyeT=ZERO
   Element( 85)%AtomicRadius=1.38d0*Angstrom2Bohr
   Element( 85)%AtomicMass=210d0
   Element( 85)%ImplicitMuffinTinRadius=0.75d0*Element( 85)%AtomicRadius
   Element( 85)%ImplicitCoreRadius=0.75d0*Element( 85)%AtomicRadius
!
   Element( 86)%AtomName='Rn'
   Element( 86)%AtomicNumber=86
   Element( 86)%NumDeepCoreElectrons=60
   Element( 86)%NumSemiCoreElectrons=18
   Element( 86)%NumValenceElectrons=8
   Element( 86)%NumCoreStates=21
   Element( 86)%NumVariations=0
   Element( 86)%DebyeT=ZERO
   Element( 86)%AtomicRadius=1.33d0*Angstrom2Bohr
   Element( 86)%AtomicMass=222d0
   Element( 86)%ImplicitMuffinTinRadius=0.75d0*Element( 86)%AtomicRadius
   Element( 86)%ImplicitCoreRadius=0.75d0*Element( 86)%AtomicRadius
!
   Element( 87)%AtomName='Fr'
   Element( 87)%AtomicNumber=87
   Element( 87)%NumDeepCoreElectrons=68
   Element( 87)%NumSemiCoreElectrons=10
   Element( 87)%NumValenceElectrons=9
   Element( 87)%NumCoreStates=21
   Element( 87)%NumVariations=0
   Element( 87)%DebyeT=ZERO
   Element( 87)%AtomicRadius=3.48d0*Angstrom2Bohr
   Element( 87)%AtomicMass=223d0
   Element( 87)%ImplicitMuffinTinRadius=0.75d0*Element( 87)%AtomicRadius
   Element( 87)%ImplicitCoreRadius=0.75d0*Element( 87)%AtomicRadius
!
   Element( 88)%AtomName='Ra'
   Element( 88)%AtomicNumber=88
   Element( 88)%NumDeepCoreElectrons=68
   Element( 88)%NumSemiCoreElectrons=10
   Element( 88)%NumValenceElectrons=10
   Element( 88)%NumCoreStates=21
   Element( 88)%NumVariations=0
   Element( 88)%DebyeT=ZERO
   Element( 88)%AtomicRadius=2.15d0*Angstrom2Bohr
   Element( 88)%AtomicMass=226d0
   Element( 88)%ImplicitMuffinTinRadius=0.75d0*Element( 88)%AtomicRadius
   Element( 88)%ImplicitCoreRadius=0.75d0*Element( 88)%AtomicRadius
!
   Element( 89)%AtomName='Ac'
   Element( 89)%AtomicNumber=89
   Element( 89)%NumDeepCoreElectrons=68
   Element( 89)%NumSemiCoreElectrons=10
   Element( 89)%NumValenceElectrons=11
   Element( 89)%NumCoreStates=21
   Element( 89)%NumVariations=0
   Element( 89)%DebyeT=ZERO
   Element( 89)%AtomicRadius=1.95d0*Angstrom2Bohr
   Element( 89)%AtomicMass=227d0
   Element( 89)%ImplicitMuffinTinRadius=0.75d0*Element( 89)%AtomicRadius
   Element( 89)%ImplicitCoreRadius=0.75d0*Element( 89)%AtomicRadius
!
   Element( 90)%AtomName='Th'
   Element( 90)%AtomicNumber=90
   Element( 90)%NumDeepCoreElectrons=68
   Element( 90)%NumSemiCoreElectrons=10
   Element( 90)%NumValenceElectrons=12
   Element( 90)%NumCoreStates=21
   Element( 90)%NumVariations=0
   Element( 90)%DebyeT=ZERO
   Element( 90)%AtomicRadius=1.80d0*Angstrom2Bohr
   Element( 90)%AtomicMass=232.04d0
   Element( 90)%ImplicitMuffinTinRadius=0.75d0*Element( 90)%AtomicRadius
   Element( 90)%ImplicitCoreRadius=0.75d0*Element( 90)%AtomicRadius
!
   Element( 91)%AtomName='Pa'
   Element( 91)%AtomicNumber=91
   Element( 91)%NumDeepCoreElectrons=68
   Element( 91)%NumSemiCoreElectrons=10
   Element( 91)%NumValenceElectrons=13
   Element( 91)%NumCoreStates=21
   Element( 91)%NumVariations=0
   Element( 91)%DebyeT=ZERO
   Element( 91)%AtomicRadius=1.80d0*Angstrom2Bohr
   Element( 91)%AtomicMass=231.04d0
   Element( 91)%ImplicitMuffinTinRadius=0.75d0*Element( 91)%AtomicRadius
   Element( 91)%ImplicitCoreRadius=0.75d0*Element( 91)%AtomicRadius
!
   Element( 92)%AtomName='U'
   Element( 92)%AtomicNumber=92
   Element( 92)%NumDeepCoreElectrons=68
   Element( 92)%NumSemiCoreElectrons=10
   Element( 92)%NumValenceElectrons=14
   Element( 92)%NumCoreStates=21
   Element( 92)%NumVariations=0
   Element( 92)%DebyeT=ZERO
   Element( 92)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 92)%AtomicMass=238.03d0
   Element( 92)%ImplicitMuffinTinRadius=0.75d0*Element( 92)%AtomicRadius
   Element( 92)%ImplicitCoreRadius=0.75d0*Element( 92)%AtomicRadius
!
   Element( 93)%AtomName='Np'
   Element( 93)%AtomicNumber=93
   Element( 93)%NumDeepCoreElectrons=68
   Element( 93)%NumSemiCoreElectrons=10
   Element( 93)%NumValenceElectrons=15
   Element( 93)%NumCoreStates=21
   Element( 93)%NumVariations=0
   Element( 93)%DebyeT=ZERO
   Element( 93)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 93)%AtomicMass=237d0
   Element( 93)%ImplicitMuffinTinRadius=0.75d0*Element( 93)%AtomicRadius
   Element( 93)%ImplicitCoreRadius=0.75d0*Element( 93)%AtomicRadius
!
   Element( 94)%AtomName='Pu'
   Element( 94)%AtomicNumber=94
   Element( 94)%NumDeepCoreElectrons=68
   Element( 94)%NumSemiCoreElectrons=10
   Element( 94)%NumValenceElectrons=16
   Element( 94)%NumCoreStates=21
   Element( 94)%NumVariations=0
   Element( 94)%DebyeT=ZERO
   Element( 94)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 94)%AtomicMass=244d0
   Element( 94)%ImplicitMuffinTinRadius=0.75d0*Element( 94)%AtomicRadius
   Element( 94)%ImplicitCoreRadius=0.75d0*Element( 94)%AtomicRadius
!
   Element( 95)%AtomName='Am'
   Element( 95)%AtomicNumber=95
   Element( 95)%NumDeepCoreElectrons=68
   Element( 95)%NumSemiCoreElectrons=10
   Element( 95)%NumValenceElectrons=17
   Element( 95)%NumCoreStates=21
   Element( 95)%NumVariations=0
   Element( 95)%DebyeT=ZERO
   Element( 95)%AtomicRadius=1.75d0*Angstrom2Bohr
   Element( 95)%AtomicMass=243d0
   Element( 95)%ImplicitMuffinTinRadius=0.75d0*Element( 95)%AtomicRadius
   Element( 95)%ImplicitCoreRadius=0.75d0*Element( 95)%AtomicRadius
!
   Element( 96)%AtomName='Cm'
   Element( 96)%AtomicNumber=96
   Element( 96)%NumDeepCoreElectrons=68
   Element( 96)%NumSemiCoreElectrons=10
   Element( 96)%NumValenceElectrons=18
   Element( 96)%NumCoreStates=21
   Element( 96)%NumVariations=0
   Element( 96)%DebyeT=ZERO
   Element( 96)%AtomicRadius=1.74d0*Angstrom2Bohr
   Element( 96)%AtomicMass=247d0
   Element( 96)%ImplicitMuffinTinRadius=0.75d0*Element( 96)%AtomicRadius
   Element( 96)%ImplicitCoreRadius=0.75d0*Element( 96)%AtomicRadius
!
   Element( 97)%AtomName='Bk'
   Element( 97)%AtomicNumber=97
   Element( 97)%NumDeepCoreElectrons=68
   Element( 97)%NumSemiCoreElectrons=10
   Element( 97)%NumValenceElectrons=19
   Element( 97)%NumCoreStates=21
   Element( 97)%NumVariations=0
   Element( 97)%DebyeT=ZERO
   Element( 97)%AtomicRadius=1.70d0*Angstrom2Bohr
   Element( 97)%AtomicMass=247d0
   Element( 97)%ImplicitMuffinTinRadius=0.75d0*Element( 97)%AtomicRadius
   Element( 97)%ImplicitCoreRadius=0.75d0*Element( 97)%AtomicRadius
!
   Element( 98)%AtomName='Cf'
   Element( 98)%AtomicNumber=98
   Element( 98)%NumDeepCoreElectrons=68
   Element( 98)%NumSemiCoreElectrons=10
   Element( 98)%NumValenceElectrons=20
   Element( 98)%NumCoreStates=21
   Element( 98)%NumVariations=0
   Element( 98)%DebyeT=ZERO
   Element( 98)%AtomicRadius=1.86d0*Angstrom2Bohr
   Element( 98)%AtomicMass=251d0
   Element( 98)%ImplicitMuffinTinRadius=0.75d0*Element( 98)%AtomicRadius
   Element( 98)%ImplicitCoreRadius=0.75d0*Element( 98)%AtomicRadius
!
   Element( 99)%AtomName='Es'
   Element( 99)%AtomicNumber=99
   Element( 99)%NumDeepCoreElectrons=68
   Element( 99)%NumSemiCoreElectrons=10
   Element( 99)%NumValenceElectrons=21
   Element( 99)%NumCoreStates=21
   Element( 99)%NumVariations=0
   Element( 99)%DebyeT=ZERO
   Element( 99)%AtomicRadius=1.86d0*Angstrom2Bohr
   Element( 99)%AtomicMass=252d0
   Element( 99)%ImplicitMuffinTinRadius=0.75d0*Element( 99)%AtomicRadius
   Element( 99)%ImplicitCoreRadius=0.75d0*Element( 99)%AtomicRadius
!
   Element(100)%AtomName='Fm'
   Element(100)%AtomicNumber=100
   Element(100)%NumDeepCoreElectrons=68
   Element(100)%NumSemiCoreElectrons=10
   Element(100)%NumValenceElectrons=22
   Element(100)%NumCoreStates=21
   Element(100)%NumVariations=0
   Element(100)%DebyeT=ZERO
   Element(100)%AtomicRadius=ZERO
   Element(100)%AtomicMass=257d0
   Element(100)%ImplicitMuffinTinRadius=0.75d0*Element(100)%AtomicRadius
   Element(100)%ImplicitCoreRadius=0.75d0*Element(100)%AtomicRadius
!
   Element(101)%AtomName='Md'
   Element(101)%AtomicNumber=101
   Element(101)%NumDeepCoreElectrons=68
   Element(101)%NumSemiCoreElectrons=10
   Element(101)%NumValenceElectrons=22
   Element(101)%NumCoreStates=21
   Element(101)%NumVariations=0
   Element(101)%DebyeT=ZERO
   Element(101)%AtomicRadius=ZERO
   Element(101)%AtomicMass=258d0
   Element(101)%ImplicitMuffinTinRadius=0.75d0*Element(101)%AtomicRadius
   Element(101)%ImplicitCoreRadius=0.75d0*Element(101)%AtomicRadius
!
   Element(102)%AtomName='No'
   Element(102)%AtomicNumber=102
   Element(102)%NumDeepCoreElectrons=68
   Element(102)%NumSemiCoreElectrons=10
   Element(102)%NumValenceElectrons=23
   Element(102)%NumCoreStates=21
   Element(102)%NumVariations=0
   Element(102)%DebyeT=ZERO
   Element(102)%AtomicRadius=ZERO
   Element(102)%AtomicMass=259d0
   Element(102)%ImplicitMuffinTinRadius=0.75d0*Element(102)%AtomicRadius
   Element(102)%ImplicitCoreRadius=0.75d0*Element(102)%AtomicRadius
!
   Element(103)%AtomName='Lr'
   Element(103)%AtomicNumber=103
   Element(103)%NumDeepCoreElectrons=68
   Element(103)%NumSemiCoreElectrons=10
   Element(103)%NumValenceElectrons=23
   Element(103)%NumCoreStates=21
   Element(103)%NumVariations=0
   Element(103)%DebyeT=ZERO
   Element(103)%AtomicRadius=ZERO
   Element(103)%AtomicMass=266d0
   Element(103)%ImplicitMuffinTinRadius=0.75d0*Element(103)%AtomicRadius
   Element(103)%ImplicitCoreRadius=0.75d0*Element(103)%AtomicRadius
!
   Initialized=.true.
!
   end subroutine initChemElement
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getAtomInfo(AtomName,ztot,zcor,zsem,zval)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   character (len=*), intent(in) :: AtomName
!
   integer (kind=IntKind), intent(out) :: ztot,zcor,zsem,zval
   integer (kind=IntKind) :: j, z
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getAtomInfo','Invalid Atom Name',AtomName)
   endif
!
   if (AtomName(1:2).eq.Element(z)%AtomName(1:2)) then
      if (AtomName == Element(z)%AtomName .or.                        &
          Element(z)%NumVariations == 0) then
         ztot=Element(z)%AtomicNumber
         zcor=Element(z)%NumDeepCoreElectrons
         zsem=Element(z)%NumSemiCoreElectrons
         zval=Element(z)%NumValenceElectrons
         return
      else
         do j = 1, Element(z)%NumVariations
            if (AtomName == Element(z)%Variation(j)%AtomName) then
               ztot=Element(z)%Variation(j)%AtomicNumber
               zcor=Element(z)%Variation(j)%NumDeepCoreElectrons
               zsem=Element(z)%Variation(j)%NumSemiCoreElectrons
               zval=Element(z)%Variation(j)%NumValenceElectrons
               return
            endif
         enddo
      endif
   endif
!  -------------------------------------------------------------------
   call ErrorHandler('getAtomInfo',                                    &
                     'Atom name is not found in the periodic table',AtomName)
!  -------------------------------------------------------------------
!
   end subroutine getAtomInfo
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getName(AtomNumber) result(name)
!  ===================================================================
   implicit none
!
   character (len=2) :: name
   character (len=7), parameter :: sname='getName'
!
   integer, intent(in) :: AtomNumber
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getName','Invalid Atom Number',AtomNumber)
   endif
!
   name=Element(AtomNumber)%AtomName(1:2)
   end function getName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isAtomName(AtomName) result(y)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
!
   logical :: y
!
   character (len=2) :: a
   character (len=7), parameter :: sname='isAtomName'
   integer (kind=IntKind) :: i, k
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   k=len_trim(adjustl(AtomName))
   if (k > MaxLenOfAtomName .or. k < 1) then
      y = .false.
      return
   endif
!
   a=trim(adjustl(AtomName))
!  if(a(1:1) < 'A' .and. a(1:1) /= '_') then
   if(a(1:1) < 'A') then
      a(1:1)=achar(iachar(a(1:1))+iachar('A')-iachar('a'))
   endif
!  if (k == 2 .and. a(2:2) >= 'A' .and. a(2:2) <= 'Z' .and. a(1:1) /= '_') then
   if (k == 2 .and. a(2:2) >= 'A' .and. a(2:2) <= 'Z') then
      a(2:2)=achar(iachar(a(2:2))-iachar('A')+iachar('a'))
   endif
!
   y = .false.
   do i=MinZtot,NumElements
      if (Element(i)%AtomName(1:2) == a(1:2)) then
         y = .true.
         exit
      endif
   enddo
!
   end function isAtomName
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getZtot(AtomName) result(z)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   character (len=2) :: a
   character (len=7), parameter :: sname='getZtot'
   integer (kind=IntKind) :: i, k, z
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   k=len_trim(adjustl(AtomName))
   if (k > MaxLenOfAtomName .or. k < 1) then
      call ErrorHandler('getZtot','Invalid Atom Name',AtomName)
   endif
!
   a=trim(adjustl(AtomName))
!  if(a(1:1) < 'A' .and. a(1:1) /= '_') then
   if(a(1:1) < 'A') then
      a(1:1)=achar(iachar(a(1:1))+iachar('A')-iachar('a'))
   endif
!  if (k == 2 .and. a(2:2) >= 'A' .and. a(2:2) <= 'Z' .and. a(1:1) /= '_') then
   if (k == 2 .and. a(2:2) >= 'A' .and. a(2:2) <= 'Z') then
      a(2:2)=achar(iachar(a(2:2))-iachar('A')+iachar('a'))
   endif
!
   z=MinZtot-1
   do i=MinZtot,NumElements
      if (Element(i)%AtomName(1:2) == a(1:2)) then
         z=Element(i)%AtomicNumber
         exit
      endif
   enddo
   if (z < MinZtot) then
      call ErrorHandler(sname,'Invalid atom name',AtomName)
   endif
!
   end function getZtot
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getZcor_a(AtomName) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   logical :: found
   integer (kind=IntKind) :: z, c, i
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getZcor','Invalid Atom Name',AtomName)
   endif
!
   found = .false.
   if (Element(z)%NumVariations == 0 .or. Element(z)%AtomName == AtomName) then
      c=Element(z)%NumDeepCoreElectrons
      found = .true.
   else
      LOOP_i: do i = 1, Element(z)%NumVariations
         if (Element(z)%Variation(i)%AtomName == AtomName) then
            c=Element(z)%Variation(i)%NumDeepCoreElectrons
            found = .true.
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   if (.not. found) then
      call ErrorHandler('getZcor','Invalid Atom Name',AtomName)
   endif
!
   end function getZcor_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getZcor_n(AtomNumber) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getZcor','Invalid Atom Number',AtomNumber)
   endif
!
   c=Element(AtomNumber)%NumDeepCoreElectrons
   end function getZcor_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getZsem_a(AtomName) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   logical :: found
   integer (kind=IntKind) :: z, c, i
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getZsem','Invalid Atom Name',AtomName)
   endif
!
   found = .false.
   if (Element(z)%NumVariations == 0 .or. Element(z)%AtomName == AtomName) then
      c=Element(z)%NumSemiCoreElectrons
      found = .true.
   else
      LOOP_i: do i = 1, Element(z)%NumVariations
         if (Element(z)%Variation(i)%AtomName == AtomName) then
            c=Element(z)%Variation(i)%NumSemiCoreElectrons
            found = .true.
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   if (.not. found) then
      call ErrorHandler('getZsem','Invalid Atom Name',AtomName)
   endif
!
   end function getZsem_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getZsem_n(AtomNumber) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getZsem','Invalid Atom Number',AtomNumber)
   endif
!
   c=Element(AtomNumber)%NumSemiCoreElectrons
   end function getZsem_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getZval_a(AtomName) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   logical :: found
   integer (kind=IntKind) :: z, c, i
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getZval','Invalid Atom Name',AtomName)
   endif
!
   found = .false.
   if (Element(z)%NumVariations == 0 .or. Element(z)%AtomName == AtomName) then
      c=Element(z)%NumValenceElectrons
      found = .true.
   else
      LOOP_i: do i = 1, Element(z)%NumVariations
         if (Element(z)%Variation(i)%AtomName == AtomName) then
            c=Element(z)%Variation(i)%NumValenceElectrons
            found = .true.
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   if (.not. found) then
      call ErrorHandler('getZval','Invalid Atom Name',AtomName)
   endif
!
   end function getZval_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getZval_n(AtomNumber) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getZval','Invalid Atom Number',AtomNumber)
   endif
!
   c=Element(AtomNumber)%NumValenceElectrons
   end function getZval_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCoreStates_a(AtomName) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   logical :: found
   integer (kind=IntKind) :: z, c, i
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getNumCoreStates','Invalid Atom Name',AtomName)
   endif
!
   found = .false.
   if (Element(z)%NumVariations == 0 .or. Element(z)%AtomName == AtomName) then
      c=Element(z)%NumCoreStates
      found = .true.
   else
      LOOP_i: do i = 1, Element(z)%NumVariations
         if (Element(z)%Variation(i)%AtomName == AtomName) then
            c=Element(z)%Variation(i)%NumCoreStates
            found = .true.
            exit LOOP_i
         endif
      enddo LOOP_i
   endif
!
   if (.not. found) then
      call ErrorHandler('getNumCoreStates','Invalid Atom Name',AtomName)
   endif
!
   end function getNumCoreStates_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumCoreStates_n(AtomNumber) result(c)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getNumCoreStates','Invalid Atom Number',AtomNumber)
   endif
!
   c=Element(AtomNumber)%NumCoreStates
   end function getNumCoreStates_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateIndex(n,l,k) result(i)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n, l, k
   integer (kind=IntKind) :: i, j
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   do j = 1, MaxNumc
      if (n == CoreState(j)%n .and. l == CoreState(j)%l .and.         &
          k == CoreState(j)%kappa) then
         i = j
         return
      endif
   enddo
!
   call ErrorHandler('getCoreStateIndex','Invalid core state config',n,l,k)
!
   end function getCoreStateIndex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateN_a(AtomName,i) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: z, c, j
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getCoreStateN','Invalid Atom Name',AtomName)
   else if (i < 1 .or. i > MaxNumc) then
      call ErrorHandler('getCoreStateN','Invalid Core State Index',i)
   else if (Element(z)%NumVariations == 0 ) then 
!           .and. i <= Element(z)%NumCoreStates) then
      c=CoreState(i)%n
      return
   else if (Element(z)%NumVariations > 0) then
      if (AtomName == Element(z)%AtomName) then
!         .and. i <= Element(z)%NumCoreStates) then
         c=CoreState(i)%n
         return
      else
         do j = 1, Element(z)%NumVariations
            if (AtomName == Element(z)%Variation(j)%AtomName) then 
!               .and. i <= Element(z)%Variation(j)%NumCoreStates) then
               c=CoreState(i)%n
               return
            endif
         enddo
      endif
   endif
!
   call ErrorHandler('getCoreStateN','Invalid atom name',AtomName)
!
   end function getCoreStateN_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateN_n(AtomNumber, i) result(c)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber, i
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getCoreStateN','Invalid Atom Number',AtomNumber)
   else if (i < 1 .or. i > Element(AtomNumber)%NumCoreStates) then
      call ErrorHandler('getCoreStateN','Invalid Core State Index',i)
   endif
!
   c=CoreState(i)%n
   end function getCoreStateN_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateL_a(AtomName,i) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: z, c, j
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getCoreStateL','Invalid Atom Name',AtomName)
   else if (i < 1 .or. i > MaxNumc) then
      call ErrorHandler('getCoreStateL','Invalid Core State Index',i)
   else if (Element(z)%NumVariations == 0 ) then
!           .and. i <= Element(z)%NumCoreStates) then
      c=CoreState(i)%l
      return
   else if (Element(z)%NumVariations > 0) then
      if (AtomName == Element(z)%AtomName) then
!         .and. i <= Element(z)%NumCoreStates) then
         c=CoreState(i)%l
         return
      else
         do j = 1, Element(z)%NumVariations
            if (AtomName == Element(z)%Variation(j)%AtomName) then
!               .and. i <= Element(z)%Variation(j)%NumCoreStates) then
               c=CoreState(i)%l
               return
            endif
         enddo
      endif
   endif
!
   call ErrorHandler('getCoreStateL','Invalid atom name',AtomName)
!
   end function getCoreStateL_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateL_n(AtomNumber,i) result(c)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber, i
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getCoreStateL','Invalid Atom Number',AtomNumber)
   else if (i < 1 .or. i > Element(AtomNumber)%NumCoreStates) then
      call ErrorHandler('getCoreStateL','Invalid Core State Index',i)
   endif
!
   c=CoreState(i)%l
   end function getCoreStateL_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateKappa_a(AtomName,i) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: z, c, j
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getCoreStateKappa','Invalid Atom Name',AtomName)
   else if (i < 1 .or. i > MaxNumc) then
      call ErrorHandler('getCoreStateKappa','Invalid Core State Index',i)
   else if (Element(z)%NumVariations == 0) then
!           .and. i <= Element(z)%NumCoreStates) then
      c=CoreState(i)%kappa
      return
   else if (Element(z)%NumVariations > 0) then
      if (AtomName == Element(z)%AtomName) then
!         .and. i <= Element(z)%NumCoreStates) then
         c=CoreState(i)%kappa
         return
      else
         do j = 1, Element(z)%NumVariations
            if (AtomName == Element(z)%Variation(j)%AtomName) then
!               .and. i <= Element(z)%Variation(j)%NumCoreStates) then
               c=CoreState(i)%kappa
               return
            endif
         enddo
      endif
   endif
!
   call ErrorHandler('getCoreStateKappa','Invalid atom name',AtomName)
!
   end function getCoreStateKappa_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateKappa_n(AtomNumber,i) result(c)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber,i
   integer (kind=IntKind) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getCoreStateKappa','Invalid Atom Number',AtomNumber)
   else if (i < 1 .or. i > Element(AtomNumber)%NumCoreStates) then
      call ErrorHandler('getCoreStateKappa','Invalid Core State Index',i)
   endif
!
   c=CoreState(i)%kappa
   end function getCoreStateKappa_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateSymbol_a(AtomName,i) result(c)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   character (len=2) :: c
   integer (kind=IntKind), intent(in) :: i
   integer (kind=IntKind) :: z, j
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('getCoreStateSymbol','Invalid Atom Name',AtomName)
   else if (i < 1 .or. i > MaxNumc) then
      call ErrorHandler('getCoreStateSymbol','Invalid Core State Index',i)
   else if (Element(z)%NumVariations == 0) then
!           .and. i <= Element(z)%NumCoreStates) then
      c=CoreState(i)%Symbol
      return
   else if (Element(z)%NumVariations > 0) then
      if (AtomName == Element(z)%AtomName) then
!         .and. i <= Element(z)%NumCoreStates) then
         c=CoreState(i)%Symbol
         return
      else
         do j = 1, Element(z)%NumVariations
            if (AtomName == Element(z)%Variation(j)%AtomName) then
!               .and. i <= Element(z)%Variation(j)%NumCoreStates) then
               c=CoreState(i)%Symbol
               return
            endif
         enddo
      endif
   endif
!
   call ErrorHandler('getCoreStateSymbol','Invalid atom name',AtomName)
!
   end function getCoreStateSymbol_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCoreStateSymbol_n(AtomNumber,i) result(c)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber, i
   character (len=2) :: c
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getCoreStateSymbol','Invalid Atom Number',AtomNumber)
   else if (i < 1 .or. i > Element(AtomNumber)%NumCoreStates) then
      call ErrorHandler('getCoreStateSymbol','Invalid Core State Index',i)
   endif
!
   c=CoreState(i)%Symbol
   end function getCoreStateSymbol_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDebyeTemperature_n(AtomNumber) result(t)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   real (kind=RealKind) :: t
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getDebyeTemperature','Invalid Atom Number',AtomNumber)
   endif
!
   t=Element(AtomNumber)%DebyeT
!
   end function getDebyeTemperature_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDebyeTemperature_a(AtomName) result(t)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind) :: z
   real (kind=RealKind) :: t
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot .or. z > NumElements) then
      call ErrorHandler('getDebyeTemperature','Invalid Atom Name',AtomName)
   endif
!
   t=Element(z)%DebyeT
!
   end function getDebyeTemperature_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicRadius_n(AtomNumber) result(r)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getAtomicRadius','Invalid Atom Number',AtomNumber)
   endif
!
   r=Element(AtomNumber)%AtomicRadius
!
   end function getAtomicRadius_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicRadius_a(AtomName) result(r)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind) :: z
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot .or. z > NumElements) then
      call ErrorHandler('getAtomicRadius','Invalid Atom Number',AtomName)
   endif
!
   r=Element(z)%AtomicRadius
!
   end function getAtomicRadius_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicMass_n(AtomNumber) result(r)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getAtomicMass','Invalid Atom Number',AtomNumber)
   endif
!
   r=Element(AtomNumber)%AtomicMass
!
   end function getAtomicMass_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAtomicMass_a(AtomName) result(r)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind) :: z
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot .or. z > NumElements) then
      call ErrorHandler('getAtomicMass','Invalid Atom Number',AtomName)
   endif
!
   r=Element(z)%AtomicMass
!
   end function getAtomicMass_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getImplicitMuffinTinRadius_n(AtomNumber) result(r)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getImplicitMuffinTinRadius','Invalid Atom Number',AtomNumber)
   endif
!
   r=Element(AtomNumber)%ImplicitMuffinTinRadius
!
   end function getImplicitMuffinTinRadius_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getImplicitMuffinTinRadius_a(AtomName) result(r)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind) :: z
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot .or. z > NumElements) then
      call ErrorHandler('getImplicitMuffinTinRadius','Invalid Atom Number',AtomName)
   endif
!
   r=Element(z)%ImplicitMuffinTinRadius
!
   end function getImplicitMuffinTinRadius_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getImplicitCoreRadius_n(AtomNumber) result(r)
!  Note: This function is not implemented for other variations of the element
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: AtomNumber
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('getImplicitCoreRadius','Invalid Atom Number',AtomNumber)
   endif
!
   r=Element(AtomNumber)%ImplicitCoreRadius
!
   end function getImplicitCoreRadius_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getImplicitCoreRadius_a(AtomName) result(r)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind) :: z
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot .or. z > NumElements) then
      call ErrorHandler('getImplicitCoreRadius','Invalid Atom Number',AtomName)
   endif
!
   r=Element(z)%ImplicitCoreRadius
!
   end function getImplicitCoreRadius_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setConfig_a(AtomName,Zc,Zs)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
   integer (kind=IntKind), intent(in), optional :: Zc, Zs
   integer (kind=IntKind) :: z
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('setConfiguration','Invalid Atom Name',AtomName)
   else if (.not.present(Zc) .and. .not.present(Zs)) then
      call WarningHandler('setConfiguration','New configuration is not provided')
      return
   endif
!
   if (present(Zc)) then
      Element(z)%NumDeepCoreElectrons = Zc
   endif
!
   if (present(Zs)) then
      Element(z)%NumSemiCoreElectrons = Zs
   endif
!
   Element(z)%NumValenceElectrons = z - Element(z)%NumSemiCoreElectrons - &
                                        Element(z)%NumDeepCoreElectrons
!
!  -------------------------------------------------------------------
!  call setNumCoreStates(z)
!  -------------------------------------------------------------------
!
   end subroutine setConfig_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setConfig_n(AtomNumber,Zc,Zs)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: AtomNumber
   integer (kind=IntKind), intent(in), optional :: Zc, Zs
   integer (kind=IntKind) :: z
!
   if (.not.Initialized) then
      call initChemElement()
   endif
!
   if (AtomNumber < MinZtot .or. AtomNumber > NumElements) then
      call ErrorHandler('setConfiguration','Invalid Atom Number',AtomNumber)
   else if (.not.present(Zc) .and. .not.present(Zs)) then
      call WarningHandler('setConfiguration','New configuration is not provided')
      return
   endif
!
   z = AtomNumber
   if (present(Zc)) then
      Element(z)%NumDeepCoreElectrons = Zc
   endif
!
   if (present(Zs)) then
      Element(z)%NumSemiCoreElectrons = Zs
   endif
!
   Element(z)%NumValenceElectrons = z - Element(z)%NumSemiCoreElectrons - &
                                        Element(z)%NumDeepCoreElectrons
!
!  -------------------------------------------------------------------
!  call setNumCoreStates(z)
!  -------------------------------------------------------------------
!
   end subroutine setConfig_n
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setNumCoreStates_a(AtomName,numc)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: AtomName
!
   integer (kind=IntKind), intent(in), optional :: numc
   integer (kind=IntKind) :: z, ne, nc, zct
!
   z=getZtot(AtomName)
   if (z < MinZtot) then
      call ErrorHandler('setNumCoreStates','Invalid Atom Name',AtomName)
   endif
!
   if (present(numc)) then
      Element(z)%NumCoreStates = numc
      ne = 0
      do nc = 1, numc
         if (CoreState(nc)%l == 0) then
            ne = ne + 2
         else
            ne = ne + 2*CoreState(nc)%l+1
         endif
      enddo     
      if (ne < Element(z)%NumDeepCoreElectrons) then
         Element(z)%NumDeepCoreElectrons = ne
         Element(z)%NumSemiCoreElectrons = 0
      else
         Element(z)%NumSemiCoreElectrons = ne - Element(z)%NumDeepCoreElectrons
      endif
      Element(z)%NumValenceElectrons = z - Element(z)%NumSemiCoreElectrons - &
                                           Element(z)%NumDeepCoreElectrons
   else  
      zct = Element(z)%NumDeepCoreElectrons + Element(z)%NumSemiCoreElectrons
      nc = 0
      ne = 0
      do while (ne < zct .and. nc < 26)
         nc = nc + 1
         if (CoreState(nc)%l == 0) then
            ne = ne + 2
         else
            ne = ne + 2*CoreState(nc)%l+1
         endif
      enddo
!
      if (ne /= zct) then
!        -------------------------------------------------------------
         call ErrorHandler('setConfiguration',                           &
                           'Inconsistency between the number of core electrons and states', &
                           zct,ne,nc)
!        -------------------------------------------------------------
      endif
      Element(z)%NumCoreStates = nc
   endif
!
   end subroutine setNumCoreStates_a
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setNumCoreStates_n(z,numc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: z
   integer (kind=IntKind), intent(in), optional :: numc
   integer (kind=IntKind) :: ne, nc, zct
!
   if (present(numc)) then
      Element(z)%NumCoreStates = numc
      ne = 0
      do nc = 1, numc
         if (CoreState(nc)%l == 0) then
            ne = ne + 2
         else
            ne = ne + 2*CoreState(nc)%l+1
         endif
      enddo
      if (ne < Element(z)%NumDeepCoreElectrons) then
         Element(z)%NumDeepCoreElectrons = ne
         Element(z)%NumSemiCoreElectrons = 0
      else
         Element(z)%NumSemiCoreElectrons = ne - Element(z)%NumDeepCoreElectrons
      endif
      Element(z)%NumValenceElectrons = z - Element(z)%NumSemiCoreElectrons - &
                                           Element(z)%NumDeepCoreElectrons
   else
      zct = Element(z)%NumDeepCoreElectrons + Element(z)%NumSemiCoreElectrons
      nc = 0
      ne = 0
      do while (ne < zct .and. nc < 26)
         nc = nc + 1
         if (CoreState(nc)%l == 0) then
            ne = ne + 2
         else
            ne = ne + 2*CoreState(nc)%l+1
         endif
      enddo
!
      if (ne /= zct) then
!        -------------------------------------------------------------
         call ErrorHandler('setConfiguration',                           &
                           'Inconsistency between the number of core electrons and states', &
                           zct,ne,nc)
!        -------------------------------------------------------------
      endif
      Element(z)%NumCoreStates = nc
   endif
!
   end subroutine setNumCoreStates_n
!  ===================================================================
!
end module ChemElementModule
