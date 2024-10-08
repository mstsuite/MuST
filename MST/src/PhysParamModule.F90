module PhysParamModule
!
!  setup oftenly used physics parameters in atomic units
!
   use KindParamModule, only : RealKind
!
   implicit none
!
   private :: RealKind
!
   real (kind=RealKind), parameter :: FineStruct   = 1.0d0/137.0359895d0
   real (kind=RealKind), parameter :: LightSpeed   = 274.0719790d0
!
!  Boltzmann = 1.38062d-16/1.60219d-12/13.6058 Ryd/K
   real (kind=RealKind), parameter :: Boltzmann    = 6.3333870d-06
!
!  BohrMagneton = 9.27410d-21/1.60219d-12/13.6058 Ryd/Gauss
   real (kind=RealKind), parameter :: BohrMagneton = 4.2543545d-10
   real (kind=RealKind), parameter :: Ryd2eV = 13.6054000d0
!
   real (kind=RealKind), parameter :: Kelvin2eV = 8.61732814974493D-05
   real (kind=RealKind), parameter :: Kelvin2Ryd = 6.33361706838587E-06
!
   real (kind=RealKind), parameter :: MassUnit2eV = 931.4941D+06
   real (kind=RealKind), parameter :: MassUnit2Ryd = MassUnit2eV/Ryd2eV
!
   real (kind=RealKind), parameter :: Bohr2Angstrom = 0.529177210903d0
   real (kind=RealKind), parameter :: Angstrom2Bohr = 1.0d0/Bohr2Angstrom
!
end module PhysParamModule
