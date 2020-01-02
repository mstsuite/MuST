module KindParamModule
!
! determine kind parameters for integer, real and complex types
!
  implicit none
!
  integer, parameter :: IntKind = kind(1)
  integer, parameter :: LongIntKind = 2*IntKind
  integer, parameter :: RealKind = kind(1.0d0) 
  integer, parameter :: CmplxKind = RealKind
  integer, parameter :: QuadRealKind = 2*RealKind
  integer, parameter :: QuadCmplxKind = 2*CmplxKind
!
end module KindParamModule
