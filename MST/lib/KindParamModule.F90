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
#ifdef NVHPC
  integer, parameter :: QuadRealKind = kind(1.0)
  integer, parameter :: QuadCmplxKind = kind(1.0)
#else
  integer, parameter :: QuadRealKind = 2*RealKind
  integer, parameter :: QuadCmplxKind = 2*CmplxKind
#endif
!
end module KindParamModule
