module MathParamModule
!
!  setup oftenly used mathematical parameters
!
   use KindParamModule, only : RealKind
   use KindParamModule, only : CmplxKind
!
   implicit none
!
   private :: RealKind, CmplxKind
!
   real (kind=RealKind), parameter :: zero  = 0.0d0
   real (kind=RealKind), parameter :: tenth = 0.1d0
   real (kind=RealKind), parameter :: eighth= 0.125d0
   real (kind=RealKind), parameter :: fifth = 0.2d0
   real (kind=RealKind), parameter :: fourth= 0.25d0
   real (kind=RealKind), parameter :: half  = 0.5d0
   real (kind=RealKind), parameter :: one   = 1.0d0
   real (kind=RealKind), parameter :: two   = 2.0d0
   real (kind=RealKind), parameter :: three = 3.0d0
   real (kind=RealKind), parameter :: four  = 4.0d0
   real (kind=RealKind), parameter :: five  = 5.0d0
   real (kind=RealKind), parameter :: six   = 6.0d0
   real (kind=RealKind), parameter :: seven = 7.0d0
   real (kind=RealKind), parameter :: eight = 8.0d0
   real (kind=RealKind), parameter :: nine  = 9.0d0
   real (kind=RealKind), parameter :: ten   = 1.0d+01
!  real (kind=RealKind), parameter :: third = 0.33333333333333d0
   real (kind=RealKind), parameter :: third = ONE/THREE
!
   real (kind=RealKind), parameter :: pi    = 3.14159265358979d0
!  real (kind=RealKind), parameter :: pi2   = 6.28318530717959d0
   real (kind=RealKind), parameter :: pi2   = TWO*PI
!  real (kind=RealKind), parameter :: pi3   = 9.42477796076938d0
   real (kind=RealKind), parameter :: pi3   = THREE*PI
!  real (kind=RealKind), parameter :: pi4   = 1.25663706143592d+01
   real (kind=RealKind), parameter :: pi4   = FOUR*PI
!
   real (kind=RealKind), parameter :: ten2m2 = 1.0d-02
   real (kind=RealKind), parameter :: ten2m3 = 1.0d-03
   real (kind=RealKind), parameter :: ten2m4 = 1.0d-04
   real (kind=RealKind), parameter :: ten2m5 = 1.0d-05
   real (kind=RealKind), parameter :: ten2m6 = 1.0d-06
   real (kind=RealKind), parameter :: ten2m7 = 1.0d-07
   real (kind=RealKind), parameter :: ten2m8 = 1.0d-08
   real (kind=RealKind), parameter :: ten2m9 = 1.0d-09
   real (kind=RealKind), parameter :: ten2m10= 1.0d-010
   real (kind=RealKind), parameter :: ten2m11= 1.0d-011
   real (kind=RealKind), parameter :: ten2m12= 1.0d-012
   real (kind=RealKind), parameter :: ten2m13= 1.0d-013
   real (kind=RealKind), parameter :: ten2m14= 1.0d-014
   real (kind=RealKind), parameter :: ten2m16= 1.0d-016
   real (kind=RealKind), parameter :: ten2m20= 1.0d-020
   real (kind=RealKind), parameter :: ten2m30= 1.0d-030
   real (kind=RealKind), parameter :: ten2p10= 1.0d+010
   real (kind=RealKind), parameter :: ten2p14= 1.0d+014
   real (kind=RealKind), parameter :: ten2p30= 1.0d+030
!
   real (kind=RealKind), parameter :: sqrt2      = 1.41421356237310d0
   real (kind=RealKind), parameter :: sqrt3      = 1.73205080756888d0
   real (kind=RealKind), parameter :: sqrt5      = 2.23606797749979d0
   real (kind=RealKind), parameter :: sqrt2_half = 0.70710678118655d0
   real (kind=RealKind), parameter :: sqrt_pi    = 1.772453850905514d0
   real (kind=RealKind), parameter :: sqrt_pio2  = 1.253314137315499d0
!
   real (kind=RealKind), parameter :: Y0    = half/sqrt_pi
   real (kind=RealKind), parameter :: Y0inv = two*sqrt_pi
!
   complex (kind=CmplxKind), parameter :: sqrtm1 = ( 0.0d0,1.0d0)
   complex (kind=CmplxKind), parameter :: czero  = ( 0.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: cone   = ( 1.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: ctwo   = ( 2.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: cthree = ( 3.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: cfour  = ( 4.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: cfive  = ( 5.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: csix   = ( 6.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: cseven = ( 7.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: ceight = ( 8.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: cnine  = ( 9.0d0,0.0d0)
   complex (kind=CmplxKind), parameter :: cten   = (10.0d0,0.0d0)
!
end module MathParamModule
