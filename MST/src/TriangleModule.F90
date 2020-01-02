module TriangleModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initTriangle,      &
          endTriangle,       &
          getNumGaussPoints, &
          getGaussPosi,      &
          getGaussWght,      &
          getArea,           &
          getCorners,        &
          genGaussPoints,    &
          getNumFinestTris
!
private
   type TriangleStruct
        real (kind=RealKind) :: area
        real (kind=RealKind) :: corner_1(3), corner_2(3), corner_3(3)
        integer (kind=IntKind) :: NumSubTriangles
        type (TriangleStruct), pointer :: SubTriangles(:)
   end type TriangleStruct
!
   type PointStruct
        real (kind=RealKind) :: position(3)
        real (kind=RealKind) :: weight
        integer (kind=IntKind) :: index
        type (PointStruct), pointer :: next
        type (PointStruct), pointer :: prev
   end type PointStruct
!
   integer (kind=IntKind) :: Method ! = 1,   4 points ( R = O(h**3) )
                                    ! = 2,   7 points ( R = O(h**4) )
                                    ! = 3,   7 points ( R = O(h**6) )
                                    ! = 4,   6 points ( R = O(h**5) )
                                    ! = 5,  12 points ( R = O(h**7) )
                                    ! = 6,  15 points ( R = O(h**8) )
                                    ! = 7,  16 points ( R = O(h**9) )
                                    ! = 8,  19 points ( R = O(h**10) )
                                    ! = 9,  25 points ( R = O(h**11) )
!
   integer (kind=IntKind) :: TriGaussPts
   integer (kind=IntKind) :: NumGaussPts
   integer (kind=IntKind) :: NumTriangles
   integer (kind=IntKind) :: NumAllLevelTris
!
   real (kind=RealKind) :: epsilon
   real (kind=RealKind) :: MyArea
   real (kind=RealKind), pointer :: wi(:), xi(:,:)
!
   type (PointStruct), pointer :: curr_point
   type (PointStruct), pointer :: first_point
!
   type (TriangleStruct), pointer :: my_triangle
!
!  *******************************************************************
!  Note: knowing the three vertices of a triangle, v1(1;3), v2(1:3), and
!        v3(1:3), the points corresponding to xi(1:3,1:N) on the triangle
!        are given by:
!           xi(1,j)*v1(1)+xi(2,j)*v2(1)+xi(3,j)*v3(1),
!           xi(1,j)*v1(2)+xi(2,j)*v2(2)+xi(3,j)*v3(2),
!           xi(1,j)*v1(3)+xi(2,j)*v2(3)+xi(3,j)*v3(3)
!        where j=1,2,...,N
!  *******************************************************************
!
!  ===========================================================================
!  M. Abramowitz and L. A. Stegun, "Handbook of Mathematical Functions"
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts1=4
   real (kind=RealKind), target :: wi1(NumGaussPts1), xi1(3,NumGaussPts1)
   data wi1/0.75000000000000d+00,                                             &
            0.08333333333333d+00,                                             &
            0.08333333333333d+00,                                             &
            0.08333333333333d+00/
   data xi1/0.33333333333333d+00, 0.33333333333333d+00, 0.33333333333333d+00, &
            1.00000000000000d+00, 0.00000000000000d+00, 0.00000000000000d+00, &
            0.00000000000000d+00, 1.00000000000000d+00, 0.00000000000000d+00, &
            0.00000000000000d+00, 0.00000000000000d+00, 1.00000000000000d+00/
!  ===========================================================================
!  M. Abramowitz and L. A. Stegun, "Handbook of Mathematical Functions"
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts2=7
   real (kind=RealKind), target :: wi2(NumGaussPts2), xi2(3,NumGaussPts2)
   data wi2/0.45000000000000d+00,                                             &
            0.05000000000000d+00,                                             &
            0.05000000000000d+00,                                             &
            0.05000000000000d+00,                                             &
            0.13333333333333d+00,                                             &
            0.13333333333333d+00,                                             &
            0.13333333333333d+00/
   data xi2/0.33333333333333d+00, 0.33333333333333d+00, 0.33333333333333d+00, &
            1.00000000000000d+00, 0.00000000000000d+00, 0.00000000000000d+00, &
            0.00000000000000d+00, 1.00000000000000d+00, 0.00000000000000d+00, &
            0.00000000000000d+00, 0.00000000000000d+00, 1.00000000000000d+00, &
            0.50000000000000d+00, 0.50000000000000d+00, 0.00000000000000d+00, &
            0.00000000000000d+00, 0.50000000000000d+00, 0.50000000000000d+00, &
            0.50000000000000d+00, 0.00000000000000d+00, 0.50000000000000d+00/
!  ===========================================================================
!  G. R. Cowper, Int. J. Num. Mech. Engng, vol 7, 405 (1973).
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts3=7
   real (kind=RealKind), target :: wi3(NumGaussPts3), xi3(3,NumGaussPts3)
   data wi3/0.22500000000000d+00,                                             &
            0.12593918054483d+00,                                             &
            0.12593918054483d+00,                                             &
            0.12593918054483d+00,                                             &
            0.13239415278851d+00,                                             &
            0.13239415278851d+00,                                             &
            0.13239415278851d+00/
   data xi3/0.33333333333333d+00, 0.33333333333333d+00, 0.33333333333333d+00, &
            0.79742698535308d+00, 0.10128650732346d+00, 0.10128650732346d+00, &
            0.10128650732346d+00, 0.79742698535308d+00, 0.10128650732346d+00, &
            0.10128650732346d+00, 0.10128650732346d+00, 0.79742698535308d+00, &
            0.59715871789771d-01, 0.47014206410511d+00, 0.47014206410511d+00, &
            0.47014206410511d+00, 0.59715871789771d-01, 0.47014206410511d+00, &
            0.47014206410511d+00, 0.47014206410511d+00, 0.59715871789771d-01/
!  ===========================================================================
!  G. R. Cowper, Int. J. Num. Mech. Engng, vol 7, 405 (1973).
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts4=6
   real (kind=RealKind), target :: wi4(NumGaussPts4), xi4(3,NumGaussPts4)
   data wi4/0.10995174365532d+00,                                             &
            0.10995174365532d+00,                                             &
            0.10995174365532d+00,                                             &
            0.22338158967801d+00,                                             &
            0.22338158967801d+00,                                             &
            0.22338158967801d+00/
   data xi4/0.81684757298046d+00, 0.09157621350977d+00, 0.09157621350977d+00, &
            0.09157621350977d+00, 0.81684757298046d+00, 0.09157621350977d+00, &
            0.09157621350977d+00, 0.09157621350977d+00, 0.81684757298046d+00, &
            0.10810301816807d+00, 0.44594849091597d+00, 0.44594849091597d+00, &
            0.44594849091597d+00, 0.10810301816807d+00, 0.44594849091597d+00, &
            0.44594849091597d+00, 0.44594849091597d+00, 0.10810301816807d+00/
!  ===========================================================================
!  G. R. Cowper, Int. J. Num. Mech. Engng, vol 7, 405 (1973).
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts5=12
   real (kind=RealKind), target :: wi5(NumGaussPts5), xi5(3,NumGaussPts5)
   data wi5/0.05084490637021d+00,                                             &
            0.05084490637021d+00,                                             &
            0.05084490637020d+00,                                             &
            0.11678627572638d+00,                                             &
            0.11678627572638d+00,                                             &
            0.11678627572638d+00,                                             &
            0.08285107561837d+00,                                             &
            0.08285107561837d+00,                                             &
            0.08285107561838d+00,                                             &
            0.08285107561837d+00,                                             &
            0.08285107561837d+00,                                             &
            0.08285107561838d+00/
   data xi5/0.87382197101700d+00, 0.06308901449150d+00, 0.06308901449150d+00, &
            0.06308901449150d+00, 0.87382197101700d+00, 0.06308901449150d+00, &
            0.06308901449150d+00, 0.06308901449150d+00, 0.87382197101700d+00, &
            0.50142650965818d+00, 0.24928674517091d+00, 0.24928674517091d+00, &
            0.24928674517091d+00, 0.50142650965818d+00, 0.24928674517091d+00, &
            0.24928674517091d+00, 0.24928674517091d+00, 0.50142650965818d+00, &
            0.63650249912140d+00, 0.31035245103378d+00, 0.05314504984482d+00, &
            0.63650249912140d+00, 0.05314504984482d+00, 0.31035245103378d+00, &
            0.31035245103378d+00, 0.63650249912140d+00, 0.05314504984482d+00, &
            0.05314504984482d+00, 0.63650249912140d+00, 0.31035245103378d+00, &
            0.31035245103378d+00, 0.05314504984482d+00, 0.63650249912140d+00, &
            0.05314504984482d+00, 0.31035245103378d+00, 0.63650249912140d+00/
!  ===========================================================================
!  M. E. Laursen and M. Gellert, Int. J. Num. Mech. Engng, vol 12, 67 (1978).
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts6=15
   real (kind=RealKind), target :: wi6(NumGaussPts6), xi6(3,NumGaussPts6)
   data wi6/0.05307780179023d+00,                                             &
            0.05307780179023d+00,                                             &
            0.05307780179023d+00,                                             &
            0.07085308369213d+00,                                             &
            0.07085308369214d+00,                                             &
            0.07085308369213d+00,                                             &
            0.07085308369214d+00,                                             &
            0.07085308369213d+00,                                             &
            0.07085308369214d+00,                                             &
            0.06927468207941d+00,                                             &
            0.06927468207942d+00,                                             &
            0.06927468207941d+00,                                             &
            0.06927468207942d+00,                                             &
            0.06927468207941d+00,                                             &
            0.06927468207942d+00/
   data xi6/0.87013897368167d+00, 0.06493051315916d+00, 0.06493051315917d+00, &
            0.06493051315916d+00, 0.87013897368167d+00, 0.06493051315917d+00, &
            0.06493051315916d+00, 0.06493051315917d+00, 0.87013897368167d+00, &
            0.28457558424917d+00, 0.51703993906933d+00, 0.19838447668150d+00, &
            0.51703993906933d+00, 0.28457558424917d+00, 0.19838447668150d+00, &
            0.28457558424917d+00, 0.19838447668150d+00, 0.51703993906933d+00, &
            0.51703993906933d+00, 0.19838447668150d+00, 0.28457558424917d+00, &
            0.19838447668150d+00, 0.28457558424917d+00, 0.51703993906933d+00, &
            0.19838447668150d+00, 0.51703993906933d+00, 0.28457558424917d+00, &
            0.31355918438493d+00, 0.04386347179237d+00, 0.64257734382270d+00, &
            0.04386347179237d+00, 0.31355918438493d+00, 0.64257734382270d+00, &
            0.31355918438493d+00, 0.64257734382270d+00, 0.04386347179237d+00, &
            0.04386347179237d+00, 0.64257734382270d+00, 0.31355918438493d+00, &
            0.64257734382270d+00, 0.31355918438493d+00, 0.04386347179237d+00, &
            0.64257734382270d+00, 0.04386347179237d+00, 0.31355918438493d+00/
!  ===========================================================================
!  M. E. Laursen and M. Gellert, Int. J. Num. Mech. Engng, vol 12, 67 (1978).
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts7=16
   real (kind=RealKind), target :: wi7(NumGaussPts7), xi7(3,NumGaussPts7)
   data wi7/0.14431560767779d+00,                                             &
            0.10321737053472d+00,                                             &
            0.10321737053472d+00,                                             &
            0.10321737053471d+00,                                             &
            0.03245849762320d+00,                                             &
            0.03245849762320d+00,                                             &
            0.03245849762319d+00,                                             &
            0.09509163426728d+00,                                             &
            0.09509163426728d+00,                                             &
            0.09509163426729d+00,                                             &
            0.02723031417444d+00,                                             &
            0.02723031417443d+00,                                             &
            0.02723031417443d+00,                                             &
            0.02723031417444d+00,                                             &
            0.02723031417443d+00,                                             &
            0.02723031417443d+00/
   data xi7/0.33333333333333d+00, 0.33333333333333d+00, 0.33333333333333d+00, &
            0.65886138449648d+00, 0.17056930775176d+00, 0.17056930775176d+00, &
            0.17056930775176d+00, 0.65886138449648d+00, 0.17056930775176d+00, &
            0.17056930775176d+00, 0.17056930775176d+00, 0.65886138449648d+00, &
            0.89890554336594d+00, 0.05054722831703d+00, 0.05054722831703d+00, &
            0.05054722831703d+00, 0.89890554336594d+00, 0.05054722831703d+00, &
            0.05054722831703d+00, 0.05054722831703d+00, 0.89890554336594d+00, &
            0.08141482341455d+00, 0.45929258829272d+00, 0.45929258829272d+00, &
            0.45929258829272d+00, 0.08141482341455d+00, 0.45929258829272d+00, &
            0.45929258829272d+00, 0.45929258829272d+00, 0.08141482341455d+00, &
            0.00839477740996d+00, 0.26311282963464d+00, 0.72849239295540d+00, &
            0.00839477740996d+00, 0.72849239295540d+00, 0.26311282963464d+00, &
            0.26311282963464d+00, 0.00839477740996d+00, 0.72849239295540d+00, &
            0.26311282963464d+00, 0.72849239295540d+00, 0.00839477740996d+00, &
            0.72849239295540d+00, 0.00839477740996d+00, 0.26311282963464d+00, &
            0.72849239295540d+00, 0.26311282963464d+00, 0.00839477740996d+00/
!  ===========================================================================
!  M. E. Laursen and M. Gellert, Int. J. Num. Mech. Engng, vol 12, 67 (1978).
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts8=19
   real (kind=RealKind), target :: wi8(NumGaussPts8), xi8(3,NumGaussPts8)
   data wi8/0.09713579628280d+00,                                             &
            0.03133470022714d+00,                                             &
            0.03133470022714d+00,                                             &
            0.03133470022714d+00,                                             &
            0.07782754100477d+00,                                             &
            0.07782754100477d+00,                                             &
            0.07782754100478d+00,                                             &
            0.07964773892721d+00,                                             &
            0.07964773892721d+00,                                             &
            0.07964773892721d+00,                                             &
            0.02557767565870d+00,                                             &
            0.02557767565870d+00,                                             &
            0.02557767565869d+00,                                             &
            0.04328353937729d+00,                                             &
            0.04328353937729d+00,                                             &
            0.04328353937729d+00,                                             &
            0.04328353937729d+00,                                             &
            0.04328353937729d+00,                                             &
            0.04328353937728d+00/
   data xi8/0.33333333333333d+00, 0.33333333333333d+00, 0.33333333333333d+00, &
            0.02063496160252d+00, 0.48968251919874d+00, 0.48968251919874d+00, &
            0.48968251919874d+00, 0.02063496160252d+00, 0.48968251919874d+00, &
            0.48968251919874d+00, 0.48968251919874d+00, 0.02063496160252d+00, &
            0.12582081701413d+00, 0.43708959149294d+00, 0.43708959149294d+00, &
            0.43708959149294d+00, 0.12582081701413d+00, 0.43708959149294d+00, &
            0.43708959149294d+00, 0.43708959149294d+00, 0.12582081701413d+00, &
            0.62359292876193d+00, 0.18820353561903d+00, 0.18820353561903d+00, &
            0.18820353561903d+00, 0.62359292876193d+00, 0.18820353561903d+00, &
            0.18820353561903d+00, 0.18820353561903d+00, 0.62359292876193d+00, &
            0.91054097321109d+00, 0.04472951339445d+00, 0.04472951339445d+00, &
            0.04472951339445d+00, 0.91054097321109d+00, 0.04472951339445d+00, &
            0.04472951339445d+00, 0.04472951339445d+00, 0.91054097321109d+00, &
            0.03683841205474d+00, 0.22196298916077d+00, 0.74119859878450d+00, &
            0.03683841205474d+00, 0.74119859878450d+00, 0.22196298916077d+00, &
            0.22196298916077d+00, 0.03683841205474d+00, 0.74119859878450d+00, &
            0.22196298916077d+00, 0.74119859878450d+00, 0.03683841205474d+00, &
            0.74119859878450d+00, 0.22196298916077d+00, 0.03683841205474d+00, &
            0.74119859878450d+00, 0.03683841205474d+00, 0.22196298916077d+00/
!  ===========================================================================
!  M. E. Laursen and M. Gellert, Int. J. Num. Mech. Engng, vol 12, 67 (1978).
!  ===========================================================================
   integer (kind=IntKind), parameter :: NumGaussPts9=25
   real (kind=RealKind), target :: wi9(NumGaussPts9), xi9(3,NumGaussPts9)
   data wi9/0.07989450474124d+00,                                             &
            0.07112380223238d+00,                                             &
            0.07112380223238d+00,                                             &
            0.07112380223237d+00,                                             &
            0.00822381869046d+00,                                             &
            0.00822381869046d+00,                                             &
            0.00822381869047d+00,                                             &
            0.04543059229617d+00,                                             &
            0.04543059229617d+00,                                             &
            0.04543059229617d+00,                                             &
            0.04543059229617d+00,                                             &
            0.04543059229617d+00,                                             &
            0.04543059229617d+00,                                             &
            0.03735985623431d+00,                                             &
            0.03735985623430d+00,                                             &
            0.03735985623431d+00,                                             &
            0.03735985623430d+00,                                             &
            0.03735985623431d+00,                                             &
            0.03735985623430d+00,                                             &
            0.03088665688456d+00,                                             &
            0.03088665688456d+00,                                             &
            0.03088665688457d+00,                                             &
            0.03088665688456d+00,                                             &
            0.03088665688456d+00,                                             &
            0.03088665688457d+00/
   data xi9/0.33333333333333d+00, 0.33333333333333d+00, 0.33333333333333d+00, &
            0.14982757879582d+00, 0.42508621060209d+00, 0.42508621060209d+00, &
            0.42508621060209d+00, 0.14982757879582d+00, 0.42508621060209d+00, &
            0.42508621060209d+00, 0.42508621060209d+00, 0.14982757879582d+00, &
            0.95338226498000d+00, 0.02330886751000d+00, 0.02330886751000d+00, &
            0.02330886751000d+00, 0.95338226498000d+00, 0.02330886751000d+00, &
            0.02330886751000d+00, 0.02330886751000d+00, 0.95338226498000d+00, &
            0.14792562620953d+00, 0.22376697357697d+00, 0.62830740021349d+00, &
            0.14792562620953d+00, 0.62830740021349d+00, 0.22376697357697d+00, &
            0.22376697357697d+00, 0.14792562620953d+00, 0.62830740021349d+00, &
            0.22376697357697d+00, 0.62830740021349d+00, 0.14792562620953d+00, &
            0.62830740021349d+00, 0.14792562620953d+00, 0.22376697357697d+00, &
            0.62830740021349d+00, 0.22376697357697d+00, 0.14792562620953d+00, &
            0.02994603195417d+00, 0.35874014186443d+00, 0.61131382618140d+00, &
            0.02994603195417d+00, 0.61131382618140d+00, 0.35874014186443d+00, &
            0.35874014186443d+00, 0.02994603195417d+00, 0.61131382618140d+00, &
            0.35874014186443d+00, 0.61131382618140d+00, 0.02994603195417d+00, &
            0.61131382618140d+00, 0.02994603195417d+00, 0.35874014186443d+00, &
            0.61131382618140d+00, 0.35874014186443d+00, 0.02994603195417d+00, &
            0.03563255958750d+00, 0.14329537042687d+00, 0.82107206998563d+00, &
            0.03563255958750d+00, 0.82107206998563d+00, 0.14329537042687d+00, &
            0.14329537042687d+00, 0.03563255958750d+00, 0.82107206998563d+00, &
            0.14329537042687d+00, 0.82107206998563d+00, 0.03563255958750d+00, &
            0.82107206998563d+00, 0.14329537042687d+00, 0.03563255958750d+00, &
            0.82107206998563d+00, 0.03563255958750d+00, 0.14329537042687d+00/
!  ===========================================================================
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initTriangle(t1,t2,t3,PSMethod)
!  ===================================================================
   use MathParamModule, only : ZERO, TEN2m8, ONE
!
   implicit   none
!
   real (kind=RealKind), intent(in) :: t1(3), t2(3), t3(3)
   integer (kind=IntKind), intent(in) :: PSMethod
!
   integer (kind=IntKind) :: i
   real (kind=RealKind) :: tw
!
   Method=PSMethod
!
   if (Method.eq.1) then
      TriGaussPts=NumGaussPts1
      wi=>wi1(1:TriGaussPts)
      xi=>xi1(1:3,1:TriGaussPts)
   else if(Method.eq.2) then
      TriGaussPts=NumGaussPts2
      wi=>wi2(1:TriGaussPts)
      xi=>xi2(1:3,1:TriGaussPts)
   else if(Method.eq.3) then
      TriGaussPts=NumGaussPts3
      wi=>wi3(1:TriGaussPts)
      xi=>xi3(1:3,1:TriGaussPts)
   else if(Method.eq.4) then
      TriGaussPts=NumGaussPts4
      wi=>wi4(1:TriGaussPts)
      xi=>xi4(1:3,1:TriGaussPts)
   else if(Method.eq.5) then
      TriGaussPts=NumGaussPts5
      wi=>wi5(1:TriGaussPts)
      xi=>xi5(1:3,1:TriGaussPts)
   else if(Method.eq.6) then
      TriGaussPts=NumGaussPts6
      wi=>wi6(1:TriGaussPts)
      xi=>xi6(1:3,1:TriGaussPts)
   else if(Method.eq.7) then
      TriGaussPts=NumGaussPts7
      wi=>wi7(1:TriGaussPts)
      xi=>xi7(1:3,1:TriGaussPts)
   else if(Method.eq.8) then
      TriGaussPts=NumGaussPts8
      wi=>wi8(1:TriGaussPts)
      xi=>xi8(1:3,1:TriGaussPts)
   else if(Method.eq.9) then
      TriGaussPts=NumGaussPts9
      wi=>wi9(1:TriGaussPts)
      xi=>xi9(1:3,1:TriGaussPts)
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initTriangle','Invalid method',Method)
!     ----------------------------------------------------------------
   endif
!
   tw = ZERO
   do i=1,TriGaussPts
      if (abs(xi(1,i)+xi(2,i)+xi(3,i)-ONE) > TEN2m8) then
!        -------------------------------------------------------------
         call ErrorHandler('initTriangle','Invalid position',Method)
!        -------------------------------------------------------------
      endif
      tw = tw+wi(i)
   enddo
   if (abs(tw-ONE) > TEN2m8) then
!     ----------------------------------------------------------------
      call ErrorHandler('initTriangle','Invalid weight',Method)
!     ----------------------------------------------------------------
   endif
!
   NumAllLevelTris=0
   allocate(my_triangle)
!  -------------------------------------------------------------------
   call setupTriangle(t1,t2,t3,my_triangle)
!  -------------------------------------------------------------------
   MyArea = my_triangle%area
   NumTriangles=0
   NumGaussPts = 0
   nullify(curr_point, first_point)
!
   end subroutine initTriangle
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupTriangle(t1,t2,t3,triangle)
!  ===================================================================
   use MathParamModule, only : HALF
   implicit   none
!
   real (kind=RealKind), intent(in) :: t1(3), t2(3), t3(3)
   real (kind=RealKind) :: v12(3), v13(3), v(3)
!
   type (TriangleStruct) :: triangle
!
   triangle%corner_1(1:3)=t1(1:3)
   triangle%corner_2(1:3)=t2(1:3)
   triangle%corner_3(1:3)=t3(1:3)
!
   v12(1:3)=t1(1:3)-t2(1:3)
   v13(1:3)=t1(1:3)-t3(1:3)
!  -------------------------------------------------------------------
   call vcross(v12,v13,v)
!  -------------------------------------------------------------------
   triangle%area=HALF*sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
   nullify(triangle%SubTriangles)
   triangle%NumSubTriangles = 0
   NumAllLevelTris = NumAllLevelTris + 1
!
   end subroutine setupTriangle
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endTriangle()
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind) :: i
!
!  -------------------------------------------------------------------
   call deleteSubTriangle(my_triangle)
!  -------------------------------------------------------------------
   deallocate(my_triangle)
!
   if (NumAllLevelTris /= 0) then
      call ErrorHandler('endTriangle','Triangles are not completely deleted')
   endif
   do i=1,NumGaussPts
      curr_point => first_point%next
      nullify(first_point%prev)
      deallocate(first_point)
      first_point => curr_point
   enddo
!  -------------------------------------------------------------------
   nullify(wi, xi)
!  -------------------------------------------------------------------
   NumGaussPts = 0
   NumTriangles = 0
   end subroutine endTriangle
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   recursive subroutine deleteSubTriangle(triangle)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind) :: i
!
   type (TriangleStruct) :: triangle
!
   do i=1,triangle%NumSubTriangles
!     ----------------------------------------------------------------
      call deleteSubTriangle( triangle%SubTriangles(i) )
!     ----------------------------------------------------------------
   enddo
   if (triangle%NumSubTriangles > 0) then
!     ----------------------------------------------------------------
      deallocate(triangle%SubTriangles)
!     ----------------------------------------------------------------
   endif
!
   NumAllLevelTris=NumAllLevelTris-1
!
   end subroutine deleteSubTriangle
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getArea() result(a)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: a
!
   a = my_triangle%area
!
   end function getArea
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumGaussPoints() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumGaussPts  
!
   end function getNumGaussPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussPosi(i) result(x)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: x(3)
!
   if (i < 1 .or. i > NumGaussPts) then
      call ErrorHandler('getGaussPointPosi','i out of range',i)
   else if (i == 1) then
      curr_point => first_point
   else
!     ----------------------------------------------------------------
      call movePointer(i)
!     ----------------------------------------------------------------
   endif
   x(1:3)=curr_point%position(1:3)
!  x(1)=xi(1,i)*Corners(1,1)+xi(2,i)*Corners(1,2)+xi(3,i)*Corners(1,3)
!  x(2)=xi(1,1)*Corners(2,1)+xi(2,i)*Corners(2,2)+xi(3,i)*Corners(2,3)
!  x(3)=xi(1,1)*Corners(3,1)+xi(2,i)*Corners(3,2)+xi(3,i)*Corners(3,3)
!
   end function getGaussPosi
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGaussWght(i) result(w)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   real (kind=RealKind) :: w
!
   if (i < 1 .or. i > NumGaussPts) then
      call ErrorHandler('getGaussPointWght','i out of range',i)
   else if (i == 1) then
      curr_point => first_point
   else
!     ----------------------------------------------------------------
      call movePointer(i)
!     ----------------------------------------------------------------
   endif
   w = curr_point%weight
!  w = wi(i)*Area
!
   end function getGaussWght
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getCorners() result(c)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: c(3,3)
!
   c(1:3,1) = my_triangle%corner_1(1:3)
   c(1:3,2) = my_triangle%corner_2(1:3)
   c(1:3,3) = my_triangle%corner_3(1:3)
!
   end function getCorners
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumFinestTris() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   n = NumTriangles
!
   end function getNumFinestTris
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine movePointer(i)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: i
   do
      if (curr_point%index == i) then
         exit
      else if (curr_point%index > i) then
         curr_point => curr_point%prev
      else
         curr_point => curr_point%next
      endif
   enddo
   end subroutine movePointer
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genGaussPoints(e)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: e
!
   epsilon = e
!
   call divtri(my_triangle)
!
   end subroutine genGaussPoints
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   recursive subroutine divtri(triangle)
!  ===================================================================
   use MathParamModule, only : ten2m8, ten2m7, ten2m6
   use MathParamModule, only : third
   use MathParamModule, only : HALF, one
!
   implicit none
!
!  *******************************************************************
!  tri_corner_1, tri_corner_2, tri_corner_3: the traingle corner points
!  triangle: structure that contains info about a triangle
!  *******************************************************************
!
   type (TriangleStruct) :: triangle  
!
   type (PointStruct), pointer :: new_point
!
   logical :: CheckPoint
!
   integer (kind=IntKind) :: i, n
!
   real (kind=RealKind) :: t1(3), t2(3), t3(3), tnew(3)
   real (kind=RealKind) :: d1, d2, d3, x, y, z
!
   if (triangle%area >= epsilon) then
      if (HALF*triangle%area < epsilon) then
!        =============================================================
!        determine the longest side of the current triangle
!        one of the sub-triangle corners is on this side
!        =============================================================
         d1=(triangle%corner_2(1)-triangle%corner_3(1))**2+               &
            (triangle%corner_2(2)-triangle%corner_3(2))**2+               &
            (triangle%corner_2(3)-triangle%corner_3(3))**2
         d2=(triangle%corner_1(1)-triangle%corner_3(1))**2+               &
            (triangle%corner_1(2)-triangle%corner_3(2))**2+               &
            (triangle%corner_1(3)-triangle%corner_3(3))**2
         d3=(triangle%corner_1(1)-triangle%corner_2(1))**2+               &
            (triangle%corner_1(2)-triangle%corner_2(2))**2+               &
            (triangle%corner_1(3)-triangle%corner_2(3))**2
         if(d1.gt.d2 .and. d1.gt.d3) then
            t1(1:3)=triangle%corner_2(1:3)
            t2(1:3)=triangle%corner_3(1:3)
            t3(1:3)=triangle%corner_1(1:3)
         else if(d2.gt.d1 .and. d2.gt.d3) then
            t1(1:3)=triangle%corner_1(1:3)
            t2(1:3)=triangle%corner_3(1:3)
            t3(1:3)=triangle%corner_2(1:3)
         else
            t1(1:3)=triangle%corner_1(1:3)
            t2(1:3)=triangle%corner_2(1:3)
            t3(1:3)=triangle%corner_3(1:3)
         endif
         tnew(1)=HALF*(t1(1)+t2(1))
         tnew(2)=HALF*(t1(2)+t2(2))
         tnew(3)=HALF*(t1(3)+t2(3))
!
         n = 2
         triangle%NumSubTriangles = n
         allocate(triangle%SubTriangles(n))
         call setupTriangle(t1,tnew,t3,triangle%SubTriangles(1))
         call setupTriangle(t2,tnew,t3,triangle%SubTriangles(2))
      else
         t1(1:3)=HALF*(triangle%corner_1(1:3)+triangle%corner_2(1:3))
         t2(1:3)=HALF*(triangle%corner_2(1:3)+triangle%corner_3(1:3))
         t3(1:3)=HALF*(triangle%corner_1(1:3)+triangle%corner_3(1:3))
!
         n = 4
         triangle%NumSubTriangles = n
         allocate(triangle%SubTriangles(n))
         call setupTriangle(triangle%corner_1,t1,t3,triangle%SubTriangles(1))
         call setupTriangle(t1,triangle%corner_2,t2,triangle%SubTriangles(2))
         call setupTriangle(t3,t2,triangle%corner_3,triangle%SubTriangles(3))
         call setupTriangle(t1,t2,t3,triangle%SubTriangles(4))
      endif
!
      do i=1,triangle%NumSubTriangles
!        =============================================================
!        Recursively call divtri to divide the triangle.
!        -------------------------------------------------------------
         call divtri(triangle%SubTriangles(i))
!        -------------------------------------------------------------
      enddo
   else
      NumTriangles=NumTriangles+1
      do i=1,TriGaussPts
         x=xi(1,i)*triangle%corner_1(1)+xi(2,i)*triangle%corner_2(1)+  &
           xi(3,i)*triangle%corner_3(1)
         y=xi(1,i)*triangle%corner_1(2)+xi(2,i)*triangle%corner_2(2)+  &
           xi(3,i)*triangle%corner_3(2)
         z=xi(1,i)*triangle%corner_1(3)+xi(2,i)*triangle%corner_2(3)+  &
           xi(3,i)*triangle%corner_3(3)
         if ((abs(xi(1,i)-ONE) < TEN2m8 .or. abs(xi(2,i)-ONE) < TEN2m8 &
             .or. abs(xi(3,i)-ONE) < TEN2m8) .and. NumGaussPts > 3) then
            CheckPoint = .true.
         else
            CheckPoint = .false.
         endif
!        -------------------------------------------------------------
         call allocPoint(x,y,z,CheckPoint,new_point)
!        -------------------------------------------------------------
         new_point%weight=new_point%weight+wi(i)*triangle%area
      enddo
   endif
!
   end subroutine divtri
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocPoint(x,y,z,CheckPoint,new_point)
!  ===================================================================
   use MathParamModule, only : ZERO, TEN2m8
!
   implicit none
!
   logical, intent(in) :: CheckPoint
!
   real (kind=RealKind), intent(in) :: x,y,z
!
   type (PointStruct), pointer :: new_point
!
!  ===================================================================
!  determine if the new point is alrealdy on the list.................
!  ===================================================================
   if (CheckPoint) then
      new_point=>first_point
      do while (new_point%index < NumGaussPts)
         if(abs(x-new_point%position(1))+abs(y-new_point%position(2))+ &
            abs(z-new_point%position(3)) < TEN2m8) then
            return
         else
            new_point=>new_point%next
         endif
      enddo
   endif
!
   NumGaussPts=NumGaussPts+1
!  -------------------------------------------------------------------
   allocate(new_point)
!  -------------------------------------------------------------------
   if (NumGaussPts == 1) then
      first_point => new_point
      nullify(new_point%prev)
   else
      curr_point%next => new_point
      new_point%prev => curr_point
   endif
!  -------------------------------------------------------------------
   nullify(new_point%next)
!  -------------------------------------------------------------------
   new_point%weight=ZERO
   new_point%index=NumGaussPts
   new_point%position(1)=x
   new_point%position(2)=y
   new_point%position(3)=z
   curr_point => new_point
!
   end subroutine allocPoint
!  ===================================================================
end module TriangleModule
