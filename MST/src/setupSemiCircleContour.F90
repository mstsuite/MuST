!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupSemiCircleContour(ne,er,eg,ew)
!  ===================================================================
!  Set up a semi-circle contour around (0.0, 0.0) with radius er by
!  generating ne number of Gaussian points, eg with weight ew, along the
!  contour.
!      Input:  ne -- the number of Gaussian quadrature points
!              er -- the radius of the contour
!      Output: eg -- the Gaussian quadrature points
!              ew -- the weight of each Gaussian quadrature point
!  *******************************************************************
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, HALF, ONE, TWO, PI, CZERO, SQRTm1, TEN2m6
!
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ne
   integer (kind=IntKind) :: ie
!
   real (kind=RealKind), intent(in) :: er
   real (kind=RealKind) :: xg(ne), wg(ne)
!
   complex (kind=CmplxKind), intent(out) :: eg(ne), ew(ne)
   complex (kind=CmplxKind) :: es, ec, check_wght
!
!  -------------------------------------------------------------------
   call gauleg(-ONE, ONE, xg, wg, ne)
!  -------------------------------------------------------------------
   ec = PI*HALF*SQRTm1
   check_wght = CZERO
   do ie = 1, ne
      es = er*exp(ec*(ONE-xg(ie)))
      eg(ie) = es
      ew(ie) = -ec*es*wg(ie)
      check_wght = check_wght + ew(ie)
   enddo
   if (abs(check_wght-TWO*er) > TEN2m6) then
      call ErrorHandler('setupSemiCircleContour','Weight check is failed', &
                        check_wght)
   endif
!
   end subroutine setupSemiCircleContour
