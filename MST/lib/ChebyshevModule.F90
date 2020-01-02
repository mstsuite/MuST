module ChebyshevModule
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, HALF, TWO, PI, SQRTM1
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
public :: ChebyshevStruct,     &
          initChebyshevSeries, &
          ChebyshevEval,       &
          calChebGrid,         &
          endChebyshevSeries
!
   interface initChebyshevSeries
      module procedure initChebyshevSeries_r,initChebyshevSeries_c
   end interface initChebyshevSeries
!
   interface ChebyshevEval
      module procedure ChebyshevEval_r, ChebyshevEval_c
   end interface ChebyshevEval
!
private
   type ChebyshevStruct
      logical :: isCmplx
      integer (kind=IntKind) :: order
      real (kind=RealKind) :: LowBound
      real (kind=RealKind) :: UpBound
      real (kind=RealKind), pointer :: ChebCoef_R(:)
      real (kind=RealKind), pointer :: ChebCoef_C(:)
      real (kind=RealKind), pointer :: ChebFunc_R(:)
      complex (kind=CmplxKind), pointer :: ChebFunc_C(:)
   end type ChebyshevStruct
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initChebyshevSeries_r(order,a,b,f,cheb_struct)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: order
   real (kind=RealKind), intent(in) :: a, b
   real (kind=RealKind), intent(in) :: f(1:order)
   type (ChebyshevStruct), target :: cheb_struct
!
   integer (kind=IntKind) :: k, j
   real (kind=RealKind) :: bma, bpa, fact, sum_j
!
   cheb_struct%isCmplx  = .false.
   cheb_struct%order    = order
   cheb_struct%LowBound = a
   cheb_struct%UpBound  = b
!
   allocate( cheb_struct%ChebCoef_R( 1:order ) )
   nullify( cheb_struct%ChebCoef_C )
   allocate( cheb_struct%ChebFunc_R( 1:order ) )
   nullify( cheb_struct%ChebFunc_C )
!
   if (a >= b) then
     call ErrorHandler("initChebyshevSeries_c","Wrong interval [a,b]", a, b)
   endif
!
   cheb_struct%LowBound = a
   cheb_struct%UpBound = b
!
   bma = HALF*( b - a )
   bpa = HALF*( b + a )
   fact = TWO/real(order,kind=RealKind)
!
   do j = 1,order
      sum_j = ZERO
      do k = 1, order
        sum_j = sum_j + f(k)*cos( PI*(j-1)*(k-HALF)/real(order,kind=RealKind) )
      enddo
      cheb_struct%ChebCoef_r(j) = fact*sum_j
   enddo
!
   end subroutine initChebyshevSeries_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initChebyshevSeries_c(order,a,b,f,cheb_struct)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: order
   real (kind=RealKind), intent(in) :: a, b
   complex (kind=CmplxKind), intent(in) :: f(1:order)
   type (ChebyshevStruct), target :: cheb_struct
!
   integer (kind=IntKind) :: k,j
   real (kind=RealKind) :: bma, bpa, fact, sum_jr, sum_jc
!
   cheb_struct%isCmplx  = .true.
   cheb_struct%order    = order
   cheb_struct%LowBound = a
   cheb_struct%UpBound  = b
!
   allocate( cheb_struct%ChebCoef_R( 1:order ) )
   allocate( cheb_struct%ChebCoef_C( 1:order ) )
   allocate( cheb_struct%ChebFunc_R( 1:order ) )
   allocate( cheb_struct%ChebFunc_C( 1:order ) )
!
   if (a >= b) then
     call ErrorHandler("initChebyshevSeries_c","Wrong interval [a,b]", a, b)
   endif
!
   cheb_struct%LowBound = a
   cheb_struct%UpBound = b
!
   bma = HALF*( b - a )
   bpa = HALF*( b + a )
   fact = TWO/real(order,kind=RealKind)
!
   do j = 1,order
      sum_jr = ZERO
      sum_jc = ZERO
      do k = 1, order
        sum_jr = sum_jr + real(f(k),kind=RealKind)*                   &
                    cos( PI*(j-1)*(k-HALF)/real(order,kind=RealKind) )
        sum_jc = sum_jc - real(f(k)*sqrtm1,kind=RealKind)*            &
                    cos( PI*(j-1)*(k-HALF)/real(order,kind=RealKind) )
      enddo
      cheb_struct%ChebCoef_r(j) = fact*sum_jr
      cheb_struct%ChebCoef_c(j) = fact*sum_jc
   enddo
!
   end subroutine initChebyshevSeries_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ChebyshevEval_r(cheb_struct,max_order,n,x,f)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), optional :: max_order
   real (kind=RealKind), intent(in) :: x(n)
   real (kind=RealKind), target :: f(1:n)
   type (ChebyshevStruct), target :: cheb_struct
!
   integer (kind=IntKind) :: j,i
   real (kind=RealKind) :: a, b, y, yy
   real (kind=RealKind) :: d, dd, temp
!
   a = cheb_struct%LowBound
   b = cheb_struct%UpBound
!
   do j = 1, n
      d  = ZERO
      dd = ZERO
!
      y  = ( TWO*x(j) - a - b ) / ( b - a )
      yy = TWO*y;
      do i = min(max_order,cheb_struct%order),2,-1
         temp = d
         d  = yy*d - dd + cheb_struct%ChebCoef_r(i)
         dd = temp
      enddo
!
      f(j) = y*d - dd + HALF*cheb_struct%ChebCoef_r(1)
   enddo
!
   end subroutine ChebyshevEval_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ChebyshevEval_c(cheb_struct,max_order,n,x,f)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind), optional :: max_order
   real (kind=RealKind), intent(in) :: x(n)
   complex (kind=CmplxKind), target :: f(1:n)
   type (ChebyshevStruct), target :: cheb_struct
!
   integer (kind=IntKind) :: j,i
   real (kind=RealKind) :: f_r, f_c
   real (kind=RealKind) :: a, b, y, yy
   real (kind=RealKind) :: d_r, dd_r, temp_r
   real (kind=RealKind) :: d_c, dd_c, temp_c
!
   a = cheb_struct%LowBound
   b = cheb_struct%UpBound
!
   do j = 1, n
      d_r = ZERO
      dd_r = ZERO
      d_c = ZERO
      dd_c = ZERO
!
      y  = ( TWO*x(j) - a - b ) / ( b - a )
      yy = TWO*y
      do i = min(max_order,cheb_struct%order),2,-1
!
         temp_r = d_r
         d_r = (yy*d_r + cheb_struct%ChebCoef_r(i)) - dd_r 
         dd_r = temp_r
!
         temp_c = d_c
         d_c = (yy*d_c + cheb_struct%ChebCoef_c(i)) - dd_c
         dd_c = temp_c
      enddo
      f_r = (y*d_r + HALF*cheb_struct%ChebCoef_r(1)) - dd_r
      f_c = (y*d_c + HALF*cheb_struct%ChebCoef_c(1)) - dd_c
!
      f(j) = cmplx( f_r, f_c, kind=CmplxKind )
   enddo
!
   end subroutine ChebyshevEval_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calChebGrid(a,b,order,r)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: order
   real (kind=RealKind), intent(in) :: a,b
   real (kind=RealKind), intent(inout) :: r(order)
!
   integer (kind=IntKind) :: i
!
   do i = 1,order
      r(i) = HALF*(b-a)*cos(PI*(i-HALF)/order)+HALF*(b+a)
   enddo
!
   end subroutine calChebGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endChebyshevSeries(cheb_st)
!  ===================================================================
   implicit none
!
   type (ChebyshevStruct), target :: cheb_st
!
   if ( associated(cheb_st%ChebCoef_R) ) deallocate( cheb_st%ChebCoef_R )
   if ( associated(cheb_st%ChebCoef_C) ) deallocate( cheb_st%ChebCoef_C )
   if ( associated(cheb_st%ChebFunc_R) ) deallocate( cheb_st%ChebFunc_R )
   if ( associated(cheb_st%ChebFunc_C) ) deallocate( cheb_st%ChebFunc_C )
!
   end subroutine endChebyshevSeries
!  ===================================================================
end module 
