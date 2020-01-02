function rtsafe(funcd,x1,x2,acc,info)
!
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : zero, half, two
!
   implicit none
!
   real (kind=RealKind) :: rtsafe
   real (kind=RealKind), intent(in) :: x1
   real (kind=RealKind), intent(in) :: x2
   real (kind=RealKind), intent(in) :: acc
!
   real (kind=RealKind) :: df, dx, dxold, f, fh, fl, temp, xh, xl
!
   integer (kind=IntKind), intent(out) :: info
   integer (kind=IntKind), parameter :: maxit = 300
   integer (kind=IntKind) :: j
!
   interface
      subroutine funcd(x,f,df)
         use KindParamModule, only : RealKind
         real (kind=RealKind), intent(in) :: x
         real (kind=RealKind), intent(out) :: f
         real (kind=RealKind), intent(out) :: df
      end subroutine funcd
   end interface
!
   info=0
!
!  -------------------------------------------------------------------
   call funcd(x1,fl,df)
   call funcd(x2,fh,df)
!  -------------------------------------------------------------------
!
   if((fl.gt.zero .and. fh.gt.zero) .or. (fl.lt.zero .and. fh.lt.zero)) then
      info=1
      return       ! 'root is not bracketed in rtsafe'
   else if(fl.eq.zero)then
      rtsafe=x1
      return
   else if(fh.eq.zero)then
      rtsafe=x2
      return
   else if(fl.lt.zero)then
      xl=x1
      xh=x2
   else
      xh=x1
      xl=x2
   endif
   rtsafe=half*(x1+x2)
   dxold=abs(x2-x1)
   dx=dxold
!  -------------------------------------------------------------------
   call funcd(rtsafe,f,df)
!  -------------------------------------------------------------------
   do j=1,maxit
      if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.zero .or.           &
         abs(two*f).gt.abs(dxold*df) ) then
         dxold=dx
         dx=half*(xh-xl)
         rtsafe=xl+dx
      else
         dxold=dx
         dx=f/df
         temp=rtsafe
         rtsafe=rtsafe-dx
      endif
      if(abs(dx).lt.acc) then
         info=2  ! 'rtsafe is not able to find a root in the range'
         return
      endif
!     ----------------------------------------------------------------
      call funcd(rtsafe,f,df)
!     ----------------------------------------------------------------
      if(abs(f).lt.acc) then
         return
      else if(f.lt.zero) then
         xl=rtsafe
      else
         xh=rtsafe
      endif
   enddo
   info=3 
   return  ! 'rtsafe exceeding maximum iterations'
end function rtsafe
