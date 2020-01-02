!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function rtflsp(func,x1,x2,xacc,info)
!     ================================================================
!
      use KindParamModule, only : IntKind, RealKind
!
      use MathParamModule, only : zero, ten2m8
!
      implicit none
!
      integer (kind=IntKind), parameter :: MAXIT=100
!
      integer (kind=IntKind), intent(out) :: info
!
      real (kind=RealKind) :: rtflsp
      real (kind=RealKind), intent (in) :: x1,x2,xacc
!
      integer j
      real (kind=RealKind) :: del,dx,f,fh,fl,swap,xh,xl,df
!
      interface
         subroutine func(x,f,df)
            use KindParamModule, only : RealKind
            real (kind=RealKind), intent(in) :: x
            real (kind=RealKind), intent(out) :: f
            real (kind=RealKind), intent(out) :: df
         end subroutine func
      end interface
!
      call func(x1,fl,df)
      call func(x2,fh,df)
!
      if (fl*fh.gt.zero) then
         info=1              !  'root must be bracketed in rtflsp'
         return
      endif
!
      info=0
      if(fl.lt.zero)then
         xl=x1
         xh=x2
      else
         xl=x2
         xh=x1
         swap=fl
         fl=fh
         fh=swap
      endif
      dx=xh-xl
      do j=1,MAXIT
         rtflsp=xl+dx*fl/(fl-fh)
         call func(rtflsp,f,df)
         if(f.lt.zero) then
            del=xl-rtflsp
            xl=rtflsp
            fl=f
         else
            del=xh-rtflsp
            xh=rtflsp
            fh=f
         endif
         dx=xh-xl
         if(abs(del).lt.xacc .or. abs(f).lt.ten2m8)return
      enddo
      info=2                 !  'rtflsp exceed maximum iterations'
      end function rtflsp
