!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function rtsec(func,x1,x2,xacc,info)
!     ================================================================
!
      use KindParamModule, only : IntKind, RealKind
!
      use MathParamModule, only : zero, ten2m8
!
      implicit none
!
      integer (kind=IntKind), parameter :: MAXIT = 100
!
      integer (kind=IntKind), intent(out) :: info
!
      real (kind=RealKind) :: rtsec
      real (kind=RealKind), intent(in) :: x1,x2,xacc
!
      integer (kind=IntKind) :: j
      real (kind=RealKind) :: dx,f,fl,swap,xl,df
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
      call func(x2,f,df)
!
      info=0
      if(abs(fl).lt.abs(f))then
         rtsec=x1
         xl=x2
         swap=fl
         fl=f
         f=swap
      else
         xl=x1
         rtsec=x2
      endif
      do j=1,MAXIT
         dx=(xl-rtsec)*f/(f-fl)
         xl=rtsec
         fl=f
         rtsec=rtsec+dx
         call func(rtsec,f,df)
         if(abs(dx).lt.xacc .or. abs(f).lt.ten2m8) return
      enddo
      info=2                 ! 'rtsec exceed maximum iterations'
      end function rtsec
