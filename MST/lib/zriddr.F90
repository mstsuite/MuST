!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function zriddr(func,x1,x2,xacc,info)
!     ================================================================
!
      use KindParamModule, only : IntKind, RealKind
!
      use MathParamModule, only : zero, half, one, ten2m8, ten2m30
!
      implicit none
!
      integer (kind=IntKind), parameter :: MAXIT=200
!
      integer (kind=IntKind), intent(out) :: info
!
      real (kind=RealKind) :: zriddr
      real (kind=RealKind), intent(in) :: x1,x2,xacc
!
      integer (kind=IntKind) :: j
      real (kind=RealKind) :: fh,fl,fm,fnew,s,xh,xl,xm,xnew,df,fac
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
      info=0
      if (abs(fl).lt.ten2m30) then
         zriddr=x1
      else if (abs(fh).lt.ten2m30) then
         zriddr=x2
      else if((fl.gt.zero.and.fh.lt.zero).or.(fl.lt.zero.and.fh.gt.zero)) then
         xl=x1
         xh=x2
         fac=sqrt(-fl*fh)
         zriddr=ten2m30
         do j=1,MAXIT
            xm=half*(xl+xh)
            call func(xm,fm,df)
            s=sqrt(fm**2-fl*fh)
            if (s.lt.ten2m30) then
               return
            endif
            xnew=xm+(xm-xl)*(sign(one,fl-fh)*fm/s)
            if (abs(xnew-zriddr).le.xacc) then
               if (abs(fm/fac).gt.ten2m8) then
                  info=3
               endif
               return
            endif
            zriddr=xnew
            call func(zriddr,fnew,df)
            if (abs(fnew/fac).le.ten2m8) return
            if(sign(fm,fnew).ne.fm) then
               xl=xm
               fl=fm
               xh=zriddr
               fh=fnew
            else if(sign(fl,fnew).ne.fl) then
               xh=zriddr
               fh=fnew
            else if(sign(fh,fnew).ne.fh) then
               xl=zriddr
               fl=fnew
            else
               info = 3    !  'never get here in zriddr'
               return
            endif
            if(abs(xh-xl).le.xacc) then
               if (abs(fnew/fac).gt.ten2m8) then
                  info=3
               endif
               return
            endif
         enddo
         info=2            !  'zriddr exceed maximum iterations'
      else
         info=1            !  'root must be bracketed in zriddr'
      endif
      return
      end function zriddr
