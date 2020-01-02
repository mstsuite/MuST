!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function zbrent(func,x1,x2,tol,info)
!     ================================================================
!
      use KindParamModule, only : IntKind, RealKind
!
      use MathParamModule, only : zero, half, one, two, three, ten2m8
!
      implicit   none
!
      integer (kind=IntKind), parameter :: itmax = 400
      integer (kind=IntKind) :: info
      integer (kind=IntKind) :: iter
!
      real (kind=RealKind) :: zbrent
      real (kind=RealKind), intent(in) :: tol,x1,x2
!
      real (kind=RealKind), parameter :: eps = ten2m8
      real (kind=RealKind) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,df
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
      info=0
!
      a=x1
      b=x2
!     ----------------------------------------------------------------
      call func(a,fa,df)
      call func(b,fb,df)
!     ----------------------------------------------------------------
!
      if((fa.gt.zero.and.fb.gt.zero).or.(fa.lt.zero.and.fb.lt.zero)) then
         info=1
         return   ! 'root must be bracketed for zbrent'
      endif
!
      c=b
      fc=fb
      do iter=1,itmax
         if((fb.gt.zero.and.fc.gt.zero).or.(fb.lt.zero.and.fc.lt.zero)) then
            c=a
            fc=fa
            d=b-a
            e=d
         endif
         if(abs(fc).lt.abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         endif
         tol1=two*eps*abs(b)+half*tol
         xm=half*(c-b)
         if(abs(xm).le.tol1 .or. fb.eq.zero)then
!           -----------------------------------------------------------
            call func(b,fb,df)
!           -----------------------------------------------------------
            if(abs(fb).lt.tol) then
               zbrent=b
               return
            endif
         endif
         if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
            s=fb/fa
            if(a.eq.c) then
               p=two*xm*s
               q=one-s
            else
               q=fa/fc
               r=fb/fc
               p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
               q=(q-one)*(r-one)*(s-one)
            endif
            if(p.gt.zero) q=-q
            p=abs(p)
            if(two*p.lt.min(three*xm*q-abs(tol1*q),abs(e*q))) then
               e=d
               d=p/q
            else
               d=xm
               e=d
            endif
         else
            d=xm
            e=d
         endif
         a=b
         fa=fb
         if(abs(d) .gt. tol1) then
            b=b+d
         else
            b=b+sign(tol1,xm)
         endif
!        -------------------------------------------------------------
         call func(b,fb,df)
!        -------------------------------------------------------------
      enddo
      info=2 ! 'zbrent exceeding maximum iterations'
      zbrent=b
!
      return
      end
