program tst_bessjy
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use BesselModule, only : SphericalBessel, SphericalNeumann
!
   implicit none
!
   integer (kind=IntKind) :: lmax
   integer (kind=IntKind) :: l
!
   real (kind=RealKind) :: x
!
   complex (kind=CmplxKind) :: xc
!
   real (kind=RealKind), allocatable :: jl(:), djl(:)
   real (kind=RealKind), allocatable :: nl(:), dnl(:)
!
   complex (kind=CmplxKind), allocatable :: jlp(:), djlp(:)
   complex (kind=CmplxKind), allocatable :: nlp(:), dnlp(:)
!
   write(6,*)' '
   write(6,*)'          ***************************************************'
   write(6,*)'          *                                                 *'
   write(6,*)'          *    TEST CODE TO CHECK THE WRONSKIAN RELATION    *'
   write(6,*)'          *                                                 *'
   write(6,*)'          *  OF SPHERICAL BESSEL FUNCTIONS: j_l(x), n_l(x)  *'
   write(6,*)'          *                                                 *'
   write(6,*)'          ***************************************************'
   write(6,*)' '
   write(6,'('' The maximum value of l (0, 1, 2,...)?  > '',$)')
   read(5,*)lmax
   write(6,'('' The Bessel function real argument x?   > '',$)')
   read(5,*)x
!
   allocate ( jl(0:lmax), djl(0:lmax), nl(0:lmax), dnl(0:lmax),  &
              jlp(0:lmax), djlp(0:lmax), nlp(0:lmax), dnlp(0:lmax) )
!
   do l=0,lmax
      call sphbes(l,x,jl(l),nl(l),djl(l),dnl(l))
   enddo
   xc=x
   call SphericalBessel(lmax,xc,jlp,djlp)
   call SphericalNeumann(lmax,xc,nlp,dnlp)
!
!  write(6,*)' '
!  do l=0,lmax
!     write(6,'(1i5,3d18.11)')l,x*x*(jl(l)*dnl(l)-djl(l)*nl(l)),      &
!                               xc*xc*(jlp(l)*dnlp(l)-djlp(l)*nlp(l))
!  enddo
!
   write(6,*)' '
   do l=0,lmax
!     write(6,'(1i5,4d18.11)')l,jl(l)-jlp(l),nl(l)-nlp(l)
      write(6,'(1i5,2d24.16)')l,jl(l)
   enddo
   write(6,*)' '
   do l=0,lmax
!     write(6,'(1i5,4d18.11)')l,djl(l)-djlp(l),dnl(l)-dnlp(l)
      write(6,'(1i5,2d24.16)')l,djl(l)
   enddo
   stop 'ok'
!
end program tst_bessjy
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine sphbes(n,x,sj,sy,sjp,syp)
! ====================================================================
!
  use KindParamModule, only : IntKind, RealKind, CmplxKind
  use MathParamModule, only : sqrt_pio2
  use MathParamModule, only : half
  use MathParamModule, only : two
  use ErrorHandlerModule,  only : ErrorHandler
!
  implicit none
!
  integer (kind=IntKind), intent(in) :: n
  real (kind=RealKind), intent(in) :: x
  real (kind=RealKind), intent(out) :: sj,sjp,sy,syp
!
  real (kind=RealKind) :: factor,order,rj,rjp,ry,ryp
!
  if (n.lt.0) then
     call ErrorHandler('SPHBES','bad argument in sphbes',n)
  else if (x.le.0.) then
     call ErrorHandler('SPHBES','bad argument in sphbes',x)
  endif
!
  order=n+half
! --------------------------------------------------------------------
  call bessjy(x,order,rj,ry,rjp,ryp)
! --------------------------------------------------------------------
  factor=sqrt_pio2/sqrt(x)
  sj=factor*rj
  sy=factor*ry
  sjp=factor*rjp-sj/(two*x)
  syp=factor*ryp-sy/(two*x)
!
  end subroutine sphbes
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine bessjy(x,xnu,rj,ry,rjp,ryp)
! ====================================================================
!
  use KindParamModule, only : IntKind, RealKind, CmplxKind
  use MathParamModule, only : zero
! use MathParamModule, only : ten2m16, ten2m30
  use MathParamModule, only : fourth
  use MathParamModule, only : half
  use MathParamModule, only : one
  use MathParamModule, only : two
  use MathParamModule, only : pi
  use ErrorHandlerModule,  only : ErrorHandler
!
  implicit none
!
  integer (kind=IntKind), parameter :: MAXIT=10000
  real (kind=Realkind), parameter :: XMIN=two
  real (kind=Realkind), parameter :: EPS=1.0d-16
  real (kind=Realkind), parameter :: FPMIN=1.0d-128
!
  integer (kind=IntKind) :: i,isign,l,nl
!
  real (kind=Realkind), intent(out) :: rj,rjp,ry,ryp,x,xnu
  real (kind=Realkind) :: a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli, &
                          dr,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,    &
                          gammi,gampl,h,p,pimu,pimu2,q,r,rjl,rjl1,rjmu,&
                          rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,  &
                          sum1,temp,w,x2,xi,xi2,xmu,xmu2
!
  if (x.le.zero) then
     call ErrorHandler('BESSJY','bad arguments in bessjy',x)
  else if(xnu.lt.zero) then
     call ErrorHandler('BESSJY','bad arguments in bessjy',xnu)
  endif
!
  if (x.lt.XMIN)then
     nl=int(xnu+half)
  else
     nl=max(0,int(xnu-x+1.5d0))
  endif
  xmu=xnu-nl
  xmu2=xmu*xmu
  xi=one/x
  xi2=two*xi
  w=xi2/pi
  isign=1
  h=xnu*xi
  if (h.lt.FPMIN) h=FPMIN
  b=xi2*xnu
  d=zero
  c=h
  do i=1,MAXIT
     b=b+xi2
     d=b-d
     if(abs(d).lt.FPMIN)d=FPMIN
     c=b-one/c
     if(abs(c).lt.FPMIN)c=FPMIN
     d=one/d
     del=c*d
     h=del*h
     if(d.lt.zero)isign=-isign
     if(abs(del-one).lt.EPS)goto 1
  enddo
  call ErrorHandler('BESSJY','x too large in bessjy; try asymptotic expansion')
1 continue
  rjl=isign*FPMIN
  rjpl=h*rjl
  rjl1=rjl
  rjp1=rjpl
  fact=xnu*xi
  do l=nl,1,-1
     rjtemp=fact*rjl+rjpl
     fact=fact-xi
     rjpl=fact*rjtemp-rjl
     rjl=rjtemp
  enddo
  if(rjl.eq.zero)rjl=EPS
  f=rjpl/rjl
  if(x.lt.XMIN) then
     x2=half*x
     pimu=pi*xmu
     if(abs(pimu).lt.EPS)then
        fact=one
     else
        fact=pimu/sin(pimu)
     endif
     d=-log(x2)
     e=xmu*d
     if(abs(e).lt.EPS)then
        fact2=one
     else
        fact2=sinh(e)/e
     endif
     call beschb(xmu,gam1,gam2,gampl,gammi)
     ff=two/pi*fact*(gam1*cosh(e)+gam2*fact2*d)
     e=exp(e)
     p=e/(gampl*pi)
     q=one/(e*pi*gammi)
     pimu2=half*pimu
     if(abs(pimu2).lt.EPS)then
        fact3=one
     else
        fact3=sin(pimu2)/pimu2
     endif
     r=pi*pimu2*fact3*fact3
     c=one
     d=-x2*x2
     sum=ff+r*q
     sum1=p
     do i=1,MAXIT
        ff=(i*ff+p+q)/(i*i-xmu2)
        c=c*d/i
        p=p/(i-xmu)
        q=q/(i+xmu)
        del=c*(ff+r*q)
        sum=sum+del
        del1=c*p-i*del
        sum1=sum1+del1
        if(abs(del).lt.(one+abs(sum))*EPS)goto 2
     enddo
     call ErrorHandler('BESSJY','bessy series failed to converge')
2    continue
     rymu=-sum
     ry1=-sum1*xi2
     rymup=xmu*xi*rymu-ry1
     rjmu=w/(rymup-f*rymu)
  else
     a=fourth-xmu2
     p=-half*xi
     q=one
     br=two*x
     bi=two
     fact=a*xi/(p*p+q*q)
     cr=br+q*fact
     ci=bi+p*fact
     den=br*br+bi*bi
     dr=br/den
     di=-bi/den
     dlr=cr*dr-ci*di
     dli=cr*di+ci*dr
     temp=p*dlr-q*dli
     q=p*dli+q*dlr
     p=temp
     do i=2,MAXIT
        a=a+2*(i-1)
        bi=bi+two
        dr=a*dr+br
        di=a*di+bi
        if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
        fact=a/(cr*cr+ci*ci)
        cr=br+cr*fact
        ci=bi-ci*fact
        if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
        den=dr*dr+di*di
        dr=dr/den
        di=-di/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
     enddo
     call ErrorHandler('BESSJY','cf2 failed in bessjy')
3    continue
     gam=(p-f)/q
     rjmu=sqrt(w/((p-f)*gam+q))
     rjmu=sign(rjmu,rjl)
     rymu=rjmu*gam
     rymup=rymu*(p+q/gam)
     ry1=xmu*xi*rymu-rymup
  endif
  fact=rjmu/rjl
  rj=rjl1*fact
  rjp=rjp1*fact
  do i=1,nl
     rytemp=(xmu+i)*xi2*ry1-rymu
     rymu=ry1
     ry1=rytemp
  enddo
  ry=rymu
  ryp=xnu*xi*rymu-ry1
!
  end subroutine bessjy
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine beschb(x,gam1,gam2,gampl,gammi)
! ====================================================================
!
  use KindParamModule, only : IntKind, RealKind, CmplxKind
  use MathParamModule, only : one
  use MathParamModule, only : eight
!
  implicit none
!
  integer (kind=IntKind), parameter :: NUSE1=7
  integer (kind=IntKind), parameter :: NUSE2=8
!
  real (kind=RealKind), intent(in) :: x
  real (kind=RealKind), intent(out) :: gam1,gam2,gammi,gampl
!
  real (kind=RealKind) :: xx,c1(7),c2(8),chebev
  SAVE c1,c2
  data c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,    &
          -3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
  data c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,  &
          -4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
  xx=eight*x*x-one
  gam1=chebev(-one,one,c1,NUSE1,xx)
  gam2=chebev(-one,one,c2,NUSE2,xx)
  gampl=gam2-x*gam1
  gammi=gam2+x*gam1
!
  end subroutine beschb
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function chebev(a,b,c,m,x)
! ====================================================================
!
  use KindParamModule, only : IntKind, RealKind, CmplxKind
  use MathParamModule, only : zero
  use MathParamModule, only : half
  use MathParamModule, only : two
  use ErrorHandlerModule,  only : ErrorHandler
!
  implicit none
!
  integer (kind=IntKind), intent(in) :: m
  integer :: j
!
  real (kind=RealKind) :: chebev
  real (kind=RealKind), intent(in) :: a,b,x,c(m)
  real (kind=RealKind) :: d,dd,sv,y,y2
!
  if ((x-a)*(x-b).gt.zero) then
     call ErrorHandler('CHEBEV','x not in range in chebev',x)
  endif
!
  d=zero
  dd=zero
  y=(two*x-a-b)/(b-a)
  y2=two*y
  do j=m,2,-1
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
  enddo
  chebev=y*d-dd+half*c(1)
  end function chebev
