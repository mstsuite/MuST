program tst_bessy
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use BesselModule, only : SphericalNeumann
!
   implicit none
!
   integer (kind=IntKind) :: lmax
   integer (kind=IntKind) :: l
!
   complex (kind=CmplxKind) :: x
!
   complex (kind=CmplxKind), allocatable :: nl(:), dnl(:)
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
   write(6,'('' The maximum value of l (0,1,...)? > '',$)')
   read(5,*)lmax
   write(6,'('' The Neumann function argument x?  > '',$)')
   read(5,*)x
!
   allocate ( nl(0:lmax), dnl(0:lmax), nlp(0:lmax), dnlp(0:lmax) )
!
   if (abs(x).lt.2.0) then
      do l=0,lmax
         call sphbes(l,x,nl(l),dnl(l))
      enddo
   else
      call SphericalNeumann(lmax,x,nl,dnl)
   endif
   call SphericalNeumann(lmax,x,nlp,dnlp)
!
   write(6,*)' '
   do l=0,lmax
      write(6,'(1i5,2d18.10)')l,nl(l)-nlp(l)
   enddo
   write(6,*)' '
   do l=0,lmax
      write(6,'(1i5,2d18.10)')l,dnl(l)-dnlp(l)
   enddo
   stop 'ok'
!
end program tst_bessy
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine sphbes(n,x,sy,syp)
! ====================================================================
!
  use KindParamModule, only : IntKind, RealKind, CmplxKind
  use MathParamModule, only : sqrt_pio2
  use MathParamModule, only : zero
  use MathParamModule, only : half
  use MathParamModule, only : two
  use ErrorHandlerModule,  only : ErrorHandler
!
  implicit none
!
  integer (kind=IntKind), intent(in) :: n
  complex (kind=CmplxKind), intent(in) :: x
  complex (kind=CmplxKind), intent(out) :: sy,syp
!
  real (kind=RealKind) :: order
!
  complex (kind=CmplxKind) :: factor,ry,ryp
!
  if (n.lt.0) then
     call ErrorHandler('SPHBES','bad argument in sphbes',n)
  else if (abs(x).eq.zero) then
     call ErrorHandler('SPHBES','bad argument in sphbes',x)
  endif
!
  order=n+half
! --------------------------------------------------------------------
  call bessy(x,order,ry,ryp)
! --------------------------------------------------------------------
  factor=sqrt_pio2/sqrt(x)
  sy=factor*ry
  syp=factor*ryp-sy/(two*x)
!
  end subroutine sphbes
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine bessy(x,xnu,ry,ryp)
! ====================================================================
!
  use KindParamModule, only : IntKind, RealKind, CmplxKind
  use CmplxFunctionModule, only : ComplexSign
  use MathParamModule, only : czero
  use MathParamModule, only : cone
  use MathParamModule, only : zero
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
  integer (kind=IntKind) :: i,l,nl,j
!
  real (kind=RealKind), intent(in) :: xnu
  complex (kind=Cmplxkind), intent(in) :: x
  complex (kind=Cmplxkind), intent(out) :: ry,ryp
!
  real (kind=RealKind) :: a,xmu,xmu2,gam1,gam2,gammi,gampl,          &
                          pimu,pimu2,fact,fact2,fact3,r,ri,bi
  complex (kind=Cmplxkind) :: br,c,d,del,del1,e,cfact,ff,p,q,        &
                          ry1,rymu,rymup,rytemp,sum,                 &
                          sum1,x2,xi,xi2,gam
!
  if(xnu.lt.zero) then
     call ErrorHandler('BESSY','bad arguments in bessjy',xnu)
  endif
!
  if (abs(x).lt.XMIN)then
     nl=int(xnu+half)
  else
     nl=max(0,int(xnu-abs(x)+1.5d0))
  endif
  xmu=xnu-nl
print *,'xmu = ',xmu
  xmu2=xmu*xmu
  xi=one/x
  xi2=two*xi
!
  if(abs(x).lt.XMIN) then
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
        cfact=cone
     else
        cfact=(exp(e)-exp(-e))/(e*two)
     endif
!    -----------------------------------------------------------------
     call beschb(xmu,gam1,gam2,gampl,gammi)
!    -----------------------------------------------------------------
     ff=two/pi*fact*(gam1*(exp(e)+exp(-e))/two+gam2*cfact*d)  ! => f0
     e=exp(e)
     p=e/(gampl*pi)                                            ! => p0
     q=cone/(e*pi*gammi)                                       ! => q0
     pimu2=half*pimu
     if(abs(pimu2).lt.EPS)then
        fact3=one
     else
        fact3=sin(pimu2)/pimu2
     endif
     r=pi*pimu2*fact3*fact3                        ! => sin^2(mu*pi/2) * 2/mu
     c=cone
     d=-x2*x2                                      ! => x*x/4
     sum=ff+r*q                                    ! sum  = c0*g0
     sum1=p                                        ! sum1 = c0*h0
     j=0
     do i=1,MAXIT
        ri=i                                       ! => k
        ff=(i*ff+p+q)/(ri*ri-xmu2)                 ! => fk
        c=c*d/ri                                   ! => ck
        p=p/(ri-xmu)                               ! => pk
        q=q/(ri+xmu)                               ! => qk
        del=c*(ff+r*q)                             ! => ck*gk
        sum=sum+del
        del1=c*p-i*del
        sum1=sum1+del1
        if(abs(del).lt.(one+abs(sum))*EPS) exit
        j=i
     enddo
     if (j.eq.MAXIT) then
        call ErrorHandler('BESSY','bessy series failed to converge')
     endif
     rymu=-sum                                     ! => Yu
     ry1=-sum1*xi2                                 ! => Yu+1
     rymup=xmu*xi*rymu-ry1                         ! => Y'u
  endif
!
! ====================================================================
! Starting from Yu and Y'u, using upward recursion relation to obtain 
! Yv and Y'v
! ====================================================================
  do i=1,nl
     rytemp=(xmu+i)*xi2*ry1-rymu
     rymu=ry1
     ry1=rytemp
  enddo
  ry=rymu
  ryp=xnu*xi*rymu-ry1
!
  end subroutine bessy
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
