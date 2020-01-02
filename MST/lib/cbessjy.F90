! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine bessjy(x,xnu,rj,ry,rjp,ryp)
! ====================================================================
!
  use MathParamModule, only : zero
  use MathParamModule, only : fourth
  use MathParamModule, only : one
  use MathParamModule, only : two
  use MathParamModule, only : pi
  use MathParamModule, only : tenm10
  use MathParamModule, only : tenm30
!
  implicit none
!
  integer (kind=IntKind), parameter :: MAXIT=10000
  real (kind=Realkind), parameter :: XMIN=two
  real (kind=Realkind), parameter :: EPS=ten2m10
  real (kind=Realkind), parameter :: FPMIN=ten2m30
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
