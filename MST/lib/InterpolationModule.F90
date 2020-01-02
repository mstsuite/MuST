module InterpolationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, HALF, ONE, FOUR, TEN2m14, TEN2p14
   use ErrorHandlerModule, only : ErrorHandler
!
public :: PolyInterp,       &
          getInterpolation, &
          FitInterp,        &
          LeastSqFitInterp
!
   interface getInterpolation
      module procedure getInterp_r, getInterp_c
   end interface
!
   interface PolyInterp
      module procedure polint_r, polint_c
   end interface
!
   interface FitInterp
      module procedure interp0, cinterp0, &
                       interp1, cinterp1
   end interface
!
   interface LeastSqFitInterp
      module procedure leastSqFit_interp, leastSqFit_cinterp, &
                       leastSqFit_interp_v, leastSqFit_cinterp_v
   end interface
!
private
!
   integer (kind=IntKind) :: j_inter, n_inter
   integer (kind=IntKind), parameter :: np_inter = 5
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterp_r(ns,xa,ya,x,err) result(y)
!  ===================================================================
   implicit none
!  *******************************************************************
!  Given arrays xa, ya of length n, and given a value x, this routine
!  returns an interpolation value y and an error estimate dy.
!  Taken from numerical recipes; Modified by Yang Wang
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: ns
!
   real (kind=RealKind), intent(in) :: x, xa(ns), ya(ns)
   real (kind=RealKind), intent(out) :: err
   real (kind=RealKind) :: y
!
   n_inter = min(ns,np_inter)
!  ns = size(xa)
!  -------------------------------------------------------------------
   call huntx(ns,xa,x,j_inter)
!  -------------------------------------------------------------------
!  j_inter=min(j_inter-(n_inter-1)/2,ns-n_inter+1)
   if (j_inter > ns-(n_inter-1)/2) then
      j_inter=ns-n_inter+1
   else if (2*j_inter+1 > n_inter) then
      j_inter=j_inter-(n_inter-1)/2
   else
      j_inter=1
   endif
!
!  -------------------------------------------------------------------
   call polint_r(n_inter,xa(j_inter:j_inter+n_inter-1),               &
                 ya(j_inter:j_inter+n_inter-1), x, y, err)
!  -------------------------------------------------------------------
!
   end function getInterp_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getInterp_c(ns,xa,ya,x,err) result(y)
!  ===================================================================
   implicit none
!  *******************************************************************
!  Given arrays xa, ya of length n, and given a value x, this routine
!  returns an interpolation value y and an error estimate dy.
!  Taken from numerical recipes; Modified by Yang Wang
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: ns
!
   real (kind=RealKind), intent(in) :: x, xa(ns)
   real (kind=RealKind), intent(out) :: err
!
   complex (kind=CmplxKind), intent(in) :: ya(ns)
   complex (kind=CmplxKind) :: y
!
   n_inter = min(ns,np_inter)
!  ns = size(xa)
!  -------------------------------------------------------------------
   call huntx(ns,xa,x,j_inter)
!  -------------------------------------------------------------------
!  j_inter=min(j_inter-(n_inter-1)/2,ns-n_inter+1)
   if (j_inter > ns-(n_inter-1)/2) then
      j_inter=ns-n_inter+1
   else if (2*j_inter+1 > n_inter) then
      j_inter=j_inter-(n_inter-1)/2
   else
      j_inter=1
   endif
!
!  -------------------------------------------------------------------
   call polint_c(n_inter,xa(j_inter:j_inter+n_inter-1),               &
                 ya(j_inter:j_inter+n_inter-1), x, y, err)
!  -------------------------------------------------------------------
!
   end function getInterp_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine polint_r(n,xa,ya,x,y,dy)
!  ===================================================================
   implicit none
!  *******************************************************************
!  Given arrays xa, ya of length n, and given a value x, this routine
!  returns an interpolation value y and an error estimate dy.
!  Taken from numerical recipes; Modified by Yang Wang
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i, m, ns
   integer (kind=IntKind), parameter :: nmax=10

   real (kind=RealKind), intent(in) :: x, xa(n), ya(n)
   real (kind=RealKind), intent(out) :: dy, y
   real (kind=RealKind) :: dif, dift, ho, hp
   real (kind=RealKind) :: den, w, c(nmax), d(nmax)
! 
!  n=size(xa)
!  if (n /= size(ya)) then
!     call ErrorHandler('polint_r','size(x) <> size(y)',n,size(ya))
   if (n > nmax) then
      call ErrorHandler('polint_r','nmax < n',nmax,n)
   else if (n < 2) then
      call ErrorHandler('polint_r','n < 2',n)
   endif
 
   ns=1
   dif=abs(x-xa(1))
   do i=1,n
      dift=abs(x-xa(i))
      if (dift.lt.dif) then
         ns=i
         dif=dift
      endif
      c(i)=ya(i)
      d(i)=ya(i)
   enddo

   y=ya(ns)
   ns=ns-1
   do m=1,n-1
      do i=1,n-m
         ho=xa(i)-x
         hp=xa(i+m)-x
         w=c(i+1)-d(i)
         den=ho-hp
         if (den.eq.ZERO) then
            call ErrorHandler('POLINT_R','den = 0')
         endif
         den=w/den
         d(i)=hp*den
         c(i)=ho*den
      enddo
      if (2*ns.lt.n-m)then
         dy=c(ns+1)
      else
         dy=d(ns)
         ns=ns-1
      endif
      y=y+dy
   enddo
   end subroutine polint_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine polint_c(n,xa,ya,x,y,err)
!  ===================================================================
   implicit none
!  *******************************************************************
!  Given arrays xa, ya of length n, and given a value x, this routine
!  returns an interpolation value y and an error estimate dy.
!  Taken from numerical recipes; Modified by Yang Wang
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i, m, ns
   integer (kind=IntKind), parameter :: nmax=10

   real (kind=RealKind), intent(in) :: x, xa(n)
   real (kind=RealKind) :: dif, dift, ho, hp
   real (kind=RealKind), intent(out) :: err
!
   complex (kind=CmplxKind), intent(in) :: ya(n)
   complex (kind=CmplxKind), intent(out) :: y
   complex (kind=CmplxKind) :: dy, den, w, c(nmax), d(nmax)
!
!  n=size(xa)
!  if (n /= size(ya)) then
!     call ErrorHandler('polint_c','size(x) <> size(y)',n,size(ya))
   if (nmax < n) then
      call ErrorHandler('polint_c','nmax < n',nmax,n)
   else if (n < 2) then
      call ErrorHandler('polint_c','n < 2',n)
   endif

   ns=1
   dif=abs(x-xa(1))
   do i=1,n
      dift=abs(x-xa(i))
      if (dift.lt.dif) then
         ns=i
         dif=dift
      endif
      c(i)=ya(i)
      d(i)=ya(i)
   enddo

   y=ya(ns)
   ns=ns-1
   do m=1,n-1
      do i=1,n-m
         ho=xa(i)-x
         hp=xa(i+m)-x
         w=c(i+1)-d(i)
         den=ho-hp
         if(abs(den).eq.ZERO) then
            call ErrorHandler('POLINT_C','den = 0')
         endif
         den=w/den
         d(i)=hp*den
         c(i)=ho*den
      enddo
      if (2*ns.lt.n-m) then
         dy=c(ns+1)
      else
         dy=d(ns)
         ns=ns-1
      endif
      y=y+dy
   enddo
   err=abs(dy)
   end subroutine polint_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine huntx(n,xx,x,jlo)
   implicit none
!
   integer (kind=IntKind), intent(inout) :: jlo
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: inc, jhi, jm
!
   real (kind=RealKind), intent(in) :: xx(:)
   real (kind=RealKind), intent(in) :: x
!
   logical :: ascnd
! 
!  *******************************************************************
!  Taken from numerical recipes. Modified by Yang Wang
!  Modified according to the Fortran 90 version
!
!  Given an array xx(1:n), and given a value x, returns a value jlo
!  such that x is between xx(jlo) and xx(jlo+1). xx must be monotonic,
!  either increasing or decreasing. jlo=0 or jlo=n is returned to
!  indicate that x is out of range. jlo on input is taken as the initial
!  guess for jlo on output.
!  *******************************************************************
!  n=size(xx,1)
   ascnd=(xx(n) >= xx(1))
   if (jlo <= 0 .or. jlo > n) then
      jlo=0
      jhi=n+1
   else
      inc=1
      if (x >= xx(jlo) .eqv. ascnd) then
         do
            jhi=jlo+inc
            if(jhi > n)then
               jhi=n+1
               exit
            else if (x < xx(jhi) .eqv. ascnd) then
               exit
            else 
               jlo=jhi
               inc=inc+inc
            endif
         enddo
      else
         jhi=jlo
         do
            jlo=jhi-inc
            if (jlo < 1) then
               jlo=0
               exit
            else if (x >= xx(jlo) .eqv. ascnd) then
               exit
            else
               jhi=jlo
               inc=inc+inc
            endif
         enddo
      endif
   endif
   do
      if (jhi-jlo <= 1) then
         if (x == xx(1)) then
            jlo=1
         else if (x == xx(n)) then
            jlo=n-1
         endif
         exit
      else
         jm=(jhi+jlo)/2
         if (x >= xx(jm) .eqv. ascnd) then
            jlo=jm
         else
            jhi=jm
         endif
      endif
   enddo
   end subroutine huntx
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine interp0(nr,r,f,rs,ps,dps)
!  ===================================================================  
!
!  --- this routine interpolates functn f on grid r
!      returns value ps and its derivative dps at rs 
!
!  modified by yw, 11/20/03
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind) :: ip1, i
!
   real (kind=RealKind), intent(in) :: r(nr),f(nr), rs
   real (kind=RealKind), intent(out) :: ps, dps
   real (kind=RealKind) :: c(4), h1, h2, c2
!
   ip1=2
   do i=2,nr-1
      if (rs > r(i)) then
         ip1=i+1
      endif
   enddo
   i=ip1-1
!  -------------------------------------------------------------------
   call fit(nr,r,f,i,c)
!  -------------------------------------------------------------------
   h1=rs-r(i)
   h2=r(ip1)-rs
   c2=ONE+c(4)*h1*h2
   ps=((c(3)-c(2))*h1+c(1))/c2
   dps=c(2)+(c(3)-c(2)+ps*c(4)*(h1-h2))/c2
   ps=ps+f(i)-c(1)+c(2)*h1
!
   end subroutine interp0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine interp1(nr,r,f,nrs,rs,ps,dps)
!  ===================================================================  
!
!  --- this routine interpolates functn f on grid r
!      returns value ps and its derivative dps at rs 
!
!  modified by yw, 11/20/03
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nrs
   integer (kind=IntKind) :: ip1, i, ir
!
   real (kind=RealKind), intent(in) :: r(nr),f(nr), rs(nrs)
   real (kind=RealKind), intent(out) :: ps(nrs), dps(nrs)
   real (kind=RealKind) :: c(4), h1, h2, c2
!
   do ir = 1,nrs
      ip1=2
      do i=2,nr-1
         if (rs(ir) > r(i)) then
            ip1=i+1
         endif
      enddo
      i=ip1-1
!     ----------------------------------------------------------------
      call fit(nr,r,f,i,c)
!     ----------------------------------------------------------------
      h1=rs(ir)-r(i)
      h2=r(ip1)-rs(ir)
      c2=ONE+c(4)*h1*h2
      ps(ir)=((c(3)-c(2))*h1+c(1))/c2
      dps(ir)=c(2)+(c(3)-c(2)+ps(ir)*c(4)*(h1-h2))/c2
      ps(ir)=ps(ir)+f(i)-c(1)+c(2)*h1
   enddo
!
   end subroutine interp1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cinterp0(nr,r,f,rs,ps,dps)
!  ===================================================================  
!
!  --- this routine interpolates functn f on grid r
!      returns value ps and its derivative dps at rs 
!
!  modified by yw, 11/20/03
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind) :: ip1, i
!
   real (kind=RealKind), intent(in) :: r(nr), rs
   real (kind=RealKind) :: c(4), h1, h2, c2, psr, psi, dpsr, dpsi, work(nr)
!
   complex (kind=CmplxKind), intent(in) :: f(nr)
   complex (kind=CmplxKind), intent(out) :: ps, dps
!
   ip1=2
   do i=2,nr-1
      if (rs > r(i)) then
         ip1=i+1
      endif
   enddo

!  Real part
   do i=1,nr
     work(i)=real(f(i),RealKind)
   enddo
   i=ip1-1
!  -------------------------------------------------------------------
   call fit(nr,r,work,i,c)
!  -------------------------------------------------------------------
   h1=rs-r(i)
   h2=r(ip1)-rs
   c2=ONE+c(4)*h1*h2
   psr=((c(3)-c(2))*h1+c(1))/c2
   dpsr=c(2)+(c(3)-c(2)+psr*c(4)*(h1-h2))/c2
   psr=psr+real(f(i),RealKind)-c(1)+c(2)*h1

!  Imaginary part
   do i=1,nr
      work(i)=aimag(f(i))
   enddo
   i=ip1-1
!  -------------------------------------------------------------------
   call fit(nr,r,work,i,c)
!  -------------------------------------------------------------------
   h1=rs-r(i)
   h2=r(ip1)-rs
   c2=ONE+c(4)*h1*h2
   psi=((c(3)-c(2))*h1+c(1))/c2
   dpsi=c(2)+(c(3)-c(2)+psi*c(4)*(h1-h2))/c2
   psi=psi+aimag(f(i))-c(1)+c(2)*h1

   ps=cmplx(psr,psi,CmplxKind)
   dps=cmplx(dpsr,dpsi,CmplxKind)
!
   end subroutine cinterp0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cinterp1(nr,r,f,nrs,rs,ps,dps)
!  ===================================================================  
!
!  --- this routine interpolates functn f on grid r
!      returns value ps and its derivative dps at rs 
!
!  modified by yw, 11/20/03
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nrs
   integer (kind=IntKind) :: ip1, i, ir
!
   real (kind=RealKind), intent(in) :: r(nr), rs(nrs)
   real (kind=RealKind) :: c(4), h1, h2, c2, psr, psi, dpsr, dpsi, work(nr)
!
   complex (kind=CmplxKind), intent(in) :: f(nr)
   complex (kind=CmplxKind), intent(out) :: ps(nrs), dps(nrs)
!
   do ir = 1,nrs
      ip1=2
      do i=2,nr-1
         if (rs(ir) > r(i)) then
            ip1=i+1
         endif
      enddo

!     Real part
      do i=1,nr
         work(i)=real(f(i),RealKind)
      enddo
      i=ip1-1
!     ----------------------------------------------------------------
      call fit(nr,r,work,i,c)
!     ----------------------------------------------------------------
      h1=rs(ir)-r(i)
      h2=r(ip1)-rs(ir)
      c2=ONE+c(4)*h1*h2
      psr=((c(3)-c(2))*h1+c(1))/c2
      dpsr=c(2)+(c(3)-c(2)+psr*c(4)*(h1-h2))/c2
      psr=psr+real(f(i),RealKind)-c(1)+c(2)*h1

!     Imaginary part
      do i=1,nr
         work(i)=aimag(f(i))
      enddo
      i=ip1-1
!     ----------------------------------------------------------------
      call fit(nr,r,work,i,c)
!     ----------------------------------------------------------------
      h1=rs(ir)-r(i)
      h2=r(ip1)-rs(ir)
      c2=ONE+c(4)*h1*h2
      psi=((c(3)-c(2))*h1+c(1))/c2
      dpsi=c(2)+(c(3)-c(2)+psi*c(4)*(h1-h2))/c2
      psi=psi+aimag(f(i))-c(1)+c(2)*h1

      ps(ir)=cmplx(psr,psi,CmplxKind)
      dps(ir)=cmplx(dpsr,dpsi,CmplxKind)
   enddo
!
   end subroutine cinterp1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine fit(nr,r,f,i,c)
!  ===================================================================
   implicit none
!  *******************************************************************
!  This routine does a four-point interpolation from i-1 to i+2
!  returns c array 4 of real which are coeficients for
!
!  f(i)-c1+c2*(r-r(i))+((c3-c2)*(r-r(i))+c1)/(1+c4*(r-r(i))*(r(i+1)-r))
!
!  It's obvious that
!
!  c3=(f(i+1)-f(i))/(r(i+1)-r(i))
!
!  c2 is adjustable and is chosen so to avoid possible singularities
!  in eqns solving for c1 and c4.
!
!  c2=((f(i+2)-f(i+1))/(r(i+2)-r(i+1))-(f(i)-f(i-1))/(r(i)-r(i-1)))
!
!  and the sign is chosen so that it is opposite of
!  (f(i+1)-f(i-1))/(r(i+1)-r(i-1))
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: nr, i
!
   real (kind=RealKind), intent(in) :: r(nr)
!
   real (kind=RealKind), intent(in) :: f(nr)
   real (kind=RealKind), intent(out) :: c(4)
!
   integer (kind=IntKind) :: ip1, i1, i2
!
   real (kind=RealKind) :: dr, h1, h2, rj
!
   real (kind=RealKind) :: drv, c0, c1, c2, c3, c4, eqn12, eqn22, gj, fj
!
   if(i < 1) then
      call ErrorHandler('FIT','i < 1',i)
   endif
!
   ip1=i+1
   if(ip1 > nr) then
      if(r(1) > r(nr)) then   ! wrap around
         ip1=1
      else
         call ErrorHandler('FIT','i > n-1',i)
      endif
   endif
!
   dr=r(ip1)-r(i)
   if(dr == ZERO) then
      call ErrorHandler('FIT','degenerate grid')
   endif
   drv=f(ip1)-f(i)
   c(3)=drv/dr
!
!  ===================================================================
!  two point interpolation
!  ===================================================================
   if(nr == 2) then
      c(1)=f(i)
      c(2)=ZERO
      c(4)=ZERO
      return
   endif
!
!  ===================================================================
!  fit the function to the two outside points
!  ===================================================================
   if(i > 1) then
      i1=i-1
   else
      i1=min(4,nr)
   endif
   if(ip1 < nr) then
      i2=ip1+1
   else
      i2=max(nr-3,1)
   endif
   if(i1 /= i2) then
      c(2)=(f(i2)-f(ip1))/(r(i2)-r(ip1))-(f(i1)-f(i))/(r(i1)-r(i))
      if(c(2)*(f(i2)-f(i1))*(r(i2)-r(i1)) > ZERO) then
         c(2)=-c(2)
      endif
!     ================================================================
!     solve the 2x2 equation for c(1) and c(4)
!     ================================================================
      h1=r(i1)-r(i)
      c0=f(i1)-f(i)
      c1=c0-c(3)*h1
      c2=c0-c(2)*h1
      h1=h1*(r(ip1)-r(i1))
      eqn12=-h1*c2
      h2=r(i2)-r(i)
      c0=f(i2)-f(i)
      c3=c0-c(3)*h2
      c4=c0-c(2)*h2
      h2=h2*(r(ip1)-r(i2))
      eqn22=-h2*c4
      gj=(c4-c2)*h1*h2
      if(gj == ZERO) then
         c(4)=ZERO
      else
         c(1)=(c1*eqn22-c3*eqn12)/gj
         c(4)=(c1*h2-c3*h1)/gj
      endif
!     ================================================================
!     If the points are on a line or nearly a line then use linear
!     interpolation
!     check this by checking to see if the denominator 1+c4*(r-r(i))*(r(i+1)-r)
!     can be zero or negative. The minimum is 1-c4*dr*dr/4
!     ================================================================
      gj=c(4)*dr*dr
      if(gj > -FOUR .and. abs(gj) > TEN2m14 .and. gj < TEN2p14) then
         c(1)=c(1)/c(4)
      else
         c(1)=f(i)
         c(4)=ZERO
      endif
   else  ! i1=i2
!     ================================================================
!     set A4=0.5*(f(i)+f(i+1)) and do a 3 point fit
!     ================================================================
      fj=HALF*(f(i)+f(ip1))
      c(1)=f(i)-fj
      c(2)=ZERO
      rj=r(i)+HALF*dr
      if(abs(r(i1)-rj) > abs(r(i2)-rj)) then
         i1=i2
      endif
      h1=r(i1)-r(i)
      c1=f(i1)-f(i)-c(3)*h1
      h1=h1*(r(ip1)-r(i1))
      if(f(i1) /= fj) then
         c(4)=-c1/((f(i1)-fj)*h1)
         gj=c(4)*dr*dr
         if(gj <= -FOUR .or. abs(gj) <= TEN2m14 .or. gj >= TEN2p14) then
            c(1)=f(i)
            c(4)=ZERO
         endif
      else
         c(1)=f(i)
         c(4)=ZERO
      endif
   endif  ! f(i1)=f(i2)
!
   end subroutine fit
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine leastSqFit_interp(nr,r,f,rs,ps,m)
!  ===================================================================
!
!  --- this routine interpolates functn f on grid r using a least square fit
!      returns value ps and its derivative dps at rs 
!
!  modified by yw, 11/20/03
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), optional :: m
   integer (kind=IntKind) :: m_poly, i
!
   real (kind=RealKind), intent(in) :: r(nr), f(nr), rs
   real (kind=RealKind), intent(out) :: ps
   real (kind=RealKind), allocatable :: c(:)
!
   if ( present(m) ) then
      m_poly = m
   else
      m_poly = 2
   endif
   allocate( c(m_poly+1) )
!  -------------------------------------------------------------------
   call leastSqFit_r(nr,r,f,m_poly,c)
!  -------------------------------------------------------------------
   ps = c(1)
   do i = 1,m_poly
      ps = ps + (rs**i)*c(i+1)
   enddo
   deallocate(c)
!
   end subroutine leastSqFit_interp
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine leastSqFit_cinterp(nr,r,f,rs,ps,m)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), optional :: m
   integer (kind=IntKind) :: m_poly, i
!
   real (kind=RealKind), intent(in) :: r(nr), rs
   complex (kind=CmplxKind), intent(in)  :: f(nr)
   complex (kind=CmplxKind), intent(out) :: ps
   complex (kind=CmplxKind), allocatable :: c(:)
!
   if ( present(m) ) then
      m_poly = m
   else
      m_poly = 3
   endif
   allocate( c(m_poly+1) )
!  -------------------------------------------------------------------
   call leastSqFit_c(nr,r,f,m_poly,c)
!  -------------------------------------------------------------------
   ps = c(1)
   do i = 1,m_poly
      ps = ps + (rs**i)*c(i+1)
   enddo
   deallocate(c)
!
   end subroutine leastSqFit_cinterp
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine leastSqFit_interp_v(nr,r,f,nrs,rs,ps,m)
!  ===================================================================
!
!  --- this routine interpolates functn f on grid r using a least square fit
!      returns value ps and its derivative dps at rs 
!
!  modified by yw, 11/20/03
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nrs
   integer (kind=IntKind), optional :: m
   integer (kind=IntKind) :: m_poly, i, irs
!
   real (kind=RealKind), intent(in) :: r(nr), f(nr), rs(nrs)
   real (kind=RealKind), intent(out) :: ps(nrs)
   real (kind=RealKind), allocatable :: c(:)
!
   if ( present(m) ) then
      m_poly = m
   else
      m_poly = 3
   endif
   allocate( c(m_poly+1) )
!  -------------------------------------------------------------------
   call leastSqFit_r(nr,r,f,m_poly,c)
!  -------------------------------------------------------------------
   do irs = 1,nrs
      ps(irs) = c(1)
      do i = 1,m_poly
         ps(irs) = ps(irs) + (rs(irs)**i)*c(i+1)
      enddo
   enddo
   deallocate(c)
!
   end subroutine leastSqFit_interp_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine leastSqFit_cinterp_v(nr,r,f,nrs,rs,ps,m)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: nr, nrs
   integer (kind=IntKind), optional :: m
   integer (kind=IntKind) :: m_poly, i, irs
!
   real (kind=RealKind), intent(in) :: r(nr), rs(nrs)
   complex (kind=CmplxKind), intent(in)  :: f(nr)
   complex (kind=CmplxKind), intent(out) :: ps(nrs)
   complex (kind=CmplxKind), allocatable :: c(:)
!
   if ( present(m) ) then
      m_poly = m
   else
      m_poly = 3
   endif
   allocate( c(m_poly+1) )
!  -------------------------------------------------------------------
   call leastSqFit_c(nr,r,f,m_poly,c)
!  -------------------------------------------------------------------
   do irs = 1,nrs
      ps(irs) = c(1)
      do i = 1,m_poly
         ps(irs) = ps(irs) + (rs(irs)**i)*c(i+1)
      enddo
   enddo
   deallocate(c)
!
   end subroutine leastSqFit_cinterp_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine leastSqFit_r(nr,r,f,m,c)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_GE
!
   implicit none
!  *******************************************************************
!  This routine does a linear square fit given the f(r) and the order
!  of the order m of the polinomial fit returns the coeficients (c) of 
!  the fitting polynom.
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: nr,m
!
   real (kind=RealKind), intent(in) :: r(1:nr)
   real (kind=RealKind), intent(in) :: f(1:nr)
!
   real (kind=RealKind), intent(inout) :: c(1:m+1)
!
   integer (kind=IntKind) :: i, j, ir
   real (kind=RealKind) :: Sum_ri0, Sum_rim, detS
   real (kind=RealKind) :: S(1:m+1,1:m+1),Sxy(1:m+1)
!
   if ( m<0 ) then
      call ErrorHandler("leastSqFit","Wrong polynomial order",m)
   endif
!
   Sxy  = ZERO
   do  i = 0,m
       Sum_ri0 = ZERO
       Sum_rim = ZERO
       do ir = 1,nr
          Sum_ri0  = Sum_ri0 + r(ir)**i
          Sum_rim  = Sum_rim + r(ir)**(2*m-i)
          Sxy(i+1) = Sxy(i+1) + f(ir)*(r(ir)**i)
       enddo
       do j = 0,i/2
          S(j+1,i-j+1) = Sum_ri0
          S(i-j+1,j+1) = Sum_ri0
          S(m-i+j+1,m+1-j) = Sum_rim
          S(m+1-j,m-i+j+1) = Sum_rim
       enddo
   enddo
!
   if ( m==0 ) then
      c(1) = Sxy(1)
   else if ( m==1 ) then
      c(2) = ( S(1,1)*Sxy(2)-Sxy(2)*S(1,2) )/( S(2,2)*S(1,1) - S(1,2)*S(2,1) )
      c(1) = ( Sxy(1) - c(2)*S(1,2) )/S(1,1)
   else
      call MtxInv_GE( m+1, S, detS )
      do i = 1,m+1
         c(i) = ZERO
         do j = 1,m+1
            c(i) = c(i) + Sxy(j)*S(j,i)
         enddo
      enddo
   endif
!
   end subroutine leastSqFit_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine leastSqFit_c(nr,r,f,m,c)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_GE
!
   implicit none
!  *******************************************************************
!  This routine does a linear square fit given the f(r) and the order
!  of the order m of the polinomial fit returns the coeficients (c) of 
!  the fitting polynom.
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: nr,m
!
   real (kind=RealKind), intent(in) :: r(1:nr)
   complex (kind=CmplxKind), intent(in) :: f(1:nr)
!
   complex (kind=CmplxKind), intent(inout) :: c(1:m+1)
!
   integer (kind=IntKind) :: i, j, ir
   complex (kind=CmplxKind) :: Sum_ri0, Sum_rim, detS
   complex (kind=CmplxKind) :: S(1:m+1,1:m+1),Sxy(1:m+1)
!
   if ( m<0 ) then
      call ErrorHandler("leastSqFit","Wrong polynomial order",m)
   endif
!
   Sxy  = CZERO
   do  i = 0,m
       Sum_ri0 = CZERO
       Sum_rim = CZERO
       do ir = 1,nr
          Sum_ri0  = Sum_ri0 + cmplx(r(ir)**i,ZERO,kind=CmplxKind)
          Sum_rim  = Sum_rim + cmplx(r(ir)**(2*m-i),ZERO,kind=CmplxKind)
          Sxy(i+1) = Sxy(i+1) + f(ir)*(r(ir)**i)
       enddo
       do j = 0,i/2
          S(j+1,i-j+1) = Sum_ri0
          S(i-j+1,j+1) = Sum_ri0
          S(m-i+j+1,m+1-j) = Sum_rim
          S(m+1-j,m-i+j+1) = Sum_rim
       enddo
   enddo
!
   if ( m==0 ) then
      c(1) = Sxy(1)
   else if ( m==1 ) then
      c(2) = ( S(1,1)*Sxy(2)-Sxy(2)*S(1,2) )/( S(2,2)*S(1,1) - S(1,2)*S(2,1) )
      c(1) = ( Sxy(1) - c(2)*S(1,2) )/S(1,1)
   else
      call MtxInv_GE( m+1, S, detS )
      do i = 1,m+1
         c(i) = CZERO
         do j = 1,m+1
            c(i) = c(i) + Sxy(j)*S(j,i)
         enddo
      enddo
   endif
!
   end subroutine leastSqFit_c
!  ===================================================================
end module InterpolationModule
