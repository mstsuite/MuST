module IntegrationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, FOURTH, THIRD, HALF
   use MathParamModule, only : ONE, TWO, THREE, FOUR, CZERO
   use MathParamModule, only : TEN2m10, TEN2m12, TEN2m14, TEN2p14
   use ErrorHandlerModule, only : ErrorHandler
!
public :: calIntegration
   interface calIntegration
      module procedure calIntegration_r, calIntegration_c
      module procedure calIntegration0_r, calIntegration0_c
      module procedure calIntegration1_r, calIntegration1_c
   end interface
!
public:: qexpIntegration
   interface qexpIntegration
      module procedure qexpr, qexpc
   end interface
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegration_r(mode,n,x,f,fint)
!  ===================================================================
   implicit none
!  *******************************************************************
!  1-d integration routine
!
!  mode = 0: integrate forward
!      /= 0: integrate backward 
!
!  n: number of mesh
!  f: integrant function
!  x: integral variable
!  fint: integral result
!     fint(x) = int {from x1 to x} * dt * f(t)
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: mode
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: n2
   integer (kind=IntKind) :: n3
   integer (kind=IntKind) :: nf
   integer (kind=IntKind) :: nstep

   real (kind=RealKind), intent(in) :: x(n)
   real (kind=RealKind), intent(in) :: f(n)
   real (kind=RealKind), intent(out) :: fint(n)

   real (kind=RealKind) :: xm
   real (kind=RealKind) :: xh
   real (kind=RealKind) :: xl
   real (kind=RealKind) :: slm
   real (kind=RealKind) :: smh
   real (kind=RealKind) :: smhsav
 
   if(mode.eq.0) then
      n1=1
      nf=n
      nstep=1
   else 
      n1=n
      nf=1
      nstep=-1
   endif
   n2=n1+nstep
   n3=n2+nstep
   nf=nf-nstep

   fint(n1)=ZERO
   xh=HALF*(x(n1)+x(n2))
!  -------------------------------------------------------------------
   call qexpIntegration(0,x(n2),x(n1),xh,x(n1),x(n2),x(n3),           &
                        f(n1),f(n2),f(n3),slm,smh)
!  -------------------------------------------------------------------
   do i=n2,nf,nstep
      smhsav=smh
      xl=HALF*(x(i-nstep)+x(i))
      xh=HALF*(x(i+nstep)+x(i))
!     ----------------------------------------------------------------
      call qexpIntegration(0,xl,x(i),xh,x(i-nstep),x(i),x(i+nstep),   &
                           f(i-nstep),f(i),f(i+nstep),slm,smh)
!     ----------------------------------------------------------------
      fint(i)=fint(i-nstep)+smhsav+slm
   enddo

   smhsav=smh
   xm=HALF*(x(nf+nstep)+x(nf))
!  -------------------------------------------------------------------
   call qexpIntegration(0,x(nf),xm,x(nf+nstep),x(nf-nstep),x(nf),x(nf+nstep), &
                        f(nf-nstep),f(nf),f(nf+nstep),slm,smh)
!  -------------------------------------------------------------------
   fint(nf+nstep)=fint(nf)+smhsav+smh
!
   end subroutine calIntegration_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegration_c(mode,n,x,f,fint)
!  ===================================================================
   implicit none
!  *******************************************************************
!  1-d integration routine
!
!  mode = 0: integrate forward
!      /= 0: integrate backward 
!
!  n: number of mesh
!  f: integrant function
!  x: integral variable
!  fint: integral result
!     fint(x) = int {from x1 to x} * dt * f(t)
!  *******************************************************************
   integer (kind=IntKind), intent(in) :: mode
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: n1
   integer (kind=IntKind) :: n2
   integer (kind=IntKind) :: n3
   integer (kind=IntKind) :: nf
   integer (kind=IntKind) :: nstep
 
   real (kind=RealKind), intent(in) :: x(n)
   real (kind=RealKind) :: xm
   real (kind=RealKind) :: xh
   real (kind=RealKind) :: xl
 
   complex (kind=CmplxKind), intent(in) :: f(n)
   complex (kind=CmplxKind), intent(out) :: fint(n)
   complex (kind=CmplxKind) :: slm
   complex (kind=CmplxKind) :: smh
   complex (kind=CmplxKind) :: smhsav
 
   if(mode.eq.0) then
      n1=1
      nf=n
      nstep=1
   else 
      n1=n
      nf=1
      nstep=-1
   endif
   n2=n1+nstep
   n3=n2+nstep
   nf=nf-nstep
 
   fint(n1)=CZERO
   xh=HALF*(x(n1)+x(n2))
!  ------------------------------------------------------------------- 
   call qexpIntegration(0,x(n2),x(n1),xh,x(n1),x(n2),x(n3),           &
                        f(n1),f(n2),f(n3),slm,smh)
!  ------------------------------------------------------------------- 
   do i=n2,nf,nstep
      smhsav=smh
      xl=HALF*(x(i-nstep)+x(i))
      xh=HALF*(x(i+nstep)+x(i))
      !---------------------------------------------------------------
      call qexpIntegration(0,xl,x(i),xh,x(i-nstep),x(i),x(i+nstep),   &
                           f(i-nstep),f(i),f(i+nstep),slm,smh)
      !---------------------------------------------------------------
      fint(i)=fint(i-nstep)+smhsav+slm
   enddo

   smhsav=smh
   xm=HALF*(x(nf+nstep)+x(nf))
!  -------------------------------------------------------------------
   call qexpIntegration(0,x(nf),xm,x(nf+nstep),x(nf-nstep),x(nf),x(nf+nstep), &
                        f(nf-nstep),f(nf),f(nf+nstep),slm,smh)
!  -------------------------------------------------------------------
   fint(nf+nstep)=fint(nf)+smhsav+smh
!
   end subroutine calIntegration_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegration0_r(nr,r,f,g0,ip0)
!  ===================================================================
   implicit none
!  *******************************************************************
!  this routine integrates real function f on a 1-d mesh r by interpolation
!  ip0: integration power:
!      g = int dr * r^ip0 * f(r)
!  the partial integrals (up to r(k)) are stored in g(k). The mesh
!  is arbitrary
!  xgz  ornl 1998
!
!  modlfied by yw psc 2003
!  *******************************************************************
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: ip0
!
   real (kind=RealKind), intent(in) :: r(nr)
   real (kind=RealKind), intent(in) :: f(nr)
   real (kind=RealKind), intent(out) :: g0
!
   integer (kind=IntKind) :: ip, ip1, ip2, i, i0, i1, i2, i3, i4, j, k
!
   real (kind=RealKind) :: pip1, pip2
   real (kind=RealKind) :: rs(4), rip1, wt, dr, h1, h2, rj
!
   real (kind=RealKind) :: cr(4), fr(4)
   real (kind=RealKind) :: gj, fj, gr, c1, c2
!
   integer (kind=IntKind), parameter :: nj=6
   integer (kind=IntKind), parameter :: njm1=nj-1
   real (kind=RealKind), parameter :: dnj=ONE/(ONE*nj)
   real (kind=RealKind), parameter :: wt0=TWO/(THREE*nj)
   real (kind=RealKind), parameter :: wt1=THREE*wt0
!
!!!ip=2*ip0+1
   ip=2*(2*ip0+1)+1
!
   pip1=ONE/real(ip+1,RealKind)
   pip2=ONE/real(ip+2,RealKind)
   i1=1
   i2=2
   i3=3
   i4=4
   i0=1
   g0 = ZERO
   if (r(1) < TEN2m10) then
      rs(1)=ZERO
!!!   fr(1)=f(1)
      fr(1)=TWO*f(1)
   else
!!!   rs(1)=sqrt(r(1))
      rs(1)=r(1)**FOURTH
!!!   fr(1)=f(1)
      fr(1)=TWO*f(1)
   endif
   do i=2,4
!!!   rs(i)=sqrt(r(i))
      rs(i)=r(i)**FOURTH
!!!   fr(i)=f(i)
      fr(i)=TWO*f(i)
   enddo
   ip1=1
   ip2=4
   do i=1,nr-1
      ip1=ip1+1
      if(i.gt.2.and.ip2.le.nr) then
!!!      rs(i4)=sqrt(r(ip2))
         rs(i4)=r(ip2)**FOURTH
!!!      fr(i4)=f(ip2)
         fr(i4)=TWO*f(ip2)
      endif
!     -------------------------------------------------------------
      call fit(4,rs,fr,i0,cr)
!     -------------------------------------------------------------
!     First integrate a linear function going through f(i) and f(ip1)
!     =============================================================
      rip1=rs(mod(i0,4)+1)
      gj=ONE
      fj=rip1
      do j=1,ip
         gj=gj*rs(i0)+fj
         fj=fj*rip1
      enddo
      fj=(gj*rs(i0)+fj)*pip2
      gj=gj*pip1
      gr=(fr(i0)-cr(3)*rs(i0))*gj+cr(3)*fj
!     =============================================================
!     Using Simpson's rule. No contribution from j=0 and j=nj.
!     =============================================================
      dr=rip1-rs(i0)
      h1=dnj*dr
      h2=ZERO
      wt=wt0
      do j=1,njm1
         h2=h2+h1
         wt=wt1-wt
         fj=wt
         rj=rs(i0)+h2
         do k=1,ip
            fj=fj*rj
         enddo
         gj=h2*(dr-h2)
         c1=(cr(3)-cr(2))*h2+cr(1)
         c2=cr(4)*gj
         gr=gr-fj*c1*c2/(ONE+c2)
      enddo
      g0 = g0 + TWO*dr*gr
      if (i > 1 .and. ip2 < nr) then
         ip2=ip2+1
         j=i1
         i1=i2
         i2=i3
         i3=i4
         i4=j
      endif
      i0=mod(i0,4)+1
   enddo
!
   end subroutine calIntegration0_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegration1_r(nr,r,f,g,ip0)
!  ===================================================================
   implicit none
!  *******************************************************************
!  this routine integrates real function f on a 1-d mesh r by interpolation
!  ip0: integration power:
!      g = int dr * r^ip0 * f(r)
!  the partial integrals (up to r(k)) are stored in g(k). The mesh
!  is arbitrary
!  xgz  ornl 1998
!
!  modlfied by yw psc 2003
!  *******************************************************************
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: ip0
!
   real (kind=RealKind), intent(in) :: r(nr)
   real (kind=RealKind), intent(in) :: f(nr)
   real (kind=RealKind), intent(out) :: g(nr)
!
   integer (kind=IntKind) :: ip, ip1, ip2, i, i0, i1, i2, i3, i4, j, k
!
   real (kind=RealKind) :: pip1, pip2
   real (kind=RealKind) :: rs(4), rip1, wt, dr, h1, h2, rj
!
   real (kind=RealKind) :: cr(4), fr(4)
   real (kind=RealKind) :: gj, fj, gr, c1, c2
!
   integer (kind=IntKind), parameter :: nj=6
   integer (kind=IntKind), parameter :: njm1=nj-1
   real (kind=RealKind), parameter :: dnj=ONE/(ONE*nj)
   real (kind=RealKind), parameter :: wt0=TWO/(THREE*nj)
   real (kind=RealKind), parameter :: wt1=THREE*wt0
!
!!!ip=2*ip0+1
   ip=2*(2*ip0+1)+1
!
   pip1=ONE/real(ip+1,RealKind)
   pip2=ONE/real(ip+2,RealKind)
   g(1)=ZERO
   i1=1
   i2=2
   i3=3
   i4=4
   i0=1
   if (r(1) < TEN2m10) then
      rs(1)=ZERO
!!!   fr(1)=f(1)
      fr(1)=TWO*f(1)
   else
!!!   rs(1)=sqrt(r(1))
      rs(1)=r(1)**FOURTH
!!!   fr(1)=f(1)
      fr(1)=TWO*f(1)
   endif
   do i=2,4
!!!   rs(i)=sqrt(r(i))
      rs(i)=r(i)**FOURTH
!!!   fr(i)=f(i)
      fr(i)=TWO*f(i)
   enddo
   ip1=1
   ip2=4
   do i=1,nr-1
      ip1=ip1+1
      if(i.gt.2.and.ip2.le.nr) then
!!!      rs(i4)=sqrt(r(ip2))
         rs(i4)=r(ip2)**FOURTH
!!!      fr(i4)=f(ip2)
         fr(i4)=TWO*f(ip2)
      endif
!     -------------------------------------------------------------
      call fit(4,rs,fr,i0,cr)
!     -------------------------------------------------------------
!     First integrate a linear function going through f(i) and f(ip1)
!     =============================================================
      rip1=rs(mod(i0,4)+1)
      gj=ONE
      fj=rip1
      do j=1,ip
         gj=gj*rs(i0)+fj
         fj=fj*rip1
      enddo
      fj=(gj*rs(i0)+fj)*pip2
      gj=gj*pip1
      gr=(fr(i0)-cr(3)*rs(i0))*gj+cr(3)*fj
!     =============================================================
!     Using Simpson's rule. No contribution from j=0 and j=nj.
!     =============================================================
      dr=rip1-rs(i0)
      h1=dnj*dr
      h2=ZERO
      wt=wt0
      do j=1,njm1
         h2=h2+h1
         wt=wt1-wt
         fj=wt
         rj=rs(i0)+h2
         do k=1,ip
            fj=fj*rj
         enddo
         gj=h2*(dr-h2)
         c1=(cr(3)-cr(2))*h2+cr(1)
         c2=cr(4)*gj
         gr=gr-fj*c1*c2/(ONE+c2)
      enddo
      g(ip1)=g(i)+TWO*dr*gr
      if (i > 1 .and. ip2 < nr) then
         ip2=ip2+1
         j=i1
         i1=i2
         i2=i3
         i3=i4
         i4=j
      endif
      i0=mod(i0,4)+1
   enddo
!
   end subroutine calIntegration1_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegration0_c(nr,r,f,g0,ip0)
!  ===================================================================
   implicit none
!  *******************************************************************
!  this routine integrates complex function f on a 1-d mesh r by
!  interpolation
!  ip0: integration power:
!      g = int dr * r^ip0 * f(r)
!  the partial integrals (up to r(k)) are stored in g(k). The mesh
!  is arbitrary
!  xgz  ornl 1998
!
!  modified by yw psc 2003, now, r is the original r mesh, not sqrt(r)
!  *******************************************************************
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: ip0
!
   real (kind=RealKind), intent(in) :: r(nr)
   complex (kind=CmplxKind), intent(in) :: f(nr)
   complex (kind=CmplxKind), intent(out) :: g0
!
   integer (kind=IntKind) :: ip, ip1, ip2, i, i0, i1, i2, i3, i4, j, k
!
   real (kind=RealKind) :: pip1, pip2
   real (kind=RealKind) :: rs(4), rip1, wt, dr, h1, h2, rj
!
   real (kind=RealKind) :: cr(4), ci(4), fr(4), fi(4)
   real (kind=RealKind) :: gj, fj, gr, gi, c1, c2
!
   integer (kind=IntKind), parameter :: nj=6
   integer (kind=IntKind), parameter :: njm1=nj-1
   real (kind=RealKind), parameter :: dnj=ONE/(ONE*nj)
   real (kind=RealKind), parameter :: wt0=TWO/(THREE*nj)
   real (kind=RealKind), parameter :: wt1=THREE*wt0
!
!!!ip=2*ip0+1
   ip=2*(2*ip0+1)+1
!
   pip1=ONE/real(ip+1,RealKind)
   pip2=ONE/real(ip+2,RealKind)
   i1=1
   i2=2
   i3=3
   i4=4
   i0=1
   if (r(1) < TEN2m10) then
      rs(1)=ZERO
!!!   fr(1)=real(f(1),RealKind)
      fr(1)=TWO*real(f(1),RealKind)
!!!   fi(1)=aimag(f(1))
      fi(1)=TWO*aimag(f(1))
   else
!!!   rs(1)=sqrt(r(1))
      rs(1)=r(1)**FOURTH
!!!   fr(1)=real(f(1),RealKind)
      fr(1)=TWO*real(f(1),RealKind)
!!!   fi(1)=aimag(f(1))
      fi(1)=TWO*aimag(f(1))
   endif
   do i=2,4
!!!   rs(i)=sqrt(r(i))
      rs(i)=r(i)**FOURTH
!!!   fr(i)=real(f(i),RealKind)
      fr(i)=TWO*real(f(i),RealKind)
!!!   fi(i)=aimag(f(i))
      fi(i)=TWO*aimag(f(i))
   enddo
   ip1=1
   ip2=4
   g0 = CZERO
   do i=1,nr-1
      ip1=ip1+1
      if (i.gt.2.and.ip2.le.nr) then
!!!      rs(i4)=sqrt(r(ip2))
         rs(i4)=r(ip2)**FOURTH
!!!      fr(i4)=real(f(ip2),RealKind)
         fr(i4)=TWO*real(f(ip2),RealKind)
!!!      fi(i4)=aimag(f(ip2))
         fi(i4)=TWO*aimag(f(ip2))
      endif
!     ----------------------------------------------------------------
      call fit(4,rs,fr,i0,cr)
      call fit(4,rs,fi,i0,ci)
!     ----------------------------------------------------------------
!     First integrate a linear function going through f(i) and f(ip1)
!     ================================================================
      rip1=rs(mod(i0,4)+1)
      gj=ONE
      fj=rip1
      do j=1,ip
         gj=gj*rs(i0)+fj
         fj=fj*rip1
      enddo
      fj=(gj*rs(i0)+fj)*pip2
      gj=gj*pip1
      gr=(fr(i0)-cr(3)*rs(i0))*gj+cr(3)*fj
      gi=(fi(i0)-ci(3)*rs(i0))*gj+ci(3)*fj
!     ================================================================
!     Using Simpson's rule. No contribution from j=0 and j=nj.
!     ================================================================
      dr=rip1-rs(i0)
      h1=dnj*dr
      h2=ZERO
      wt=wt0
      do j=1,njm1
         h2=h2+h1
         wt=wt1-wt
         fj=wt
         rj=rs(i0)+h2
         do k=1,ip
            fj=fj*rj
         enddo
         gj=h2*(dr-h2)
         c1=(cr(3)-cr(2))*h2+cr(1)
         c2=cr(4)*gj
         gr=gr-fj*c1*c2/(ONE+c2)
         c1=(ci(3)-ci(2))*h2+ci(1)
         c2=ci(4)*gj
         gi=gi-fj*c1*c2/(ONE+c2)
      enddo
      g0 = g0+TWO*dr*cmplx(gr,gi,CmplxKind)
      if (i > 1 .and. ip2 < nr) then
         ip2=ip2+1
         j=i1
         i1=i2
         i2=i3
         i3=i4
         i4=j
      endif
      i0=mod(i0,4)+1
   enddo
!
   end subroutine calIntegration0_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calIntegration1_c(nr,r,f,g,ip0)
!  ===================================================================
   implicit none
!  *******************************************************************
!  this routine integrates complex function f on a 1-d mesh r by
!  interpolation
!  ip0: integration power:
!      g = int dr * r^ip0 * f(r)
!  the partial integrals (up to r(k)) are stored in g(k). The mesh
!  is arbitrary
!  xgz  ornl 1998
!
!  modified by yw psc 2003, now, r is the original r mesh, not sqrt(r)
!  *******************************************************************
!
   integer (kind=IntKind), intent(in) :: nr
   integer (kind=IntKind), intent(in) :: ip0
!
   real (kind=RealKind), intent(in) :: r(nr)
   complex (kind=CmplxKind), intent(in) :: f(nr)
   complex (kind=CmplxKind), intent(out) :: g(nr)
!
   integer (kind=IntKind) :: ip, ip1, ip2, i, i0, i1, i2, i3, i4, j, k
!
   real (kind=RealKind) :: pip1, pip2
   real (kind=RealKind) :: rs(4), rip1, wt, dr, h1, h2, rj
!
   real (kind=RealKind) :: cr(4), ci(4), fr(4), fi(4)
   real (kind=RealKind) :: gj, fj, gr, gi, c1, c2
!
   integer (kind=IntKind), parameter :: nj=6
   integer (kind=IntKind), parameter :: njm1=nj-1
   real (kind=RealKind), parameter :: dnj=ONE/(ONE*nj)
   real (kind=RealKind), parameter :: wt0=TWO/(THREE*nj)
   real (kind=RealKind), parameter :: wt1=THREE*wt0
!
!!!ip=2*ip0+1
   ip=2*(2*ip0+1)+1
!
   pip1=ONE/real(ip+1,RealKind)
   pip2=ONE/real(ip+2,RealKind)
   i1=1
   i2=2
   i3=3
   i4=4
   i0=1
   if (r(1) < TEN2m10) then
      rs(1)=ZERO
!!!   fr(1)=real(f(1),RealKind)
      fr(1)=TWO*real(f(1),RealKind)
!!!   fi(1)=aimag(f(1))
      fi(1)=TWO*aimag(f(1))
   else
!!!   rs(1)=sqrt(r(1))
      rs(1)=r(1)**FOURTH
!!!   fr(1)=real(f(1),RealKind)
      fr(1)=TWO*real(f(1),RealKind)
!!!   fi(1)=aimag(f(1))
      fi(1)=TWO*aimag(f(1))
   endif
   do i=2,4
!!!   rs(i)=sqrt(r(i))
      rs(i)=r(i)**FOURTH
!!!   fr(i)=real(f(i),RealKind)
      fr(i)=TWO*real(f(i),RealKind)
!!!   fi(i)=aimag(f(i))
      fi(i)=TWO*aimag(f(i))
   enddo
   ip1=1
   ip2=4
   g(1)=CZERO
   do i=1,nr-1
      ip1=ip1+1
      if (i.gt.2.and.ip2.le.nr) then
!!!      rs(i4)=sqrt(r(ip2))
         rs(i4)=r(ip2)**FOURTH
!!!      fr(i4)=real(f(ip2),RealKind)
         fr(i4)=TWO*real(f(ip2),RealKind)
!!!      fi(i4)=aimag(f(ip2))
         fi(i4)=TWO*aimag(f(ip2))
      endif
!     ----------------------------------------------------------------
      call fit(4,rs,fr,i0,cr)
      call fit(4,rs,fi,i0,ci)
!     ----------------------------------------------------------------
!     First integrate a linear function going through f(i) and f(ip1)
!     ================================================================
      rip1=rs(mod(i0,4)+1)
      gj=ONE
      fj=rip1
      do j=1,ip
         gj=gj*rs(i0)+fj
         fj=fj*rip1
      enddo
      fj=(gj*rs(i0)+fj)*pip2
      gj=gj*pip1
      gr=(fr(i0)-cr(3)*rs(i0))*gj+cr(3)*fj
      gi=(fi(i0)-ci(3)*rs(i0))*gj+ci(3)*fj
!        =============================================================
!        Using Simpson's rule. No contribution from j=0 and j=nj.
!        =============================================================
      dr=rip1-rs(i0)
      h1=dnj*dr
      h2=ZERO
      wt=wt0
      do j=1,njm1
         h2=h2+h1
         wt=wt1-wt
         fj=wt
         rj=rs(i0)+h2
         do k=1,ip
            fj=fj*rj
         enddo
         gj=h2*(dr-h2)
         c1=(cr(3)-cr(2))*h2+cr(1)
         c2=cr(4)*gj
         gr=gr-fj*c1*c2/(ONE+c2)
         c1=(ci(3)-ci(2))*h2+ci(1)
         c2=ci(4)*gj
         gi=gi-fj*c1*c2/(ONE+c2)
      enddo
      g(ip1) = g(i)+TWO*dr*cmplx(gr,gi,CmplxKind)
      if (i > 1 .and. ip2 < nr) then
         ip2=ip2+1
         j=i1
         i1=i2
         i2=i3
         i3=i4
         i4=j
      endif
      i0=mod(i0,4)+1
   enddo
!
   end subroutine calIntegration1_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qexpr(ip,xl,xm,xh,x1,x2,x3,y1,y2,y3,slm,smh)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ip
!
   real (kind=RealKind), intent(in) :: xl,xm,xh,x1,x2,x3
   real (kind=RealKind), intent(in) :: y1,y2,y3
   real (kind=RealKind), intent(out) :: slm,smh
!
   real (kind=RealKind) :: xl2,xm2,xh2,xl3,xm3,xh3
   real (kind=RealKind) :: pinv,pinv2,enxl,enxm,enxh
   real (kind=RealKind) :: xenxl,xenxm,xenxh,x2enxl,x2enxm,x2enxh
   real (kind=RealKind) :: b,c,ap,bp,app,bpp
!
   xl2=xl*xl
   xm2=xm*xm
   xh2=xh*xh
   xl3=xl2*xl
   xm3=xm2*xm
   xh3=xh2*xh
!
   b=(y2-y1)/(x2-x1)
   c=( (y3-y1)/(x3-x1) -b )/(x3-x2)
   ap=y1-b*x1 + c*x1*x2
   bp=b -c*(x1+x2)
!
   if (ip.ne.0) then
      pinv=one/ip
      pinv2=pinv*pinv
      enxl=exp(ip*xl)
      enxm=exp(ip*xm)
      enxh=exp(ip*xh)
      xenxl=xl*enxl
      xenxm=xm*enxm
      xenxh=xh*enxh
      x2enxl=xl2*enxl
      x2enxm=xm2*enxm
      x2enxh=xh2*enxh
      app=ap-bp*pinv+two*c*pinv2
      bpp=bp-two*c*pinv
      slm=(app*(enxm-enxl)+bpp*(xenxm-xenxl)+c*(x2enxm-x2enxl))*pinv
      smh=(app*(enxh-enxm)+bpp*(xenxh-xenxm)+c*(x2enxh-x2enxm))*pinv
   else
      slm=ap*(xm-xl)+half*bp*(xm2-xl2)+third*c*(xm3-xl3)
      smh=ap*(xh-xm)+half*bp*(xh2-xm2)+third*c*(xh3-xm3)
   endif
!
   end subroutine qexpr
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qexpc(ip,xl,xm,xh,x1,x2,x3,y1,y2,y3,slm,smh)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: ip
!
   real (kind=RealKind), intent(in) :: xl,xm,xh,x1,x2,x3
!
   real (kind=RealKind) :: xl2,xm2,xh2,xl3,xm3,xh3
   real (kind=RealKind) :: pinv,pinv2,enxl,enxm,enxh
   real (kind=RealKind) :: xenxl,xenxm,xenxh,x2enxl,x2enxm,x2enxh
!
   complex (kind=CmplxKind), intent(in) :: y1,y2,y3
   complex (kind=CmplxKind), intent(out) :: slm,smh
!
   complex (kind=CmplxKind) :: b,c,ap,bp,app,bpp
!
   xl2=xl*xl
   xm2=xm*xm
   xh2=xh*xh
   xl3=xl2*xl
   xm3=xm2*xm
   xh3=xh2*xh
!
   b=(y2-y1)/(x2-x1)
   c=( (y3-y1)/(x3-x1) -b )/(x3-x2)
   ap=y1-b*x1 + c*x1*x2
   bp=b -c*(x1+x2)
!
   if (ip.ne.0) then
      pinv=one/ip
      pinv2=pinv*pinv
      enxl=exp(ip*xl)
      enxm=exp(ip*xm)
      enxh=exp(ip*xh)
      xenxl=xl*enxl
      xenxm=xm*enxm
      xenxh=xh*enxh
      x2enxl=xl2*enxl
      x2enxm=xm2*enxm
      x2enxh=xh2*enxh
      app=ap-bp*pinv+two*c*pinv2
      bpp=bp-two*c*pinv
      slm=(app*(enxm-enxl)+bpp*(xenxm-xenxl)+c*(x2enxm-x2enxl))*pinv
      smh=(app*(enxh-enxm)+bpp*(xenxh-xenxm)+c*(x2enxh-x2enxm))*pinv
   else
      slm=ap*(xm-xl)+half*bp*(xm2-xl2)+third*c*(xm3-xl3)
      smh=ap*(xh-xm)+half*bp*(xh2-xm2)+third*c*(xh3-xm3)
   endif
!
   end subroutine qexpc
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
!!    if(abs(gj) <= TEN2m14) then
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
!!    if(gj > -FOUR .and. gj /= ZERO .and. gj < TEN2p14) then
      if(gj > -FOUR .and. abs(gj) > TEN2m14 .and. gj < TEN2p14) then
         c(1)=c(1)/c(4)
      else
         c(1)=f(i)
         c(4)=ZERO
      endif
   else  ! i1=i2
!     ================================================================
!     set A4=HALF*(f(i)+f(i+1)) and do a 3 point fit
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
!!       if(gj <= -FOUR .or. gj == ZERO .or. gj >= TEN2p14) then
         if(gj <= -FOUR .or. abs(gj) <= TEN2m14 .or. gj >= TEN2p14) then
            c(1)=f(i)
            c(4)=ZERO
         endif
      else
         c(1)=f(i)
         c(4)=ZERO
      endif
   endif  ! f(i1)=f(i2)
   end subroutine fit
!  ===================================================================
end module IntegrationModule
