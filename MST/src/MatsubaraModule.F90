module MatsubaraModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use PublicParamDefinitionsModule, only : NicholsonPoints, GaussianPoints
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
!
public :: initMatsubara,        &
          endMatsubara,         &
          getNumMatsubaraPoles, &
          calMatsubaraPoles
!
private
!
   integer (kind=IntKind) :: NumPoles
   integer (kind=IntKind) :: PoleType
   integer (kind=IntKind) :: print_level
!
   real (kind=RealKind) :: Chempot
   real (kind=RealKind) :: Temperature
   real (kind=RealKind) :: FD_Width
!
   logical :: Initialized = .false.
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMatsubara(eGrid,Temp,eWidth,NumPs,iprint)
!  ===================================================================
   use MathParamModule, only : ZERO, ONE, TWO, TEN2m6
   use PhysParamModule, only : Boltzmann
   implicit none
!
   integer (kind=IntKind), intent(in) :: eGrid
   integer (kind=IntKind), intent(in), optional :: NumPs, iprint
   integer (kind=IntKind) :: nume
   integer (kind=IntKind), parameter :: nkbts = 5
!
   real (kind=RealKind), intent(in) :: Temp
   real (kind=RealKind), intent(in), optional :: eWidth
   real (kind=RealKind) :: w, sigma, kbt
!
   if (Temp < ZERO) then
!     ----------------------------------------------------------------
      call ErrorHandler('initMatsubara','Temperature < 0',Temp)
!     ----------------------------------------------------------------
   else if (present(NumPs)) then
      if (NumPs < 1) then
!        -------------------------------------------------------------
         call ErrorHandler('initMatsubara','Number energy points < 1',NumPs)
!        -------------------------------------------------------------
      endif
   else if (present(iprint)) then
      print_level = iprint
   else
      print_level = -1
   endif
!
   Temperature = Temp
!
   if (eGrid == NicholsonPoints) then
      if (present(eWidth)) then
         if (eWidth < TEN2m6) then
!           ----------------------------------------------------------
            call ErrorHandler('initMatsubara',                        &
                              'Energy band width is too small',eWidth)
!           ----------------------------------------------------------
         endif
         w = eWidth
      else
         w = ONE
      endif
      if (present(NumPs)) then
         NumPoles = NumPs
         sigma=TWO*NumPoles*kbt
      else
         kbt = Temperature*Boltzmann
         NumPoles = 0
         sigma=TWO*NumPoles*kbt
         do while (TWO*sigma < w+nkbts*kbt)
            NumPoles=NumPoles+1
            sigma=TWO*NumPoles*kbt
         enddo
      endif
      FD_Width = TWO*sigma+nkbts*kbt
   else if (eGrid == GaussianPoints) then
      if (present(NumPs)) then
         NumPoles = NumPs
      else
         NumPoles = 30
      endif
      FD_Width = ONE  ! For now, this number is set arbitrarily in this case
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initMatsubara','Invalid e-grid type',eGrid)
!     ----------------------------------------------------------------
   endif
!
   PoleType = eGrid
   Initialized = .true.
!
   end subroutine initMatsubara
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMatsubara()
!  ===================================================================
   implicit none
!
   NumPoles = 0
   PoleType = -1
   Chempot = 0.0d0
   Initialized = .false.
!
   end subroutine endMatsubara
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumMatsubaraPoles() result(n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
      call ErrorHandler('getNumMatsubaraPoles','MatsubaraModule is not initialized')
   endif
!
   n = NumPoles
!
   end function getNumMatsubaraPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calMatsubaraPoles(ebot,etop,epole,weight,core_top)
!  ===================================================================
   use MathParamModule, only : CZERO
   implicit none
!
   real (kind=RealKind), intent(in) :: ebot, etop
   real (kind=RealKind), intent(in), optional :: core_top
   real (kind=RealKind) :: etopcor
!
   complex (kind=CmplxKind), intent(out) :: epole(:), weight(:)
!
   if (.not.Initialized) then
      call ErrorHandler('calMatsubaraPoles','MatsubaraModule is not initialized')
   else if (size(epole) < NumPoles) then
      call ErrorHandler('calMatsubaraPoles','size of epole array < NumPoles', &
                        size(epole),NumPoles)
   else if (size(weight) < NumPoles) then
      call ErrorHandler('calMatsubaraPoles','size of weight array < NumPoles',&
                        size(weight),NumPoles)
   endif
!
   if (present(core_top)) then
      etopcor = core_top
   else
      etopcor = -2.5d0
   endif
!
   epole = CZERO; weight = CZERO
!
   if (PoleType == NicholsonPoints) then
!     ----------------------------------------------------------------
      call calNicholsonPoles(etopcor,ebot,etop,NumPoles,epole,weight, &
                             print_level,'none')
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call calPolyFermi(NumPoles,epole,weight,mu=etop)
!     ----------------------------------------------------------------
   endif
!
   end subroutine calMatsubaraPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calNicholsonPoles(etopcor,ebot,chempot,nume,epole,weight, &
                                iprint,istop)
!  ===================================================================
   use MathParamModule, only : ZERO, ONE, TWO, PI, SQRTm1, CZERO, TEN2m8
   use PhysParamModule, only : Boltzmann
!
   implicit   none
!
   character (len=*), intent(in) :: istop
   character (len=25), parameter :: sname = 'calNocholsonPoles'
!
   integer (kind=IntKind), intent(in) :: nume
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind) :: ie, je, nmatp
!
   real (kind=RealKind), intent(in) :: etopcor
   real (kind=RealKind), intent(in) :: ebot
   real (kind=RealKind), intent(in) :: chempot
!
   real (kind=RealKind) :: kbt, sigma
!
   complex (kind=CmplxKind), intent(out) :: epole(nume)
   complex (kind=CmplxKind), intent(out) :: weight(nume)
!
   complex (kind=CmplxKind) :: exptht
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
!  ==================================================================
!  constructs semi-circle energy contour and a set of Matsubara
!  frequencies i*PI*kb*t*(2*m+1), where m=1,2,3,.....................
!  this is for Don Nicholson et al's approximate Fermi function......
!  with sigma = 2*N*kB*T (see: D.M.C. Nicholson et al. PRB50, 14686(1994)) 
!  ------------------------------------------------------------------
!
!                                 |
!                               x | x Im(e)
!                           x     |     x       Matsubara poles
!                       x         |         x
!                     x           |           x
!                   x             |             x
!                  x              |              x
!          -+------+--------------|--------------+----------
!        etopcor  ebot            |            chempot   Real(e)
!                  |<----------2*sigma---------->|
!
!  -------------------------------------------------------------------
!  Description of the variables:
!  -------------------------------------------------------------------
!  Boltzmann   : the Boltzmann constant
!  chempot     : Chemical potential
!  etopcor     : upper limit of the core states energy. It is used for 
!                checking if the width of the generated contour.
!  ebot        : bottom of the valence band on the real axis.
!  nume        : total number of generated Matsubara energies.
!  iprint      : print level control (>0 for print of energies).
!
!  -------------------------------------------------------------------
!  returns   : nume, (epole(ie), weight(ie),ie=1,nume)
!  *******************************************************************
!
!  ===================================================================
!  convert temperature from Kelvin to Rydbergs.
!  ===================================================================
   kbt=Temperature*Boltzmann
   sigma=TWO*nume*kbt
   nmatp=2*nume
!
!  ==================================================================
   if (iprint >= 0) then
      write(6,'(/,a)')' ***********************************'
      write(6,'(  a)')' *        calNicholsonPoles        *'
      write(6,'(a,/)')' ***********************************'
      write(6,'('' ebot                :'',d12.5)') ebot
      write(6,'('' chempot             :'',d12.5)') chempot
      write(6,'('' T (Kelvin)          :'',d12.3)') temperature
      write(6,'('' T (Rydberg)         :'',d12.5)') kbt
   endif
!
   epole = CZERO
   weight = CZERO
!
   if (abs(temperature) < TEN2m8) then
      call WarningHandler(sname,'Temperature is too small',temperature)
      return
   endif
!
!  ===================================================================
!  make sure that bottom of contour is well above uppermost core st.
   if (etopcor+0.1 > chempot-FD_Width) then
      call ErrorHandler(sname,'etopcor+0.1 > chempot-FD_Width',       &
                        etopcor,chempot-FD_Width)
   endif
!
!  ===================================================================
   if(iprint >= 0) then
      write(6,'('' No. Nicholson poles :'',i5)')nume
      write(6,'('' Sigma               :'',f13.8)') sigma
      write(6,'('' f[ebot]             :'',d15.8)')                   &
                   ONE/(ONE+((ebot-chempot+sigma)/sigma)**nmatp)
      write(6,'('' Chempot-2*sigma     :'',f13.8)') chempot-2*sigma
      write(6,'('' f[etopcor]          :'',d15.8)')                   &
                   ONE/(ONE+((etopcor-chempot+sigma)/sigma)**nmatp)
   endif
!
!  ===================================================================
!  generate poles and residue weights: epole(ie), weight(ie)..........
!    epole corresponds to the pole of fP(z) = 1/{1+[(z-mu+sigma)/sigma]^{2N}}
!    where sigma = 2*N*kB*T (see: D.M.C. Nicholson et al. PRB50, 14686(1994)) 
!    and weight is the residule of fP at the pole. They are given by
!       epole = -sigma + sigma*exp[i*(2*ie-1)*PI/(2*N)], ie = 1,2,...,N
!      weight = -2*PI*i*sigma*exp[i*(2*ie-1)*PI/(2*N)]/(2*N)
!             = -2*PI*i*kB*T*exp[i*(2*ie-1)*PI/(2*N)]
!  To derive the residual, the following algebra is used:
!         a^n - b^n = (a-b) * \sum_{j=0}^{n-1}( a^j * b^{n-1-j} )
!  ===================================================================
   do ie=nume,1,-1
      je=nume+1-ie
      exptht=exp(SQRTm1*(2*ie-1)*PI/nmatp)
      epole(je)=chempot-sigma*(ONE-exp(SQRTm1*(2*ie-1)*PI/nmatp))
      weight(je)=-TWO*PI*kbt*SQRTm1*exptht
   enddo
!
!  ===================================================================
!  printout if needed..............................................
   if (iprint > 0) then
      write(6,'(/)')
      do ie=1,nume
         write(6,'(12x,''n,e,de:'',i5,4e12.5)')ie,epole(ie),weight(ie)
      enddo
      do ie=1,nume
         write(6,'(12x,''n,e,f[real(e)]:'',i5,3e12.5)') ie,epole(ie), &
               ONE/(ONE+((dble(epole(ie))-chempot+sigma)/sigma)**nmatp)
      enddo
   endif
!
!  ===================================================================
   if (nocaseCompare(istop,sname)) then
!     ----------------------------------------------------------------
      call StopHandler(sname)
!     ----------------------------------------------------------------
   endif
!
   end subroutine calNicholsonPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calPolyFermi(ne,xg,wg,nterms,mu)
!  ===================================================================
!  This subroutine is adapted from code written by X-G Zhang and his student.
!  It is modified and integrated into MST package by Yang Wang
!
!  mu - The chemical potential
!  xg - Gaussian points in Ryd units (output)
!  wg - Gaussian quadrature weights in Ryd units (output)
!  ne - number of Gaussian points (input)
!  ******************************************************************* 
!  subroutine polyfermi(Temp,xg,wg,ne,nterms,mu)   
!
   use MathParamModule, only : PI, SQRTm1, CONE, PI2, ZERO
   use PhysParamModule, only : Boltzmann
!
   implicit none
!
   integer (kind=IntKind), parameter :: ntm=100000
   integer (kind=IntKind), parameter :: nem=100
!
   complex (kind=CmplxKind) :: z,s,ephi,emphi,b,c,w
!
   real (kind=RealKind) :: phi(ntm),x(ntm),wr(ntm),pl(ntm,3),sum(0:nem)
   real (kind=RealKind) :: al(ntm),bl(ntm)
   real (kind=RealKind) :: pe(0:nem,ntm),xe(ntm),we(ntm)
   real (kind=RealKind) :: gamma, r, q, p, a, sum1, sum2, p1, p2, x1, x2
!
   integer (kind=IntKind) :: ipvt(nem)
   integer (kind=IntKind) :: n, i, ilm1, ilm2, il, nz, l, info
!
   integer (kind=IntKind), intent(in) :: ne
   integer (kind=IntKind), intent(in), optional :: nterms
   real (kind=RealKind), intent(in), optional :: mu
   real (kind=RealKind) :: T
   complex (kind=CmplxKind), intent(out) :: xg(ne),wg(ne)
!      
!  Temperature in Kelvin
   T = Temperature*Boltzmann  ! Change to Rydberg units. Add by Yang Wang
   gamma=3.d0-sqrt(8.d0)
!
   if (ne.gt.nem) then
      call ErrorHandler('calPolyFermi','ne > nem',ne,nem)
   endif
!
!!!n=2+int(0.5+0.1d0/T)   ! ywg
   if (present(nterms)) then
      n = nterms
   else 
      n=max(2+int(0.5+0.1d0/T),ne+1)   ! ywg
   endif
!
   if (n.gt.ntm .or. ne.gt.n) then
      call ErrorHandler('calPolyFermi','n > ntm or n <= ne',n,ntm,ne)
   endif
!
   pe = ZERO
   sum(0)=0.d0
   sum1=0.d0
   r=(1.d0+gamma)/(2.d0*n)
   q=(1.d0-gamma)/(2.d0*n)
   p=(1.d0-gamma*gamma)/(8.d0*n)
   a=r*r
   do i=1,n
      phi(i)=(i-0.5d0)*pi/n
      ephi=exp(sqrtm1*phi(i))
      emphi=CONE/ephi
      b=r*emphi+q*ephi
      c=emphi*emphi-1.d0
      s=sqrt(b*b-a*c)
      z=ephi*(-b+s)/(0.5d0*a)
      w=-sqrtm1*(1.d0+0.5d0*z*r)*(1.d0-z*q)/(1.d0-z*p)*T
      if (ne.eq.n) then
         xg(i)=z*T
         wg(i)=w
      else ! ne < n, since code stops if ne > n
         wr(i)=1.d0     ! changed
         x(i)=phi(i)    ! changed
         pl(i,1)=1.d0
         sum(0)=sum(0)+wr(i)
         sum1=sum1+wr(i)*x(i)
      endif
   enddo ! i
!      
   if (ne.lt.n) then
      sum(1)=0.d0
      sum2=0.d0
      do i=1,n
         al(1)=sum1/sum(0)
         bl(1)=0.d0
         pl(i,2)=x(i)-al(1)
         p2=pl(i,2)*pl(i,2)
         sum(1)=sum(1)+wr(i)*p2
         sum2=sum2+wr(i)*x(i)*p2
      enddo
!
      ilm1=1
      il=2
      do l=2,ne
         sum(l)=0.d0
         sum1=sum2
         sum2=0.d0
         ilm2=ilm1
         ilm1=il
         il=mod(il,3)+1
         do i=1,n
            al(l)=sum1/sum(l-1)
            bl(l)=sum(l-1)/sum(l-2)
            pl(i,il)=(x(i)-al(l))*pl(i,ilm1)-bl(l)*pl(i,ilm2)
            p2=pl(i,il)*pl(i,il)
            sum(l)=sum(l)+wr(i)*p2
            sum2=sum2+wr(i)*x(i)*p2
         enddo
      enddo ! l
!
      nz=0
      do i=1,n
         if (pl(i,il).eq.0.d0) then
            nz=nz+1
            pe(0,nz)=1.d0
            xe(nz)=x(i)
            pe(1,nz)=xe(nz)-al(1)
            do l=2,ne
               pe(l,nz)=(xe(nz)-al(l))*pe(l-1,nz)-bl(l)*pe(l-2,nz)
            enddo ! l
         else if (i.gt.1) then
            if (pl(i,il)*pl(i-1,il).lt.0.d0) then
               nz=nz+1
               pe(0,nz)=1.d0
               x1=x(i-1)
               x2=x(i)
               p1=pl(i-1,il)
               p2=pl(i,il)
!
               do while (abs(x1-x2).gt.1.d-14)
                  xe(nz)=x1-p1*(x1-x2)/(p1-p2)
                  pe(1,nz)=xe(nz)-al(1)
                  do l=2,ne
                     pe(l,nz)=(xe(nz)-al(l))*pe(l-1,nz)-bl(l)*pe(l-2,nz)
                  enddo ! l
                  if (abs(p1).lt.abs(p2)) then
                     x2=xe(nz)
                     p2=pe(ne,nz)
                  else
                     x1=xe(nz)
                     p1=pe(ne,nz)
                  endif
               enddo
            endif
         endif
      enddo  ! i
!     we(nz)=0.d0
      we = ZERO
      we(1)=sum(0)
!
      call dgesv(ne,1,pe(0,1),nem+1,ipvt,we,nem,info)
!
      do nz=1,ne
         ephi=exp(dcmplx(0.d0,xe(nz)))   !changed
         emphi=conjg(ephi)
         b=r*emphi+q*ephi
         c=emphi*emphi-1.d0
         s=sqrt(b*b-a*c)
         z=ephi*(-b+s)/(0.5d0*a)
         w=-SQRTm1*(1.d0+0.5d0*z*r)*(1.d0-z*q)/(1.d0-z*p)
         xg(nz)=z*T
         wg(nz)=w*we(nz)*T    ! changed
      enddo ! nz
   endif ! ne
!
!  ===================================================================
!  The following lines are added by Yang Wang
!  ===================================================================
   wg = wg*PI2  ! added by Yang Wang
   if (present(mu)) then
      xg = xg + mu
   endif
!
   end subroutine calPolyFermi
!  ==================================================================
end module MatsubaraModule
