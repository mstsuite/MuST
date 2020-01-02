!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calNicholsonPoles(etopcor,ebot,chempot,temperature,     &
                                nume,epole,weight,iprint,istop)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, TWO, PI, SQRTm1, CZERO, TEN2m8
   use PhysParamModule, only : Boltzmann
   use ErrorHandlerModule, only : StopHandler, ErrorHandler, WarningHandler
!
   implicit   none
!
   character (len=*), intent(in) :: istop
   character (len=25), parameter :: sname = 'calNocholsonPoles'
!
   integer (kind=IntKind), intent(out) :: nume
   integer (kind=IntKind), intent(in) :: iprint
!
   integer (kind=IntKind) :: ie, je, nmatp
   integer (kind=IntKind), parameter :: nkbts = 5
!
   real (kind=RealKind), intent(in) :: etopcor
   real (kind=RealKind), intent(in) :: ebot
   real (kind=RealKind), intent(in) :: temperature
   real (kind=RealKind), intent(in) :: chempot
!
   real (kind=RealKind) :: kbt, sigma
!
   complex (kind=CmplxKind), intent(out) :: epole(:)
   complex (kind=CmplxKind), intent(out) :: weight(:)
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
!  temperature : temperature in K
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
      nume = 0
      return
   endif
!
   sigma=ZERO
   nume=0
   sigma=TWO*nume*kbt
   do while(TWO*sigma < chempot-ebot+nkbts*kbt)
      nume=nume+1
      sigma=TWO*nume*kbt
   enddo
   if (nume > size(epole)) then
      call ErrorHandler(sname,'nume > size of epole',nume,size(epole))
   endif
!
   nmatp=2*nume
!
!  ===================================================================
!  make sure that bottom of contour is well above uppermost core st.
   if (etopcor+0.1 > chempot-TWO*sigma-nkbts*kbt) then
      call ErrorHandler(sname,'etopcor+0.1 > ebot-2*nkbts*kbt',       &
                        etopcor,ebot-TWO*nkbts*kbt)
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
