program testExchCorr
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use ScfDataModule, only : n_spin_pola, excorr_name
!
   use MathParamModule, only : ZERO
!
   use ExchCorrFunctionalModule, only : initExchCorrFunctional,   &
                                        endExchCorrFunctional,    &
                                        isLDAFunctional,          &
                                        isGGAFunctional,          &
                                        isMGGAFunctional,         &
                                        isHybridFunctional,       &
                                        calSphExchangeCorrelation,&
                                        calExchangeCorrelation,   &
                                        getExchCorrPot, getExchCorrEnDen
!
   implicit   none
!
   integer (kind=IntKind) :: i, is
   integer (kind=IntKind) :: iprint = 0
!
   real (kind=RealKind) :: rho(5) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
   real (kind=RealKind) :: der_rho(5) = (/0.1, 0.2, 0.3, 0.4, 0.0/)
   real (kind=RealKind) :: mag(5) = (/0.1, 0.2, 0.3, 0.4, 0.0/)
   real (kind=RealKind) :: der_mag(5) = (/0.3, 0.4, 0.5, 0.0, 0.0/)
   real (kind=RealKind) :: vxc(5,2), exc(5)
   real (kind=RealKind), pointer :: p_vxc(:), p_exc(:)
!
!  -------------------------------------------------------------------
   call startProcess()
!  -------------------------------------------------------------------
!
!  ===================================================================
!  initialize exchange-correlation scheme, potential generation, total
!  energy calculation, and charge distribution book-keeping
!  -------------------------------------------------------------------
   call initExchCorrFunctional(n_spin_pola,excorr_name,iprint)
!  -------------------------------------------------------------------
   if (iprint >= 0) then
      if (isLDAFunctional()) then
         write(6,'(/,2x,a)')'========================================='
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: LDA'
         write(6,'(2x,a,/)')'========================================='
      else if (isGGAFunctional()) then
         write(6,'(/,2x,a)')'========================================='
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: GGA'
         write(6,'(2x,a,/)')'========================================='
      else if (isHybridFunctional()) then
         write(6,'(/,2x,a)')'============================================'
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: Hybrid'
         write(6,'(2x,a,/)')'============================================'
      else if (isMGGAFunctional()) then
         write(6,'(/,2x,a)')'=========================================='
         write(6,'(2x,a)')  'Exchange-Correlation Functional Type: MGGA'
         write(6,'(2x,a,/)')'=========================================='
      else
         call ErrorHandler('testExchCorr','Unknown exchange-correlation functional type')
      endif
   endif
!
   write(6,'(/)')
   do i = 1, 5
!     ----------------------------------------------------------------
      call calSphExchangeCorrelation(rho(i),der_rho(i),mag(i),der_mag(i))
!     ----------------------------------------------------------------
      do is = 1, n_spin_pola
!        -------------------------------------------------------------
         vxc(i,is) = getExchCorrPot(is)
         exc(i) = getExchCorrEnDen()
!        -------------------------------------------------------------
      enddo
   enddo
   do is = 1, n_spin_pola
      do i = 1, 5
         write(6,'(a,2i3,2x,2f15.8)')'i, is, vxc,exc = ',i,is,vxc(i,is),exc(i)
      enddo
   enddo
   write(6,'(/)')
!
!  -------------------------------------------------------------------
   call calSphExchangeCorrelation(5,rho,der_rho,mag,der_mag)
!  -------------------------------------------------------------------
   do is = 1, n_spin_pola
!     ----------------------------------------------------------------
      p_vxc => getExchCorrPot(5,is)
      p_exc => getExchCorrEnDen(5)
!     ----------------------------------------------------------------
      do i = 1, 5
         write(6,'(a,2i3,2x,4f15.8)')'i, is, rho, mom, vxc,exc = ',i,is, &
                                      rho(i),mag(i),p_vxc(i),p_exc(i)
      enddo 
   enddo
!
!  -------------------------------------------------------------------
   call endExchCorrFunctional()
   call finishProcess()
!  -------------------------------------------------------------------
!
end program testExchCorr
