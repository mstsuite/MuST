!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupMixingScheme(LocalNumAtoms, n_spin_pola)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use ScfDataModule, only : isSimpleMixing, isDGAMixing, isBroydenMixing
   use ScfDataModule, only : isPotentialMixing, isChargeMixing
!
   use AtomModule, only : getMixingParam4Rho
   use AtomModule, only : getMixingParam4Pot
   use AtomModule, only : getMixingParam4Mom
   use AtomModule, only : getLocalNumSpecies
!
   use MixingModule, only : setSimpleMixing, setDGAMixing, setBroydenMixing
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms, n_spin_pola
!
   integer (kind=IntKind) :: id, is, ia, num_mix
!
   real (kind=RealKind) :: alpha_mix
!
   is = 1
!  do is = 1,n_spin_pola
      num_mix = 0
      do id=1,LocalNumAtoms
         if (isPotentialMixing()) then
!           ----------------------------------------------------------
            alpha_mix = getMixingParam4Pot(id)
!           ----------------------------------------------------------
         else if (isChargeMixing()) then
            if (is == 1) then
!              -------------------------------------------------------
               alpha_mix = getMixingParam4Rho(id)
!              -------------------------------------------------------
            else
!              -------------------------------------------------------
               alpha_mix = getMixingParam4Mom(id)
!              -------------------------------------------------------
            endif
         else
!           ----------------------------------------------------------
            call ErrorHandler('setupMixingScheme','Unkown mixing quantity')
!           ----------------------------------------------------------
         endif
!
         do ia = 1, getLocalNumSpecies(id)
            num_mix = num_mix + 1
            if (isSimpleMixing()) then
!              -------------------------------------------------------
               call setSimpleMixing(is,num_mix,alpha_mix)
!              -------------------------------------------------------
            else if (isDGAMixing()) then
!              -------------------------------------------------------
               call setDGAMixing(is,num_mix,alpha_mix)
!              -------------------------------------------------------
            else if (isBroydenMixing()) then
!              -------------------------------------------------------
               call setBroydenMixing(is,num_mix,alpha_mix)
!              -------------------------------------------------------
            else
!              -------------------------------------------------------
               call ErrorHandler('main','Unkown mixing scheme')
!              -------------------------------------------------------
            endif
         enddo
      enddo
!  enddo
!
!  if (n_spin_pola == 2 .and. isPotentialMixing()) then
!     if (isSimpleMixing()) then
!        -------------------------------------------------------------
!        call setSimpleMixing(n_spin_pola+1,1,alpha_mix)
!        -------------------------------------------------------------
!     else if (isDGAMixing()) then
!        -------------------------------------------------------------
!        call setDGAMixing(n_spin_pola+1,1,alpha_mix)
!        -------------------------------------------------------------
!     else if (isBroydenMixing()) then
!        -------------------------------------------------------------
!        call setBroydenMixing(n_spin_pola+1,1,alpha_mix)
!        -------------------------------------------------------------
!     endif
!  endif
!
   end subroutine setupMixingScheme
