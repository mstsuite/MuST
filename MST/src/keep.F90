!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine keep(sdstep,kscf,max_rms, Efermi, TotalEnergy, PV3)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use SystemModule, only : getNumAtoms
   use SystemModule, only : getSystemID, getSystemTitle, getNumAtomTypes
   use SystemModule, only : getNumAtomsOfType, getAtomTypeName
!
   use SystemVolumeModule, only : getSystemVolume
!
   use ScfDataModule, only : isSimpleMixing, isDGAMixing, isBroydenMixing
   use ScfDataModule, only : isPotentialMixing, isChargeMixing
   use ScfDataModule, only : n_spin_pola
!
   use AtomModule, only : getMixingParam4Rho, getMixingParam4Pot
!
   use ChargeDistributionModule, only : getAverageMoment
!
   use BookKeepingModule, only : isBookKeepingNew, &
                                 writeHeadLine,    &
                                 insertColumn,     &
                                 writeRow,         &
                                 getHeadLineValue, &
                                 getLastColumnValue
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: kscf, sdstep
!
   integer (kind=IntKind), save :: E_offset
   integer (kind=IntKind), save :: LastIter = -1
   integer (kind=intKind) :: i, n, m
!
   real (kind=RealKind), intent(in) :: max_rms(4)
   real (kind=RealKind), intent(in) :: Efermi
   real (kind=RealKind), intent(in) :: TotalEnergy, PV3
!
   real (kind=RealKind) :: mom
!
   character (len=80) :: string_tmp, text, myunit
   character (len=8) :: snum
   character (len=30) :: ver
!
   interface
      function getTokenPosition(k,s,n) result(p)
         character (len=*), intent(in) :: s
         integer, intent(in) :: k
         integer, intent(out), optional :: n
         integer :: p
      end function getTokenPosition
   end interface
!
!  TotalEnergy = getEnergyPerAtom()
!  PV3 = getPressurePerAtom()
!
   if (isBookKeepingNew()) then
#include "git_version.h"
      i = getTokenPosition(6,myunit)
      ver = trim(myunit(i:))
!     ----------------------------------------------------------------
      call writeHeadLine('Source Code Version',ver)
!     ----------------------------------------------------------------
      call writeHeadLine('System Description',getSystemTitle())
!     ----------------------------------------------------------------
      write(string_tmp,'(i6)')getNumAtoms()
      call writeHeadLine('Number of Atoms in Unit Cell',string_tmp)
!     ----------------------------------------------------------------
      write(string_tmp,'(i6)')getNumAtomTypes()
      call writeHeadLine('Number of Atomic Species',string_tmp)
!     ----------------------------------------------------------------
      write(snum,'(i8)')getNumAtomsOfType(1); snum = adjustl(snum)
      write(string_tmp,'(''  '',a,''('',a,'')'')')trim(snum),         &
                                                  trim(getAtomTypeName(1))
      n = 4+len_trim(snum)+len_trim(getAtomTypeName(1))
      do i = 2, getNumAtomTypes()
         write(snum,'(i8)')getNumAtomsOfType(i); snum = adjustl(snum)
         m = 4+len_trim(snum)+len_trim(getAtomTypeName(i))
         write(string_tmp(n+1:n+m),'('', '',a,''('',a,'')'')')        &
               trim(snum), trim(getAtomTypeName(i))
         n = n + m
      enddo
      call writeHeadLine('Number of Atoms in Each Type',string_tmp)
!     ----------------------------------------------------------------
      write(string_tmp,'(d15.8)')getSystemVolume()
      call writeHeadLine('Unit Cell Volume (au^3)',string_tmp)
!     ----------------------------------------------------------------
      write(string_tmp,'(f12.5)')getSystemVolume()/real(getNumAtoms(),RealKind)
      call writeHeadLine('Average Atomic Volume (au^3)',string_tmp)
!     ----------------------------------------------------------------
      E_offset = int(TotalEnergy,kind=IntKind)
      write(string_tmp,'(i6)')E_offset
      call writeHeadLine('Energy Offset',string_tmp)
!     ----------------------------------------------------------------
      LastIter = 0
   else if (LastIter < 0) then
      text = getHeadLineValue('Energy Offset')
      read(text,'(i6)')E_offset
      text = getLastColumnValue('  Iter')
      read(text,'(i7)')LastIter
   endif
!
!  -------------------------------------------------------------------
   write(string_tmp,'(i7)')LastIter+kscf
   call insertColumn('Iter',string_tmp)
!  -------------------------------------------------------------------
   if (abs(TotalEnergy-E_offset) >= 1000.0) then
      write(string_tmp,'(f10.4,''    '')')TotalEnergy-E_offset
   else if (abs(TotalEnergy-E_offset) >= 100.0) then
      write(string_tmp,'(f10.5,''    '')')TotalEnergy-E_offset
   else
      write(string_tmp,'(f10.6,''    '')')TotalEnergy-E_offset
   endif
   call insertColumn('Energy',string_tmp)
!  -------------------------------------------------------------------
   if (abs(PV3) >= 1000.0) then
      write(string_tmp,'(f9.3,''    '')')PV3
   else if (abs(PV3) >= 100.0) then
      write(string_tmp,'(f9.4,''    '')')PV3
   else
      write(string_tmp,'(f9.5,''    '')')PV3
   endif
   call insertColumn('3PV',string_tmp)
!  -------------------------------------------------------------------
   write(string_tmp,'(f8.5,''    '')')Efermi
   call insertColumn('Efermi',string_tmp)
!  -------------------------------------------------------------------
   write(string_tmp,'(e10.3,''    '')')max_rms(1)
   call insertColumn('Rms_rho',string_tmp)
!  -------------------------------------------------------------------
   write(string_tmp,'(e10.3,''    '')')max_rms(2)
   call insertColumn('Rms_pot',string_tmp)
!  -------------------------------------------------------------------
   if ( isChargeMixing() ) then
      string_tmp = ' rho'
   else
      string_tmp = ' pot'
   endif
   call insertColumn('Mix',string_tmp)
!  -------------------------------------------------------------------
   if (isSimpleMixing()) then
      string_tmp = 'S'
   else if (isDGAMixing()) then
      string_tmp = 'A'
   else if (isBroydenMixing()) then
      string_tmp = 'B'
   endif
   call insertColumn('Alg',string_tmp)
!  -------------------------------------------------------------------
   if ( isChargeMixing() ) then
      if (getMixingParam4Rho(1) >= 0.99999999) then
         write(string_tmp,'(f7.3,''    '')')getMixingParam4Rho(1)
      else
         write(string_tmp,'(f7.4,''    '')')getMixingParam4Rho(1)
      endif
   else
      if (getMixingParam4Pot(1) >= 0.99999999) then
         write(string_tmp,'(f7.3,''    '')')getMixingParam4Pot(1)
      else
         write(string_tmp,'(f7.4,''    '')')getMixingParam4Pot(1)
      endif
   endif
   call insertColumn('Alpha',string_tmp)
   if (n_spin_pola==2) then
      mom = getAverageMoment()
      write(string_tmp,'(f7.4,''    '')') mom
      call insertColumn('Mom',string_tmp)
   endif
!  -------------------------------------------------------------------
   if (kscf == 1) then
!     ----------------------------------------------------------------
      call writeRow('**',2)
!     ----------------------------------------------------------------
   else if (sdstep==1) then
!     ----------------------------------------------------------------
      call writeRow('*T',2)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call writeRow('  ',2)
!     ----------------------------------------------------------------
   endif
!
   end subroutine keep
