!  ==================================================================
   subroutine printSystemMovie(iscf,isd,fu,EnergyOffset)
!  ==================================================================
   use KindParamModule, only: IntKind, RealKind
   use MathParamModule, only : ZERO
   use ChemElementModule, only : MaxLenOfAtomName
   use MPPModule, only : MyPE
!
   use ScfDataModule, only : printScfData, ntstep, &
                             n_spin_pola, n_spin_cant, NumSS_IntEs
   use SystemModule, only : printSystem, getAtomPosition, getAtomName, &
                     getNumAtoms, getForce, getRMSInfo, getAtomEnergy, &
                     getMomentDirection, getMomentDirectionOld, &
                     getConstrainField, getSiteLIZ
   use PotentialTypeModule, only : printPotentialType
   use SystemVolumeModule, only : getAtomicVPVolume, &
                                  getAtomicInscribedVolume
   use ChargeDistributionModule, only : getGlobalMTSphereElectronTable, &
                                        getGlobalVPCellElectronTable,   &
                                        getGlobalMTSphereMomentTable,   &
                                        getGlobalVPCellMomentTable,     &
                                        getGlobalExchangeEnergyTable
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: iscf, isd, fu, EnergyOffset
!
   logical :: file_exist
!
   integer (kind=IntKind) :: i, is, nspin, GlobalNumAtoms
   integer (kind=IntKind), pointer :: siteLIZ(:)
!
   character (len=MaxLenOfAtomName), pointer :: atname(:)
   character (len=160) :: fname
!
   real (kind=RealKind), pointer :: atpos(:,:), enpres(:,:), force(:,:), rms(:,:)
   real (kind=RealKind), pointer :: exc_en(:)
   real (kind=RealKind), pointer :: vol_VP(:), vol_MT(:)
   real (kind=RealKind), pointer :: evec(:,:), evec_old(:,:), bcon(:,:)
   real (kind=RealKind), pointer :: q_VP(:), q_MT(:), mom_MT(:), mom_VP(:)
!
   GlobalNumAtoms = getNumAtoms()

   vol_VP => getAtomicVPVolume()
   vol_MT => getAtomicInscribedVolume()
   siteLIZ=> getSiteLIZ()
   atpos  => getAtomPosition()
   atname => getAtomName()
   enpres => getAtomEnergy()
   force  => getForce()
   rms    => getRMSInfo()

   if (n_spin_cant==2) then
      evec_old   => getMomentDirectionOld()
      evec   => getMomentDirection()
   endif
   bcon   => getConstrainField()
   q_VP   => getGlobalVPCellElectronTable()
   q_MT   => getGlobalMTSphereElectronTable()

   if ( n_spin_pola==2 ) then
      exc_en => getGlobalExchangeEnergyTable()
      mom_MT => getGlobalMTSphereMomentTable()
      mom_VP => getGlobalVPCellMomentTable()
   endif
!

   write(fu,'(a)')   '# ****************************************************'
   write(fu,'(a,i4)')'# SCF Iteration: ',iscf
   if (ntstep>1) then
      write(fu,'(a,i4)')'# S-D Time Step: ',isd
   endif
   write(fu,'(a)')   '# ****************************************************'
   if (n_spin_cant==2) then
      do i = 1, GlobalNumAtoms
         write(fu,'(1x,a3,1x,i7,17(1x,f10.6),1x,i4,1x,e16.9,1x,f10.6,1x,e16.9,4(1x,f10.6),4(1x,e9.2))')  &
                  atname(i), i, atpos(1:3,i), force(1:3,i), evec_old(1:3,i), evec(1:3,i), &
                  bcon(1:3,i), vol_MT(i), vol_VP(i), siteLIZ(i), enpres(1,i)-EnergyOffset,&
                  exc_en(i),&
                  enpres(2,i), q_MT(i), q_VP(i), mom_MT(i), mom_VP(i), rms(1:4,i)
      enddo
   else if (n_spin_pola==2) then
      do i = 1, GlobalNumAtoms
         write(fu,'(1x,a3,1x,i7,8(1x,f10.6),1x,i4,1x,e16.9,1x,f10.6,1x,e13.6,4(1x,f10.6),2(1x,e9.2))') &
                  atname(i), i, atpos(1:3,i), force(1:3,i), vol_MT(i), vol_VP(i),       &
                  siteLIZ(i), enpres(1,i)-EnergyOffset, &
                  exc_en(i),&
                  enpres(2,i), q_MT(i), q_VP(i),  &
                  mom_MT(i), mom_VP(i), rms(1:2,i)
      enddo
   else
      do i = 1, GlobalNumAtoms
         write(fu,'(1x,a3,i7,8(1x,f10.6),1x,i4,1x,e16.9,1x,e13.6,2(1x,f10.6),2(1x,e9.2))') &
                  atname(i), i, atpos(1:3,i), force(1:3,i), vol_MT(i), vol_VP(i),    &
                  siteLIZ(i), enpres(1,i)-EnergyOffset, enpres(2,i),                 &
                  q_MT(i), q_VP(i), rms(1:2,i)
      enddo
   endif

!  ==================================================================
   end subroutine printSystemMovie
