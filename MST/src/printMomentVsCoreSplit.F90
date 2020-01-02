subroutine printMomentVsCoreSplit(iter,funit)
   use KindParamModule, only : IntKind, RealKind
!
   use SystemModule, only : getAtomName, getNumAtoms
!
   use CoreStatesModule, only : getCoreSplitTable, &
                                getCoreNumStatesTable, getCoreDescriptionTable
!
   use ChargeDistributionModule, only : getGlobalMTSphereMomentTable, &
                                        getGlobalVPCellMomentTable
!
   implicit none
!
   integer (kind=IntKind), intent(in), optional :: iter, funit
   integer (kind=IntKind) :: ig, ic, nc, fu, GlobalNumAtoms
   integer (kind=IntKind), pointer :: CoreNumStatesTable(:)
!
   logical :: FileExist
!
   character (len=2), pointer :: CoreDescriptionTable(:,:)
!
   real (kind=RealKind), pointer :: Mvp_Table(:), Mmt_Table(:), CoreSplitTable(:,:)
!
   if (present(funit)) then
      fu = funit
   else
      fu = 6
   endif
!
   GlobalNumAtoms = getNumAtoms()
!
   if (fu /= 6) then
      inquire(file='MomentVsCoreSplit',exist=FileExist)
      if (FileExist .and. present(iter)) then
         open(unit=fu,file='MomentVsCoreSplit',form='formatted',     &
              status='old',position='append')
      else
         open(unit=fu,file='MomentVsCoreSplit',form='formatted',     &
              status='unknown')
         write(fu,'(/,80(''-''))')
         write(fu,'(/,20x,a)')'**************************************'
         write(fu,'( 20x,a )')'* Output from printMomentVsCoreSplit *'
         write(fu,'(20x,a,/)')'**************************************'
      endif
   endif
!
   if (present(iter)) then
      write(fu,'(/,a,i5)')'# ITERATION :',iter
   endif
!
   CoreNumStatesTable => getCoreNumStatesTable()
   CoreSplitTable => getCoreSplitTable()
   CoreDescriptionTable => getCoreDescriptionTable()
   Mmt_Table => getGlobalMTSphereMomentTable()
   MvP_table => getGlobalVPCellMomentTable()
!
   ig = maxloc(CoreNumStatesTable,1)
   nc = CoreNumStatesTable(ig)
   write(fu,'(35(''=''),$)')
   ic=1
   do while (ic <= nc) 
      write(fu,'(11(''=''),$)')
      if (CoreDescriptionTable(ic,ig)(2:2)=='s') then
         ic = ic + 1
      else
         ic = ic + 2
      endif
   enddo
   write(fu,'(/,a,$)')' Atom   Index      Mvp        Mmt'
   ic = 1
   do while (ic <= nc) 
      write(fu,'(9x,a,$)')CoreDescriptionTable(ic,ig)
      if (CoreDescriptionTable(ic,ig)(2:2)=='s') then
         ic = ic + 1
      else
         ic = ic + 2
      endif
   enddo
   write(fu,'(a)')''
!
   do ig = 1, GlobalNumAtoms
      write(fu,'(2x,a3,2x,i6,2(2x,f9.5),$)')getAtomName(ig), ig, Mvp_Table(ig),Mmt_Table(ig)
      nc = CoreNumStatesTable(ig)
      ic = 1
      do while (ic <= nc)
         write(fu,'(2x,f9.5,$)')CoreSplitTable(ic,ig)
         if (CoreDescriptionTable(ic,ig)(2:2)=='s') then
            ic = ic + 1
         else
            ic = ic + 2
         endif
      enddo
      write(fu,'(a)')''
   enddo
!
   if (fu == 6) then
      write(fu,'(35(''=''),$)')
      ig = maxloc(CoreNumStatesTable,1)
      nc = CoreNumStatesTable(ig)
      ic=1
      do while (ic <= nc) 
         write(fu,'(11(''=''),$)')
         if (CoreDescriptionTable(ic,ig)(2:2)=='s') then
            ic = ic + 1
         else
            ic = ic + 2
         endif
      enddo
      write(fu,'(/)')
   else
      close(unit=fu)
   endif
!
end subroutine printMomentVsCoreSplit
