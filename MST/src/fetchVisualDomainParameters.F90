!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine fetchVisualDomainParameters(na,nb,nc,vcell,v0)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO
!
   use ScfDataModule, only : TableID
!
   use InputModule, only : getKeyValue, getNumKeyValues, isKeyExisting
!
   use ErrorHandlerModule, only : ErrorHandler
!
   implicit none
!
   character (len=80) :: svalue
!
   integer (kind=IntKind), intent(out) :: na, nb, nc
   integer (kind=IntKind) :: GridDim, nabc(3)
!
   real (kind=RealKind), intent(out) :: vcell(3,3), v0(3)
   real (kind=RealKind) :: a0_vgrid
!
   if ( getKeyValue(TableID,'Visual Grid Type (0<D<4)', svalue) == 0 ) then
      read(svalue,*)GridDim
   else
      call ErrorHandler('fetchVisualDomainParameters','Visual Grid Type (0<D<4)','Not exist')
   endif
   if (GridDim < 1 .or. GridDim > 3) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid visual grid type input',GridDim)
   endif
!
   if ( getKeyValue(TableID,'Origin Grid Vector', svalue) == 0 ) then
      read(svalue,*)v0
   else
      call ErrorHandler('fetchVisualDomainParameters','Origin Grid Vector','Not exist')
   endif
!
   vcell = ZERO
   if (getKeyValue(TableID,'Grid Vector 1',3,vcell(:,1)) /= 0) then
      call ErrorHandler('fetchVisualDomainParameters','Grid Vector 1','Not exist')
   endif
   if (GridDim > 1) then
      if (getKeyValue(TableID,'Grid Vector 2',3,vcell(:,2)) /= 0) then
         call ErrorHandler('fetchVisualDomainParameters','Grid Vector 2','Not exist')
      endif
   endif
   if (GridDim > 2) then
      if (getKeyValue(TableID,'Grid Vector 3',3,vcell(:,3)) /= 0) then
         call ErrorHandler('fetchVisualDomainParameters','Grid Vector 3','Not exist')
      endif
   endif
!
   if ( getKeyValue(TableID,'Grid Points', svalue) == 0 ) then
      read(svalue,*)na, nb, nc
   else
      call ErrorHandler('fetchVisualDomainParameters','Grid Points','Not exist')
   endif
!
   if ( getKeyValue(TableID,'Grid Scale', a0_vgrid) /= 0 ) then
      call ErrorHandler('fetchVisualDomainParameters','Grid Scale','Not exist')
   endif
   if ( a0_vgrid > ZERO ) then
      v0  = a0_vgrid*v0
      vcell = a0_vgrid*vcell
   endif
!
   end subroutine fetchVisualDomainParameters
!  ===================================================================
