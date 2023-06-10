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
   integer (kind=IntKind) :: GridDim, nabc(3), rstatus
!
   real (kind=RealKind), intent(out) :: vcell(3,3), v0(3)
   real (kind=RealKind) :: a0_vgrid
!
   rstatus = getKeyValue(TableID,'Visual Grid Type (0<D<4)',GridDim)
   if (rstatus /= 0) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid input for Visual Grid Type (0<D<4)')
   else if (GridDim < 1 .or. GridDim > 3) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid visual grid type input',GridDim)
   endif
!
   rstatus = getKeyValue(TableID,'Visual Grid Origin Vector',3,v0)
   if (rstatus /= 0) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid input for Visual Origin Vector')
   endif
!
   vcell = ZERO
   rstatus = getKeyValue(TableID,'Visual Grid Vector 1',3,vcell(:,1))
   if (rstatus /= 0) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid input for Visual Grid Vector 1')
   endif
   if (GridDim > 1) then
      rstatus = getKeyValue(TableID,'Visual Grid Vector 2',3,vcell(:,2))
      if (rstatus /= 0) then
         call ErrorHandler('fetchVisualDomainParameters','Invalid input for Visual Grid Vector 2')
      endif
   endif
   if (GridDim > 2) then
      rstatus = getKeyValue(TableID,'Visual Grid Vector 3',3,vcell(:,3))
      if (rstatus /= 0) then
         call ErrorHandler('fetchVisualDomainParameters','Invalid input for Visual Grid Vector 3')
      endif
   endif
!
   rstatus = getKeyValue(TableID,'Visual Grid Points',GridDim,nabc)
   if (rstatus /= 0) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid input for Visual Grid Points')
   endif
   if (GridDim == 1) then
      na = nabc(1)
      nb = 1; nc = 1
   else if (GridDim == 2) then
      na = nabc(1)
      nb = nabc(2)
      nc = 1
   else
      na = nabc(1)
      nb = nabc(2)
      nc = nabc(3)
   endif
   if (na < 1 .or. nb < 1 .or. nc < 1) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid Visual Grid Points',na,nb,nc)
   endif
!
   rstatus = getKeyValue(TableID,'Visual Grid Scale',a0_vgrid)
   if (rstatus /= 0) then
      call ErrorHandler('fetchVisualDomainParameters','Invalid input for Visual Grid Scale')
   endif
   if ( a0_vgrid > ZERO ) then
      v0  = a0_vgrid*v0
      vcell = a0_vgrid*vcell
   else
      call ErrorHandler('fetchVisualDomainParameters','Invalid Visual Grid Scale',a0_vgrid)
   endif
!
   end subroutine fetchVisualDomainParameters
!  ===================================================================
