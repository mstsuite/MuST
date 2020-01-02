!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateMixCmplxValues(LocalNumAtoms,n_spin_pola,CmplxArrayList)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
!
   use MathParamModule, only : ZERO, HALF
!
   use PublicTypeDefinitionsModule, only : MixListCmplxStruct
!
   use ScfDataModule, only : isPotentialMixing, isChargeMixing
!
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
!
   use DataServiceCenterModule, only : getDataStorage,                &
                                       ComplexType, ComplexMark,      &
                                       createDataStorage,             &
                                       isDataStorageExisting
!
   use PotentialModule, only : getOldPotential => getPotential
   use PotentialModule, only : setOldVdif => setVdif
!
   use AtomModule, only : getLocalNumSpecies
!
   use ChargeDensityModule, only : isChargeComponentZero, getRhoLmax
   use ChargeDensityModule, only : getChargeComponentFlag, getChargeDensity
   use ChargeDensityModule, only : getMomentDensity
   use ChargeDensityModule, only : setNewVPCharge => setVPCharge
   use ChargeDensityModule, only : setNewVPMom => setVPMomSize
!
   use PotentialGenerationModule, only : getNewPotential => getPotential
   use PotentialGenerationModule, only : setNewVdif => setVdif
   use PotentialGenerationModule, only : getPotLmax
   use PotentialGenerationModule, only : getPotComponentFlag
!
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms,n_spin_pola
!
   type (MixListCmplxStruct), target :: CmplxArrayList
!
   integer (kind=IntKind) :: id, nr, is, data_size, data_size_ns, ia
   integer (kind=IntKind) :: lmax, jl, jl_nonZero, ind_jl, num_items, p_item
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind) :: factor, chg, mom
   complex (kind=CmplxKind), pointer :: pStore_old(:,:), pStore_new(:,:)
   complex (kind=CmplxKind), pointer :: ptmp1(:), ptmp2(:)
!
   type (MixListCmplxStruct), pointer :: p_CAL
!
   data_size = 0
   num_items = 0
   do id = 1,LocalNumAtoms
      nr = getNumRmesh(id)
      if ( isPotentialMixing() ) then
         lmax = getPotLmax(id)
         flag_jl => getPotComponentFlag(id)
         jl_nonZero = 0
         do jl = 1, (lmax+1)*(lmax+2)/2
            if ( flag_jl(jl) /= 0 ) then
               jl_nonZero = jl_nonZero + 1
            endif
         enddo
!         jl_nonZero = (lmax+1)*(lmax+2)/2
         data_size = max(data_size,nr*jl_nonZero)
      else
         lmax = getRhoLmax(id)
         flag_jl => getChargeComponentFlag(id)
         jl_nonZero = 0
         do jl = 1, (lmax+1)*(lmax+2)/2
            if ( flag_jl(jl) /= 0 ) then
               jl_nonZero = jl_nonZero + 1
            endif
         enddo
!         jl_nonZero = (lmax+1)*(lmax+2)/2
         data_size = max(data_size,nr*jl_nonZero)
      endif
      num_items = num_items + getLocalNumSpecies(id)
   enddo
   data_size_ns = data_size*n_spin_pola+n_spin_pola
   if (.not.isDataStorageExisting('MixingVectorOld')) then
      call ErrorHandler('updateMixCmplxValues','MixingVectorOld not defined')
   endif
   if (.not.isDataStorageExisting('MixingVectorNew')) then
      call ErrorHandler('updateMixCmplxValues','MixingVectorNew not defined')
   endif
!
   pStore_old => getDataStorage('MixingVectorOld',data_size_ns,num_items,ComplexMark)
   pStore_new => getDataStorage('MixingVectorNew',data_size_ns,num_items,ComplexMark)
   p_CAL => CmplxArrayList
   p_item = 0
   do id = 1, LocalNumAtoms
      nr = getNumRmesh(id)
      flag_jl => getPotComponentFlag(id)
      do ia = 1, getLocalNumSpecies(id)
         p_item = p_item + 1
         p_CAL%vector_old => pStore_old(1:data_size_ns,p_item)
         p_CAL%vector_new => pStore_new(1:data_size_ns,p_item)
         do is = 1, n_spin_pola
            ind_jl = data_size*(is -1)+is-1
            if (isPotentialMixing()) then
               lmax = getPotLmax(id)
               do jl = 1, (lmax+1)*(lmax+2)/2
                  if ( flag_jl(jl) /= 0 ) then
                     ptmp2 => getNewPotential("Total",id,ia,is,jl)
                     ptmp2(1:nr) = p_CAL%vector_new(ind_jl+1:ind_jl+nr)
                     ind_jl = ind_jl+ nr
                  endif
               enddo
               if (is == 1) then
                  call setOldVdif( real(p_CAL%vector_new(ind_jl+1),kind=RealKind) )
                  call setNewVdif( real(p_CAL%vector_new(ind_jl+1),kind=RealKind) )
               endif
            else
               lmax = getRhoLmax(id)
               factor = real(3-is*2,kind=RealKind)
               flag_jl => getChargeComponentFlag(id)
               do jl = 1, (lmax+1)*(lmax+2)/2
                  if ( flag_jl(jl) /=0  ) then
                     ptmp1 => getChargeDensity('TotalNew', id, ia, jl)
                     if ( n_spin_pola==2 ) then
                        ptmp2 => getMomentDensity('TotalNew', id, ia, jl)
                     endif
                     if ( n_spin_pola==2 ) then
                        if ( is == 1 ) then
                           ptmp1(1:nr) = p_CAL%vector_new(ind_jl+1:ind_jl+nr)
                        else
                           ptmp1(1:nr) = ptmp1+p_CAL%vector_new(nr*(is-1)+1:nr*is)
                        endif
                        if ( is == 1 ) then
                           ptmp2(1:nr) = p_CAL%vector_new(ind_jl+1:ind_jl+nr)
                        else
                          ptmp2(1:nr) = ptmp2-p_CAL%vector_new(ind_jl+1:ind_jl+nr)
                        endif
                     else
                        ptmp1(1:nr) = p_CAL%vector_new(ind_jl+1:ind_jl+nr)
                     endif
                     ind_jl = ind_jl+ nr
                  endif
               enddo
               if (n_spin_pola==2) then
                  if ( is == 1 ) then
                     chg = real(p_CAL%vector_new(ind_jl+1),kind=RealKind)
                     mom = real(p_CAL%vector_new(ind_jl+1),kind=RealKind)
                  else
                     chg = chg + real(p_CAL%vector_new(ind_jl+1),kind=RealKind)
                     mom = mom - real(p_CAL%vector_new(ind_jl+1),kind=RealKind)
                     call setNewVPCharge(chg,id,ia)
                     call setNewVPMom(mom,id,ia)
                  endif
               else
                  chg = real(p_CAL%vector_new(ind_jl+1),kind=RealKind)
                  call setNewVPCharge(chg,id,ia)
               endif
            endif
         enddo
         p_CAL => p_CAL%next
      enddo
   enddo
!
   end subroutine updateMixCmplxValues
