!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine updateMixRealValues(LocalNumAtoms,n_spin_pola,RealArrayList)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
!
   use MathParamModule, only : ZERO, HALF
!
   use ScfDataModule, only : isPotentialMixing, isChargeMixing
!
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
!
   use DataServiceCenterModule, only : getDataStorage, RealType, RealMark, &
                                       createDataStorage, isDataStorageExisting
!
   use PotentialModule, only : getOldSphPotr => getSphPotr
   use PotentialModule, only : getOldVdif => getVdif
   use PotentialModule, only : setOldVdif => setVdif
!
   use AtomModule, only : getLocalNumSpecies
!
   use PotentialGenerationModule, only : getNewSphPotr => getSphPotr
   use PotentialGenerationModule, only : setNewVdif => setVdif
!
   use ChargeDensityModule, only : getSphChargeDensity, &
                                   getSphMomentDensity
!
   use PublicTypeDefinitionsModule, only : MixListRealStruct
!
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms,n_spin_pola
   integer (kind=IntKind) :: id, ia, nr, is, data_size, num_items, p_item
!
   real (kind=RealKind), pointer :: pStore_old(:,:), pStore_new(:,:)
   real (kind=RealKind), pointer :: ptmp1(:), ptmp2(:)
!
   type (MixListRealStruct), target :: RealArrayList
   type (MixListRealStruct), pointer :: p_RAL
!
   data_size = 0; num_items = 0
   do id = 1,LocalNumAtoms
      nr = getNumRmesh(id)+1
      data_size = max(data_size,nr)
      num_items = num_items + getLocalNumSpecies(id)
   enddo
   data_size = data_size*n_spin_pola
   if (.not.isDataStorageExisting("MixingVectorOld")) then
      call ErrorHandler('updateMixRealValues','MixingVectorOld not defined')
   endif
   if (.not.isDataStorageExisting("MixingVectorNew")) then
      call ErrorHandler('updateMixRealValues','MixingVectorNew not defined')
   endif
!
   pStore_old => getDataStorage("MixingVectorOld",data_size,num_items,RealMark)
   pStore_new => getDataStorage("MixingVectorNew",data_size,num_items,RealMark)
   p_RAL => RealArrayList
   p_item = 0
   do id = 1, LocalNumAtoms
      nr = getNumRmesh(id)+1
      do ia = 1, getLocalNumSpecies(id)
         p_item = p_item + 1
         p_RAL%vector_old => pStore_old(1:nr*n_spin_pola,p_item)
         p_RAL%vector_new => pStore_new(1:nr*n_spin_pola,p_item)
         do is = 1, n_spin_pola
            if (isPotentialMixing()) then
               ptmp2 => getNewSphPotr(id,ia,is)
               ptmp2 = p_RAL%vector_new(nr*(is-1)+1:nr*is-1)
!              call setNewVdif(p_RAL%vector_new(nr*is))
!              if (id == 1) then
!                 write(6,*)'Mix Vr: spin =',is
!                 call writeMatrix('Mix Vr',ptmp2,nr-1)
!              endif
               if (is == 1) then
                  call setOldVdif(p_RAL%vector_new(nr))
                  call setNewVdif(p_RAL%vector_new(nr))
               endif
            else
               ptmp1 => getSphChargeDensity('TotalNew', id, ia)
               if ( n_spin_pola==2 ) then
                  ptmp2 => getSphMomentDensity('TotalNew', id, ia)
               endif
               if ( n_spin_pola==2 ) then
                  if ( is == 1 ) then
                     ptmp1 = p_RAL%vector_new(nr*(is-1)+1:nr*is)
                  else
                     ptmp1 = ptmp1+p_RAL%vector_new(nr*(is-1)+1:nr*is)
                  endif
                  if ( is == 1 ) then
                     ptmp2 = p_RAL%vector_new(nr*(is-1)+1:nr*is)
                  else
                     ptmp2 = ptmp2-p_RAL%vector_new(nr*(is-1)+1:nr*is)
                  endif
               else
                  ptmp1 = p_RAL%vector_new(nr*(is-1)+1:nr*is)
               endif
            endif
         enddo
         p_RAL => p_RAL%next
      enddo
   enddo
!
   end subroutine updateMixRealValues
