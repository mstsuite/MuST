program tst_Data
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use DataServiceCenterModule, only : initDataServiceCenter, endDataServiceCenter
   use DataServiceCenterModule, only : RealType, IntegerType, ComplexType
   use DataServiceCenterModule, only : RealMark, IntegerMark, ComplexMark
   use DataServiceCenterModule, only : createDataStorage, getDataStorage
   use DataServiceCenterModule, only : getDataStorageType, getDataStorageSize
!
   implicit none
!
   integer (kind=IntKind), pointer :: num_atoms
   integer (kind=IntKind), pointer :: NumAtoms
   integer (kind=IntKind) :: i,ir
!
   real (kind=RealKind), pointer :: pot1(:,:)
   real (kind=RealKind), pointer :: pot2(:,:)
   real (kind=RealKind), pointer :: pot3(:)
   real (kind=RealKind), pointer :: pot4(:)
!
   call initDataServiceCenter()
!
   call createDataStorage('NumberOfAtoms',1,IntegerType)
   call createDataStorage('Potential',50,RealType)
!
!  ===================================================================
!
   num_atoms => getDataStorage('NumberOfAtoms',IntegerMark)
   num_atoms = 5
!
   pot1 => getDataStorage("Potential",10,5,RealMark)
   do i=1,num_atoms
      do ir = 1,10
          pot1(ir,i) = 1.0d0*ir+10.0d0*i
      enddo
   enddo
   nullify( pot1 )
!
!  ===================================================================
!
   NumAtoms => getDataStorage('NumberOfAtoms',IntegerMark)
   print *,'NumAtoms = ',NumAtoms
!
!  pot3 => getDataStorage("Potential",50,RealMark)
!  do i=1,num_atoms
!     do ir = 1,10
!         pot3(ir+(i-1)*10) = 1.0d0*ir+10.0d0*i
!     enddo
!  enddo
!  nullify( pot3 )
!
   pot2 => getDataStorage("Potential",10,5,RealMark)
   do i=1,num_atoms
      do ir = 1,10
          print *,'ir,i,pot2(ir,i)',ir,', ',i,', ',pot2(ir,i)
      enddo
   enddo
   nullify( pot2 )
!
!  pot4 => getDataStorage("Potential",50,RealMark)
!  do i=1,num_atoms
!     do ir = 1,10
!         print *,'ir,i,pot4(ir,i)',ir,', ',i,', ',pot4(ir+(i-1)*10)
!     enddo
!  enddo
!
!  ===================================================================
!
   call endDataServiceCenter()
!
end program tst_Data
