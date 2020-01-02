module PotentialTypeModule
   use KindParamModule, only : IntKind
   use ErrorHandlerModule, only : ErrorHandler
!  
public :: initPotentialType,         &
          endPotentialType,          &
          isMuffinTinPotential,      &
          isASAPotential,            & 
          isMuffinTinASAPotential,   &
          isSphericalPotential,      &
          setIsSphericalPotential,   &
          isFullPotential,           &
          isTestPotential,           &
          isEmptyLatticePotential,   &
          isMathieuPotential,        &
          isCoulombPotential,        &
          isMuffinTinTestPotential,  &
          isMuffinTinFullPotential,  &
          printPotentialType
!
private
!
   logical :: Initialized = .false.
   logical :: isSphPot = .false.
!
   integer (kind=IntKind) :: pot_type
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initPotentialType(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : checkParameter
   implicit none
!
   integer (kind=IntKind), intent(in) :: pt
!
   if (checkParameter('Pot_Type',pt)) then
      pot_type = pt
   endif
   Initialized = .true.
!
   end subroutine initPotentialType
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endPotentialType()
!  ===================================================================
   implicit none
!
   pot_type = -1
   Initialized = .false.
!
   end subroutine endPotentialType
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printPotentialType(fu_in)
!  ===================================================================
   use PublicParamDefinitionsModule, only : getParameterInfo
   implicit none
!
   integer(kind=IntKind), optional, intent(in) :: fu_in
   integer(kind=IntKind) :: fu
!
   if (present(fu_in)) then
      fu = fu_in
   else
      fu = 6
      write(fu,'(/,80(''-''))')
      write(fu,'(/,21x,a)')'**************************************'
      write(fu,'( 21x,a )')'*   Output from printPotentialType   *'
      write(fu,'(21x,a,/)')'**************************************'
      write(fu,'(/,80(''=''))')
   endif
!
   write(fu,'(''# Potential Type          :'',$)')
   write(fu,*)trim(getParameterInfo('Pot_Type',pot_type))
   if (fu==6) then
      write(fu,'(/,80(''=''))')
   endif
!
   end subroutine printPotentialType
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isSphericalPotential() result(pt)
!  ===================================================================
   implicit none
   logical :: pt
!
   if (.not.Initialized) then
      call ErrorHandler('isSphericalPotential','module is not initialized')
   endif
!
   pt = isSphPot
!
   end function isSphericalPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setIsSphericalPotential(pt)
!  ===================================================================
   implicit none
!
   logical, intent(in) :: pt
!
   if (.not.Initialized) then
      call ErrorHandler('isSphericalPotential','module is not initialized')
   endif
!
   isSphPot = pt
!
   end subroutine setIsSphericalPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isMuffinTinPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : MuffinTin, MuffinTinASA
   implicit none
   logical :: pt
!
   if (.not.Initialized) then
      call ErrorHandler('isMuffinTinPotential','module is not initialized')
   endif
!
   if (pot_type == MuffinTin .or. pot_type == MuffinTinASA) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isMuffinTinPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isASAPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : ASA
   implicit none
   logical :: pt
!
   if (.not.Initialized) then
      call ErrorHandler('isASAPotential','module is not initialized')
   endif
!
   if (pot_type == ASA) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isASAPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isMuffinTinASAPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : MuffinTinASA
   implicit none
   logical :: pt
!
   if (.not.Initialized) then
      call ErrorHandler('isMuffinTinASAPotential','module is not initialized')
   endif
!
   if (pot_type == MuffinTinASA) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isMuffinTinASAPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isFullPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : FullPot, EmptyLattice,    &
                                            Mathieu, MuffinTinFullPot
   implicit none
   logical :: pt
!
   if (.not.Initialized) then
      call ErrorHandler('isFullPotential','module is not initialized')
   endif
!
   if (pot_type == FullPot .or. pot_type == EmptyLattice .or.         &
       pot_type == Mathieu .or. pot_type == MuffinTinFullPot) then
       pt = .true.
   else
       pt = .false.
   endif
   end function isFullPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isTestPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : EmptyLattice, Mathieu,    &
                                            MuffinTinTest, Coulomb
   implicit none
   logical :: pt
!
   if (pot_type == EmptyLattice) then
       pt = .true.
   else if (pot_type == Mathieu) then
       pt = .true.
   else if (pot_type == MuffinTinTest) then
       pt = .true.
   else if (pot_type == Coulomb) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isTestPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isEmptyLatticePotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : EmptyLattice
   implicit none
   logical :: pt
!
   if (pot_type == EmptyLattice) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isEmptyLatticePotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isMathieuPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : Mathieu
   implicit none
   logical :: pt
!
   if (pot_type == Mathieu) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isMathieuPotential
!  ===================================================================
!
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isCoulombPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : Coulomb
   implicit none
   logical :: pt
!
   if (pot_type == Coulomb) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isCoulombPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isMuffinTinTestPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : MuffinTinTest
   implicit none
   logical :: pt
!
   if (pot_type == MuffinTinTest) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isMuffinTinTestPotential
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isMuffinTinFullPotential() result(pt)
!  ===================================================================
   use PublicParamDefinitionsModule, only : MuffinTinFullPot
   implicit none
   logical :: pt
!
   if (pot_type == MuffinTinFullPot) then
       pt = .true.
   else
       pt = .false.
   endif
!
   end function isMuffinTinFullPotential
!  ===================================================================
end module PotentialTypeModule
