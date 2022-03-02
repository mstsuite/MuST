module KuboDataModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicParamDefinitionsModule, only : MaxLenFileName
!
public :: initKuboData,                 &
          endKuboData,                  &
          useStepFunctionForSigma,     &
          useCubicSymmetryForSigma,    &
          includeVertexCorrections,    &
          printTildeMatrices,          &
          printKuboData
!
public
   integer (kind=IntKind) :: TableID

!  Conductivity Parameters
   integer (kind=IntKind), private :: use_sf = 1
   integer (kind=IntKind), private :: use_csymm = 1
   integer (kind=IntKind), private :: is_ef_rp = 0
   integer (kind=IntKind), private :: print_tilde = 0
   integer (kind=IntKind), private :: vertex_corr = 1
   real (kind=RealKind), private :: imag_part = 0.001
   real (kind=RealKind), private :: sigma_real_part = 0.5
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initKuboData(tbl_id)
!  ===================================================================
   use MathParamModule, only : TEN2m6, TEN2m8, ZERO, TEN2m4, ONE
   use PublicParamDefinitionsModule, only : ASA, MuffinTin, MuffinTinASA
   use InputModule, only : getKeyValue
   implicit none
!
   integer (kind=IntKind), intent(in) :: tbl_id
!
   integer (kind=IntKind) :: rstatus, n
   character(len=120) :: svalue
!
   character (len=50) :: s50
!
   real (kind=RealKind) :: rp(3)
!
   TableID = tbl_id
!
   rstatus = getKeyValue(tbl_id,'Fermi Energy Imaginary Part',imag_part)
   rstatus = getKeyValue(tbl_id,'Use Different Fermi Energy', is_ef_rp)
   rstatus = getKeyValue(tbl_id,'Fermi Energy Real Part', sigma_real_part)
   rstatus = getKeyValue(tbl_id,'Integrate Upto Muffin Tin', use_sf)
   rstatus = getKeyValue(tbl_id,'Use Cubic Symmetry', use_csymm)
   rstatus = getKeyValue(tbl_id,'Vertex Corrections', vertex_corr)
   rstatus = getKeyValue(tbl_id,'Print Tilde Matrices', print_tilde)

   end subroutine initKuboData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endKuboData()
!  ===================================================================
   implicit none
   
   end subroutine endKuboData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printKuboData(fu_in)
!  ===================================================================
   use PublicParamDefinitionsModule, only : MuffinTin, ASA, MuffinTinASA, &
                                            MuffinTinTest
   implicit none
!
   integer (kind=IntKind), optional :: fu_in
   integer (kind=IntKind) :: fu
!
   if (present(fu_in)) then
      fu = fu_in
   else
      fu = 6
   endif
!
   end subroutine printKuboData
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function useStepFunctionForSigma() result(md)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: md

   md = use_sf
   
   end function useStepFunctionForSigma
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function useCubicSymmetryForSigma() result(md)
!  ===================================================================
   implicit none
   logical :: md

   if (use_csymm == 0) then
     md = .false.
   else
     md = .true.
   endif
   
   end function useCubicSymmetryForSigma
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function includeVertexCorrections() result(md)
!  ===================================================================
   implicit none
   logical :: md

   if (vertex_corr == 0) then
     md = .false.
   else
     md = .true.
   endif

   end function includeVertexCorrections
!  ===================================================================
!
!  *******************************************************************
!  
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
   function getFermiEnergyImagPart() result(del)
!  ===================================================================
   implicit none

   real (kind=RealKind) :: del

   del = imag_part

   end function getFermiEnergyImagPart
!  ===================================================================
!
!  *******************************************************************
!  
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
   function isFermiEnergyRealPart() result(use_ef_rp)
!  ===================================================================
   implicit none

   logical :: use_ef_rp


   if (is_ef_rp == 0) then
     use_ef_rp = .false.
   else if (is_ef_rp == 1) then
     use_ef_rp = .true.
   endif   

   end function isFermiEnergyRealPart
!  ===================================================================
!  
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
   function getFermiEnergyRealPart() result(alt_ef)
!  ===================================================================
   implicit none

   real (kind=RealKind) :: alt_ef

   alt_ef = sigma_real_part

   end function getFermiEnergyRealPart
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function printTildeMatrices() result(print_mat)
!  ===================================================================
   implicit none

   logical :: print_mat

   if (print_tilde == 1) then
     print_mat = .true.
   else if (print_tilde == 0) then
     print_mat = .false.
   else
     call ErrorHandler('printTildeMatrices', &
             'Incorrect Value (choose 0 or 1)', print_mat)
   endif

   end function printTildeMatrices
!  ==================================================================
end module KuboDataModule
