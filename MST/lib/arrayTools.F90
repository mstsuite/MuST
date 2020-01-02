!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasCmplx2Real(a,n) result(b)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : WarningHandler, ErrorHandler
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: m
!
   real (kind=RealKind), target :: a(:)
   real (kind=RealKind), pointer :: b(:)
!
   m = 2*n
   if (size(a) < m) then
      call ErrorHandler('aliasCmplx2Real','The source array size is less that n',n)
   endif
   b => a(1:m)
!
   end function aliasCmplx2Real
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasReal2Cmplx(a,n) result(b)
!  ===================================================================
   use KindParamModule, only : IntKind, CmplxKind
   use ErrorHandlerModule, only : WarningHandler, ErrorHandler
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
   integer (kind=IntKind) :: m
!
   complex (kind=CmplxKind), target :: a(:)
   complex (kind=CmplxKind), pointer :: b(:)
!
   m = n/2
   if (size(a) < m) then
      call ErrorHandler('aliasReal2Cmplx','The source array size is less that n',n)
   else if (mod(n,2) /= 0) then
      call WarningHandler('aliasReal2Cmplx','The source array size is not an even number',n)
   endif
   b => a(1:m)
!
   end function aliasReal2Cmplx
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine copyReal2Cmplx(a,b,n)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
!
   real(kind=RealKind), target :: a(2*((n-1)/2+1))
   complex(kind=CmplxKind), intent(out) :: b(n)
!
   b = transfer(a,b)
!
   end subroutine copyReal2Cmplx
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine copyCmplx2Real(a,b,n)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
!
   complex(kind=CmplxKind), target :: a(n/2)
   real(kind=RealKind), intent(out) :: b(n)
!
   if (mod(n,2) == 0) then
      b = transfer(a,b)
   else
      call ErrorHandler('copyCmplx2Real','The target array size is not an even number',n)
   endif
!
   end subroutine copyCmplx2Real
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray1_l(array,n1)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1
   logical, target  :: array(n1)
!
   logical, pointer :: parray(:)
!
   parray => array(1:n1)
!
   end function aliasArray1_l
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray2_l(array,n1,n2)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2
   logical, target  :: array(n1,n2)
!
   logical, pointer :: parray(:,:)
!
   parray => array(1:n1,1:n2)
!
   end function aliasArray2_l
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray3_l(array,n1,n2,n3)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3
   logical, target  :: array(n1,n2,n3)
!
   logical, pointer :: parray(:,:,:)
!
   parray => array(1:n1,1:n2,1:n3)
!
   end function aliasArray3_l
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray4_l(array,n1,n2,n3,n4)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3,n4
   logical, target  :: array(n1,n2,n3,n4)
!
   logical, pointer :: parray(:,:,:,:)
!
   parray => array(1:n1,1:n2,1:n3,1:n4)
!
   end function aliasArray4_l
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray1_i(array,n1)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1
   integer(kind=IntKind), target  :: array(n1)
!
   integer(kind=IntKind), pointer :: parray(:)
!
   parray => array(1:n1)
!
   end function aliasArray1_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray2_i(array,n1,n2)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2
   integer(kind=IntKind), target  :: array(n1,n2)
!
   integer(kind=IntKind), pointer :: parray(:,:)
!
   parray => array(1:n1,1:n2)
!
   end function aliasArray2_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray3_i(array,n1,n2,n3)       result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3
   integer(kind=IntKind), target  :: array(n1,n2,n3)
!
   integer(kind=IntKind), pointer :: parray(:,:,:)
!
   parray => array(1:n1,1:n2,1:n3)
!
   end function aliasArray3_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray4_i(array,n1,n2,n3,n4)       result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3,n4
   integer(kind=IntKind), target  :: array(n1,n2,n3,n4)
!
   integer(kind=IntKind), pointer :: parray(:,:,:,:)
!
   parray => array(1:n1,1:n2,1:n3,1:n4)
!
   end function aliasArray4_i
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray1_r(array,n1)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1
   real(kind=RealKind), target  :: array(n1)
!
   real(kind=RealKind), pointer :: parray(:)
!
   parray => array(1:n1)
!
   end function aliasArray1_r
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray2_r(array,n1,n2)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2
   real(kind=RealKind), target  :: array(n1,n2)
!
   real(kind=RealKind), pointer :: parray(:,:)
!
   parray => array(1:n1,1:n2)
!
   end function aliasArray2_r
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray3_r(array,n1,n2,n3)       result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3
   real(kind=RealKind), target  :: array(n1,n2,n3)
!
   real(kind=RealKind), pointer :: parray(:,:,:)
!
   parray => array(1:n1,1:n2,1:n3)
!
   end function aliasArray3_r
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray4_r(array,n1,n2,n3,n4)       result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3,n4
   real(kind=RealKind), target  :: array(n1,n2,n3,n4)
!
   real(kind=RealKind), pointer :: parray(:,:,:,:)
!
   parray => array(1:n1,1:n2,1:n3,1:n4)
!
   end function aliasArray4_r
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray1_c(array,n1)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, CmplxKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1
   complex(kind=CmplxKind), target :: array(n1)
!
   complex(kind=CmplxKind), pointer :: parray(:)
!
   parray => array(1:n1)
!
   end function aliasArray1_c
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray2_c(array,n1,n2)              result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, CmplxKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2
   complex(kind=CmplxKind), target :: array(n1,n2)
!
   complex(kind=CmplxKind), pointer :: parray(:,:)
!
   parray => array(1:n1,1:n2)
!
   end function aliasArray2_c
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray3_c(array,n1,n2,n3)       result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, CmplxKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3
   complex(kind=CmplxKind), target :: array(n1,n2,n3)
!
   complex(kind=CmplxKind), pointer :: parray(:,:,:)
!
   parray => array(1:n1,1:n2,1:n3)
!
   end function aliasArray3_c
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray4_c(array,n1,n2,n3,n4)       result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind, CmplxKind
   implicit none
!
   integer(kind=IntKind), intent(in) :: n1,n2,n3,n4
   complex(kind=CmplxKind), target :: array(n1,n2,n3,n4)
!
   complex(kind=CmplxKind), pointer :: parray(:,:,:,:)
!
   parray => array(1:n1,1:n2,1:n3,1:n4)
!
   end function aliasArray4_c
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray1_s(array,n1)                    result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: n1
   character (len=1), intent(in), target :: array(n1)
!
   character (len=1), pointer :: parray(:)
!
   parray => array(1:n1)
!
   end function aliasArray1_s
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray2_s(array,n1,n2)                result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: n1, n2
!
   character (len=1), intent(in), target :: array(n1,n2)
!
   character (len=1), pointer :: parray(:,:)
!
   parray => array(1:n1,1:n2)
!
   end function aliasArray2_s
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray3_s(array,n1,n2,n3)             result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: n1, n2, n3
!
   character (len=1), intent(in), target :: array(n1,n2,n3)
!
   character (len=1), pointer :: parray(:,:,:)
!
   parray => array(1:n1,1:n2,1:n3)
!
   end function aliasArray3_s
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function aliasArray4_s(array,n1,n2,n3,n4)             result(parray)
!  ===================================================================
   use KindParamModule, only : IntKind
   implicit none
!
   integer (kind=IntKind), intent(in) :: n1, n2, n3, n4
!
   character (len=1), intent(in), target :: array(n1,n2,n3,n4)
!
   character (len=1), pointer :: parray(:,:,:,:)
!
   parray => array(1:n1,1:n2,1:n3,1:n4)
!
   end function aliasArray4_s
!  ===================================================================
