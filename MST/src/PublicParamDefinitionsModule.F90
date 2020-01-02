Module PublicParamDefinitionsModule
   use KindParamModule, only : IntKind
!
public :: checkParameter,  &
          getParameterInfo
!
   character (len=12) :: StandardInputFile= 'i_mst_stdin'
   integer (kind=IntKind), parameter :: MaxLenFileName = 100
!
!  ===================================================================
!  E-contour  type parameters:
!  ===================================================================
   character (len=*), parameter :: ec_pname = 'E_contour'
   integer (kind=IntKind), parameter :: HalfCircle=0
   integer (kind=IntKind), parameter :: RectBox=1
   integer (kind=IntKind), parameter :: HorizLine=2
   integer (kind=IntKind), parameter :: VertLine=3
   integer (kind=IntKind), parameter :: ButterFly=4
   integer (kind=IntKind), parameter :: MatsubaraPoles=5
   integer (kind=IntKind), parameter :: ec_min = 0
   integer (kind=IntKind), parameter :: ec_max = 5
   character (len=20), parameter :: &
                       ec_info(ec_min:ec_max) = (/'Half circle    ',  &
                                                  'Rectangluar box',  &
                                                  'Horizontal line',  &
                                                  'Vertical line  ',  &
                                                  'ButterFly shape',  &
                                                  'Matsubara poles'/)
!
!  ===================================================================
!  E-grid  type parameters:
!  ===================================================================
   character (len=*), parameter :: eg_pname = 'E_grid'
   integer (kind=IntKind), parameter :: EqualInterval=0
   integer (kind=IntKind), parameter :: GaussianPoints=1
   integer (kind=IntKind), parameter :: LogInterval=2
   integer (kind=IntKind), parameter :: NicholsonPoints=3
   integer (kind=IntKind), parameter :: eg_min = 0
   integer (kind=IntKind), parameter :: eg_max = 3
   character (len=20), parameter :: &
                       eg_info(eg_min:eg_max) = (/'Uniform mesh       ', &
                                                  'Gaussian point mesh', &
                                                  'Logrithm mesh      ', &
                                                  'Nicholson pole mesh'/)
!
!  ===================================================================
!  Print DOS instruction parameters:
!  ===================================================================
   character (len=*), parameter :: pd_pname = 'Print_DOS'
   integer (kind=IntKind), parameter :: PrintDOSswitchOff=0
   integer (kind=IntKind), parameter :: PrintDOSswitchOn=1 ! The input value > 0 to specify an atom
   integer (kind=IntKind), parameter :: PrintEachAtomDOS=-1
   integer (kind=IntKind), parameter :: pd_min = -1
   integer (kind=IntKind), parameter :: pd_max = 1
   character (len=80), parameter :: &
      pd_info(pd_min:pd_max) = (/'Calculate DOS and print the DOS of each atom and system average       ', &
                                 'DOS is not calculated and printed                                     ', &
                                 'Calculate DOS and print the DOS of a specified atom and system average'/)
!
!  ===================================================================
!  Potential type parameters:                
!  ===================================================================
   character (len=*), parameter :: pt_pname = 'Pot_type'
   integer (kind=IntKind), parameter :: MuffinTin     = 0
   integer (kind=IntKind), parameter :: ASA           = 1
   integer (kind=IntKind), parameter :: MuffinTinASA  = 2
   integer (kind=IntKind), parameter :: FullPot       = 3
   integer (kind=IntKind), parameter :: MuffinTinTest = 4
   integer (kind=IntKind), parameter :: EmptyLattice  = 5
   integer (kind=IntKind), parameter :: Mathieu       = 6
   integer (kind=IntKind), parameter :: Coulomb       = 7
   integer (kind=IntKind), parameter :: MuffinTinFullPot = 8
   integer (kind=IntKind), parameter :: pt_min = 0
   integer (kind=IntKind), parameter :: pt_max = 8
   character (len=50), parameter :: &
                       pt_info(pt_min:pt_max) = (/'Muffin-Tin                               ', &
                                                  'ASA                                      ', &
                                                  'Muffin-Tin ASA                           ', &
                                                  'Full-Potential                           ', &
                                                  'Muffin-Tin Potential Test                ', &
                                                  'Full-Potential and Empty Lattice Test    ', &
                                                  'Full-Potential and Mathieu Potential Test', &
                                                  'Coulomb Potential Test                   ', &
                                                  'Full-Potential inside Muffin-Tin         '/)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function checkParameter(pname,pvalue) result(t)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   implicit none
!
   character (len=*), intent(in) :: pname
!
   integer (kind=IntKind), intent(in) :: pvalue
!
   logical :: t
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   if ( nocaseCompare(pname,ec_pname) ) then
      if (pvalue < ec_min .or. pvalue > ec_max) then
         call ErrorHandler(pname,'Undefined parameter value: ',pvalue)
      endif
   else if ( nocaseCompare(pname,eg_pname) ) then
      if (pvalue < eg_min .or. pvalue > eg_max) then
         call ErrorHandler(pname,'Undefined parameter value: ',pvalue)
      endif
   else if ( nocaseCompare(pname,pd_pname) ) then
      if (pvalue < pd_min .or. pvalue > pd_max) then
         call ErrorHandler(pname,'Undefined parameter value: ',pvalue)
      endif
   else if ( nocaseCompare(pname,pt_pname) ) then
      if (pvalue < pt_min .or. pvalue > pt_max) then
         call ErrorHandler(pname,'Undefined parameter value: ',pvalue)
      endif
   else
      call ErrorHandler('checkParameter','Undefined parameter type: ',pvalue)
   endif
!
   t = .true.
!
   end function checkParameter
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getParameterInfo(pname,pvalue) result(s)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   implicit none
!
   character (len=*), intent(in) :: pname
!
   integer (kind=IntKind), intent(in) :: pvalue
!
   character (len=80) :: s
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   s = ' '
!
   if ( nocaseCompare(pname,ec_pname) ) then
      s = ec_info(pvalue)
   else if ( nocaseCompare(pname,eg_pname) ) then
      s = eg_info(pvalue)
   else if ( nocaseCompare(pname,pd_pname) ) then
      s = pd_info(pvalue)
   else if ( nocaseCompare(pname,pt_pname) ) then
      s = pt_info(pvalue)
   else
      call ErrorHandler('getParameterInfo','Undefined parameter type: ',pname)
   endif
!
   end function getParameterInfo
end Module PublicParamDefinitionsModule
