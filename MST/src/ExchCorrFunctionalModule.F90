module ExchCorrFunctionalModule
!
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, TWO, THREE, THIRD, PI4, TEN2m6, &
                               TEN2m8, TEN2m10, HALF
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler,        &
                                  StopHandler
!
#ifdef LIBXC
   use xc_f90_types_m
   use xc_f90_lib_m
!  use xc_f03_lib_m
! #define xc_f90_pointer_t xc_f03_func_t
! #define xc_f90_func_init xc_f03_func_init
! #define xc_f90_info_number xc_f03_info_number
! #define xc_f90_info_family xc_f03_info_family
! #define xc_f90_functional_get_number xc_f03_functional_get_number
#endif
!
public :: initExchCorrFunctional,         &
          endExchCorrFunctional,          &
          isLDAFunctional,                &
          isGGAFunctional,                &
          isMGGAFunctional,               &
          isHybridFunctional,             &
          getExchCorrPot,                 &
          getExchCorrEnDen,               &
          calSphExchangeCorrelation,      &
          calExchangeCorrelation
!
   interface getExchCorrPot
      module procedure getExchCorrPot_s, getExchCorrPot_v
   end interface
!
   interface getExchCorrEnDen
      module procedure getExchCorrEnDen_s, getExchCorrEnDen_v
   end interface
!
   interface calSphExchangeCorrelation
      module procedure calSphExchangeCorrelation_s, calSphExchangeCorrelation_v
   end interface
!
   interface calExchangeCorrelation
      module procedure calExchangeCorrelation_s, calExchangeCorrelation_v
   end interface
!
private 
   logical :: Initialized = .false.
!
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: LegacyFunctionalID
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: FunctionalID_X
   integer (kind=IntKind) :: FunctionalID_C
   integer (kind=IntKind) :: FunctionalID_XC
!
   character (len=50) :: FunctionalName_X
   character (len=50) :: FunctionalName_C
   character (len=50) :: FunctionalName_XC
   character (len=15) :: FunctionalType
   character (len=3), parameter :: GGA = 'GGA'
   character (len=3), parameter :: LDA = 'LDA'
   character (len=10), parameter :: HybridGGA = 'Hybrid GGA'
   character (len=11), parameter :: HybridMGGA = 'Hybrid MGGA'
   character (len=4), parameter :: MGGA = 'MGGA'
   character (len=6), parameter :: UnkownXC = 'Unkown'
   character (len=6), parameter :: Null_String = ' '
!
   real (kind=RealKind) ::    Vexc_s(2)
   real (kind=RealKind) ::    Eexc_s
   real (kind=RealKind), allocatable, target :: Vexc_v(:,:)
   real (kind=RealKind), allocatable, target :: Eexc_v(:)
!
   integer (kind=IntKind) :: NumFunctionals = 0
   integer (kind=IntKind) :: nV, nE
!
   logical :: Exchange = .false.
   logical :: Correlation = .false.
   logical :: ECin1 = .false.
!
#ifdef LIBXC
   TYPE(xc_f90_pointer_t), target :: xc_func_1
   TYPE(xc_f90_pointer_t), target :: xc_func_2
   TYPE(xc_f90_pointer_t), target :: xc_info_1
   TYPE(xc_f90_pointer_t), target :: xc_info_2
   real (kind=RealKind), parameter :: energy_units_conv = TWO ! From Hartree to Rydberg
#else
   real (kind=RealKind), parameter :: energy_units_conv = ONE
#endif
!
contains
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initExchCorrFunctional( pola, excorr_name, iprint )
!  ===================================================================
   use StringModule, only : initString, endString, getNumTokens, readToken
!
   implicit none
!
   character (len=*), intent(in) :: excorr_name
   integer (kind=IntKind), intent(in) :: pola, iprint
#ifdef LIBXC
   integer (kind=IntKind) :: vmajor, vminor, vmicro, fid1, fid2, i, n
   character (len=40) :: fname1, fname2, fname
   character (len=120) :: func_name, func_kind, func_family, ref
   TYPE(xc_f90_pointer_t), pointer :: xc_func_p
   TYPE(xc_f90_pointer_t), pointer :: xc_info_p
#endif
!
   interface
      function isNumber(s) result(t)
         character (len=*), intent(in) :: s
         logical :: t
      end function isNumber
   end interface
!
   if (pola < 1 .or. pola > 2) then
      call ErrorHandler( "initExchCorrFunctional", "Invalid spin index", pola )
   endif
!
   n_spin_pola = pola
   print_level = iprint
!
   FunctionalType = Null_String
   FunctionalName_X  = Null_String
   FunctionalName_C  = Null_String
   FunctionalName_XC = Null_String
   FunctionalID_X = -1
   FunctionalID_C = -1
   FunctionalID_XC = -1
!
#ifdef LIBXC
!  -------------------------------------------------------------------
   call xc_f90_version(vmajor, vminor, vmicro)
!  -------------------------------------------------------------------
   if (print_level >= 0) then
      write(6,'(/,3x,68("*"))')
      write(6,'(3x,a,66x,a)')"*","*"
      write(6,'(3x,a,3x,a,3x,a)')"*",                                 &
         "LibXC is used for evaluating exchange-correlation functional","*"
      write(6,'(3x,a,22x,a,i1,".",i1,".",i1,23x,a)')"*",              &
         "LibXC version:  ", vmajor, vminor, vmicro, "*"
      write(6,'(3x,a,66x,a)')"*","*"
      write(6,'(3x,68("*"))')
   endif
#endif
!
   if (isNumber(excorr_name)) then
      read(excorr_name,*) LegacyFunctionalID
   else
      LegacyFunctionalID = -1
   endif
!
   if ( LegacyFunctionalID == 0 ) then
      FunctionalName_X = 'LDA_X'
      FunctionalName_C = 'LDA_C_VBH'
      FunctionalType = LDA
      NumFunctionals = 2
#ifdef LIBXC
      FunctionalID_X   = XC_LDA_X
      FunctionalID_C   = XC_LDA_C_vBH
#else
      FunctionalID_X   = 1
      FunctionalID_C   = 17
#endif
   else if ( LegacyFunctionalID == 1 ) then
      FunctionalName_X = 'LDA_X'
      FunctionalName_C = 'LDA_C_VWN'
      FunctionalType = LDA
      NumFunctionals = 2
#ifdef LIBXC
      FunctionalID_X   = XC_LDA_X
      FunctionalID_C   = XC_LDA_C_VWN_1
#else
      FunctionalID_X   = 1
      FunctionalID_X   = 28
#endif
   else if ( LegacyFunctionalID == 2 ) then
      FunctionalName_X = 'LDA_X'
      FunctionalName_C = 'LDA_C_PZ'
      FunctionalType = LDA
      NumFunctionals = 2
#ifdef LIBXC
      FunctionalID_X  = XC_LDA_X
      FunctionalID_C  = XC_LDA_C_PZ
#else
      FunctionalID_X   = 1
      FunctionalID_C   = 9
      call ErrorHandler( "initExchCorrFunctional", "This functional requires enabling LibXC", LegacyFunctionalID )
#endif
   else if ( LegacyFunctionalID == 3 ) then
      FunctionalName_X = 'GGA_X_PW91'
      FunctionalName_C = 'GGA_C_PW91'
      FunctionalType = GGA
      NumFunctionals = 2
#ifdef LIBXC
      FunctionalID_X  = XC_GGA_X_PW91
      FunctionalID_C  = XC_GGA_C_PW91
#else
      FunctionalID_X  = 109
      FunctionalID_C  = 134
      call ErrorHandler( "initExchCorrFunctional", "This functional requires enabling LibXC", LegacyFunctionalID )
#endif
   else if ( LegacyFunctionalID == 4 ) then
      FunctionalName_X = 'GGA_X_PBE'
      FunctionalName_C = 'GGA_C_PBE'
      FunctionalType = GGA
      NumFunctionals = 2
#ifdef LIBXC
      FunctionalID_X   = XC_GGA_X_PBE
      FunctionalID_C   = XC_GGA_C_PBE
#else
      FunctionalID_X   = 101
      FunctionalID_C   = 130
      call ErrorHandler( "initExchCorrFunctional", "This functional requires enabling LibXC", LegacyFunctionalID )
#endif
   else
#ifdef LIBXC
      call initString(excorr_name,separator="+")
      NumFunctionals = getNumTokens()
      if (NumFunctionals == 1) then
         call readToken(1,fname1)
      else if (NumFunctionals == 2) then
         call readToken(1,fname1)
         call readToken(2,fname2)
      else
         call ErrorHandler( "initExchCorrFunctional", "Invalid functional name(s)", excorr_name )
      endif
      call endString()
#else
      call ErrorHandler( "initExchCorrFunctional", "Functional ID is invalid or LibXC is required", LegacyFunctionalID )
#endif
   endif
!
#ifdef LIBXC
   if ( LegacyFunctionalID >= 0 ) then
      if (NumFunctionals == 1) then
         fname1 = FunctionalName_XC
      else
         fname1 = FunctionalName_X
         fname2 = FunctionalName_C
      endif
   endif
!
   if (NumFunctionals == 1) then
      fid1 = xc_f90_functional_get_number(fname1)
      if (n_spin_pola == 1) then
!        -------------------------------------------------------------
         call xc_f90_func_init(xc_func_1, xc_info_1, fid1, XC_UNPOLARIZED)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call xc_f90_func_init(xc_func_1, xc_info_1, fid1, XC_POLARIZED)
!        -------------------------------------------------------------
      endif
      if (xc_f90_info_kind(xc_info_1) == XC_EXCHANGE_CORRELATION) then
         call xc_f90_info_name(xc_info_1,FunctionalName_XC)
         FunctionalID_XC = xc_f90_info_number(xc_info_1)
         ECin1 = .true.
      else
         call ErrorHandler( "initExchCorrFunctional",                 &
                            "The functional does not include both exchange and correlation", fname1 )
      endif
   else
      fid1 = xc_f90_functional_get_number(fname1)
      fid2 = xc_f90_functional_get_number(fname2)
      if (n_spin_pola == 1) then
!        -------------------------------------------------------------
         call xc_f90_func_init(xc_func_1, xc_info_1, fid1, XC_UNPOLARIZED)
!        -------------------------------------------------------------
         call xc_f90_func_init(xc_func_2, xc_info_2, fid2, XC_UNPOLARIZED)
!        -------------------------------------------------------------
      else
!        -------------------------------------------------------------
         call xc_f90_func_init(xc_func_1, xc_info_1, fid1, XC_POLARIZED)
!        -------------------------------------------------------------
         call xc_f90_func_init(xc_func_2, xc_info_2, fid2, XC_POLARIZED)
!        -------------------------------------------------------------
      endif
      if (xc_f90_info_kind(xc_info_1) == XC_EXCHANGE) then
         call xc_f90_info_name(xc_info_1,FunctionalName_X)
         FunctionalID_X = xc_f90_info_number(xc_info_1)
         Exchange = .true.
      else if (xc_f90_info_kind(xc_info_1) == XC_CORRELATION) then
         call xc_f90_info_name(xc_info_1,FunctionalName_C)
         FunctionalID_C = xc_f90_info_number(xc_info_1)
         Correlation = .true.
      else
         call ErrorHandler( "initExchCorrFunctional",                 &
                            "The functional is neither exchange nor correlation", fname1 )
      endif
      if (xc_f90_info_kind(xc_info_2) == XC_EXCHANGE) then
         call xc_f90_info_name(xc_info_2,FunctionalName_X)
         FunctionalID_X = xc_f90_info_number(xc_info_2)
         Exchange = .true.
      else if (xc_f90_info_kind(xc_info_2) == XC_CORRELATION) then
         call xc_f90_info_name(xc_info_2,FunctionalName_C)
         FunctionalID_C = xc_f90_info_number(xc_info_2)
         Correlation = .true.
      else
         call ErrorHandler( "initExchCorrFunctional",                 &
                            "The functional is neither exchange nor correlation", fname2 )
      endif
      if (.not.Exchange .or. .not.Correlation) then
         call ErrorHandler( "initExchCorrFunctional",                 &
                            "The functional does not include both exchange and correlation", excorr_name )
      else
         ECin1 = .false.
      endif
   endif
!
!  ===================================================================
!  Determine the functional type via xc_info_1.
!  ===================================================================
   select case (xc_f90_info_family(xc_info_1))
      case (XC_FAMILY_LDA);
         FunctionalType = LDA
      case (XC_FAMILY_GGA);
         FunctionalType = GGA
      case (XC_FAMILY_HYB_GGA);
         FunctionalType = HybridGGA
      case (XC_FAMILY_MGGA);
         FunctionalType = MGGA
      case (XC_FAMILY_HYB_MGGA);
         FunctionalType = HybridMGGA
      case default;
         FunctionalType = UnkownXC
   end select
!
!  ===================================================================
!  The following piece of code is copied from program xcinfo in the online
!  libxc manual
!  -------------------------------------------------------------------
   if (print_level >= 0) then
      do n = 1, NumFunctionals
         if (n == 1) then
            xc_info_p => xc_info_1
            fname = fname1
         else
            xc_info_p => xc_info_2
            fname = fname2
         endif
         select case(xc_f90_info_kind(xc_info_p))
            case (XC_EXCHANGE)
               write(func_kind, '(a)') "an exchange functional"
            case (XC_CORRELATION)
               write(func_kind, '(a)') "a correlation functional"
            case (XC_EXCHANGE_CORRELATION)
               write(func_kind, '(a)') "an exchange-correlation functional"
            case (XC_KINETIC)
               write(func_kind, '(a)') "a kinetic energy functional"
            case default
               write(func_kind, '(a)') "of unknown kind"
         end select
         select case (xc_f90_info_family(xc_info_p))
            case (XC_FAMILY_LDA);
               write(func_family,'(a)') "LDA"
            case (XC_FAMILY_GGA);
               write(func_family,'(a)') "GGA"
            case (XC_FAMILY_HYB_GGA);
               write(func_family,'(a)') "Hybrid GGA"
            case (XC_FAMILY_MGGA);
               write(func_family,'(a)') "MGGA"
            case (XC_FAMILY_HYB_MGGA);
               write(func_family,'(a)') "Hybrid MGGA"
            case default;
               write(func_family,'(a)') "unknown"
         end select
         call xc_f90_info_name(xc_info_p,func_name)
         write(6,'(/,"Functional name: ",a)')fname
         write(6,'("This functional ''", a, "'' is ",a,               &
        &          ", it belongs to the ''",a,                        &
        &          "'' family and is defined in the reference(s):")') &
               trim(func_name), trim(func_kind), trim(func_family)
         i = 0
         call xc_f90_info_refs(xc_info_p, i, ref)
         do while(i >= 0)
            write(6,'(a,i1,2a)') "[", i, "] ", trim(ref)
            call xc_f90_info_refs(xc_info_p, i, ref)
         end do
      enddo
   endif
#endif
!
   nV = 0
   nE = 0
   Initialized = .true.
!
   end subroutine initExchCorrFunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endExchCorrFunctional()
!  ===================================================================
   implicit none
!
   if ( allocated(Vexc_v) ) then
      deallocate( Vexc_v )
   endif
!
   if ( allocated(Eexc_v) ) then
      deallocate( Eexc_v )
   endif
!
   nV = 0
   nE = 0
!
#ifdef LIBXC
   if (ECin1) then
!     ----------------------------------------------------------------
      call xc_f90_func_end(xc_func_1)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call xc_f90_func_end(xc_func_1)
!     ----------------------------------------------------------------
      call xc_f90_func_end(xc_func_2)
!     ----------------------------------------------------------------
   endif
#endif
!
   FunctionalName_X  = Null_String
   FunctionalName_C  = Null_String
   FunctionalName_XC = Null_String
   FunctionalID_X = -1
   FunctionalID_C = -1
   FunctionalID_XC = -1
   FunctionalType  = Null_String
!
   NumFunctionals = 0
   Exchange = .false.
   Correlation = .false.
   ECin1 = .false.
   Initialized = .false.
!
   end subroutine endExchCorrFunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrEnDen_s() result(E)
!  ===================================================================
   implicit none
!
   real (kind=RealKind) :: E
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrEnDen_s',                         &
                        'ExchCorrFunctional is not initialized')
   endif
!
   E = Eexc_s
!
   end function getExchCorrEnDen_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrEnDen_v(n_Rpts) result(pE)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_Rpts
!
   real (kind=RealKind), pointer :: pE(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrEnDen_v',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( n_Rpts < 1 .or. n_Rpts > nE ) then
      call ErrorHandler("getExchCorrEnDen_v","Num of points out of range",n_RptS)
   endif
!
   pE => Eexc_v(1:n_Rpts)
!
   end function getExchCorrEnDen_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrPot_s(is) result(V)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   real (kind=RealKind) :: V
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrPot_s',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( is < 1 .or. is > n_spin_pola ) then
      call ErrorHandler("getExchCorrPot_s","Wrong spin index",is)
   endif
!
   V = Vexc_s(is)
!
   end function getExchCorrPot_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getExchCorrPot_v(n_Rpts,is) result(pV)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is, n_Rpts
!
   real (kind=RealKind), pointer :: pV(:)
!
   if (.not.Initialized) then
      call ErrorHandler('getExchCorrPot_v',                         &
                        'ExchCorrFunctional is not initialized')
   else if ( is < 1 .or. is > n_spin_pola ) then
      call ErrorHandler("getExchCorrPot_v","Invalid spin index",is)
   else if ( n_Rpts < 1 .or. n_Rpts > nV ) then
      call ErrorHandler("getExchCorrPot_v","n_Rpts out of range",n_Rpts)
   endif
!
   pV => Vexc_v(1:n_Rpts,is)
!
   end function getExchCorrPot_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSphExchangeCorrelation_s(rho_den,der_rho_den,        &
                                          mag_den,der_mag_den)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: is , n
!
   real (kind=RealKind), intent(in) :: rho_den
   real (kind=RealKind), intent(in), optional :: der_rho_den
   real (kind=RealKind), intent(in), optional :: mag_den
   real (kind=RealKind), intent(in), optional :: der_mag_den
!
#ifdef LIBXC
   real (kind=RealKind) :: rho(2), exc(1), vxc(2)
   real (kind=RealKind) :: sigma(3), vsig(3)
   TYPE(xc_f90_pointer_t), pointer :: xc_func_p
#endif
!
   if (.not.Initialized) then
      call ErrorHandler('calSphExchangeCorrelation_s',                &
                        'ExchCorrFunctional is not initialized')
   else if (n_spin_pola == 2 .and. .not.present(mag_den)) then
      call ErrorHandler('calSphExchangeCorrelation_s',                &
                        'magnetic moment density is required')
   endif
!
   if (LegacyFunctionalID == 0 .or. LegacyFunctionalID == 1) then
      if (n_spin_pola == 1) then
!        ------------------------------------------------------------
         call calExchCorr_s( rho_den )
!        ------------------------------------------------------------
      else
         do is = 1, n_spin_pola
!           ---------------------------------------------------------
            call calExchCorr_s( rho_den, mag_den, is )
!           ---------------------------------------------------------
         enddo
      endif
#ifdef LIBXC
   else 
      if (n_spin_pola == 1) then
         rho(1) = rho_den
      else
         rho(1) = HALF*(rho_den+mag_den)
         rho(2) = HALF*(rho_den-mag_den)
      endif
      Vexc_s = ZERO; Eexc_s = ZERO
      do n = 1, NumFunctionals
         if (n == 1) then
            xc_func_p => xc_func_1
         else
            xc_func_p => xc_func_2
         endif 
         if (FunctionalType == LDA) then
!           ---------------------------------------------------------
            call xc_f90_lda_exc_vxc(xc_func_p, 1, rho(1), exc(1), vxc(1))
!           ---------------------------------------------------------
         else if (FunctionalType == GGA) then
            sigma = ZERO
            if (n_spin_pola == 1 .and. present(der_rho_den)) then
               sigma(1) = der_rho_den*der_rho_den
            else if (n_spin_pola == 2 .and.                           &
                     present(der_rho_den) .and.  present(der_mag_den)) then
               sigma(1) = (HALF*(der_rho_den+der_mag_den))**2
               sigma(2) = HALF*(der_rho_den+der_mag_den)*HALF*(der_rho_den-der_mag_den)
               sigma(3) = HALF*(der_rho_den-der_mag_den)*HALF*(der_rho_den-der_mag_den)
            else
!              ------------------------------------------------------
               call ErrorHandler('calSphExchangeCorrelation_s','Density derivative is required')
!              ------------------------------------------------------
            endif
!           ---------------------------------------------------------
            call xc_f90_gga_exc_vxc(xc_func_p, 1, rho(1), sigma(1), exc(1), vxc(1), vsig(1))
!           ---------------------------------------------------------
         else
            call ErrorHandler('calSphExchangeCorrelation_s',          &
                 'Evauation of the functionals other LDA and GGA has not been implemented yet.')
         endif
         do is = 1, n_spin_pola
            Vexc_s(is) = Vexc_s(is) + vxc(is)
         enddo
         Eexc_s = Eexc_s + exc(1)
      enddo
      Vexc_s = energy_units_conv*Vexc_s
      Eexc_s = energy_units_conv*Eexc_s
#else
   else
      call ErrorHandler( "calExchangeCorrelation_ssp",                &
                         "The LDA functional is not implemented.",LegacyFunctionalID )
#endif
   endif
!
   end subroutine calSphExchangeCorrelation_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSphExchangeCorrelation_v(n_Rpts,rho_den,der_rho_den, &
                                          mag_den,der_mag_den)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_Rpts
   integer (kind=IntKind) :: is, i, n, js
!
   real (kind=RealKind), intent(in) :: rho_den(n_Rpts)
   real (kind=RealKind), intent(in), optional :: der_rho_den(n_Rpts)
   real (kind=RealKind), intent(in), optional :: mag_den(n_Rpts)
   real (kind=RealKind), intent(in), optional :: der_mag_den(n_Rpts)
!
#ifdef LIBXC
   real (kind=RealKind), allocatable :: vxc(:), exc(:), rho(:)
   real (kind=RealKind), allocatable :: sigma(:), vsig(:)
   TYPE(xc_f90_pointer_t), pointer :: xc_func_p
#endif
!
   if (.not.Initialized) then
      call ErrorHandler('calSphExchangeCorrelation_v',                &
                        'ExchCorrFunctional is not initialized')
   else if (n_spin_pola == 2 .and. .not.present(mag_den)) then
      call ErrorHandler('calSphExchangeCorrelation_vsp',              &
                        'magnetic moment density is required')
   endif
!
   if ( allocated(Eexc_v) .and. nE < n_Rpts ) then
      deallocate( Eexc_v )
   endif
   if ( .not.allocated(Eexc_v) ) then
      allocate( Eexc_v(n_Rpts) )
      nE = n_Rpts
   endif
!
   if ( allocated(Vexc_v) .and. nV < n_Rpts ) then
      deallocate( Vexc_v )
   endif
   if ( .not.allocated(Vexc_v) ) then
      allocate( Vexc_v(n_Rpts,n_spin_pola) )
      nV = n_Rpts
   endif
!
   if (LegacyFunctionalID == 0 .or. LegacyFunctionalID == 1) then
      if (n_spin_pola == 1) then
!        ------------------------------------------------------------
         call calExchCorr_v( n_Rpts, rho_den )
!        ------------------------------------------------------------
      else
         do is = 1, n_spin_pola
!           ---------------------------------------------------------
            call calExchCorr_v( n_Rpts, rho_den, mag_den, is )
!           ---------------------------------------------------------
         enddo
      endif
#ifdef LIBXC
   else 
      allocate( rho(n_Rpts*n_spin_pola) )
      allocate( vxc(n_Rpts*n_spin_pola), exc(n_Rpts) )
      if (FunctionalType == GGA) then
         if (n_spin_pola == 1) then
            allocate(sigma(n_Rpts), vsig(n_Rpts))
         else
            allocate(sigma(3*n_Rpts), vsig(3*n_Rpts))
         endif
      endif
      if (n_spin_pola == 1) then
         rho = rho_den
      else
         do i = 1, n_Rpts
            rho(2*i-1) = HALF*(rho_den(i)+mag_den(i))
            rho(2*i)   = HALF*(rho_den(i)-mag_den(i))
         enddo
      endif
      Vexc_v = ZERO; Eexc_v = ZERO
      do n = 1, NumFunctionals
         if (n == 1) then
            xc_func_p => xc_func_1
         else
            xc_func_p => xc_func_2
         endif 
         if (FunctionalType == LDA) then
!           ---------------------------------------------------------
            call xc_f90_lda_exc_vxc(xc_func_p, n_Rpts, rho(1), exc(1), vxc(1))
!           ---------------------------------------------------------
         else if (FunctionalType == GGA) then
            if (n_spin_pola == 1 .and. present(der_rho_den)) then
               do i = 1, n_Rpts
                  sigma(i) = der_rho_den(i)*der_rho_den(i)
               enddo
            else if (n_spin_pola == 2 .and.                           &
                     present(der_rho_den) .and.  present(der_mag_den)) then
               do i = 1, n_Rpts
                  sigma(3*i-2) = HALF*(der_rho_den(i)+der_mag_den(i))*HALF*(der_rho_den(i)+der_mag_den(i))
                  sigma(3*i-1) = HALF*(der_rho_den(i)+der_mag_den(i))*HALF*(der_rho_den(i)-der_mag_den(i))
                  sigma(3*i  ) = HALF*(der_rho_den(i)-der_mag_den(i))*HALF*(der_rho_den(i)-der_mag_den(i))
               enddo
            else
!              ------------------------------------------------------
               call ErrorHandler('calSphExchangeCorrelation_v','Density derivative is required')
!              ------------------------------------------------------
            endif
!           ---------------------------------------------------------
            call xc_f90_gga_exc_vxc(xc_func_p, n_Rpts, rho(1), sigma(1), exc(1), vxc(1), vsig(1))
!           ---------------------------------------------------------
         else
            call ErrorHandler('calSphExchangeCorrelation_s',         &
                 'Evauation of the functionals other LDA and GGA has not been implemented yet.')
         endif
         do is = 1, n_spin_pola
            js = n_spin_pola - is
            do i = 1, n_Rpts
               Vexc_v(i,is) = Vexc_v(i,is) + vxc(n_spin_pola*i-js)
            enddo
         enddo
         do i = 1, n_Rpts
            Eexc_v(i) = Eexc_v(i) + exc(i)
          enddo
      enddo
      Vexc_v = energy_units_conv*Vexc_v
      Eexc_v = energy_units_conv*Eexc_v
      deallocate( rho, vxc, exc )
      if (FunctionalType == GGA) then
         deallocate(sigma,vsig)
      endif
#else
   else
      call ErrorHandler( "calSphExchangeCorrelation_v",               &
                         "The LDA functional is not implemented.",    &
                         LegacyFunctionalID )
#endif
   endif
!
   end subroutine calSphExchangeCorrelation_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchangeCorrelation_s(rho_den,grad_rho_den,          &
                                       mag_den,grad_mag_den)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: is, n 
!
   real (kind=RealKind), intent(in) :: rho_den
   real (kind=RealKind), intent(in), optional :: grad_rho_den(3)
   real (kind=RealKind), intent(in), optional :: mag_den
   real (kind=RealKind), intent(in), optional :: grad_mag_den(3)
!
#ifdef LIBXC
   real (kind=RealKind) :: rho(2), exc(1), vxc(2)
   real (kind=RealKind) :: sigma(3), vsig(3)
   TYPE(xc_f90_pointer_t), pointer :: xc_func_p
#endif
!
   if (.not.Initialized) then
      call ErrorHandler('calExchangeCorrelation_s',                  &
                        'ExchCorrFunctional is not initialized')
   else if (n_spin_pola == 2 .and. .not.present(mag_den)) then
      call ErrorHandler('calExchangeCorrelation_s',                  &
                        'magnetic moment density is required')
   endif
!
   if (LegacyFunctionalID == 0 .or. LegacyFunctionalID == 1) then
      if (n_spin_pola == 1) then
!        ------------------------------------------------------------
         call calExchCorr_s( rho_den )
!        ------------------------------------------------------------
      else
         do is = 1, n_spin_pola
!           ---------------------------------------------------------
            call calExchCorr_s( rho_den, mag_den, is )
!           ---------------------------------------------------------
         enddo
      endif
#ifdef LIBXC
   else 
      if (n_spin_pola == 1) then
         rho(1) = rho_den
      else
         rho(1) = HALF*(rho_den+mag_den)
         rho(2) = HALF*(rho_den-mag_den)
      endif
      Vexc_s = ZERO; Eexc_s = ZERO
      do n = 1, NumFunctionals
         if (n == 1) then
            xc_func_p => xc_func_1
         else
            xc_func_p => xc_func_2
         endif 
         if (FunctionalType == LDA) then
!           ---------------------------------------------------------
            call xc_f90_lda_exc_vxc(xc_func_p, 1, rho(1), exc(1), vxc(1))
!           ---------------------------------------------------------
         else if (FunctionalType == GGA) then
            if (n_spin_pola == 1 .and. present(grad_rho_den)) then
               sigma(1) = grad_rho_den(1)**2+grad_rho_den(2)**2+grad_rho_den(3)**2
            else if (n_spin_pola == 2 .and.                           &
                     present(grad_rho_den) .and.  present(grad_mag_den)) then
               sigma(1) = HALF*(grad_rho_den(1)+grad_mag_den(1))*HALF*(grad_rho_den(1)+grad_mag_den(1)) + &
                          HALF*(grad_rho_den(2)+grad_mag_den(2))*HALF*(grad_rho_den(2)+grad_mag_den(2)) + &
                          HALF*(grad_rho_den(3)+grad_mag_den(3))*HALF*(grad_rho_den(3)+grad_mag_den(3))
               sigma(2) = HALF*(grad_rho_den(1)+grad_mag_den(1))*HALF*(grad_rho_den(1)-grad_mag_den(1)) + &
                          HALF*(grad_rho_den(2)+grad_mag_den(2))*HALF*(grad_rho_den(2)-grad_mag_den(2)) + &
                          HALF*(grad_rho_den(3)+grad_mag_den(3))*HALF*(grad_rho_den(3)-grad_mag_den(3))
               sigma(3) = HALF*(grad_rho_den(1)-grad_mag_den(1))*HALF*(grad_rho_den(1)-grad_mag_den(1)) + &
                          HALF*(grad_rho_den(2)-grad_mag_den(2))*HALF*(grad_rho_den(2)-grad_mag_den(2)) + &
                          HALF*(grad_rho_den(3)-grad_mag_den(3))*HALF*(grad_rho_den(3)-grad_mag_den(3))
            else
!              ------------------------------------------------------
               call ErrorHandler('calExchangeCorrelation_s','Density gradient is required')
!              ------------------------------------------------------
            endif
!           ---------------------------------------------------------
            call xc_f90_gga_exc_vxc(xc_func_p, 1, rho(1), sigma(1), exc(1), vxc(1), vsig(1))
!           ---------------------------------------------------------
         else
            call ErrorHandler('calExchangeCorrelation_s',            &
                 'Evauation of the functionals other LDA and GGA has not been implemented yet.')
         endif
         do is = 1, n_spin_pola
            Vexc_s(is) = Vexc_s(is) + vxc(is)
         enddo
         Eexc_s = Eexc_s + exc(1)
      enddo
      Vexc_s = energy_units_conv*Vexc_s
      Eexc_s = energy_units_conv*Eexc_s
#else
   else
      call ErrorHandler( "calExchangeCorrelation_s",                 &
                         "The LDA functional is not implemented.",   &
                         LegacyFunctionalID )
#endif
   endif
!
   end subroutine calExchangeCorrelation_s
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchangeCorrelation_v( n_Rpts, rho_den, grad_rho_den,&
                                        mag_den, grad_mag_den )
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n_Rpts
   integer (kind=IntKind) :: is, i, n, js
!
   real (kind=RealKind), intent(in) :: rho_den(n_Rpts)
   real (kind=RealKind), intent(in), optional :: grad_rho_den(3,n_Rpts)
   real (kind=RealKind), intent(in), optional :: mag_den(n_Rpts)
   real (kind=RealKind), intent(in), optional :: grad_mag_den(3,n_Rpts)
!
#ifdef LIBXC
   real (kind=RealKind), allocatable :: vxc(:), exc(:), rho(:)
   real (kind=RealKind), allocatable :: sigma(:), vsig(:)
   TYPE(xc_f90_pointer_t), pointer :: xc_func_p
#endif
!
   if (.not.Initialized) then
      call ErrorHandler('calExchangeCorrelation_v',                   &
                        'ExchCorrFunctional is not initialized')
   else if (n_spin_pola == 2 .and. .not.present(mag_den)) then
      call ErrorHandler('calExchangeCorrelation_v',                   &
                        'magnetic moment density is required')
   endif
!
   if ( allocated(Eexc_v) .and. nE < n_Rpts ) then
      deallocate( Eexc_v )
   endif
   if ( .not.allocated(Eexc_v) ) then
      allocate( Eexc_v(n_Rpts) )
      nE = n_Rpts
   endif
!
   if ( allocated(Vexc_v) .and. nV < n_Rpts ) then
      deallocate( Vexc_v )
   endif
   if ( .not.allocated(Vexc_v) ) then
      allocate( Vexc_v(n_Rpts,n_spin_pola) )
      nV = n_Rpts
   endif
!
   if (LegacyFunctionalID == 0 .or. LegacyFunctionalID == 1) then
      if (n_spin_pola == 1) then
!        ------------------------------------------------------------
         call calExchCorr_v( n_Rpts, rho_den )
!        ------------------------------------------------------------
      else
         do is = 1, n_spin_pola
!           ---------------------------------------------------------
            call calExchCorr_v( n_Rpts, rho_den, mag_den, is )
!           ---------------------------------------------------------
         enddo
      endif
#ifdef LIBXC
   else 
      allocate( rho(n_Rpts*n_spin_pola) )
      allocate( vxc(n_Rpts*n_spin_pola), exc(n_Rpts) )
      if (n_spin_pola == 1) then
         rho = rho_den
      else
         do i = 1, n_Rpts
            rho(2*i-1) = HALF*(rho_den(i)+mag_den(i))
            rho(2*i)   = HALF*(rho_den(i)-mag_den(i))
         enddo
      endif
      Vexc_v = ZERO; Eexc_v = ZERO
      do n = 1, NumFunctionals
         if (n == 1) then
            xc_func_p => xc_func_1
         else
            xc_func_p => xc_func_2
         endif 
         if (FunctionalType == LDA) then
!           ---------------------------------------------------------
            call xc_f90_lda_exc_vxc(xc_func_p, n_Rpts, rho(1), exc(1), vxc(1))
!           ---------------------------------------------------------
         else if (FunctionalType == GGA) then
            if (n_spin_pola == 1 .and. present(grad_rho_den)) then
               allocate(sigma(n_Rpts), vsig(n_Rpts))
               do i = 1, n_Rpts
                  sigma(i) = grad_rho_den(1,i)**2+grad_rho_den(2,i)**2+grad_rho_den(3,i)**2
               enddo
            else if (n_spin_pola == 2 .and.                           &
                     present(grad_rho_den) .and.  present(grad_mag_den)) then
               allocate(sigma(3*n_Rpts), vsig(3*n_Rpts))
               do i = 1, n_Rpts
                  sigma(3*i-2) = HALF*(grad_rho_den(1,i)+grad_mag_den(1,i))*HALF*(grad_rho_den(1,i)+grad_mag_den(1,i)) + &
                                 HALF*(grad_rho_den(2,i)+grad_mag_den(2,i))*HALF*(grad_rho_den(2,i)+grad_mag_den(2,i)) + &
                                 HALF*(grad_rho_den(3,i)+grad_mag_den(3,i))*HALF*(grad_rho_den(3,i)+grad_mag_den(3,i))
                  sigma(3*i-1) = HALF*(grad_rho_den(1,i)+grad_mag_den(1,i))*HALF*(grad_rho_den(1,i)-grad_mag_den(1,i)) + &
                                 HALF*(grad_rho_den(2,i)+grad_mag_den(2,i))*HALF*(grad_rho_den(2,i)-grad_mag_den(2,i)) + &
                                 HALF*(grad_rho_den(3,i)+grad_mag_den(3,i))*HALF*(grad_rho_den(3,i)-grad_mag_den(3,i))
                  sigma(3*i)   = HALF*(grad_rho_den(1,i)-grad_mag_den(1,i))*HALF*(grad_rho_den(1,i)-grad_mag_den(1,i)) + &
                                 HALF*(grad_rho_den(2,i)-grad_mag_den(2,i))*HALF*(grad_rho_den(2,i)-grad_mag_den(2,i)) + &
                                 HALF*(grad_rho_den(3,i)-grad_mag_den(3,i))*HALF*(grad_rho_den(3,i)-grad_mag_den(3,i))
               enddo
            else
!              ------------------------------------------------------
               call ErrorHandler('calExchangeCorrelation_v','Density gradient is required')
!              ------------------------------------------------------
            endif
!           ---------------------------------------------------------
            call xc_f90_gga_exc_vxc(xc_func_p, n_Rpts, rho(1), sigma(1), exc(1), vxc(1), vsig(1))
!           ---------------------------------------------------------
            deallocate(sigma,vsig)
         else
!           ---------------------------------------------------------
            call ErrorHandler('calExchangeCorrelation_v',            &
                 'Evauation of the functionals other LDA and GGA has not been implemented yet.')
!           ---------------------------------------------------------
         endif
         do is = 1, n_spin_pola
            js = n_spin_pola - is
            do i = 1, n_Rpts
               Vexc_v(i,is) = Vexc_v(i,is) + vxc(n_spin_pola*i-js)
            enddo
         enddo
         do i = 1, n_Rpts
            Eexc_v(i) = Eexc_v(i) + exc(i)
         enddo
      enddo
      Vexc_v = energy_units_conv*Vexc_v
      Eexc_v = energy_units_conv*Eexc_v
      deallocate( rho, vxc, exc )
#else
   else
      call ErrorHandler( "calExchangeCorrelation_v",                  &
                         "The LDA functional is not implemented.",    &
                         LegacyFunctionalID )
#endif
   endif
!
   end subroutine calExchangeCorrelation_v
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchCorr_s( rho_den, mag_den, js )
!  ===================================================================
!
   implicit   none
!
   integer (kind=IntKind), intent(in), optional :: js
   integer (kind=IntKind) :: is
!
   real (kind=RealKind), intent(in) :: rho_den
   real (kind=RealKind), intent(in), optional :: mag_den
!
   real (kind=RealKind) :: dz
   real (kind=RealKind) :: r_s
   real (kind=RealKind) :: sp
!
!  ===================================================================
!  exchange-correlation potential and density at a point
!  ===================================================================
   if (present(js)) then
      is = js
   else
      is = 1
   endif
!
   sp = THREE-TWO*is 
   if ( rho_den > ZERO ) then
      if (present(mag_den)) then
         dz = mag_den/rho_den
      else
         dz = ZERO
      endif
      if (abs(dz) > ONE+TEN2m10) then
!        =============================================================
!        call ErrorHandler('calExchCorr','abs(dz) > 1',dz,.true.)
! modified here to help the vaccume case 05/29/17
         call WarningHandler('calExchCorr','abs(dz) > 1',dz,.true.)
         if ( dz > ONE ) then
            dz = ONE
         else
            dz = -ONE
         endif
!     else if ( dz > ONE ) then
!        dz = ONE
!     else if (dz < -ONE) then
!        dz = -ONE
!        =============================================================
      endif
      r_s = (THREE/(PI4*rho_den))**THIRD
!     ----------------------------------------------------------------
      call LDAfunctional( r_s, dz, sp, Eexc_s, Vexc_s(is) )
!     ----------------------------------------------------------------
   else
      Vexc_s(is) = ZERO
      Eexc_s = ZERO
   endif
!
   end subroutine calExchCorr_s 
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calExchCorr_v( n_Rpts, rho_den, mag_den, js )
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in), optional :: js
   integer (kind=IntKind) :: is
   integer (kind=IntKind), intent(in) :: n_Rpts
!
   real (kind=RealKind), intent(in) :: rho_den(n_Rpts)
   real (kind=RealKind), intent(in), optional :: mag_den(n_Rpts)
!
   integer (kind=IntKind) :: ir
!
   real (kind=RealKind) :: dz
   real (kind=RealKind) :: r_s
   real (kind=RealKind) :: sp
!
!  ===================================================================
!  exchange-correlation potential and density for an array of points
!  ===================================================================
   if (present(js)) then
      is = js
   else
      is = 1
   endif
!
   sp = THREE-TWO*is 
!
   do ir = 1,n_Rpts
      if ( rho_den(ir) > ZERO ) then
         if (present(mag_den)) then
            dz = mag_den(ir)/rho_den(ir)
         else
            dz = ZERO
         endif
         if (abs(dz) > ONE+TEN2m10) then
!           ==========================================================
!           call ErrorHandler('calExchCorr','abs(dz) > 1',dz,.true.)
! modified here to help the vaccume case 05/29/17
            call WarningHandler('calExchCorr','abs(dz) > 1',dz,.true.)
            if ( dz > ONE ) then
               dz = ONE
            else
               dz = -ONE
            endif
!        else if ( dz > ONE ) then
!           dz = ONE
!        else if (dz < -ONE) then
!           dz = -ONE
!           ==========================================================
         endif
         r_s = (THREE/(PI4*rho_den(ir)))**THIRD
!        -------------------------------------------------------------
         call LDAfunctional( r_s, dz, sp, Eexc_v(ir), Vexc_v(ir,is) )
!        -------------------------------------------------------------
      else
         Vexc_v(ir,is) = ZERO
         Eexc_v(ir) = ZERO
      endif
   enddo
!
   end subroutine calExchCorr_v 
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine LDAfunctional(r_s, dz, sp, exchg, vxchg)
!  ===================================================================
!
   implicit none
!
   integer (kind=IntKind) ::  i
!
   real (kind=RealKind), intent(in) ::  r_s
   real (kind=RealKind), intent(in) ::  dz
   real (kind=RealKind), intent(in) ::  sp
   real (kind=RealKind), intent(out) ::  exchg
   real (kind=RealKind), intent(out) ::  vxchg
!
   real (kind=RealKind) ::  vxx(2)
   real (kind=RealKind) ::  g(3)
   real (kind=RealKind) ::  dg(3)
   real (kind=RealKind) ::  tbq(3)
   real (kind=RealKind) ::  tbxq(3)
   real (kind=RealKind) ::  bxx(3)
   real (kind=RealKind) ::  q(3)
   real (kind=RealKind) ::  bb(3)
   real (kind=RealKind) ::  cx(3)
   real (kind=RealKind) ::  fm
   real (kind=RealKind) ::  fdz
   real (kind=RealKind) ::  ex
   real (kind=RealKind) ::  exf
   real (kind=RealKind) ::  xp
   real (kind=RealKind) ::  xf
   real (kind=RealKind) ::  gp
   real (kind=RealKind) ::  gf
   real (kind=RealKind) ::  exc
   real (kind=RealKind) ::  excf
   real (kind=RealKind) ::  dedz
   real (kind=RealKind) ::  gpp
   real (kind=RealKind) ::  gfp
   real (kind=RealKind) ::  depd
   real (kind=RealKind) ::  defd
   real (kind=RealKind) ::  decd
   real (kind=RealKind) ::  bfc
   real (kind=RealKind) ::  zp1
   real (kind=RealKind) ::  zm1
   real (kind=RealKind) ::  xr
   real (kind=RealKind) ::  pex
   real (kind=RealKind) ::  xrsq
   real (kind=RealKind) ::  qi
   real (kind=RealKind) ::  txb
   real (kind=RealKind) ::  fx
   real (kind=RealKind) ::  arct
   real (kind=RealKind) ::  dxs
   real (kind=RealKind) ::  vcc
   real (kind=RealKind) ::  facc
   real (kind=RealKind) ::  ecp
   real (kind=RealKind) ::  zp3
   real (kind=RealKind) ::  zm3
   real (kind=RealKind) ::  zp3m3
   real (kind=RealKind) ::  fx1
   real (kind=RealKind) ::  z4
   real (kind=RealKind) ::  fz
   real (kind=RealKind) ::  beta
   real (kind=RealKind) ::  ec
   real (kind=RealKind) ::  f3ex
!
!  ==================================================================
!  data for von Barth-Hedin
!  ==================================================================
   real (kind=RealKind), parameter ::  ccp = 0.0450d+00
   real (kind=RealKind), parameter ::  rp = 21.0d+00
   real (kind=RealKind), parameter ::  ccf = 0.02250d+00
   real (kind=RealKind), parameter ::  rf = 52.9166820d+00
!
!  ==================================================================
!  data for Vosko-Wilks-Nusair
!  ==================================================================
!
   real (kind=RealKind), parameter ::  a(3) = (/-0.0337740d+00,      &
                                       0.06218140d+00,0.03109070d+00/)
   real (kind=RealKind), parameter ::  b(3) = (/1.131070d+00,        &
                                       3.727440d+00,7.060420d+00/)
   real (kind=RealKind), parameter ::  c(3) = (/13.00450d+00,        &
                                       12.93520d+00,18.05780d+00/)
   real (kind=RealKind), parameter ::  x0(3) = (/-0.00475840d+00,    &
                                       -0.104980d+00,-0.32500d+00/)
   real (kind=RealKind), parameter ::  cst = 1.923661050d+00
   real (kind=RealKind), parameter ::  aip = 0.916330590d+00
   real (kind=RealKind), parameter ::  fnot = 1.709920950d+00
   real (kind=RealKind), parameter ::  bip = 0.259921050d+00
   real (kind=RealKind), parameter ::  for3 = ONE + THIRD
   real (kind=RealKind), parameter ::  thrd = THIRD
!
   if ( LegacyFunctionalID == 0 ) then  
!     ================================================================
!     von Barth-Hedin  exch-corr potential
!     j. phys. c5,1629(1972)
!     ================================================================
      fm  = 2.0d+00**(4.0d+00/3.0d+00)-2.0d+00
      fdz = ((1.0d+00+dz)**(4.0d+00/3.0d+00)                         &
           +(1.0d+00-dz)**(4.0d+00/3.0d+00)-2.0d+00)/fm
      ex  = -0.916330d+00/r_s
      exf = ex*2.0d+00**0.333333330d+00
      xp  = r_s/rp
      xf  = r_s/rf
      gp = (1.0d+00+xp**3)*log(1.0d+00+1.0d+00/xp)                   &
          -xp*xp +xp/2.0d+00 - 0.333333330d+00
      gf = (1.0d+00+xf**3)*log(1.0d+00+1.0d+00/xf)                   &
          -xf*xf +xf/2.0d+00 - 0.333333330d+00
      exc  = ex-ccp*gp
      excf = exf-ccf*gf
      dedz = (4.0d+00/3.0d+00)*(excf-exc)                            &
            *((1.0d+00+dz)**(1.0d+00/3.0d+00)                        &
            -(1.0d+00-dz)**(1.0d+00/3.0d+00))/fm
      gpp = 3.0d+00*xp*xp*log(1.0d+00+1.0d+00/xp)-1.0d+00/xp         &
           +1.50d+00-3.0d+00*xp
      gfp = 3.0d+00*xf*xf*log(1.0d+00+1.0d+00/xf)-1.0d+00/xf         &
           +1.50d+00-3.0d+00*xf
      depd = -ex/r_s-ccp/rp*gpp
      defd = -exf/r_s-ccf/rf*gfp
      decd = depd+(defd-depd)*fdz
!     ================================================================
!     exchange-correlation energy
!     ================================================================
      exchg = exc + (excf-exc)*fdz
!     ================================================================
!     exchange-correlation potential
!     ================================================================
      vxchg = exc+(excf-exc)*fdz-r_s*decd/3.0d+00                     &
              +sp*(1.0d+00-sp*dz)*dedz
!
   else if ( LegacyFunctionalID == 1 ) then
!     ================================================================
!     Vosko-Wilks-Nusair exch-corr potential
!     From G.S. Painter : Phys. Rev. B24 4264,1981
!     ================================================================
!
!     ================================================================
!     generate constant coefficients for the parameterization (v-w-n)
!     ================================================================
      do i = 1,3
         cx(i)  = x0(i)**2 + b(i)*x0(i) + c(i)
         bfc    = 4.0d+00*c(i) - b(i)**2.0d+00
         q(i)   = sqrt(bfc)
         bxx(i) = b(i)*x0(i)/cx(i)
         tbq(i) = 2.0d+00*b(i)/q(i)
         tbxq(i)= tbq(i) + 4.0d+00*x0(i)/q(i)
         bb(i)  = 4.0d+00*b(i)*( 1 - x0(i)*(b(i) + 2.0d+00*x0(i))/cx(i) )
      enddo
!
      zp1  = 1.0d+00 + dz
      zm1  = 1.0d+00 - dz
      xr   = sqrt(r_s)
      pex  = -aip/r_s
      xrsq = r_s
!
!     ===============================================================
!     generate g(i)=alpha,epsilon fct.s
!     and their derivatives dg(i)
!     1=alpha(spin stiffness)  2=ecp  3=ecf
!     ===============================================================
      do i = 1,3
         qi   = q(i)
         txb  = 2.0d+00*xr + b(i)
         fx   = xrsq + xr*b(i) + c(i)
         arct = atan2(qi,txb)
         dxs  = (xr-x0(i))**2/fx
         g(i) = a(i)*( log(xrsq/fx) + tbq(i)*arct-bxx(i)*(log(dxs)  &
                +tbxq(i)*arct) )
         dg(i) = a(i)*( 2.0d+00/xr - txb/fx                         &
                -bxx(i)*(2.0d+00/(xr-x0(i))-txb/fx)                 &
                -bb(i)/(qi**2 + txb**2) )
      enddo
!
      ecp = g(2)
      zp3 = zp1**thrd
      zm3 = zm1**thrd
      zp3m3 = zp3-zm3
!     ===============================================================
!     part of last term in vx   eq(13)
!     ===============================================================
      fx1  = .50d+00*for3*pex*zp3m3
      z4   = dz**4
      fz   = cst*(zp1**for3 + zm1**for3 - 2.0d+00)
      beta = fnot*( g(3)-g(2) )/g(1) -1.0d+00
      ec   = ecp + fz*g(1)*( 1.0d+00 + z4*beta )/fnot
      ex   = pex*( 1.0d+00 + fz*bip )
      f3ex = for3*ex
!     ===============================================================
!     echange-correlation energy
!     ===============================================================
      exchg = ec + ex
!     ===============================================================
!     exchange potential
!     ===============================================================
      vxx(1) = f3ex + fx1*zm1
      vxx(2) = f3ex - fx1*zp1
!     ===============================================================
!     correlation potential
!     ===============================================================
      vcc = ec - xr*( (1.0d+00-z4*fz)*dg(2) + z4*fz*dg(3)           &
           +(1.0d+00 - z4)*fz*dg(1)/fnot )/6.0d+00
!
      facc = 4.0d+00*g(1)*( dz**3*fz*beta            &
            +( 1.0d+00 + beta*z4 )*zp3m3/(6.0d+00*bip) )/fnot
!
!     ===============================================================
!     exch-corr. potential for each spin as called in newpot
!     ===============================================================
!
      if ( sp >= 0 ) then 
         vxchg = vcc + zm1*facc + vxx(1)
      else
         vxchg = vcc - zp1*facc + vxx(2)
      endif
   else
      call ErrorHandler('LDAfunctional','Functional is not implemented yet', &
                        LegacyFunctionalID)
   endif
!
   end subroutine LDAfunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLDAFunctional() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if (FunctionalType == LDA) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isLDAFunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isGGAFunctional() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if (FunctionalType == GGA .or. FunctionalType == MGGA .or.         &
       FunctionalType == HybridGGA .or. FunctionalType == HybridMGGA) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isGGAFunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isMGGAFunctional() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if (FunctionalType == MGGA) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isMGGAFunctional
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isHybridFunctional() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if (FunctionalType == HybridGGA .or. FunctionalType == HybridMGGA) then
      y = .true.
   else
      y = .false.
   endif
!
   end function isHybridFunctional
!  ===================================================================
end module ExchCorrFunctionalModule
