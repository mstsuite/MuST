!  ********************************************************************
!  * MODULE NAME    : OutputModule                                    *
!  *                                                                  *
!  * VERSION NUMBER : 1.0                                             *
!  * LAST MODIFIED  : MARCH 09, 2004                                  *
!  *                                                                  *
!  * DESCRIPTION: MODULE for determining standard output level,       *
!  *              standard output PE, and standard output file unit.  *
!  *                                                                  *
!  * EXTERNAL MODULE DEPENDENCE                                       *
!  * ================================================================ *
!  *    KindParamModule                                               *
!  *    ErrorHandlerModule                                            *
!  *    InputModule                                                   *
!  *                                                                  *
!  * PUBLIC FUNCTIONS                                                 *
!  * ================================================================ *
!  *    subroutine initOutput(tbl_id)                                 *
!  *    Purpose: initialize the module for output determination,      *
!  *             and open output file if necessary                    *
!  *    Input:   tbl_id  = integer, input data table ID               *
!  *    Output:  none                                                 *
!  *    ============================================================= *
!  *    subroutine endOutput()                                        *
!  *    Purpose: clean the memory allocated within the module,        *
!  *             and close output file if necessary                   *
!  *    Input:   none                                                 *
!  *    Output:  none                                                 *
!  *    ============================================================= *
!  *    function getStandardOutputLevel() result(p_level)             *
!  *    Purpose: returns the print information level                  *
!  *    Input:   none                                                 *
!  *    Output:  p_level = integer, print information level ( >= -1 ) *
!  *                                                                  *
!  ********************************************************************
module OutputModule
   use KindParamModule, only : IntKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use InputModule, only : getKeyValue
!
public :: initOutput,               &
          endOutput,                &
          isStandardOutToScreen,    &
          isOutputAtomBased,        &
          getStandardOutputLevel,   &
          getDensityPrintFlag,      &
          getDensityPrintFormat,    &
          getDensityPrintFile,      &
          getPotentialPrintFlag,    &
          getPotentialPrintFormat,  &
          getPotentialPrintFile,    &
          getWavePrintFlag,         &
          getWavePrintFormat,       &
          getWavePrintFile

!
   interface getStandardOutputLevel
      module procedure getStandardOutputLevel_0, getStandardOutputLevel_1
   end interface
!
   character (len=70), public :: cpath
!
private
!
   logical :: Initialized = .false.
   logical :: OutToScreen = .true.
   logical :: AtomBasedPrintDefined = .false.
!
   character (len=80) :: densityPrintFile
   character (len=80) :: potentialPrintFile
   character (len=80) :: wavePrintFile
!
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: CPU_based_print_level
   integer (kind=IntKind) :: densityPrintFlag, densityPrintFormat
   integer (kind=IntKind) :: potentialPrintFlag, potentialPrintFormat
   integer (kind=IntKind) :: wavePrintFlag, wavePrintFormat
   integer (kind=IntKind), allocatable :: ATOM_based_print_level(:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initOutput(tbl_id,nlvl)
!  ===================================================================
   use StringModule, only : initString, endString, getNumTokens
!
   use MPPModule, only : MyPE, NumPEs
!
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms
!
   use SystemModule, only : getNumAtoms
!
   implicit none
!
   character (len=1)  :: output_to_screen
   character (len=6)  :: cnode
   character (len=40)  :: systemid
!   character (len=70)  :: cpath
   character (len=70)  :: str_pe, str_atom, str_lvl
   character (len=120) :: opfile
!
   integer (kind=IntKind), intent(in) :: tbl_id
   integer (kind=IntKind), intent(in), optional :: nlvl
   integer (kind=IntKind) :: num_o_pe, num_o_atom, num_o_lvl, i, j, id
   integer (kind=IntKind) :: rstatus, iflag
   integer (kind=IntKind), allocatable :: o_pe(:), o_atom(:), p_lvl(:)
   integer (kind=IntKind), parameter :: offset = 10**(len(cnode)-1)
   integer (kind=IntKind), parameter :: undefined_lvl = -1010
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Output Level (>= -1)',str_lvl)
   call initString(str_lvl)
   num_o_lvl = getNumTokens()
   call endString()
!  -------------------------------------------------------------------
   if (num_o_lvl <= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler('initOutput','No data entry for Output Level')
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  The following code is added to enable the atom based print_level
!  -------------------------------------------------------------------
   LocalNumAtoms = getLocalNumAtoms()
!  -------------------------------------------------------------------
   allocate(ATOM_based_print_level(LocalNumAtoms))
   ATOM_based_print_level(1:LocalNumAtoms) = -1
   densityPrintFlag = -1
   densityPrintFormat = -1
   potentialPrintFlag = -1
   potentialPrintFormat = -1
   wavePrintFlag = -1
   wavePrintFormat = -1
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Output Atom ID (>= 0)',str_atom)
!  -------------------------------------------------------------------
!
   num_o_atom = 0
   if (len_trim(str_atom) > 0) then
!     ----------------------------------------------------------------
      call initString(str_atom)
!     num_o_atom = min(getNumTokens(), NumPEs)     ! debugged
      num_o_atom = min(getNumTokens(), getNumAtoms())
      call endString()
!     ----------------------------------------------------------------
      allocate(o_atom(1:num_o_atom), p_lvl(1:num_o_lvl))
      read(str_atom,*)o_atom(1:num_o_atom)
      read(str_lvl,*)p_lvl(1:num_o_lvl)
      LOOP_i: do i=1,num_o_atom
         if (o_atom(i) < 0) then
            if (i <= num_o_lvl) then
               j = i
            else
               j = num_o_lvl
            endif
            ATOM_based_print_level(1) = p_lvl(j)
            ATOM_based_print_level(2:LocalNumAtoms) = ATOM_based_print_level(1)
            exit LOOP_i
         else if (o_atom(i) == 0) then
            ATOM_based_print_level = undefined_lvl
         else  ! if (o_atom(i) > 0) then
            if (i <= num_o_lvl) then
               j = i
            else
               j = num_o_lvl
            endif
            LOOP_id: do id = 1, LocalNumAtoms
               if (o_atom(i) ==getGlobalIndex(id)) then
                  ATOM_based_print_level(id) = p_lvl(j)
                  exit LOOP_id
               endif
            enddo LOOP_id
         endif
      enddo LOOP_i
      deallocate(o_atom, p_lvl)
      if ( maxval(ATOM_based_print_level)>=0 ) then
         AtomBasedPrintDefined = .true.
      else
         AtomBasedPrintDefined = .false.
      endif
   else
      AtomBasedPrintDefined = .false.
   endif
!
!  ===================================================================
!  Determine the CPU_based print level
!  ===================================================================
   CPU_based_print_level = -1
!
   if ( getKeyValue(tbl_id, 'Output Electron Density ID (>= -1)',iflag) == 0 ) then
      densityPrintFlag = iflag
   endif
   if ( getKeyValue(tbl_id, 'Output Density Format', iflag) == 0 ) then
      densityPrintFormat = iflag
   endif
   if ( getKeyValue(tbl_id, 'Output Potential ID (>= -1)', iflag) == 0 ) then
      potentialPrintFlag = iflag
   endif
   if ( getKeyValue(tbl_id, 'Output Wave Function ID (>= -1)', iflag) == 0 ) then
      wavePrintFlag = iflag
   endif
!
!  temporary
!
!  rstatus = getKeyValue(tbl_id, 'Output Potential Format', potentialPrintFormat)
   potentialPrintFormat = densityPrintFormat
   wavePrintFormat = densityPrintFormat
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Output Proc. ID (>= -1)',str_pe)
!  -------------------------------------------------------------------
!
   if (len_trim(str_pe) > 0) then
!     ----------------------------------------------------------------
      call initString(str_pe)
      num_o_pe = min(getNumTokens(), NumPEs)
      call endString()
!     ----------------------------------------------------------------
!
      if (num_o_pe <= 0) then
!        -------------------------------------------------------------
         call ErrorHandler('initOutput','No data entry for Output Proc. ID')
!        -------------------------------------------------------------
      endif
!
      allocate(o_pe(1:num_o_pe), p_lvl(1:num_o_pe))
      read(str_pe,*)o_pe(1:num_o_pe)
      read(str_lvl,*)p_lvl(1:num_o_pe)
      do i=1,num_o_pe
         if (i <= num_o_lvl) then
            j = i
         else
            j = num_o_lvl
         endif
         if (o_pe(i) == MyPE .or. o_pe(i) == -1) then
            CPU_based_print_level = p_lvl(j)
            exit
         endif
      enddo
      deallocate(o_pe, p_lvl)
      if (.not.AtomBasedPrintDefined .and. CPU_based_print_level>=0 .and. &
          minval(ATOM_based_print_level) == undefined_lvl ) then
         ATOM_based_print_level = CPU_based_print_level
!        AtomBasedPrintDefined = .true.
      endif
!   else if (AtomBasedPrintDefined) then
!      do id = 1, LocalNumAtoms
!         CPU_based_print_level = max(CPU_based_print_level,           &
!                                     ATOM_based_print_level(id))
!      enddo
   else if (MyPE == 0) then
      CPU_based_print_level = 0
   endif
!
   if ( CPU_based_print_level < 0 ) then
      if ( maxval(ATOM_based_print_level) > -1 ) then
         CPU_based_print_level = 0
      endif
   else if (minval(ATOM_based_print_level) > undefined_lvl .and.      &
            maxval(ATOM_based_print_level) < 0) then
      CPU_based_print_level = -1
   endif
!
   if (present(nlvl)) then
      CPU_based_print_level = min(CPU_based_print_level,nlvl)
   endif
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Output to Screen (y/n)',output_to_screen)
!  -------------------------------------------------------------------
   if (output_to_screen == 'n' .or. output_to_screen == 'N') then
      OutToScreen = .false.
   else if (output_to_screen == 'y' .or. output_to_screen == 'Y') then
      OutToScreen = .true.
   else
!     ----------------------------------------------------------------
      call ErrorHandler('initOutput','Unknown output_to_screen value', &
                        output_to_screen)
!     ----------------------------------------------------------------
   endif
!
   write(cnode,'(i6)')MyPE+offset
   if (MyPE >= offset) then
!     ----------------------------------------------------------------
      call WarningHandler('initOutput','Processor logical index is too long')
!     ----------------------------------------------------------------
   endif
   cnode(:1)='n'
!
!  -------------------------------------------------------------------
   rstatus = getKeyValue(tbl_id,'Text Identification',systemid)
   rstatus = getKeyValue(tbl_id,'Current File Path',cpath)
!  -------------------------------------------------------------------
!
   if (.not.OutToScreen .and. CPU_based_print_level >= 0) then
!     ================================================================
!     construct filename for standard output..........................
!     default to screen ::  if(output_to_screen.eq.'y')...............
!     otherwise to file=opfile........................................
!     ================================================================
      opfile=trim(cpath)//'o_'//trim(cnode)//'_'//trim(systemid)
!     ----------------------------------------------------------------
      open(unit=6,file=opfile,form='formatted',status='unknown')
!     ----------------------------------------------------------------
   endif
!
   if (densityPrintFlag >= 0) then
!     ================================================================
!     construct filename for electronic density output
!     ----------------------------------------------------------------
      densityPrintFile='Charge_'//trim(cnode)
      if (densityPrintFormat == 0) then
         densityPrintFile=trim(densityPrintFile)//'.dat'
      else if ( densityPrintFormat == 1 ) then
         densityPrintFile=trim(densityPrintFile)//'.vtk'
      endif
   endif
!
   if (potentialPrintFlag >= 0) then
!     ================================================================
!     construct filename for electronic density output
!     ----------------------------------------------------------------
      potentialPrintFile='Pot_'//trim(cnode)
      if (potentialPrintFormat == 0) then
         potentialPrintFile=trim(potentialPrintFile)//'.dat'
      else if ( potentialPrintFormat == 1 ) then
         potentialPrintFile=trim(potentialPrintFile)//'.vtk'
      endif
   endif
!
   if (wavePrintFlag >= 0) then
!     ================================================================
!     construct filename for electronic density output
!     ----------------------------------------------------------------
      wavePrintFile='WaveFunction_'//trim(cnode)
      if (wavePrintFormat == 0) then
         wavePrintFile=trim(wavePrintFile)//'.dat'
      else if ( wavePrintFormat == 1 ) then
         wavePrintFile=trim(wavePrintFile)//'.vtk'
      endif
   endif
!
   Initialized = .true.
!
   end subroutine initOutput
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endOutput()
!  ===================================================================
   implicit none
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
!
   if (.not.OutToScreen .and. CPU_based_print_level >= 0) then
      close(6)
   endif
!
   if (allocated(ATOM_based_print_level)) then
      deallocate(ATOM_based_print_level)
   endif
!
   Initialized = .false.
!
   end subroutine endOutput
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isStandardOutToScreen() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('isStandardOutToScreen','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
!
   y = OutToScreen
!
   end function isStandardOutToScreen
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isOutputAtomBased() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('isOutputAtomBased','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
!
   y = AtomBasedPrintDefined
!
   end function isOutputAtomBased
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStandardOutputLevel_0() result(p_level)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: p_level
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getOutputLevel','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
!
   p_level = CPU_based_print_level
!
   end function getStandardOutputLevel_0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getStandardOutputLevel_1(local_index) result(p_level)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: local_index
   integer (kind=IntKind) :: p_level
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getOutputLevel','OutputModule is not initialized')
!     ----------------------------------------------------------------
   else if (local_index < 1 .or. local_index > LocalNumAtoms) then
!     ----------------------------------------------------------------
      call ErrorHandler('getOutputLevel','Invalid local atom index', &
                        local_index)
!     ----------------------------------------------------------------
   endif
!
   if ( AtomBasedPrintDefined ) then
      p_level = ATOM_based_print_level(local_index)
   else
      p_level = CPU_based_print_level
   endif
!
   end function getStandardOutputLevel_1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDensityPrintFlag() result(denFlag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: denFlag
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   denFlag = densityPrintFlag
!
   end function getDensityPrintFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDensityPrintFormat() result(denFormat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: denFormat
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   denFormat = densityPrintFormat
   end function getDensityPrintFormat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getDensityPrintFile(lmax, rho_type) result(denFile)
!  ===================================================================
   implicit none
!
   character(len=*), optional :: rho_type
   integer (kind=IntKind), intent(in) :: lmax
!
   character (len=120) :: denFile
   character (len=3) :: Lchar
   integer (kind=IntKind) :: offset
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   if ( lmax<0 ) then
      if (present(rho_type)) then
         denFile = trim(cpath)//rho_type(1:5)//trim(densityPrintFile)
      else
         denFile = trim(cpath)//trim(densityPrintFile)
      endif
   else
      offset = 100+lmax
      write(Lchar,'(i3)') offset
      Lchar(1:1) = "L"
      if ( present(rho_type) ) then
         denFile = trim(cpath)//rho_type(1:5)//Lchar//"_"//trim(densityPrintFile)
      else
         denFile = trim(cpath)//Lchar//"_"//trim(densityPrintFile)
      endif
   endif
!
   end function getDensityPrintFile
!  ===================================================================
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotentialPrintFlag() result(potFlag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: potFlag
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   potFlag = potentialPrintFlag
!
   end function getPotentialPrintFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotentialPrintFormat() result(potFormat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: potFormat
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   potFormat = potentialPrintFormat
   end function getPotentialPrintFormat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getPotentialPrintFile(lmax,pot_type) result(potFile)
!  ===================================================================
   implicit none
!
   character(len=*), optional :: pot_type
   integer (kind=IntKind), intent(in) :: lmax
!
   character (len=120) :: potFile
   character (len=3) :: Lchar
   integer (kind=IntKind) :: offset
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   if ( lmax<0 ) then
      if (present(pot_type)) then
         potFile = trim(cpath)//pot_type(1:5)//trim(potentialPrintFile)
      else
         potFile = trim(cpath)//trim(potentialPrintFile)
      endif
   else
      offset = 100+lmax
      write(Lchar,'(i3)') offset
      Lchar(1:1) = "L"
      if (present(pot_type)) then
         potFile = trim(cpath)//pot_type(1:5)//Lchar//"_"//trim(potentialPrintFile)
      else
         potFile = trim(cpath)//Lchar//"_"//trim(potentialPrintFile)
      endif
   endif
!
   end function getPotentialPrintFile
!  ===================================================================
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWavePrintFlag() result(waveFlag)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: waveFlag
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   waveFlag = wavePrintFlag
!
   end function getWavePrintFlag
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWavePrintFormat() result(waveFormat)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: waveFormat
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   waveFormat = wavePrintFormat
   end function getWavePrintFormat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getWavePrintFile(lmax,wave_type) result(waveFile)
!  ===================================================================
   implicit none
!
   character(len=*), optional :: wave_type
   integer (kind=IntKind), intent(in) :: lmax
!
   character (len=120) :: waveFile
   character (len=3) :: Lchar
   integer (kind=IntKind) :: offset
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('endOutput','OutputModule is not initialized')
!     ----------------------------------------------------------------
   endif
   if ( lmax<0 ) then
      if (present(wave_type)) then
         waveFile = trim(cpath)//wave_type(1:3)//trim(wavePrintFile)
      else
         waveFile = trim(cpath)//trim(wavePrintFile)
      endif
   else
      offset = 100+lmax
      write(Lchar,'(i3)') offset
      Lchar(1:1) = "L"
      if (present(wave_type)) then
         waveFile = trim(cpath)//wave_type(1:3)//Lchar//"_"//trim(wavePrintFile)
      else
         waveFile = trim(cpath)//Lchar//"_"//trim(wavePrintFile)
      endif
   endif
!
   end function getWavePrintFile
!  ===================================================================
end module OutputModule
