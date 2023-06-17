!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDataOnGrid( gp, value_name, value_type, denOnGrid, iprint )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : PI2
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use PublicTypeDefinitionsModule, only : UniformGridStruct
!
   use PhysParamModule, only : Bohr2Angstrom
!
   use Uniform3DGridModule, only : getUniform3DGrid, getGridPosition
!
   use OutputModule, only : getDensityPrintFlag,     &
                            getDensityPrintFile,     &
                            getDensityPrintFormat,   &
                            getPotentialPrintFlag,   &
                            getPotentialPrintFile,   &
                            getPotentialPrintFormat, &
                            getWavePrintFlag,        &
                            getWavePrintFile,        &
                            getWavePrintFormat
!
   use SystemModule, only : getNumAtoms, getBravaisLattice, getAtomName
   use SystemModule, only : getAtomPosition, getNumVacancies
!
   implicit none
!
!  ==================================================================
!  Variables
!  ==================================================================
!
!  character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: value_name
   character(len=*), intent(in) :: value_type
!
   type (UniformGridStruct), intent(in) :: gp
!
   integer (kind=IntKind), intent(in) :: iprint
!
   real (kind=RealKind), intent(in), target :: denOnGrid(:)
   real (kind=RealKind), pointer :: p_denOnGrid(:,:,:)
   real (kind=RealKind) :: bravais(3,3)
!
   logical :: isCharge
   logical :: isPotential
   logical :: isWave
!
   character (len=90) :: fname
!
   integer (kind=IntKind) :: i, j, k, ig, na
   integer (kind=IntKind) :: printFormat, funit
!
   real (kind=RealKind) :: rg(3), sqrt_rg
!
   interface
      function nocaseCompare(s1,s2) result(t)
         character (len=*), intent(in) :: s1
         character (len=*), intent(in) :: s2
         logical :: t
      end function nocaseCompare
   end interface
!
   interface
      function aliasArray3_r(a,n1,n2,n3) result(p)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer(kind=IntKind), intent(in) :: n1,n2,n3
         real(kind=RealKind), target :: a(n1,n2,n3)
         real(kind=RealKind), pointer :: p(:,:,:)
      end function aliasArray3_r
   end interface
!
!  gp => getUniform3DGrid(grid_name)
   isCharge    = .false.
   isPotential = .false.
   isWave      = .false.
!
!  ==================================================================
!  Body
!  ==================================================================
   if ( nocaseCompare(value_name,'Charge') ) then
      isCharge = .true.
   else if ( nocaseCompare(value_name,'Potential') ) then
      isPotential=.true.
   else if ( nocaseCompare(value_name,'WaveFunct') ) then
      isWave=.true.
   else
      call ErrorHandler("printDensityOnGrid", "Unknown value name",value_name)
   endif
!
   if ( isCharge ) then
      printFormat = getDensityPrintFormat()
      fname = getDensityPrintFile(iprint,trim(value_type))
   else if ( isPotential ) then
      printFormat = getPotentialPrintFormat()
      fname = getPotentialPrintFile(iprint,trim(value_type))
   else if ( isWave ) then
      printFormat = getWavePrintFormat()
      fname = getWavePrintFile(iprint,trim(value_type))
   endif
!
   funit = 11
   open(unit=funit, file=trim(fname), status='unknown')
   if ( printFormat == 0 ) then
      write(funit, *) "# FORMAT: x    y    z    value"
!!    ig = 0
!!    do k = gp%gstart(3), gp%gend(3)
!!       do j = gp%gstart(2), gp%gend(2)
!!          do i = gp%gstart(1), gp%gend(1)
!!             ig = ig + 1
!!             rg = gp%grid_step_a*(i-1) + gp%grid_step_b*(j-1) + gp%grid_step_c*(k-1) + &
!!                  gp%vec_origin
!!             sqrt_rg = sqrt(rg(1)**2+rg(2)**2+rg(3)**2)
!              WRITE(funit,'(5f15.8,2i5)') rg, sqrt_rg, &
!                   denOnGrid(ig), gp%AtomOnGrid%NumLocalAtoms, gp%AtomOnGrid%NumGridPointsOnCellBound
!!             WRITE(funit,'(4f15.8)') rg, denOnGrid(ig)
!!          enddo
!!       enddo
!!    enddo
      do ig = 1, gp%ng
         rg = getGridPosition(gp,ig)
         WRITE(funit,'(4f15.8)') rg, denOnGrid(ig)
      enddo
!
   else if (printFormat == 1) then
      write(funit, '(''# vtk DataFile Version 1.0'')')
      write(funit, '(''value on the grid example'')')
      write(funit, '(''ASCII'')')
      write(funit, *)
      write(funit, '(''DATASET STRUCTURED_GRID'')')
      write(funit, '(''DIMENSIONS '' , 3i5)') gp%nga, gp%ngb, gp%ngc
      write(funit, '(''POINTS '', i10, '' double'')') gp%ng
!
      do k = gp%gstart(3), gp%gend(3)
         do j = gp%gstart(2), gp%gend(2)
            do i = gp%gstart(1), gp%gend(1)
               rg = gp%grid_step_a*(i-1) + gp%grid_step_b*(j-1) + gp%grid_step_c*(k-1) + &
                    gp%vec_origin
               WRITE(funit,'(3f10.6)') rg
            enddo
         enddo
      enddo
!
      write(funit, '(a, i10)') "POINT_DATA ", gp%ng
      write(funit, '(a)')      "SCALARS value double 1"
      write(funit, '(a)')      "LOOKUP_TABLE default"
!
!!    do ig = 1, gp%NumLocalGridPoints
      do ig = 1, gp%ng
         write (funit, '(f15.8)') denOnGrid(ig)
      enddo
!
   else if (printFormat == 2) then ! writing the data is xsf format
!
      p_denOnGrid =>  aliasArray3_r(denOnGrid,gp%nga,gp%ngb,gp%ngc)
      bravais = getBravaisLattice() 
      write(funit, '(a)')'# This XSF file is produced by MuST'
      write(funit, '(a)')'# The quantities in this file are in the unit of Angstroms'
      write(funit, '(a)')'CRYSTAL'
      write(funit, '(a)')'PRIMVEC'
      write(funit, '(3(2x,f15.11))')bravais(1:3,1)*Bohr2Angstrom
      write(funit, '(3(2x,f15.11))')bravais(1:3,2)*Bohr2Angstrom
      write(funit, '(3(2x,f15.11))')bravais(1:3,3)*Bohr2Angstrom
      write(funit, '(a)')'PRIMCOORD'
      na = getNumAtoms(); i = 1
      if (getNumVacancies() > 0) then
         write(funit, '(a)')    '# Including vacancies:'
         write(funit, '(a,2i5)')'#',na, i
         write(funit, '(a)')    '# Excluding vacancies:'
      endif
      write(funit, '(2i5)')na-getNumVacancies(), i
      do i = 1, na
         if (getAtomName(i) /= 'Va') then
            write(funit,'(2x,a2,3(2x,f15.11))')getAtomName(i),getAtomPosition(i)*Bohr2Angstrom
         endif
      enddo
      if (getNumVacancies() > 0) then
         do i = 1, na
            if (getAtomName(i) == 'Va') then
               write(funit,'(a,3(2x,f15.11))')'# Va',getAtomPosition(i)*Bohr2Angstrom
            endif
         enddo
      endif
      write(funit,'(a)')'BEGIN_BLOCK_DATAGRID_3D'
      write(funit,'(2x,a)')value_name//'_density'
      write(funit,'(2x,a)')'BEGIN_DATAGRID_3D'//'_'//value_name//'_density'
      write(funit,'(2x,3i5)')gp%nga,gp%ngb,gp%ngc
      write(funit,'(2x,3f10.5)')gp%vec_origin(1:3)*Bohr2Angstrom
      write(funit,'(2x,3f10.5)')gp%cell(1:3,1)*Bohr2Angstrom
      write(funit,'(2x,3f10.5)')gp%cell(1:3,2)*Bohr2Angstrom
      write(funit,'(2x,3f10.5)')gp%cell(1:3,3)*Bohr2Angstrom
!     write(funit,'(2x,3f10.5)')gp%grid_step_a(1:3)*Bohr2Angstrom
!     write(funit,'(2x,3f10.5)')gp%grid_step_b(1:3)*Bohr2Angstrom
!     write(funit,'(2x,3f10.5)')gp%grid_step_c(1:3)*Bohr2Angstrom
      do k = 1, gp%ngc
         do j = 1, gp%ngb
            write(funit,'(2x,5F12.5)')(p_denOnGrid(i,j,k)/Bohr2Angstrom**3,i=1,gp%nga)
         enddo
!        if (k < gp%ngc) then
!           write(funit,'(a)')' '
!        endif
      enddo
      write(funit,'(2x,a)')'END_DATAGRID_3D'
      write(funit,'(a)')'END_BLOCK_DATAGRID_3D'
   endif
   call FlushFile(funit)
   close(funit)
!
   end subroutine printDataOnGrid
