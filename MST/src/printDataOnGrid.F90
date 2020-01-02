!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printDataOnGrid( grid_name, value_name, value_type, denOnGrid, lmax )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
!
   use PublicTypeDefinitionsModule, only : UniformGridStruct
!
   use Uniform3DGridModule, only : getUniform3DGrid
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
   implicit none
!
!  ==================================================================
!  Variables
!  ==================================================================
!
   character(len=*), intent(in) :: grid_name
   character(len=*), intent(in) :: value_name
   character(len=*), intent(in) :: value_type
!
   type (UniformGridStruct), pointer :: gp
!
   integer (kind=IntKind), intent(in) :: lmax
!
   real (kind=RealKind), intent(in) :: denOnGrid(:)
!
   logical :: isCharge
   logical :: isPotential
   logical :: isWave
!
   character (len=90) :: fname
!
   integer (kind=IntKind) :: i, j, k, ig
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
   gp => getUniform3DGrid(grid_name)
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
      fname = getDensityPrintFile(lmax,trim(value_type))
   else if ( isPotential ) then
      printFormat = getPotentialPrintFormat()
      fname = getPotentialPrintFile(lmax,trim(value_type))
   else if ( isWave ) then
      printFormat = getWavePrintFormat()
      fname = getWavePrintFile(lmax,trim(value_type))
   endif
!
   funit = 11
   open(unit=funit, file=trim(fname), status='unknown')
   if ( printFormat == 0 ) then
!
      write(funit, *) "# xyz sqrt(r^2) format"
      ig = 0
      do k = gp%gstart(3), gp%gend(3)
         do j = gp%gstart(2), gp%gend(2)
            do i = gp%gstart(1), gp%gend(1)
               ig = ig + 1
               rg = gp%grid_step_a*(i-1) + gp%grid_step_b*(j-1) + gp%grid_step_c*(k-1) + &
                    gp%vec_origin
               sqrt_rg = sqrt(rg(1)**2+rg(2)**2+rg(3)**2)
               WRITE(funit,'(5f15.8,2i5)') rg, sqrt_rg, &
                    denOnGrid(ig), gp%AtomOnGrid%NumLocalAtoms, gp%AtomOnGrid%NumGridPointsOnCellBound
            enddo
         enddo
      enddo
!
   else if (printFormat == 1) then
!
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
      do ig = 1, gp%NumLocalGridPoints
         write (funit, '(f15.8)') denOnGrid(ig)
      enddo
!
   endif
   call FlushFile(funit)
   close(funit)
!
   end subroutine printDataOnGrid
