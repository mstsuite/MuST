!  ==================================================================
   subroutine driverSystemMovie(initMovie,TotalEnergy,fu,EnergyOffset)
!  ==================================================================
   use KindParamModule, only: IntKind, RealKind
!
   use PotentialTypeModule, only : printPotentialType
   use ScfDataModule, only : printScfData
   use SystemModule, only : printSystem
   use ScfDataModule, only : n_spin_cant, n_spin_pola
!
   implicit none
!
   logical, intent(in) ::initMovie
!
   real (kind=RealKind), intent(in) :: TotalEnergy
!
   integer (kind=IntKind), intent(out)  :: EnergyOffset
   integer (kind=IntKind), intent(out)  :: fu
!
   logical :: file_exist
!
   character (len=80) :: text, text0
!
   integer (kind=IntKind) :: i
!
   character (len=160) :: fname
!
   fu = 220
   if ( initMovie ) then
      fname = 'System_movie.dat'
      inquire(file=fname,exist=file_exist)
      if (file_exist) then
         open(unit=fu,file=fname,form='formatted',status='old',position='rewind')
         do i=1,31
            read(fu,*) text
         enddo
         if (text(1:14)=="#Energy Offset") then
            read(text,'(a27,i8)') text0,EnergyOffset
            write(6,*) "driverSystemMovie:: EnergyOffset",EnergyOffset
         else
            EnergyOffset = int(TotalEnergy,kind=IntKind)
         endif
         close(fu)
         open(unit=fu,file=fname,form='formatted',status='old',position='append')
      else
         EnergyOffset = int(TotalEnergy,kind=IntKind)
         open(unit=fu,file=fname,form='formatted',status='unknown')
         call print_version(fu)
         call printPotentialType(fu)
         call printScfData(fu)
         call printSystem(fu)
         write(fu,'(a,i8)') "# Energy Offset           :",EnergyOffset
         if ( n_spin_cant==2 ) then
            write(fu,'(a,i5)') "# Column Legend : ", 31
            write(fu,'(a)') "#         1. Type"
            write(fu,'(a)') "#         2. Index"
            write(fu,'(a)') "#         3. Position-X"
            write(fu,'(a)') "#         4. Position-Y"
            write(fu,'(a)') "#         5. Position-Z"
            write(fu,'(a)') "#         6. Force-X"
            write(fu,'(a)') "#         7. Force-Y"
            write(fu,'(a)') "#         8. Force-Z"
            write(fu,'(a)') "#         9. EvecIn-X"
            write(fu,'(a)') "#        10. EvecIn-Y"
            write(fu,'(a)') "#        11. EvecIn-Z"
            write(fu,'(a)') "#        12. EvecOut-X"
            write(fu,'(a)') "#        13. EvecOut-Y"
            write(fu,'(a)') "#        14. EvecOut-Z"
            write(fu,'(a)') "#        15. Bcon-X"
            write(fu,'(a)') "#        16. Bcon-Y"
            write(fu,'(a)') "#        17. Bcon-Z"
            write(fu,'(a)') "#        18. Volume - MT"
            write(fu,'(a)') "#        19. Volume - VP"
            write(fu,'(a)') "#        20. LIZ Size"
            write(fu,'(a)') "#        21. Energy"
            write(fu,'(a)') "#        22. Exc_energy"
            write(fu,'(a)') "#        23. Pressure"
            write(fu,'(a)') "#        24. Charge-MT"
            write(fu,'(a)') "#        25. Charge-VP"
            write(fu,'(a)') "#        26. Moment-MT"
            write(fu,'(a)') "#        27. Moment-VP"
            write(fu,'(a)') "#        28. RMS-Charge"
            write(fu,'(a)') "#        29. RMS-Potential"
            write(fu,'(a)') "#        30. RMS-Evec"
            write(fu,'(a)') "#        31. RMS-Bcon"
            write(fu,'(a,$)') '#'
            write(fu,'(a)')repeat('-',335)
            write(fu,'(a1,a,$)') '#'," Type "
            write(fu,'(a,$)') " Indx "
            write(fu,'(a,$)') "     X     "
            write(fu,'(a,$)') "     Y     "
            write(fu,'(a,$)') "     Z     "
            write(fu,'(a,$)') "  Force-X  "
            write(fu,'(a,$)') "  Force-Y  "
            write(fu,'(a,$)') "  Force-Z  "
            write(fu,'(a,$)') "  EvecIn-X "
            write(fu,'(a,$)') "  EvecIn-Y "
            write(fu,'(a,$)') "  EvecIn-Z "
            write(fu,'(a,$)') " EvecOut-X "
            write(fu,'(a,$)') " EvecOut-Y "
            write(fu,'(a,$)') " EvecOut-Z "
            write(fu,'(a,$)') "   Bcon-X  "
            write(fu,'(a,$)') "   Bcon-Y  "
            write(fu,'(a,$)') "   Bcon-Z  "
            write(fu,'(a,$)') "   Vol-MT  "
            write(fu,'(a,$)') "   Vol-VP  "
            write(fu,'(a,$)') " LIZ "
            write(fu,'(a,$)') " Energy     "
            write(fu,'(a,$)') "      Exch_energy "
            write(fu,'(a,$)') "     Pressure "
            write(fu,'(a,$)') "     Charge-MT "
            write(fu,'(a,$)') " Charge-VP "
            write(fu,'(a,$)') "  Moment-MT "
            write(fu,'(a,$)') " Moment-VP "
            write(fu,'(a,$)') "RMS-Rho "
            write(fu,'(a,$)') "  RMS-Pot "
            write(fu,'(a,$)') " RMS-Evec "
            write(fu,'(a)')   " RMS-Bcon "
            write(fu,'(a,$)') '#'
            write(fu,'(a)') repeat('-',335)
         else if ( n_spin_pola==2 ) then
            write(fu,'(a,i5)') "# Column Legend : ", 20
            write(fu,'(a)') "#         1. Type"
            write(fu,'(a)') "#         2. Index"
            write(fu,'(a)') "#         3. Position-X"
            write(fu,'(a)') "#         4. Position-Y"
            write(fu,'(a)') "#         5. Position-Z"
            write(fu,'(a)') "#         6. Force-X"
            write(fu,'(a)') "#         7. Force-Y"
            write(fu,'(a)') "#         8. Force-Z"
            write(fu,'(a)') "#         9. Volume - MT"
            write(fu,'(a)') "#        10. Volume - VP"
            write(fu,'(a)') "#        11. LIZ Size"
            write(fu,'(a)') "#        12. Energy"
            write(fu,'(a)') "#        13. Exch_energy"
            write(fu,'(a)') "#        14. Pressure"
            write(fu,'(a)') "#        15. Charge-MT"
            write(fu,'(a)') "#        16. Charge-VP"
            write(fu,'(a)') "#        17. Moment-MT"
            write(fu,'(a)') "#        18. Moment-VP"
            write(fu,'(a)') "#        19. RMS-Charge"
            write(fu,'(a)') "#        20. RMS-Potential"
            write(fu,'(a,$)') '#'
            write(fu,'(a)')repeat('-',216)
            write(fu,'(a1,a,$)') '#'," Type "
            write(fu,'(a,$)') " Indx "
            write(fu,'(a,$)') "     X     "
            write(fu,'(a,$)') "     Y     "
            write(fu,'(a,$)') "     Z     "
            write(fu,'(a,$)') "  Force-X  "
            write(fu,'(a,$)') "  Force-Y  "
            write(fu,'(a,$)') "  Force-Z  "
            write(fu,'(a,$)') "  Vol-MT   "
            write(fu,'(a,$)') "  Vol-VP   "
            write(fu,'(a,$)') " LIZ "
            write(fu,'(a,$)') " Energy     "
            write(fu,'(a,$)') "      Exch_energy "
            write(fu,'(a,$)') "     Pressure "
            write(fu,'(a,$)') "     Charge-MT "
            write(fu,'(a,$)') " Charge-VP "
            write(fu,'(a,$)') "  Moment-MT "
            write(fu,'(a,$)') " Moment-VP "
            write(fu,'(a,$)') "RMS-Rho "
            write(fu,'(a)')   "  RMS-Pot "
            write(fu,'(a,$)') '#'
            write(fu,'(a)')repeat('-',216)
         else
            write(fu,'(a,i5)') "# Column Legend : ", 17
            write(fu,'(a)') "#         1. Type"
            write(fu,'(a)') "#         2. Index"
            write(fu,'(a)') "#         3. Position-X"
            write(fu,'(a)') "#         4. Position-Y"
            write(fu,'(a)') "#         5. Position-Z"
            write(fu,'(a)') "#         6. Force-X"
            write(fu,'(a)') "#         7. Force-Y"
            write(fu,'(a)') "#         8. Force-Z"
            write(fu,'(a)') "#         9. Volume - MT"
            write(fu,'(a)') "#        10. Volume - VP"
            write(fu,'(a)') "#        11. LIZ Size"
            write(fu,'(a)') "#        12. Energy"
            write(fu,'(a)') "#        13. Pressure"
            write(fu,'(a)') "#        14. Charge-MT"
            write(fu,'(a)') "#        15. Charge-VP"
            write(fu,'(a)') "#        16. RMS-Charge"
            write(fu,'(a)') "#        17. RMS-Potential"
            write(fu,'(a,$)') '#'
            write(fu,'(a)')repeat('-',179)
            write(fu,'(a1,a,$)') '#'," Type "
            write(fu,'(a,$)') " Indx "
            write(fu,'(a,$)') "     X     "
            write(fu,'(a,$)') "     Y     "
            write(fu,'(a,$)') "     Z     "
            write(fu,'(a,$)') "  Force-X  "
            write(fu,'(a,$)') "  Force-Y  "
            write(fu,'(a,$)') "  Force-Z  "
            write(fu,'(a,$)') "  Vol-MT   "
            write(fu,'(a,$)') "  Vol-VP   "
            write(fu,'(a,$)') " LIZ "
            write(fu,'(a,$)') "      Energy     "
            write(fu,'(a,$)') "     Pressure "
            write(fu,'(a,$)') " Charge-MT "
            write(fu,'(a,$)') " Charge-VP "
            write(fu,'(a,$)') "  RMS-Rho "
            write(fu,'(a)')   "  RMS-Pot "
            write(fu,'(a,$)') '#'
            write(fu,'(a)')repeat('-',179)
         endif
      endif
   else
      close(fu)
   endif
!  ==================================================================
   end subroutine driverSystemMovie
