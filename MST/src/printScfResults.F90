!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printScfResults(GlobalNumAtoms,LocalNumAtoms,              &
                              node_print_level,atom_print_level,n_write, &
                              movie,iscf,nscf,Converged)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use MPPModule, only : syncAllPEs, MyPE
!
   use ScfDataModule, only : n_spin_pola, n_spin_cant
!
   use PotentialTypeModule, only : isFullPotential, isMuffinTinFullPotential
!
   use OutputModule, only : getDensityPrintFlag, getPotentialPrintFlag
!
   use ValenceDensityModule, only : printRho_L
!
   use ChargeDensityModule, only : printChargeDensity,        &
                                   printChargeDensity_L,      &
                                   printMomentDensity_L
!
   use PotentialGenerationModule, only : printNewPot_L => printPot_L
   use PotentialGenerationModule, only : printMadelungShiftTable, printPotentialGeneration
!
   use ChargeDistributionModule, only : printChargeDistribution
!
   implicit none
!
   logical, intent(in) :: Converged
!
   character (len=12) :: anm
!
   integer (kind=IntKind), intent(in) :: GlobalNumAtoms
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: node_print_level
   integer (kind=IntKind), intent(in) :: atom_print_level(LocalNumAtoms)
   integer (kind=IntKind), intent(inout) :: n_write
   integer (kind=IntKind), intent(in) :: movie, iscf, nscf
   integer (kind=IntKind) :: id, funit
!
   call syncAllPEs()
!
   if (GlobalNumAtoms > 100) then
      funit = 21     ! if tooo many atoms, charge table is written to a file
   else
      funit = 6
   endif
!
   if ( node_print_level >= 0 ) then
!     ----------------------------------------------------------------
      call printPotentialGeneration()
!     ----------------------------------------------------------------
      call printChargeDistribution(iscf,funit,n_spin_cant)
!     ----------------------------------------------------------------
   endif
!
   n_write = n_write + 1
!
   if ( n_write==movie .or. iscf==nscf .or. Converged ) then
      if ( getDensityPrintFlag()>=1) then
         call printDensityOnGrid(LocalNumAtoms,node_print_level)
      endif
!
      if ( getPotentialPrintFlag()>=1) then
         call printPotentialOnGrid(LocalNumAtoms,node_print_level)
      endif
!
      if ( node_print_level >= 0 ) then
!        -------------------------------------------------------------
         call printMadelungShiftTable(iscf,funit)
!        -------------------------------------------------------------
         if (n_spin_pola == 2) then
!           ----------------------------------------------------------
            call printMomentVsCoreSplit(iscf,funit)
!           ----------------------------------------------------------
         endif
      endif
!     if (n_spin_cant == 2) then
!        -------------------------------------------------------------
!        mom_table => getGlobalVPCellMomentTable()
!        -------------------------------------------------------------
!        call updateSystemMovie(mom_table)
!        -------------------------------------------------------------
!        call writeMomentMovie(iscf,itstep)
!        -------------------------------------------------------------
!     endif
      n_write = 0
   endif
!
!  output the charge density and magnetic moment (added by xianglin)
!  The following codes need some cleaning....
!
   if ( isFullPotential() ) then
      if ( getDensityPrintFlag()>=2 .and. (iscf==nscf .or. Converged) ) then
         do id = 1,LocalNumAtoms
!           ----------------------------------------------------------------
            call printChargeDensity_L(id,"TotalNew")
            call printChargeDensity_L(id,"Valence")
!           ----------------------------------------------------------------
            if ( n_spin_pola==2 ) then
!              -------------------------------------------------------------
               call printMomentDensity_L(id,"Valence")
               call printMomentDensity_L(id,"TotalNew")
!              -------------------------------------------------------------
            endif
         enddo
      endif
      if ( node_print_level >=1 ) then
         do id = 1,LocalNumAtoms
            if ( atom_print_level(id) >= 0 ) then
               call printChargeDensity_L(id,"Valence")
               call printChargeDensity_L(id,"TotalNew")
               if ( n_spin_pola==2 ) then
                  call printMomentDensity_L(id,"Valence")
                  call printMomentDensity_L(id,"TotalNew")
               endif
               if ( isFullPotential().and.atom_print_level(id)>=1) then
                  call printChargeDensity_L(id,"Valence",1)
                  if (.not.isMuffinTinFullPotential()) then
                     call printChargeDensity_L(id,"Pseudo")
                     call printChargeDensity_L(id,"Pseudo",1)
                  endif
               endif
            endif
         enddo
      endif
   endif
!
   if (iscf < 10) then
      write(anm,'(a,i1,a)')'scf',iscf,'post_'
   else if (iscf < 100) then
      write(anm,'(a,i2,a)')'scf',iscf,'post_'
   else if (iscf < 1000) then
      write(anm,'(a,i3,a)')'scf',iscf,'post_'
   else
      write(anm,'(a,i4,a)')'scf',iscf,'post_'
   endif
!
   if (node_print_level >= 1) then
      do id = 1,LocalNumAtoms
!        -------------------------------------------------------------------
         call printRho_L(id,aux_name=anm)
!        -------------------------------------------------------------------
      enddo
   endif
!
   if ( isFullPotential() ) then
      if ( getPotentialPrintFlag() >= 2 .and. node_print_level >= 1 ) then
         do id = 1,LocalNumAtoms
!           call printNewPot_L(id,-1,"Coulomb",aux_name=anm)
            call printNewPot_L(id,-1,"Total",aux_name=anm)
!           call printNewPot_L(id,-1,"Exchg",aux_name=anm)
            if ( atom_print_level(id) >= 0 ) then
!              call printNewPot_L(id,-1,"Total",1,aux_name=anm)
!              call printNewPot_L(id,-1,"Coulomb",1,aux_name=anm)
!              call printNewPot_L(id,-1,"Exchg",1,aux_name=anm)
!              ----------------------------------------
!              Components of Coulomb potential
!              ----------------------------------------
!              call printNewPot_L(id,-1,"Tilda",aux_name=anm)
!              call printNewPot_L(id,-1,"Madelung",aux_name=anm)
!              if (.not.isMuffinTinFullPotential()) then
!                 call printNewPot_L(id,-1,"Pseudo",aux_name=anm)
!              endif
!              call printNewPot_L(id,-1,"Tilda",1,aux_name=anm)
!              call printNewPot_L(id,-1,"Madelung",1,aux_name=anm)
!              if (.not.isMuffinTinFullPotential()) then
!                 call printNewPot_L(id,-1,"Pseudo",1,aux_name=anm)
!              endif
            endif
         enddo
      endif
   endif
!
   end subroutine printScfResults
!  ===================================================================
