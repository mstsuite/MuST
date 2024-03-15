!  *******************************************************************
!  Subroutine for calculating the superconducting transition temperature
!  of metals using Gaspari-Gyorffy formula,
!  Input:
!         kmax_phi_max  -- The size of the 1st dimension of partial DOS and
!                          phase shift arrays
!         kmax_phi      -- Total number of the (l,m) indexes = (lmax+1)^2
!         n_spin_pola   -- Spin size = 1 or 2
!         atomic_number -- Atomic number
!         ss_dos_mt     -- The DOS (within muffin-tin sphere) of the atom
!         dos_mt        -- The DOS (within muffin-tin sphere) of the atom in crystal
!         dos_ws        -- The DOS (within atomic cell) of the atom in crystal
!         phase_shift   -- The partial phase shift of the atom
!         ss_pdos_mt    -- The partial DOS (within muffin-tin sphere) of the atom
!         partial_dos_mt-- The partial DOS (within muffin-tin sphere) of the atom
!                          in crystal
!         iprint        -- Printing switch = 0 or 1
!
!  Other parameters: For 1 <= kl <= kmax_phi, lofk(kl) gives angular momentum
!                    quantum number l, and mofk(kl) gives magnetic quantum
!                    number m
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine gaspari_gyorffy_formula(LocalNumAtoms,n_spin_pola,efermi,  &
                                           PartialDOS,iprint)
!  ===========================================================================
        use KindParamModule, only : IntKind, RealKind, CmplxKind
        use MathParamModule, only : ZERO, HALF, ONE, TWO, PI, TEN2m6, CONE, SQRTm1
        use PhysParamModule, only : Boltzmann, Ryd2eV, Bohr2Angstrom,         &
                                    LightSpeed, MassUnit2Ryd, Kelvin2Ryd, MassUnit2eV
        use PublicTypeDefinitionsModule, only : PDOSStruct
        use IntegerFactorsModule, only : lofk, mofk
        use InputModule, only : getKeyValue, getKeyLabelIndexValue
        use ChemElementModule, only : getDebyeTemperature, getAtomicMass, MaxLenOfAtomName
        use ChemElementModule, only : getName
        use MPPModule, only : MyPE
        use AtomModule, only : getLocalNumSpecies, getLocalAtomicNumber,      &
                               getLocalSpeciesContent, getLocalAtomName
        use GroupCommModule, only : GlobalSumInGroup, getGroupID, getMyPEinGroup
        use GroupCommModule, only : isGroupExisting
        use SystemModule, only : getNumAtomTypes, getNumAtoms, getAtomType, getAtomTypeName
        use SystemModule, only : getAtomicNumber
        use ErrorHandlerModule, only : ErrorHandler
        use Atom2ProcModule, only : getGlobalIndex
        use ScfDataModule, only : isKKRCPA
!
        implicit none
!
        integer (kind=IntKind), intent(in) :: LocalNumAtoms
        integer (kind=IntKind), intent(in) :: n_spin_pola
        integer (kind=IntKind), intent(in) :: iprint
        integer (kind=IntKind) :: atomic_number, aGID, bGID
        integer (kind=IntKind) :: id, ia, is, l, m, kl, klp1, kmax_phi, lmax_kkr
        integer (kind=IntKind) :: kmax_kkr, kmax_kkr_max
        integer (kind=IntKind) :: NumAtomTypes, NumAtoms, ig
!
        real (kind=RealKind), intent(in) :: efermi
!
        type (PDOSStruct), intent(in) :: PartialDOS(LocalNumAtoms)
!
        real (kind=RealKind), pointer :: dos_mt(:)
        real (kind=RealKind), pointer :: dos_ws(:)
        real (kind=RealKind), pointer :: phase_shift(:,:)
        real (kind=RealKind), pointer :: ss_pdos_mt(:,:)
        real (kind=RealKind), pointer :: partial_dos_mt(:,:)
!
        real (kind=RealKind) :: I_Total,sinTerm,Inum,Iden
        real (kind=RealKind), allocatable :: Phonon_Freq2(:)
        real (kind=RealKind) :: EPC, EPCnum,EPCden, EPCnumUP, EPCnumDOWN
        real (kind=RealKind) :: AtomMass
        real (kind=RealKind) :: Coulomb,Coulombnum,Coulombden,mu_star
        real (kind=RealKind) :: SuperTemp,SuperTempExp
        real (kind=RealKind) :: SumI,Iconst,IconstUp,IconstDown,SumIUP,SumIDOWN
        real (kind=RealKind) :: DebyeTemp
        real (kind=RealKind) :: species_content
        real (kind=RealKind) :: total_dos, dos_per_spin, total_dos_mt, dos_mt_per_spin, Ef
        real (kind=RealKind) :: spdos, cpdos,pps, cfac, val_au, av_phonon_freq, eta
        real (kind=RealKind) :: nr(0:4)  ! assuming lmax <= 4
!
        character (len=MaxLenOfAtomName), pointer :: AtomTypeName(:)
        integer (kind=IntKind), pointer :: AtomType(:)
        integer (kind=IntKind), allocatable :: ind_array(:)
        real (kind=RealKind), allocatable :: val_array(:)
!
        complex (kind=CmplxKind) :: kappa
        complex (kind=CmplxKind), pointer :: t_mat(:,:)
        complex (kind=CmplxKind), allocatable :: phase_shift_aver(:,:)
!
        NumAtoms = getNumAtoms()
        NumAtomTypes = getNumAtomTypes()
        AtomTypeName => getAtomTypeName()
        AtomType => getAtomType()
!
        allocate(ind_array(NumAtoms), val_array(NumAtoms))
        allocate(Phonon_Freq2(NumAtoms))
!
        if (n_spin_pola == 1) then
           cfac = HALF
        else
           cfac = ONE
        endif
!
        if (MyPE == 0 .and. iprint >= 0) then
           write(6,*), ' '
           write(6,*), ' '
           write(6,*), "***************************************************"
           write(6,*),         'Output of Gaspari-Gyorffy Formula'
           write(6,*), 'Calculate the Superconducting Transition Teperature'
           write(6,*), "***************************************************"
           write(6,*), 'n_spin_pola: ', n_spin_pola
           do id = 1, LocalNumAtoms
              write(6,*), 'id: ', id
              write(6,*), 'kmax_phi: ', PartialDOS(id)%kmax_phi
              do ia = 1, getLocalNumSpecies(id)
                 write(6,*), 'ia: ', ia
                 write(6,*), 'atomic number: ', getLocalAtomicNumber(id,ia)
                 if (n_spin_pola == 1) then
                    write(6,*), 'dos_ws per spin: (St/Ryd)', PartialDOS(id)%dos_ws(1,ia)*HALF
                    write(6,*), 'dos_ws per spin: (St/eV)', PartialDOS(id)%dos_ws(1,ia)*HALF/Ryd2eV
                 else
                    write(6,*), 'dos_ws for spin-up: (St/Ryd)', PartialDOS(id)%dos_ws(1,ia)
                    write(6,*), 'dos_ws for spin-up: (St/eV)', PartialDOS(id)%dos_ws(1,ia)/Ryd2eV
                    write(6,*), 'dos_ws for spin-down: (St/Ryd)', PartialDOS(id)%dos_ws(2,ia)
                    write(6,*), 'dos_ws for spin-down: (St/eV)', PartialDOS(id)%dos_ws(2,ia)/Ryd2eV
                 endif
              enddo
           enddo
           write(6,*), 'Efermi: ', efermi
           !do kl=1,kmax_phi
           !   write(6,*), 'kl: ', kl, ' l: ', lofk(kl), ' m: ', mofk(kl)
           !   write(6,*), 'phaseshift: ', phase_shift(kl,n_spin_pola)
           !   write(6,*), 'dos: ', partial_dos_mt(kl,n_spin_pola)
           !   write(6,*), 'ss_dos: ', ss_pdos_mt(kl,n_spin_pola)
           !enddo
        endif
!
!       Determine the Debye Temperature and phonon frequency
!       ==============================================================
        if (getKeyValue(1,'Average of phonon frequency (1/sec)',av_phonon_freq,default_param=.false.) == 0) then
           DebyeTemp = av_phonon_freq/Boltzmann
           write(6,*), 'Debye 1'
        else if (getKeyValue(1,'Average of phonon frequency (K)',av_phonon_freq,default_param=.false.) == 0) then
           DebyeTemp = av_phonon_freq
           write(6,*), 'Debye 2'
        else
           if (getKeyValue(1,'Debye Temperature (K)',DebyeTemp, default_param=.false.) /= 0) then
              write(6,*), 'Debye 3'
              if (NumAtomTypes == 1)  then
                 write(6,*), 'NumAtomTypes is 1'
                 atomic_number = getLocalAtomicNumber(1,1)
                 DebyeTemp = getDebyeTemperature(atomic_number)
              !!! I think this else statement is not properly indented
              else
                 call ErrorHandler('Gaspari-Gyorffy-Formula','Debye temperature is missing from input')
              endif
           endif
        endif
!
!       Determine the average phonon frequency squared
!       ==============================================================
        ind_array = 0; val_array = ZERO; Phonon_Freq2 = ZERO
        if (getKeyLabelIndexValue(1,'Atomic mass times <omega^2> (eV/Anst^2)', &
                                  NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
           do ig = 1, NumAtoms
              val_au = val_array(ig)*Bohr2Angstrom**2/Ryd2eV ! This is M*<omega^2> in atomic units
              atomic_number = getAtomicNumber(ig)
              if (MyPE == 0) then
                 write(6,*) getName(atomic_number), ': M*<omega^2> = ',val_array(ig),'(eV/Anst^2)'
              endif
              AtomMass=getAtomicMass(atomic_number)*MassUnit2Ryd/LightSpeed**2 ! Atomic Mass in Ryd.
              Phonon_Freq2(ig) = val_au/AtomMass
           enddo
        else if (getKeyLabelIndexValue(1,'Atomic mass times <omega^2> (Ryd/BohrRad^2)', &
                                  NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
           do ig = 1, NumAtoms
              atomic_number = getAtomicNumber(ig)
              if (MyPE == 0) then
                 write(6,*) getName(atomic_number), ': M*<omega^2> = ',val_array(ig),'(Ryd/BohrRad^2)'
              endif
              AtomMass=getAtomicMass(atomic_number)*MassUnit2Ryd/LightSpeed**2 ! Atomic Mass in Ryd.
              Phonon_Freq2(ig) = val_array(ig)/AtomMass
           enddo
        else
           if (getKeyLabelIndexValue(1,'Average of phonon frequency squared (1/sec^2)', &
                                     NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
              do ig = 1, NumAtoms
                 atomic_number = getAtomicNumber(ig)
                 if (MyPE == 0) then
                    write(6,*) getName(atomic_number), ': <omega^2> = ',val_array(ig),'(1/sec^2)'
                 endif
                 Phonon_Freq2(ig) = val_array(ig)
              enddo
           else if (getKeyLabelIndexValue(1,'Average of phonon frequency squared (K^2)', &
                                          NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
              do ig = 1, NumAtoms
                 atomic_number = getAtomicNumber(ig)
                 if (MyPE == 0) then
                    write(6,*) getName(atomic_number), ': <omega^2> = ',val_array(ig),'(K^2)'
                 endif
                 Phonon_Freq2(ig) = val_array(ig)*Kelvin2Ryd**2
              enddo
           else
              Phonon_Freq2 = HALF*(DebyeTemp*Kelvin2Ryd)**2
           endif
        endif
        do ig = 1, NumAtoms
           if (Phonon_Freq2(ig) < TEN2m6) then
              call ErrorHandler('Gaspari-Gyorffy-Formula','Invalid value of <omega^2>',Phonon_Freq2(ig))
           endif
        enddo
!
        deallocate(ind_array, val_array)
!
        if (MyPE == 0) then
           write(6,*), 'Writing Debye Temp'
           write(6,*), 'Debye Temp (K): ', DebyeTemp
           do ig = 1, NumAtoms  ! Loop atomic sites
           !!! Want to use ia to get specific atomic number
              do ia = 1, NumAtomTypes
                 write(6,*), 'gettting atomic number', ig, ia
                 atomic_number = getAtomicNumber(ig,ia)
              enddo
              write(6,'(1x,a,a,d15.8,a)')getName(atomic_number),': Average of Phonon Frequency Square =', &
                                         Phonon_Freq2(ig),' (Ryd^2)'
!             write(6,'(40x,a,d15.8,a)')'=',Phonon_Freq2(ig)/Kelvin2Ryd**2,' (K^2)'
           enddo
        endif
!
        aGID = getGroupID("Unit Cell")
        if (isGroupExisting('E-K Plane')) then
           bGID = getGroupID('E-K Plane')
        else
           bGID = getGroupID('Energy Mesh')
        endif
!
!       ==============================================================
!       Calculate total dos
!       ==============================================================
        total_dos = ZERO
        total_dos_mt = ZERO
        do id = 1, LocalNumAtoms  ! Loop atomic sites
           do ia = 1, getLocalNumSpecies(id)  ! Loop over atomic species
              dos_ws => PartialDOS(id)%dos_ws(:,ia)
              dos_mt => PartialDOS(id)%dos_mt(:,ia)
              species_content = getLocalSpeciesContent(id,ia)
              do is = 1, n_spin_pola
!                =====================================================
!                For now, I am using the muffin-tin DOS. It needs to 
!                be updated.
!                =====================================================
                 total_dos = total_dos + species_content*dos_ws(is)
                 total_dos_mt = total_dos_mt + species_content*dos_mt(is)
              enddo
           enddo
        enddo
        call GlobalSumInGroup(aGID,total_dos)
        dos_per_spin = HALF*total_dos  ! states/Ryd/spin (has both spin up and down)
        dos_mt_per_spin = HALF*total_dos_mt  ! states/Ryd/spin
        if (MyPE == 0) then
           write(6,'(/,1x,a,f10.5,a)')'WS-Volume DOS of Unit cell per spin =', dos_per_spin,' (states/Ryd/spin))'
           write(6,'(  1x,a,f10.5,a)')'                                    =', dos_per_spin/Ryd2eV,' (states/eV/spin)'
           write(6,'(  1x,a,f10.5,a)')'MT-Volume DOS of Unit cell per spin =', dos_mt_per_spin,' (states/Ryd/spin))'
           write(6,'(  1x,a,f10.5,a)')'                                    =', dos_mt_per_spin/Ryd2eV,' (states/eV/spin)'
        endif
!
!       Calculate Lamda (or EPC) ...........
!       ==============================================================
        kmax_kkr_max = 0
        do id = 1, LocalNumAtoms
           kmax_kkr_max = max(kmax_kkr_max,PartialDOS(id)%kmax_kkr)
        enddo
        if (isKKRCPA()) then
           allocate(phase_shift_aver(kmax_kkr_max,n_spin_pola))
        endif
        kappa = sqrt(efermi)
        nr = ZERO
        EPC = ZERO
        do id = 1, LocalNumAtoms  ! Loop atomic sites
           if (isKKRCPA()) then
              do is = 1, n_spin_pola
                 t_mat => PartialDOS(id)%t_aver(:,:,is)
                 do kl = 1, PartialDOS(id)%kmax_kkr
                     phase_shift_aver(kl,is) = atan(kappa*t_mat(kl,kl)/(-CONE+SQRTm1*kappa*t_mat(kl,kl)))
                 enddo
              enddo
           endif
           kmax_phi = PartialDOS(id)%kmax_phi
           kmax_kkr = PartialDOS(id)%kmax_kkr
           lmax_kkr = lofk(kmax_phi)
           ig = getGlobalIndex(id)
           do ia = 1, getLocalNumSpecies(id)  ! Loop over atomic species
                                              ! at each atomic site. The
                                              ! number of species is usually 1,
                                              ! except for the KKR-CPA case.
              species_content = getLocalSpeciesContent(id,ia)
              atomic_number = getLocalAtomicNumber(id,ia)
              phase_shift => PartialDOS(id)%phase_shift(:,:,ia)
              ss_pdos_mt => PartialDOS(id)%ss_pdos_mt(:,:,ia)
              partial_dos_mt => PartialDOS(id)%partial_dos_mt(:,:,ia)
!
              if (getMyPEinGroup(bGID) == 0) then
                 write(6,'(/,1x,a,a)')'For species: ',getLocalAtomName(id,ia)
              endif
!
              SumI = ZERO
              SumIUP = ZERO
              SumIDOWN = ZERO
              do is = 1, n_spin_pola
                 if (n_spin_pola == 2 .and. getMyPEinGroup(bGID) == 0) then
                    write(6,'(1x,a,i4)')'spin index = ',is
                 endif
                 do l=0,lmax_kkr
                    cpdos = ZERO
                    spdos = ZERO
                    do m = -l, l
                       kl = (l+1)**2-l+m
                       cpdos = cpdos + partial_dos_mt(kl,is)
                       spdos = spdos + ss_pdos_mt(kl,is)
                    enddo
                    nr(l) = cpdos/spdos
                    cpdos = cfac*cpdos
                    spdos = cfac*spdos
                    if (getMyPEinGroup(bGID) == 0) then
                       kl = (l+1)**2-l
                       if (phase_shift(kl,is) > PI*HALF) then
                          pps = phase_shift(kl,is) - PI
                       else if (phase_shift(kl,is) < -PI*HALF) then
                          pps = phase_shift(kl,is) + PI
                       else
                          pps = phase_shift(kl,is)
                       endif
                       if (l == 0) then
                          write(6,'(1x,5(a,f12.5))')'s-state: phase shift =',pps, &
                                                    ', partial DOS =',cpdos, ', DOS ratio =',nr(l), ', f_l = ', &
                                                    cpdos/(PartialDOS(id)%dos_ws(is,ia))
                       else if (l == 1) then
                          write(6,'(1x,5(a,f12.5))')'p-state: phase shift =',pps, &
                                                    ', partial DOS =',cpdos, ', DOS ratio =',nr(l), ', f_l = ', &
                                                    cpdos/(PartialDOS(id)%dos_ws(is,ia))
                       else if (l == 2) then
                          write(6,'(1x,5(a,f12.5))')'d-state: phase shift =',pps, &
                                                    ', partial DOS =',cpdos, ', DOS ratio =',nr(l),', f_l = ', &
                                                    cpdos/(PartialDOS(id)%dos_ws(is,ia))
                       else if (l == 3) then
                          write(6,'(1x,5(a,f12.5))')'f-state: phase shift =',pps, &
                                                    ', partial DOS =',cpdos, ', DOS ratio =',nr(l),', f_l = ', &
                                                    cpdos/(PartialDOS(id)%dos_ws(is,ia))
                       endif
                    endif
                 enddo
                 do l=0,lmax_kkr-1 ! sum over l
                    klp1 = (l+2)**2-l-1
                    kl = (l+1)**2-l
                    sinTerm = sin(phase_shift(klp1,is)-phase_shift(kl,is))
                    if (n_spin_pola == 1) then
                       SumI = SumI + TWO*(l+1)*sinTerm**2*nr(l+1)*nr(l)
                    else
                       if (is == 1) then
                       !  SumIUP = SumIUP + TWO*(l+1)*sinTerm**2*nr(l+1)*nr(l)
                          SumIUP = SumIUP + (l+1)*sinTerm**2*nr(l+1)*nr(l)
                       else
                       !  SumIDOWN = SumIDOWN + TWO*(l+1)*sinTerm**2*nr(l+1)*nr(l)
                          SumIDOWN = SumIDOWN + (l+1)*sinTerm**2*nr(l+1)*nr(l)
                       endif
                    endif
                 enddo
              enddo
!
              if (n_spin_pola == 1) then
                 Iconst = efermi/PI**2/dos_per_spin**2 ! Ryd^3
                 I_Total = Iconst*SumI ! Ryd^2/au^2
              else
                 IconstUp = efermi/PI**2/PartialDOS(id)%dos_ws(1,ia)**2 ! Ryd^3
                 IconstDown = efermi/PI**2/PartialDOS(id)%dos_ws(2,ia)**2 ! Ryd^3
              endif
              Iconst = efermi/PI**2/dos_per_spin**2 ! Ryd^3
!
              AtomMass=getAtomicMass(atomic_number)*MassUnit2Ryd/LightSpeed**2
              EPCden = AtomMass*Phonon_Freq2(ig) ! Ryd/au^2
!
              if (n_spin_pola == 1) then
                 EPCnum = dos_per_spin*I_Total  ! Ryd/au^2
                 eta = EPCnum
              else
             !!  EPCnumUP = PartialDOS(id)%dos_ws(1,ia)*IconstUp*SumIUP
             !!  EPCnumDOWN = PartialDOS(id)%dos_ws(2,ia)*IconstDown*SumIDOWN
                 EPCnumUP = dos_per_spin*Iconst*SumIUP
                 EPCnumDOWN = dos_per_spin*Iconst*SumIDOWN
                 eta = EPCnumUp+EPCnumDown
              endif
!
              EPC = EPC + getLocalSpeciesContent(id,ia)*eta/EPCden ! unitless
!
              if (getMyPEinGroup(bGID) == 0) then
                 write(6,'(1x,a,f12.5)')'AtomicMass (Ryd/c^2) =', AtomMass
                 write(6,'(1x,a,f12.5)')'M<Omega^2> (Ryd/au^2)=', EPCden
                 write(6,'(1x,a,f12.5)')'M<Omega^2> (eV/A^2)  =', EPCden*Ryd2eV/Bohr2Angstrom**2
                 if (n_spin_pola == 1) then
                    write(6,'(1x,a,f12.5)')'SumI                 =', SumI
                    write(6,'(1x,a,f12.5)')'I (eV^2/A^2)         =', I_Total*(Ryd2eV/Bohr2Angstrom)**2
                    write(6,'(1x,a,f12.5)')'I (Ryd^2/au^2)       =', I_Total
                    write(6,'(1x,a,f12.5)')'eta (Ryd/au^2)       =', EPCnum
                    write(6,'(1x,a,f12.5)')'eta (eV/A^2)         =', EPCnum*Ryd2eV/Bohr2Angstrom**2
                    write(6,'(1x,a,f12.5)')'eta/(M<Omega^2>      =', EPCnum/EPCden
                 else
                    write(6,'(1x,a,f12.5)')'eta up (Ryd/au^2)    =', EPCnumUP
                    write(6,'(1x,a,f12.5)')'eta up (eV/A^2)      =', EPCnumUP*Ryd2eV/Bohr2Angstrom**2
                    write(6,'(1x,a,f12.5)')'eta down (Ryd/au^2)  =', EPCnumDOWN
                    write(6,'(1x,a,f12.5)')'eta down (eV/A^2)    =', EPCnumDOWN*Ryd2eV/Bohr2Angstrom**2
                    write(6,'(1x,a,f12.5)')'eta (Ryd/au^2)       =', eta
                    write(6,'(1x,a,f12.5)')'eta (eV/A^2)         =', eta*Ryd2eV/Bohr2Angstrom**2
                    write(6,'(1x,a,f12.5)')'eta/(M<Omega^2>      =', eta/EPCden
                 endif
                 write(6,'(1x,a,f12.5)')'lamda                =', eta/EPCden
              endif
           enddo  ! Loop over ia
        enddo  ! Loop over id
        call GlobalSumInGroup(aGID, EPC)
        if (MyPE == 0) then
           write(6,'(/,1x,a)')'For the system ...'
           write(6,'(1x,a,f12.5)')'lamda (EPC)                  =', EPC
        endif
!
!       Calculate mu ...
!       ==============================================================
        Coulombnum = 0.26D0*total_dos/Ryd2eV
        Coulombden = ONE+total_dos/Ryd2eV
        Coulomb = Coulombnum/Coulombden
        if (MyPE == 0) then
           write(6,'(1x,a,f12.5)')'mu* (Coulomb pseudopotential)=', Coulomb
        endif
        if (getKeyValue(1,'mu* (e-e interaction constant)',mu_star, default_param=.false.) == 0) then
           Coulomb = mu_star
           if (MyPE == 0) then
              write(6,'(1x,a,f12.5)')'Using the input mu* (Coulomb pseudopotential): ', Coulomb
           endif
        endif
!
!       if (EPC-Coulomb*(ONE+0.62D0*EPC) < ZERO) then
!          if (MyPE == 0) then
!             write(6,*), 'Warning: EPC-Coulomb*(ONE+0.62D0*EPC) < 0'
!             write(6,*), 'EPC/(ONE+0.62D0*EPC) =', EPC/(ONE+0.62D0*EPC)
!          endif
!       endif
!
!       Calculate Tc ...
!       ==============================================================
        if (EPC > Coulomb) then
           if (MyPE == 0) then
              write(6,'(1x,a,f12.5)')'EPC-Coulomb*(ONE+0.62D0*EPC) =', EPC-Coulomb*(ONE+0.62D0*EPC)
           endif
           SuperTempExp = -1.04D0*(ONE+EPC)/(EPC-Coulomb*(ONE+0.62D0*EPC))
           SuperTemp=DebyeTemp*exp(SuperTempExp)/1.45D0
        else
           if (MyPE == 0) then
              write(6,'(1x,a)') 'EPC <= Coulomb'
           endif
           SuperTemp=ZERO
        endif
        if (MyPE == 0) then
           write(6,'(/,1x,a)')'************************************************'
           write(6,'(1x,a,f12.5)')'Superconducting Transition Temp (K):',SuperTemp
           write(6,'(1x,a,/)')'************************************************'
        endif
!
        deallocate(Phonon_Freq2)
!
        if (isKKRCPA()) then
           deallocate(phase_shift_aver)
        endif
!
!       ==============================================================
!       Using the partial phase shift and DOS data published in PRB 15, 4221
!       (1977) to check against the calculated Tc published in the paper
!       --------------------------------------------------------------
        call testPRB15_4221_1977()
!       --------------------------------------------------------------
!
        end subroutine gaspari_gyorffy_formula
!   ==================================================================
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine testPRB15_4221_1977()
!   ==================================================================
!
!       Test the data published in PRB 15, 4221 (1977)
!
!       **************************************************************
        use KindParamModule, only : IntKind, RealKind
        use MathParamModule, only : ZERO, HALF, ONE, TWO, PI
        use PhysParamModule, only : Boltzmann, Ryd2eV, Bohr2Angstrom,         &
                                    LightSpeed, MassUnit2Ryd, Kelvin2Ryd
        use MPPModule, only : MyPE
!
        integer (kind=IntKind) :: atomic_number
!
        real (kind=RealKind) :: nr(0:4)  ! assuming lmax <= 4
        real (kind=RealKind) :: ps(0:4)  ! assuming lmax <= 4
        real (kind=RealKind) :: DebyeTemp, Ef, dos_per_spin
        real (kind=RealKind) :: EPCnum,EPC,Coulomb,SuperTemp
!
        if (MyPE == 0) then
           write(6,'(a)')  '............................................'
           write(6,'(/,a)')'Check against the data for Nb, Ti, and V published in PRB 15, 4221 (1977)'
           write(6,'(a)')  '............................................'
        endif
!       ==============================================================
!       For Nb: Input from Table II and Table IV
!       --------------------------------------------------------------
        atomic_number = 41
        ps(0) = -0.932; ps(1) = -0.363; ps(2) = 1.142; ps(3) = 0.007
        nr(0) = 0.807; nr(1) = 2.927; nr(2) = 0.681; nr(3) = 3.845
        DebyeTemp = 275.0D0; Ef = 0.676; dos_per_spin = 9.71
!
        call calculateMcMillanTc(atomic_number,ps,nr,DebyeTemp,Ef,    &
                                 dos_per_spin,EPCnum,EPC,Coulomb,SuperTemp)
!
        if (MyPE == 0) then
           write(6,'(/,a)')'For Niobium'
           write(6,'(a,f12.5)')'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
           write(6,'(a,f12.5)')'lamda = ',EPC
           write(6,'(a,f12.5)')'Mu* (Coulomb pseudopotential): ', Coulomb
           write(6,'(a,f12.5)')'Superconducting Transition Temp: ', SuperTemp
        endif
!       ==============================================================
!       For Ti: Input from Table II and Table IV
!       --------------------------------------------------------------
        atomic_number = 22
        ps(0) = -0.641; ps(1) = -0.157; ps(2) = 0.754; ps(3) = 0.003
        nr(0) = 0.414; nr(1) = 3.101; nr(2) = 0.788; nr(3) = 6.382
        DebyeTemp = 420.0D0; Ef = 0.588; dos_per_spin = 12.38
!
        call calculateMcMillanTc(atomic_number,ps,nr,DebyeTemp,Ef,    &
                                 dos_per_spin,EPCnum,EPC,Coulomb,SuperTemp)
!
        if (MyPE == 0) then
           write(6,'(/,a)')'For Titanium'
           write(6,'(a,f12.5)')'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
           write(6,'(a,f12.5)')'lamda = ',EPC
           write(6,'(a,f12.5)')'Mu* (Coulomb pseudopotential): ', Coulomb
           write(6,'(a,f12.5)')'Superconducting Transition Temp: ', SuperTemp
        endif
!       ==============================================================
!       For V: Input from Table II and Table IV
!       --------------------------------------------------------------
        atomic_number = 23
        ps(0) = -0.691; ps(1) = -0.173; ps(2) = 1.030; ps(3) = 0.003
        nr(0) = 0.550; nr(1) = 2.842; nr(2) = 0.634; nr(3) = 6.935
        DebyeTemp = 380.0D0; Ef = 0.675; dos_per_spin = 12.70
!
        call calculateMcMillanTc(atomic_number,ps,nr,DebyeTemp,Ef,    &
                                 dos_per_spin,EPCnum,EPC,Coulomb,SuperTemp)
!
        if (MyPE == 0) then
           write(6,'(/,a)')'For Vanadium'
           write(6,'(a,f12.5)')'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
           write(6,'(a,f12.5)')'lamda = ',EPC
           write(6,'(a,f12.5)')'Mu* (Coulomb pseudopotential): ', Coulomb
           write(6,'(a,f12.5)')'Superconducting Transition Temp: ', SuperTemp
        endif
!       ==============================================================
!       For K: Input from Table II and Table IV
!       --------------------------------------------------------------
        atomic_number = 19
        ps(0) = -0.165; ps(1) = -0.026; ps(2) = 0.033; ps(3) = 0.001
        nr(0) = 1.200; nr(1) = 1.178; nr(2) = 1.037; nr(3) = 0.934
        DebyeTemp = 91.0D0; Ef = 0.159; dos_per_spin = 5.26
!
        call calculateMcMillanTc(atomic_number,ps,nr,DebyeTemp,Ef,    &
                                 dos_per_spin,EPCnum,EPC,Coulomb,SuperTemp)
!
        if (MyPE == 0) then
           write(6,'(/,a)')'For Potassium'
           write(6,'(a,f12.5)')'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
           write(6,'(a,f12.5)')'lamda = ',EPC
           write(6,'(a,f12.5)')'Mu* (Coulomb pseudopotential): ', Coulomb
           if (SuperTemp < 100000.0D0) then
              write(6,'(a,f12.5)')'Superconducting Transition Temp: ', SuperTemp
           else
              write(6,'(a,d12.5)')'Superconducting Transition Temp: ', SuperTemp
           endif
        endif
!
    end subroutine testPRB15_4221_1977
!   ==================================================================
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calculateMcMillanTc(anum,ps,nr,thetaD,Ef,dos_per_spin, &
                                   eta,lamda,mu,Tc)
!   ==================================================================
    use KindParamModule, only : IntKind, RealKind
    use ChemElementModule, only : getAtomicMass
    use MathParamModule, only : ZERO, HALF, ONE, TWO, PI
    use PhysParamModule, only : Boltzmann, Ryd2eV, Bohr2Angstrom,     &
                                LightSpeed, MassUnit2Ryd, Kelvin2Ryd
    implicit none
!
    real (kind=RealKind), intent(in) :: ps(0:3),nr(0:3)
    real (kind=RealKind), intent(in) :: thetaD,Ef,dos_per_spin
    real (kind=RealKind), intent(out) :: eta,lamda,mu,Tc
!
    real (kind=RealKind) :: Phonon_Freq2, Iconst, SumI, sinTerm, I_Total
    real (kind=RealKind) :: AtomMass, EPCden, total_dos
    real (kind=RealKind) :: Coulombnum, Coulombden, SuperTempExp
!
    integer (kind=IntKind), intent(in) :: anum
    integer (kind=IntKind) :: l
!
    Phonon_Freq2 = HALF*(thetaD*Kelvin2Ryd)**2
    Iconst = Ef/PI**2/dos_per_spin**2 ! Ryd^3
    total_dos = TWO*dos_per_spin ! dos_per_spin = DOS per spin
    AtomMass=getAtomicMass(anum)*MassUnit2Ryd/LightSpeed**2
!
    SumI = ZERO
    do l=0,2 ! sum over l
       sinTerm = sin(ps(l+1)-ps(l))
       SumI = SumI + TWO*(l+1)*sinTerm**2*nr(l+1)*nr(l)
    enddo
    I_Total = Iconst*SumI ! Ryd^2/au^2
!
    eta = dos_per_spin*I_Total  ! Ryd/au^2
    EPCden = AtomMass*Phonon_Freq2 ! Ryd/au^2
    lamda = eta/EPCden ! unitless
!
    Coulombnum = 0.26D0*total_dos/Ryd2eV
    Coulombden = ONE+total_dos/Ryd2eV
    mu = Coulombnum/Coulombden
!
    SuperTempExp = -1.04D0*(ONE+lamda)/(lamda-mu*(ONE+0.62D0*lamda))
    if (mu < lamda) then
       Tc=thetaD*exp(SuperTempExp)/1.45D0
    else
       Tc = ZERO
    endif
!
    end subroutine calculateMcMillanTc
