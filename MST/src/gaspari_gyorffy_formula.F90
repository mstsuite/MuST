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
        use KindParamModule, only : IntKind, RealKind
        use MathParamModule, only : ZERO, HALF, ONE, TWO, PI
        use PhysParamModule, only : Boltzmann, Ryd2eV, Bohr2Angstrom,         &
                                    LightSpeed, MassUnit2Ryd, Kelvin2Ryd
        use PublicTypeDefinitionsModule, only : PDOSStruct
        use IntegerFactorsModule, only : lofk, mofk
        use InputModule, only : getKeyValue
        use ChemElementModule, only : getDebyeTemperature, getAtomicMass
        use MPPModule, only : myPE
        use AtomModule, only : getLocalNumSpecies, getLocalAtomicNumber,      &
                               getLocalSpeciesContent, getLocalAtomName
        use GroupCommModule, only : GlobalSumInGroup, getGroupID, getMyPEinGroup
        use GroupCommModule, only : isGroupExisting
        use SystemModule, only : getNumAtomTypes, getNumAtoms
        use ErrorHandlerModule, only : ErrorHandler
        implicit none
!
        integer (kind=IntKind), intent(in) :: LocalNumAtoms
        integer (kind=IntKind), intent(in) :: n_spin_pola
        integer (kind=IntKind), intent(in) :: iprint
        integer (kind=IntKind) :: atomic_number, aGID, bGID
        integer (kind=IntKind) :: id, ia, is, l, m, kl, klp1, kmax_phi, lmax_kkr
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
        real (kind=RealKind) :: Phonon_Freq2
        real (kind=RealKind) :: EPC, EPCnum,EPCden
        real (kind=RealKind) :: AtomMass
        real (kind=RealKind) :: Coulomb,Coulombnum,Coulombden
        real (kind=RealKind) :: SuperTemp,SuperTempExp
        real (kind=RealKind) :: SumI,Iconst
        real (kind=RealKind) :: DebyeTemp
        real (kind=RealKind) :: species_content
        real (kind=RealKind) :: total_dos, dos_per_spin, total_dos_mt, dos_mt_per_spin, Ef
        real (kind=RealKind) :: spdos, cpdos
        real (kind=RealKind) :: nr(0:4)  ! assuming lmax <= 4
        real (kind=RealKind) :: ps(0:4)  ! assuming lmax <= 4
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
                 write(6,*), 'dos_ws: ', PartialDOS(id)%dos_ws(:,ia)
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
        if (getKeyValue(1,'Debye Temperature (K)',DebyeTemp, default_param=.false.) /= 0) then
           if (getNumAtomTypes() == 1)  then
              atomic_number = getLocalAtomicNumber(1,1)
              DebyeTemp = getDebyeTemperature(atomic_number)
           else
              call ErrorHandler('Gaspari-Gyorffy-Formula',           &
                               'Debye temperature is missing from input')
           endif
        endif
        Phonon_Freq2 = HALF*DebyeTemp**2
        if (MyPE == 0) then
           write(6,*), 'Debye Temp (K): ', DebyeTemp
           write(6,*), 'Phonon Frequency (K): ', sqrt(Phonon_Freq2)
        endif
        Phonon_Freq2 = Phonon_Freq2*Kelvin2Ryd*Kelvin2Ryd
!
        aGID = getGroupID("Unit Cell")
        if (isGroupExisting('E-K Plane')) then
           bGID = getGroupID('E-K Plane')
        else
           bGID = getGroupID('Energy Mesh')
        endif
!
!       Calculate total dos
        total_dos = ZERO
        total_dos_mt = ZERO
        if (getNumAtoms() > 1) then
           do id = 1, LocalNumAtoms  ! Loop atomic sites
              do ia = 1, getLocalNumSpecies(id)  ! Loop over atomic species
                 dos_ws => PartialDOS(id)%dos_ws(:,ia)
                 dos_mt => PartialDOS(id)%dos_mt(:,ia)
                 species_content = getLocalSpeciesContent(id,ia)
                 do is = 1, n_spin_pola
!                   ==================================================
!                   For now, I am using the muffin-tin DOS. It needs to 
!                   be updated.
!                   ==================================================
                    total_dos = total_dos + species_content*dos_ws(is)
                    total_dos_mt = total_dos_mt + species_content*dos_mt(is)
                 enddo
              enddo
           enddo
           call GlobalSumInGroup(aGID,total_dos)
        else
           dos_ws => PartialDOS(1)%dos_ws(:,1)
           dos_mt => PartialDOS(1)%dos_mt(:,1)
           do is = 1, n_spin_pola
!             ========================================================
!             For now, I am using the muffin-tin DOS. It needs to 
!             be updated.
!             ========================================================
              total_dos = total_dos + dos_ws(is)
              total_dos_mt = total_dos_mt + dos_mt(is)
           enddo
        endif
        dos_per_spin = HALF*total_dos  ! states/Ryd/spin
        dos_mt_per_spin = HALF*total_dos_mt  ! states/Ryd/spin
        if (MyPE == 0) then
           write(6,'(1x,a,f12.5)'), 'Atomic DOS (states/Ryd/spin): ', dos_per_spin
           write(6,'(1x,a,f12.5)'), 'Atomic DOS (states/eV/spin): ', dos_per_spin/Ryd2eV
           write(6,'(1x,a,f12.5)'), 'Muffin-tin DOS (states/Ryd/spin): ', dos_mt_per_spin
           write(6,'(1x,a,f12.5)'), 'Muffin-tin DOS (states/eV/spin): ', dos_mt_per_spin/Ryd2eV
        endif
!
        Iconst = efermi/PI**2/dos_per_spin**2 ! Ryd^3
!
!       Calculate Lamda (or EPC) ...........
        nr = ZERO
        EPC = ZERO
        do id = 1, LocalNumAtoms  ! Loop atomic sites
           kmax_phi = PartialDOS(id)%kmax_phi
           lmax_kkr = lofk(kmax_phi)
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
              SumI = ZERO
              do is = 1, n_spin_pola
                 do l=0,lmax_kkr
                    cpdos = ZERO
                    spdos = ZERO
                    do m = -l, l
                       kl = (l+1)**2-l+m
                       cpdos = cpdos + partial_dos_mt(kl,is)
                       spdos = spdos + ss_pdos_mt(kl,is)
                    enddo
!                   nr(l) = partial_dos_mt(kl,is)/ss_pdos_mt(kl,is)
                    nr(l) = cpdos/spdos
                 enddo
                 if (MyPE == 0) then
                    write(6,'(1x,2(a,f12.5))')'s-state: phase shift =',phase_shift(1,is), ', DOS ratio =',nr(0)
                    write(6,'(1x,2(a,f12.5))')'p-state: phase shift =',phase_shift(3,is), ', DOS ratio =',nr(1)
                    write(6,'(1x,2(a,f12.5))')'d-state: phase shift =',phase_shift(7,is), ', DOS ratio =',nr(2) 
                    write(6,'(1x,2(a,f12.5))')'f-state: phase shift =',phase_shift(13,is),', DOS ratio =',nr(3)
                 endif
                 do l=0,lmax_kkr-1 ! sum over l
                    klp1 = (l+2)**2-l-1
                    kl = (l+1)**2-l
                    sinTerm = sin(phase_shift(klp1,is)-phase_shift(kl,is))
                    SumI = SumI + TWO*(l+1)*sinTerm**2*nr(l+1)*nr(l)
                 enddo
              enddo
              I_Total = Iconst*SumI ! Ryd^2/au^2
!
              AtomMass=getAtomicMass(atomic_number)*MassUnit2Ryd/LightSpeed**2
              EPCnum = dos_per_spin*I_Total  ! Ryd/au^2
              EPCden = AtomMass*Phonon_Freq2 ! Ryd/au^2
              EPC = EPC + EPCnum/EPCden ! unitless
              if (getMyPEinGroup(bGID) == 0) then
                 write(6,'(/,1x,a,a)')'For species: ',getLocalAtomName(id,ia)
                 write(6,*), 'SumI: ', SumI
                 write(6,*), 'I (eV^2/A^2): ', I_Total*(Ryd2eV/Bohr2Angstrom)**2
                 write(6,*), 'AtomicMass (Ryd/c^2): ', AtomMass
                 write(6,*), 'eta (Ryd/au^2): ', EPCnum
                 write(6,*), 'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
                 write(6,*), 'M<Omega^2> (Ryd/au^2): ', EPCden
                 write(6,*), 'M<Omega^2> (eV/A^2): ', EPCden*Ryd2eV/Bohr2Angstrom**2
                 write(6,*), 'eta/(M<Omega^2>: ', EPCnum/EPCden
              endif
           enddo  ! Loop over ia
        enddo  ! Loop over id
        call GlobalSumInGroup(aGID, EPC)
        if (MyPE == 0) then
           write(6,*), 'lamda (EPC): ', EPC
        endif
!
!       Calculate mu ...
        Coulombnum = 0.26D0*total_dos/Ryd2eV
        Coulombden = ONE+total_dos/Ryd2eV
        Coulomb = Coulombnum/Coulombden
        if (MyPE == 0) then
           write(6,*), 'mu* (Coulomb pseudopotential): ', Coulomb
        endif
!
!       Calculate Tc ...
        SuperTempExp = -1.04D0*(ONE+EPC)/(EPC-Coulomb*(ONE+0.62D0*EPC))
        SuperTemp=DebyeTemp*exp(SuperTempExp)/1.45D0
        if (MyPE == 0) then
           write(6,*), 'Superconducting Transition Temp: ', SuperTemp
        endif
!
!       ==============================================================
!       Test the data published in PRB 15, 4221 (1977)
!       ==============================================================
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
        call testPRB15_4221_1977(atomic_number,ps,nr,DebyeTemp,Ef,    &
                                 dos_per_spin,EPCnum,EPC,Coulomb,SuperTemp)
!
        if (MyPE == 0) then
           write(6,'(/,a)')'For Niobium'
           write(6,'(a,f12.5)')'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
           write(6,'(a,f12.5)')'lamda = ',EPC
           write(6,'(a,f12.5)')'mu* (Coulomb pseudopotential): ', Coulomb
           write(6,'(a,f12.5)')'Superconducting Transition Temp (K): ', SuperTemp
        endif
!       ==============================================================
!       For Ti: Input from Table II and Table IV
!       --------------------------------------------------------------
        atomic_number = 22
        ps(0) = -0.641; ps(1) = -0.157; ps(2) = 0.754; ps(3) = 0.003
        nr(0) = 0.414; nr(1) = 3.101; nr(2) = 0.788; nr(3) = 6.382
        DebyeTemp = 420.0D0; Ef = 0.588; dos_per_spin = 12.38
!
        call testPRB15_4221_1977(atomic_number,ps,nr,DebyeTemp,Ef,    &
                                 dos_per_spin,EPCnum,EPC,Coulomb,SuperTemp)
!
        if (MyPE == 0) then
           write(6,'(/,a)')'For Titanium'
           write(6,'(a,f12.5)')'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
           write(6,'(a,f12.5)')'lamda = ',EPC
           write(6,'(a,f12.5)')'mu* (Coulomb pseudopotential): ', Coulomb
           write(6,'(a,f12.5)')'Superconducting Transition Temp (K): ', SuperTemp
        endif
!       ==============================================================
!       For V: Input from Table II and Table IV
!       --------------------------------------------------------------
        atomic_number = 23
        ps(0) = -0.691; ps(1) = -0.173; ps(2) = 1.030; ps(3) = 0.003
        nr(0) = 0.550; nr(1) = 2.842; nr(2) = 0.634; nr(3) = 6.935
        DebyeTemp = 380.0D0; Ef = 0.675; dos_per_spin = 12.70
!
        call testPRB15_4221_1977(atomic_number,ps,nr,DebyeTemp,Ef,    &
                                 dos_per_spin,EPCnum,EPC,Coulomb,SuperTemp)
!
        if (MyPE == 0) then
           write(6,'(/,a)')'For Vanadium'
           write(6,'(a,f12.5)')'eta (eV/A^2): ', EPCnum*Ryd2eV/Bohr2Angstrom**2
           write(6,'(a,f12.5)')'lamda = ',EPC
           write(6,'(a,f12.5)')'mu* (Coulomb pseudopotential): ', Coulomb
           write(6,'(a,f12.5)')'Superconducting Transition Temp (K): ', SuperTemp
        endif
!       ==============================================================
!
        end subroutine gaspari_gyorffy_formula
!   ==================================================================
!
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine testPRB15_4221_1977(anum,ps,nr,thetaD,Ef,dos_per_spin, &
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
    Tc=thetaD*exp(SuperTempExp)/1.45D0
!
    end subroutine testPRB15_4221_1977
