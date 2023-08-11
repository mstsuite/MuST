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
        use PublicTypeDefinitionsModule, only : PDOSStruct
        use IntegerFactorsModule, only : lofk, mofk
        use InputModule, only : getKeyValue
        use ChemElementModule, only : getDebyeTemperature, getAtomicMass
        use MPPModule, only : myPE
        use AtomModule, only : getLocalNumSpecies, getLocalAtomicNumber,      &
                               getLocalSpeciesContent
        implicit none
!
        integer (kind=IntKind), intent(in) :: LocalNumAtoms
        integer (kind=IntKind), intent(in) :: n_spin_pola
        integer (kind=IntKind), intent(in) :: iprint
        integer (kind=IntKind) :: atomic_number
        integer (kind=IntKind) :: id, ia, l, kl, kmax_phi
!
        real (kind=RealKind), intent(in) :: efermi
!
        type (PDOSStruct), intent(in) :: PartialDOS(LocalNumAtoms)
!
        real (kind=RealKind), pointer :: ss_dos_mt(:)
        real (kind=RealKind), pointer :: dos_mt(:)
        real (kind=RealKind), pointer :: phase_shift(:,:)
        real (kind=RealKind), pointer :: ss_pdos_mt(:,:)
        real (kind=RealKind), pointer :: partial_dos_mt(:,:)
!
        real (kind=RealKind) :: I_Total,sinTerm,Inum,Iden
        real (kind=RealKind) :: Phonon_Freq,kB,hbar
        real (kind=RealKind) :: EPC, EPCnum,EPCden
        real (kind=RealKind) :: AtomMass
        real (kind=RealKind) :: Coulomb,Coulombnum,Coulombden
        real (kind=RealKind) :: SuperTemp,SuperTempExp
        real (kind=RealKind) :: SumI,Iconst1,Iconst2,Pi
        real (kind=RealKind) :: DebyeTemp
        real (kind=RealKind) :: species_content
!
        if (myPE == 0 .and. iprint >= 0) then
        write(6,*), ' '
        write(6,*), ' '
        write(6,*), "***************************************************"
        write(6,*),         'Output of Gaspari-Gyorffy Formula'
        write(6,*), 'Calculate the Superconducting Transition Teperature'
        write(6,*), "***************************************************"
        write(6,*), 'n_spin_pola: ', n_spin_pola
        do id = 1, LocalNumAtoms
           write(6,*), 'kmax_phi: ', PartialDOS(id)%kmax_phi
           do ia = 1, getLocalNumSpecies(id)
        write(6,*), 'atomic number: ', getLocalAtomicNumber(id,ia)
        write(6,*), 'ss_dos_mt: ', PartialDOS(id)%ss_dos_mt(:,ia)
        write(6,*), 'dos_mt: ', PartialDOS(id)%dos_mt(:,ia)
           enddo
        enddo
        write(6,*), 'Efermi: ', efermi
        !do kl=1,kmax_phi
        !write(6,*), 'kl: ', kl, ' l: ', lofk(kl), ' m: ', mofk(kl)
        !write(6,*), 'phaseshift: ', phase_shift(kl,n_spin_pola)
        !write(6,*), 'dos: ', partial_dos_mt(kl,n_spin_pola)
        !write(6,*), 'ss_dos: ', ss_pdos_mt(kl,n_spin_pola)
        !enddo
        write(6,*), ' '
        endif
!
        SumI=0
        open(unit=1,file="ss_pdos_mt")
        open(unit=2,file="partial_dos_mt")
        open(unit=3,file="phase_shift")
!
        do id = 1, LocalNumAtoms  ! Loop atomic sites
           kmax_phi = PartialDOS(id)%kmax_phi
           do ia = 1, getLocalNumSpecies(id)  ! Loop over atomic species
                                              ! at each atomic site. The
                                              ! number of species is usually 1,
                                              ! except for the KKR-CPA case.
              species_content = getLocalSpeciesContent(id,ia)
              atomic_number = getLocalAtomicNumber(id,ia)
              ss_dos_mt => PartialDOS(id)%ss_dos_mt(:,ia)
              dos_mt => PartialDOS(id)%dos_mt(:,ia)
              phase_shift => PartialDOS(id)%phase_shift(:,:,ia)
              ss_pdos_mt => PartialDOS(id)%ss_pdos_mt(:,:,ia)
              partial_dos_mt => PartialDOS(id)%partial_dos_mt(:,:,ia)

        do l=0,lofk(kmax_phi)-1 ! sum over l
            sinTerm = sin(phase_shift((l+2)**2,n_spin_pola)-phase_shift((l+1)**2,n_spin_pola))
            Inum = (2*(l+1)*sinTerm**2*partial_dos_mt((l+2)**2-l-1,n_spin_pola)* &
                                partial_dos_mt((l+1)**2-l,n_spin_pola))
            Iden = ss_pdos_mt((l+2)**2,n_spin_pola)*ss_pdos_mt((l+1)**2,n_spin_pola)
            SumI = SumI + Inum/Iden
        enddo
        Pi=3.1415926535893D0
        kB= 8.617E-5 ! eV/K
        hbar=6.582E-16 ! eVs
        Iconst1 = (2*0.511*10**6)/3892729 ! 2m_e/hbar^2 => 1/(eV*A^2)
        Iconst2 = (efermi*13.606D0**3)/(Pi**2*(dos_mt(n_spin_pola))**2) ! eV^3
        I_Total = Iconst1*Iconst2*SumI ! eV^2/A^2
        write(6,*), 'I (eV^2/A^2): ', I_Total
        write(6,*), 'I (Ry^2/au^2): ', I_Total*(0.073)**2*(0.529)**2
        ! Calculate lambda
        AtomMass=getAtomicMass(atomic_number)
        if (getKeyValue(1,'Debye Temperature',DebyeTemp) /= 0) then
            DebyeTemp = getDebyeTemperature(atomic_number)
        endif
        Phonon_Freq = (1D0/2D0)*(kB*DebyeTemp)**2/(hbar**2) ! 1/s^2
        EPCnum = dos_mt(n_spin_pola)/13.606D0*I_Total ! eV/A^2
        write(6,*), 'eta (eV/A^2): ', EPCnum
        EPCden = (AtomMass*1.66E-27*6.25E18)*(Phonon_Freq/1E20) ! eV/A^2
        EPC = EPCnum/EPCden ! unitless
        write(6,*), 'EPC: ', EPC
        ! Calculate mu
        Coulombnum = 0.13D0*dos_mt(n_spin_pola)/13.606D0
        Coulombden = 1D0+dos_mt(n_spin_pola)/13.606D0
        Coulomb = Coulombnum/Coulombden
        write(6,*), 'Coulomb: ', Coulomb
        write(6,*), 'AtomicMass: ', AtomMass
        write(6,*), 'Debye Temp: ', DebyeTemp
        write(6,*), 'Phonon Frequency: ', Phonon_Freq
        SuperTempExp = -1.04*(1D0+EPC)/(EPC-Coulomb*(1D0+0.62D0*EPC))
        SuperTemp=(DebyeTemp/1.45D0)*exp(SuperTempExp)

           enddo  ! Loop over ia
        enddo  ! Loop over id
        write(6,*), 'Superconducting Transition Temp: ', SuperTemp
        end subroutine gaspari_gyorffy_formula
