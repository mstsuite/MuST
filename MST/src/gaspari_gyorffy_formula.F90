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
        use SystemModule, only : getAtomicNumber, getNumAlloyElements, getAlloyElementContent
        use SystemModule, only : getAtomName
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
        real (kind=RealKind), pointer :: ss_pDOS_aver(:,:)
        real (kind=RealKind), pointer :: pDOS_aver(:,:)
!
        real (kind=RealKind) :: I_Total,sinTerm,Inum,Iden
        real (kind=RealKind), allocatable :: PhononFreq2(:), AtomMass(:)
        real (kind=RealKind) :: EPC, EPC_pc, EPCnum(2), EPCden
        real (kind=RealKind) :: Coulomb,Coulombnum,Coulombden,mu_star
        real (kind=RealKind) :: SuperTemp,SuperTempExp
        real (kind=RealKind) :: SumI(2),Iconst
!
        real (kind=RealKind) :: DebyeTemp, M_aver, TD
        real (kind=RealKind) :: species_content
        real (kind=RealKind) :: total_dos, dos_per_spin, total_dos_mt, dos_mt_per_spin, Ef
        real (kind=RealKind) :: spdos, cpdos,pps, cfac, val_au, av_phonon_freq, eta, sfac
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
        allocate(PhononFreq2(NumAtoms), AtomMass(NumAtoms))
!
        if (n_spin_pola == 1) then
           cfac = HALF
           sfac = TWO
        else
           cfac = ONE
           sfac = ONE
        endif
!
        if (MyPE == 0) then
           write(6,*), ' '
           write(6,*), ' '
           write(6,*), "***************************************************"
           write(6,*),         'Output of Gaspari-Gyorffy Formula'
           write(6,*), 'Calculate the Superconducting Transition Teperature'
           write(6,*), "***************************************************"
           write(6,'(1x,a,t30,a,i3)')'n_spin_pola',':',n_spin_pola
           do id = 1, LocalNumAtoms
              write(6,'(1x,a,t30,a,i3)')'id',':',id
              write(6,'(1x,a,t30,a,i3)')'kmax_phi',':',PartialDOS(id)%kmax_phi
              write(6,'(1x,a,t30,a,i3)')'Number of species on site',':',getLocalNumSpecies(id)
              do ia = 1, getLocalNumSpecies(id)
                 if (getLocalNumSpecies(id) > 1) then
                    write(6,'(1x,a,t30,a,i3,a)'), 'ia',':',ia,' ------------'
                 endif
                 write(6,'(1x,a,t30,a,i3)')'Atomic number',':',getLocalAtomicNumber(id,ia)
                 if (n_spin_pola == 1) then
                    write(6,'(1x,a,t30,a,d15.8,a)')'dos_ws per spin',':',PartialDOS(id)%dos_ws(1,ia)*HALF,' (St/Ryd)'
                    write(6,'(1x,a,t30,a,d15.8,a)')'dos_ws per spin',':',PartialDOS(id)%dos_ws(1,ia)*HALF/Ryd2eV,' (St/eV)'
                 else
                    write(6,'(1x,a,t30,a,d15.8,a)')'dos_ws for spin-up',':',PartialDOS(id)%dos_ws(1,ia),' (St/Ryd)'
                    write(6,'(1x,a,t30,a,d15.8,a)')'dos_ws for spin-up',':',PartialDOS(id)%dos_ws(1,ia)/Ryd2eV,' (St/eV)'
                    write(6,'(1x,a,t30,a,d15.8,a)')'dos_ws for spin-down',':',PartialDOS(id)%dos_ws(2,ia),' (St/Ryd)'
                    write(6,'(1x,a,t30,a,d15.8,a)')'dos_ws for spin-down',':',PartialDOS(id)%dos_ws(2,ia)/Ryd2eV,' (St/eV)'
                 endif
              enddo
           enddo
           write(6,'(/,1x,a,t30,a,f12.8)')'Efermi',':',efermi
           !do kl=1,kmax_phi
           !   write(6,*), 'kl: ', kl, ' l: ', lofk(kl), ' m: ', mofk(kl)
           !   write(6,*), 'phaseshift: ', phase_shift(kl,n_spin_pola)
           !   write(6,*), 'dos: ', partial_dos_mt(kl,n_spin_pola)
           !   write(6,*), 'ss_dos: ', ss_pdos_mt(kl,n_spin_pola)
           !enddo
        endif
!
!       Determine the atomic mass
!       ==============================================================
        do ig = 1, NumAtoms
           TD = ZERO
           M_aver = ZERO
           do ia = 1,  getNumAlloyElements(ig)  ! Loop over atomic species
              atomic_number = getAtomicNumber(ig,ia)
              species_content = getAlloyElementContent(ig,ia)
              M_aver = M_aver + species_content*getAtomicMass(atomic_number)
           enddo
           AtomMass(ig) = M_aver*MassUnit2Ryd/LightSpeed**2 ! In atomic units
        enddo
!
!       Determine the Debye Temperature
!       ==============================================================
        if (getKeyValue(1,'Average of phonon frequency (1/sec)',av_phonon_freq,default_param=.false.) == 0) then
           DebyeTemp = av_phonon_freq/Boltzmann
!          write(6,*), 'Debye 1'
        else if (getKeyValue(1,'Average of phonon frequency (K)',av_phonon_freq,default_param=.false.) == 0) then
           DebyeTemp = av_phonon_freq
!          write(6,*), 'Debye 2'
        else if (GetKeyValue(1,'Debye Temperature (K)',DebyeTemp, default_param=.false.) /= 0) then
!          ===========================================================
!          In case the Debye temperature is not provided from the input
!          ===========================================================
           if (NumAtoms == 1) then  ! For a simple metal or a random
                                    ! alloy on a simple lattice case
!             write(6,*), 'Debye 3'
!             ========================================================
!             The Debye temperature is to be taken from the internal data base.
!             In particular, if system is an alloy, it will be averaged over the 
!             constitutional species, The reliability of this approach
!             to the alloy's Debye temperature is questionable, though
!             ========================================================
              DebyeTemp = ZERO
              do ia = 1, getLocalNumSpecies(1)  ! Loop over atomic species
                 atomic_number = getLocalAtomicNumber(1,ia)
                 species_content = getLocalSpeciesContent(1,ia)
                 TD = getDebyeTemperature(atomic_number)
                 if (TD < ONE) then
                    write(6,'(3a)')'WARNING:: The Debye temp for species ',getName(atomic_number), &
                                   ' is not available in database.'
                 else
                    DebyeTemp = DebyeTemp + species_content*TD
                 endif
              enddo
           else ! For the multi-sublattice case
!             ========================================================
!             The Debye temperature needs to be given from the input
!             ========================================================
              call ErrorHandler('Gaspari-Gyorffy-Formula',            &
                                'In supe cell case, Debye Temperature (K) is required')
           endif
        endif
!
!       Determine the average phonon frequency squared
!       ==============================================================
        ind_array = 0; val_array = ZERO; PhononFreq2 = ZERO
        if (getKeyLabelIndexValue(1,'Atomic mass times <omega^2> (eV/Anst^2)', &
                                  NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
           do ig = 1, NumAtoms
              if (MyPE == 0) then
                 atomic_number = getAtomicNumber(ig)
                 write(6,*) getName(atomic_number), ': M*<omega^2> = ',val_array(ig),'(eV/Anst^2)'
              endif
              val_au = val_array(ig)*Bohr2Angstrom**2/Ryd2eV ! This is M*<omega^2> in atomic units
              PhononFreq2(ig) = val_au/AtomMass(ig)
           enddo
        else if (getKeyLabelIndexValue(1,'Atomic mass times <omega^2> (Ryd/BohrRad^2)', &
                                       NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
           do ig = 1, NumAtoms
              if (MyPE == 0) then
                 atomic_number = getAtomicNumber(ig)
                 write(6,*) getName(atomic_number), ': M*<omega^2> = ',val_array(ig),'(Ryd/BohrRad^2)'
              endif
              PhononFreq2(ig) = val_array(ig)/AtomMass(ig)
           enddo
        else if (getKeyLabelIndexValue(1,'Average of phonon frequency squared (1/sec^2)', &
                                       NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
           do ig = 1, NumAtoms
              if (MyPE == 0) then
                 atomic_number = getAtomicNumber(ig)
                 write(6,*) getName(atomic_number), ': <omega^2> = ',val_array(ig),'(1/sec^2)'
              endif
              PhononFreq2(ig) = val_array(ig)
           enddo
        else if (getKeyLabelIndexValue(1,'Average of phonon frequency squared (K^2)', &
                                       NumAtomTypes,AtomTypeName,NumAtoms,AtomType,ind_array,val_array) == 0) then
           do ig = 1, NumAtoms
              if (MyPE == 0) then
                 atomic_number = getAtomicNumber(ig)
                 write(6,*) getName(atomic_number), ': <omega^2> = ',val_array(ig),'(K^2)'
              endif
              PhononFreq2(ig) = val_array(ig)*Kelvin2Ryd**2
           enddo
        else
!          ===========================================================
!          In case where the average of phonon frequency squared is not
!          given from the input, we take the Debye temperature from the 
!          internal data base and use it to calculate the average of
!          the phonon frequency squared.
!          ===========================================================
           if (NumAtoms == 1) then
              PhononFreq2(1) = HALF*(DebyeTemp*Kelvin2Ryd)**2
           else
              do ig = 1, NumAtoms
                 TD = ZERO
                 do ia = 1,  getNumAlloyElements(ig)  ! Loop over atomic species
                    atomic_number = getAtomicNumber(ig,ia)
                    species_content = getAlloyElementContent(ig,ia)
                    TD = TD + species_content*getDebyeTemperature(atomic_number)
                 enddo
                 PhononFreq2(ig) = HALF*(TD*Kelvin2Ryd)**2
              enddo
           endif
        endif
        do ig = 1, NumAtoms
           if (PhononFreq2(ig) < TEN2m6) then
              call ErrorHandler('Gaspari-Gyorffy-Formula','Invalid value of <omega^2>',PhononFreq2(ig))
           endif
        enddo
!
        deallocate(ind_array, val_array)
!
        if (MyPE == 0) then
           write(6,'(/,1x,a,t30,a,f8.3,/)')'Debye Temp (K)',':',DebyeTemp
           do ig = 1, NumAtoms  ! Loop atomic sites
              if (getNumAlloyElements(ig) == 1) then
                 write(6,'(1x,a3,a,d15.8,a)')getAtomName(ig,1),': Average of Phonon Frequency Square =', &
                                             PhononFreq2(ig),' (Ryd^2)'
              else
                 write(6,'(1x,a3,a,d15.8,a)')'CPA',': Average of Phonon Frequency Square =', &
                                             PhononFreq2(ig),' (Ryd^2)'
              endif
              write(6,'(41x,a,d15.8,a)')'=',PhononFreq2(ig)/Kelvin2Ryd**2,' (K^2)'
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

        Iconst = efermi/PI**2/dos_per_spin**2 ! Ryd^3

        if (MyPE == 0) then
           write(6,'(/,1x,a,f10.5,a)')'WS-Volume DOS of Unit cell per spin =', dos_per_spin,' (states/Ryd/spin))'
           write(6,'(  1x,a,f10.5,a)')'                                    =', dos_per_spin/Ryd2eV,' (states/eV/spin)'
           write(6,'(  1x,a,f10.5,a)')'MT-Volume DOS of Unit cell per spin =', dos_mt_per_spin,' (states/Ryd/spin))'
           write(6,'(  1x,a,f10.5,a)')'                                    =', dos_mt_per_spin/Ryd2eV,' (states/eV/spin)'
        endif
!
        kmax_kkr_max = 0
        do id = 1, LocalNumAtoms
           kmax_kkr_max = max(kmax_kkr_max,PartialDOS(id)%kmax_kkr)
        enddo

        kappa = sqrt(efermi)

!       Calculate Lamda (or EPC) ...........
!       ==============================================================
        nr = ZERO
        EPC = ZERO
        do id = 1, LocalNumAtoms  ! Loop atomic sites
           kmax_phi = PartialDOS(id)%kmax_phi
           kmax_kkr = PartialDOS(id)%kmax_kkr
           lmax_kkr = lofk(kmax_phi)

           ig = getGlobalIndex(id) 
!
           EPCden = AtomMass(ig)*PhononFreq2(ig) ! Ryd/au^2

           do ia = 1, getLocalNumSpecies(id)  ! Loop over atomic species
                                              ! at each atomic site. The
                                              ! number of species is usually 1,
                                              ! except for the KKR-CPA case.
              species_content = getLocalSpeciesContent(id,ia)
              phase_shift => PartialDOS(id)%phase_shift(:,:,ia)
              ss_pdos_mt => PartialDOS(id)%ss_pdos_mt(:,:,ia)
              partial_dos_mt => PartialDOS(id)%partial_dos_mt(:,:,ia)
!
              if (getMyPEinGroup(bGID) == 0) then
                 write(6,'(/,1x,a,a)')'For species: ',getLocalAtomName(id,ia)
              endif
!
              SumI = ZERO
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
                    SumI(is) = SumI(is) + sfac*(l+1)*sinTerm**2*nr(l+1)*nr(l)
                 enddo
              enddo ! loop over is
!
              eta = ZERO
              do is = 1, n_spin_pola
                 EPCnum(is) = dos_per_spin*Iconst*SumI(is)
                 eta = eta + EPCnum(is)  ! Ryd/au^2
              enddo
!
              EPC = EPC + getLocalSpeciesContent(id,ia)*eta/EPCden ! unitless
!
              if (getMyPEinGroup(bGID) == 0) then
                 write(6,'(/)')
                 write(6,'(1x,a,f12.5)')'AtomicMass (Ryd/c^2)    =', AtomMass(ig)
                 write(6,'(1x,a,f12.5)')'M<Omega^2> (Ryd/au^2)   =', EPCden
                 write(6,'(1x,a,f12.5)')'M<Omega^2> (eV/A^2)     =', EPCden*Ryd2eV/Bohr2Angstrom**2
                 do is = 1, n_spin_pola
                    if (n_spin_pola == 1) then
                       write(6,'(/)')
                    else
                       write(6,'(/,1x,a,i2)')'Spin index              =', is
                    endif
                    write(6,'(1x,a,f12.5)')'SumI                    =', SumI(is)
                    write(6,'(1x,a,f12.5)')'I (eV^2/A^2)            =', Iconst*SumI(is)*(Ryd2eV/Bohr2Angstrom)**2
                    write(6,'(1x,a,f12.5)')'I (Ryd^2/au^2)          =', Iconst*SumI(is)
                    write(6,'(1x,a,f12.5)')'EPCnum(is) (Ryd/au^2)   =', EPCnum(is)
                    write(6,'(1x,a,f12.5)')'EPCnum(is) (eV/A^2)     =', EPCnum(is)*Ryd2eV/Bohr2Angstrom**2
                    write(6,'(1x,a,f12.5)')'EPCnum(is)/(M<Omega^2>  =', EPCnum(is)/EPCden
                 enddo
                 write(6,'(1x,a,f12.5)')'eta (Ryd/au^2)          =', eta
                 write(6,'(1x,a,f12.5)')'eta (eV/A^2)            =', eta*Ryd2eV/Bohr2Angstrom**2
                 write(6,'(1x,a,f12.5)')'lamda = eta/(M<Omega^2> =', eta/EPCden
              endif
           enddo  ! Loop over ia
        enddo  ! Loop over id

        call GlobalSumInGroup(aGID, EPC)

        if (MyPE == 0) then
           write(6,'(/,1x,a)')'For the system ...'
           write(6,'(1x,a,t32,a,f12.5)')'Total lamda (EPC)/unit cell','=',EPC
        endif
!
!       ==============================================================
!       In the KKR-CPA case, we also consider a different approach to
!       calculation of averaged EPC (or lambda). In this approach,
!       which was first suggested by P. Chatterje, t_cpa is used to 
!       calculate the averaged partial phase shift.
!       ==============================================================
        if (isKKRCPA()) then
           if (MyPE == 0) then
              write(6,'(/,1x,a)')'An alternative approach to the AVERAGED ETA AND LAMBDA'
           endif
           allocate(phase_shift_aver(kmax_kkr_max,n_spin_pola))
           allocate(ss_pDOS_aver(kmax_kkr_max,n_spin_pola))
           allocate(pDOS_aver(kmax_kkr_max,n_spin_pola))

           nr = ZERO
           EPC_pc = ZERO
           do id = 1, LocalNumAtoms  ! Loop atomic sites
              kmax_phi = PartialDOS(id)%kmax_phi
              kmax_kkr = PartialDOS(id)%kmax_kkr
              lmax_kkr = lofk(kmax_phi)

              ig = getGlobalIndex(id)
!
              EPCden = AtomMass(ig)*PhononFreq2(ig) ! Ryd/au^2

              pDOS_aver = ZERO
              ss_pDOS_aver = ZERO
              do ia = 1, getLocalNumSpecies(id)
                 partial_dos_mt => PartialDOS(id)%partial_dos_mt(:,:,ia)
                 ss_pdos_mt => PartialDOS(id)%ss_pdos_mt(:,:,ia)
                 species_content = getLocalSpeciesContent(id,ia)
                 do is = 1, n_spin_pola
                    do kl = 1, PartialDOS(id)%kmax_kkr
                       pDOS_aver(kl,is) = pDOS_aver(kl,is) + species_content*partial_dos_mt(kl,is)
                       ss_pDOS_aver(kl,is) = ss_pDOS_aver(kl,is) + species_content*ss_pdos_mt(kl,is)
                    enddo
                 enddo
              enddo

              do is = 1, n_spin_pola
                 t_mat => PartialDOS(id)%t_aver(:,:,is)
                 do kl = 1, PartialDOS(id)%kmax_kkr
                    phase_shift_aver(kl,is) = atan(kappa*t_mat(kl,kl)/(-CONE+SQRTm1*kappa*t_mat(kl,kl)))
                 enddo
              enddo

              SumI = ZERO
              do is = 1, n_spin_pola
                 do l=0,lmax_kkr
                    cpdos = ZERO
                    spdos = ZERO
                    do m = -l, l
                       kl = (l+1)**2-l+m
                       cpdos = cpdos + pDOS_aver(kl,is)
                       spdos = spdos + ss_pDOS_aver(kl,is)
                    enddo
                    nr(l) = cpdos/spdos
                 enddo
                 do l=0,lmax_kkr-1 ! sum over l
                    klp1 = (l+2)**2-l-1
                    kl = (l+1)**2-l
                    sinTerm = abs(sin(phase_shift_aver(klp1,is)-phase_shift_aver(kl,is)))
                    SumI(is) = SumI(is) + sfac*(l+1)*sinTerm**2*nr(l+1)*nr(l)
                 enddo
              enddo
!
              eta = ZERO
              do is = 1, n_spin_pola
                 EPCnum(is) = dos_per_spin*Iconst*SumI(is)
                 eta = eta + EPCnum(is)  ! Ryd/au^2
              enddo
!
              EPC_pc = EPC_pc + eta/EPCden ! unitless
!
              if (getMyPEinGroup(bGID) == 0) then
                 do is = 1, n_spin_pola
                    if (n_spin_pola == 1) then
                       write(6,'(/)')
                    else
                       write(6,'(/,1x,a,i2)')'Spin index              =', is
                    endif
                    write(6,'(1x,a,f12.5)')'SumI                    =', SumI(is)
                    write(6,'(1x,a,f12.5)')'I (eV^2/A^2)            =', Iconst*SumI(is)*(Ryd2eV/Bohr2Angstrom)**2
                    write(6,'(1x,a,f12.5)')'I (Ryd^2/au^2)          =', Iconst*SumI(is)
                    write(6,'(1x,a,f12.5)')'EPCnum(is) (Ryd/au^2)   =', EPCnum(is)
                    write(6,'(1x,a,f12.5)')'EPCnum(is) (eV/A^2)     =', EPCnum(is)*Ryd2eV/Bohr2Angstrom**2
                    write(6,'(1x,a,f12.5)')'EPCnum(is)/(M<Omega^2>  =', EPCnum(is)/EPCden
                 enddo
                 write(6,'(1x,a,f12.5)')'eta (Ryd/au^2)          =', eta
                 write(6,'(1x,a,f12.5)')'eta (eV/A^2)            =', eta*Ryd2eV/Bohr2Angstrom**2
                 write(6,'(1x,a,f12.5)')'lamda = eta/(M<Omega^2> =', eta/EPCden
              endif
           enddo  ! Loop over id

           call GlobalSumInGroup(aGID, EPC_pc)

           if (MyPE == 0) then
              write(6,'(/,1x,a)')'For the system ...'
              write(6,'(1x,a,t32,a,f12.5)')'Total lamda (EPC)/unit cell','=', EPC_pc
           endif
           deallocate(phase_shift_aver, ss_pDOS_aver, pDOS_aver)
        endif

!       Calculate mu ...
!       ==============================================================
        Coulombnum = 0.26D0*total_dos/Ryd2eV
        Coulombden = ONE+total_dos/Ryd2eV
        Coulomb = Coulombnum/Coulombden
        if (getKeyValue(1,'mu* (e-e interaction constant)',mu_star, default_param=.false.) == 0) then
           Coulomb = mu_star
           if (MyPE == 0) then
              write(6,'(1x,a)')'Using the input mu* (Coulomb pseudopotential) ......'
           endif
        endif
        if (MyPE == 0) then
           write(6,'(/,1x,a,t32,a,f12.5)')'mu* (Coulomb pseudopotential)','=', Coulomb
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
              write(6,'(1x,a,t32,a,f12.5)')'EPC-Coulomb*(ONE+0.62D0*EPC)','=', EPC-Coulomb*(ONE+0.62D0*EPC)
           endif
           SuperTempExp = -1.04D0*(ONE+EPC)/(EPC-Coulomb*(ONE+0.62D0*EPC))
           SuperTemp=DebyeTemp*exp(SuperTempExp)/1.45D0
        else
           if (MyPE == 0) then
              write(6,'(1x,a)') 'EPC <= Coulomb'
           endif
           SuperTemp=ZERO
        endif
!
        if (MyPE == 0) then
           write(6,'(/,1x,a)')'************************************************'
           write(6,'(1x,a,f12.5)')'Superconducting Transition Temp (K):',SuperTemp
           write(6,'(1x,a,/)')'************************************************'
        endif
!
        if (isKKRCPA()) then
           if (EPC_pc > Coulomb) then
              SuperTempExp = -1.04D0*(ONE+EPC_pc)/(EPC_pc-Coulomb*(ONE+0.62D0*EPC_pc))
              SuperTemp=DebyeTemp*exp(SuperTempExp)/1.45D0
           else
              SuperTemp=ZERO
           endif
           if (MyPE == 0) then
              write(6,'(/,1x,a)')'With an alternative approach to the GG formula applied to random alloys ...'
              write(6,'(/,1x,a)')'************************************************'
              write(6,'(1x,a,f12.5)')'Superconducting Transition Temp (K):',SuperTemp
              write(6,'(1x,a,/)')'************************************************'
           endif
        endif
!
        deallocate(PhononFreq2, AtomMass)
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
    real (kind=RealKind) :: PhononFreq2, Iconst, SumI, sinTerm, I_Total
    real (kind=RealKind) :: AtomMass, EPCden, total_dos
    real (kind=RealKind) :: Coulombnum, Coulombden, SuperTempExp
!
    integer (kind=IntKind), intent(in) :: anum
    integer (kind=IntKind) :: l
!
    PhononFreq2 = HALF*(thetaD*Kelvin2Ryd)**2
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
    EPCden = AtomMass*PhononFreq2 ! Ryd/au^2
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
