!  *******************************************************************
!  Sunroutine for calculating the superconducting transition temperature
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
!                    quantum number l, and mofk(kl) gives magnetic quantum number m
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine gaspari_gyorffy_formula(kmax_phi_max,kmax_phi,n_spin_pola, &
                                      atomic_number,                     &
                                      ss_dos_mt,dos_mt,phase_shift,      &
                                      ss_pdos_mt,partial_dos_mt,iprint)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   use IntegerFactorsModule, only : lofk, mofk
   use ChemElementModule, only : getDebyeTemperature
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax_phi_max, kmax_phi, n_spin_pola
   integer (kind=IntKind), intent(in) :: atomic_number, iprint
!
   real (kind=RealKind), intent(in) :: ss_dos_mt(n_spin_pola)
   real (kind=RealKind), intent(in) :: dos_mt(n_spin_pola)
   real (kind=RealKind), intent(in) :: phase_shift(kmax_phi_max,n_spin_pola)
   real (kind=RealKind), intent(in) :: ss_pdos_mt(kmax_phi_max,n_spin_pola)
   real (kind=RealKind), intent(in) :: partial_dos_mt(kmax_phi_max,n_spin_pola)
!
!
   end subroutine gaspari_gyorffy_formula
