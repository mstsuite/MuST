!  =============================================================================
!  Table of Input parameters and their default value
!  =============================================================================
   n = n + 1; Keys(n) = 'Current File Path';                     Values(n) = './'
!  n = n + 1; Keys(n) = 'Current File Name';                     Values(n) = 'i_lsms_D'
!  n = n + 1; Keys(n) = 'Info Table File Name';                  Values(n) = 'info_table_D'
   n = n + 1; Keys(n) = 'No. Iterations (> 0)';                  Values(n) = '60'
   n = n + 1; Keys(n) = 'Method of SCF Calculation';             Values(n) = '2'
!                     = -2: Single Site; = -1: ScreenKKR-LSMS; = 0: Screen-KKR; = 1: LSMS; = 2: KKR; = 3: KKR-CPA
   n = n + 1; Keys(n) = 'Output to Screen (y/n)';                Values(n) = 'n'
   n = n + 1; Keys(n) = 'Output Level (>= -1)';                  Values(n) = '0'
   n = n + 1; Keys(n) = 'Output Proc. ID (>= -1)';               Values(n) = '0'
!                     = -1: All processors write data to individual output file
!                     >= 0: the logical indexes (separated by space or ",") of the processors that write to unit 6
   n = n + 1; Keys(n) = 'Output Atom ID (>= 0)';                 Values(n) = '0' ! This format is to allow back-compatable 
   n = n + 1; Keys(n) = 'Output Atom ID (>= -1)';                Values(n) = '0'
!                     = -1: The results for each atom in the unit cell are written to the output
!                     =  0: The results for the atoms mapped on the output procesor are written to the output
!                     >= 1: the global indexes (separated by space and/or ",") of the atoms whose info are written 
!                           to the output
   n = n + 1; keys(n) = 'Generate System Movie';                 Values(n) = '1'
!                     = 0 - No system movie output(Default)
!                     = 1 - generate movie
   n = n + 1; Keys(n) = 'Stop-at Routine Name';                  Values(n) = 'main'
   n = n + 1; Keys(n) = 'No. Iter for Each Pot. Write';          Values(n) = '5'
   n = n + 1; Keys(n) = 'No. Iter for Each Movie';               Values(n) = '0'
!                     = 0 - No system movie output(Default)
!                    >= 1 - A system movie is written at the end of every this number of SCF iterations
   n = n + 1; Keys(n) = 'Calc. Harris Energy (H.E.)';            Values(n) = '0'
!                     = 0: Do Not Calc. H.E.; = 1: H.E. at Updated mu; = 2: H.E. at Fixed mu
   n = n + 1; Keys(n) = 'No. Gauss Pts. along r';                Values(n) = '80'
   n = n + 1; Keys(n) = 'No. Gauss Pts. along theta';            Values(n) = '60'
   n = n + 1; Keys(n) = 'Valence Band Bottom Est.';              Values(n) = '-0.5'
   n = n + 1; Keys(n) = 'Temperature Parameter (K)';             Values(n) = '0.000'
   n = n + 1; Keys(n) = 'DOS Run ID';                            Values(n) = '0'
!                     =  0: DOS along the real energy axis will not be calculated;
!                     >  0: An integer representing the global index of an atom, to indicate the DOS to be calculated
!                           for the atom;
!                     = -1: DOS will be calculated for all the atoms.
   n = n + 1; Keys(n) = 'Uniform Grid Parameters';               Values(n) = '64  64  64'
!                     > 0: Three intergers used to define the grid numbers along three Bravais lattice vector directions
!             Note: Uniform grid is used to generate the new potential. These numbers must be an interger power of 2,
!                   and will affect the accuracy of the Poisson solver.
   n = n + 1; Keys(n) = 'Visual Grid Type (0<D<4)';              Values(n) = '3'
!                     = 1: Define a line (1D); = 2: Define a plane (2D); = 3: Define a rectilinear region (3D)
   n = n + 1; Keys(n) = 'Visual Grid Scale';                     Values(n) = '1.00'
   n = n + 1; Keys(n) = 'Visual Grid Origin Vector';             Values(n) = '0.0 0.0 0.0'
   n = n + 1; Keys(n) = 'Visual Grid Vector 1';                  Values(n) = '10.00  0.00  0.00'
   n = n + 1; Keys(n) = 'Visual Grid Vector 2';                  Values(n) = '0.00 10.00  0.00'
   n = n + 1; Keys(n) = 'Visual Grid Vector 3';                  Values(n) = '0.00  0.00 10.00'
!             Note: Vectors defining the grid boundaries
   n = n + 1; Keys(n) = 'Visual Grid Points';                    Values(n) = '10 10 10'
!             Note: Visual grid will use only first D parameters for Grid Vectors and Grid Points
   n = n + 1; Keys(n) = 'Visual Line Vector';                    Values(n) = '1.00  0.00  0.00'
   n = n + 1; Keys(n) = 'Visual Line Points';                    Values(n) = '10'
   n = n + 1; Keys(n) = 'Output Electron Density ID (>= -1)';    Values(n) = '-1'
!                     = -1: do no print; = 0:  print brief data; = 1: print data on visual grid
   n = n + 1; Keys(n) = 'Output Density Format';                 Values(n) = '2'
!                     = 0: x y z rho format; 
!                     = 1: Legacy ".vtk"  format for ParaView visualization
!                     = 2: ".xsf" format for VESTA
   n = n + 1; Keys(n) = 'Energy (Ryd) Tol (> 0)';                Values(n) = '0.000001'
   n = n + 1; Keys(n) = 'Potential Tol (> 0)';                   Values(n) = '0.0000001'
   n = n + 1; Keys(n) = 'Fermi Energy Tol (> 0)';                Values(n) = '0.0000001'
   n = n + 1; Keys(n) = 'SuperLU Tol (> 0)';                     Values(n) = '0.0000001'
   n = n + 1; Keys(n) = 'K-space Check Tol (> 0)';               Values(n) = '0.0000001'
   n = n + 1; Keys(n) = 'Other RMS Tol (> 0)';                   Values(n) = '0.0000001'
   n = n + 1; Keys(n) = 'No. Atoms in System (> 0)';             Values(n) = '1'
   n = n + 1; Keys(n) = 'Atomic Position File Name';             Values(n) = 'position.dat'
   n = n + 1; Keys(n) = 'Text Identification';                   Values(n) = 'D'
   n = n + 1; Keys(n) = 'Alloy System Description';              Values(n) = 'D, Default structure'
   n = n + 1; Keys(n) = 'Val. Electron Rel (>= 0)';              Values(n) = '0'
!                     = 0: Non-relativitic; = 1: Scalar-relativitic; = 2: Full-relativitic
   n = n + 1; Keys(n) = 'Core Electron Rel (>= 0)';              Values(n) = '0'
!                     = 0: Non-relativitic; = 1: Full-relativitic
   n = n + 1; Keys(n) = 'Additional Electrons';                  Values(n) = '0.0'
!             Note: The additional number of electrons in a unit cell. This is an arbitrary number
!                   less, equal, or greater than 0.0. This parameter exists for testing purposes.
!                   It might be useful for some adventure. By default this paramater is set to zero.
   n = n + 1; Keys(n) = 'Charge Symmetry (>=0)';                 Values(n) = '1'
!                     = 0: No symmetry imposed; = 1: Symmetry impossed from spherical scattering calculation
   n = n + 1; Keys(n) = 'Single Site Solver (>= 0)';             Values(n) = '2'
!                     = 0: Spherical Schrodinger; = 1: Spherical Dirac; 2: Full Schrodinger; 3: Full Dirac
   n = n + 1; Keys(n) = 'Single Site Solution Method (>=-1)';    Values(n) = '2'
!             Note:: This parameter defines how the site potential is used in
!                    the iterative coupled integral solver. The solver always enforces
!                    a truncated potential if the irregular solutions are used. Default = 0.
!                     = -1: Untruncated potential in solver, surface integration
!                     =  0: Untruncated potential in solver, Untruncated potential in volume integration
!                     =  1: Untruncated potential in solver, Truncated potential in volume integration
!                     =  2: Truncated potential in solver
   n = n + 1; Keys(n) = 'Irregular Solutions (>=0)';             Values(n) = '1'
!                     = 0: Solver does NOT USE the IrrSolulions; = 1: Solver does USE the IrrSolulions
   n = n + 1; Keys(n) = 'Pole Search Step (>0.0)';               Values(n) = '0.010'
!             Note: Defaults to 0.01 for unspecified value, <= 0.0 or > 0.4
   n = n + 1; Keys(n) = 'Bound State Contour Integration Radius (>0.0)'; Values(n) = '0.001'
!             Note: This parameter defines the radius of the semi-circle contour integrating around the bound state
!                   energy to get the bound state density.
   n = n + 1; Keys(n) = 'Resonance State Max Width (>0.0)';      Values(n) = '0.050'
!             Note: If a pole is in the 4th quadrant and its distance from the real axis is less this value, 
!                   the pole is considered as a resonance, and the imaginary part is the width of the resonance 
   n = n + 1; Keys(n) = 'Resonance State Contour Integration Radius (>0.0)'; Values(n) = '0.002'
!             Note: This parameter defines the radius of the semi-circle contour around the resonance energy if
!                   its width is less than this value.
   n = n + 1; Keys(n) = 'No. Gauss Pts. along Resonance State Contour'; Values(n) = '5'
!             Note: This parameter specifies the number of Gaussian Quadrature points used for the integration
!                   along the semi-circle contour around the resonance energy.
   n = n + 1; Keys(n) = 'No. Gauss Pts. along Bound State Contour'; Values(n) = '5'
!             Note: This parameter specifies the number of Gaussian Quadrature points used for the integration
!                   along the semi-circle contour around the bound energy.
   n = n + 1; Keys(n) = 'Solutions Lmax Cutoff';                 Values(n) = '-1'
!                     = l<0: uses the lmax of the wave functions defined in info_table
!                     = l: finite, the coupled integral solver is used up to this l, while beyond it 
!                          uses only the radial solutions.
   n = n + 1; Keys(n) = 'Compute Phase Shifts (>=0)';            Values(n) = '0'
!                     = 0: Does not print out phase shifts; = 1: generates phase shifts and stops the code after it
   n = n + 1; Keys(n) = 'SS Lmax Potential Solver';              Values(n) = '-1'
!             Note: Defines the potential lmax used in the single site solver
!                   If it is less than 0 then defaults to the info_table lmax, otherwise
!                   the minimum of this lmax and potential lmax is used.
   n = n + 1; Keys(n) = 'Potential Type (>= 0)';                 Values(n) = '3'
!                     = 0: Muffin-tin;       = 1: ASA;               = 2: Muffin-tin ASA; 
!                     = 3: Full;             = 4: Muffin-Tin Test;   = 5: Empty Lattice;
!                     = 6: Mathieu Potential
   n = n + 1; Keys(n) = 'Exch-Corr. LDA Type (>= 0)';            Values(n) = '0'
!             Note: The input can be either one of the following numbers or, e.g.,
!                   LDA_X+LDA_C_HL for Hedin-Lundqvist LDA functional, or GGA_X_PBE+GGA_C_PBE
!                   for PBE GGA function, etc.. The naming convention here follows the definition given in LibXC.
!                     = 0: Barth-Hedin; = 1: Vosk-Wilk-Nusair; = 2: Perdew-Zunger; = 3: Perdew-Wang GGA; = 4: PBE
   n = n + 1; Keys(n) = 'LDA Improvement Scheme';                Values(n) = '0'
!                     = 0: No improvement; = 1: LDA + U; = 2: LDA + SIC; = 3: LDA + DFMT
   n = n + 1; Keys(n) = 'LDA+U Parameter File Name';             Values(n) = 'UJ.dat'
!             Note: If the name is 'None', the data will be obtained from other way
   n = n + 1; Keys(n) = 'Moment Direction File Name';            Values(n) = 'Evec_input.dat'
   n = n + 1; Keys(n) = 'Spin Index Param (>= 1)';               Values(n) = '1'
!                     = 1: No Spin; = 2: Spin-polarized; = 3: Spin-canted
   n = n + 1; Keys(n) = 'Interstitial Electron Spin';            Values(n) = '1'
!                     = 1: No Spin; = 2: Spin-polarized
   n = n + 1; Keys(n) = 'Canted Moment Torque Coef.';            Values(n) = '0.0'
   n = n + 1; Keys(n) = 'Calculate J_ij (y/n)';                  Values(n) = 'n'
   n = n + 1; Keys(n) = 'Read E-mesh from emeshs.inp';           Values(n) = '0'
!                     = 0. No; 1: Yes. In this case, the following parameters for the energy contour have no effect
   n = n + 1; Keys(n) = 'Contour Type (>= 0)';                   Values(n) = '0'
!                     = 0: Semi-circle; = 1: Rectangle Box; = 2: Horizontal Line; = 3: Vertical Line
   n = n + 1; Keys(n) = 'Number of Contours (> 0)';              Values(n) = '1'
   n = n + 1; Keys(n) = 'Energy Grid Type (>= 0)';               Values(n) = '1'
!                     = 0: Equal Interval; = 1: Gaussian Points; = 2: Log Interval; = 3: Nicholson Points
   n = n + 1; Keys(n) = 'No. Energy Grids';                      Values(n) = '30'
   n = n + 1; Keys(n) = 'No. Extra Energy Points';               Values(n) = '5'
   n = n + 1; Keys(n) = 'Offset Energy Point';                   Values(n) = '0'
   n = n + 1; Keys(n) = 'Real Axis Bottom, erbot';               Values(n) = '-0.40'
   n = n + 1; Keys(n) = 'Real Axis Top, ertop';                  Values(n) = '0.00'
   n = n + 1; Keys(n) = 'Imag Axis Bottom, eibot';               Values(n) = '0.001'
   n = n + 1; Keys(n) = 'Imag Axis Top, eitop';                  Values(n) = '1.000'
   n = n + 1; Keys(n) = 'Iterate Fermi energy';                  Values(n) = '0'
!                     = 1: On; = 0: Off
   n = n + 1; Keys(n) = 'SS Real Axis Int. Method';              Values(n) = '0'
!                     = 0: Unimesh; = 1: Adaptive; = 2: Uniform; = 3: Gaussian Quadrature; = 4: Romberg Method
   n = n + 1; Keys(n) = 'SS Real Axis Int. Points';              Values(n) = '300'
   n = n + 1; Keys(n) = 'T-matrix inversion (>= 0)';             Values(n) = '2'
   n = n + 1; Keys(n) = 'M-matrix inversion (>= 0)';             Values(n) = '10'
!                     = 0: LU method; = 2: QMR method
   n = n + 1; Keys(n) = 'No. Spin-dynamics Time Steps (>= 0)';   Values(n) = '1'
   n = n + 1; Keys(n) = 'Spin-dynamics Time Step';               Values(n) = '1.000'
   n = n + 1; Keys(n) = 'Mixing quantity type';                  Values(n) = '1'
!                     = 0: Charge mixing; = 1: Potential mixing
   n = n + 1; Keys(n) = 'Mixing algorithm';                      Values(n) = '2'
!                     = 0: Simple mixing; = 1: DGA mixing; = 2: Broyden mixing
   n = n + 1; Keys(n) = 'Lloyd correction';                      Values(n) = '0'
!                     = 0: Off; = 1: On
   n = n + 1; Keys(n) = 'Lloyd mode';                            Values(n) = '1'
!                     = 1: first method (to be defined by Markus); = 2: second method (to be defined by Markus)
   n = n + 1; Keys(n) = 'K-space Solver Method';                 Values(n) = '0'
!                     = 0: SuperLU; = 1: SuperLU with direct solve checking; = 2: Direct ( 1CPU only )
   n = n + 1; Keys(n) = 'Read K-mesh from kmeshs.inp';           Values(n) = '0'
!                     = 0: No; = 1: Yes. In this case, the following parameters for k-points generation have no effect
   n = n + 1; Keys(n) = 'Scheme to Generate K (>=0)';            Values(n) = '0'
!                     = 0: Special K-points; = 1: Tetrahedron; = 2: Direction
   n = n + 1; Keys(n) = 'No. K Meshs in IBZ (> 0)';              Values(n) = '1'
!             Note: This is not the number of k-points. This is the number of sets of k-points.
   n = n + 1; Keys(n) = 'Kx, Ky, Kz Division (> 0)';             Values(n) = '16  16  16'
   n = n + 1; Keys(n) = 'Symmetrize BZ Integration';             Values(n) = '1'
!                     = 0: No; = 1: Yes; = 2: Yes (Equiv. points)
   n = n + 1; Keys(n) = 'Large sphere radius (a.u.)';            Values(n) = '1000.0'
   n = n + 1; Keys(n) = 'Default Potential Input File Name';     Values(n) = 'D_fp_v'
   n = n + 1; Keys(n) = 'Default Potential Input File Form';     Values(n) = '1'
!                     = 0: ASCII Format; = 1: XDR Format; = 2: HDF Format; = 3: Machine Dependent Binary
   n = n + 1; Keys(n) = 'Default Potential Output File Name';    Values(n) = 'D_fp_w'
   n = n + 1; Keys(n) = 'Default Potential Output File Form';    Values(n) = '1'
!                     = 0: ASCII Format; = 1: XDR Format; = 2: HDF Format; = 3: Machine Dependent Binary
   n = n + 1; Keys(n) = 'Default Moment Direction';              Values(n) = '0.00  0.00  1.00'
   n = n + 1; Keys(n) = 'Default Constrain Field';               Values(n) = '0.00  0.00  0.00'
   n = n + 1; Keys(n) = 'Default Lmax-T matrix';                 Values(n) = '4'
!             Note: Lmax-T matrix  - controls the KKR size
   n = n + 1; Keys(n) = 'Default Lmax-Wave Func';                Values(n) = '4'
!             Note: Lmax-Wave Func - expansion of the solutions of Schrodinger eq. (Lmax= Lmax_T+n*Lmax_Pot, where n is
!                   the number of iterations of the iterative integral equations solver). This is set internally
!                   if is less than 0, otherwise takes the value specified in this file >= Lmax_T).
   n = n + 1; Keys(n) = 'Default Lmax-Potential';                Values(n) = '8'
!             Note: Lmax-Potential - expansion of the LDA potential to be used scf ( 0<= Lmax-Pot <= 2*Lmax-T ).
   n = n + 1; Keys(n) = 'Default Lmax-Trunc Pot';                Values(n) = '12'
   n = n + 1; Keys(n) = 'Default Lmax-Charge Den';               Values(n) = '8'
!             Note: Lmax-Charge    - expansion of the charge density (0<=Lmax-Pot<=2*Lmax-T).
   n = n + 1; Keys(n) = 'Default Lmax-Step Func';                Values(n) = '16'
!             Note: Lmax-Step Func - step function expansion used to converge
!                   the volume integration method of single site solver (0<= Lmax-Step < Infinity).
   n = n + 1; Keys(n) = 'Default LIZ # Neighbors';               Values(n) = '350'
!             Note: The neighbors number represents the maximum number of neighbors of a site. It can be less than
!                   the specified number if the neighbors on the next shell adds up over its maximum set.
   n = n + 1; Keys(n) = 'Default LIZ # NN Shells';               Values(n) = '16'
   n = n + 1; Keys(n) = 'Default LIZ Shell Lmax';                Values(n) = '4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4'
   n = n + 1; Keys(n) = 'Default LIZ Cutoff Radius';             Values(n) = '25.0'
   n = n + 1; Keys(n) = 'Default Rho  Mix Param.';               Values(n) = '0.1000'
   n = n + 1; Keys(n) = 'Default Pot  Mix Param.';               Values(n) = '0.1000'
   n = n + 1; Keys(n) = 'Default Mom  Mix Param.';               Values(n) = '0.1500'
   n = n + 1; Keys(n) = 'Default Chg  Mix Param.';               Values(n) = '1.0000'
   n = n + 1; Keys(n) = 'Default Evec Mix Param.';               Values(n) = '0.0000'
   n = n + 1; Keys(n) = 'Default Maximum Core Radius';           Values(n) = '0.0d0'
   n = n + 1; Keys(n) = 'Default Maximum Muffin-tin Radius';     Values(n) = '0.0d0'
   n = n + 1; Keys(n) = 'Default No. Rad Points ndivin ';        Values(n) = '1001'
!                     = 0: Not specified; > 0: Speciflied. Note: 0 < r(j)=exp(j*hin) <= rmt, j=1,2,...,ndivin
   n = n + 1; Keys(n) = 'Default No. Rad Points ndivout';        Values(n) = '0'
!                     = 0: Not specified; > 0: Speciflied. Note: rmt < r(j)=exp(j*hout) <= rmax, j=1,2,...,ndivout
   n = n + 1; Keys(n) = 'Default Integer Factor nmult';          Values(n) = '1'
!                     = 0: Not specified; > 0: Speciflied. Note: r(j) = exp(j*hin), hin = nmult*hout
   n = n + 1; Keys(n) = 'Default Radial Grid Exponential Step';  Values(n) = '0.01'
!                     = 0.0: Not specified; 
!                     > 0.0: 0.005 - 0.02 is recommended. Note: r(j) = exp(j*hin), here hin is the 
!                            exponential step. In this case, you may want to leave ndivin unspecified
   n = n + 1; Keys(n) = 'Default Pseudo Charge Radius';          Values(n) = '0.9'
!             Note: This is the ratio of the pseudo charge radius and the muffin-tin radius
   n = n + 1; Keys(n) = 'Default Screen Pot.';                   Values(n) = '3.0'
   n = n + 1; Keys(n) = 'Default Lmax-Screen';                   Values(n) = '3'
   n = n + 1; Keys(n) = 'Default Rcut-Screen';                   Values(n) = '4.8'
   n = n + 1; Keys(n) = 'Local SIC';                             Values(n) = '0'
   n = n + 1; Keys(n) = 'Default Mixing Parameter';              Values(n) = '0.1000'
   n = n + 1; Keys(n) = 'Frozen-Core Calculation';               Values(n) = '0'
!             Note: = 0: Not a frozen core calculation; > 0: a SCF iteration value beyond which the frozon core applies
   n = n + 1; Keys(n) = 'Frozen-Core File Name';                 Values(n) = ' '
!             Note: This file will be read if a file name is given
   n = n + 1; Keys(n) = 'Maximum Effective Medium Iterations';   Values(n) = '40'           
!             Note: If this value is set to be 0, an ATA (instead of CPA) calculation will be performed.
   n = n + 1; Keys(n) = 'Effective Medium Mixing Scheme';        Values(n) = '2'
!                     = 0: Simple mixing; = 1: Anderson mixing; = 2: Broyden mixing; = 3: Anderson Mixing by Messina group
   n = n + 1; Keys(n) = 'Effective Medium Mixing Parameters';    Values(n) = '0.1000  0.01'
!             Note: The first mixing value is for the energy points in standard mixing mode; the second mixing value
!                   is for the energy points in conservative mixing mode
   n = n + 1; Keys(n) = 'Effective Medium Mixing eSwitch Value'; Values(n) = '0.003'
!             Note: If Re[E] > 0 and Im[E] < eSwitch, the effective medium iteration is switched to the conservative mode
   n = n + 1; Keys(n) = 'Number of iterations with aggressive mixing'; Values(n) = '25'
!             Note: This is the number of CPA iterations for each aggressive mixing scheme, e.g., Broyden and Anderson 
   n = n + 1; Keys(n) = 'Maximum number of mixing scheme changes'; Values(n) = '10'
!             Note: This is the maximum number of mixing scheme changes
   n = n + 1; Keys(n) = 'Effective Medium T-matrix Tol (> 0)';   Values(n) = '0.0000001'
   n = n + 1; Keys(n) = 'Default Core Radius';                   Values(n) = '0'
!                     =-1: Using the circumscribed sphere radius for full-potential calculations
!                     = 0: Using the inscribed sphere radius for full-potential calculations,
!                          the muffin-tin radius for muffin-tin calculations, or the ASA radius for ASA calculations.
!                     = 1: Using the implicit core radius defined in ChemElementModule
!                     = A specific real value (> 0.0, in atomic units). However, if this value < 1.0, it will be treated
!                       as a multiplier, and the core radius is this value times the inscribed sphere radius.
   n = n + 1; Keys(n) = 'Default Muffin-tin Radius';             Values(n) = '0'
!                     = 0: Using the inscribed sphere radius
!                     = 1: Using the implicit muffin-tin radius defined in ChemElementModule
!                     = A specific real value (> 0.0, in atomic units). However, if this value < 1.0, it will be treated
!                       as a multiplier, and the muffin-tin radius is this value times the inscribed sphere radius.
   n = n + 1; Keys(n) = 'Default Radical Plane Ratio';           Values(n) = '1.000'
!  Note: "Ratio" here is the radical plane distance from atom ralative to a 
!        universal value applied to the entire system. This universal value does 
!        not have to be specified, so long as all the atoms in the system have 
!        the same denominator when choosing this value. If all atoms have the same
!        radical plane ratio (RPR), a radical plane is set at half way between two
!        atoms. If each atom has its own RPR, which is in effect treated as a 
!        "weight" for the atom when setting up a radical plane. Specifically, a 
!        radical plane between atom i and atom j, separated by distance D, will be 
!        set at distance from atom i at D*RPR(i)/(RPR(i)+RPR(j)), and at distance 
!        from atom j at D*RPR(j)/(RPR(i)+RPR(j))
   n = n + 1; Keys(n) = 'Uniform Grid Origin';                   Values(n) = '0'
!                     = 0: Using unit cell corner             
!                     = 1: Using Unit cell center
   n = n + 1; Keys(n) = 'Uniform Grid Origin Vector';            Values(n) = '0.0  0.0  0.0'
   n = n + 1; Keys(n) = 'Core States Normalization Range';       Values(n) = '0'
!                     = 0: Up to bounding sphere radius
!                     = 1: Up to infinity
   n = n + 1; Keys(n) = 'Mixing Parameter for Finding Ef';       Values(n) = '0.5'
!        Note: At each SCF step, the calculated fermi energy will be mixing with the old fermi energy
!              to determine the new fermi energy for the next SCF iteration
   n = n + 1; Keys(n) = 'Mixing Switch for Finding Ef';          Values(n) = '0.01'
!        Note: If the difference between the calculated fermi energy and the old fermi energy is
!              greater than this switching value, mixing will be engaged. otherwise, no fermi energy
!              mixing is applied. One may keep the fermi energy fixed by setting this switching value 
!              to 0.0 and also setting fermi energy mixing parameter to 0.0.
   n = n + 1; Keys(n) = 'Renormalize Green function';            Values(n) = '0'
!        Note: The Green function may be renormalized using the the integrated DOS value
!              calculated from multiple scattering Lloyd formula and/or single scattering
!              Krein formula.
!                     = 0: Do NOT renormalize the Green function
!                     = 1: Renormalize the Green function
   n = n + 1; Keys(n) = 'Full-potential Semi-core';              Values(n) = '0'
!                     = 0: Do NOT treat semi-core states with full-potential;
!                     = 1: Treat semi-core states with full-potential.
   n = n + 1; Keys(n) = 'Include CPA/SRO Charge Correction';     Values(n) = '0'
!                     = 0: Do NOT include charge correction term in potential for the CPA calculation;
!                     = 1: Include the charge correction term in potential/energy for the CPA calculation;
   n = n + 1; Keys(n) = 'Use Linear Relation';                   Values(n) = '0'
!                     = 0: DO NOT use the linear q-V relation to determine the charge correction term;
!                     = 1: Use the linear q-V relation to determine the charge correction term;
   n = n + 1; Keys(n) = 'Tolerance for setting up polyhedra';    Values(n) = '0.0000005d0'
!        Note: This is the tolerance value used in setting up atomic polyhedra in the unit cell.
!              If one gets error message for determining the corners or edges of a polyhedron, it is
!              likely that this tolerance value needs to be adjusted.
   n = n + 1; Keys(n) = 'Ewald parameter for KKR';               Values(n) = '0.5d0'
   n = n + 1; Keys(n) = 'Calculate Superconducting Tc';          Values(n) = '0'
!                     = 0: No (Default)
!                     = 1: Yes
   n = n + 1; Keys(n) = 'mu* (e-e interaction constant)';        Values(n) = '0'
!                     = 0: Bennemann and Garland formula (Default)
!                     = 1: Modified BG formula by Papaconstantopoulos
!                     = A real positive number (< 1.0)
   n = n + 1; Keys(n) = 'Average of phonon frequency (1/sec)';   Values(n) = '0.0'
!        Note; the input format is a real positive value
   n = n + 1; Keys(n) = 'Average of phonon frequency (K)';  Values(n) = '0.0'
!        Note; the input format is a real positive value
   n = n + 1; Keys(n) = 'Average of phonon frequency squared (1/sec^2)';  Values(n) = '0.0'
!        Note; the input format is an atomic index value (corresponding to the atom in the position data), or 
!              an chemical element symbol, followed by a real positve value
   n = n + 1; Keys(n) = 'Average of phonon frequency squared (K^2)'; Values(n) = '0.0'
!        Note; the input format is an atomic index value (corresponding to the atom in the position data), or 
!              an chemical element symbol, followed by a real positve value
   n = n + 1; Keys(n) = 'Atomic mass times <omega^2> (eV/Anst^2)';  Values(n) = '0.0'
!        Note; the input format is an atomic index value (corresponding to the atom in the position data), or 
!              an chemical element symbol, followed by a real positve value
   n = n + 1; Keys(n) = 'Atomic mass times <omega^2> (Ryd/BohrRad^2)';    Values(n) = '0.0'
!        Note; the input format is an atomic index value (corresponding to the atom in the position data), or 
!              an chemical element symbol, followed by a real positve value
   n = n + 1; Keys(n) = 'Debye Temperature (K)';                 Values(n) = '0'
!                     = 0: getting from the internal database
!                     = A real positive number provided by user
   n = n + 1; Keys(n) = 'Negative charge density tolerance';     Values(n) = '0.00001'
!        Note: this is a tolerence value for rho(r) < 0: if -tol < rho(r) < 0, we will set rho(r) = 0.
   n = n + 1; Keys(n) = 'Imaginary energy shift';     Values(n) = '0.001'
!        Note: this is a tolerence value for rho(r) < 0: if -tol < rho(r) < 0, we will set rho(r) = 0.
   n = n + 1; Keys(n) = 'Perform DFT+DMFT Calculation';          Values(n) = '0'
!                     = 0: No;
!                     = 1: Yes.
   n = n + 1; Keys(n) = 'Default Local Orbital Labels';          Values(n) = 'd'
!        Note: The values are: s, p, d, and/or f, separated by ' ', ',', or ';'
   n = n + 1; Keys(n) = 'Default Local s-Orbital Energy';        Values(n) = '0.2'
   n = n + 1; Keys(n) = 'Default Local p-Orbital Energy';        Values(n) = '0.3'
   n = n + 1; Keys(n) = 'Default Local d-Orbital Energy';        Values(n) = '0.5'
   n = n + 1; Keys(n) = 'Default Local f-Orbital Energy';        Values(n) = '0.4'
   n = n + 1; Keys(n) = 'Default Local Orbital Energy';          Values(n) = '0.5'
!        Note: The orbital energy value is real and is in Rydberg unit
!  =============================================================================
!  Additional input parameters can be added to this table using the 
!  following format:
!    n = n + 1; Keys(n) = 'Key here';  Values(n) = 'Value here'
!  =============================================================================
