============================================
Input Parameters
============================================

This document provides a comprehensive list of input parameters for the MuST code. 
It includes default values, allowed ranges, and technical logic extracted from 
the source input files.

Each parameter is defined in the main input file (typically start with ``i_``, for example, ``i_mst``) using the format:

   ``Key :: Value``

.. contents:: Table of Contents
   :depth: 2

General Execution & File I/O
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Current File Path
     - ./
     - —
   * - Output to Screen (y/n)
     - n
     - y or n
   * - Output Level (>= -1)
     - 0
     - Must be ≥ -1
   * - Output Proc. ID (>= -1)
     - 0
     - -1: All processors write to individual output files; ≥ 0: logical processor indices (space/comma separated) that write to unit 6
   * - Output Atom ID (>= -1)
     - 0
     - -1: Output all atoms in unit cell; 0: Output atoms mapped to output processor; ≥ 1: global atom indices (space/comma separated)


SCF Calculations
===============================
.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - No. Iterations (> 0)
     - 60
     - Must be > 0
   * - Method of SCF Calculation
     - 2
     - 1: LSMS; 2: KKR; 3: KKR-CPA
   * - No. Gauss Pts. along r
     - 80
     - — 
   * - No. Gauss Pts. along theta
     - 60
     - — 
   * - Valence Band Bottom Est.
     - -0.5
     - — 
   * - Temperature Parameter (K)
     - 0.000
     - — 
   * - DOS Run ID
     - 0
     - 0: DOS along real energy axis will not be calculated; >0: integer representing global atom index for DOS calculation; -1: DOS calculated for all atoms
   * - Ewald parameter for KKR
     - 0.5d0
     - Only needed in KKR/KKR-CPA calculations 

Tolerance Parameters
===========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Energy (Ryd) Tol (> 0)
     - 0.000001
     - Must be > 0
   * - Potential Tol (> 0)
     - 0.0000001
     - Must be > 0
   * - Fermi Energy Tol (> 0)
     - 0.0000001
     - Must be > 0
   * - Other RMS Tol (> 0)
     - 0.0000001
     - Must be > 0
   * - Tolerance for setting up polyhedra
     - 0.0000005d0
     - Note: Adjust if errors occur when determining corners/edges of polyhedra
   * - Negative charge density tolerance
     - 0.00001
     - Note: If -tol < rho(r) < 0, set rho(r) = 0

System
=============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - No. Atoms in System (> 0)
     - 1
     - Must be > 0
   * - Atomic Position File Name
     - position.dat
     - — 
   * - Text Identification
     - D
     - — 
   * - Alloy System Description
     - D, Default structure
     - — 
   * - Default Radical Plane Ratio
     - 1.000
     - Note: Ratio of radical plane distance relative to universal value applied to system

Full-potential Method Parameters
===================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Single Site Solver (>= 0)
     - 2
     - 0: Spherical Schrodinger; 1: Spherical Dirac; 2: Full Schrodinger; 3: Full Dirac
   * - Single Site Solution Method (>= -1)
     - 2
     - -1: Untruncated potential in solver, surface integration; 0: Untruncated potential in solver, Untruncated potential in volume integration; 1: Untruncated potential in solver, Truncated potential in volume integration; 2: Truncated potential in solver. Note: Defines how site potential is used in the iterative coupled integral solver; truncated potential is enforced if irregular solutions are used; default = 0
   * - Irregular Solutions (>=0)
     - 1
     - 0: Solver does NOT use the Irregular Solutions; 1: Solver does use the Irregular Solutions
   * - Pole Search Step (>0.0)
     - 0.010
     - Note: Defaults to 0.01 for unspecified value; valid range <=0.0 or >0.4
   * - Bound State Contour Integration Radius (>0.0)
     - 0.001
     - Note: Radius of the semi-circle contour around the bound state energy for density integration
   * - Resonance State Max Width (>0.0)
     - 0.050
     - Note: If a pole in 4th quadrant has distance from real axis < this value, considered as a resonance; imaginary part is width
   * - Resonance State Contour Integration Radius (>0.0)
     - 0.002
     - Note: Radius of semi-circle contour around resonance energy if its width < this value
   * - No. Gauss Pts. along Resonance State Contour
     - 5
     - Note: Number of Gaussian quadrature points for integration along semi-circle contour around resonance energy
   * - No. Gauss Pts. along Bound State Contour
     - 5
     - Note: Number of Gaussian quadrature points for integration along semi-circle contour around bound energy
   * - Solutions Lmax Cutoff
     - -1
     - l<0: uses lmax of wave functions in info_table; l≥0: coupled integral solver used up to this l, beyond it only radial solutions used
   * - Compute Phase Shifts (>=0)
     - 0
     - 0: Does not print phase shifts; 1: generates phase shifts and stops code
   * - SS Lmax Potential Solver
     - -1
     - Note: Defines potential lmax in single site solver; if <0, defaults to info_table lmax; otherwise minimum of this lmax and potential lmax is used
   * - Uniform Grid Parameters
     - 64 64 64
     - >0: Three integers defining grid numbers along Bravais lattice vectors; must be powers of 2; affects Poisson solver accuracy
   * - Default Lmax-Wave Func
     - 4
     - Note: Expansion of solutions of Schrödinger equation; internally computed if <0, otherwise takes specified value ≥ Lmax_T
   * - Default Lmax-Potential
     - 8
     - Note: Expansion of LDA potential used in SCF (0 ≤ Lmax-Pot ≤ 2*Lmax-T)
   * - Default Lmax-Trunc Pot
     - 12
     - — 
   * - Default Lmax-Charge Den
     - 8
     - Note: Expansion of charge density (0 ≤ Lmax-Charge ≤ 2*Lmax-T)
   * - Default Pseudo Charge Radius
     - 0.9
     - Note: Ratio of pseudo charge radius to muffin-tin radius
   * - Uniform Grid Origin
     - 0
     - 0: unit cell corner; 1: unit cell center
   * - Uniform Grid Origin Vector
     - 0.0 0.0 0.0
     - — 
   * - Full-potential Semi-core
     - 0
     - 0: do NOT treat semi-core with full-potential; 1: treat semi-core with full-potential

Electronic Structure & Potential
================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Val. Electron Rel (>= 0)
     - 0
     - 0: Non-relativistic; 1: Scalar-relativistic; 2: Full-relativistic
   * - Charge Symmetry (>=0)
     - 1
     - 0: No symmetry imposed; 1: Symmetry imposed from spherical scattering calculation
   * - Potential Type (>= 0)
     - 0
     - 0: Muffin-tin; 1: ASA; 2: Muffin-tin ASA; 3: Full; 4: Muffin-Tin Test; 5: Empty Lattice; 6: Mathieu Potential
   * - Exch-Corr. LDA Type (>= 0)
     - 0
     - 0: Barth-Hedin; 1: Vosk-Wilk-Nusair; 2: Perdew-Zunger; 3: Perdew-Wang GGA; 4: PBE  
       Note: Can also use LibXC-style functional names like LDA_X+LDA_C_HL or GGA_X_PBE+GGA_C_PBE

Spin and Magnetic-Moment Parameters
=======================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Moment Direction File Name
     - Evec_input.dat
     - Require only for spin-canted case 
   * - Spin Index Param (>= 1)
     - 1
     - 1: No Spin; 2: Spin-polarized; 3: Spin-canted
   * - Interstitial Electron Spin
     - 1
     - 1: No Spin; 2: Spin-polarized
   * - Default Moment Direction
     - 0.00 0.00 1.00
     - — 
   * - Default Constrain Field
     - 0.00 0.00 0.00
     - — 

Energy Contour Parameters
==========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Read E-mesh from emeshs.inp
     - 0
     - 0: No; 1: Yes. If yes, following parameters for energy contour have no effect
   * - Contour Type (>= 0)
     - 0
     - 0: Semi-circle; 1: Rectangle Box; 2: Horizontal Line; 3: Vertical Line
   * - Number of Contours (> 0)
     - 1
     - Must be > 0
   * - Energy Grid Type (>= 0)
     - 1
     - 0: Equal Interval; 1: Gaussian Points; 2: Log Interval; 3: Nicholson Points
   * - No. Energy Grids
     - 30
     - — 
   * - No. Extra Energy Points
     - 5
     - — 
   * - Offset Energy Point
     - 0
     - — 
   * - Real Axis Bottom, erbot
     - -0.40
     - — 
   * - Real Axis Top, ertop
     - 0.00
     - — 
   * - Imag Axis Bottom, eibot
     - 0.001
     - — 
   * - Imag Axis Top, eitop
     - 1.000
     - — 
   * - Iterate Fermi energy
     - 0
     - 0: Off; 1: On
   * - SS Real Axis Int. Method
     - 0
     - 0: Unimesh; 1: Adaptive; 2: Uniform; 3: Gaussian Quadrature; 4: Romberg Method
   * - SS Real Axis Int. Points
     - 300
     - — 
   * - Imaginary energy shift
     - 0.001
     - Note: If -tol < rho(r) < 0, set rho(r) = 0 (DOS calculation)

SCF Mixing Schemes and Parameters
==================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Mixing quantity type
     - 1
     - 0: Charge mixing; 1: Potential mixing
   * - Mixing algorithm
     - 2
     - 0: Simple mixing; 1: DGA mixing; 2: Broyden mixing
   * - Lloyd correction
     - 0
     - 0: Off; 1: On
   * - Lloyd mode
     - 1
     - 1: first method (to be defined by Markus); 2: second method (to be defined by Markus)
   * - Default Rho Mix Param.
     - 0.1000
     - — 
   * - Default Pot Mix Param.
     - 0.1000
     - — 
   * - Default Mom Mix Param.
     - 0.1500
     - — 
   * - Default Chg Mix Param.
     - 1.0000
     - — 
   * - Default Evec Mix Param.
     - 0.0000
     - — 
   * - Default Mixing Parameter
     - 0.1000
     - Note: Default value for mixing
   * - Mixing Parameter for Finding Ef
     - 0.5
     - Note: New Fermi energy = mixing of calculated and old Fermi energies at each SCF step
   * - Mixing Switch for Finding Ef
     - 0.01
     - Note: If difference > switch value, mixing is applied; 0 disables Fermi energy mixing

Kpoint Parameters (Required only for KKR and KKR-CPA)
=======================================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - K-space Solver Method
     - 0
     - 0: SuperLU; 1: SuperLU with direct solve checking; 2: Direct (1 CPU only)
   * - Read K-mesh from kmeshs.inp
     - 0
     - 0: No; 1: Yes. If yes, following parameters for k-points generation have no effect
   * - Scheme to Generate K (>=0)
     - 0
     - 0: Special K-points; 1: Tetrahedron; 2: Direction
   * - No. K Meshs in IBZ (> 0)
     - 1
     - Note: Number of sets of k-points, not total k-points
   * - Kx, Ky, Kz Division (> 0)
     - 16 16 16
     - — 
   * - Symmetrize BZ Integration
     - 1
     - 0: No; 1: Yes; 2: Yes (Equivalent points)

Potential and LMAX Parameters
================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Large sphere radius (a.u.)
     - 1000.0
     - — 
   * - Default Potential Input File Name
     - D_fp_v
     - — 
   * - Default Potential Input File Form
     - 1
     - 0: ASCII Format; 1: XDR Format; 2: HDF Format; 3: Machine Dependent Binary
   * - Default Potential Output File Name
     - D_fp_w
     - — 
   * - Default Potential Output File Form
     - 1
     - 0: ASCII Format; 1: XDR Format; 2: HDF Format; 3: Machine Dependent Binary
   * - Default Lmax-T matrix
     - 4
     - Note: Controls KKR size
   * - Default Lmax-Step Func
     - 16
     - Note: Step function expansion used for volume integration in single site solver (0 ≤ Lmax-Step < Infinity)

Local Interaction Zone (LIZ) Parameters (Required only for LSMS)
=================================================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Default LIZ # Neighbors
     - 350
     - Note: Maximum number of neighbors of a site; actual number may be less if the next shell exceeds maximum
   * - Default LIZ # NN Shells
     - 16
     - — 
   * - Default LIZ Shell Lmax
     - 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
     - — 
   * - Default LIZ Cutoff Radius
     - 25.0
     - — 

Radial Grid Parameters
=======================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Default No. Rad Points ndivin
     - 1001
     - 0: Not specified; >0: Specified. Note: 0 < r(j)=exp(j*hin) ≤ rmt, j=1,2,...,ndivin
   * - Default No. Rad Points ndivout
     - 0
     - 0: Not specified; >0: Specified. Note: rmt < r(j)=exp(j*hout) ≤ rmax, j=1,2,...,ndivout
   * - Default Integer Factor nmult
     - 1
     - 0: Not specified; >0: Specified. Note: r(j)=exp(j*hin), hin = nmult*hout
   * - Default Radial Grid Exponential Step
     - 0.01
     - 0.0: Not specified; >0.0: recommended 0.005–0.02. Note: r(j)=exp(j*hin), hin = exponential step
   * - Default Muffin-tin Radius
     - 0
     - 0: inscribed sphere radius; 1: implicit muffin-tin radius; >0: specific value in a.u.

Core States Solver Parameters
========================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Frozen-Core Calculation
     - 0
     - Note: 0: Not a frozen core calculation; >0: SCF iteration beyond which frozen core applies
   * - Frozen-Core File Name
     - ' '
     - Note: File to read if a frozen-core file name is given
   * - Core Electron Rel (>= 0)
     - 0
     - 0: Non-relativistic; 1: Full-relativistic
   * - Default Maximum Core Radius
     - 0.0d0
     - — 
   * - Default Core Radius
     - 0
     - -1: circumscribed sphere radius; 0: inscribed sphere / muffin-tin / ASA radius; 1: implicit core radius; >0: specific value in a.u.
   * - Core States Normalization Range
     - 0
     - 0: up to bounding sphere radius; 1: up to infinity

Effective Medium Parameters
=================================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Maximum Effective Medium Iterations
     - 40
     - Note: 0 → ATA (instead of CPA) calculation will be performed
   * - Effective Medium Mixing Scheme
     - 2
     - 0: Simple mixing; 1: Anderson mixing; 2: Broyden mixing; 3: Anderson Mixing by Messina group
   * - Effective Medium Mixing Parameters
     - 0.1000 0.01
     - Note: First value = energy points in standard mixing; second value = energy points in conservative mode
   * - Effective Medium Mixing eSwitch Value
     - 0.003
     - Note: If Re[E] > 0 and Im[E] < eSwitch, iteration switches to conservative mode
   * - Number of iterations with aggressive mixing
     - 25
     - Note: Number of CPA iterations per aggressive mixing scheme
   * - Maximum number of mixing scheme changes
     - 10
     - Note: Maximum number of mixing scheme changes
   * - Effective Medium T-matrix Tol (>0)
     - 0.0000001
     - — 
   * - Include CPA/SRO Charge Correction
     - 0
     - 0: do NOT include correction; 1: include correction term in potential/energy for CPA
   * - Use Linear Relation
     - 0
     - 0: do NOT use linear q-V relation; 1: use linear q-V relation for charge correction


MPI-GPU Parameters
==========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Maximum MPI tasks per GPU for KKR Matrix Inverse
     - 28
     - Note: If the number of MPI tasks exceeds this, GPU acceleration of KKR matrix inverse is disabled
   * - Maximum MPI tasks per GPU for KKR Matrix Calculation
     - 4
     - Note: If MPI tasks exceed this, KKR matrix is calculated on CPU and copied to GPU for inversion
   * - Maximum MPI tasks per GPU for Gij Matrix Calculation
     - 4
     - Note: If MPI tasks exceed this, Gij matrix is calculated on CPU and copied to GPU for inversion

Under development Methods
========================================


Spin Dynamics Parameters
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Canted Moment Torque Coef.
     - 0.0
     - — 
   * - Calculate J_ij (y/n)
     - n
     - y or n
   * - No. Spin-dynamics Time Steps (>= 0)
     - 1
     - — 
   * - Spin-dynamics Time Step
     - 1.000
     - — 

Screened KKR Parameters
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - SuperLU Tol (> 0)
     - 0.0000001
     - Must be > 0
   * - K-space Check Tol (> 0)
     - 0.0000001
     - Must be > 0
   * - Default Screen Pot.
     - 3.0
     - Note: Default value of the screened potential
   * - Default Lmax-Screen
     - 3
     - Note: Maximum angular momentum for screening
   * - Default Rcut-Screen
     - 4.8
     - Note: Cutoff radius for screening (a.u.)


DFT+DMFT Parameters
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Perform DFT+DMFT Calculation
     - 0
     - 0: No; 1: Yes
   * - Default Local Orbital Labels
     - d
     - Note: Allowed values are s, p, d, f, separated by ' '(space), ','(comma), or ';'(semicolon)
   * - Default Local s-Orbital Energy
     - 0.2
     - Note: Real value in Rydberg units
   * - Default Local p-Orbital Energy
     - 0.3
     - Note: Real value in Rydberg units
   * - Default Local d-Orbital Energy
     - 0.5
     - Note: Real value in Rydberg units
   * - Default Local f-Orbital Energy
     - 0.4
     - Note: Real value in Rydberg units
   * - Default Local Orbital Energy
     - 0.5
     - Note: Real value in Rydberg units
   * - Initial Fermi Energy Setting
     - 0
     - 0: Set to average of Fermi energy from input potential
       1: Set to average of Max and Min of Fermi energy from input potential
       >0: User-specified real number as initial Fermi energy

Beyond DFT
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Key
     - Default Value
     - Options / Notes
   * - Calculate Superconducting Tc
     - 0
     - 0: No (Default); 1: Yes
   * - mu* (e-e interaction constant)
     - 0
     - 0: Bennemann and Garland formula (Default); 0,A: Modified BG formula with factor A (<1.0)
   * - Average of phonon frequency (1/sec)
     - 0.0
     - Note: Input format is a real positive value
   * - Average of phonon frequency (K)
     - 0.0
     - Note: Input format is a real positive value
   * - Average of phonon frequency squared (1/sec^2)
     - 0.0
     - Note: Input format: atomic index or chemical element symbol, followed by a real positive value
   * - Average of phonon frequency squared (K^2)
     - 0.0
     - Note: Input format: atomic index or chemical element symbol, followed by a real positive value
   * - Atomic mass times <omega^2> (eV/Anst^2)
     - 0.0
     - Note: Input format: atomic index or chemical element symbol, followed by a real positive value
   * - Atomic mass times <omega^2> (Ryd/BohrRad^2)
     - 0.0
     - Note: Input format: atomic index or chemical element symbol, followed by a real positive value
   * - Debye Temperature (K)
     - 0
     - 0: get from internal database; >0: user-provided real positive value
   * - LDA Improvement Scheme
     - 0
     - 0: No improvement; 1: LDA + U; 2: LDA + SIC; 3: LDA + DFMT
   * - LDA+U Parameter File Name
     - UJ.dat
     - Note: If name is 'None', data will be obtained by other methods
   * - Local SIC
     - 0
     - Note: 0: Off; 1: On (local self-interaction correction)



Visualization Parameters
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Generate System Movie
     - 0
     - 0: No system movie output (default); 1: generate movie
   * - Stop-at Routine Name
     - main
     - —
   * - No. Iter for Each Pot. Write
     - 5
     - —
   * - No. Iter for Each Movie
     - 0
     - 0: No system movie output (default); ≥ 1: a system movie is written at the end of every this number of SCF iterations
   * - Visual Grid Type (0<D<4)
     - 3
     - 1: Define a line (1D); 2: Define a plane (2D); 3: Define a rectilinear region (3D)
   * - Visual Grid Scale
     - 1.00
     - —
   * - Visual Grid Origin Vector
     - 0.0 0.0 0.0
     - —
   * - Visual Grid Vector 1
     - 10.00 0.00 0.00
     - Note: Vector defining grid boundaries
   * - Visual Grid Vector 2
     - 0.00 10.00 0.00
     - Note: Vector defining grid boundaries
   * - Visual Grid Vector 3
     - 0.00 0.00 10.00
     - Note: Vector defining grid boundaries
   * - Visual Grid Points
     - 10 10 10
     - Note: Only first D parameters are used for Grid Vectors and Grid Points
   * - Visual Line Vector
     - 1.00 0.00 0.00
     - —
   * - Visual Line Points
     - 10
     - —
   * - Output Electron Density ID (>= -1)
     - -1
     - -1: do not print; 0: print brief data; 1: print data on visual grid
   * - Output Density Format
     - 2
     - 0: x y z rho format; 1: Legacy ".vtk" format for ParaView; 2: ".xsf" format for VESTA

Under development Parameters
----------------------------------
.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Calc. Harris Energy (H.E.)
     - 0
     - 0: Do not calculate H.E.; 1: H.E. at Updated μ; 2: H.E. at Fixed μ
   * - Additional Electrons
     - 0.0
     - Note: Additional number of electrons in a unit cell; can be <, =, or > 0; mainly for testing; default = 0
   * - T-matrix inversion (>= 0)
     - 2
     - — 
   * - M-matrix inversion (>= 0)
     - 10
     - 0: LU method; 2: QMR method
   * - Renormalize Green function
     - 0
     - 0: do NOT renormalize; 1: renormalize using integrated DOS from Lloyd/Krein formula


Atomic Position File (position.dat)
===================================

The ``position.dat`` file defines the structural framework and chemical composition of the system.

Structure Definition
--------------------
* **Lattice Constant**: The first line defines the global scaling factor (e.g., 5.53).
* **Bravais Vectors**: Three lines defining the primitive lattice vectors.
* **Site Position**: Cartesian coordinates (x, y, z) for each atom in the cell.

Crystal Structure: CuZn (BCC)
--------------------------------------

.. code-block:: text

    5.53
    # Bravais lattice
        1.00000000000    0.00000000000    0.0000000000
        0.00000000000    1.00000000000    0.0000000000
        0.00000000000    0.00000000000    1.0000000000
    # Atomic position
    Cu  0.00000000000    0.00000000000    0.00000000000
    Zn  0.50000000000    0.50000000000    0.50000000000

Crystal Structure: Cu-Zn Random Alloy (FCC, CPA)
---------------------------------------------------

.. code-block:: text

    5.53
    # Bravais lattice
         0.50000000000    0.50000000000   -0.50000000000
         0.50000000000   -0.50000000000    0.50000000000
        -0.50000000000    0.50000000000    0.50000000000
    # Atomic position
    CPA  0.00000000000    0.00000000000    0.00000000000  Cu  0.50000  Zn  0.50000
    
Notes
-----
- For random alloy calculations, use 'CPA' as the virtual atom name.
- Coordinates are Cartesian (x, y, z) of the virtual atom.
- The atomic species and their site concentrations follow the coordinates.
