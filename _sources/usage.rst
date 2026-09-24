********************
Package Description
********************

Available Scientific Packages
=============================

- **MST**: Perform KKR, LSMS, single-site, and Cluster Averaged KKR-CPA.
- **lsms**: Perform LSMS and Wang-Landau LSMS. This package is optimized for extreme performance on petascale and exascale systems.
- **KUBO**: Perform first-principles electrical conductivity calculations.

User Support Folders
====================

- **potentials**: Contains the starting potential for selected elements.
- **architecture**: Contains preset makefile parameters ("architecture files") for a wide variety of computer systems.
- **docs**: Contains installation instructions, license information, and user guides.
- **external**: Contains external libraries required or optionally required by MuST, e.g., FFTW, Lua, P3DFFT, and LibXC.
- **tutorials**: Contains hands-on exercises and training materials.
- **ase_must**: Provides Atomic Simulation Environment (ASE) support for MuST.

Executables under ``bin/``
##########################

1. ``mst``
==========

The corresponding main source code is ``MST/src/mst2.F90``. It performs *ab initio*
electronic structure calculations for 3-D structures.

**Main features:**

- Linear scaling calculation based on LSMS method
- Calculations based on KKR, KKR-CPA, or LSMS method
- Muffin-tin potential or full-potential
- Non-relativistic, scalar-relativistic, or relativistic
- Non-spin polarized, spin-polarized, or spin-canted
- Special k-points method for Brillouin zone integration
- LDA or GGA exchange-correlation potentials

**Input files:**

- ``i_*``: Main input file containing control parameters for SCF calculations,
  system definition, atomic positions, and potential file names
- Position data file: File name specified in ``i_*``
- Potential file(s): File name and format specified in input
- ``info_*`` (obsolete): Atom-based control parameters (specified in ``i_*``)
- ``kmeshs.inp`` (optional): Used for testing purposes
- ``emeshs.inp`` (optional): Used for testing purposes
- ``Evec_*`` (optional): Used for spin-canted calculations

**Output files:**

- ``o_n*``: Main output file (may not be created if output is redirected to screen)
- ``k_n*``: Contains total energy and Fermi energy for each SCF iteration
- New potential file: Name and format defined in ``i_*`` or ``info_*``

**Execution:**

.. code-block:: bash

   mpirun -np <num_cores> $(MuST_PATH)/bin/mst2 < i_file

Example input files are available under ``MST/sample/``.

.. note::

   Unless modified in the architecture file, the executable name is ``mst2``.

**Command Line Options:**

.. list-table::
   :header-rows: 1
   :widths: 8 30 45 17

   * - Short Option
     - Long Option(s)
     - Description
     - Default Value

   * - ``-wp``
     - ``--write-procs``, ``--write_procs``
     - Number of Writing Processes
     - ``0``

   * - ``-rp``
     - ``--read-procs``, ``--read_procs``
     - Number of Reading Processes
     - ``0``

   * - ``-app``
     - ``--atoms-per-proc``, ``--atoms_per_proc``
     - Number of Atoms Per Process
     - ``0``

   * - ``-cpu``
     - ``--cpu-only``, ``--cpu_only``
     - Run on CPU without Acceleration
     - ``-1``

   * - ``-abm``
     - ``--acc-bigmat``, ``--acc_bigmat``
     - Accelerate KKR Matrix Calculation
     - ``-1``

   * - ``-cbm``
     - ``--cpu-bigmat``, ``--cpu_bigmat``
     - KKR Matrix Calculation on CPU
     - ``-1``

   * - ``-agm``
     - ``--acc-gijmat``, ``--acc_gijmat``
     - Accelerate Gij Matrix Calculation
     - ``-1``

   * - ``-cgm``
     - ``--cpu-gijmat``, ``--cpu_gijmat``
     - Gij Matrix Calculation on CPU
     - ``-1``

   * - ``-tl``
     - ``--timing-lsms``, ``--timing_lsms``
     - Timing the LSMS Matrix Calculations
     - ``-1``

   * - ``-nbc``
     - ``--no-blocking-constr``, ``--no_blocking_constr``
     - Constructing the LSMS Matrix without using blocks
     - ``-1``

   * - ``-pb``
     - ``--print-blocking``, ``--print_blocking``
     - Print Blocking Details in LSMS Matrix Inverse
     - ``-1``

   * - ``-nbi``
     - ``--no-blocking-inverse``, ``--no_blocking_inverse``
     - The LSMS Matrix Inverse is Performed without Blocking
       (*reserved -- accepted but not yet implemented, see note below*)
     - ``-1``

   * - ``-re``
     - ``--report-energy``, ``--report_energy``
     - Report GPU Energy Consuption
     - ``-1``

   * - ``-ptm``
     - ``--print-tau``, ``--print_tau``
     - Print the Tau Matrix
     - ``-1``

   * - ``-em``
     - ``--emul-scheme``, ``--emul_scheme``
     - Emulation Scheme
     - ``0``

   * - ``-emp``
     - ``--emul-param``, ``--emul_param``
     - Emulation Precision Parameters
     - ``12,14,16``

   * - ``-cz``
     - ``--contour-zone``, ``--contour_zone``
     - Energy Contour Zone Parameters
     - ``0.000,0.010``

.. note::

   The switch-style options above default to ``-1``, which is the sentinel for
   *not supplied*: ``getCmdLineOption`` returns ``1`` for such a key and ``0``
   once the switch is given, so a test of the form
   ``getCmdLineOption('Run on CPU without Acceleration') == 0`` is true exactly
   when the user passed ``-cpu``/``--cpu-only``. An option whose default is not
   ``-1`` therefore reads as *always supplied*.

   ``-wp``, ``-rp`` and ``-app`` take an integer argument; ``-em`` takes an
   integer scheme selector (``0`` none, ``1`` Ozaki 1, ``2`` Ozaki 2); ``-emp``
   and ``-cz`` take comma-separated lists. The remaining options are switches
   that take no argument.

   The acceleration options (``-cpu``, ``-abm``, ``-cbm``, ``-agm``, ``-cgm``,
   ``-re``) have an effect only in a build configured with GPU support.
   ``-cpu``/``--cpu-only`` takes precedence over the others and over the
   ``MUST_*_GPU`` environment variables, disabling every accelerated path:
   the KKR matrix (``ClusterMatrixModule``), the radial density interpolation
   (``DensityOnGridModule``), the pseudo-potential back-projection
   (``PseudoPotBackProjModule``) and the single-site radial march
   (``SSMarchModule``).

.. warning::

   ``-nbi``/``--no-blocking-inverse`` is **reserved and currently has no
   effect**. The option is parsed and stored, but no source file queries its key
   (``'The LSMS Matrix Inverse is Performed without Blocking'``), so the CPU
   inverse always takes the blocked path through ``invertMatrixBlock`` in
   ``MST/src/MatrixBlockInversionModule.F90``. Supplying the flag is harmless
   but does nothing. This is the counterpart of ``-nbc``, which *is* wired in
   and does select a non-blocked construction of the KKR matrix.

---

2. ``lsms``
===========

Main source: ``lsms/src/Main/lsms.cpp``. Performs *ab initio*, linear-scaling
electronic structure calculations for 3-D systems.

**Main features:**

- Linear scaling based on LSMS method
- Muffin-tin potential
- Non-relativistic or scalar-relativistic
- Non-spin polarized, spin-polarized, or spin-canted
- LDA or GGA exchange-correlation

**Execution:**

.. code-block:: bash

   mpirun -np <num_cores> $(MuST_PATH)/bin/lsms < i_file

Example inputs: ``lsms/Test/``

---

3. ``wl-lsms``
==============

Main source: ``lsms/src/Main/wl_lsms.cpp``. Performs Wang–Landau Monte Carlo
simulations using LSMS-generated energy data.

**Main features:**

- Wang–Landau Monte Carlo method
- Drives linear-scaling *ab initio* energy calculations

**Execution:**

.. code-block:: bash

   mpirun -np <num_cores> $(MuST_PATH)/bin/wl-lsms < i_file

Example inputs: ``lsms/Test/``

---

4. ``genap``
============

Utility code (``MST/util/generateAtomPosition.F90``) for generating unit cell
samples of ordered compounds or disordered alloys.

**Execution:**

.. code-block:: bash

   $(MuST_PATH)/bin/genap

Input is provided interactively via the terminal.

---

5. ``measureVoronoi``
=====================

Utility (``MST/util/measureVoronoi.F90``) for computing geometric properties
of Voronoi polyhedra for atoms in a unit cell.

**Execution:**

.. code-block:: bash

   mpirun -np <num_cores> $(MuST_PATH)/bin/measureVoronoi < i_file

.. note::

   The input file is the same as used for ``mst2``.

---

6. ``murn``
===========

Utility (``MST/util/murn_new.F90``) for determining ground-state properties
(lattice constant, volume, bulk modulus) from energy vs. volume data.

**Execution:**

.. code-block:: bash

   $(MuST_PATH)/bin/murn < input_file

Example input: ``MST/sample/Co/a0/inp_murn``

---

7. ``newa``
===========

Utility (``MST/util/newa.F``) for generating initial atomic potentials.

**Input file:**

- ``*_a_in``: Specifies atom type, spin, output names, etc.

**Output files:**

- ``*_a_out``: Standard output file
- ``*_a_pot``: Generated potential file

**Execution:**

.. code-block:: bash

   $(MuST_PATH)/bin/newa < input_file

Example: ``MST/sample/Mg/Atom/Mg_a_in``

---

8. ``newss``
============

Utility (``MST/util/newss.F``) for generating starting potentials for
KKR/KKR-CPA/LSMS calculations.

**Input files:**

- ``*_ss_in``: Defines lattice, structure, and potential parameters
- ``*_a_pot``: Generated using ``newa``

**Output files:**

- ``*_ss_out``: Main output data
- ``*_ss_k``: Energy and RMS per SCF iteration
- ``*_ss_pot``: Generated starting potential

**Execution:**

.. code-block:: bash

   $(MuST_PATH)/bin/newss < input_file

Example: ``MST/sample/Mg/Atom/Mg_ss_in``
