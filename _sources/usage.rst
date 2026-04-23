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

   * - Short Option
     - Long Option(s)
     - Description
     - Default Value

   * - -wp
     - --write-procs, --write_procs
     - Number of Writing Processes
     - 0

   * - -rp
     - --read-procs, --read_procs
     - Number of Reading Processes
     - 0

   * - -app
     - --atoms-per-proc, --atoms_per_proc
     - Number of Atoms Per Process
     - 0

   * - -cpu
     - --cpu-only, --cpu_only
     - Run on CPU without Acceleration
     - -1

   * - -abm
     - --acc-bigmat, --acc_bigmat
     - Accelerate KKR Matrix Calculation
     - -1

   * - -tl
     - --timing-lsms, --timing_lsms
     - Timing the LSMS Matrix Calculations
     - -1

   * - -pb
     - --print-blocking, --print_blocking
     - Print Blocking Details in LSMS Matrix Inverse
     - -1

   * - -nb
     - --no-blocking, --no_blocking
     - Perform LSMS Matrix Inverse without Blocking
     - -1

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
