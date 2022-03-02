Usage
=====

Executables under bin/:
=======================
1. mst:
The corresponding main source code is MST/src/mst2.F90. It performs ab
initio electronic structure calculations for 3-d structures.
Main features:
   * Linear scaling calculation based on LSMS method
   * Calculations based on KKR, KKR-CPA, or LSMS method
   * Muffin-tin potential or Full-potential
   * Non-relativistic, Scalar-relativistic, or Relativistic
   * Non-spin polarized, Spin-polarized, or Spin-canted
   * Special k-points method for BZ integration
   * LDA or GGA for the exchange-correlation potentials
Input files:
   * i_* file: The main input file which contains the controling parameters
     for running the SCF calculation and the parameters defining the system, and
     position data and potential file names
   * position data file: The actual file name is specified in the i_* file.
   * potential file(s): The actual file name and format
   * info_* file (obsolete): The actual file name is specifiled in the i_* file. It
     contains atom based controling parameters
   * kmeshs.inp (optional): This is an optional input file, only used for the testing
     purpose.
   * emeshs.inp (optional): This is an optional input file, only used for the testing
     purpose.
   * Evec_* (optional): This is an optional input file, only used in the spin-canted
     calculation case.
Output files:
   * o_n* file: This file contains the major output information. Note that
     it will not be created if the output is instructed (in the i_* input file)
     to be printed out to the screen.
   * k_n* file: This file contains a brief information of the total energy
     and the Fermi energy from each SCF iteration
   * new potential file: The actual file name and format is specified in the
     i_* or info_* file.
Execution:
   mpirun -np number_of_CPU_cores $(MuST_PATH)/bin/mst2 < i_file
Example input files for various structures can be found under MST/sample/.
Note: Unless otherwise changed name in archiecture file, the executable name is called
"mst2" by default.

2. lsms
The corresponding main source code is lsms/src/Main/lsms.cpp. It performs ab
initio, linear scaling, electronic structure calculations for 3-d structures.
Main features:
   * Linear scaling calculation based on LSMS method
   * Muffin-tin potential
   * Non-relativistic, or Scalar-relativistic
   * Non-spin polarized, Spin-polarized, or Spin-canted
   * LDA or GGA for the exchange-correlation potentials
Execution:
   mpirun -np number_of_CPU_cores $(MuST_PATH)/bin/lsms < i_file
Example input files for various structures can be found under lsms/Test/.

3. wl-lsms
The corresponding main source code is lsms/src/Main/wl_lsms.cpp. It performs
Wang-Landau Monte-Carlo simulation of random unit cell samples with energy data
obtained from LSMS electronic structure calculation.
Main features:
   * Wang-Landau Monte-Carlo simulation method
   * Driving linear scaling ab initio calculation of the energy data for the unit
     cell samples
Execution:
   mpirun -np number_of_CPU_cores $(MuST_PATH)/bin/wl-lsms < i_file
Example input files for various structures can be found under lsms/Test/.

4. genap:
A utility code (main: MST/util/generateAtomPosition.F90) for generating unit cell sample of
ordered compounds or disordered alloys (with random distribution or short-range order)
Execution:
   $(MuST_PATH)/bin/genap
The input data can be taken at the prompt on computer screen.

5. measureVoronoi
A utility code (main: MST/util/measureVoronoi.F90) for determining the geometric properties of
voronoi polyhedra generated for each atom in a unit cell sample.
Execution:
   mpirun -np number_of_CPU_cores $(MuST_PATH)/bin/measureVoronoi < i_file
Note, the input file, i_file, is the same as the one used for running bin/mst2.

6. murn
A utility code (main: MST/util/murn_new.F90) for determining the ground state properties
(lattice constant, unit cell volume, and bulk modulus) of a structure with given data for
 energy versus volume (or lattice constant).
Execution:
   $(MuST_PATH)/bin/murn < input_file
An example input file for murn, inp_murn, can be found under MST/sample/Co/a0/.

7. newa:
A utility code (main: MST/util/newa.F) for generating an initial atomic potential
Input file:
!   _a_in: input file specifying the atom type, spin information, output file name, etc
Output files:
   *_a_out: standard file, whose name is specified in the input file
   *_a_pot: potential file, whose name is specified in the input file
Execution:
   $(MuST_PATH)/bin/newa < input_file
An example input file for newa, Mg_a_in, for generating Mg atom potential can be found under
MST/sample/Mg/Atom/.

8. newss:
A utility code (main: MST/util/newss.F) for generating an initial potential for the KKR/KKR-CPA/LSMS
based electronic structure calculations.
Input files:
   *_ss_in: input file specifying lattice constant, crystal structure, potential file name, etc.
   *_a_pot: potential file generated from newa
Output files:
   *_ss_out: contains major ouput data
   *_ss_k:   contains a brief information of the total energy and the rms from each SCF iteration
   *_ss_pot: the starting potential for the KKR/LSMS calculation
Execution:
   $(MuST_PATH)/bin/newss < input_file
An example input file for newss, Mg_ss_in, for generating Mg starting potential for KKR/KKR-CPA/LSMS
can be found under MST/sample/Mg/Atom/.
