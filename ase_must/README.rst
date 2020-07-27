======================
ASE Interface for MuST
======================

Introduction
------------
`Atomic Simulation Environment (ASE)`_ is a Python library which makes it easy to set up, steer, and analyze atomistic simulations.
It interfaces with a number of ab initio calculators such as VASP, Quantum Espresso, Gaussian etc. An ASE interface for the MuST package
has been included locally with this code. To use this interface, you will need the ASE package installed in your system.

.. _Atomic Simulation Environment (ASE): https://wiki.fysik.dtu.dk/ase/about.html

Setting Up
----------

To use this interface, build the MuST package, `install ASE`_ and then add the path to ``MuST/ase_must`` folder of your MuST installation to your :code:`PYTHONPATH` variable.

.. _install ASE: https://wiki.fysik.dtu.dk/ase/install.html

Getting Started
---------------
It is highly recommended to get familiar with basics of ASE before using this interface. In particular, you should go through
ASE's guide for its `Atoms object`_.

.. _Atoms object: https://wiki.fysik.dtu.dk/ase/ase/atoms.html

The following guide with walk you through generating starting potentials and setting up a KKR and KKR-CPA calculation using this ASE interface.
Please note that at this point, only potential energy has been implemented in this interface. You cannot calculate forces yet.

Generating Starting Potentials
------------------------------
Every calculation requires starting potentials. You can generate single site starting potential for each element
present in the :code:`Atoms` object using
:code:`generate_starting_potentials` function.

::

    from ase_must.must import generate_starting_potentials

    atoms = bulk('Al', a=4.05)
    generate_starting_potentials(atoms, crystal_type=1, a=4.05) # crystal_type=1 for FCC, 2 for BCC

You will see the output files generated in your working directory.

KKR Calculation
---------------
Use ``El_ss_pot`` (El is the element),
as starting potential file name in :code:`MuST` calculator to perform a KKR calculation::

    from ase_must.must import MuST

    calc = MuST(default_in_pot='Al_ss_pot', # Starting potential file name
           default_out_pot='Al',        # Output potential file name
           pot_in_form=0,               # Starting potential file format
           method=2,                    # Method of calculation (2=KKR)
           potential_type=0,            # Potential type (0=Muffin-tin)
           xc=4,                        # XC functional (4=PBE)
           spin=1,                      # Number of spins
           kpts=(12,12,12),             # K-points grid
           bzsym=1)                     # Symmetrize BZ Integration

    # Attach calculator to the atoms object
    atoms.set_calculator(calc)

    # Trigger KKR calculation using .get_potential_energy()
    energy = atoms.get_potential_energy()
    print(energy, 'eV')

KKR-CPA Calculation
-------------------
To perform a CPA calculation, you need to pass the CPA formulation using ``info`` attribute of :code:`Atoms`
object.

* CPA sites should be defined as a dictionary: ``{'CPA': [list_of_sites]}``
* Each element in the ``[list_of_sites]`` is a dictionary ``{'index': int, 'Symbol1': fraction, 'Symbol2': fraction ...}``
* Here, ``index`` should match the index from :code:`Atoms` object. That's where it's position is taken
  from. Length of this list should match the number of atoms in :code:`Atoms` object.
* The cell dimensions and positions are still taken from the :code:`Atoms` object. Use ``method=3``
  in :code:`MuST <MuST.ase_must.must.MuST>` calculator to perform a KKR-CPA calculation

::

    atoms = bulk('Au', 'fcc', a=4.2)

    # Add CPA sites
    atoms.info = {'CPA': [{'index': 0, 'Au': 0.5, 'Pd': 0.5}]}

    # Generate starting potentials
    generate_starting_potentials(atoms, crystal_type=1, a=4.2, cpa=True)

    # Create a MuST calculator object
    calc = MuST(default_in_pot='Au_ss_pot Pd_ss_pot',
           default_out_pot='AuPd',
           pot_in_form=0,
           nscf=90,
           method=3,                    # Method of calculation (3=KKR-CPA)
           potential_type=0,
           xc=4,
           spin=1,
           kpts=(12,12,12),
           bzsym=1)

    # Attach calculator to the atoms object
    atoms.set_calculator(calc)

    # Trigger KKR-CPA calculation using .get_potential_energy()
    energy = atoms.get_potential_energy()
    print(energy, 'eV')

If your atoms object has multiple atoms, you can define multiple CPA sites in the same way::

    atoms.info = {'CPA': [{'index': 0, 'Al': 0.5, 'Ni': 0.5}, {'index': 1, 'Al': 0.5, 'Ni': 0.5}]}

Parallelization
_______________
MuST calculator uses single core by default. You can pass the :code:`ntasks` keyword while creating the MuST object to
define the number of cores you want to use for your calculation. Please note that it will run :code:`mpirun -np <ntasks> mst2 < i_new`
command when you do so. ``i_new`` is the filename used for the input file that :code:`MuST` object creates.

Parameter Keywords
------------------
Keyword related to each :code:`MuST`
parameter is mentioned in this table:

==========================      =====================================
keyword                          parameter
==========================      =====================================
``default_in_ot``               Default Potential Input File Name
``in_pot``                      Potential Input File Name
``default_out_pt``              Default Potential Output File Name
``nscf``                        No. Iterations (> 0)
``method``                      Method of SCF Calculation
``out_to_scr``                  Output to Screen (y/n)
``out_level``                   Output Level (>= -1)
``out_proc_id``                 Output Proc. ID (>= -1)
``out_atom_id``                 Output Atom ID (>= -1)
``generate_movie``              Generate System Movie
``stop_rout_name``              Stop-at Routine Name
``write_pot_niter``             No. Iter for Each Pot. Write
``movie_niter``                 No. Iter for Each Movie
``calc_harris_energy``          Calc. Harris Energy (H.E.)
``ngauss_r``                    No. Gauss Pts. along r
``ngauss_theta``                No. Gauss Pts. along theta
``vband_bot_est``               Valence Band Bottom Est.
``temperature``                 Temperature Parameter (K)
``dos_id``                      DOS Run ID
``uniform_grid``                Uniform Grid Parameters
``visual_grid_type``            Visual Grid Type (0<D<4)
``grid_scale``                  Grid Scale
``grid_origin``                 Origin Grid Vector
``grid_1``                      Grid Vector 1
``grid_2``                      Grid Vector 2
``grid_3``                      Grid Vector 3
``grid_pts``                    Grid Points
``e_density_out_id``            Output Electron Density ID (>= -1)
``density_format``              Output Density Format
``etol``                        Energy (Ryd) Tol (> 0)
``ptol``                        Potential Tol (> 0)
``fetol``                       Fermi Energy Tol (> 0)
``slu_tol``                     SuperLU Tol (> 0)
``ktol``                        K-space Check Tol (> 0)
``rms_tol``                     Other RMS Tol (> 0)
``val_e_rel``                   Val. Electron Rel (>= 0)
``core_e_rel``                  Core Electron Rel (>= 0)
``add_electrons``               Additional Electrons
``charge_sym``                  Charge Symmetry (>=0)
``ss_solver``                   Single Site Solver (>= 0)
``ss_method``                   Single Site Solution Method (>=-1)
``irreg_sols``                  Irregular Solutions (>=0)
``pole_step``                   Pole Search Step (>0.0)
``sol_lmax_cutoff``             Solutions Lmax Cutoff
``compute_phase_shifts``        Compute Phase Shifts (>=0)
``lmax_solver``                 SS Lmax Potential Solver
``potential_type``              Potential Type (>= 0)
``xc``                          Exch-Corr. LDA Type (>= 0)
``lda_improve_scheme``          LDA Improvement Scheme
``lda_file_name``               LDA+U Parameter File Name
``moment_dir_file``             Moment Direction File Name
``spin``                        Spin Index Param (>= 1)
``int_espin``                   Interstitial Electron Spin
``canted_torque_coef``          Canted Moment Torque Coef.
``calc_j_ij``                   Calculate J_ij (y/n)
``read_mesh``                   Read E-mesh from emeshs.inp
``contour_type``                Contour Type (>= 0)
``n_contours``                  Number of Contours (> 0)
``egrid_type``                  Energy Grid Type (>= 0)
``n_egrids``                    No. Energy Grids
``extra_energy_pts``            No. Extra Energy Points
``offset_energy_pt``            Offset Energy Point
``erbot``                       Real Axis Bottom erbot
``ertop``                       Real Axis Top ertop
``eibot``                       Imag Axis Bottom eibot
``eitop``                       Imag Axis Top eitop
``iterate_fermi_energy``        Iterate Fermi energy
``real_axis_method``            SS Real Axis Int. Method
``real_axis_points``            SS Real Axis Int. Points
``t_inversion``                 T-matrix inversion (>= 0)
``m_inversion``                 M-matrix inversion (>= 0)
``n_time_steps``                No. Spin-dynamics Time Steps (>= 0)
``time_step``                   Spin-dynamics Time Step
``mix_quantity``                Mixing quantity type
``mix_algo``                    Mixing algorithm
``lloyd_correction``            Lloyd correction
``lloyd_mode``                  Lloyd mode
``k_solver``                    K-space Solver Method
``read_kmesh``                  Read K-mesh from kmeshs.inp
``k_scheme``                    Scheme to Generate K (>=0)
``n_kmesh``                     No. K Meshs in IBZ (> 0)
``kpts``                        Kx, Ky, Kz Division (> 0)
``bzsym``                       Symmetrize BZ Integration
``large_sphere_radius``         Large sphere radius (a.u.)
``pot_in_form``                 Default Potential Input File Form
``pot_out_form``                Default Potential Output File Form
``moment_direction``            Default Moment Direction
``constrain_field``             Default Constrain Field
``lmax_T``                      Default Lmax-T matrix
``lmax_wave_func``              Default Lmax-Wave Func
``lmax_pot``                    Default Lmax-Potential
``lmax_trunc_pot``              Default Lmax-Trunc Pot
``lmax_charge_den``             Default Lmax-Charge Den
``lmax_step_func``              Default Lmax-Step Func
``liz_neighbors``               Default LIZ # Neighbors
``liz_nn_shells``               Default LIZ # NN Shells
``liz_shell_lmax``              Default LIZ Shell Lmax
``liz_cutoff``                  Default LIZ Cutoff Radius
``rho_mix_param``               Default Rho  Mix Param.
``pot_mix_param``               Default Pot  Mix Param.
``mom_mix_param``               Default Mom  Mix Param.
``chg_mix_param``               Default Chg  Mix Param.
``evec_mix_param``              Default Evec Mix Param.
``max_core_radius``             Default Maximum Core Radius
``max_mt_radius``               Default Maximum Muffin-tin Radius
``ndivin``                      Default No. Rad Points ndivin
``ndivout``                     Default No. Rad Points ndivout
``nmult``                       Default Integer Factor nmult
``pseudo_charge_radius``        Default Pseudo Charge Radius
``screen_pot``                  Default Screen Pot.
``lmax_screen``                 Default Lmax-Screen
``rcut_screen``                 Default Rcut-Screen
``local_sic``                   Local SIC
``mix_param``                   Default Mixing Parameter
``frozen_core_calc``            Frozen-Core Calculation
``frozen_core_file``            Frozen-Core File Name
``em_iter``                     Maximum Effective Medium Iterations
``em_scheme``                   Effective Medium Mixing Scheme
``em_mix_param``                Effective Medium Mixing Parameters
``em_eswitch``                  Effective Medium Mixing eSwitch Value
``em_tmatrix_tol``              Effective Medium T-matrix Tol (> 0)
``core_radius``                 Default Core Radius
``mt_radius``                   Default Muffin-tin Radius
``radical_plane_ratio``         Default Radical Plane Ratio
==========================      =====================================

Best place to look for default value of each parameter is the DefaultParameters.h_ file in MuST source code.

.. _DefaultParameters.h: https://github.com/mstsuite/MuST/blob/master/MST/src/DefaultParameters.h
