Equation of State using Birch-Murnaghan EOS
===========================================

Dataset Download
-------------------

The complete dataset (inputs, scripts, and example workflow) is available here:

`Download FeNi3 EOS tutorial files <_static/FeNi3_EV.tar.gz>`_

The archive contains all files required to reproduce the calculation for the FeNi\ :sub:`3` alloy in the L1\ :sub:`2` phase.

Overview
-------------

This tutorial describes how to compute the equation of state (EOS) by fitting total energy vs. volume data using the Birch-Murnaghan EOS.

The goal is to determine:

- Equilibrium volume :math:`V_0`
- Minimum energy :math:`E_0`
- Bulk modulus :math:`B_0`
- Pressure derivative :math:`B_0'`

Required Files
-------------------

.. list-table::
   :header-rows: 1

   * - File
     - Description
     - Usage
   * - i_mst
     - Main input file
     - Controls the simulation
   * - Evec_input.dat
     - Moment direction data
     - 
   * - Fe_mt_v
     - Fe pseudopotential
     - 
   * - Ni_mt_v
     - Ni pseudopotential
     - 
   * - position.dat
     - Structure file
     - Atomic positions
   * - volume.dat
     - Volume scaling factors
     - Scales lattice constant
   * - job.sh
     - Job submission script
     - Adjust for cluster and MPI
   * - job-submit
     - Batch execution script
     - Run multiple volume calculations
   * - extract.sh
     - Data extraction script
     - Generates energy-volume data
   * - eos.py
     - EOS fitting script
     - Computes equation of state

Step 1: Extract Dataset
---------------------------

.. code-block:: bash

   tar -xzf FeNi3_EV.tar.gz

Step 2: Prepare Volume Scaling
---------------------------------

Edit ``volume.dat``:

.. code-block:: text

   0.94
   0.96
   0.98
   1.00
   1.02
   1.04
   1.06

The lattice constant scales as:

:math:`a = a_0 \times s`

Step 3: Configure Input File
-------------------------------

Edit ``i_mst``:

- Set pseudopotentials
- Ensure SCF calculation
- Use tight convergence

Step 4: Configure Job Script
--------------------------------

.. code-block:: bash

   mpirun -np 32 ./mst.x > output.log

MPI tasks should satisfy:

:math:`N_{\text{MPI}} = k \times N_{\text{atoms}}`

Step 5: Run Calculations
----------------------------

.. code-block:: bash

   ./job-submit

Step 6: Extract Energy-Volume Data
-------------------------------------

.. code-block:: bash

   ./extract.sh

Output:

.. code-block:: text

   Volume, Energy
   V_1, E_1
   V_2, E_2

Step 7: Fit Birch-Murnaghan EOS
---------------------------------

.. code-block:: bash

   python eos.py

The EOS is given by:

:math:`E(V) = E_0 + \frac{9 V_0 B_0}{16} \left[ \left( \left( \frac{V_0}{V} \right)^{2/3} - 1 \right)^3 B_0' + \left( \left( \frac{V_0}{V} \right)^{2/3} - 1 \right)^2 \left( 6 - 4 \left( \frac{V_0}{V} \right)^{2/3} \right) \right]`

.. figure:: _static/FeNi3_EV.png
   :align: center
   :width: 60%

   Energy vs. volume curve for FeNi\ :sub:`3` in the L1\ :sub:`2` phase. The discrete points correspond to calculated total energies at different volumes, while the smooth curve represents the Birch-Murnaghan EOS fit. The minimum of the curve gives the equilibrium volume :math:`V_0`.

Outputs:

- :math:`V_0`
- :math:`E_0`
- :math:`B_0`
- :math:`B_0'`

Workflow Summary
-------------------

.. code-block:: text

   FeNi3_EV.tar.gz -> extract -> run jobs -> extract.sh -> result.csv -> eos.py
