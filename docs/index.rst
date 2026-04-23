Welcome to the MuST Documentation!
===================================

.. image:: _static/MuST_logo.png
   :align: center
   :width: 200px
   :height: 200px

**MuST** is a research project supported by the National Science Foundation (NSF) to build 
a public *ab initio* electronic structure calculation software package with petascale 
and beyond computing capability, for first-principles studies of quantum phenomena in 
disordered materials.  

The MuST package is now (as of January 1, 2020) free to download on GitHub 
(`https://github.com/mstsuite <https://github.com/mstsuite>`_) under a BSD 3-clause license.

MuST is developed based on full-potential multiple scattering theory (also known as the 
Korringa-Kohn-Rostoker method) using the Green’s function approach. It builds upon decades 
of development of research codes led by Malcolm Stocks and his postdocs and students in 
the Theory Group of the Metals and Ceramics Division (later Materials Science and 
Technology Division) at Oak Ridge National Laboratory.  

The original research codes include:

- Korringa-Kohn-Rostoker Coherent Potential Approximation (KKR-CPA): A highly efficient 
  *ab initio* method for the study of random alloys.
- Locally Self-consistent Multiple Scattering (LSMS) method: A linear-scaling *ab initio* 
  code capable of treating extremely large disordered systems using the largest parallel 
  supercomputers available.

It was originally suggested by Mark Jarrell, and demonstrated by model calculations, that 
strong disorder and localization effects can be studied within the LSMS formalism using 
cluster embedding in an effective medium with the Typical Medium Dynamical Cluster 
Approximation (TMDCA), enabling a scalable approach for first-principles studies of quantum 
materials.

The ultimate goal of the MuST project is to provide a computational framework for investigating 
quantum phase transitions and electron localization in the presence of disorder in real 
materials, and to enable computational studies of local chemical correlation effects on 
magnetic structure, phase stability, and mechanical properties of solid-state materials 
with complex structures.

The starting point of the MuST package is the integration of two research codes: 
``LSMS`` (formerly LSMS3) and ``MST`` (formerly MST2), both originally based on the
legacy LSMS-1 code developed in the mid-1990s at Oak Ridge National Laboratory.

LSMS
----

The ``LSMS`` code, maintained by Markus Eisenbach, is primarily written in C++. It 
consists of muffin-tin LSMS with an interface for Monte-Carlo simulation drivers. 
``LSMS`` is one of the baseline benchmark codes for DOE COREL systems and has been 
selected as a CAAR project for exascale computing on the Frontier system. It demonstrates 
nearly ideal linear scaling with 96% parallel efficiency on the Titan machine at ORNL.

MST
---

The ``MST`` code, maintained by Yang Wang, is mainly written in Fortran 90. It focuses 
on physics capabilities and serves as a platform for implementing and testing full-potential 
multiple scattering theory and its numerical algorithms. It includes LSMS, KKR, and KKR-CPA codes 
and supports:

1. Muffin-tin and full-potential calculations  
2. Non-relativistic, scalar-relativistic, and fully-relativistic approaches  
3. Non-spin-polarized, spin-polarized, and spin-canted *ab initio* electronic structure calculations  

KUBO
----

The ``KUBO`` code, maintained by Vishnu Raghuraman, is mainly written in Fortran 90. It 
implements the Kubo-Greenwood formula within the KKR-CPA framework and calculates the 
electrical conductivity of random alloys.

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   usage
   theory
   input
   tutorials
   contrib
   acknowledge
   cite
