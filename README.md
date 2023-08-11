<img src=MuST_logo.png></img> <br>

<b><i>MuST</i></b> (<strong><i>Mu</i></strong>ltiple <b><i>S</i></b>cattering <b><i>T</i></b>heory) is an ab initio electronic structure calculation software suite, with petascale and beyond computing capability, for the first principles study of quantum phenomena in disordered materials. <br>
It is capable of performing
<ul>
  <li> KKR for ordered structures
  <li> KKR-CPA for random structures (with/without short range chemical order)
  <li> LSMS calculations for large systems
  <li> Kubo-Greenwood method for residual resistivity calculation
  <li> ...and many more upcoming features!
</ul>
This repository is actively developed and maintained - please check for regular updates! <br><br> 

<p float="left">
  <a href='https://must.readthedocs.io/en/latest/?badge=latest'>
    <img src='https://readthedocs.org/projects/must/badge/?version=latest' alt='Documentation Status' />
</a>
<img src=https://img.shields.io/github/license/mstsuite/MuST> 
<a href="https://github.com/mstsuite/MuST/wiki">
  <img src=https://img.shields.io/badge/-MuST%20Wiki-blue> </a>
<a href="https://www.youtube.com/channel/UCvlVeAb_m4kvBa-3_q43TIQ">
  <img src=https://img.shields.io/badge/-MuST%20Youtube%20Channel-red></a>
 </p>

## User Guide
All the relevant information and instructions are provided in the <a href="https://must.readthedocs.io/en/latest/index.html#">documentation</a>

For using MuST without compilations, please check out the Docker image <a href="https://hub.docker.com/r/liangstein/must">here</a>

## Available Scientific Packages
<ul>
  <li> <b>MST</b>: Perform KKR, LSMS, single-site and Cluster Averaged KKR-CPA.
  <li> <b>lsms</b>: Peform LSMS and Wang-Landau LSMS. This package is built for extreme performance on petascale/exascale systems.
  <li> <b> KUBO </b>: Perform first-principles electrical conductivity calculation.
</ul>

### User Support Folders
<ul>
<li> <b>Potentials</b> folder contains the starting potential for selected elements.
<li> <b>architecture</b> folder contains preset makefile parameters ("architecture files") for a wide variety of computer systems
<li> <b>docs</b> folder contains install instructions, license information, and users guide.
  <li> <b>external</b> folder contains external libraries required or optionally required by MuST, e.g., FFTW, Lua, P3DFFT, and LibXC.
  <li> <b>Tutorials</b> folder contains hands-on exercises and training materials.
<li> <b>ase_must</b> folder provides Atomic Simulation Environment (ASE) support for MuST.
</ul>

## Selected references
### KKR Method/Multiple Scattering Theory
* J. Korringa, _On the calculation of the energy of a Bloch wave in a metal_, Physica **13**, 392 (1947).

* W. Kohn and N. Rostoker, _Solution of the Schrodinger equation in periodic lattices with an application to metallic Lithium_, Phys. Rev. **94**, 1111 (1954).

* J. S. Faulkner and G. M. Stocks, _Calculating properties with the coherent-potential approximation_, Phys. Rev. B **21**, 3222 (1980).

* A. Gonis, _Green functions for ordered and disordered systems_, North-Holland Amsterdam, 1992

* A. Gonis and W. H. Butler, _Multiple Scattering in Solids_, (Graduate Texts in Contemporary Physics), Springer 1999.

* J. Zabloudil, R. Hammerling, L. Szunyogh, and P. Weinberger, _Electron Scattering in Solid Matter: A Theoretical and Computational Treatise_, (Springer Series in Solid-State Sciences), Springer 2004.

* H. Ebert, D. Kodderitzsch and J. Minar, _Calculating condensed matter properties using the KKR-Green's function method - recent developments and applications_, Rep. Prog. Phys. **74**, 096501 (2011).

* J.S. Faulkner, G.M. Stocks, and Y. Wang, _Multiple Scattering Theory: Electronic Structure of Solids_, IOP Publishing Ltd. 2019.

### KKR-CPA Method
* P. Soven, _Coherent-Potential Model of Substitutional Disordered Alloys_, Phys. Rev. **156**, 809 (1967).

* B. Velicky, S. Kirkpatrick, and H. Ehrenreich, _Single-Site Approximations in the Electronic Theory of Simple Binary Alloys_, Phys. Rev. **175**, 747 (1968).

* B. Gyorffy, _Coherent-Potential Approximation for a Nonoverlapping-Muffin-Tin-Potential Model of Random Substitutional Alloys_, Phys. Rev. B **5**, 2382 (1972).

* G. Stocks, W. Temmerman, and B. Gyorffy, _Complete Solution of the Korringa-Kohn-Rostoker Coherent-Potential-Approximation Equations: Cu-Ni Alloys_, Phys. Rev. Lett. **41**, 339 (1978).

* J. S. Faulkner and G. M. Stocks, _Calculating properties with the coherent-potential approximation_, Phys. Rev. B **21**, 3222 (1980).

* G. M. Stocks and H. Z. Winter, _Self-consistent-field-Korringa-Kohn-Rostoker-coherent-potential approximation for random alloys_, Z. Physik B-Condensed Matter **46**, 95 (1982).
* 

## Citation
If you publish results obtained using LSMS we ask that you cite the following publications:

* Y. Wang, G. M. Stocks, W. A. Shelton, D. M. C. Nicholson, W. M. Temmerman, and Z. Szotek. _Order-n multiple scattering approach to electronic structure calculations_. Phys. Rev. Lett. **75**, 2867 (1995).

and if the GPU accelerated version was used, please cite additionally:

* M. Eisenbach, J. Larkin, J. Lutjens, S. Rennich, and J. H. Rogers. _GPU acceleration of the locally selfconsistent multiple scattering code for first principles calculation of the ground state and statistical physics of materials_. Computer Physics Communications **211**, 2 (2017).

and for calculations using Monte-Carlo simulations:

* M. Eisenbach, C.-G. Zhou, D. M. C. Nicholson, G. Brown, J. Larkin, and T. C. Schulthess. _A Scalable Method for Ab Initio Computation of Free Energies in Nanoscale Systems_. Proceedings of the Conference on High Performance Computing Networking, Storage and Analysis, ACM, New York, 64 (2009)
