# MuST
MuST is an ab initio electronic structure calculation software package, with petascale and beyond computing capability, for the first principles study of quantum phenomena in disordered materials. It is capable of performing KKR, KKR-CPA, and LSMS calculations for ordered or disordered structures.

## A sketch of MuST package
In the top directory of MuST, there are following files and directories: 

DESCRIPTION: A brief description of MuST project

GUIDE: A users guide explaining how to use the package

INSTALL: An instruction for how to build and install MuST

LICENSE: License information

Makefile: the makfile for building MuST

architecture/: contains architecture files for some selected systems and their environments

bin/: contains exectables for running KKR, KKR-CPA, LSMS, and WL-LSMS calculations.

Documentation/: repository for storing instructions, license information, and users guide.

external/: contains external libraries required or optionally required by MuST, e.g., FFTW, Lua, P3DFFT, and LibXC.

lsms/: contains LSMS and WL-LSMS codes targeted to extreme performance on petascale/exascale systems.

MST/: contains MST packages targeted to physics development and capabilities, e.g. FP/MT KKR/LSMS/KKR-CPA, etc.

KUBO/: contains KUBO package for first-principles electrical conductivity calculation

Potentials/: contains the starting potential for selected elements.

Tutorials/: contains hands-on exercises and training materials.


## Selected references on KKR Method/Multiple Scattering Theory
* J. Korringa, _On the calculation of the energy of a Bloch wave in a metal_, Physica **13**, 392 (1947).

* W. Kohn and N. Rostoker, _Solution of the Schrodinger equation in periodic lattices with an application to metallic Lithium_, Phys. Rev. **94**, 1111 (1954).

* J. S. Faulkner and G. M. Stocks, _Calculating properties with the coherent-potential approximation_, Phys. Rev. B **21**, 3222 (1980).

* A. Gonis, _Green functions for ordered and disordered systems_, North-Holland Amsterdam, 1992

* A. Gonis and W. H. Butler, _Multiple Scattering in Solids_, (Graduate Texts in Contemporary Physics), Springer 1999.

* J. Zabloudil, R. Hammerling, L. Szunyogh, and P. Weinberger, _Electron Scattering in Solid Matter: A Theoretical and Computational Treatise_, (Springer Series in Solid-State Sciences), Springer 2004.

* H. Ebert, D. Kodderitzsch and J. Minar, _Calculating condensed matter properties using the KKR-Green's function method - recent developments and applications_, Rep. Prog. Phys. **74**, 096501 (2011).

* J.S. Faulkner, G.M. Stocks, and Y. Wang, _Multiple Scattering Theory: Electronic Structure of Solids_, IOP Publishing Ltd. 2019.

## Selected references on KKR-CPA Method
* P. Soven, _Coherent-Potential Model of Substitutional Disordered Alloys_, Phys. Rev. **156**, 809 (1967).

* B. Velicky, S. Kirkpatrick, and H. Ehrenreich, _Single-Site Approximations in the Electronic Theory of Simple Binary Alloys_, Phys. Rev. **175**, 747 (1968).

* B. Gyorffy, _Coherent-Potential Approximation for a Nonoverlapping-Muffin-Tin-Potential Model of Random Substitutional Alloys_, Phys. Rev. B **5**, 2382 (1972).

* G. Stocks, W. Temmerman, and B. Gyorffy, _Complete Solution of the Korringa-Kohn-Rostoker Coherent-Potential-Approximation Equations: Cu-Ni Alloys_, Phys. Rev. Lett. **41**, 339 (1978).

* J. S. Faulkner and G. M. Stocks, _Calculating properties with the coherent-potential approximation_, Phys. Rev. B **21**, 3222 (1980).

* G. M. Stocks and H. Z. Winter, _Self-consistent-field-Korringa-Kohn-Rostoker-coherent-potential approximation for random alloys_, Z. Physik B-Condensed Matter **46**, 95 (1982).

## Selected references on LSMS method
If you publish results obtained using LSMS we ask that you cite the following publications:

* Y. Wang, G. M. Stocks, W. A. Shelton, D. M. C. Nicholson, W. M. Temmerman, and Z. Szotek. _Order-n multiple scattering approach to electronic structure calculations_. Phys. Rev. Lett. **75**, 2867 (1995).

and if the GPU accelerated version was used, please cite additionally:

* M. Eisenbach, J. Larkin, J. Lutjens, S. Rennich, and J. H. Rogers. _GPU acceleration of the locally selfconsistent multiple scattering code for first principles calculation of the ground state and statistical physics of materials_. Computer Physics Communications **211**, 2 (2017).

and for calculations using Monte-Carlo simulations:

* M. Eisenbach, C.-G. Zhou, D. M. C. Nicholson, G. Brown, J. Larkin, and T. C. Schulthess. _A Scalable Method for Ab Initio Computation of Free Energies in Nanoscale Systems_. Proceedings of the Conference on High Performance Computing Networking, Storage and Analysis, ACM, New York, 64 (2009)
