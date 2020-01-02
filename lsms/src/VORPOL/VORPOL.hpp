/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_VORPOL_H
#define LSMS_VORPOL_H

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "Array3d.hpp"

class VoronoiPolyhedra {
public:
/*
// VORPOL data [this is only used from routines imported from LSMS_1.9 .. VORPOL and replaces the common block /common_step/
      integer    ncrit
      real*8     gwwylm(ipngaussr*(iprcrit-1))
      real*8     grwylm(ipngaussr*(iprcrit-1))
      include 'atom_param.h'
      complex*16 wylm((2*iplmax+1)*(iplmax+1),ipngaussr,iprcrit-1)

c This common block is defined in subroutine setup_vorpol. It
c contains the step function and the r-mesh for interstitial integration
c ncrit = number of critical r points
c grwylm = array containing the r-mesh for interstitial integration
c gwwylm = array containing the weights for interstitial integration
c wylm = array containing the step function for each r-mesh point
 common /common_step/ wylm,gwwylm,grwylm,ncrit
*/
  int ncrit;
  Real rInscribedSphere; // inscribed sphere
  Real omegaInt; // interstitial volume
  Complex dipint[6];
  Array3d<Complex> wylm;
  Matrix<Real> gwwylm, grwylm;
};

extern "C"
{
  void setup_vorpol_(int *my_atom,int *num_atoms,
                     Real *atom_position_1,
                     Real *atom_position_2,Real *atom_position_3,
                     Real *system_bravais,
                     int *lmax,Real *clm,int *ngaussq,int *ngaussr,
                     Real *rmt,Real *omegint,Complex *dipint,Real *rad,
                     int *ipvp,int *ipnode, int *ipcorn, int *ipedge,int *iprcrit,
                     Real *gwwylm, Real*grwylm,
                     int *ncrit, Complex *wylm,
                     int *iprint,char *istop,int istopl_len);
}

#endif
