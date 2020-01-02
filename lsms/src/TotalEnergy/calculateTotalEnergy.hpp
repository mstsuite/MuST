#ifndef CALCULATE_TOTAL_ENERGY_HPP
#define CALCULATE_TOTAL_ENERGY_HPP

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Real.hpp"

void calculateTotalEnergy(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local, CrystalParameters &crystal); 

extern "C"
{
  void janake_(Real *vrold, Real *vrnew, Real *rhotot, Real *rhonew, Real *corden,
               Real *rr, Real *rins, Real *rmt, int *jmt, int *jws, int *komp, Real *atcon,
               Real *ztotss_dum, Real *atvol, Real *vx, Real *enxc, Real *evalsum, 
               Real *ecorv, Real *esemv, Real *etot, Real *press, Real *rspin,
               int *iprpts, int *iprint, char *istop, int istopl_len);
}

#endif
