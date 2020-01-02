#ifndef LSMS_ENERGYCONTOURINTEGRATION_H
#define LSMS_ENERGYCONTOURINTEGRATION_H

#include "Complex.hpp"
#include "SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

// typedef enum {EnergyGridBox=1, EnergyGridGauss=2} EnergyGridType;

void energyContourIntegration(LSMSCommunication &comm,LSMSSystemParameters &lsms, LocalTypeInfo &local);

extern "C"
{
  void congauss_(double *ebot,double *etop,double *eibot,Complex *egrd,Complex *dele1,int *npts,int *nume,
                 double *pi,int * ipepts, int *iprint,char *istop, int istop_len);
}

#endif
