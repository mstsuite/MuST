#ifndef LSMS_POTENTIALIO_H
#define LSMS_POTENTIALIO_H

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

int loadPotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);
int writePotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);
void initialAtomSetup(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);

int loadAlloyBank(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyMixingDesc &alloyDesc, AlloyAtomBank &alloyBank);

#endif
