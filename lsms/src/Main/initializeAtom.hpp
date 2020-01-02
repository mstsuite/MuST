#ifndef LSMS_INITIALIZEATOMS_H
#define LSMS_INITIALIZEATOMS_H
#include <stdio.h>
#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "SingleSite/AtomData.hpp"

void initializeAtom(AtomData &a);
int initializeNewPotentials(LSMSCommunication &comm,LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);
int initializeNewAlloyBank(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyMixingDesc &alloyDesc, AlloyAtomBank &alloyBank);

#endif
