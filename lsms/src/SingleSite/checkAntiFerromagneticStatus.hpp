#ifndef CHECK_ANTIFERROMAGNETIC_STATUS_H
#define CHECK_ANTIFERROMAGNETIC_STATUS_H

#include "AtomData.hpp"
#include "Main/SystemParameters.hpp"


void checkIfSpinHasFlipped(LSMSSystemParameters &lsms, AtomData &a);

void swapCoreStateEnergies(AtomData &a);

#endif
