#include "checkAntiFerromagneticStatus.hpp"


void checkIfSpinHasFlipped(LSMSSystemParameters &lsms, AtomData &a)
{
  if (lsms.n_spin_cant == 1 && lsms.n_spin_pola == 2)
  {
    Real spinFlipDirection = (a.xvalws[0] - a.xvalws[1]) * (a.xvalwsNew[0] - a.xvalwsNew[1]);

    if (spinFlipDirection > 0.0) 
      a.spinFlipped = false;
    else
      a.spinFlipped = true;
  }

}


void swapCoreStateEnergies(AtomData &a)
{

  Real tmp;
  for (int coreEnergyLevel=0; coreEnergyLevel<a.numc; coreEnergyLevel++)
  {
    tmp = a.ec(coreEnergyLevel, 0);
    a.ec(coreEnergyLevel, 0) = a.ec(coreEnergyLevel, 1);
    a.ec(coreEnergyLevel, 1) = tmp;
  }

}
