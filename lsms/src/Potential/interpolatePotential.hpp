#ifndef INTERPOLATE_POTENTIAL_HPP
#define INTERPOLATE_POTENTIAL_HPP

#include "Main/SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"

extern "C"
{
  void interp_(Real *r, Real *f, int *nr, Real *rs, Real *ps, Real *dps, int *deriv);
}

void interpolatePotential(LSMSSystemParameters &lsms, AtomData &atom);

#endif
