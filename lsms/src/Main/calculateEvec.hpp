#ifndef CALCULATEEVEC_HPP
#define CALCULATEEVEC_HPP
 
#include <cmath>
#include "Main/SystemParameters.hpp"
//#include "SingleSite/AtomData.hpp"

void calculateEvec(LSMSSystemParameters &lsms, LocalTypeInfo &local);

void mixEvec(LSMSSystemParameters &lsms, LocalTypeInfo &local, Real alpev);

#endif
