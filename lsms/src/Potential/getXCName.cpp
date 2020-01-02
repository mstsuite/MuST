#include "Main/SystemParameters.hpp"
#include "Potential/libxcInterface.hpp"
#include "getXCName.hpp"
#include <string>

bool getXCName(LSMSSystemParameters &lsms, std::string &name)
{
  if(lsms.xcFunctional[0]==0) // built in functionals
  {
    switch(lsms.xcFunctional[1])
    {
      case 1: name="von Barth-Hedin (LSMS_1)"; return true;
      case 2: name="Vosko-Wilk-Nusair (LSMS_1)"; return true;
    }
    name="Illegal Exchange-Correlation Functional (built in)!"; return false;
  } else if(lsms.xcFunctional[0]==1) { // libxc functionals
#ifdef USE_LIBXC
    name="";
    
    if(lsms.libxcFunctional.numFunctionals==0) name="none";
    for(int i=0; i<lsms.libxcFunctional.numFunctionals; i++)
    {
      name.append(lsms.libxcFunctional.functional[i].info->name);
      if(i!=lsms.libxcFunctional.numFunctionals-1) name.append(" + ");
    }
    name.append(" (libxc)");
    return true;
#else
    name="Illegal Exchange-Correlation Functional (LSMS not linked to libXC)!"; return false;
#endif
  } else { // unknown functional!!
    name="Illegal Exchange-Correlation Functional!"; return false;
  }
  return false;
}
