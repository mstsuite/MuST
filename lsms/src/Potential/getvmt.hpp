#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "SingleSite/AtomData.hpp"
#include "Real.hpp"

//getvmt.f in LSMS 1
void getvmt(LSMSSystemParameters lsms, AtomData &atom, CrystalParameters &crystal, Real qsub[], int &mytype, Real &vmt, Real &vmt1, Real &u0);

