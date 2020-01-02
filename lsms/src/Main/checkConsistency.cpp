#include "checkConsistency.hpp"

// Check that average integrated valence DOS = average valence
// (from mufind_c.f)
void checkIntegratedValenceDOS()
{

}



// Check that the charges are internally consistent
// (from lsms_main.f in LSMS 1.9, after getchg())
// ad-hoc now, needs amendment
void checkCharge(LSMSSystemParameters &lsms, AtomData &atom)
{
  if(lsms.global.iprint>=-1)
  {
    Real q = atom.xvalmt[0] + (lsms.n_spin_pola - 1) * atom.xvalmt[1];
    if (std::abs(q - atom.qvalmt) > 0.000005)
      printf("checkCharge [%d] :: Trouble! MT q = %16.8f qvalmt = %16.8f\n",lsms.commRank, q, atom.qvalmt); 
  }
}
