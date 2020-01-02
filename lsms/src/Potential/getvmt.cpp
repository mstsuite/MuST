#include "getvmt.hpp"

// getvmt.f in LSMS 1
// vshift still needs to be calculated somewhere else first! (now local and =0.0)

void getvmt(LSMSSystemParameters lsms, AtomData &atom, CrystalParameters &crystal, Real qsub[], int &mytype, Real &vmt, Real &vmt1, Real &u0)
{
/*
  ================================================================
  Calculate muffin-tin zero potential and its contribution to the
  Coulomb energy..................................................
  note: lattice constant factor is included in madmat.............
  ================================================================
*/
  vmt = 0.0;       // contribution to site-independent constant MT potential
                   // (yet to be averaged and divided by omegint)
  vmt1 = 0.0;      // site-dependent constant MT potential
                   // (in the case of ASA-MT, vmt is the shift inside rins,
                   // and vmt1 is the shift for the whole ASA sphere.)
  u0 = 0.0;        // contribution to the total energy

  for(int i=0; i<crystal.num_types; i++)
  {
    vmt1 += atom.madelungMatrix[i] * qsub[i];
    u0 += atom.madelungMatrix[i] * qsub[i] * qsub[mytype];
  }

  vmt1 *= 2.0;
/*
  =============================================================
  Calculate vmt, the site-independent electro-static
  contribution to the muffin-tin zero potential
  =============================================================
*/
  Real dq_mt = atom.qtotmt - atom.ztotss;
  Real rSphere;

  switch (lsms.mtasa)
  {
    case 0:                               // Muffin-Tin case
      //rSphere = atom.rmt;
      rSphere = atom.rInscribed;
      break;
    default:                              // ASA or ASA-MT cases
      rSphere = atom.rws;
  }

  switch (lsms.mtasa)
  {
    case 1:                               // ASA case
      vmt = 2.0 * dq_mt / rSphere;
      break;
    
    default:                              // Muffin-Tin case
      Real surfaceAreaMT = 4.0 * M_PI * rSphere * rSphere;
      vmt = vmt1 * atom.omegaMT + surfaceAreaMT * (atom.omegaMT * atom.rhoInt / 5.0 + qsub[mytype]);
      u0 += atom.rhoInt * atom.omegaMT * (-6.0/5.0 * atom.rhoInt * atom.omegaMT + 3.0 * qsub[mytype]) / rSphere;
      vmt1 += surfaceAreaMT * atom.rhoInt;
  
      if (lsms.mtasa >= 2)                // MT-ASA case
      {
        // additional contribution to vmt1 from variation of qCorrection
        // Delta Q_i, the charge due to shape correction
        Real qCorrection = atom.qtotws - atom.qtotmt - atom.rhoInt*(atom.omegaWS - atom.omegaMT);
        vmt1 -= 2.0 * qCorrection / rSphere;
        u0 += qCorrection * (2.0 * dq_mt + qCorrection) / rSphere;
      }
  }

  Real vshift = 0.0;
  vmt1 -= vshift;

return;

}
