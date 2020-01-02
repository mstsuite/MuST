#ifndef LSMS_SHIFT_POTENTIALS
#define LSMS_SHIFT_POTENTIALS

class PotentialShifter {

public:

  bool vSpinShiftFlag {false};
  double minShift {0.0};
  double maxShift {0.0};

  void resetPotentials(LocalTypeInfo &local) {
    for (int i=0; i<local.num_local; i++) 
      vr0[i] = local.atom[i].vr;
  }

  void resize(int n) { vr0.resize(n); }

  void applyShifts(LocalTypeInfo &local)
  {
    for (int i=0; i<local.num_local; i++)
    {
      for (int ir=0; ir<local.atom[i].r_mesh.size(); ir++)
      {
//        Real deltaPotential=local.atom[i].vSpinShift*local.atom[i].r_mesh[ir]*
//                            (local.atom[i].rhotot(ir,0)-local.atom[i].rhotot(ir,1));
        Real deltaPotential = 0.5 * local.atom[i].vSpinShift * local.atom[i].r_mesh[ir];
        local.atom[i].vr(ir,0) = vr0[i](ir,0) - deltaPotential;
        local.atom[i].vr(ir,1) = vr0[i](ir,1) + deltaPotential;
      }
    }
  }

  void restorePotentials(LocalTypeInfo &local) {
    for (int i=0; i<local.num_local; i++)
      local.atom[i].vr=vr0[i];
  }               

  std::vector<Matrix<Real> > vr0; // the unshifted potential

};

#endif
