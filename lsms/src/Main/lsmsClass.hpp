// lsms class to encapsulate a version of LSMS_1.9 for use in gWL etc.

#ifndef LSMS_CLASS_H
#define LSMS_CLASS_H

#include <mpi.h>

#include <vector>

#include "Real.hpp"
#include "../Potential/PotentialShifter.hpp"
#include "mixing.hpp"

class LSMS {

public:

  //LSMS(MPI_Comm comm, const char * i_lsms, const char * out_prefix, PotentialShifter &potentialShifter);
  //LSMS(MPI_Comm comm, const char * i_lsms, const char * out_prefix);
  LSMS(MPI_Comm _comm, const char* i_lsms, const char* out_prefix, int my_group = 0);
  ~LSMS();

  int myWalkerID {0};

  int version() { return LSMS_version; }
  int numSpins() { return lsms.num_atoms; }
  bool vSpinShiftFlag() { return potentialShifter.vSpinShiftFlag; }
  double potentialMinShift() { return potentialShifter.minShift; }
  double potentialMaxShift() { return potentialShifter.maxShift; }

  void setEvec(std::vector<std::vector<Real>> &);
  void setEvec(Real *);
  void setEvecAndSpinPotentialShift(Real *);
  void setEvecAndSpinPotentialShift4(Real *);

  void getEvec(std::vector<std::vector<Real>> &);
  void getEvec(Real *);
  void getMag(std::vector<std::vector<Real>> &);
  void getMag(Real *);

  void setOccupancies(int*);
  void getOccupancies(int*);
  void getAlloyInfo(AlloyMixingDesc&, int**);

  Real oneStepEnergy(Real *eb);
  Real oneStepEnergy()
  {
    Real eband;
    return oneStepEnergy(&eband);
  }
  Real oneStepEnergy(std::vector<std::vector<Real>> &ev)
  {
    setEvec(ev);
    return oneStepEnergy();
  }

  Real multiStepEnergy();

  Real scfEnergy(Real *eb);
  Real scfEnergy()
  {
    Real eb;
    return scfEnergy(&eb);
  }

  void setEfTol(Real e) { efTol = e; }
  void setRmsTol(Real r) { rmsTolerance = r; }
  Real getEf(void) { return lsms.chempot; }
  void setEnergyTol(Real e) { energyTolerance = e; }
  void writePot(char *name);

  Real energyDifference;

  long energyLoopCount;

  void savePotentials(std::vector<Matrix<Real> > &vrs)
  {
    if(vrs.size()!=local.num_local) vrs.resize(local.num_local);
    for(int i=0; i<local.num_local; i++) vrs[i]=local.atom[i].vr;
  }

  void restorePotentials(std::vector<Matrix<Real> > &vrs)
  {
    for(int i=0; i<local.num_local; i++) local.atom[i].vr=vrs[i];
    mixing->prepare(comm,lsms,local.atom);
  }

  void replaceAtom(AtomData& currAtom, AtomData& newAtom);

private:

  char prefix[256];
  int LSMS_version;
  int max_num_local;
  LSMSSystemParameters lsms;
  LSMSCommunication comm;
  CrystalParameters crystal;
  LocalTypeInfo local;
  MixingParameters mix;

  PotentialShifter potentialShifter;

  Mixing *mixing;

  Real efTol;
  Real energyTolerance;
  Real rmsTolerance;

  // retain a bank of atoms (per species type) that can be used to load 
  // a guess (or frozen) potential when site occupancy changes.
  AlloyMixingDesc alloyDesc;
  AlloyAtomBank   alloyBank;

};

#endif

