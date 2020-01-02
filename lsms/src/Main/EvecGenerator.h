// mode: -*- c++ -*-

// Header for the evec generator class to be used in lsms.cc
// also example implementations:
// ConstantEvecGenerator
// RandomEvecGenerator

#ifndef EVEC_GENERATOR_H
#define EVEC_GENERATOR_H

#include <stdio.h>
#include <vector>
#include "random_evec.h"

typedef enum {SpinMove, OccupancyMove} MoveChoice_t;

/*
YingWai's note  (Dec 5, 13):
1. The split of the original generateEvec into "determineAcceptance, 
   updateHistogram and generateEvec" is only implemented in 
   WangLandau.h and WangLandau2d.h
2. generateEvec is renamed to updateHistogram in ExhautiveIsing.h.
   Nothing has been changed there otherwise.
   Main program (wl_lsms.cpp) was changed in a way where ExhaustiveIsing 
   was NOT considered. Use with caution!
*/

class EvecGenerator
{
 public:

  virtual double getDOSRatio(int instance, double energy)
  { return 0.0; }

// YingWai: do we really need the function overloading for determineAcceptance?   (Dec 4, 13)
  virtual bool determineAcceptance(int instance, double energy)
  { return false; }

  virtual bool determineAcceptance(int instance, double energy, double magnetization)
  {
    return determineAcceptance(instance, energy);
  }

  virtual bool updateHistogram(int instance, double *evecs, bool accepted)
  { return false; }
/*
  virtual bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool accepted)
  {
    return updateHistogram(instance, evecs, accepted);
  }
*/
  virtual bool updateHistogram(int instance, double *evecs, double *potentialShifts, bool accepted)
  {
    return updateHistogram(instance, evecs, accepted);
  }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, int check)
  { return false; }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, int check)
  { 
    return updateHistogramFromRE(instance, evecs, energy, check);
  }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, double *potentialShifts, int check)
  {
    return updateHistogramFromRE(instance, evecs, energy, check);
  }

  virtual bool updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, double *potentialShifts, int check)
  {
    return updateHistogramFromRE(instance, evecs, energy, check);
  }
/*
// YingWai: keep these function overloadings for updateHistogram for the moment,
//          but might want to get rid of them later when code is stable  (Dec 4, 13)

  virtual bool updateHistogram(int instance, double *evecs, double energy)
  { return false; }
  
  virtual bool updateHistogram(int instance, double *evecs, double energy, bool *accepted)
  {
    *accepted = true;
    return updateHistogram(instance, evecs, energy);
  }

  virtual bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool *accepted)
  {
    *accepted = true;
    return updateHistogram(instance, evecs, energy);
  }
*/

  virtual void generatePotentialShift(int instance, double *potentialShifts, bool accepted) {};

  virtual void generateEvec(int instance, double *evecs, bool accepted) {};

  virtual void initializeEvecAndPotentialShift(int instance, double *evecs, double *potentialShifts) {};

  virtual void initializeEvec(int instance, double *evecs)
  { 
    //generateEvec(instance, evecs, 0.0); 
    bool accepted {false};
    generateEvec(instance, evecs, accepted); 
  }

  virtual void generateUnsampledEvec(int instance, double *evecs, double energy) 
  { 
    //generateEvec(instance, evecs, energy);
    bool accepted {false};
    generateEvec(instance, evecs, accepted); 
  }

  virtual void startSampling(bool isspin = true, bool isocc = false) {;}

  virtual void writeState(const char *name) {;}

  void setVerbosity(int v) { verbosity = v; }

  int verbosity;

  // dummy functions to be used for occupancy variables
  // these routines were added to the Evec class because they share Monte Carlo logic
  virtual MoveChoice_t selectMoveType(bool isSpinSim, bool isOccSim)  { return SpinMove; }
  virtual void setAlloyClasses(const AlloyMixingDesc&, int* siteclass) {}
  virtual void generateOccupancies(int instance, int *occ, bool acceptedOcc) {}
  virtual void generateUnsampledOcc(int inst, int *occ) {}
  virtual void initializeOccupancies(int inst, int *occ) {}
  virtual void updateLog(int instance, double *evecs, int * occ, double energy, bool accepted, MoveChoice_t MoveChoice,
      bool isspin, bool isocc) {}
};


class ConstantEvecGenerator : public EvecGenerator
{
 public:

  ConstantEvecGenerator(size_t num_spins, double evec[3], int nw = 0)
  {
    printf("Constructor of ConstantEvecGenerator called.\n");
    n_spins = num_spins;
    setEvec(evec);
    n_walker = nw;
    if (nw > 0)
    {
      walker_step.resize(nw);
      for (int i=0; i<nw; i++) walker_step[i] = 0;
    }
  }

  void setEvec(double evec[3])
  {
    ev[0] = evec[0];
    ev[1] = evec[1];
    ev[2] = evec[2];
  }

  void generatePotentialShift(int instance, double *potentialShifts, bool accepted)
  {
    for (size_t i=0; i<n_spins; i++)
      potentialShifts[i] = 0.0;
  }

  void generateEvec(int instance, double *evecs, bool accepted)
  {
    for (size_t j=0; j<3*n_spins; j+=3)
    {
      evecs[j] = ev[0];
      evecs[j+1] = ev[1];
      evecs[j+2] = ev[2];
    }
    if (instance < n_walker) walker_step[instance]++;
    //return false;
  }

  void initializeEvecAndPotentialShift(int instance, double *evecs, double *potentialShifts)
  {
    initializeEvec(instance, evecs);
    generatePotentialShift(instance, potentialShifts, true);
  }

  void initializeEvec(int instance, double *evecs)
  {
    for(size_t j=0; j<3*n_spins; j+=3)
    {
      evecs[j] = ev[0];
      evecs[j+1] = ev[1];
      evecs[j+2] = ev[2];
    }
  }

 private:

  size_t n_spins;
  double ev[3];
  int n_walker;
  std::vector<int> walker_step;

};


class RandomEvecGenerator : public EvecGenerator
{
 public:

  RandomEvecGenerator(size_t num_spins)
  { n_spins = num_spins; }

  void generateEvec(int inst, double *evecs, double energy)
  {
    for(size_t j=0; j<3*n_spins; j+=3)
      random_evec(&evecs[j]);
    //return false;
  }

  void initializeEvec(int instance, double *evecs)
  { generateEvec(instance, evecs, 0.0); }

 private:
  size_t n_spins;

};

#endif
