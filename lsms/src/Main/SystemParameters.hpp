/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_SYSTEM_PARAM_H
#define LSMS_SYSTEM_PARAM_H
#include <stdio.h>
#include <string.h>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "LSMSMode.hpp"

#include "SingleSite/AtomData.hpp"

#include "Misc/Indices.hpp"


const int numFunctionalIndices=3;
#include "Potential/libxcInterface.hpp"

enum LinearSolver : int {defaultSolver = 0, oldSolver = 1, newSolver = 2};

class LSMSGlobals {
public:
  void setIstop(const char *c){strncpy(istop,c,32); for(int i=strlen(c); i<32; i++) istop[i]=' ';}
  bool checkIstop(const char *c){return (strncmp(istop,c,32)==0);}
  int iprpts,ipcore;
  int iprint;
  int print_node,default_iprint;
  char istop[32];
  LinearSolver linearSolver;

// for GPU only
  int GPUThreads;
};

class EnergyContourParameters {
public:
  int grid, npts;
  Real ebot,etop,eitop,eibot;
// Grouping of energies for single site solver
  int maxGroupSize;
  int groupSize() {int nume=npts; if(grid==0) nume=1; else if(grid==2) nume++; return std::min(nume,maxGroupSize);}
};

enum Relativity : int {none=0, scalar=1, full=2};

class LSMSSystemParameters {
public:
  char systemid[80];
  char title[80];
  LSMSMode lsmsMode;
  int rank;
  char potential_file_in[128];
  char potential_file_out[128];
  int pot_in_type,pot_out_type;
  char alloy_file_in[128];
  char alloy_file_out[128]; 
  int alloy_in_type,alloy_out_type;
  char infoEvecFileIn[128];
  char infoEvecFileOut[128];
  char localAtomDataFile[128];

  int mixing; // combines LSMS_1's mix_quant & mix_algor : -1 don't mix. mix_quant=mixing%4; mix_algor=mixing>>2;
              // mix_quant  0: charge, 1: potential
              // mix_algor  0: simple (linear) mixing; 1: broyden
              // charge, simple: 0; broyden: 4
              // potential, simple: 1; boyden: 5
  Real alphaDV; // mixing parameter for density or potential
  Real rmsTolerance; // rms Convergence criterion
  int num_atoms;
  int nspin;
  Relativity relativity;
  int nrelc,nrelv;
  int n_spin_cant;
  int n_spin_pola;
  int nscf;
  int writeSteps;
  int mtasa;
  int xcFunctional[numFunctionalIndices]; // specifies the energy functional to use in the calculations
                            // the first entry specifies the set of DFTs to use and the following numbers specify
                            // the actual functional:
                            // densityFunctional[0]=0: build in functionals
                            //     densityFunctional[1]=
                            //               1: von barth-hedin, j. phys. c5,1629(1972)
                            //               2: vosko--wilk-nusair, from g.s. painter, phys. rev. b24 4264(1981)
                            // densityFunctional[0]=1: functionals from libxc 
  //char *xcName;
  LibxcInterface libxcFunctional;
  int vSpinShiftFlag;      // if !=0 : shift the spin up and down potentials according to atom.vSpinShift
                           // this is used in WL-LSMS with moment magnitude fluctuations
  //double vSpinShift_min;   // vSpinShift_min, vSpinShift_max define the range for atom.vSpinShift
  //double vSpinShift_max;
  int fixRMT; // n_fix_mt from LSMS_1:
              //   0 -> set rmt to calculated inscribed sphere in atomic volume
              //   1 -> set it to rmt read in from potential file
  Real temperature;
  Real clight;
  int maxlmax;
  LSMSGlobals global;
  AngularMomentumIndices angularMomentumIndices;
  EnergyContourParameters energyContour;

// no. of Gaussian points for volume integration
  int ngaussr,ngaussq;
// prefered block size for zblock_lu: 0 use the default
  int zblockLUSize;

// Properties of the whole system:
  Real chempot;                // Chemical potential
  Real zvaltss;                // Total valence charge
  Real volumeTotal;            // Total cell volume
  Real volumeNorm;             // Volume renormalization factor
  Real volumeInterstitial;     // Total interstitial volume
  Real u0;                     // Contribution of the Muffin-tin zero potential to the Coulomb energy
  Real totalEnergy;            // Total energy
  //Real pressure;               // Pressure

// repeat the MPI rank from comm for reporting purposes
  int commRank;

  Real adjustContourBottom;    // if >0.0. set ebot to largestCorestate + adjustContourBottom
  Real largestCorestate;       // maximum of the core levels
};

extern const char *potentialTypeName[];

class AtomType {
public:
  AtomType() : pot_in_idx(-1), store_id(-1), forceZeroMoment(0) {}
  char name[4];
  int lmax,Z,Zc,Zs,Zv;
  int forceZeroMoment;
  int first_instance, number_of_instances;
  Real rsteps[4];
  Real rLIZ, rad;
  int node,local_id;
  int store_id;   // position in tmatStore
  int pot_in_idx;
  Real conc;
  int alloy_class;
};

class CrystalParameters {
public:
  int maxlmax;
  CrystalParameters() : bravais(3,3) {}
  void resize(size_t n) {type.resize(n); position.resize(3,n); evecs.resize(3,n);}
  void resizeTypes(size_t n) {types.resize(n);}
  Matrix<Real> bravais;
  Real omega; // bravais lattice volume
  int num_atoms,num_types;
  Matrix<Real> position,evecs;
  std::vector<int> type;
  std::vector<AtomType> types;
};

class LocalTypeInfo {
public:
  //LocalTypeInfo() : num_local(0), atom(0) {}
  //~LocalTypeInfo() { //if(num_local) delete[] atom;}
  // void setNumLocal(int n) {num_local=n; atom = new AtomData[n]; global_id.resize(n);}
  void setNumLocal(int n)
  {
    num_local=n; atom.resize(n); global_id.resize(n); n_per_type.resize(n);
    for(int i=0; i<num_local; i++) atom[i].reset();
  }
  void setGlobalId(int rank,CrystalParameters &crystal)
  {
    for(int i=0; i<crystal.num_types; i++)
    {
      if(rank==crystal.types[i].node)
      {
        global_id[crystal.types[i].local_id]=i;
        n_per_type[crystal.types[i].local_id]=crystal.types[i].number_of_instances;
      }
    }
  }
  void setMaxPts(int n) {for(int i=0; i<num_local; i++) atom[i].resizePotential(n);}
  void setMaxCore(int n) {for(int i=0; i<num_local; i++) atom[i].resizeCore(n);}
  int maxNrmat(void) {int v=0; for(int i=0; i<num_local; i++) if(atom[i].nrmat>v) v=atom[i].nrmat; return v;}
  int maxjws(void) {int m=0; for(int i=0; i<num_local; i++) if(atom[i].jws>m) m=atom[i].jws; return m;}
  int num_local;
  std::vector<int> global_id;
  std::vector<AtomData> atom;
  std::vector<int> n_per_type;

  int lDimTmatStore,blkSizeTmatStore;
  Matrix<Complex> tmatStore;
  std::vector<int> tmatStoreGlobalIdx;

  Real qrms[2];
  Real vrms[2];
};

// for Wang-Landau for metallic alloys
typedef std::vector< std::vector<AtomType> > AlloyMixingDesc;
typedef std::vector< std::vector<AtomData> > AlloyAtomBank;
// first index is alloy class (i.e. atomic types that can mix)
// second index is atomic component within a class


void printLSMSGlobals(FILE *f,LSMSSystemParameters &lsms);
void printLSMSSystemParameters(FILE *f,LSMSSystemParameters &lsms);
void printCrystalParameters(FILE *f, CrystalParameters &crystal);
void printAlloyParameters(FILE *f, AlloyMixingDesc &alloyDesc);
void printLocalTypeInfo(FILE *f, LocalTypeInfo &local);
void printLIZInfo(FILE * f, AtomData &atom);

#endif
