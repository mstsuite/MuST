#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "Madelung/Madelung.hpp"
#include "VORPOL/VORPOL.hpp"
#include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "calculateChemPot.hpp"

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

void initLSMSLuaInterface(lua_State *L);
int readInput(lua_State *L, LSMSSystemParameters &lsms, CrystalParameters &crystal);
void buildLIZandCommLists(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                          CrystalParameters &crystal, LocalTypeInfo &local);
void setupVorpol(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local,
                 SphericalHarmonicsCoeficients &shc);

int main(int argc, char *argv[])
{
  LSMSSystemParameters lsms;
  LSMSCommunication comm;
  CrystalParameters crystal;
  LocalTypeInfo local;

  Real eband;


  lua_State *L=lua_open();
  luaL_openlibs(L);
  initLSMSLuaInterface(L);

  initializeCommunication(comm);

  lsms.global.iprpts=1051;
  lsms.global.ipcore=15;
  lsms.global.setIstop("main");
  lsms.global.iprint=-1;
  lsms.ngaussr=10;
  lsms.ngaussq=40;
  if(comm.rank==0) lsms.global.iprint=0;

  if(comm.rank==0)
  {
    printf("LSMS_3: Program started\n");
    printf("Using %d MPI processes\n",comm.size);
#ifdef _OPENMP
    printf("Using %d OpenMP threads\n",omp_get_max_threads());
#endif
    printf("Reading input file '%s'\n",argv[1]);

    if(luaL_loadfile(L, argv[1]) || lua_pcall(L,0,0,0))
    {
      fprintf(stderr,"!! Cannot run input file!!\n");
      exit(1);
    }

    if(readInput(L,lsms,crystal))
    {
      fprintf(stderr,"!! Something wrong in input file!!\n");
      exit(1);
    }
  }

  communicateParameters(comm,lsms,crystal);
  // printf("maxlmax=%d\n",lsms.maxlmax);

  local.setNumLocal(distributeTypes(crystal, comm));
  local.setGlobalId(comm.rank,crystal);

  lsms.angularMomentumIndices.init(2*crystal.maxlmax);
  sphericalHarmonicsCoeficients.init(2*crystal.maxlmax);

  gauntCoeficients.init(lsms,lsms.angularMomentumIndices,sphericalHarmonicsCoeficients);
  iFactors.init(lsms,crystal.maxlmax);

  buildLIZandCommLists(comm, lsms, crystal, local);

// initialize the potential accelerators (GPU)
// we need to know the max. size of the kkr matrix to invert: lsms.n_spin_cant*local.maxNrmat()
// which is only available after building the LIZ

  acceleratorInitialize(lsms.n_spin_cant*local.maxNrmat());
  

// set maximal number of radial grid points and core states if reading from bigcell file
  local.setMaxPts(lsms.global.iprpts);
  local.setMaxCore(lsms.global.ipcore);

  if(lsms.global.iprint>=2)
  {
    printLSMSSystemParameters(stdout,lsms);
    printCrystalParameters(stdout,crystal);
  }
  if(lsms.global.iprint>=1)
  {
    fprintf(stdout,"LIZ for atom 0 on this node\n");
    printLIZInfo(stdout,local.atom[0]);
    printCommunicationInfo(stdout, comm);
  }
 
 
  loadPotentials(comm,lsms,crystal,local);

  setupVorpol(lsms,crystal,local,sphericalHarmonicsCoeficients);

// need to calculate madelung matrices
  calculateMadelungMatrices(lsms,crystal,local);

  if(lsms.global.iprint>=1)
  {
    printLocalTypeInfo(stdout,local);
  }

  calculateCoreStates(comm,lsms,local);
  if(lsms.global.iprint>=0)
    printf("Finished calculateCoreStates(...)\n");

// -----------------------------------------------------------------------------
//                                 MAIN SCF LOOP
// -----------------------------------------------------------------------------

    energyContourIntegration(comm,lsms,local);
    calculateChemPot(comm,lsms,local,eband);

// -----------------------------------------------------------------------------

  acceleratorFinalize();
  finalizeCommunication();
  lua_close(L);
  return 0;
}
