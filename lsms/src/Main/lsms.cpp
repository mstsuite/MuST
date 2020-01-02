/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

// #include <fenv.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define USE_PAPI 1
#ifdef USE_PAPI
#include <papi.h>
#endif

#ifdef USE_GPTL
#include "gptl.h"
#endif

#include <hdf5.h>

#include "lua.hpp"
//#include "lua.h"
//#include "lauxlib.h"
//#include "lualib.h"

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
#include "calculateDensities.hpp"
#include "mixing.hpp"
#include "calculateEvec.hpp"
#include "Potential/calculateChargesPotential.hpp"
#include "Potential/interpolatePotential.hpp"
#include "Potential/PotentialShifter.hpp"
#include "TotalEnergy/calculateTotalEnergy.hpp"
#include "SingleSite/checkAntiFerromagneticStatus.hpp"

#include "Misc/readLastLine.hpp"

#include "writeInfoEvec.cpp"
#include "write_restart.hpp"

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
#include "Accelerator/DeviceStorage.hpp"
// void * deviceStorage;
DeviceStorage *deviceStorage;
#endif
#ifdef BUILDKKRMATRIX_GPU
#include "Accelerator/buildKKRMatrix_gpu.hpp"

std::vector<DeviceConstants> deviceConstants;
// void *allocateDConst(void);
// void freeDConst(void *);
#endif
// std::vector<void *> deviceConstants;
// std::vector<void *> deviceStorage;


void initLSMSLuaInterface(lua_State *L);
int readInput(lua_State *L, LSMSSystemParameters &lsms, CrystalParameters &crystal, MixingParameters &mix, PotentialShifter &potentialShifter,
     AlloyMixingDesc &alloyDesc);
void buildLIZandCommLists(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                          CrystalParameters &crystal, LocalTypeInfo &local);
void setupVorpol(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local,
                 SphericalHarmonicsCoeficients &shc);

void calculateVolumes(LSMSCommunication &comm, LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);

/*
static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}
*/


int main(int argc, char *argv[])
{
  LSMSSystemParameters lsms;
  LSMSCommunication comm;
  CrystalParameters crystal;
  LocalTypeInfo local;
  MixingParameters mix;
  PotentialShifter potentialShifter;
  AlloyMixingDesc alloyDesc;
  AlloyAtomBank alloyBank;

  char inputFileName[128];

  Real eband;

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  initLSMSLuaInterface(L);

  // feenableexcept(FE_INVALID);

#ifdef USE_GPTL
  GPTLinitialize();
#endif
  initializeCommunication(comm);
  H5open();

  // set input file name (default 'i_lsms')
  strncpy(inputFileName, "i_lsms", 10);
  if (argc > 1)
    strncpy(inputFileName, argv[1], 120);

  lsms.global.iprpts = 1051;
  lsms.global.ipcore = 30;
  lsms.global.setIstop("main");
  lsms.global.iprint = 0;
  lsms.global.default_iprint = -1;
  lsms.global.print_node = 0;
  lsms.ngaussr = 10;
  lsms.ngaussq = 40;
  lsms.vSpinShiftFlag = 0;
#ifdef _OPENMP
  lsms.global.GPUThreads = std::min(12, omp_get_max_threads());
#else
  lsms.global.GPUThreads = 1;
#endif
  if (comm.rank == 0)
    lsms.global.iprint = 0;

  if (comm.rank == 0)
  {
    printf("LSMS_3: Program started\n");
    printf("Using %d MPI processes\n", comm.size);
#ifdef _OPENMP
    printf("Using %d OpenMP threads\n", omp_get_max_threads());
#endif
    acceleratorPrint();
#ifdef BUILDKKRMATRIX_GPU
    printf("Using GPU to build KKR matrix.\n");
#endif
    printf("Reading input file '%s'\n", inputFileName);
    fflush(stdout);

    if (luaL_loadfile(L, inputFileName) || lua_pcall(L,0,0,0))
    {
      fprintf(stderr, "!! Cannot run input file!!\n");
      exit(1);
    }

    printf("Loaded input file!\n");
    fflush(stdout);

    if (readInput(L, lsms, crystal, mix, potentialShifter, alloyDesc))
    {
      fprintf(stderr, "!! Something wrong in input file!!\n");
      exit(1);
    }

    printf("System information:\n");
    printf("===================\n");
    printf("Number of atoms        : %10d\n", crystal.num_atoms);
    printf("Number of atomic types : %10d\n", crystal.num_types);
    switch (lsms.mtasa)
    {
      case 1:
        printf("Performing Atomic Sphere Approximation (ASA) calculation\n");
        break;
      case 2:
        printf("Performing Atomic Sphere Approximation + Muffin-Tin (ASA-MT) calculation\n");
        break;
      default:
        printf("Performing Muffin-Tin (MT) calculation\n");
    }
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif
 
  communicateParameters(comm, lsms, crystal, mix, alloyDesc);
  if (comm.rank != lsms.global.print_node)
    lsms.global.iprint = lsms.global.default_iprint;
  // printf("maxlmax=%d\n",lsms.maxlmax);
  if(comm.rank == 0)
  {
    printf("communicated Parameters.\n");
    fflush(stdout);
  }

  local.setNumLocal(distributeTypes(crystal, comm));
  local.setGlobalId(comm.rank, crystal);

  if(comm.rank == 0)
  {
    printf("set global ids.\n");
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

  // set up exchange correlation functionals
  if (lsms.xcFunctional[0] == 1)         // use libxc functional
    lsms.libxcFunctional.init(lsms.n_spin_pola, lsms.xcFunctional); 

  lsms.angularMomentumIndices.init(2*crystal.maxlmax);
  sphericalHarmonicsCoeficients.init(2*crystal.maxlmax);

  gauntCoeficients.init(lsms, lsms.angularMomentumIndices, sphericalHarmonicsCoeficients);
  iFactors.init(lsms, crystal.maxlmax);

  double timeBuildLIZandCommList = MPI_Wtime();
  if (lsms.global.iprint >= 0)
  {
    printf("building the LIZ and Communication lists [buildLIZandCommLists]\n");
    fflush(stdout);
  }
  buildLIZandCommLists(comm, lsms, crystal, local);
  timeBuildLIZandCommList = MPI_Wtime() - timeBuildLIZandCommList;
  if (lsms.global.iprint >= 0)
  {
    printf("time for buildLIZandCommLists [num_local=%d]: %lf sec\n",
           local.num_local, timeBuildLIZandCommList);
    fflush(stdout);
  }

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

// initialize the potential accelerators (GPU)
// we need to know the max. size of the kkr matrix to invert: lsms.n_spin_cant*local.maxNrmat()
// which is only available after building the LIZ

  acceleratorInitialize(lsms.n_spin_cant*local.maxNrmat(), lsms.global.GPUThreads);
  local.tmatStore.pinMemory();
#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
  // deviceStorage = allocateDStore();
  deviceStorage = new DeviceStorage;
#endif
#ifdef BUILDKKRMATRIX_GPU
  deviceConstants.resize(local.num_local);
  // for(int i=0; i<local.num_local; i++) deviceConstants[i]=allocateDConst();
#endif

  for (int i=0; i<local.num_local; i++)
    local.atom[i].pmat_m.resize(lsms.energyContour.groupSize());  

// set maximal number of radial grid points and core states if reading from bigcell file
  local.setMaxPts(lsms.global.iprpts);
  local.setMaxCore(lsms.global.ipcore);

  if (lsms.global.iprint >= 0) printLSMSGlobals(stdout, lsms);
  if (lsms.global.iprint >= 0) printLSMSSystemParameters(stdout, lsms);
  if (lsms.global.iprint >= 1) printCrystalParameters(stdout, crystal);
  if (lsms.global.iprint >= 0) printAlloyParameters(stdout, alloyDesc);
  if (lsms.global.iprint >= 0)
  {
    fprintf(stdout,"LIZ for atom 0 on this node\n");
    printLIZInfo(stdout, local.atom[0]);
    if(local.atom[0].forceZeroMoment)
      fprintf(stdout,"\nMagnetic moment of atom 0 forced to be zero!\n\n");
  }
  if (lsms.global.iprint >= 1)
  {
    printCommunicationInfo(stdout, comm);
  }
  fflush(stdout);

//  initialAtomSetup(comm,lsms,crystal,local);

// the next line is a hack for initialization of potentials from scratch to work.


#ifdef LSMS_DEBUG
  if(lsms.global.iprint >= 0)
  {
    printf("Entering the Voronoi construction BEFORE loading the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  /* if(lsms.pot_in_type < 0) */ setupVorpol(lsms, crystal, local, sphericalHarmonicsCoeficients);

#ifdef LSMS_DEBUG
  if(lsms.global.iprint >= 0)
  {
    printf("Entering the LOADING of the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  loadPotentials(comm, lsms, crystal, local);

  if ( alloyDesc.size() > 0 )
  {
    if(lsms.global.iprint >= 0)
    {
      printf("Entering the LOADING of the alloy banks.\n");
      fflush(stdout);
    }
    loadAlloyBank(comm,lsms,alloyDesc,alloyBank); 
  }

// for testing purposes:
//  std::vector<Matrix<Real> > vrs;
//  vrs.resize(local.num_local);
//  for(int i=0; i<local.num_local; i++) vrs[i]=local.atom[i].vr;
// -------------------------------------

#ifdef LSMS_DEBUG
  if(lsms.global.iprint >= 0)
  {
    printf("Entering the Voronoi construction AFTER loading the potentials.\n");
    fflush(stdout);
  }
  MPI_Barrier(comm.comm);
#endif

  setupVorpol(lsms, crystal, local, sphericalHarmonicsCoeficients);

#ifdef LSMS_DEBUG
  MPI_Barrier(comm.comm);
#endif

// Generate new grids after new rmt is defined
  for (int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].generateNewMesh)
      interpolatePotential(lsms, local.atom[i]);
  }

  calculateVolumes(comm, lsms, crystal, local);

//  loadPotentials(comm,lsms,crystal,local);

// initialize Mixing
  Mixing *mixing;
  setupMixing(mix, mixing, lsms.global.iprint);

// need to calculate madelung matrices
  calculateMadelungMatrices(lsms, crystal, local);

  if (lsms.global.iprint >= 1)
  {
    printLocalTypeInfo(stdout, local);
  }

  calculateCoreStates(comm, lsms, local);
  if (lsms.global.iprint >= 0)
    printf("Finished calculateCoreStates(...)\n");

// check that vrs have not changed ...
//  bool vr_check=false;
//  for(int i=0; i<local.num_local; i++)
//  {
//    vr_check=true;
//    for(int j=0; j<vrs[i].n_row();j++)
//      for(int k=0; k<vrs[i].n_col(); k++)
//        vr_check=vr_check && (vrs[i](j,k)==local.atom[i].vr(j,k));
//    if(!vr_check)
//      printf("Potential %d has changed!!\n",i);
//  }
//  printf("Potentials have been checked\n");
// --------------------------------------------

// meis: Test if Potential for atoms 0 and 1 are the same
/*
  if(local.num_local>1)
  {
    bool vr_check=true;
    for(int j=0; j<local.atom[0].vr.n_row();j++)
      for(int k=0; k<local.atom[0].vr.n_col(); k++)
        vr_check=vr_check && (local.atom[0].vr(j,k)==local.atom[1].vr(j,k));
    if(!vr_check)
      printf("Potentials 0 and 1 are different!!\n");
    printf("Potentials have been checked\n");
  }
*/

  if (lsms.n_spin_cant > 1)
  {
    for (int i=0; i<local.num_local; i++)
      local.atom[i].get_b_basis();
  }
  else
  {
    for (int i=0; i<local.num_local; i++)
      local.atom[i].reset_b_basis();
  }

  mixing -> prepare(comm, lsms, local.atom);

#ifdef USE_PAPI
  #define NUM_PAPI_EVENTS 2
  int hw_counters = PAPI_num_counters();
  if (hw_counters > NUM_PAPI_EVENTS)
    hw_counters = NUM_PAPI_EVENTS;
  int papi_events[NUM_PAPI_EVENTS]; 
  char *papi_event_name[] = {"PAPI_TOT_INS", "PAPI_FP_OPS"};
  // get events from names:
  for (int i=0; i<NUM_PAPI_EVENTS; i++)
  {
    if (PAPI_event_name_to_code(papi_event_name[i], &papi_events[i]) != PAPI_OK)
      if (hw_counters > i)
        hw_counters = i;
  }
  long long papi_values[NUM_PAPI_EVENTS+4];
  if (hw_counters > NUM_PAPI_EVENTS)
    hw_counters = NUM_PAPI_EVENTS;
  long long papi_real_cyc_0 = PAPI_get_real_cyc();
  long long papi_real_usec_0 = PAPI_get_real_usec();
  long long papi_virt_cyc_0 = PAPI_get_virt_cyc();
  long long papi_virt_usec_0 = PAPI_get_virt_usec();
  PAPI_start_counters(papi_events, hw_counters);
#endif

// -----------------------------------------------------------------------------
//                                 MAIN SCF LOOP
// -----------------------------------------------------------------------------

  bool converged = false;

  if (lsms.global.iprint >= 0)
    printf("Total number of iterations:%d\n", lsms.nscf);

  double timeScfLoop = MPI_Wtime();
  double timeCalcChemPot = 0.0;

  int iterationStart = 0;
  int potentialWriteCounter = 0;

  FILE *kFile = NULL;
  if (comm.rank == 0)
  {
    iterationStart = readNextIterationNumber("k.out");
    kFile = fopen("k.out","a");
  }

  int iteration;
  for (iteration=0; iteration<lsms.nscf && !converged; iteration++)
  {
    if (lsms.global.iprint >= 0)
      printf("SCF iteration %d:\n", iteration);

    // Calculate band energy
    energyContourIntegration(comm, lsms, local);
    double dTimeCCP = MPI_Wtime();
    // if(!lsms.global.checkIstop("buildKKRMatrix"))

    // Calculate chemical potential 
    calculateChemPot(comm, lsms, local, eband);
    dTimeCCP = MPI_Wtime() - dTimeCCP;
    timeCalcChemPot += dTimeCCP;

    // Calculate magnetic moments for each site and check if spin has flipped
    calculateEvec(lsms, local);
    mixEvec(lsms, local, 0.0);
    for (int i=0; i<local.num_local; i++) {
      local.atom[i].newConstraint();
      checkIfSpinHasFlipped(lsms, local.atom[i]);
    }

    // Calculate charge densities, potentials, and total energy
    calculateAllLocalChargeDensities(lsms, local);
    calculateChargesPotential(comm, lsms, local, crystal, 0);
    checkAllLocalCharges(lsms, local);
    calculateTotalEnergy(comm, lsms, local, crystal);

    // Mix charge density
    mixing -> updateChargeDensity(comm, lsms, local.atom);
 
    // Recalculate core states
    // - swap core state energies for different spin channels first if spin has flipped
    //   (from LSMS 1: lsms_main.f:2101-2116)
    for (int i=0; i<local.num_local; i++) {
      if (local.atom[i].spinFlipped)
      {
        checkIfSpinHasFlipped(lsms, local.atom[i]);
        if (!local.atom[i].spinFlipped)
          swapCoreStateEnergies(local.atom[i]);
      }
    }
    calculateCoreStates(comm, lsms, local);

    // If charge is mixed, recalculate potential and mix (need a flag for this from input)
    calculateChargesPotential(comm, lsms, local, crystal, 1);
    mixing -> updatePotential(comm, lsms, local.atom);

    // Real rms = 0.5 * (local.qrms[0] + local.qrms[1]);
    Real rms = 0.0;
    for(int i=0; i<local.num_local; i++)
      rms = std::max(rms, 0.5*(local.atom[i].qrms[0]+local.atom[i].qrms[1]));
    globalMax(comm, rms);
    
// check for convergence
    converged = rms < lsms.rmsTolerance;
    /*
    converged = true;
    for (int i=0; i<local.num_local; i++)
    {
      converged = converged
                && (0.5*(local.atom[i].qrms[0]+local.atom[i].qrms[1])<lsms.rmsTolerance);
    }
    globalAnd(comm, converged);
    */

    if (comm.rank == 0)
    {
      printf("Band Energy = %lf Ry %10s", eband, "");
      printf("Fermi Energy = %lf Ry\n", lsms.chempot);
      printf("Total Energy = %lf Ry\n", lsms.totalEnergy);
      printf("RMS = %lg\n",rms);
      if(lsms.global.iprint > 0)
      {
        printf("  qrms[0] = %lg   qrms[1] = %lg\n",local.qrms[0], local.qrms[1]);
        printf("  local.atom[i]:\n");
        for (int i=0; i<local.num_local; i++)
        {
          printf("  %d : qrms[0] = %lg   qrms[1] = %lg\n",i,local.atom[i].qrms[0], local.atom[i].qrms[1]);
          printf("  %d : vrms[0] = %lg   vrms[1] = %lg\n",i,local.atom[i].vrms[0], local.atom[i].vrms[1]);
        }
      }
    }

    if (kFile != NULL)
    {
      fprintf(kFile,"%4d %20.12lf %12.6lf %12.6lf  %14.10lf\n",
              iterationStart+iteration, lsms.totalEnergy, lsms.chempot, local.atom[0].mtotws, rms);
      fflush(kFile);
    }

    // Recalculate core states for new potential if we are performing scf calculations
    calculateCoreStates(comm, lsms, local);

    // Periodically write the new potential for scf calculations 
    potentialWriteCounter++;
    if ((lsms.pot_out_type >= 0 && potentialWriteCounter >= lsms.writeSteps)
        || converged)
    {
      if (comm.rank == 0) std::cout << "Writing new potentials and restart file.\n";
      writePotentials(comm, lsms, crystal, local);
      potentialWriteCounter = 0;
      if (comm.rank == 0)
      { 
        writeRestart("i_lsms.restart", lsms, crystal, mix, potentialShifter, alloyDesc);
      }
    }

  }

 timeScfLoop = MPI_Wtime() - timeScfLoop;
  
  writeInfoEvec(comm, lsms, crystal, local, eband, lsms.infoEvecFileOut);
  if(lsms.localAtomDataFile[0]!=0)
    writeLocalAtomData(comm, lsms, crystal, local, eband, lsms.localAtomDataFile);

  if (kFile != NULL)
    fclose(kFile);

 

// -----------------------------------------------------------------------------

#ifdef USE_PAPI
  PAPI_stop_counters(papi_values,hw_counters);
  papi_values[hw_counters  ] = PAPI_get_real_cyc()-papi_real_cyc_0;
  papi_values[hw_counters+1] = PAPI_get_real_usec()-papi_real_usec_0;
  papi_values[hw_counters+2] = PAPI_get_virt_cyc()-papi_virt_cyc_0;
  papi_values[hw_counters+3] = PAPI_get_virt_usec()-papi_virt_usec_0;
  long long accumulated_counters[NUM_PAPI_EVENTS+4];
  MPI_Reduce(papi_values,accumulated_counters,hw_counters+4,
             MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  if (comm.rank == 0)
  {
    for (int i=0; i<hw_counters; i++)
      std::cout<<"Accumulated: "<<(papi_event_name[i])<<" = "<<(accumulated_counters[i])<<"\n";
    std::cout<<"PAPI accumulated real cycles : "<<(accumulated_counters[hw_counters])<<"\n";
    std::cout<<"PAPI accumulated user cycles : "<<(accumulated_counters[hw_counters+2])<<"\n";
    double gflops_papi = ((double)accumulated_counters[1])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_double = ((double)accumulated_counters[2])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gflops_hw_single = ((double)accumulated_counters[3])/
      (1000.0*(double)papi_values[hw_counters+1]);
    double gips = ((double)accumulated_counters[0])/(1000.0*(double)papi_values[hw_counters+1]);
    std::cout<<"PAPI_FP_OPS real GFLOP/s : "<<(gflops_papi)<<"\n";
    std::cout<<"PAPI hw double real GFLOP/s : "<<(gflops_hw_double)<<"\n";
    std::cout<<"PAPI hw single real GFLOP/s : "<<(gflops_hw_single)<<"\n";
    std::cout<<"PAPI real GINST/s : "<<(gips)<<"\n";
    std::cout<<"Time (s) : " << (double)papi_values[hw_counters+1] << "\n";
  }
#endif

  if (lsms.pot_out_type >= 0)
  {
    if (comm.rank == 0) std::cout << "Writing new potentials.\n";
    writePotentials(comm, lsms, crystal, local);
    if (comm.rank == 0)
    {
      std::cout << "Writing restart file.\n";
      writeRestart("i_lsms.restart", lsms, crystal, mix, potentialShifter, alloyDesc);
    }
  }

  long long fomScale = calculateFomScale(comm, local);

  if (comm.rank == 0)
  {
    printf("Band Energy = %.15lf Ry\n", eband);
    printf("Fermi Energy = %.15lf Ry\n", lsms.chempot);
    printf("Total Energy = %.15lf Ry\n", lsms.totalEnergy);
    printf("timeScfLoop[rank==0] = %lf sec\n", timeScfLoop);
    printf("     number of iteration:%d\n",iteration);
    printf("timeScfLoop/iteration = %lf sec\n", timeScfLoop / (double)iteration);
    // printf(".../lsms.nscf = %lf sec\n", timeScfLoop / (double)lsms.nscf);
    printf("timeCalcChemPot[rank==0]/iteration = %lf sec\n", timeCalcChemPot / (double)iteration);
    // printf("timeCalcChemPot[rank==0]/lsms.nscf = %lf sec\n", timeCalcChemPot / (double)lsms.nscf);
    printf("timeBuildLIZandCommList[rank==0]: %lf sec\n",
           timeBuildLIZandCommList);
    // fom = [ \sum_#atoms (LIZ * (lmax+1)^2)^3 ] / time per iteration
    //     = [ \sum_#atoms (LIZ * (lmax+1)^2)^3 ] * lsms.nscf / timeScfLoop
    // fom_e = fom * energy contour points
    // fomScale = \sum_#atoms (LIZ * (lmax+1)^2)^3
    // energyContourPoints
    long energyContourPoints = 1;
    if(lsms.energyContour.grid==2)
    {
      energyContourPoints = lsms.energyContour.npts+1;
    }
    printf("FOM Scale = %lf\n",(double)fomScale);
    printf("Energy Contour Points = %ld\n",energyContourPoints);
    printf("FOM = %lg/sec\n",fomScale * (double)iteration / timeScfLoop);
    // printf("FOM = %lg/sec\n",fomScale * (double)lsms.nscf / timeScfLoop);
    printf("FOM * energyContourPoints = = %lg/sec\n",
            (double)energyContourPoints * (double)fomScale * (double)iteration / timeScfLoop);
    //         (double)energyContourPoints * (double)fomScale * (double)lsms.nscf / timeScfLoop);
  }

  local.tmatStore.unpinMemory();

#ifdef BUILDKKRMATRIX_GPU
  // for(int i=0; i<local.num_local; i++) freeDConst(deviceConstants[i]);
#endif

#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
  // freeDStore(deviceStorage);
  delete deviceStorage;
#endif
#ifdef BUILDKKRMATRIX_GPU
  deviceConstants.clear();
#endif

  acceleratorFinalize();
  delete mixing;

#ifdef USE_GPTL
  GPTLpr(comm.rank);
#endif

  H5close();
  finalizeCommunication();
  lua_close(L);
  return 0;
}
