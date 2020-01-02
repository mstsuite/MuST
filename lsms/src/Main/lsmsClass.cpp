// lsms class to encapsulate a version of LSMS_1.9 for use in gWL etc.

#include <mpi.h>
#include <iostream>
#include <vector>

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
#include "Potential/calculateChargesPotential.hpp"
#include "Potential/interpolatePotential.hpp"
#include "Potential/PotentialShifter.hpp"
#include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "calculateChemPot.hpp"
#include "calculateDensities.hpp"
#include "calculateEvec.hpp"
#include "TotalEnergy/calculateTotalEnergy.hpp"
#include "SingleSite/checkAntiFerromagneticStatus.hpp"

#include "lsmsClass.hpp"

#ifdef _OPENMP
#include <omp.h>
#else
#ifndef LSMS_DUMMY_OPENMP
#define LSMS_DUMMY_OPENMP
inline int omp_get_max_threads() {return 1;}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
#endif
#endif

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
#include "Accelerator/DeviceStorage.hpp"
// void *deviceStorage;
DeviceStorage *deviceStorage;
#endif

#ifdef BUILDKKRMATRIX_GPU
// void *allocateDStore(void);
// void freeDStore(void *);
// void *allocateDConst(void);
// void freeDConst(void *);
#include "Accelerator/buildKKRMatrix_gpu.hpp"

std::vector<DeviceConstants> deviceConstants;
//std::vector<void *> deviceConstants;
#endif

void initLSMSLuaInterface(lua_State *L);
int readInput(lua_State *L, LSMSSystemParameters &lsms, CrystalParameters &crystal, MixingParameters &mix, PotentialShifter &potentialshifter,
       AlloyMixingDesc &alloyDesc);
void buildLIZandCommLists(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                          CrystalParameters &crystal, LocalTypeInfo &local);
void setupVorpol(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local,
                 SphericalHarmonicsCoeficients &shc);
void calculateVolumes(LSMSCommunication &comm, LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);

// Constructor
//LSMS::LSMS(MPI_Comm _comm, const char* i_lsms, const char* out_prefix, PotentialShifter &ps)
//LSMS::LSMS(MPI_Comm _comm, const char* i_lsms, const char* out_prefix)
LSMS::LSMS(MPI_Comm _comm, const char* i_lsms, const char* out_prefix, int my_group)
{

  myWalkerID = my_group;

  char h[] = "NO OUTPUT";

  if (out_prefix == NULL)
    strncpy(prefix, h, 255);
  else
    strncpy(prefix, out_prefix, 255);

  energyLoopCount = 0;

/*
  if (rank == 0)
    std::cout<<"initializing lsms::lsms with i_lsms="<<i_lsms<<" and out_prefix="
	       <<prefix<<std::endl;
*/

  efTol = 1.0e-6;
  energyTolerance = 1.0e-8;
  // rmsTolerance = 1.0e-6;

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  initLSMSLuaInterface(L);

  initializeCommunication(comm, _comm);

  lsms.global.setIstop("main");
  lsms.global.iprpts         = 1051;
  lsms.global.ipcore         = 30;
  lsms.global.iprint         = -1;
  lsms.global.default_iprint = -1;
  lsms.global.print_node     = 0;
  lsms.global.GPUThreads     = 1;
#ifdef _OPENMP
  lsms.global.GPUThreads     = std::min(16,omp_get_max_threads());
#else
  lsms.global.GPUThreads     = 1;
#endif
  lsms.ngaussr               = 10;
  lsms.ngaussq               = 40;
  // if (comm.rank == comm.size-1) lsms.global.iprint = 0;
  //lsms.vSpinShiftFlag = 0;

  if (comm.rank == 0)
  {
    if (strncmp(prefix,"0_",2) == 0)
    {
      printf("LSMS_3: Program started\n");
      printf("Using %d MPI processes\n", comm.size);
      printf("Using %d OpenMP threads\n", omp_get_max_threads());
      acceleratorPrint();
      printf("Reading input file '%s'\n", i_lsms);
    }

    if(luaL_loadfile(L, i_lsms) || lua_pcall(L,0,0,0))
    {
      fprintf(stderr,"!! Cannot run input file '%s'!!\n", i_lsms);
      exit(1);
    }

    if(readInput(L,lsms,crystal,mix,potentialShifter,alloyDesc))
    {
      fprintf(stderr, "!! Something wrong in input file!!\n");
      exit(1);
    }
  }

  communicateParameters(comm, lsms, crystal, mix, alloyDesc);
  communicatePotentialShiftParameters(comm, potentialShifter);
  if (comm.rank != lsms.global.print_node)
    lsms.global.iprint = lsms.global.default_iprint;

  local.setNumLocal(distributeTypes(crystal, comm));
  max_num_local = local.num_local;
  globalMax(comm, max_num_local);
  local.setGlobalId(comm.rank, crystal);

  // set up exchange correlation functionals
  if(lsms.xcFunctional[0] == 1)  // libxc functional
    lsms.libxcFunctional.init(lsms.n_spin_pola, lsms.xcFunctional);

  lsms.angularMomentumIndices.init(2 * crystal.maxlmax);
  sphericalHarmonicsCoeficients.init(2 * crystal.maxlmax);

  gauntCoeficients.init(lsms, lsms.angularMomentumIndices, sphericalHarmonicsCoeficients);
  iFactors.init(lsms, crystal.maxlmax);

  buildLIZandCommLists(comm, lsms, crystal, local);

  // initialize the potential accelerators (GPU)
  // we need to know the max. size of the kkr matrix to invert: lsms.n_spin_cant*local.maxNrmat()
  // which is only available after building the LIZ

  acceleratorInitialize(lsms.n_spin_cant * local.maxNrmat(), lsms.global.GPUThreads);
  local.tmatStore.pinMemory();
#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
  // deviceStorage = allocateDStore();
  deviceStorage = new DeviceStorage;
#endif
#ifdef BUILDKKRMATRIX_GPU
  deviceConstants.resize(local.num_local);
  // for(int i=0; i<local.num_local; i++) deviceConstants[i] = allocateDConst();
#endif

  for(int i=0; i<local.num_local; i++)
    local.atom[i].pmat_m.resize(lsms.energyContour.groupSize());

  // set maximal number of radial grid points and core states if reading from bigcell file
  local.setMaxPts(lsms.global.iprpts);
  local.setMaxCore(lsms.global.ipcore);

  if (lsms.global.iprint >= 0)
  {
    printLSMSSystemParameters(stdout, lsms);
    printCrystalParameters(stdout, crystal);
    printAlloyParameters(stdout,alloyDesc);
  }
  if (lsms.global.iprint >= 1)
  {
    fprintf(stdout,"LIZ for atom 0 on this node\n");
    printLIZInfo(stdout, local.atom[0]);
    printCommunicationInfo(stdout, comm);
  }
 
  loadPotentials(comm, lsms, crystal, local);
  setupVorpol(lsms, crystal, local, sphericalHarmonicsCoeficients);

  // for Wang-Landau for alloys
  if ( alloyDesc.size() > 0 ) 
    loadAlloyBank(comm,lsms,alloyDesc,alloyBank); 

  // Generate new grids after new rmt is defined
  for (int i=0; i<local.num_local; i++)
  {
    if (local.atom[i].generateNewMesh)
      interpolatePotential(lsms, local.atom[i]);
  }

  calculateVolumes(comm, lsms, crystal, local);

  // need to calculate madelung matrices
  calculateMadelungMatrices(lsms, crystal, local);

  if (lsms.global.iprint >= 0)
    printLocalTypeInfo(stdout, local);

  // initialize Mixing
  setupMixing(mix, mixing, lsms.global.iprint);

  // set and copy potentialShifter
  potentialShifter.resize(local.num_local);
  potentialShifter.resetPotentials(local);

  calculateCoreStates(comm, lsms, local);
  if (lsms.global.iprint >= 0)
    printf("Finished calculateCoreStates(...)\n");

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

  // do the same for alloy potentials
  if ( lsms.n_spin_cant > 1 ) {
    for (int i = 0; i < alloyBank.size(); i++)
      for (int j = 0; j < alloyBank[i].size(); j++)
        alloyBank[i][j].get_b_basis();
  } 
  else {
    for (int i = 0; i < alloyBank.size(); i++)
      for (int j = 0; j < alloyBank[i].size(); j++)
        alloyBank[i][j].reset_b_basis();
  }

  mixing -> prepare(comm, lsms, local.atom);


  LSMS_version = 3000;
  //  if (comm.rank == 0)
  //    std::cout<<prefix<<" lsms::lsms - LSMS version is "<<LSMS_version<<std::endl; 

  lua_close(L);
}


// Destructor
LSMS::~LSMS()
{
  local.tmatStore.unpinMemory();
#ifdef BUILDKKRMATRIX_GPU
  // for (int i=0; i<local.num_local; i++)
  //   freeDConst(deviceConstants[i]);
#endif
#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
  // freeDStore(deviceStorage);
  delete deviceStorage;
#endif
#ifdef BUILDKKRMATRIX_GPU
  deviceConstants.clear();
#endif
  acceleratorFinalize();
  // finalizeCommunication();
}


void LSMS::setEvec(std::vector<std::vector<Real> > &ev)
{
  MPI_Status status;
  //distribute the evecs to the respective processes
  if (comm.rank == 0)
  {
    std::vector<MPI_Request> request(crystal.num_types);
    int n_req = 0;
    for (int p=0; p<crystal.num_types; p++)
    {
      int l = crystal.types[p].local_id;
      int n = crystal.types[p].node;
      if (n == 0)
      {
        local.atom[l].evec[0] = (ev[p])[0];
        local.atom[l].evec[1] = (ev[p])[1];
        local.atom[l].evec[2] = (ev[p])[2];
      }
      else
        MPI_Isend(&(ev[p])[0], 3, MPI_DOUBLE, n, l, comm.comm, &request[n_req++]);
    }
    for (int i=0; i<n_req; i++)
      MPI_Wait(&request[i],&status);
  }
  else
  {
    std::vector<MPI_Request> request(local.num_local);
    for (int p=0; p<local.num_local; p++)
      MPI_Irecv(&local.atom[p].evec[0], 3, MPI_DOUBLE, 0, p, comm.comm, &request[p]);
    for (int i=0; i<local.num_local; i++)
      MPI_Wait(&request[i], &status);
  }

  // make sure evecs are normalized
  for (int p=0; p<local.num_local; p++)
  {
    local.atom[p].reset();
    Real m = std::sqrt(local.atom[p].evec[0] * local.atom[p].evec[0] +
                       local.atom[p].evec[1] * local.atom[p].evec[1] +
		       local.atom[p].evec[2] * local.atom[p].evec[2]);

    local.atom[p].evec[0] /= m;
    local.atom[p].evec[1] /= m;
    local.atom[p].evec[2] /= m;
    local.atom[p].evecNew[0] = local.atom[p].evec[0];
    local.atom[p].evecNew[1] = local.atom[p].evec[1];
    local.atom[p].evecNew[2] = local.atom[p].evec[2];
    local.atom[p].get_b_basis();
  }
}


void LSMS::setEvec(Real *ev)
{
  MPI_Status status;
  //distribute the evecs to the respective processes
  if (comm.rank == 0)
  {
    std::vector<MPI_Request> request(crystal.num_types);
    int n_req = 0;
    for (int p=0; p<crystal.num_types; p++)
    {
      int l = crystal.types[p].local_id;
      int n = crystal.types[p].node;
      if (n == 0)
      {
        local.atom[l].evec[0] = ev[3*p+0];
        local.atom[l].evec[1] = ev[3*p+1];
        local.atom[l].evec[2] = ev[3*p+2];
      } 
      else
        MPI_Isend(&ev[3*p], 3, MPI_DOUBLE, n, l, comm.comm, &request[n_req++]);
    }
    for (int i=0; i<n_req; i++)
      MPI_Wait(&request[i], &status);
  }
  else
  {
    std::vector<MPI_Request> request(local.num_local);
    for (int p=0; p<local.num_local; p++)
      MPI_Irecv(&local.atom[p].evec[0], 3, MPI_DOUBLE, 0, p, comm.comm, &request[p]);
    for (int i=0; i<local.num_local; i++)
      MPI_Wait(&request[i], &status);
  }

  // make sure evecs are normalized
  for(int p=0; p<local.num_local; p++)
  {
    local.atom[p].reset();
    Real m = std::sqrt(local.atom[p].evec[0] * local.atom[p].evec[0] +
                       local.atom[p].evec[1] * local.atom[p].evec[1] +
		       local.atom[p].evec[2] * local.atom[p].evec[2]);

    local.atom[p].evec[0] /= m;
    local.atom[p].evec[1] /= m;
    local.atom[p].evec[2] /= m;
    local.atom[p].evecNew[0] = local.atom[p].evec[0];
    local.atom[p].evecNew[1] = local.atom[p].evec[1];
    local.atom[p].evecNew[2] = local.atom[p].evec[2];
    local.atom[p].get_b_basis();
  }

}


void LSMS::setEvecAndSpinPotentialShift(Real *ev)
{
  MPI_Status status;
  potentialShifter.vSpinShiftFlag = true;

  //distribute the evecs to the respective processes
  if (comm.rank == 0)
  {
    std::vector<MPI_Request> request(crystal.num_types);
    int n_req = 0;
    for(int p=0; p<crystal.num_types; p++)
    {
      int l = crystal.types[p].local_id;
      int n = crystal.types[p].node;
      if (n == 0)
      {
        local.atom[l].evec[0] = ev[3*p+0];
        local.atom[l].evec[1] = ev[3*p+1];
        local.atom[l].evec[2] = ev[3*p+2];
      }
      else {
        MPI_Isend(&ev[3*p], 3, MPI_DOUBLE, n, l, comm.comm, &request[n_req++]);
      }
    }
    for(int i=0; i<n_req; i++)
      MPI_Wait(&request[i], &status);
  } 
  else
  {
    std::vector<MPI_Request> request(local.num_local);
    for(int p=0; p<local.num_local; p++)
      MPI_Irecv(&local.atom[p].evec[0], 3, MPI_DOUBLE, 0, p, comm.comm, &request[p]);
    for (int i=0; i<local.num_local; i++)
      MPI_Wait(&request[i], &status);
  }

  for (int p=0; p<local.num_local; p++)
  {
    local.atom[p].reset();
    Real m = std::sqrt(local.atom[p].evec[0] * local.atom[p].evec[0] +
                       local.atom[p].evec[1] * local.atom[p].evec[1] +
		       local.atom[p].evec[2] * local.atom[p].evec[2]);

    local.atom[p].vSpinShift = m - 1.0;

    local.atom[p].evec[0] /= m;
    local.atom[p].evec[1] /= m;
    local.atom[p].evec[2] /= m;
    local.atom[p].evecNew[0] = local.atom[p].evec[0];
    local.atom[p].evecNew[1] = local.atom[p].evec[1];
    local.atom[p].evecNew[2] = local.atom[p].evec[2];
    local.atom[p].evecOut[0] = local.atom[p].evec[0];
    local.atom[p].evecOut[1] = local.atom[p].evec[1];
    local.atom[p].evecOut[2] = local.atom[p].evec[2];
    local.atom[p].get_b_basis();
  }

}

void LSMS::setEvecAndSpinPotentialShift4(Real *evAndSpinShifts)
{
  MPI_Status status;
  potentialShifter.vSpinShiftFlag = true;

  //distribute the evecs to the respective processes
  if (comm.rank == 0)
  {
    std::vector<MPI_Request> request(crystal.num_types);
    int n_req = 0;
    for (int p=0; p<crystal.num_types; p++)
    {
      int l = crystal.types[p].local_id;
      int n = crystal.types[p].node;
      if (n == 0)
      {
        local.atom[l].evec[0] = evAndSpinShifts[4*p+0];
        local.atom[l].evec[1] = evAndSpinShifts[4*p+1];
        local.atom[l].evec[2] = evAndSpinShifts[4*p+2];
        local.atom[l].vSpinShift = evAndSpinShifts[4*p+3];
      }
      else
        MPI_Isend(&evAndSpinShifts[4*p], 4, MPI_DOUBLE, n, l, comm.comm, &request[n_req++]);
    }
    for (int i=0; i<n_req; i++)
      MPI_Wait(&request[i], &status);
  }
  else
  {
    std::vector<MPI_Request> request(local.num_local);
    Real recvBuffer[4*local.num_local];

    for (int p=0; p<local.num_local; p++)
      MPI_Irecv(&recvBuffer[4*p], 4, MPI_DOUBLE, 0, p, comm.comm, &request[p]);

    for(int p=0; p<local.num_local; p++)
    {
      MPI_Wait(&request[p], &status);
    
      local.atom[p].evec[0] = recvBuffer[4*p+0];
      local.atom[p].evec[1] = recvBuffer[4*p+1];
      local.atom[p].evec[2] = recvBuffer[4*p+2];
      local.atom[p].vSpinShift = recvBuffer[4*p+3];
    }
  }

  for (int p=0; p<local.num_local; p++)
  {
    local.atom[p].reset();
    Real m = std::sqrt(local.atom[p].evec[0] * local.atom[p].evec[0]+
                       local.atom[p].evec[1] * local.atom[p].evec[1]+
                       local.atom[p].evec[2] * local.atom[p].evec[2]);

    local.atom[p].evec[0] /= m;
    local.atom[p].evec[1] /= m;
    local.atom[p].evec[2] /= m;
    local.atom[p].evecNew[0] = local.atom[p].evec[0];
    local.atom[p].evecNew[1] = local.atom[p].evec[1];
    local.atom[p].evecNew[2] = local.atom[p].evec[2];
    local.atom[p].evecOut[0] = local.atom[p].evec[0];
    local.atom[p].evecOut[1] = local.atom[p].evec[1];
    local.atom[p].evecOut[2] = local.atom[p].evec[2];
    local.atom[p].get_b_basis();
  }

}


void LSMS::getEvec(std::vector<std::vector<Real> > &ev)
{
  Array3d<Real> r_buf(4,max_num_local,comm.size);
  Matrix<Real> s_buf(4,max_num_local);

  for (int i=0; i<max_num_local; i++)
    s_buf(0,i) = -1.0;

  for (int i=0; i<local.num_local; i++)
  {
    s_buf(0,i) = Real(local.global_id[i]);
    s_buf(1,i) = local.atom[i].evec[0];
    s_buf(2,i) = local.atom[i].evec[1];
    s_buf(3,i) = local.atom[i].evec[2];
  }

  MPI_Gather(&s_buf(0,0),   4*max_num_local, MPI_DOUBLE,
             &r_buf(0,0,0), 4*max_num_local, MPI_DOUBLE, 0, comm.comm);

  if (comm.rank == 0)
  {
    for (int p=0; p<comm.size; p++)
    {
      for (int i=0; i<max_num_local && r_buf(0,i,p)>=0; i++)
      {
        int j = r_buf(0,i,p);
        (ev[j])[0] = r_buf(1,i,p);
        (ev[j])[1] = r_buf(2,i,p);
        (ev[j])[2] = r_buf(3,i,p);
      }
    }
  }
}


void LSMS::getEvec(Real *ev)
{
  Array3d<Real> r_buf(4,max_num_local,comm.size);
  Matrix<Real> s_buf(4,max_num_local);

  for (int i=0; i<max_num_local; i++)
    s_buf(0,i) = -1.0;

  for (int i=0; i<local.num_local; i++)
  {
    s_buf(0,i) = Real(local.global_id[i]);
    s_buf(1,i) = local.atom[i].evec[0];
    s_buf(2,i) = local.atom[i].evec[1];
    s_buf(3,i) = local.atom[i].evec[2];
  }

  MPI_Gather(&s_buf(0,0),   4*max_num_local, MPI_DOUBLE,
             &r_buf(0,0,0), 4*max_num_local, MPI_DOUBLE, 0, comm.comm);

  if (comm.rank == 0)
  {
    for (int p=0; p<comm.size; p++)
    {
      for (int i=0; i<max_num_local && r_buf(0,i,p)>=0; i++)
      {
        int j     = r_buf(0,i,p);
        ev[3*j+0] = r_buf(1,i,p);
        ev[3*j+1] = r_buf(2,i,p);
        ev[3*j+2] = r_buf(3,i,p);
      }
    }
  }

}


void LSMS::getMag(std::vector<std::vector<Real> > &ev)
{
  Array3d<Real> r_buf(4, max_num_local, comm.size);
  Matrix<Real> s_buf(4, max_num_local);
  Real mag;

  for (int i=0; i<max_num_local; i++)
    s_buf(0,i) = -1.0;

  for (int i=0; i<local.num_local; i++)
  {
    s_buf(0,i) = Real(local.global_id[i]);
    s_buf(1,i) = local.atom[i].evec[0]*mag;
    s_buf(2,i) = local.atom[i].evec[1]*mag;
    s_buf(3,i) = local.atom[i].evec[2]*mag;
    mag        = local.atom[i].mtotws;
/*
    s_buf(1,i)=local.atom[i].dosckint[1] + local.atom[i].evec[0] * local.atom[i].mcpsc_mt;
    s_buf(2,i)=local.atom[i].dosckint[2] + local.atom[i].evec[1] * local.atom[i].mcpsc_mt;
    s_buf(3,i)=local.atom[i].dosckint[3] + local.atom[i].evec[2] * local.atom[i].mcpsc_mt;
*/

  }

  MPI_Gather(&s_buf(0,0),   4*max_num_local, MPI_DOUBLE,
             &r_buf(0,0,0), 4*max_num_local, MPI_DOUBLE, 0, comm.comm);

  if (comm.rank == 0)
  {
    for (int p=0; p<comm.size; p++)
    {
      for (int i=0; i<max_num_local && r_buf(0,i,p)>=0; i++)
      {
        int j      = r_buf(0,i,p);
        (ev[j])[0] = r_buf(1,i,p);
        (ev[j])[1] = r_buf(2,i,p);
        (ev[j])[2] = r_buf(3,i,p);
      }
    }
  }
}


void LSMS::getMag(Real *ev)
{
  Real r_buf[4*max_num_local*comm.size];
  Real mag;
  Real s_buf[4*max_num_local];

  for (int i=0; i<max_num_local; i++)
    s_buf[4*i] = -1.0;

  for (int i=0; i<local.num_local; i++)
  {
    s_buf[4*i]   = Real(local.global_id[i]);
    s_buf[1+4*i] = local.atom[i].evec[0]*mag;
    s_buf[2+4*i] = local.atom[i].evec[1]*mag;
    s_buf[3+4*i] = local.atom[i].evec[2]*mag;
    mag          = local.atom[i].mtotws;
  }

  MPI_Gather(s_buf, 4*max_num_local, MPI_DOUBLE,
             r_buf, 4*max_num_local, MPI_DOUBLE, 0, comm.comm);

  if (comm.rank == 0)
  {
    for (int p=0; p<comm.size; p++)
    {
      // for(int i=0; i<max_num_local && r_buf[4*(i+max_num_local*p)]>=0; i++)
      for (int i=0; i<max_num_local; i++)
      {
        int j = int(r_buf[4*(i+max_num_local*p)]);
        if (j >= 0)
        {
          ev[3*j+0] = r_buf[1 + 4*(i+max_num_local*p)];
          ev[3*j+1] = r_buf[2 + 4*(i+max_num_local*p)];
          ev[3*j+2] = r_buf[3 + 4*(i+max_num_local*p)];
        }
      }
    }
  }
}

void LSMS::replaceAtom(AtomData &currAtom, AtomData &newAtom) {

  // defined in a separate routine in case more care required
  currAtom = newAtom;
  currAtom.vrNew = currAtom.vr;
  currAtom.rhoNew = currAtom.rhotot;
  currAtom.resetLocalDensities(); 
}

void LSMS::setOccupancies(int *occ) {

  MPI_Status status;
 
  // distribute occupancies to respective processes
  if( comm.rank == 0 ) {

    // for every atomic site
    std::vector<MPI_Request> request(crystal.num_types);
    int n_req = 0;
    for(int p = 0; p < crystal.num_types; p++)
    {
      int owner = crystal.types[p].node;
      int l = crystal.types[p].local_id;
      int ac = crystal.types[p].alloy_class;

      // if root node owns site
      if( owner == 0 ) {

        // if occupancy changed, redefine atom by pulling from bank
        if( local.atom[l].ztotss != alloyBank[ac][occ[p]].ztotss ) {
          replaceAtom( local.atom[l], alloyBank[ac][occ[p]]);
          crystal.types[p].Z = alloyDesc[ac][occ[p]].Z;
          crystal.types[p].Zc = alloyDesc[ac][occ[p]].Zc;
          crystal.types[p].Zs = alloyDesc[ac][occ[p]].Zs;
          crystal.types[p].Zv = alloyDesc[ac][occ[p]].Zv; 
        }
      } 

      // if not owned by root, send site occupancy
      else {   
        MPI_Isend(&occ[p],1,MPI_INTEGER,owner,l,comm.comm,&request[n_req++]);
      }
    }

    // wait for all nodes to recieve
    for(int i = 0; i < n_req; i++)
      MPI_Wait(&request[i], &status);

  } 
  else {

    // if not root node, pick up local occupancies
    std::vector<MPI_Request> request(local.num_local);
    std::vector<int> xi(local.num_local);

    for(int p = 0; p < local.num_local; p++) 
      MPI_Irecv(&xi[p],1,MPI_INTEGER,0,p,comm.comm,&request[p]);

    // wait for all recieves to complete
    for(int i = 0; i < local.num_local; i++)
      MPI_Wait(&request[i], &status);

    /* dprintf("Walker %d, LSMS %d: Recieved local occupancies: ",myWalkerID,comm.rank);
    for(int p = 0; p < local.num_local; p++)
      dprintf("%d",xi[p]);
    dprintf("\n"); */

    // if a site occupancy changed, pull from alloy bank
    for(int p = 0; p < local.num_local; p++) {
      int ac = local.atom[p].alloy_class;
      // dprintf("Walker %d, LSMS %d: Local Atom %d: Alloy Class %d\n",myWalkerID,comm.rank,p,ac);
      if( local.atom[p].ztotss != alloyBank[ac][xi[p]].ztotss ) {
        replaceAtom( local.atom[p], alloyBank[ac][xi[p]] );
        int g = local.global_id[p];
        crystal.types[g].Z = alloyDesc[ac][xi[p]].Z;
        crystal.types[g].Zc = alloyDesc[ac][xi[p]].Zc;
        crystal.types[g].Zs = alloyDesc[ac][xi[p]].Zs;
        crystal.types[g].Zv = alloyDesc[ac][xi[p]].Zv;
      }
    }
  }

  // alloy bank potentiial is from a different mesh
  for (int i=0; i<local.num_local; i++) {
    if( lsms.fixRMT==0 ) {
      local.atom[i].rmt=local.atom[i].rInscribed;
      local.atom[i].generateNewMesh = true;
    }

    if(local.atom[i].generateNewMesh) 
      interpolatePotential(lsms, local.atom[i]);
  }

  // for debugging purposes, reset potentials
  for (int i=0; i<local.num_local; i++) {
    local.atom[i].vrNew = local.atom[i].vr;
    local.atom[i].rhoNew = local.atom[i].rhotot;
    local.atom[i].resetLocalDensities(); 
  }

  // re-calculate global properties: chempot, zvaltss, ...
  lsms.zvaltss=0.0;
  lsms.chempot=0.0;
  for(int i=0; i<local.num_local; i++)
  {
    lsms.zvaltss+=local.atom[i].zvalss*Real(local.n_per_type[i]);
    lsms.chempot+=local.atom[i].efermi*Real(local.n_per_type[i]);
  }
  Real fspace[2];
  fspace[0]=lsms.zvaltss;
  fspace[1]=lsms.chempot;
  
  globalSum(comm,fspace,2);
  
  lsms.zvaltss=fspace[0]; // /Real(lsms.num_atoms);
  lsms.chempot=fspace[1]/Real(lsms.num_atoms);

  // recalculate core states
  calculateCoreStates(comm,lsms,local);

  // reset mixing
  mixing->prepare(comm,lsms,local.atom);
}

void LSMS::getOccupancies(int *occ_out) {
 // not implemented
}

void LSMS::getAlloyInfo(AlloyMixingDesc &alloyDesc, int **siteclass) {

  printf("crystal.num_types = %d\n", crystal.num_types);

  alloyDesc = this->alloyDesc;
  *siteclass = (int*) malloc(sizeof(int)*crystal.num_types);

  for(int i = 0; i < crystal.num_types; i++)
    (*siteclass)[i] = crystal.types[i].alloy_class;
}


// oneStepEnergy calculates the frozen potential energy without converging the Fermi energy
Real LSMS::oneStepEnergy(Real *eb)
{
  Real eband;

  calculateCoreStates(comm,lsms,local);
  energyContourIntegration(comm,lsms,local);
  calculateChemPot(comm,lsms,local,eband);
  calculateEvec(lsms,local);
  mixEvec(lsms,local,0.0);
  calculateAllLocalChargeDensities(lsms,local);
  calculateLocalCharges(lsms, local, 0);
  checkAllLocalCharges(lsms, local);

  energyLoopCount++;

  *eb=eband;
  return eband;
}


Real LSMS::multiStepEnergy()
{
  Real eband,ef;
  int iterationCount=1;

  ef=lsms.chempot;

  if(potentialShifter.vSpinShiftFlag)
  {
    potentialShifter.applyShifts(local);
  }

  calculateCoreStates(comm,lsms,local);

  energyContourIntegration(comm,lsms,local);
  energyLoopCount++;

  calculateChemPot(comm,lsms,local,eband);
  while(std::abs(ef-lsms.chempot)>efTol && iterationCount<lsms.nscf)
  {
    ef=lsms.chempot;
    iterationCount++;
    energyContourIntegration(comm,lsms,local);
    calculateChemPot(comm,lsms,local,eband);
    energyLoopCount++;
  }

  calculateEvec(lsms,local);
  mixEvec(lsms,local,0.0);
  calculateAllLocalChargeDensities(lsms,local);
  calculateLocalCharges(lsms,local,0);
  checkAllLocalCharges(lsms, local);
  
  // calculate the Zeeman contribution from the spin shift and adjust the band energy acordingly
  if(1==0)
  if(potentialShifter.vSpinShiftFlag)
  {
    calculateEvec(lsms,local);
    mixEvec(lsms,local,0.0);
    calculateAllLocalChargeDensities(lsms,local);
    calculateLocalCharges(lsms,local,0);
    checkAllLocalCharges(lsms, local);
    Real eZeeman=0.0;
    for(int i=0; i<local.num_local; i++)
    {
      Real s=1.0;
      printf(" local.atom[%d].qvalws = %lf\n",i,local.atom[i].qvalws);
      printf(" local.atom[%d].mvalws = %lf\n",i,local.atom[i].mvalws);
      printf(" local.atom[%d].xvalws = %lf  %lf\n",i,local.atom[i].xvalwsNew[0],local.atom[i].xvalwsNew[1]);
      printf(" local.atom[%d].evec = %lf  %lf  %lf\n",i,local.atom[i].evec[0],local.atom[i].evec[1],local.atom[i].evec[2]);
      printf(" local.atom[%d].evecNew = %lf  %lf  %lf\n",i,local.atom[i].evecNew[0],local.atom[i].evecNew[1],local.atom[i].evecNew[2]);
      if (local.atom[i].spinFlipped)
        printf(" local.atom[%d].spinFlipped = true\n", i);
      else
        printf(" local.atom[%d].spinFlipped = false\n", i);
        
/*
      if(local.atom[i].evec[0]*local.atom[i].evecNew[0]+
         local.atom[i].evec[1]*local.atom[i].evecNew[1]+
         local.atom[i].evec[2]*local.atom[i].evecNew[2] < 0.0) s=-1.0;
      eZeeman+=-s*2.0*(local.atom[i].xvalwsNew[0]*local.atom[i].vSpinShift
                   +local.atom[i].xvalwsNew[1]*local.atom[i].vSpinShift);
*/
      eZeeman+=-s*(local.atom[i].xvalwsNew[0]*local.atom[i].vSpinShift
                   +local.atom[i].xvalwsNew[1]*local.atom[i].vSpinShift);
    }
    globalSum(comm,eZeeman);
    // if(lsms.global.iprint>=0)
    {

      printf("LSMS::multiStepEnergy(): eZeeman=%lf\n",eZeeman);
    }
    eband-=eZeeman;
  }

  if(lsms.global.iprint>=0)
  {
    if(iterationCount<lsms.nscf) printf("LSMS::multiStepEnergy() converged in %d steps.\n",iterationCount);
    else printf("LSMS::multiStepEnergy() did not converge in %d steps.\n",lsms.nscf);
  }
  return eband;
}


Real LSMS::scfEnergy(Real *eb)
{
  FILE *kFile = NULL;
  Real oldTotalEnergy = 0.0;
  
  Real eband;
  Real eZeeman;

  int iterationCount = 0;
  if (lsms.global.iprint >= 0)
    printf("Total number of iterations:%d\n", lsms.nscf);

  // double timeScfLoop     = MPI_Wtime();
  // double timeCalcChemPot = 0.0;

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

  if (potentialShifter.vSpinShiftFlag)
    potentialShifter.applyShifts(local);

  calculateCoreStates(comm,lsms,local);

  bool converged = false;

  // FILE *kFile = fopen("k.out","w");

  for (iterationCount = 0; iterationCount < lsms.nscf && !converged; iterationCount++)
  {
    //if (lsms.global.iprint >= 0)
    //  printf("SCF iteration %d:\n", iterationCount);
    Real rms=0.5*(local.qrms[0]+local.qrms[1]);
    // if(comm.rank==0) printf("Walker %d: On SCF iteration %d. Last RMS = %17.15f\n", myWalkerID, iterationCount, rms);

    energyContourIntegration(comm, lsms, local);

    //if(potentialShifter.vSpinShiftFlag)
    //  potentialShifter.restorePotentials(local);

    // double dTimeCCP = MPI_Wtime();
    // if(!lsms.global.checkIstop("buildKKRMatrix"))
    calculateChemPot(comm, lsms, local, eband);
    // dTimeCCP = MPI_Wtime() - dTimeCCP;
    // timeCalcChemPot += dTimeCCP;
    calculateEvec(lsms, local);
    mixEvec(lsms, local, 0.0);
    for (int i=0; i<local.num_local; i++)
    {
      local.atom[i].newConstraint();
      checkIfSpinHasFlipped(lsms, local.atom[i]);
    }

    calculateAllLocalChargeDensities(lsms, local);
    calculateChargesPotential(comm, lsms, local, crystal, 0);
    checkAllLocalCharges(lsms, local);
    calculateTotalEnergy(comm, lsms, local, crystal);

    // calculate the Zeeman contribution from the spin shift and adjust the band energy accordingly
    eZeeman = 0.0;
    if (1 == 0)
      if (potentialShifter.vSpinShiftFlag)
      {
        calculateLocalCharges(lsms,local,0);
        eZeeman = 0.0;
        for (int i=0; i<local.num_local; i++)
        {
          Real s = 1.0;
          /*
          if(local.atom[i].evec[0]*local.atom[i].evecNew[0]+
             local.atom[i].evec[1]*local.atom[i].evecNew[1]+
             local.atom[i].evec[2]*local.atom[i].evecNew[2] < 0.0) s=-1.0;
          eZeeman += -s*2.0*(local.atom[i].xvalwsNew[0]*local.atom[i].vSpinShift
                     +local.atom[i].xvalwsNew[1]*local.atom[i].vSpinShift);
          */
          eZeeman += 0.5*s*(local.atom[i].xvalwsNew[0]*local.atom[i].vSpinShift
                     +local.atom[i].xvalwsNew[1]*local.atom[i].vSpinShift);
        }
        globalSum(comm, eZeeman);
        if (lsms.global.iprint >= 0)
          printf("LSMS::scfEnergy(): eZeeman=%lf\n",eZeeman);
      }
  
    if (potentialShifter.vSpinShiftFlag)
      potentialShifter.restorePotentials(local);

    mixing -> updateChargeDensity(comm, lsms, local.atom);

    // LSMS 1: lsms_main.f:2101-2116
    for (int i=0; i<local.num_local; i++) {
      if (local.atom[i].spinFlipped)
      {   
        checkIfSpinHasFlipped(lsms, local.atom[i]);
        if (!local.atom[i].spinFlipped)
          swapCoreStateEnergies(local.atom[i]);
      }   
    }

    calculateCoreStates(comm, lsms, local);

    // If charge is mixed, recalculate the potential  (need a flag for this from input)
    calculateChargesPotential(comm,lsms,local,crystal,1);
    mixing -> updatePotential(comm,lsms,local.atom);

    if (potentialShifter.vSpinShiftFlag)
      potentialShifter.resetPotentials(local);

    rms = 0.0;
    for(int i=0; i<local.num_local; i++)
      rms = std::max(rms, 0.5*(local.atom[i].qrms[0]+local.atom[i].qrms[1]));
    globalMax(comm, rms);

    if (potentialShifter.vSpinShiftFlag)
      lsms.totalEnergy -= eZeeman;

    energyDifference = oldTotalEnergy - lsms.totalEnergy;
    if (lsms.global.iprint >= 0 && comm.rank == 0)
    {
      printf("Band Energy = %lf Ry %10s", eband, "");
      printf("Fermi Energy = %lf Ry\n", lsms.chempot);
      printf("Total Energy = %lf Ry\n", lsms.totalEnergy);
      printf("Energy Change = %lg Ry\n", energyDifference);
      printf("RMS = %lf\n",rms);
    }

    if (kFile != NULL)
      fprintf(kFile,"%3d %20.12lf %12.6lf   %12.6lf\n",iterationCount,lsms.totalEnergy,lsms.chempot,rms);

    if (potentialShifter.vSpinShiftFlag)
      potentialShifter.applyShifts(local);

    // calculate core states for new potential if we are performing scf calculations
    calculateCoreStates(comm,lsms,local);

// check for convergence
    converged = rms < lsms.rmsTolerance;
    /*
    converged = true;
    for (int i=0; i<local.num_local; i++)
    {
      converged = converged
                && (0.5*(local.qrms[0]+local.qrms[1])<lsms.rmsTolerance);
    }
    globalAnd(comm, converged);
    */

//    if(std::abs(energyDifference)<energyTolerance && rms<lsms.rmsTolerance)
//      converged=true;

    energyLoopCount++;

    oldTotalEnergy = lsms.totalEnergy;
  }

  if (potentialShifter.vSpinShiftFlag)
    potentialShifter.restorePotentials(local);

  //calculateChemPot(comm,lsms,local,*eband);
  *eb = eband;

  // if(comm.rank==0) dprintf("Walker %d: Fermi Energy = %lf Ry.\n", myWalkerID, lsms.chempot);

  return lsms.totalEnergy;
}


void LSMS::writePot(char* name)
{
  printf("*******    lsms::writePot not implemented yet! *********\n");
}

