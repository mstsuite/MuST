// main driver for LSMS_3 class

#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "VORPOL/VORPOL.hpp"
#include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "calculateChemPot.hpp"
#include "lsmsClass.hpp"
#include "ReplicaExchangeWL.hpp"
#include "EvecGenerator.h"
#include "WangLandau.h"
// #include "WangLandau_withoutKernel.h"
#include "ExhaustiveIsing.h"
#include "WangLandau2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//#define USE_PAPI 1

#ifdef USE_PAPI
#include <papi.h>
#endif

#define R_VALUE_OFFSET 2

int main(int argc, char *argv[])
{

  // MPI information
  int size, rank, world_rank, my_group {-1};
  int *lsms_rank0;            // an array to store the world_rank of the group leaders
  // size has to be size_lsms * num_lsms
  MPI_Comm local_comm;
  MPI_Status status;


  // Simulation info
  int num_lsms {1};              // number of parallel LSMS instances
  int size_lsms {-1};            // number of atoms in a lsms instance
  int num_window {1};            // number of energy windows
  int num_lsms_per_window {1};   // number of LSMS instances per windows
  int num_steps {1};             // number of energy calculations
  int MCstepCount {0};           // count the Monte Carlo steps executed
  int REstepCount {0};           // count the replica exchange steps executed
  int exchangeFrequency {10};    // number of MC steps between replica exchange

  double max_time;               // maximum walltime for this run in seconds
  bool restrict_time {false};    // was the maximum time specified?
  bool restrict_steps {false};   // or the max. numer of steps?

  double walltime_0, walltime;
  double restartWriteFrequency = 30.0 * 60.0;
  double nextWriteTime = restartWriteFrequency;

  double magnetization {0.0};
  double energy_accumulator {0.0};    // accumulates the energy to calculate the mean
  int energies_accumulated {0};

  typedef enum {Constant, Random, WangLandau_1d, ExhaustiveIsing, WangLandau_2d} EvecGenerationMode;
  typedef enum {MagneticMoment, MagneticMomentZ, MagneticMomentX, MagneticMomentY} SecondDimension;

  EvecGenerationMode evecGenerationMode = Constant;
  SecondDimension second_dimension = MagneticMoment;
  double ev0[3] {0.0, 0.0, 1.0};

  //Ying Wai: should return_moments_flag be initialized to false?
  bool return_moments_flag {true}; // true -> return all magnetic moments from lsms run at each step.
  bool generator_needs_moment {false};

  typedef enum {OneStepEnergy, MultiStepEnergy, ScfEnergy} EnergyCalculationMode;
  EnergyCalculationMode energyCalculationMode = OneStepEnergy;
  int energyIndex {1}; // index for the return value to use for the MC step (0: total energy, 1: band energy)


  // Input/output filenames
  char prefix[40] {};
  char i_lsms_name[64];
  char gWL_in_name[64], gWL_out_name[64];
  char mode_name[64];
  char energy_calculation_name[64];
  char stupid[37];

  char step_out_name[64];
  char wl_step_out_name[128];
  char *wl_stepf = nullptr;
  bool step_out_flag = false;
  std::ofstream step_out_file;

  sprintf(i_lsms_name, "i_lsms");
  gWL_in_name[0] = gWL_out_name[0] = 0;
  mode_name[0] = 0;
  energy_calculation_name[0] = 0;


  // check command line arguments
  for (int i=0; i<argc; i++)
  {
    //if (!strcmp("-num_lsms", argv[i])) num_lsms = atoi(argv[++i]);
    if (!strcmp("-size_lsms", argv[i])) size_lsms = atoi(argv[++i]);
    if (!strcmp("-num_window", argv[i])) num_window = atoi(argv[++i]);
    if (!strcmp("-num_lsms_per_window", argv[i])) num_lsms_per_window = atoi(argv[++i]);
    if (!strcmp("-num_steps", argv[i])) {
      num_steps = atoi(argv[++i]);
      restrict_steps = true;
    }
    if (!strcmp("-exchange_frequency", argv[i])) exchangeFrequency = atoi(argv[++i]);
    if (!strcmp("-walltime", argv[i])) {
      max_time = 60.0 * atof(argv[++i]);
      restrict_time = true;
    }
    if (!strcmp("-i", argv[i])) strncpy(i_lsms_name, argv[++i], 64);
    if (!strcmp("-random_dir", argv[i])) evecGenerationMode = Random;
    if (!strcmp("-step_out", argv[i])) {
      strncpy(step_out_name, argv[++i],64);
      step_out_flag = true;
      return_moments_flag = true;
    }
    if (!strcmp("-wl_out", argv[i])) strncpy(gWL_out_name, argv[++i], 64);
    if (!strcmp("-wl_in", argv[i])) strncpy(gWL_in_name, argv[++i], 64);
    if (!strcmp("-mode", argv[i])) strncpy(mode_name, argv[++i], 64);
    if (!strcmp("-energy_calculation", argv[i])) strncpy(energy_calculation_name, argv[++i], 64);
  }

  num_lsms = num_window * num_lsms_per_window;
  lsms_rank0 = (int *) malloc(sizeof(int) * (num_lsms+1));

  if (!(restrict_steps || restrict_time)) restrict_steps = true;

  if (mode_name[0] != 0)
  {
    if (!strcmp("constant", mode_name)) evecGenerationMode = Constant;
    if (!strcmp("random", mode_name)) evecGenerationMode = Random;
    if (!strcmp("1d", mode_name)) evecGenerationMode = WangLandau_1d;
    if (!strcmp("ising", mode_name)) evecGenerationMode = ExhaustiveIsing;
    if (!strcmp("2d", mode_name)) evecGenerationMode = WangLandau_2d;
    if (!strcmp("2d-m", mode_name)) 
    {
      evecGenerationMode = WangLandau_2d;
      second_dimension = MagneticMoment;
    }
    if (!strcmp("2d-x", mode_name)) 
    {
      evecGenerationMode = WangLandau_2d;
      second_dimension = MagneticMomentX;
    }
    if (!strcmp("2d-y", mode_name)) 
    {
      evecGenerationMode = WangLandau_2d;
      second_dimension = MagneticMomentY;
    }
    if (!strcmp("2d-z", mode_name)) 
    {
      evecGenerationMode = WangLandau_2d;
      second_dimension = MagneticMomentZ;
    }
  }

  // make sure 'return_moments_flag' is set correctly
  switch (evecGenerationMode)
  {
    case Constant : break;
    case Random : break;
    case WangLandau_1d :
      generator_needs_moment = true;
      return_moments_flag = true;
      break;
    case ExhaustiveIsing : break;
    case WangLandau_2d :
      generator_needs_moment = true;
      return_moments_flag = true;
      break;
    default:
      std::cout << " ERROR: UNKNOWN EVEC GENERATION MODE\n";
      exit(1);
  }
  //Ying Wai: this line is redundant
  //if (generator_needs_moment) return_moments_flag = true;

  if (energy_calculation_name[0] != 0)
  {
    if(energy_calculation_name[0] == 'o') 
    {
      energyCalculationMode = OneStepEnergy;
      energyIndex = 1; 
    }
    if(energy_calculation_name[0] == 'm')
    {
      energyCalculationMode = MultiStepEnergy;
      energyIndex = 1;
    }
    if(energy_calculation_name[0] == 's')
    {
      energyCalculationMode = ScfEnergy;
      energyIndex = 0;
    }
  }

#ifdef USE_PAPI
#define NUM_PAPI_EVENTS 4
  int hw_counters = PAPI_num_counters();
  if(hw_counters>NUM_PAPI_EVENTS) hw_counters=NUM_PAPI_EVENTS;
  int papi_events[NUM_PAPI_EVENTS]; // = {PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_FP_OPS,PAPI_VEC_INS};
  char *papi_event_name[] = {"PAPI_TOT_INS","PAPI_FP_OPS",
                             "RETIRED_SSE_OPERATIONS:DOUBLE_ADD_SUB_OPS:DOUBLE_MUL_OPS:DOUBLE_DIV_OPS:OP_TYPE",
                             "RETIRED_SSE_OPERATIONS:SINGLE_ADD_SUB_OPS:SINGLE_MUL_OPS:SINGLE_DIV_OPS:OP_TYPE"};
  // "RETIRED_INSTRUCTIONS",
  // "RETIRED_MMX_AND_FP_INSTRUCTIONS:PACKED_SSE_AND_SSE2",
  // "RETIRED_SSE_OPERATIONS:DOUBLE_ADD_SUB_OPS:DOUBLE_MUL_OPS:DOUBLE_DIV_OPS:1",
  // "RETIRED_SSE_OPERATIONS:SINGLE_ADD_SUB_OPS:SINGLE_MUL_OPS:SINGLE_DIV_OPS:1"
  // get events from names:
  for(int i=0; i<NUM_PAPI_EVENTS; i++)
  {
    if(PAPI_event_name_to_code(papi_event_name[i],&papi_events[i]) != PAPI_OK)
    {
      // printline("Error in obtaining PAPI event code for: "+ttos(papi_event_name[i]),
      //           std::cerr,parameters.myrankWorld);
      // printline("Skipping all following events",
      //           std::cerr,parameters.myrankWorld);
      if(hw_counters>i) hw_counters=i;
    }
  }
  long long papi_values[NUM_PAPI_EVENTS+4];
  // printline("PAPI: "+ttos(hw_counters)+" counters available",std::cout,parameters.myrankWorld);
  if(hw_counters>NUM_PAPI_EVENTS) hw_counters=NUM_PAPI_EVENTS;
  long long papi_real_cyc_0 = PAPI_get_real_cyc();
  long long papi_real_usec_0 = PAPI_get_real_usec();
  long long papi_virt_cyc_0 = PAPI_get_virt_cyc();
  long long papi_virt_usec_0 = PAPI_get_virt_usec();
  PAPI_start_counters(papi_events,hw_counters);
#endif

#ifndef SVN_REV
#define SVN_REV "unknown"
#endif


  // initialize MPI:
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  world_rank = rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  walltime_0 = MPI_Wtime();

  // Rank 0 prints out simulation info 
  if(world_rank == 0)
  {
    std::cout << "LSMS_3" << std::endl;
    std::cout << " SVN revision " << SVN_REV << std::endl << std::endl;
#ifdef USE_PAPI
    std::cout << " Using Papi counters" << std::endl << std::endl; 
#endif
    std::cout << " Size of LSMS instances = " << size_lsms << " atoms\n";
    std::cout << " Number of LSMS instances = " << num_lsms << std::endl;
    std::cout << " LSMS Energy calculated using ";

    switch (energyCalculationMode)
    {
      case OneStepEnergy:
        std::cout << "oneStepEnergy [frozen potential band energy]" << std::endl;
        break;
      case MultiStepEnergy:
        std::cout << "multiStepEnergy [frozen potential band energy with converged Fermi energy]" << std::endl;
        break;
      case ScfEnergy:
        std::cout << "scfEnergy [self-consistent total energy]" << std::endl;
        break;
      default:
        std::cout << "UNKNOWN ENERGY CALCULATION METHOD" << std::endl;
        exit(1);
    }

    if (restrict_steps)
      std::cout << " Number of gWL steps = " << num_steps << std::endl;
    if (restrict_time)
      std::cout << " Maximum walltime = " << max_time << "s\n";

    //std::cout << " Processor alignment (process allocation quantization) = " << align << std::endl;

    switch (evecGenerationMode)
    {
      case Constant :
        std::cout << " Constant moments direction along "
                  << ev0[0] << " " << ev0[1] << " " << ev0[2] << std::endl;
        break;
      case Random :
        std::cout << " Random distribution of moments (no Wang-Landau)" << std::endl;
        break;
      case WangLandau_1d :
        std::cout << " Wang-Landau for one continuous variable (energy)" << std::endl;
        break;
      case ExhaustiveIsing :
        std::cout << " Exhaustive Ising sampling" << std::endl;
        break;
      case WangLandau_2d :
        std::cout << " Wang-Landau for two continuous variable (energy, ";
        switch (second_dimension)
        {
          case MagneticMoment  : std::cout << "magnitude of magnetization)"; break;
          case MagneticMomentX : std::cout<<"x component of magnetization)"; break;
          case MagneticMomentY : std::cout<<"y component of magnetization)"; break;
          case MagneticMomentZ : std::cout<<"z component of magnetization)"; break;
        }
        std::cout << std::endl;
        break;
      default:
        std::cout << " ERROR: UNKNOWN EVEC GENERATION MODE\n";
        exit(1);
    }

    if (step_out_flag)
      std::cout << " Step output written to: " << step_out_name << std::endl;
    std::cout << std::endl;

    if (step_out_flag && (evecGenerationMode == WangLandau_1d))
    {
      // step_out_flag = false;
      snprintf(wl_step_out_name, 127, "wl1d_%s", step_out_name);
      wl_stepf = wl_step_out_name;
    }

    if (step_out_flag)
    {
      step_out_file.open(step_out_name);
      step_out_file << "#";
      for(int i=0; i<argc; i++) step_out_file << " " << argv[i];
      step_out_file << std::endl << size_lsms << std::endl;
    }
  }


  if (num_lsms == 1)
  {
    LSMS lsms_calc(MPI_COMM_WORLD, i_lsms_name, "1_");
      
    if (world_rank == 0)
    {
      std::cout << "executing LSMS(C++) for " << lsms_calc.numSpins() << " atoms\n";
      std::cout << "  LSMS version = " << lsms_calc.version() << std::endl;
    }

    if (energyCalculationMode == OneStepEnergy)
      std::cout << "one step Energy = " << lsms_calc.oneStepEnergy() << std::endl;
    else if (energyCalculationMode == MultiStepEnergy)
      std::cout << "multi-step Energy = " << lsms_calc.multiStepEnergy() << std::endl;
    else if (energyCalculationMode == ScfEnergy)
      std::cout << "self-consistent Energy = " << lsms_calc.scfEnergy() << std::endl;
    else
    {
      std:: cout << "ERROR: Unknown energy calculation mode for lsms_calc in wl-lsms main!\n";
      MPI_Abort(MPI_COMM_WORLD, 5);
    }
  }
  else
  {
    // build the LSMS communicators
    int comm_size = (size - (size % num_lsms)) / num_lsms;
    my_group = (world_rank - (world_rank % comm_size)) / comm_size;
    for (int i=0; i<num_lsms; i++)
      lsms_rank0[i] = i * comm_size;

    MPI_Comm_split(MPI_COMM_WORLD, my_group, 0, &local_comm);
    MPI_Comm_rank(local_comm, &rank);

    //std::cout << "world_rank " << world_rank << " -> group " << my_group 
    //          << " (local_rank " << rank << " ) \n";

    // build the WL communicators
    MPI_Comm WLwalkersComm;
    MPI_Group MPI_GROUP_WORLD, WLwalkersGroup;
    
    MPI_Comm_group (MPI_COMM_WORLD, &MPI_GROUP_WORLD);
    MPI_Group_incl (MPI_GROUP_WORLD, num_lsms, lsms_rank0, &WLwalkersGroup);
    MPI_Comm_create (MPI_COMM_WORLD, WLwalkersGroup, &WLwalkersComm);

    // Prepare seeds for RNG
    std::seed_seq seq {MPI_Wtime()};
    std::vector<unsigned> seeds(num_lsms);
    seq.generate(seeds.begin(), seeds.end());

    // Initialize REWL class
    REWL *rewl {};
    if (rank == 0) 
      rewl = new REWL(WLwalkersComm, num_window, num_lsms_per_window, seeds[my_group]);

    // now we get ready to do some calculations...

    double energy {0.0};
    double bandEnergy {0.0};  
    double** evecs {};
    int op {5};
    //int i_values[10];

    evecs = (double **) malloc(sizeof(double *));
    evecs[0] = (double *) malloc(sizeof(double) * 3 * size_lsms);
    for (int i = 0; i < 3*size_lsms; i++)
      evecs[0][i] = 0.0;

    snprintf(prefix, 38, "%d_", my_group);
    LSMS lsms_calc(local_comm, i_lsms_name, prefix);

    char *wl_inf {};
    char *wl_outf {};
    // YingWai: To read from the corresponding WLrestart file, 
    // the name of wl_inf needs to be reconstructed for different walker (using my_group)
    // Apr 30, 14
    if (gWL_in_name) wl_inf = gWL_in_name;
    if (gWL_out_name) wl_outf = gWL_out_name;
  
    EvecGenerator *generator {};
  
    if (world_rank == 0)
    {
      snprintf(prefix, 38, "Group %4d: ", my_group);
      std::cout << prefix << "executing LSMS(C++) for " << lsms_calc.numSpins() << " atoms\n"    ;
      std::cout << prefix << "  LSMS version = " << lsms_calc.version() << std::endl;
    }
  
    if (rank == 0)
    {
      // Initialize the correct evec generator
      switch (evecGenerationMode)
      {   
        case Random :
             generator = new RandomEvecGenerator(size_lsms);
             break;
        case Constant :
             generator = new ConstantEvecGenerator(size_lsms, ev0, num_lsms);
             break;
        case WangLandau_1d :
             generator = new WL1dEvecGenerator<std::mt19937>(size_lsms, 1, num_window, num_lsms_per_window, my_group, evecs, wl_inf, wl_outf, wl_stepf);
             break;
        case ExhaustiveIsing :
             generator = new ExhaustiveIsing1dEvecGenerator(size_lsms, num_lsms, evecs, wl_inf, wl_outf);
             break;
        case WangLandau_2d :
             generator = new WL2dEvecGenerator<std::mt19937>(size_lsms, num_lsms, evecs, wl_inf, wl_outf, wl_stepf);
             break;
        default :
             std::cerr << "REWL Main :: The code should never arrive here: UNKNOWN EVEC GENERATION MODE\n";
             exit(1);
      }
  
      // Generate the initial configuration
      std::cout << "REWL Main :: starting main calculation in group " << my_group << std::endl;    
  
      generator -> initializeEvec(my_group, evecs[0]);

      // Ying Wai: do we want to put it here?
      // see if the energy of the state is within the energy range
      // if not, generate a new state until it is       (Feb 4, 14)

      generator -> startSampling();
    
    }
    
    op = 5;
    MPI_Bcast(&op, 1, MPI_INT, 0, local_comm);
    /* 
       Recognized opcodes:
       5: calculate energy
          recognized energy calculation modes:
           OneStepEnergy :
               calclulate frozen potential band energy in one step (don't converge Ef)
               use only if the Fermi energy will not change due to MC steps!
               The only method available in LSMS_1.9
           MultiStepEnergy :
               calculate frozen potential band energy after converging Fermi energy
               This should be the new default method if the Fermi energy doesn't change
               multiStepEnergy only performs one step and should be equivalent to oneStepEnergy
               The tolerance for Ef convergence can be set with LSMS::setEfTol(Real).
               The default tolerance is set in the LSMS::LSMS constructor (currently 1.0e-6).
               The maximum number of steps is read from the LSMS input file 'nscf' parameter.
           ScfEnergy :
               calculate the selfconsistent total energy.
               The maximum number of steps is read from the LSMS input file 'nscf' parameter.
        
       10: get number of sites
    */
  
    bool more_work = true;
    while (more_work)
    {
 
      switch (op)
      {
        case 5 :
          
          // rank 0 sends evecs to the rest of the group
          lsms_calc.setEvec(evecs[0]);
  
          // calculate energy collectively
          if (energyCalculationMode == OneStepEnergy)
            energy = lsms_calc.oneStepEnergy(&bandEnergy);
          else if (energyCalculationMode == MultiStepEnergy)
            bandEnergy = energy = lsms_calc.multiStepEnergy();
          else if (energyCalculationMode == ScfEnergy)
            energy = lsms_calc.scfEnergy(&bandEnergy);
          else
          {
            std::cout << "ERROR: Unknown energy calculation mode for lsms_calc in wl-lsms main!\n";
            MPI_Abort(MPI_COMM_WORLD, 5);
          }
      
          // In case the magnetization has changed, get the new one back
          if (return_moments_flag)
            lsms_calc.getMag(evecs[0]);
   
          break;
  
        case 10 :
          //  i_values[0] = lsms_calc.numSpins();
          //  MPI_Send(i_values, 10, MPI_INT, 0, 1010, MPI_COMM_WORLD);  
          break;
  
        default :
          // printf("world rank %d: received exit\n",world_rank); 
          more_work = false;
         
      }
  
      //====== Up to this point, energy, bandEnergy and evecs are updated ======//
  
      if (rank == 0) 
      {
        // newConfigAcceptance:
        // true  - the proposed configuration is accepted
        // false - the proposed configuration is rejected, update with the old config.
        bool newConfigAcceptance {false};

        // these are for checking and debugging, not for real calculations 
        energy_accumulator += energy;
        energies_accumulated++;
  
        if (generator_needs_moment)
        {
          double m0 {0.0}, m1 {0.0}, m2 {0.0};
          for (int i = 0; i < 3*size_lsms; i += 3)
          {
            m0 += evecs[0][i];
            m1 += evecs[0][i+1];
            m2 += evecs[0][i+2];
          }
          switch (second_dimension)
          {
            case  MagneticMoment : magnetization = std::sqrt(m0*m0 + m1*m1 + m2*m2); break;
            case  MagneticMomentX : magnetization = m0; break;
            case  MagneticMomentY : magnetization = m1; break;
            case  MagneticMomentZ : magnetization = m2; break;
          }

          // Determine if the configuration is accepted
          newConfigAcceptance = generator -> determineAcceptance(0, energy, magnetization);
          
          // Update histogram, change gamma and check flatness
          //if (generator -> updateHistogram(0, evecs[0], newConfigAcceptance))
          //  more_work = false;
        }
        else 
        {
          // Determine if the configuration is accepted
          newConfigAcceptance = generator -> determineAcceptance(0, energy);

          // Update histogram, change gamma and check flatness
          //if (generator -> updateHistogram(0, evecs[0], newConfigAcceptance))
          //  more_work = false;
        }
/*
        //YingWai's check on the evecs
        FILE* yingwai;
        char checkFile[50];
        sprintf(checkFile, "checkEvec%05d.dat", my_group);
        yingwai = fopen(checkFile, "a");
        fprintf(yingwai, "Acceptance = ");
        if (newConfigAcceptance)
          fprintf(yingwai, "Y\n");
        else
          fprintf(yingwai, "N\n");
        for (int i = 0; i < 3*size_lsms; i++)
          fprintf(yingwai, "%15.8f\n", evecs[0][i]);
        fclose(yingwai);
*/
        // Update histogram, change gamma and check flatness
        if (generator -> updateHistogram(0, evecs[0], newConfigAcceptance))
          more_work = false;
        
        MCstepCount++;

        // Replica-exchange
        MPI_Barrier(WLwalkersComm);       // necessary?

        bool replicaExchangeAcceptance {false};
        if (MCstepCount % exchangeFrequency == 0) 
        {
          Real myDOSRatio {0.0};

          Real energyTemp = energy;
          rewl -> assignSwapPartner();
          rewl -> swapEnergy(energyTemp);

          // Determines if the energy received is within my energy range
          // and get the DOS ratio
          myDOSRatio = generator -> getDOSRatio(0, energyTemp);

          // Exchange DOS ratios and calculate acceptance probability
          replicaExchangeAcceptance = rewl -> determineAcceptance(myDOSRatio);

          // Exchange configurations
          if (replicaExchangeAcceptance)
          {
            energy = energyTemp;

            //printf("YW's check. Walker %5d: Before exchange.\n", my_group);
            //for (int yw=0; yw<3*size_lsms; yw++)
            //  printf("YW. Walker %5d: i = %d, evecs = %10.5e\n", my_group, yw, evecs[0][yw]);
            rewl -> swapConfig(evecs[0], 3*size_lsms);

            // calculate magnetization   (or do it before determineAcceptance?)
            if (generator_needs_moment)
            {
              double m0 {0.0}, m1 {0.0}, m2 {0.0};
              for(int i = 0; i < 3*size_lsms; i += 3)
              {
                m0 += evecs[0][i];
                m1 += evecs[0][i+1];
                m2 += evecs[0][i+2];
              }
              switch (second_dimension)
              {
                case  MagneticMoment : magnetization = std::sqrt(m0*m0 + m1*m1 + m2*m2); break;
                case  MagneticMomentX : magnetization = m0; break;
                case  MagneticMomentY : magnetization = m1; break;
                case  MagneticMomentZ : magnetization = m2; break;
              }
              if (generator -> updateHistogramFromRE(0, evecs[0], energy, magnetization, my_group))
                more_work = false;
            }
            else {
              if (generator -> updateHistogramFromRE(0, evecs[0], energy, my_group))
                more_work = false;
            }
          }
          REstepCount++;
        }

        // Prepare a new configuration
        if (MCstepCount % exchangeFrequency == 0 && replicaExchangeAcceptance) 
          generator -> generateEvec(0, evecs[0], replicaExchangeAcceptance);
        else
          generator -> generateEvec(0, evecs[0], newConfigAcceptance);
        num_steps -= 1;

        if (restrict_steps && num_steps <= 0) more_work = false;
        if (restrict_steps) std::cout << "      " << num_steps << " steps remaining\n";
        walltime = MPI_Wtime() - walltime_0;
        if (restrict_time && walltime >= max_time) more_work = false;
        if (restrict_time) std::cout << "      " << max_time - walltime << " seconds remaining\n";
 
        // Determine if more work will be done; this will be broadcasted to all later
        if (!more_work) op = 2;
      
        if (step_out_flag && newConfigAcceptance)
        {
          step_out_file << "# iteration " << energies_accumulated << std::endl;
          step_out_file.precision(15);
          step_out_file << energies_accumulated << std::endl;
          step_out_file << energy << "  " << bandEnergy << std::endl;
          for (int j=0; j<3*size_lsms; j+=3)
            step_out_file << evecs[0][j] << "  " << evecs[0][j+1] << "  " 
                          << evecs[0][j+2] << std::endl;
        }
  
        // write restart file every restartWriteFrequency seconds
        if (walltime > nextWriteTime)
        {
          if (rank == 0)
          {
            char restartFile[51];
            sprintf(restartFile, "WLrestart%05d.jsn", my_group);
            generator -> writeState(restartFile);
            //generator -> writeState("WLrestart.jsn");
            nextWriteTime += restartWriteFrequency;
          }
        }
 
      }
 
      MPI_Bcast(&op, 1, MPI_INT, 0, local_comm);
      if (op == 2) more_work = false;
  
    }

    // YingWai: or can be replaced by a collective call?        (Dec 18, 13)
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
      char restartFile[51];
      sprintf(restartFile, "WLrestart%05d.jsn", my_group);
      generator -> writeState(restartFile);

      delete generator;
      MPI_Comm_free(&WLwalkersComm);
      MPI_Group_free(&WLwalkersGroup);
    }
      
    free(evecs[0]);
    free(evecs);
  
    // Freeing WL communication related comm.
    //MPI_Comm_free(&WLwalkersComm);
    //MPI_Group_free(&WLwalkersGroup);
    MPI_Group_free(&MPI_GROUP_WORLD);
  }

  // Wrapping up:

  if (num_lsms > 1)
  {
    // make sure everyone arrives here:
    MPI_Bcast (stupid, 37, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Comm_free(&local_comm);
  }

  if (world_rank == 0)
  {
    if (step_out_flag)
    {
      step_out_file << "# end\n-1\n"
                    << energy_accumulator / double(energies_accumulated) << std::endl;
      step_out_file.close();
    }
    std::cout << "Finished all scheduled calculations. Freeing resources.\n";
    std::cout<<"Energy mean = "<<energy_accumulator/double(energies_accumulated)<<"Ry\n";

    double walltime = MPI_Wtime() - walltime_0;
    std::cout << " WL-LSMS finished in " << walltime << " seconds.\n";
    std::cout << " Monte-Carlo steps / walltime = "
              << double(MCstepCount) / walltime << "/sec\n";
  }

#ifdef USE_PAPI
  PAPI_stop_counters(papi_values,hw_counters);
  papi_values[hw_counters  ] = PAPI_get_real_cyc()-papi_real_cyc_0;    // real time counter value in clock cycles
  papi_values[hw_counters+1] = PAPI_get_real_usec()-papi_real_usec_0;  // real time counter value in microseconds
  papi_values[hw_counters+2] = PAPI_get_virt_cyc()-papi_virt_cyc_0;    // virtual time counter value in clock cycles
  papi_values[hw_counters+3] = PAPI_get_virt_usec()-papi_virt_usec_0;  // virtual time counter values in microseconds
  long long accumulated_counters[NUM_PAPI_EVENTS+4];
/*
  for(int i=0; i<hw_counters; i++)
  {
  printline(ttos(papi_event_name[i])+" = "+ttos(papi_values[i]),
  std::cout,parameters.myrankWorld);
  }
  printline("PAPI real cycles : "+ttos(papi_values[hw_counters]),
  std::cout,parameters.myrankWorld);
  printline("PAPI real usecs : "+ttos(papi_values[hw_counters+1]),
  std::cout,parameters.myrankWorld);
  printline("PAPI user cycles : "+ttos(papi_values[hw_counters+2]),
  std::cout,parameters.myrankWorld);
  printline("PAPI user usecs : "+ttos(papi_values[hw_counters+3]),
  std::cout,parameters.myrankWorld);
*/
  MPI_Reduce(papi_values,accumulated_counters,hw_counters+4,
             MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  if(world_rank==0)
  {
    for(int i=0; i<hw_counters; i++)
    {
      std::cout<<"Accumulated: "<<(papi_event_name[i])<<" = "<<(accumulated_counters[i])<<"\n";
    }
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
  }
#endif

  MPI_Finalize();
  return 0;
}
