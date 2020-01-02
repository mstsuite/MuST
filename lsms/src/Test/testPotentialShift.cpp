// driver for LSMS_3 class to test Potential shift

#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <fenv.h>

#include <iostream>
#include <fstream>
#include "../Main/SystemParameters.hpp"
#include "../Main/PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/CoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "VORPOL/VORPOL.hpp"
// #include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
// #include "calculateChemPot.hpp"
#include "../Main/lsmsClass.hpp"
#include "../Potential/PotentialShifter.hpp"

#define R_VALUE_OFFSET 2

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
  int size, rank, worldRank;
  int size_lsms; // number of atoms in a lsms instance
  int num_steps; // number of energy calculations
  double shift_min, shift_max;

  double magnetization;
  double magnetization0;

  double walltime_0,walltime;

  MPI_Status status;

  char i_lsms_name[64];

  char step_out_name[64];
  FILE *stepOutFile;

  double ev0[3];

  std::vector<Matrix<Real> > savedPotentials;

  bool return_moments_flag=false; // true-> return all magnetic moments from lsms run at each step.
  bool generator_needs_moment=false;
  bool reset_potentials=false;

  typedef enum {OneStepEnergy, MultiStepEnergy, ScfEnergy} EnergyCalculationMode;
  EnergyCalculationMode energyCalculationMode = MultiStepEnergy;

  // feenableexcept(FE_INVALID);

  ev0[0]=ev0[1]=0.0; ev0[2]=1.0;
  // size has to be align + size_lsms*num_lsms
  size_lsms=-1;
  num_steps=1;

  shift_min=0.0;
  shift_max=0.0;

  sprintf(i_lsms_name,"i_lsms");
  sprintf(step_out_name,"step.out");

  // check command line arguments
  for(int i=0; i<argc; i++)
    {
      if(!strcmp("-num_steps",argv[i])) {num_steps=atoi(argv[++i]);}
      if(!strcmp("-shift_min",argv[i])) {shift_min=atof(argv[++i]);}
      if(!strcmp("-shift_max",argv[i])) {shift_max=atof(argv[++i]);}
      if(!strcmp("-i",argv[i])) strncpy(i_lsms_name,argv[++i],64);
      if(!strcmp("-step_out",argv[i]))
                 {strncpy(step_out_name,argv[++i],64);}
      if(!strcmp("-reset_potentials",argv[i])) {reset_potentials=true;}
    }

  // initialize MPI:
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  worldRank=rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  walltime_0 = MPI_Wtime();

#ifndef SVN_REV
#define SVN_REV "unknown"
#endif

 
  double energy;
  double bandEnergy;
  double *evec,*m;
  int op;

  //YingWai's trial
  //PotentialShifter potentialShifter;
  //LSMS lsms_calc(MPI_COMM_WORLD, i_lsms_name, "", potentialShifter);
  LSMS lsms_calc(MPI_COMM_WORLD, i_lsms_name, "");

  size_lsms=lsms_calc.numSpins();

  lsms_calc.setEnergyTol(1.0e-6);
  lsms_calc.setRmsTol(1.0e-5);

 if(rank==0)
  {
    std::cout<<"LSMS_3"<<std::endl;
    std::cout<<" SVN revision "<<SVN_REV<<std::endl<<std::endl;
    std::cout<<"executing LSMS(C++) for "<<lsms_calc.numSpins()<<" atoms\n";
    std::cout<<"  LSMS version = "<<lsms_calc.version()<<std::endl;
    std::cout<<" Number of steps = "<<num_steps<<std::endl;
    std::cout<<" Step output written to: "<<step_out_name<<std::endl;
    if(reset_potentials) std::cout<<"Potentials will be reset at each step."<<std::endl;
    std::cout<<std::endl;


    stepOutFile=fopen(step_out_name,"w");
  }

  lsms_calc.savePotentials(savedPotentials);

  evec=(double *)malloc(sizeof(double)*3*size_lsms);
  m=(double *)malloc(sizeof(double)*3*size_lsms);

// loop over shifts:
  double shift_delta=shift_min;
  if(num_steps>1)
    shift_delta=(shift_max-shift_min)/(double(num_steps-1));

// calculated unconstraine magnetization:
  // initialize evec
  for(int i=0; i<3*size_lsms; i+=3)
  { evec[i]=0.0; evec[i+1]=0.0;  evec[i+2]=1.0;}
  lsms_calc.setEvecAndSpinPotentialShift(evec);

  // bandEnergy=energy=lsms_calc.multiStepEnergy();
  energy=lsms_calc.scfEnergy(&bandEnergy);

  lsms_calc.getMag(m);
  if(rank==0)
  {
    magnetization0=0.0;
    for(int i=0; i<3*size_lsms; i+=3)
    {
      magnetization0+=std::sqrt(m[i]*m[i]+m[i+1]*m[i+1]+m[i+2]*m[i+2]);
    }
  }


  for(int step=0; step<num_steps; step++)
  {
    double shift=shift_min+double(step)*shift_delta;

    printf("========================================\n");
    printf("  Step #%d : shift=%6.4lf\n",step,shift);
// initialize evec
    for(int i=0; i<3*size_lsms; i+=3)
    { evec[i]=0.0; evec[i+1]=0.0;  evec[i+2]=1.0+shift;}

    if(reset_potentials) lsms_calc.restorePotentials(savedPotentials);
    lsms_calc.setEvecAndSpinPotentialShift(evec);

    // energy=lsms_calc.oneStepEnergy(&bandEnergy);
    // bandEnergy=energy=lsms_calc.multiStepEnergy();
    energy=lsms_calc.scfEnergy(&bandEnergy);

    lsms_calc.getMag(m);
    if(rank==0)
    {
      magnetization=0.0;
      for(int i=0; i<3*size_lsms; i+=3)
      {
        magnetization+=std::sqrt(m[i]*m[i]+m[i+1]*m[i+1]+m[i+2]*m[i+2]);
      }

// for frozen potential a correction in the energy is needed:
      if(1==0)
        energy=energy-shift*(magnetization0-magnetization);

      magnetization=magnetization/double(size_lsms);
      energy=energy/double(size_lsms);
      fprintf(stepOutFile,"%4d %6.4lf %15.10lf %15.10lf %15.10lf %15.10lf\n",step,shift,energy,magnetization,lsms_calc.getEf(),
              lsms_calc.energyDifference);
      fflush(stepOutFile);
      printf("  Energy = %15.10lf Moment=%15.10lf Fermi Energy = %15.10lf\n",energy,magnetization,lsms_calc.getEf());
    }
  }

  free(evec);
  free(m);

  if(rank==0)
  {
    fclose(stepOutFile);
    std::cout<<" LSMS finished in "<<MPI_Wtime() - walltime_0<<" seconds.\n";
  }

  MPI_Finalize();
  return 0;
}
