// -*- mode: c++ -*-
#ifndef LSMS_WANG_LANDAU_H
#define LSMS_WANG_LANDAU_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
// we use BOOST for the random number generator
// #include <boost/random.hpp>
#include <random>
#include <ctime>
#include "../../mjson/json.h"
#include "EvecGenerator.h"
#include "Graph1dMoments.hpp"
#include "../Potential/PotentialShifter.hpp"

void inline performGlobalUpdate(Graph1dMoments<double,double> &g, double kappa, double lambda, double omega)
{
  for(int i=0; i<g.getN(); i++)
    if(g[i]>omega) g[i]+=kappa*std::exp(-lambda/(g[i]-omega));
}

class StatesWriter
{
public:
  StatesWriter(const char *filename=NULL)
  {
    if(filename==NULL) writeFlag=false;
    else {writeFlag=true; of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(17);}
  }
  ~StatesWriter() {if(writeFlag) of.close();}
  void writeHeader(double lnF, int numWalk, int numSpin, double **spins)
  {
    if(writeFlag)
      {
	of<<lnF<<" "<<numWalk<<" "<<numSpin<<std::endl;
	for(int i=0; i<numWalk; i++)
	  {
	    of<<i;
	    for(int j=0; j<3*numSpin; j++) of<<" "<<spins[i][j];
	    of<<std::endl;
	  }
      }
  }
  void writeChange(int iWalk, int numRet, int ispin, double *ev, double E)
  {
    if(writeFlag)
      of<<iWalk<<" "<<numRet<<" "<<ispin<<" "<<ev[0]<<" "<<ev[1]<<" "<<ev[2]<<" "<<E<<std::endl;
  }
  void newFile(const char *filename=NULL)
  {
    if(writeFlag) of.close();
    writeFlag=false;
    if(filename!=NULL){
      writeFlag=true;
      of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(8);
    }
  }
private:
  bool writeFlag;
  std::ofstream of;
};

template<class RNG = std::mt19937>
// template<class RNG = boost::mt19937>
class WL1dEvecGenerator : public EvecGenerator
{
 public:

  WL1dEvecGenerator(int num_spins, int num_instances, int numWindows, 
                    int numInstancesPerWindow, int groupID, double** ev_p,
                    PotentialShifter &potentialShifter, double *potentialShifts_p,
                    const char *init_file_name=NULL, const char *out_file_name=NULL,
                    const char *states_filename=NULL);

  double getDOSRatio(int instance, double energy);

  bool determineAcceptance(int instance, double energy);
  bool determineAcceptance(int instance, double energy, double magnetization);

  bool updateHistogram(int instance, double *evecs, bool accepted);
  //bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool accepted);
  bool updateHistogram(int instance, double *evecs, double *potentialShifts, bool accepted);

  bool updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, int check);
  bool updateHistogramFromRE(int instance, double *evecs, double energy, int check);
  bool updateHistogramFromRE(int instance, double *evecs, double energy, double *potentialShifts, int check);
  bool updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, double *potentialShifts, int check);
  
  void generateEvec(int instance, double *evecs, bool accepted);
  //void generateEvec(int instance, double *evecs, double energy);

  void generatePotentialShift(int instance, double *potentialShifts, bool accepted);

  void initializeEvecAndPotentialShift(int instance, double *evecs, double *potentialShifts);

  void initializeEvec(int instance, double *evecs);

  void generateUnsampledEvec(int instance, double *evecs, double energy)
  {
    initializeEvec(instance, evecs); 
    //return false;
  }

  void startSampling(void)
  { sw.writeHeader(gamma, n_walkers, n_spins, evecs_pointer); }

  void writeState(const char *name);
  void writeDos(const char *name);

 private:
  int n_walkers;
  int n_spins;
  int walkerID;
  double ** evecs_pointer;
  double *potentialShifts_pointer;
  int n_initialized_from_file;

  std::string dos_out_name;

  int stepsSinceLastHistogramUpdate;
  int numberOfUpdatesSinceLastBoost;
  int cycleCount;
  int modificationFactorChanges;

  // Random number generator and distribution:
  RNG rng;
  // boost::uniform_real<double> rnd;            //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd;    //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd01;  //(0.0,1.0);
  std::uniform_real_distribution<double> rnd11;  //(-1.0,1.0);
  std::uniform_real_distribution<double> rnd0pi; //(0.0,2.0*M_PI);

  /*
  // Histogramm and dos:
  double xMin, xMax, interval;
  int nX;
  double *dos; // std::vector<double> dos;
  int *histo; // std::vector<int> histo;
  int hMinimum;
  */

  double hMinimum;
  Graph1dMoments<double,double> dos, histo;
  Kernel1d<double,double> dosKernel, histoKernel, nullKernel;
  KernelType kernelType;

  unsigned long accept, reject, acceptSinceLastChange;
  double flatnessCriterion;

  double gamma, gammaFinal;
  int flipPerUpdate, updateCycle;

  // instance specific:
  std::vector<long> ref0, ref1;
  std::vector<double> position, magnetizationAtPosition;
  std::vector<bool> out;
  std::vector<int> lastChange;
  std::vector<int> lastChangePotentialShift;
  std::vector<int> lastAccepted;
  std::vector<int> lastAcceptedPotentialShiftIndex;
  std::vector<double> lastAcceptedEnergy;
  std::vector<double> lastAcceptedEvec;
  std::vector<double> oldSpin;  // oldSpin[instance*3 + {x=0, y=1, z=2}]
  std::vector<double> oldPotentialShift;
  std::vector<double> lastAcceptedPotentialShift;

  std::vector<int> numRetentions;

  char *statesFile;
  StatesWriter sw;

  int changeMode;
  bool histogramUpdateMode;
  int updatesPerBin;

  struct {double kappa, lambda, omega; int frequency, changes;} globalUpdate;

#ifdef ISING
    void inline random_evec_1(double ev[3])
  {
    ev[0]=ev[1]=0.0;
    ev[2]=1.0;
    if(rng()%2 == 0) ev[2]=-ev[2];
  }
#else
  void inline random_evec_1(double ev[3])
  {
/*
    double x,y,z;
    do {
      x = rnd(rng);
      y = rnd(rng);
    } while(x*x+y*y>1);
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r;
    y *= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z;
    r=1.0/sqrt(x*x+y*y+z*z);
    ev[0]=x*r; ev[1]=y*r; ev[2]=z*r;
*/

    double u = rnd11(rng);
    double v = rnd11(rng);

    while (u*u+v*v >= 1.0)
    {
      u = rnd11(rng);
      v = rnd11(rng);
    }

    ev[0]=2.0*u*std::sqrt(1.0-u*u-v*v);
    ev[1]=2.0*v*std::sqrt(1.0-u*u-v*v);
    ev[2]=1.0-2.0*(u*u+v*v);
  }
#endif

#ifdef ISING    
  void inline random_evec(double ev[3])
  {
    ev[2]=-ev[2];
  }
#else
  void inline random_evec(double ev[3])
  {
/*
    double x, y, z;
    do {
      x = rnd(rng); y = rnd(rng);
    } while(x*x+y*y>1); 
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r; y*= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z; 
    // Project out the parallel component;
    r = x*ev[0] + y*ev[1] + z*ev[2];
    x -= r*ev[0]; y -= r*ev[1]; z -= r*ev[2];
    r = x*x + y*y + z*z;
    double t = 1-0.3*rnd(rng);
    ev[0] *= t; ev[1] *= t; ev[2] *= t;
    r = sqrt((1-t*t)/r);
    ev[0] += x*r; ev[1] += y*r; ev[2] += z*r;
    r=1.0/sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
    ev[0]*=r; ev[1]*=r; ev[2]*=r;
*/
    
      
    double theta = rnd0pi(rng);
    double u = rnd11(rng);

    ev[0]=std::sqrt(1.0-u*u) * std::cos(theta);
    ev[1]=std::sqrt(1.0-u*u) * std::sin(theta);
    ev[2]=u;
  }
#endif

  double minVxShift {0.0}, maxVxShift {0.0}, rangeVxShift {0.0};

  double inline randomPotentialShift()
  {
    return minVxShift + rangeVxShift * rnd01(rng);
  }
};


//Constructor
template<class RNG>
WL1dEvecGenerator<RNG>::WL1dEvecGenerator(int num_spins, int num_instances, int numWindows,
                                          int numInstancesPerWindow, int groupID, double** ev_p,
                                          PotentialShifter &potentialShifter, double *potentialShifts_p,
					  const char *init_file_name, const char *out_file_name, 
                                          const char *states_filename)
: sw(states_filename), rnd01(0.0,1.0), rnd11(-1.0,1.0), rnd0pi(0.0,2.0*M_PI)
{

  if (potentialShifter.vSpinShiftFlag) {
    minVxShift = potentialShifter.minShift;
    maxVxShift = potentialShifter.maxShift;
    rangeVxShift = maxVxShift - minVxShift;
  }
  potentialShifts_pointer = potentialShifts_p;

  verbosity = 3;

  changeMode = 4;

  globalUpdate.frequency = 0;
  globalUpdate.changes = 0;
  globalUpdate.kappa = 1.0;
  globalUpdate.lambda = 1.0;
  globalUpdate.omega = 0.5;

  histogramUpdateMode = false;
  updatesPerBin = 100;

  long nX = -1;
  double interval = 0.01;
  double kernelWidth = 0.1;
  double xMin = -std::numeric_limits<double>::max();
  double xMax = 1.0;

  double xMinForEachWindow[numWindows];
  double xMaxForEachWindow[numWindows];
  double overlap = 0.75;
  double xMinFullRange = -std::numeric_limits<double>::max();
  double xMaxFullRange = 1.0;
  bool readLocalEnergyRange = false;
  bool readEnergyRangeForRestart = false;
  bool readPositions = false;

  // Ying Wai's Note  (Feb 5, 14):
  // * seed_seq is a sequence needed for generating random number seeds for different walkers
  // * initialized with system time here, but will be overridden if seed is provided in input file
  // * initialization from input file is recommended
  std::seed_seq seq {time(nullptr)};

  n_spins = num_spins;
  n_walkers = num_instances;
  walkerID = groupID;
  n_initialized_from_file = 0;
  evecs_pointer = ev_p;

  ref0.resize(n_walkers);
  for (int i=0; i<n_walkers; i++)
    ref0[i] = -1;         // ref0[i] = HUGE;
  ref1.resize(n_walkers);
  position.resize(n_walkers);
  magnetizationAtPosition.resize(n_walkers);
  out.resize(n_walkers);

  lastChange.resize(n_walkers);
  lastAccepted.resize(n_walkers);
  lastAcceptedEnergy.resize(n_walkers);
  lastAcceptedEvec.resize(3*n_walkers);

  oldSpin.resize(3*n_walkers);
  oldPotentialShift.resize(n_walkers);

  lastChangePotentialShift.resize(n_walkers);
  lastAcceptedPotentialShiftIndex.resize(n_walkers);
  lastAcceptedPotentialShift.resize(n_walkers);

  for (int i=0; i<n_walkers; i++)
  {
    lastAccepted[i] = -2;
    lastAcceptedPotentialShiftIndex[i] = -2;
    lastAcceptedPotentialShift[i] = 0.0;
    //YingWai's trial fix    (Sep 2, 14)
    initializeEvecAndPotentialShift(i, evecs_pointer[i], potentialShifts_pointer);
    //initializeEvec(i,evecs_pointer[i]);
  }

  for (int i=0; i<3*n_walkers; i++)
    lastAcceptedEvec[i] = 0.0;

  statesFile = NULL;
  if (states_filename != NULL)
  {
    statesFile = (char*) malloc(sizeof(char)*(1+strlen(states_filename)));
    strcpy(statesFile, states_filename);
  }

  numRetentions.resize(n_walkers);
  for(int i=0; i<n_walkers; i++) numRetentions[i] = 0;

  /*
  nX = -1;
  xMin = -HUGE; xMax= 1.0; interval = 0.01; // (xMax-xMin)/double(nX);
  */
  dos_out_name = "dos1d.out";
  stepsSinceLastHistogramUpdate = 0;
  numberOfUpdatesSinceLastBoost = 0;
  modificationFactorChanges = 0;
  cycleCount = 0;
  hMinimum = 1;     //   10
  acceptSinceLastChange = accept = reject = 0;
  gammaFinal = 1.e-6;
  flipPerUpdate = 1; //  100
  updateCycle = 5*num_instances;  // 1000
  gamma = 1.0;
  flatnessCriterion = 0.65;
  updatesPerBin = 100;

  kernelType = None;

  // special processing flags:
  int clearHistogram = 0;
  int setFirstWalkerToFM = 0;

  // dos = NULL;   // dos.resize(nX);
  // histo = NULL; // histo.resize(nX);
  dos.setDeltaAndClear(interval);
  histo.setDeltaAndClear(interval);

  if(init_file_name != NULL && init_file_name[0] != 0)
  {
//    if (walkerID == 0) {     // this is the processor with world_rank = 0

    std::string label, value;

    dos_out_name = std::string(init_file_name) + ".out";
    if(out_file_name != NULL && out_file_name[0] != 0)
      dos_out_name = out_file_name;
    std::ifstream inp(init_file_name);
    std::ostringstream buff;

    std::string line;
    while(std::getline(inp,line)) 
      buff << line << std::endl;

    inp.close();

    std::string fileString = buff.str();
    const char* fileChars  = fileString.c_str();
    json_t *json_root = NULL;

    json_parse_document(&json_root, (char *)fileChars);

    if(json_root == NULL || json_root->type != JSON_OBJECT)
    {
      std::ostringstream message;
      std::cerr << "WL1dEvecGenerator Constructor :: (" << init_file_name << ") parsing failed (bad format)\n";
      exit(1);
    }
  
    for(json_t *it = json_root->child; it != NULL; it=it->next)
    {
      std::string label = it->text;

      // Specifying energy ranges for a fresh start:
      // Choice 1: specify xMinFullRange, xMaxFullRange and overlap, subwindows will be calculated automatically
      // Choice 2: specify xMinForEachWindow and xMaxForEachWindow (arrays), this overrides Choice 1
      if(label=="xMinFullRange") xMinFullRange = atof(it->child->text);
      else if(label=="xMaxFullRange") xMaxFullRange = atof(it->child->text);
      else if(label=="overlap") overlap = atof(it->child->text);
      else if(label=="xMinForEachWindow")
      {
        readLocalEnergyRange = true;
        json_t *a = it->child;
        int j = 0;
        for (json_t *i = a->child; i != NULL; i = i->next)
          xMinForEachWindow[j++] = atof(i->text);
        if (j != numWindows) {
          std::cout << "ERROR #(xMin) " << j << " != numWindows "<< numWindows << std::endl;
          exit(1);
        }
      }
      else if(label=="xMaxForEachWindow")
      {
        readLocalEnergyRange = true;
        json_t *a = it->child;
        int j = 0;
        for (json_t *i = a->child; i != NULL; i = i->next)
          xMaxForEachWindow[j++] = atof(i->text);
        if (j != numWindows) {
          std::cout << "ERROR #(xMin) " << j << " != numWindows "<< numWindows << std::endl;
          exit(1);
        }
      }
      // Specifying energy ranges for a restart:
      // xMin and xMax specify the range for the subwindow
      else if(label=="xMin") 
      {
        readEnergyRangeForRestart = true;
        xMin = atof(it->child->text);
      }
      else if(label=="xMax")
      {
        readEnergyRangeForRestart = true;
        xMax = atof(it->child->text);
      }

      else if(label=="interval") interval = atof(it->child->text);
      else if(label=="kernelWidth") kernelWidth = atof(it->child->text);
      else if(label=="kernelType")
      {
        std::string strValue(it->child->text);
        kernelType = getKernelType(strValue);
      }
      else if(label=="gamma") gamma = atof(it->child->text);
      else if(label=="gammaFinal") gammaFinal = atof(it->child->text);
      else if(label=="nX")
      {
        nX = atoi(it->child->text);
        if (readEnergyRangeForRestart) {
	  dos.setRangeAndClear(xMin,xMax,nX);
	  histo.setRangeAndClear(xMin,xMax,nX);
        }
	/*
        if(dos!=NULL) free(dos);
        if(histo!=NULL) free(histo);
        dos=(double *)calloc(nX,sizeof(double));
        histo=(int *)calloc(nX,sizeof(int));
	*/
      }
      else if(label=="flipPerUpdate") flipPerUpdate=atoi(it->child->text);
      else if(label=="updateCycle") updateCycle=atoi(it->child->text);
      else if(label=="cycleCount") cycleCount=atoi(it->child->text);
      else if(label=="changeMode") changeMode=atoi(it->child->text);
      else if(label=="flatnessCriterion") flatnessCriterion=atof(it->child->text);
      else if(label=="histogramMinimum") hMinimum=atof(it->child->text);
      else if(label=="updatesPerBin") updatesPerBin=atoi(it->child->text);
      else if(label=="globalUpdate.frequency") globalUpdate.frequency=atoi(it->child->text);
      else if(label=="globalUpdate.changes") globalUpdate.changes=atoi(it->child->text);
      else if(label=="globalUpdate.kappa") globalUpdate.kappa=atof(it->child->text);
      else if(label=="globalUpdate.lambda") globalUpdate.lambda=atof(it->child->text);
      else if(label=="globalUpdate.omega") globalUpdate.omega=atof(it->child->text);
      else if(label=="seed") {
        //rng.seed(atoi(it->child->text));
        seq = {atoi(it->child->text)};
      }
      else if(label=="accept") accept=atol(it->child->text);
      else if(label=="acceptSinceLastChange") acceptSinceLastChange=atol(it->child->text);
      else if(label=="reject") reject=atol(it->child->text);
//*
      else if(label=="rngState")
      {
        std::string strValue(it->child->text);
        std::stringstream strStream(strValue, std::stringstream::in);
        strStream>>rng;
      }
//*/
      else if(label=="dos")
      {
        json_t *a = it->child;
        int j = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
          dos[j++] = atof(i->text);
        if(j != dos.getN())
        {
          std::cout << "ERROR #(dos) "<< j << " != nX " << dos.getN() << std::endl;
          exit(1);
        }
      }
      else if(label=="histo")
      {
        json_t *a = it->child;
        int j = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
          histo[j++] = atof(i->text);
        if(j != histo.getN())
        {
          std::cout << "ERROR #(histo) "<< j << " != nX " << histo.getN() << std::endl;
          exit(1);
        }
      }
      else if(label=="moments")
      {
        json_t *a = it->child->child;
        int k = atoi(a->text);
        dos.setNumberOfMoments(k);
        a = a->next;
        int jj = 0;
        for(json_t *j=a->child; j!=NULL; j=j->next)
        {
          dos.setNumberOfSamplesAtIdx(jj, atoi(j->text));
          jj++;
        }
        json_t *i=a->next;
        for(int kk=0; kk<k; kk++)
        {
          int jj = 0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
          {
            dos.setMomentAtIdx(jj, kk, atof(j->text));
            jj++;
          }
          i = i->next;
        }
      }
      else if(label=="moments.k")
      {
        dos.setNumberOfMoments(atoi(it->child->text));
      }
      else if(label=="evecs")
      {
        json_t *a = it->child;
        n_initialized_from_file = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized_from_file < n_walkers)
          {
            int k = 0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
              evecs_pointer[n_initialized_from_file][k++] = atof(j->text);
            // initialize oldSpin and lastChange to point to site 0
            lastChange[n_initialized_from_file] = 0;
            oldSpin[  3*n_initialized_from_file] = evecs_pointer[n_initialized_from_file][0];
            oldSpin[1+3*n_initialized_from_file] = evecs_pointer[n_initialized_from_file][1];
            oldSpin[2+3*n_initialized_from_file] = evecs_pointer[n_initialized_from_file][2];
            n_initialized_from_file++;
            if(k != 3*n_spins) {
              std::cout << "ERROR #(evecs) " << k << " != 3*n_spins " << 3*n_spins << std::endl;
              exit(1);
            }
          }
        }
      }
      else if(label=="potentialShifts")
      {
        json_t *a = it->child;
        int j = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
          potentialShifts_pointer[j++] = atof(i->text);
        if (j != n_spins) {
          std::cout << "ERROR #(potentialShifts) " << j << " != n_spins " << n_spins << std::endl;
          exit(1);
        }
      }
      else if(label=="oldSpin")
      {
        json_t *a = it->child;
        int n_initialized = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized < n_walkers)
          {
            int k = 0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
              oldSpin[n_initialized*3+k++] = atof(j->text);
            n_initialized++;
            if(k != 3) {
              std::cout << "ERROR #(oldSpin) " << k << " != 3" << std::endl;
              exit(1);
            }
          }
        }
      }
      else if(label=="lastChange")
      {
        json_t *a = it->child;
        int j = 0;
        int n_initialized = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized < n_walkers)
          {
            lastChange[j++] = atoi(i->text);
            n_initialized++;
          }
        }
      }
      else if(label=="lastChangePotentialShift")
      {
        json_t *a = it->child;
        int j = 0;
        int n_initialized = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized < n_walkers)
          {
            lastChangePotentialShift[j++] = atoi(i->text);
            n_initialized++;
          }
        }
      }
      else if(label=="position")
      {
        json_t *a = it->child;
        int j = 0;
        int n_initialized = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized < n_walkers)
          {
            position[j++] = atof(i->text);
            n_initialized++;
          }
        }
        readPositions = true;
      }
      else if(label=="magnetizationAtPosition")
      {
        json_t *a = it->child;
        int j = 0;
        int n_initialized = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized < n_walkers)
          {
            magnetizationAtPosition[j++] = atof(i->text);
            n_initialized++;
          }
        }
//        readMagnetizationAtPositions=true;
      }
      else if(label=="oldPotentialShift")
      {
        json_t *a = it->child;
        int j = 0;
        int n_initialized = 0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized < n_walkers)
          {
            oldPotentialShift[j++] = atof(i->text);
            n_initialized++;
          }
        }
      }
      else if(label=="stepsSinceLastHistogramUpdate") stepsSinceLastHistogramUpdate=atoi(it->child->text);
      else if(label=="numberOfUpdatesSinceLastBoost") numberOfUpdatesSinceLastBoost=atoi(it->child->text);
      else if(label=="modificationFactorChanges") modificationFactorChanges=atoi(it->child->text);
      else if(label=="clearHistogram") clearHistogram=atoi(it->child->text);
      else if(label=="setFirstWalkerToFM") setFirstWalkerToFM=atoi(it->child->text);
      else std::cout<<"WARNING: unknown label: "<<label<<std::endl;
    }

    json_free_value(&json_root);

    // YingWai: need to communicate all the read-in data to different group leaders    (Dec 17, 13)
    // Just for testing to see if communication is OK:
    //int yingwai = 817;
    //MPI_Bcast(&yingwai, 1, MPI_INT, 0, MPI_COMM_WORLD);

//  }
//  else {
    // Other group leaders receive from world_rank = 0
    // Just for testing to see if communication is OK:
    //MPI_Bcast(&yingwai, 1, MPI_INT, 0, MPI_COMM_WORLD);
//  }

  }// finished reading in from file

  // Seeding random number generators
  std::vector<unsigned> rngSeed(numWindows * numInstancesPerWindow);
  //std::vector<unsigned> rngSeed(n_walkers + 2);
  seq.generate(rngSeed.begin(), rngSeed.end());
  rng.seed(rngSeed[walkerID]);

  // Determining energy ranges:
  if (!readEnergyRangeForRestart) {
    // Choice 1 follow-up: Each instance calculates its own xMin and xMax
    if (!readLocalEnergyRange) {
      double energySubwindowWidth = (xMaxFullRange - xMinFullRange) / (1.0 + double(numWindows - 1) * (1.0 - overlap));
      xMin = xMinFullRange + floor(groupID / numInstancesPerWindow) * (1.0 - overlap) * energySubwindowWidth;
      xMax = xMin + energySubwindowWidth;
    }
    // Choice 2 follow-up: Each instance takes the corresponding read-in values 
    else {
      xMin = xMinForEachWindow[groupID];
      xMax = xMaxForEachWindow[groupID];
    }
  }
  if (nX < 0)
    nX = (xMax - xMin) / interval;
  else if (nX > 0)
    interval = (xMax - xMin) / double(nX);

  dos.setRangeAndClear(xMin, xMax, nX);
  histo.setRangeAndClear(xMin, xMax, nX); 

/*

  // Check if energy ranges, dos and energy are read in correctly
  printf("WL1dEvecGenerator Constructor :: setting energy range.
          GroupID = %5d: [%10.3f, %10.3f]\n", groupID, xMin, xMax);
  printf("WL1dEvecGenerator Constructor :: GroupID: %5d:
          # of elements in DOS = %10d\n", groupID, dos.getN());
  printf("WL1dEvecGenerator Constructor :: GroupID: %5d: 
          initial energy (position[0]) = %20.8f\n", position[0]);
*/

  // set xMax and xMin or interval depending on nX:
  //if(nX > 0)
  //{
  //  interval = (xMax-xMin) / double(nX);
  //  dos.setRange(xMin,xMax);
  //  histo.setRange(xMin,xMax);
  //}
  //else
  //{
  //  dos.setDeltaAndClear(interval);
  //  histo.setDeltaAndClear(interval);
  //}

  dosKernel.setWidthAndClear(interval,kernelWidth);
  histoKernel.setWidthAndClear(interval,kernelWidth);
  nullKernel.setWidthAndClear(interval,kernelWidth);
  initKernel(kernelType,dosKernel);
  dosKernel.scale(gamma/dosKernel(0.0));
  initKernel(kernelType,histoKernel);
  histoKernel.scale(1.0/histoKernel(0.0));
  initKernel(kernelType,nullKernel);
  nullKernel.scale(0.0);
//  histoKernel.scale(interval);

  if(readPositions)
    for(int i=0; i<n_walkers; i++)
      ref0[i] = dos.idx_withinRange(position[i]);

  if(clearHistogram != 0)
  {
    std::cout << "Clearing Wang-Landau Histogram.\n";
    histo.clear();
  }

  if(setFirstWalkerToFM != 0)
  {
    std::cout << "Setting first walker to FM.\n";
    for(int i=0; i<n_spins; i++)
    {
      evecs_pointer[0][  3*i] = 0.0;
      evecs_pointer[0][1+3*i] = 0.0;
      evecs_pointer[0][2+3*i] = 1.0;
    }
// initialize oldSpin and lastChange to point to site 0
    lastChange[0] = 0;
    oldSpin[0] = evecs_pointer[0][0];
    oldSpin[1] = evecs_pointer[0][1];
    oldSpin[2] = evecs_pointer[0][2];

// initialize oldPotentialShift and lastChangePotentialShift to point to site 0
    lastChangePotentialShift[0] = 0;
    oldPotentialShift[0] = 0.0;
  }

  if(out_file_name != NULL && out_file_name[0] != 0)
    dos_out_name = out_file_name;
  std::cout << "Wang-Landau output will be written to: " << dos_out_name << std::endl;

}


template<class RNG>
void WL1dEvecGenerator<RNG>::initializeEvecAndPotentialShift(int inst, double *evecs, double *potentialShifts)
{
  for (size_t j=0; j<n_spins; j++)
    potentialShifts[j] = 0.0;

  initializeEvec(inst, evecs);
}

// YingWai: antiferromagnetic case should be set for walkers in the last window, not inst = 1
// later ...     (Apr 30, 14)
template<class RNG>
void WL1dEvecGenerator<RNG>::initializeEvec(int inst, double *evecs)
{
  bool firstFerromagnetic = true;
  bool secondAntiferromagnetic = false;
  //if (inst >= n_initialized_from_file)
  //{
    if(firstFerromagnetic)
    {
      if(inst == 0)
        for(size_t j=0; j<3*n_spins; j+=3)
	{
          evecs[j] = 0.0;
          evecs[j+1] = 0.0;
          evecs[j+2] = 1.0;
        }
      else if (inst == 1)
        for(size_t j=0; j<3*n_spins; j+=3)
        {
          if(secondAntiferromagnetic)
            evecs[j+2] = (j/3)%2 ? 1.0 : -1.0 ;
          else
            evecs[j+2] = 1.0;
          evecs[j+1] = 0.0;
          evecs[j] = 0.0;
        }
      else
	for(size_t j=0; j<3*n_spins; j+=3)
	  random_evec_1(&evecs[j]);
    }
    else
      for(size_t j=0; j<3*n_spins; j+=3)
        random_evec_1(&evecs[j]);
  //}

  lastChange[0] = 0;
  oldSpin[0] = evecs[0];
  oldSpin[1] = evecs[1];
  oldSpin[2] = evecs[2];
  
  out[inst] = true;
}


template<class RNG>
void WL1dEvecGenerator<RNG>::writeState(const char* name)
{
  // if(!syncronizeGraphs(dos,histo))
  if(dos.getN()!=histo.getN())
  {
    std::cout<<"Histogramm size dosn't match DOS! Clearing histogramm!\n";
    histo.setRangeAndClear(dos.getMinX(),dos.getMaxX(),dos.getN());
  }
  std::ofstream ofile(name);
  if(ofile)
  {
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(8);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << ",\n";
    ofile<<"\"xMax\" : " << dos.getMaxX() << ",\n";
    ofile<<"\"nX\" : " << dos.getN() << ",\n";
    if(kernelType!=None)
      ofile<<"\"kernelWidth\" : " << dosKernel.getWidth() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma << ",\n";
    ofile<<"\"accept\" : " << accept <<",\n";
    ofile<<"\"acceptSinceLastChange\" : " << acceptSinceLastChange <<",\n";
    ofile<<"\"reject\" : " << reject <<",\n";
    ofile<<"\"gammaFinal\" : " << gammaFinal << ",\n";

    ofile<<"\"changeMode\" : "<< changeMode <<",\n";
    ofile<<"\"flatnessCriterion\" : "<< flatnessCriterion <<",\n";
    ofile<<"\"histogramMinimum\" : "<< hMinimum <<",\n";
    ofile<<"\"updatesPerBin\" : "<< updatesPerBin <<",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << ",\n";
    ofile<<"\"updateCycle\" : " << updateCycle << ",\n";

    ofile<<"\"globalUpdate.frequency\" : "<<  globalUpdate.frequency << ",\n";
    ofile<<"\"globalUpdate.changes\" : "<<  globalUpdate.changes << ",\n";
    ofile<<"\"globalUpdate.kappa\" : "<<  globalUpdate.kappa << ",\n";
    ofile<<"\"globalUpdate.lambda\" : "<<  globalUpdate.lambda << ",\n";
    ofile<<"\"globalUpdate.omega\" : "<<  globalUpdate.omega << ",\n";

    ofile<<"\"dos\" : ["<<std::endl;
    for(int i=0; i<dos.getN(); i++) ofile<<dos[i]<<((i==dos.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    ofile<<"\"histo\" : ["<<std::endl;
    for(int i=0; i<histo.getN(); i++) ofile<<histo[i]<<((i==histo.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    if(dos.getNumberOfMoments()>0)
    {
      ofile<<"\"moments\" : ["<<std::endl;
      ofile<<dos.getNumberOfMoments()<<", ["<<std::endl;
      for(int i=0; i<dos.getN(); i++) ofile<<dos.getNumberOfSamplesAtIdx(i)<<((i==dos.getN()-1)?"\n":",\n");
      ofile<<"], [\n";
      for(int j=0; j<dos.getNumberOfMoments(); j++)
      {
        for(int i=0; i<dos.getN(); i++) ofile<<dos.getMomentAtIdx(i,j)<<((i==dos.getN()-1)?"\n":",\n");
        if(j!=dos.getNumberOfMoments()-1) ofile<<"], [\n";
      }
      ofile<<"] ],\n";
    }
    ofile<<"\"rngState\" : \""<<rng<<"\",\n";
    ofile<<"\"evecs\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[\n";
      for(int j=0; j<3*n_spins; j+=3)
        ofile<<evecs_pointer[i][j]<<", "<<evecs_pointer[i][j+1]
             <<", "<<evecs_pointer[i][j+2]<<((j==3*n_spins-3)?"\n":",\n");
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    // YingWai: should be written out only when potential shift != 0   (Sep 2, 14) 
    ofile<<"\"potentialShifts\" : [\n";
    for(int j=0; j<n_spins; j++)
      ofile<<potentialShifts_pointer[j]<<((j==n_spins-1)?"\n":",\n");
    ofile<<"],\n";

    ofile<<"\"oldSpin\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[ "<<oldSpin[3*i]<<", "<<oldSpin[3*i+1]
             <<", "<<oldSpin[3*i+2];
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";

    // YingWai: should be written out only when potential shift != 0   (Sep 2, 14) 
    ofile<<"\"oldPotentialShift\" : [\n";
    for(int i=0; i<n_walkers; i++)
      ofile<< oldPotentialShift[i]<<((i==n_walkers-1)?"\n":",\n");
    ofile<<"],\n";

    ofile<<"\"lastChange\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< lastChange[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";

    // YingWai: should be written out only when potential shift != 0   (Sep 2, 14) 
    ofile<<"\"lastChangePotentialShift\" : [\n";
    for(int i=0; i<n_walkers; i++)
      ofile<< lastChangePotentialShift[i]<<((i==n_walkers-1)?"\n":",\n");
    ofile<<"],\n";

    ofile<<"\"position\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< position[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"magnetizationAtPosition\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< magnetizationAtPosition[i]<<((i==n_walkers-1)?"\n":",\n");
    }
    ofile<<"],\n";
    ofile<<"\"stepsSinceLastHistogramUpdate\" : " << stepsSinceLastHistogramUpdate << ",\n";
    ofile<<"\"numberOfUpdatesSinceLastBoost\" : " << numberOfUpdatesSinceLastBoost << ",\n";
    ofile<<"\"modificationFactorChanges\" : " << modificationFactorChanges << ",\n";
    ofile<<"\"cycleCount\" : " << cycleCount << "\n";
    ofile<<"}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
} 

template<class RNG>
void WL1dEvecGenerator<RNG>::writeDos(const char* name)
{
  std::ofstream ofile(name);
  // write dos;
  // we are using JSON as our file format
  if(ofile)
  {
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << "," <<std::endl;
    ofile<<"\"xMax\" : " << dos.getMaxX() << "," <<std::endl;
    ofile<<"\"nX\" : " << dos.getN() << "," <<std::endl;
    ofile<<"\"kernelWidth\" : " << dosKernel.getWidth() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma <<"," <<std::endl;
    ofile<<"\"gammaFinal\" : " << gammaFinal << "," <<std::endl;
    ofile<<"\"globalUpdate.changes\" : "<<  globalUpdate.changes << ",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << "," <<std::endl;
    ofile<<"\"updateCycle\" : " << updateCycle << "," <<std::endl;
    ofile<<"\"dos\" : ["<<std::endl;
    for(int i=0; i<dos.getN(); i++) ofile<<dos[i]<<((i==dos.getN()-1)?"\n":",\n");
    ofile<<"],\n";
    if(dos.getNumberOfMoments()>0)
    {
      ofile<<"\"moments\" : ["<<std::endl;
      ofile<<dos.getNumberOfMoments()<<", ["<<std::endl;
      for(int i=0; i<dos.getN(); i++) ofile<<dos.getNumberOfSamplesAtIdx(i)<<((i==dos.getN()-1)?"\n":",\n");
      ofile<<"], [\n";
      for(int j=0; j<dos.getNumberOfMoments(); j++)
      {
        for(int i=0; i<dos.getN(); i++) ofile<<dos.getMomentAtIdx(i,j)<<((i==dos.getN()-1)?"\n":",\n");
        if(j!=dos.getNumberOfMoments()-1) ofile<<"], [\n";
      }
      ofile<<"] ],\n";
    }
    ofile<<"\"histo\" : ["<<std::endl;
    for(int i=0; i<histo.getN(); i++) ofile<<histo[i]<<((i==histo.getN()-1)?"\n":",\n");
    ofile<<"]\n}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
}


// YingWai: Should deal with energy in the same bin the same way as in determineAcceptance
// ... later     (Feb 13, 14)
template<class RNG>
double WL1dEvecGenerator<RNG>::getDOSRatio(int instance, double energy)
{
  double DOSRatio = 0.0;

// YingWai: addition for debugging ref0     (Dec 4, 14)
  if (ref0[instance] < -1 || ref0[instance] > dos.getN()-1) {
    printf("ref0 error in getDOSRatio. ref0[%d] = %ld\n", instance, ref0[instance]);
    //exit(1);
  }

  // energy between xMin and xMax ?
  ref1[instance] = dos.idx_withinRange(energy);

  // YingWai: should this be updated?    (Feb 13, 14) 
  stepsSinceLastHistogramUpdate++; // counter[0]++

  if (ref1[instance] < 0 || ref1[instance] > dos.getN()-1) {
    out[instance] = true;
    DOSRatio = 0.0;
  }
  else {
    out[instance] = false;
    if (ref0[instance] < 0)
      DOSRatio = 1.0;
    else
      DOSRatio = exp(dos[ref0[instance]] - dos[ref1[instance]]);
  }

  return DOSRatio;
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::determineAcceptance(int instance, double energy)
{
  return determineAcceptance(instance, energy, 0.0);
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::determineAcceptance(int instance, double energy, double magnetization)
{

  bool accept_step (false);
  // +++++++++ Added by Odbadrakh and Don on Aug 11, 2010
  // +++ If a flip of spin results in energy in the same bin as previous accepted 
  // +++ state, it is accepted all the time. We change it by accepting if the
  // +++ new state has lower density of states according to a linear spline of DOS
  
  double dos_differ {0.0};
  double energy_differ {0.0};
  double to_go_or_not {0.0};

  // +++++++++ end of modification. Follow the new variables to see the actual changes 

  // YingWai: addition for debugging ref0     (Dec 4, 14)
  if (ref0[instance] < -1 || ref0[instance] > dos.getN()-1) {
    printf("ref0 error in determineAcceptance. ref0[%d] = %ld\n", instance, ref0[instance]);
    //exit(1);
  }

  // energy between xMin and xMax ?
  ref1[instance] = dos.idx_withinRange(energy);

  // dos.test();

  stepsSinceLastHistogramUpdate++; // counter[0]++

  // record initial energy for step out
  //if (lastAccepted[instance] == -2)
  //{
  //  position[instance] = lastAcceptedEnergy[instance] = energy;
  //  lastAccepted[instance] = -1;
  //}
 
  numRetentions[instance]++;
  //ref1[instance] = grid;

  if (ref1[instance] < 0 || ref1[instance] > dos.getN()-1)  // energy out of range, reject the move
  {
    out[instance] = true;
    accept_step = false;
    // YingWai: what to set for numRetentions[instance] ?
    // probably no changes are needed...   (Jan 27, 14)
  }
  else {

    out[instance] = false;
    if (ref0[instance] < 0)    // no energy has been visited before
    {
      accept_step = true;
      //ref0[instance] = ref1[instance];
      //position[instance] = energy;
      //magnetizationAtPosition[instance] = magnetization;
    }
    else 
    {
      if (abs(ref1[instance] - ref0[instance]) < 1)  // same energy bin
      {
// Actual change made by Odbadrakh, Aug 30, 2010
        dos_differ = dos[ref0[instance]-1] - dos[ref0[instance]];
        energy_differ = energy - lastAcceptedEnergy[instance];
        to_go_or_not = dos_differ * energy_differ;
        if (to_go_or_not >= 0.0)        // accepts all downhill changes
        {
          accept_step = true;
          //ref0[instance] = ref1[instance];
          //position[instance] = energy;
          //magnetizationAtPosition[instance] = magnetization;
        } 
        else                            // uphill moves
        {
          if(rnd(rng) < exp(to_go_or_not/dos.getDelta()))
          { //std::cout<<" dos "<<instance<<"  weight "<<dos.getDelta()<<std::endl;
            accept_step = true;
            //ref0[instance] = ref1[instance];
            //position[instance] = energy;
            //magnetizationAtPosition[instance] = magnetization;
          }
          else
          {
            accept_step = false;
          }
        }
      }
      else
      { 
        if (dos[ref1[instance]] <= dos[ref0[instance]])
        {
          accept_step = true;
          //ref0[instance] = ref1[instance];
          //position[instance] = energy;
          //magnetizationAtPosition[instance] = magnetization;
        }
        else
        {
          if(rnd(rng) < exp(dos[ref0[instance]] - dos[ref1[instance]]))
          {
            accept_step = true;
          //ref0[instance] = ref1[instance];
          //position[instance] = energy;
          //magnetizationAtPosition[instance] = magnetization;
          }
          else 
          {
            accept_step = false;
          }
        }
      
      }

    }
  }

  if (accept_step) {
    ref0[instance] = ref1[instance];
    position[instance] = energy;
    magnetizationAtPosition[instance] = magnetization;
  }

// End of change made on Aug 30, 2010
  if (verbosity > 2)
    std::cout.precision(10);
    std::cout << "WangLandau 1d EvecGenerator step "
              << modificationFactorChanges << " : " << numberOfUpdatesSinceLastBoost << " : "
              << stepsSinceLastHistogramUpdate << " nX = " << dos.getN()
	      << " [" << dos.getMinX() << ", " << dos.getMaxX() << "] "
              << (accept_step ? "accepted" : "rejected")
              << " Energy = " << energy << ", Instance " << instance << std::endl;

  return accept_step;

}


/*
template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, double energy, bool *accepted)
{
  return updateHistogram(instance, evecs, energy, 0.0, accepted);
}
*/

template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, double *potentialShifts, bool accepted)
{

  if (accepted) {
    lastAcceptedPotentialShiftIndex[instance] = lastChangePotentialShift[instance];
    lastAcceptedPotentialShift[instance] = potentialShifts[lastChangePotentialShift[instance]];
  }

  return updateHistogram(instance, evecs, accepted);
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, bool accepted)
{

  if (lastAccepted[instance] == -2 && out[instance])      // initial state outside energy range
    return false;

// Update histogram and DOS
  if(stepsSinceLastHistogramUpdate >= flipPerUpdate)
  {
    stepsSinceLastHistogramUpdate = 0;       // counter[0]
    numberOfUpdatesSinceLastBoost++;         // counter[1]
    cycleCount++;

    // modify the DoS
    //if(!out[instance])
    //{

      // addKernel(dos,dosKernel,energy);
      if (!histogramUpdateMode)             // standard Wang-Landau
      {
        addKernel(dos, dosKernel, position[instance]);
        dos.addMomentsAtIdx(dos.idx_withinRange(position[instance]), magnetizationAtPosition[instance]);
        // if(accepted) addKernel(histo,histoKernel,position[instance]);
        addKernel(histo, histoKernel, position[instance]);
        // addKernel(dos, dosKernel, energy);
        // addKernel(histo, histoKernel, energy);
      }
      else
      {
        // addKernel(dos, nullKernel, position[instance]);
        if (accepted) addKernel(histo, histoKernel, position[instance]);
      }
    //} 
    //else
    //{
    //  std::cerr << "ATTENTION: We should never reach this place in WL1dEvecGenerator!\n";
    //  exit(1);
    //}
  }

// 1/t algorithm, change gamma at every step
// Reference: Phys. Rev. E 75, 046701 (2007).
  if(changeMode & (8+16+32))
  {
    long n = accept + reject;
    double dn;

    if ((changeMode &  (8+16)) == 8) n = accept;
    else if ((changeMode &  (8+16)) == 16) n = reject;
    else if ((changeMode &  (8+16)) == 8+16) n = accept + reject;

    if (changeMode & 32) dn = double(n) / double(n_walkers);
    else dn = double(n);

    gamma = 1.0 / dn;
    if (gamma > 1.0) gamma = 1.0;

    initKernel(kernelType, dosKernel);
    dosKernel.scale(gamma);
  }

// 1. write configuration
// 2. check histogram flatness
// 3. perform global update
  if (cycleCount >= updateCycle)
  {
    cycleCount = 0;
    // syncronizeGraphs(dos, histo);
    if (dos.getN() != histo.getN())
    {
      std::cout << "Histogramm size dosn't match DOS! Clearing histogram!\n";
      histo.setRangeAndClear(dos.getMinX(), dos.getMaxX(), dos.getN());
    }

    char stateFile[51];
    sprintf(stateFile, "WL1d.state%05d", walkerID);
    writeState(stateFile);

    if (!histogramUpdateMode)        // Normal Wang-Landau
    {
      // calculate minimum nonzero histogram
      // we look only at the histogram inside the energy interval that was actually sampled if we use kernel updates
      double hMin, hMax, hMean;

      if (kernelType == None)
      {
        hMean = histo.getMeanY();
        histo.getMinMaxY(hMin,hMax);
      }
      else 
      {
        hMean = histo.getMeanYInInterval(dos.getMinX() + dosKernel.getWidth(),
                                         dos.getMaxX() - dosKernel.getWidth());
        histo.getMinMaxYInInterval(dos.getMinX() + dosKernel.getWidth(),
                                   dos.getMaxX() - dosKernel.getWidth(),
                                   hMin, hMax);
      }

      // double currentFlatnessCriterion = double(hMin-hMax) / hMean;
      double currentFlatnessCriterion = double(hMin) / hMean;
      
      std::cout << "# acceptance ratio = " << double(accept) / double(accept+reject) << "\n";
      std::cout << "# current flatness = " << currentFlatnessCriterion << (changeMode & 4 ? " *":"") << "\n";
      std::cout << "# current histogram minimum = " << hMin << (changeMode & 2 ? " *":"") << "\n";
      std::cout << "# average accepted steps/bin since last gamma change = "
                << double(acceptSinceLastChange) / double(histo.getN())
                << (changeMode & 1 ? " *":"") << "\n";

      if (changeMode != 0 && changeMode < 8)
      {
        // perform global update
        if (globalUpdate.frequency > 0 && (globalUpdate.frequency*histo.getN()) < acceptSinceLastChange) 
        {
          char fn[256];
          snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
          writeDos(fn);
          globalUpdate.changes++;
          std::cout<<"# global update "<<globalUpdate.changes<<std::endl;
          performGlobalUpdate(dos,globalUpdate.kappa,globalUpdate.lambda,globalUpdate.omega);
          histo.clear();
          acceptSinceLastChange=0;
        }
        else if (((acceptSinceLastChange>=histo.getN()*updatesPerBin || !(changeMode &1)) &&
                  (hMin >= hMinimum || !(changeMode & 2)) &&
                  (currentFlatnessCriterion>flatnessCriterion || !(changeMode & 4)) ))
        {
          std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " is finished.\n";
          modificationFactorChanges++; // counter[2]

          // write dos;
          // we are using JSON as our file format
          char fn[256];
          snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
          writeDos(fn);
          // writeDos(dos_out_name.data());

          // clear the Histogram
          histo.clear();
          acceptSinceLastChange = 0;
 
          // change gamma
          dosKernel.scale(0.5);
          gamma = 0.5 * gamma;
          std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " begins.\n";
          if(statesFile != NULL)
	  {
            // char fn[256];
            snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
            sw.newFile(fn);
            sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
          }
        }
      }
    } 
    else            // histogram update mode != 0
    {
      if (acceptSinceLastChange >= histo.getN() * updatesPerBin)
      {
        acceptSinceLastChange = 0;
        for (int ix=0; ix<dos.getN(); ix++)
          dos[ix] += gamma * histo[ix];
        histo.clear();
        gamma = 0.5 * gamma;

        std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " begins.\n";
        if (statesFile != NULL)
        {
          char fn[256];
          snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
          sw.newFile(fn);
          sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
        }
      }
    }
  }

  if (accepted) 
  {
    sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);

    lastAccepted[instance] = lastChange[instance];
    lastAcceptedEnergy[instance] = position[instance];
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]];

    accept++;
    acceptSinceLastChange++;
  }
  else reject++;

  if (gamma < gammaFinal) return true;
  else return false;

}


template<class RNG>
void WL1dEvecGenerator<RNG>::generatePotentialShift(int instance, double *potentialShifts, bool accepted)
{

  if (accepted) {
    lastChangePotentialShift[instance] = int(rnd(rng)*n_spins);
    oldPotentialShift[instance] = potentialShifts[lastChangePotentialShift[instance]];
    potentialShifts[lastChangePotentialShift[instance]] = randomPotentialShift();
  }
  else {
    potentialShifts[lastChangePotentialShift[instance]] = oldPotentialShift[instance];
    lastChangePotentialShift[instance] = int(rnd(rng)*n_spins);
    oldPotentialShift[instance] = potentialShifts[lastChangePotentialShift[instance]];
    potentialShifts[lastChangePotentialShift[instance]] = randomPotentialShift();
  }

}


template<class RNG>
void WL1dEvecGenerator<RNG>::generateEvec(int instance, double *evecs, bool accepted)
{
  if (accepted)
  {
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
    numRetentions[instance] = 0;
  }
  else 
  {
    evecs[  3*lastChange[instance]] = oldSpin[  3*instance];
    evecs[1+3*lastChange[instance]] = oldSpin[1+3*instance];
    evecs[2+3*lastChange[instance]] = oldSpin[2+3*instance];
    lastChange[instance] = int(rnd(rng)*n_spins);
    oldSpin[  3*instance] = evecs[  3*lastChange[instance]];
    oldSpin[1+3*instance] = evecs[1+3*lastChange[instance]];
    oldSpin[2+3*instance] = evecs[2+3*lastChange[instance]];
    random_evec(&evecs[3*lastChange[instance]]);
  }

}


// The following four are for the replica exchange to call to update histogram
template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogramFromRE(int instance, double *evecs, double energy, double *potentialShifts, int check)
{
  return updateHistogramFromRE(instance, evecs, energy, 0.0, potentialShifts, check);
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, double *potentialShifts, int check)
{

  lastAcceptedPotentialShiftIndex[instance] = lastChangePotentialShift[instance];
  lastAcceptedPotentialShift[instance] = potentialShifts[lastChangePotentialShift[instance]];
  return updateHistogramFromRE(instance, evecs, energy, magnetization, check);
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogramFromRE(int instance, double *evecs, double energy, int check)
{
  return updateHistogramFromRE(instance, evecs, energy, 0.0, check);
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::updateHistogramFromRE(int instance, double *evecs, double energy, double magnetization, int check)
{

//  if (lastAccepted[instance] == -2 && out[instance])      // initial state outside energy range
//    return false;

// Update histogram and DOS
  //if(stepsSinceLastHistogramUpdate >= flipPerUpdate)
  //{
    // YingWai: are these correct to be kept here?    (Feb 15, 14)
    stepsSinceLastHistogramUpdate = 0;       // counter[0]
    numberOfUpdatesSinceLastBoost++;         // counter[1]
    cycleCount++;

    ref0[instance] = ref1[instance];
    position[instance] = energy;
    magnetizationAtPosition[instance] = magnetization;

    //YingWai's addition for ref0, ref1 debugging   (Dec 5, 14)
    printf("ref1 check in updateHistogramFromRE. ref1[%d] = %ld, ref0 = %ld\n", instance, ref1[instance], ref0[instance]);

    // modify the DoS
      // addKernel(dos,dosKernel,energy);
      if (!histogramUpdateMode)             // standard Wang-Landau
      {
        addKernel(dos, dosKernel, position[instance]);
        dos.addMomentsAtIdx(dos.idx_withinRange(position[instance]), magnetizationAtPosition[instance]);
        // if(accepted) addKernel(histo,histoKernel,position[instance]);
        addKernel(histo, histoKernel, position[instance]);
        // addKernel(dos, dosKernel, energy);
        // addKernel(histo, histoKernel, energy);
      }
      else
      {
        // addKernel(dos, nullKernel, position[instance]);
        addKernel(histo, histoKernel, position[instance]);
      }
  //}

// 1/t algorithm, change gamma at every step
// Reference: Phys. Rev. E 75, 046701 (2007).
  if(changeMode & (8+16+32))
  {
    long n = accept + reject;
    double dn;

    if (changeMode &  (8+16) == 8) n = accept;
    else if (changeMode &  (8+16) == 16) n = reject;
    else if (changeMode &  (8+16) == 8+16) n = accept + reject;

    if (changeMode & 32) dn = double(n) / double(n_walkers);
    else dn = double(n);

    gamma = 1.0 / dn;
    if (gamma > 1.0) gamma = 1.0;

    initKernel(kernelType, dosKernel);
    dosKernel.scale(gamma);
  }

// 1. write configuration
// 2. check histogram flatness
// 3. perform global update
  if (cycleCount >= updateCycle)
  {
    cycleCount = 0;
    // syncronizeGraphs(dos, histo);
    if (dos.getN() != histo.getN())
    {
      std::cout << "Histogramm size dosn't match DOS! Clearing histogram!\n";
      histo.setRangeAndClear(dos.getMinX(), dos.getMaxX(), dos.getN());
    }

    char stateFile[51];
    sprintf(stateFile, "WL1d.state%05d", walkerID);
    writeState(stateFile);

    if (!histogramUpdateMode)        // Normal Wang-Landau
    {
      // calculate minimum nonzero histogram
      // we look only at the histogram inside the energy interval that was actually sampled if we use kernel updates
      double hMin, hMax, hMean;

      if (kernelType == None)
      {
        hMean = histo.getMeanY();
        histo.getMinMaxY(hMin,hMax);
      }
      else 
      {
        hMean = histo.getMeanYInInterval(dos.getMinX() + dosKernel.getWidth(),
                                         dos.getMaxX() - dosKernel.getWidth());
        histo.getMinMaxYInInterval(dos.getMinX() + dosKernel.getWidth(),
                                   dos.getMaxX() - dosKernel.getWidth(),
                                   hMin, hMax);
      }

      // double currentFlatnessCriterion = double(hMin-hMax) / hMean;
      double currentFlatnessCriterion = double(hMin) / hMean;
      
      std::cout << "# acceptance ratio = " << double(accept) / double(accept+reject) << "\n";
      std::cout << "# current flatness = " << currentFlatnessCriterion << (changeMode & 4 ? " *":"") << "\n";
      std::cout << "# current histogram minimum = " << hMin << (changeMode & 2 ? " *":"") << "\n";
      std::cout << "# average accepted steps/bin since last gamma change = "
                << double(acceptSinceLastChange) / double(histo.getN())
                << (changeMode & 1 ? " *":"") << "\n";

      if (changeMode != 0 && changeMode < 8)
      {
        // perform global update
        if (globalUpdate.frequency > 0 && (globalUpdate.frequency*histo.getN()) < acceptSinceLastChange) 
        {
          char fn[256];
          snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
          writeDos(fn);
          globalUpdate.changes++;
          std::cout<<"# global update "<<globalUpdate.changes<<std::endl;
          performGlobalUpdate(dos,globalUpdate.kappa,globalUpdate.lambda,globalUpdate.omega);
          histo.clear();
          acceptSinceLastChange=0;
        }
        else if (((acceptSinceLastChange>=histo.getN()*updatesPerBin || !(changeMode &1)) &&
                  (hMin >= hMinimum || !(changeMode & 2)) &&
                  (currentFlatnessCriterion>flatnessCriterion || !(changeMode & 4)) ))
        {
          std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " is finished.\n";
          modificationFactorChanges++; // counter[2]

          // write dos;
          // we are using JSON as our file format
          char fn[256];
          snprintf(fn,255,"Dos_Global_%02d_Changes_%03d.jsn",globalUpdate.changes,modificationFactorChanges);
          writeDos(fn);
          // writeDos(dos_out_name.data());

          // clear the Histogram
          histo.clear();
          acceptSinceLastChange = 0;
 
          // change gamma
          dosKernel.scale(0.5);
          gamma = 0.5 * gamma;
          std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " begins.\n";
          if(statesFile != NULL)
	  {
            // char fn[256];
            snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
            sw.newFile(fn);
            sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
          }
        }
      }
    } 
    else            // histogram update mode != 0
    {
      if (acceptSinceLastChange >= histo.getN() * updatesPerBin)
      {
        acceptSinceLastChange = 0;
        for (int ix=0; ix<dos.getN(); ix++)
          dos[ix] += gamma * histo[ix];
        histo.clear();
        gamma = 0.5 * gamma;

        std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " begins.\n";
        if (statesFile != NULL)
        {
          char fn[256];
          snprintf(fn,255,"%s_%02d",statesFile,modificationFactorChanges);
          sw.newFile(fn);
          sw.writeHeader(gamma,n_walkers,n_spins,evecs_pointer);
        }
      }
    }
  }

  sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);
  lastAccepted[instance] = -1;
  lastAcceptedEnergy[instance] = position[instance];
  lastAcceptedEvec[  3*instance] = 0.0;
  lastAcceptedEvec[1+3*instance] = 0.0;
  lastAcceptedEvec[2+3*instance] = 0.0;
  // YingWai: Commented out the following 2 lines
  // so that RE statistics not to be confused with normal WL's   (Feb 15, 14)
  //accept++;
  //acceptSinceLastChange++;

  if (gamma < gammaFinal) return true;
  else return false;

}




#endif

