// -*- mode: c++ -*-
#ifndef LSMS_WANG_LANDAU_H
#define LSMS_WANG_LANDAU_H

#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// we use BOOST for the random number generator
// #include <boost/random.hpp>
#include <random>
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
  void writeHeader(double lnF, int numWalk, int numSpin, double **spins, 
      bool isSpinSim = true, bool isOccSim = false, int **occup = NULL)
  {

    if( !writeFlag ) return;

    if( isSpinSim ) {
      of << lnF << " " << numWalk << " " << numSpin << std::endl;
      for(int i = 0; i < numWalk; i++) {
        of << i;
        for(int j = 0; j < 3*numSpin; j++) 
          of << " " << spins[i][j];
        of << std::endl;
      }
    }

    if( isOccSim ) return;
  }
  void writeChange(int iWalk, int numRet, int ispin, double *ev, double E)
  {
    if(writeFlag)
      of<<iWalk<<" "<<numRet<<" "<<ispin<<" "<<ev[0]<<" "<<ev[1]<<" "<<ev[2]<<" "<<E<<std::endl;
  }
  void writeChangeOcc(int iWalk, int numRet, int i, int j, int occ_i, int occ_j, double E)
  {
    if(!writeFlag) return;
      of<<iWalk<<" "<<numRet<<" ("<<i<<"<-->"<<j<<") "<<occ_i<<" "<< occ_j <<" "<<E<<std::endl;
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
  WL1dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
                    PotentialShifter &potentialShifter,
                    const char *init_file_name=NULL, const char *out_file_name=NULL,
                    const char *states_filename=NULL, int ** occ_p = NULL);

  bool determineAcceptance(int instance, double energy);
  bool determineAcceptance(int instance, double energy, double magnetization);

  bool updateHistogram(int instance, double *evecs, bool accepted);
/*
  bool updateHistogram(int instance, double *evecs, double energy, bool *accepted);
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool *accepted);
  bool updateHistogram(int instance, double *evecs, double energy) {bool h; return updateHistogram(instance, evecs, energy, &h);}
  bool updateHistogram(int instance, double *evecs) {std::cerr<<"Need energy for WL1dEvecGenerator\n"; exit(1);}
*/

  void generateEvec(int instance, double *evecs, bool accepted);
  //void generateEvec(int instance, double *evecs, double energy);
  void generatePotentialShift(int instance, double *potentialShifts, bool accepted);

  void initializeEvecAndPotentialShift(int inst, double *evecs, double *potentialShifts);
  void initializeEvec(int instance, double *evecs);

  void generateUnsampledEvec(int instance, double *evecs, double energy)
  {
    initializeEvec(instance, evecs); 
    //return false;
  }

  void startSampling(bool isSpinSim = true, bool isOccSim = false)
  { sw.writeHeader(gamma, n_walkers, n_spins, evecs_pointer, isSpinSim, isOccSim, occupancy_ptr); }

  void writeState(const char *name);
  void writeDos(const char *name);

  // Wang-Landau for occupancy variables
  MoveChoice_t selectMoveType(bool isSpinSim, bool isOccSim);
  void setAlloyClasses(const AlloyMixingDesc&, int* siteclass);
  void generateOccupancies(int instance, int *occ, bool acceptedOcc);
  void generateUnsampledOcc(int inst, int *occ);
  void initializeOccupancies(int inst, int *occ);
  void updateLog(int instance, double *evecs, int * occ, double energy, bool accepted, MoveChoice_t MoveChoice,
      bool isspin, bool isocc);

 private:
  int n_walkers;
  int n_spins;
  double ** evecs_pointer;
  int n_initialized_from_file;
  int n_initialized_from_file_occ;

  std::string dos_out_name;

  int stepsSinceLastHistogramUpdate;
  int numberOfUpdatesSinceLastBoost;
  int cycleCount;
  int modificationFactorChanges;

  // Random number generator and distribution:
  RNG rng;
  // boost::uniform_real<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd01; //(0.0,1.0);

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

  std::vector<int> visited;
  bool maskedHistogram;
  
  unsigned long accept, reject, acceptSinceLastChange;
  double flatnessCriterion;

  double gamma, gammaFinal;
  int flipPerUpdate, updateCycle;

  bool fixEnergyWindow;
  
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

  // changes to accomodate alloying (i.e. variable site occupancies)
  std::vector< std::pair<int,int> > lastSwapOcc;
  std::vector< std::pair<int,int> > lastAcceptedSwap;
  std::vector< std::pair<int,int> > lastAcceptedOcc;
  AlloyMixingDesc alloyDesc;
  std::vector<int> siteAlloyClass;
  std::vector<int> numSites_per_AlloyClass;
  int ** occupancy_ptr;

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
    
    /*  
    ev[2]=1.0-2.0*rnd(rng);
    // ev[2]=rnd11(rng);
    double phi=2.0*M_PI*rnd(rng);
    // double phi=rnd0pi(rng);
    double cos_theta=sqrt(1-ev[2]*ev[2]);
    ev[0]=cos_theta*cos(phi);
    ev[1]=cos_theta*sin(phi);
    */  
  }
#endif

  double minVxShift {0.0}, maxVxShift {0.0}, rangeVxShift {0.0};

  double inline randomPotentialShift()
  {
    return minVxShift + rangeVxShift * rnd01(rng);
  }
};


template<class RNG>
WL1dEvecGenerator<RNG>::WL1dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
                                          PotentialShifter &potentialShifter,
            const char *init_file_name, const char *out_file_name, 
                                          const char *states_filename, int **occ_p)
: sw(states_filename)
{
  if (potentialShifter.vSpinShiftFlag) {
    minVxShift = potentialShifter.minShift;
    maxVxShift = potentialShifter.maxShift;
    rangeVxShift = maxVxShift - minVxShift;
  }

  verbosity=0;

  changeMode = 4;

  globalUpdate.frequency=0;
  globalUpdate.changes=0;
  globalUpdate.kappa=1.0;
  globalUpdate.lambda=1.0;
  globalUpdate.omega=0.5;

  maskedHistogram=false;
  visited.clear();
  
  long nX=-1;
  histogramUpdateMode=false;
  updatesPerBin=100;
  double interval=0.01;
  double kernelWidth=0.1;
  double xMin=-std::numeric_limits<double>::max();// double xMin=-HUGE;
  double xMax=1.0;
  bool readPositions=false;
  n_spins=num_spins;
  n_walkers = num_instances;
  n_initialized_from_file = 0;
  n_initialized_from_file_occ = 0;
  evecs_pointer = ev_p;
  occupancy_ptr = occ_p;
  ref0.resize(n_walkers);  for(int i=0; i<n_walkers; i++) ref0[i]=-1; // ref0[i]=HUGE;
  ref1.resize(n_walkers);
  position.resize(n_walkers);
  magnetizationAtPosition.resize(n_walkers);
  out.resize(n_walkers);
  lastChange.resize(n_walkers);
  lastAccepted.resize(n_walkers);
  lastAcceptedEnergy.resize(n_walkers);
  lastAcceptedEvec.resize(3*n_walkers);

  oldPotentialShift.resize(n_walkers);

  lastChangePotentialShift.resize(n_walkers);
  lastAcceptedPotentialShiftIndex.resize(n_walkers);
  lastAcceptedPotentialShift.resize(n_walkers);

  lastSwapOcc.resize(n_walkers);
  lastAcceptedSwap.resize(n_walkers);
  lastAcceptedOcc.resize(n_walkers);
  for(int i=0; i<n_walkers; i++)
    lastSwapOcc[i].first = -1;

  fixEnergyWindow = false;
  
  // fprintf(stderr, "WARNING: Fixing random number seed!\n");
  // rng.seed(1);

  for(int i=0; i<n_walkers; i++)
  {
    lastAccepted[i]=-2;
    lastAcceptedPotentialShiftIndex[i] = -2;
    lastAcceptedPotentialShift[i] = 0.0;
  }
  for(int i=0; i<3*n_walkers; i++)  lastAcceptedEvec[i]=0.0;

  for(int i=0; i<n_walkers; i++) initializeEvec(i,evecs_pointer[i]);

  oldSpin.resize(3*n_walkers);

  statesFile=NULL;
  if(states_filename!=NULL)
    {
      statesFile=(char *)malloc(sizeof(char)*(1+strlen(states_filename)));
      strcpy(statesFile,states_filename);
    }

  numRetentions.resize(n_walkers);
  for(int i=0; i<n_walkers; i++) numRetentions[i]=0;

  /*
  nX = -1;
  xMin = -HUGE; xMax= 1.0; interval = 0.01; // (xMax-xMin)/double(nX);
  */
  dos_out_name="dos1d.out";
  stepsSinceLastHistogramUpdate=0;
  numberOfUpdatesSinceLastBoost=0;
  modificationFactorChanges=0;
  cycleCount=0;
  hMinimum= 1;     //   10
  acceptSinceLastChange=accept=reject=0;
  gammaFinal=1.e-6;
  flipPerUpdate=1; //  100
  updateCycle= 5*num_instances;  // 1000
  gamma = 1.0;
  flatnessCriterion = 0.75;
  updatesPerBin =100;

  kernelType=None;

  // special processing flags:
  int clearHistogram = 0;
  int setFirstWalkerToFM = 0;

  // dos = NULL;   // dos.resize(nX);
  // histo = NULL; // histo.resize(nX);
  dos.setDeltaAndClear(interval);
  histo.setDeltaAndClear(interval);

  if(init_file_name!=NULL && init_file_name[0]!=0)
  {
    std::string label, value;

    std::cout<<"Reading Wang-Landau configuration from: "<<init_file_name<<std::endl;

    dos_out_name=std::string(init_file_name)+".out";
    if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
    std::ifstream inp(init_file_name);
    std::ostringstream buff;

    std::string line;
    while(std::getline(inp,line)) 
      buff << line << std::endl;

    inp.close();

    std::string fileString = buff.str();
    const char* fileChars  = fileString.c_str();
    json_t *json_root=NULL;

    json_parse_document(&json_root, (char *)fileChars);

    if(json_root == NULL || json_root->type != JSON_OBJECT)
    {
      std::ostringstream message;
      std::cerr << "In WL1dEvecGenerator(" << init_file_name << ") parsing failed (bad format)\n";
      exit(1);
    }
  
    for(json_t *it = json_root->child; it != NULL; it=it->next)
    {
      std::string label = it->text;
      if(label=="xMin") xMin=atof(it->child->text);
      else if(label=="xMax") xMax=atof(it->child->text);
      else if(label=="fixEnergyWindow")
      {
        if( atoi(it->child->text) != 0 ) fixEnergyWindow = true;
      } 
      else if(label=="interval") interval=atof(it->child->text);
      else if(label=="kernelWidth") kernelWidth=atof(it->child->text);
      else if(label=="kernelType")
      {
        std::string strValue(it->child->text);
        kernelType=getKernelType(strValue);
      }
      else if(label=="gamma") gamma=atof(it->child->text);
      else if(label=="gammaFinal") gammaFinal=atof(it->child->text);
      else if(label=="nX")
      {
        nX=atoi(it->child->text);
  dos.setRangeAndClear(xMin,xMax,nX);
  histo.setRangeAndClear(xMin,xMax,nX);
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
      else if(label=="seed") rng.seed(atoi(it->child->text));
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
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          dos[j++]=atof(i->text);
        }
        if(j!=dos.getN()) {std::cout<<"ERROR #(dos) "<<j<<" != nX "<<dos.getN()<<std::endl; exit(1);}
      }
      else if(label=="histo")
      {
        json_t *a = it->child;
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          histo[j++]=atof(i->text);
        }
        if(j!=histo.getN()) {std::cout<<"ERROR #(histo) "<<j<<" != nX "<<histo.getN()<<std::endl; exit(1);}
      }
      else if(label=="moments")
      {
        json_t *a = it->child->child;
        int k = atoi(a->text);
        dos.setNumberOfMoments(k);
        a=a->next;
        int jj=0;
        for(json_t *j=a->child; j!=NULL; j=j->next)
        {
          dos.setNumberOfSamplesAtIdx(jj,atoi(j->text));
          jj++;
        }
        json_t *i=a->next;
        for(int kk=0; kk<k; kk++)
        {
          int jj=0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
          {
            dos.setMomentAtIdx(jj,kk,atof(j->text));
            jj++;
          }
          i=i->next;
        }
      }
      else if(label=="moments.k")
      {
        dos.setNumberOfMoments(atoi(it->child->text));
      }
      // additional read for occupancy scheme
      // mirrors 'evecs' read
      else if( label == "occupancies" ) {
        json_t *a = it->child;
        n_initialized_from_file_occ=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized_from_file_occ<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              occupancy_ptr[n_initialized_from_file_occ][k++]=atoi(j->text);
            }
            // initialize lastSwapOcc to -1 to indicate no previous swap
            lastSwapOcc[n_initialized_from_file_occ].first = -1;
           n_initialized_from_file_occ++;
           if(k != n_spins) {std::cout<<"ERROR #(occs) "<<k<<" != n_spins "<<n_spins<<std::endl; exit(1);}
          }
        } 
      }
      else if( label == "lastSwapOcc" ) {
        json_t *a = it->child;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              if( k == 0 ) lastSwapOcc[n_initialized].first  = atoi(j->text);
              if( k == 1 ) lastSwapOcc[n_initialized].second = atoi(j->text);
              k++;
            }
            n_initialized++;
            if(k!=2) {std::cout<<"ERROR #(lastSwapOcc) "<<k<<" != 2"<<std::endl; exit(1);}
          }
        } 
      }
      //
      else if(label=="evecs")
      {
        json_t *a = it->child;
        n_initialized_from_file=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized_from_file<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              evecs_pointer[n_initialized_from_file][k++]=atof(j->text);
            }
// initialize oldSpin and lastChange to point to site 0
            lastChange[n_initialized_from_file]=0;
            oldSpin[  3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][0];
            oldSpin[1+3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][1];
            oldSpin[2+3*n_initialized_from_file]=evecs_pointer[n_initialized_from_file][2];
//
            n_initialized_from_file++;
            if(k!=3*n_spins) {std::cout<<"ERROR #(evecs) "<<k<<" != 3*n_spins "<<3*n_spins<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="oldSpin")
      {
        json_t *a = it->child;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            int k=0;
            for(json_t *j=i->child; j!=NULL; j=j->next)
            {
              oldSpin[n_initialized*3+k++]=atof(j->text);
            }
            n_initialized++;
            if(k!=3) {std::cout<<"ERROR #(oldSpin) "<<k<<" != 3"<<std::endl; exit(1);}
          }
        }
      }
      else if(label=="lastChange")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers) lastChange[j++]=atoi(i->text);
          n_initialized++;
        }
      }
      else if(label=="position")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            position[j++]=atof(i->text);
            n_initialized++;
          }
        }
        readPositions=true;
      }
      else if(label=="magnetizationAtPosition")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            magnetizationAtPosition[j++]=atof(i->text);
            n_initialized++;
          }
        }
//        readMagnetizationAtPositions=true;
      }
      else if(label=="stepsSinceLastHistogramUpdate") stepsSinceLastHistogramUpdate=atoi(it->child->text);
      else if(label=="numberOfUpdatesSinceLastBoost") numberOfUpdatesSinceLastBoost=atoi(it->child->text);
      else if(label=="modificationFactorChanges") modificationFactorChanges=atoi(it->child->text);
      else if(label=="clearHistogram") clearHistogram=atoi(it->child->text);
      else if(label=="setFirstWalkerToFM") setFirstWalkerToFM=atoi(it->child->text);
      else if(label=="maskedHistogram") maskedHistogram=(atoi(it->child->text)!=0);
      else if(label=="visited")
      {
        maskedHistogram=true;
        json_t *a = it->child;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          visited.push_back(atoi(i->text));
        }
      }
      else std::cout<<"WARNING: unknown label: "<<label<<std::endl;
    }

    // set xMax and xMin or interval depending on nX:
    if(nX>0)
    {
      interval=(xMax-xMin)/double(nX);
      dos.setRange(xMin,xMax); histo.setRange(xMin,xMax);
    } else {
      dos.setDeltaAndClear(interval); histo.setDeltaAndClear(interval);
    }

    json_free_value(&json_root);
  }

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
    {
      ref0[i]=dos.idx(position[i]);
      lastAccepted[i] = -1;
      lastAcceptedEnergy[i] = position[i];
    }
      
  if(clearHistogram!=0)
  {
    std::cout<<"Clearing Wang-Landau Histogram.\n";
    histo.clear();
  }

  if(setFirstWalkerToFM!=0)
  { if(setFirstWalkerToFM<0)
    {
      std::cout<<"Setting all walkers to FM.\n";
      setFirstWalkerToFM=n_walkers;
    } // else {
      std::cout<<"Setting first "
               << std::min(setFirstWalkerToFM,n_walkers)
               << " walkers to FM.\n";
    // }
    n_initialized_from_file=std::max(setFirstWalkerToFM,n_initialized_from_file);
    for(int n=0; n<std::min(setFirstWalkerToFM,n_walkers); n++)
    {
      for(int i=0; i<n_spins; i++)
      {
        evecs_pointer[n][  3*i]=0.0;
        evecs_pointer[n][1+3*i]=0.0;
        evecs_pointer[n][2+3*i]=1.0;
      }
// initialize oldSpin and lastChange to point to site 0
     lastChange[n]=0;
     oldSpin[3*n+0]=evecs_pointer[n][0];
     oldSpin[3*n+1]=evecs_pointer[n][1];
     oldSpin[3*n+2]=evecs_pointer[n][2];

// initialize oldPotentialShift and lastChangePotentialShift to point to site 0
      lastChangePotentialShift[n] = 0;
      oldPotentialShift[n] = 0.0;
    }
  }

  if(maskedHistogram && !fixEnergyWindow)
  {
    std::cout<<"masked Histogram requires fixed Energy Window!\n";
    exit(2);
  }
  if(maskedHistogram && kernelType!=None)
  {
    std::cout<<"masked Histogram only works with kernel type None!.\n";
    exit(2);
  }
  if(maskedHistogram && (visited.size()!=nX))
  {
    std::cout<<"visited mask resized.\n";
    visited.resize(nX);
    for(size_t i=0; i<visited.size(); i++) 
     if(histo[i]>=1.0) visited[i]=1; else visited[i]=0;
  }
    
  if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
  std::cout<<"Wang-Landau output will be written to: "<<dos_out_name<<std::endl;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::initializeEvecAndPotentialShift(int inst, double *evecs, double *potentialShifts)
{
  for (size_t j=0; j<n_spins; j++)
    potentialShifts[j] = 0.0;

  initializeEvec(inst, evecs);
}

template<class RNG>
void WL1dEvecGenerator<RNG>::initializeEvec(int inst, double *evecs)
{
  bool firstFerromagnetic=true;
  if(inst>=n_initialized_from_file)
    {
      if(firstFerromagnetic)
        if(inst==0)
    for(size_t j=0; j<3*n_spins; j+=3)
      {evecs[j]=0.0; evecs[j+1]=0.0; evecs[j+2]=1.0;}
        else if(inst==1)
          for(size_t j=0; j<3*n_spins; j+=3)
            {evecs[j+2]= (j/3)%2 ? 1.0 : -1.0 ; evecs[j+1]=0.0; evecs[j]=0.0;}
        else
    for(size_t j=0; j<3*n_spins; j+=3)
      random_evec_1(&evecs[j]);
      else
        for(size_t j=0; j<3*n_spins; j+=3)
          random_evec_1(&evecs[j]);
    }

  out[inst]=false;
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
    ofile.precision(15);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << ",\n";
    ofile<<"\"xMax\" : " << dos.getMaxX() << ",\n";
    if( fixEnergyWindow )
      ofile<<"\"fixEnergyWindow\" : 1,\n";
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
    if(maskedHistogram)
    {
      ofile<<"\"maskedHistogram\" : 1,\n";
      ofile<<"\"visited\" : ["<<std::endl;
      for(int i=0; i<visited.size(); i++) ofile<<visited[i]<<((i==visited.size()-1)?"\n":",\n");
      ofile<<"],\n";
    }
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
    // additional output for occupancies
    ofile << "\"occupancies\" : [\n";
    for(int i = 0; i < n_walkers; i++) {
      ofile << "[ ";
      int j = 0, stride = 20;
      while( j < n_spins ) {
        ofile << occupancy_ptr[i][j];
        if( j+1 != n_spins ) ofile << ", ";
        if( j+1 % stride == 0 ) ofile << "\n";
        j++;
      }
      ofile << " ]";
      if( i+1 != n_walkers ) ofile << ", ";
      ofile << "\n";
    }
    ofile << "],\n";
    ofile << "\"lastSwapOcc\" : [\n";
    for(int i = 0; i < n_walkers; i++) {
      ofile << "[ ";
      ofile << lastSwapOcc[i].first;
      ofile << ", ";
      ofile << lastSwapOcc[i].second;
      ofile << " ]";
      if( i+1 != n_walkers ) ofile << ", ";
      ofile << "\n";
    }
    ofile << "],\n";
    //
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
    ofile<<"\"oldSpin\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<<"[ "<<oldSpin[3*i]<<", "<<oldSpin[3*i+1]
             <<", "<<oldSpin[3*i+2];
      ofile<<((i==n_walkers-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"lastChange\" : [\n";
    for(int i=0; i<n_walkers; i++)
    {
      ofile<< lastChange[i]<<((i==n_walkers-1)?"\n":",\n");
    }
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
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(15);
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
    if(maskedHistogram)
    {
      ofile<<"\"maskedHistogram\" : 1,\n";
      ofile<<"\"visited\" : ["<<std::endl;
      for(int i=0; i<visited.size(); i++) ofile<<visited[i]<<((i==visited.size()-1)?"\n":",\n");
      ofile<<"],\n";
    }
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::determineAcceptance(int instance, double energy)
{
  return determineAcceptance(instance, energy, 0.0);
}


template<class RNG>
bool WL1dEvecGenerator<RNG>::determineAcceptance(int instance, double energy, double magnetization)
{

  bool accept_step;
  // +++++++++ Added by Odbadrakh and Don on Aug 11, 2010
  // +++ If a flip of spin results in energy in the same bin as previous accepted 
  // +++ state, it is accepted all the time. We change it by accepting if the
  // +++ new state has lower density of states according to a linear spline of DOS
  
  double dos_differ;
  double energy_differ;
  double to_go_or_not;

  // +++++++++ end of modification. Follow the new variables to see the actual changes 
  // energy between xMin and xMax ?
  int grid = dos.idx(energy);

  // dos.test();

  stepsSinceLastHistogramUpdate++; // counter[0]++

  // record initial energy for step out
  if (lastAccepted[instance] == -2)
  {
    position[instance] = lastAcceptedEnergy[instance] = energy;
    lastAccepted[instance] = -1;
  }
  
  if (grid < 0 || grid > dos.getN()-1)
  {

    // guard for fixed window
    if( fixEnergyWindow ) {

      // default action is reject
      accept_step = false;

      // not in window yet
      if( ref0[instance] < 0 || ref0[instance] > dos.getN()-1 ) { 
        out[instance] = true;
        if( grid < 0 && energy - position[instance] > 0.0) accept_step = true;
        if( grid > dos.getN()-1 && energy - position[instance] < 0.0) accept_step = true;
      }

      // already inside window
      else
        out[instance] = false;

      return accept_step;
    }

    dos.extendTo(energy - dosKernel.getWidth()); histo.extendTo(energy - histoKernel.getWidth());
    dos.extendTo(energy + dosKernel.getWidth()); histo.extendTo(energy + histoKernel.getWidth());
    grid = dos.idx(energy);
    // we need to adjust the grid point of the walker positions (ref0) for all walkers
    for (long ii=0; ii<n_walkers; ii++)
      if (ref0[ii] >= 0) ref0[ii] = dos.idx(position[ii]);
  }
  numRetentions[instance]++;
  out[instance] = false;
  ref1[instance] = grid;
  if (ref0[instance] < 0 || ref0[instance] > dos.getN()-1 )
  {
      accept_step = true;
      ref0[instance] = ref1[instance];
      position[instance] = energy;
      magnetizationAtPosition[instance] = magnetization;
  }
  else 
  {
    if (abs(ref1[instance] - ref0[instance]) < 1 ) 
    {
// Actual change made by Odbadrakh, Aug 30, 2010
      dos_differ = dos[ref0[instance]-1] - dos[ref0[instance]];
      energy_differ = energy - lastAcceptedEnergy[instance];
      to_go_or_not = dos_differ * energy_differ;
      if (to_go_or_not >= 0.0)        // accepts all downhill changes
      {
        accept_step = true;
        ref0[instance] = ref1[instance];
        position[instance] = energy;
        magnetizationAtPosition[instance] = magnetization;
      } 
      else       // uphill moves
      {
        if(rnd(rng) < exp(to_go_or_not/dos.getDelta()))
        { //std::cout<<" dos "<<instance<<"  weight "<<dos.getDelta()<<std::endl;
          accept_step = true;
          ref0[instance] = ref1[instance];
          position[instance] = energy;
          magnetizationAtPosition[instance] = magnetization;
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
        ref0[instance] = ref1[instance];
        position[instance] = energy;
        magnetizationAtPosition[instance] = magnetization;
      }
      else
      {
        if(rnd(rng) < exp(dos[ref0[instance]] - dos[ref1[instance]]))
        {
          accept_step = true;
          ref0[instance] = ref1[instance];
          position[instance] = energy;
          magnetizationAtPosition[instance] = magnetization;
        }
        else 
        {
          accept_step = false;
        }
      }
      
    }

 }

// End of change made on Aug 30, 2010
  if (verbosity > 2)
    std::cout << "WangLandau 1d EvecGenerator step "
              << modificationFactorChanges << ":" << numberOfUpdatesSinceLastBoost << ":"
              << stepsSinceLastHistogramUpdate << " nX=" << dos.getN()
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
bool WL1dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, bool accepted)
{

// Update histogram and DOS
  if(stepsSinceLastHistogramUpdate >= flipPerUpdate)
  {
    stepsSinceLastHistogramUpdate = 0;       // counter[0]
    numberOfUpdatesSinceLastBoost++;         // counter[1]
    cycleCount++;
    // modify the DoS
    if(!out[instance])
    {
      // addKernel(dos,dosKernel,energy);
      if (!histogramUpdateMode)             // standard Wang-Landau
      {
        addKernel(dos, dosKernel, position[instance]);
        dos.addMomentsAtIdx(dos.idx(position[instance]), magnetizationAtPosition[instance]);
        // if(accepted) addKernel(histo,histoKernel,position[instance]);
        addKernel(histo, histoKernel, position[instance]);
        if(maskedHistogram) visited[histo.idx(position[instance])]=1;
        // addKernel(dos, dosKernel, energy);
        // addKernel(histo, histoKernel, energy);
      }
      else
      {
        // addKernel(dos, nullKernel, position[instance]);
        if(accepted) addKernel(histo, histoKernel, position[instance]);
        if(maskedHistogram) visited[histo.idx(position[instance])]=1;
      }
    } 
    else
    {
      if(!fixEnergyWindow) {
        std::cerr << "ATTENTION: We should never reach this place in WL1dEvecGenerator!\n";
        exit(1);
      }
    }
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
      std::cout << "Histogramm size dosn't match DOS! Clearing histogramm!\n";
      histo.setRangeAndClear(dos.getMinX(), dos.getMaxX(), dos.getN());
    }

    writeState("WL1d.state");

    if (!histogramUpdateMode)        // Normal Wang-Landau
    {
      // calculate minimum nonzero histogram
      // we look only at the histogram inside the energy interval that was actually sampled if we use kernel updates
      double hMin, hMax, hMean;

      if (kernelType == None)
      {
        if(maskedHistogram)
        {
          hMean = histo.getMeanYMasked(visited);
          histo.getMinMaxYMasked(hMin,hMax,visited);
        }
        else
        {
          hMean = histo.getMeanY();
          histo.getMinMaxY(hMin,hMax);
        }
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
      
      std::cout << "# acceptence ratio = " << double(accept) / double(accept+reject) << "\n";
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
    // Suffian: Moved to routine updateLog(...)
    /* sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);
    lastAccepted[instance] = lastChange[instance];
    lastAcceptedEnergy[instance] = position[instance];
    //lastAcceptedEnergy[instance] = energy;
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]]; */
    accept++;
    acceptSinceLastChange++;
  }
  else reject++;

  if (gamma < gammaFinal) return true;
  else return false;

}

template<class RNG>
void WL1dEvecGenerator<RNG>::updateLog(int instance, double *evecs, int * occ, double energy, bool accepted, MoveChoice_t MoveChoice,
    bool isSpinSim, bool isOccupancySim) {

  // use static variables to keep track of the last move type
  static bool firstrun = true;
  static std::vector<MoveChoice_t> lastAcceptedMove;
  static FILE* filep;

  // initialize static variables
  if( firstrun ) {
    lastAcceptedMove.resize(n_walkers);

    if( isOccupancySim ) {
      filep = fopen("occ_vs_energies","w");
      fprintf(filep,"# site occupancies, energy\n");
    }
    firstrun = false;
  }

  // Suffian: my preference is to print all configurations and energy
  if( isOccupancySim ) {
    for(int i = 0; i < n_spins; i++)
      fprintf(filep, "%1d", occ[i]);
    fprintf(filep, "   %20.15f", energy);
    fprintf(filep, "   %5d", instance);
    if( accepted ) fprintf(filep, "   accepted");
    fprintf(filep, "\n");
    fflush(filep);
  }

  // return unless state changed
  if( !accepted ) return;

  // print output corresponding to previous state
  if( lastAcceptedMove[instance] == SpinMove )
    sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);
  else
    sw.writeChangeOcc(instance, numRetentions[instance], 
        lastAcceptedSwap[instance].first,
        lastAcceptedSwap[instance].second,
        lastAcceptedOcc[instance].first, 
        lastAcceptedOcc[instance].second, 
        lastAcceptedEnergy[instance]); 

  // remember the new configuration
  lastAcceptedEnergy[instance] = position[instance];
  if( MoveChoice == SpinMove ) {
    lastAccepted[instance] = lastChange[instance];
    //lastAcceptedEnergy[instance] = energy;
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]]; 
  }
  else {
    lastAcceptedSwap[instance].first  = lastSwapOcc[instance].first;
    lastAcceptedSwap[instance].second = lastSwapOcc[instance].second;
    lastAcceptedOcc[instance].first  = occ[ lastSwapOcc[instance].first ];
    lastAcceptedOcc[instance].second = occ[ lastSwapOcc[instance].second ];
  }

  lastAcceptedMove[instance] = MoveChoice;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::generatePotentialShift(int instance, double *potentialShifts, bool accepted)
{

  if (accepted)
  {
    lastChangePotentialShift[instance] = int(rnd(rng)*n_spins);
    oldPotentialShift[instance] = potentialShifts[lastChangePotentialShift[instance]];
    potentialShifts[lastChangePotentialShift[instance]] = randomPotentialShift();
  } else {
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

// additional routines for Wang-Landau for alloying
// swapping atoms preserves concentration

template<class RNG>
MoveChoice_t WL1dEvecGenerator<RNG>::selectMoveType(bool isSpinSim, bool isOccSim) {

  // choose between spin and occupancy equally
  // is this a good choice??
  if( isSpinSim && isOccSim ) {
    int c = int(rnd(rng)*2);
    if( c == 0 ) 
      return SpinMove;
    else
      return OccupancyMove;
  }

  if( isSpinSim )
    return SpinMove;
  if( isOccSim )
    return OccupancyMove;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::setAlloyClasses(const AlloyMixingDesc& alloyDesc, int *siteclass) {
  this->alloyDesc = alloyDesc;

  siteAlloyClass.resize(n_spins);
  for(int i = 0; i < n_spins; i++)
    siteAlloyClass[i] = siteclass[i];

  // remember number of sites per alloy class 
  int nclasses = alloyDesc.size();
  numSites_per_AlloyClass.resize(nclasses);
  int *nsites = numSites_per_AlloyClass.data();
  for(int i = 0; i < nclasses; i++)
    nsites[i] = 0;
  for(int i = 0; i < n_spins; i++)
    nsites[ siteAlloyClass[i] ]++;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::generateOccupancies(int instance, int *occ, bool acceptedOcc) {

  int temp, atom_i, atom_j;

  // if previous move exists and not accepted 
  // revert occupancies
  if( !acceptedOcc && lastSwapOcc[instance].first >= 0 ) {
    atom_i = lastSwapOcc[instance].first;
    atom_j = lastSwapOcc[instance].second;

    temp = occ[atom_i];
    occ[atom_i] = occ[atom_j];
    occ[atom_j] = temp;

    /*
    printf("Master: Reverted occupancies:");
    for(int i = 0; i < n_spins; i++)
      printf("%d ",occ[i]);
    printf("\n"); */
  }
  else
    numRetentions[instance] = 0;

  // pick an alloy class
  int ac = int(rnd(rng)*alloyDesc.size());

  // pick two distinct atoms within mixing class
  // here we assume there are few alloy classes
  do{ atom_i = int(rnd(rng)*n_spins); } 
    while( siteAlloyClass[atom_i] != ac );
  do{ atom_j = int(rnd(rng)*n_spins); } 
    while( siteAlloyClass[atom_j] != ac || occ[atom_j] == occ[atom_i] );

  // swap atoms
  temp = occ[atom_i];
  occ[atom_i] = occ[atom_j];
  occ[atom_j] = temp;

  // remember which atoms were swapped
  lastSwapOcc[instance].first  = atom_i;
  lastSwapOcc[instance].second = atom_j;
}

template<class RNG>
void WL1dEvecGenerator<RNG>::generateUnsampledOcc(int inst, int *occ) {

  // determine number atoms of each component type given concentration
  std::vector< std::vector<int> > atomsLeft;
  atomsLeft.resize(alloyDesc.size());
  for(int i = 0; i < atomsLeft.size(); i++) {
    atomsLeft[i].resize(alloyDesc[i].size());
    for(int j = 0; j < atomsLeft[i].size(); j++) {

      // how many sites for this alloy class and component?
      double nfrac = alloyDesc[i][j].conc * numSites_per_AlloyClass[i];
      double error = fmod(nfrac, 1.0);

      // ensure concentrations commensurate with number of sites
      const double tol = 1.e-06;
      if( tol < error && error < 1.0-tol ) {
        printf("error: alloy concentrations are incommensurate with supercell\n");
        printf("alloy class = %d, component = %d, concentration = %f\n", 
            i+1, j+1, alloyDesc[i][j].conc);
        printf("desired sites = %f\n", nfrac);
        exit(1);
      }

      atomsLeft[i][j] = int(nfrac + 0.5);
    }
  }

  // distribute atoms across sites
  for(int i = 0; i < n_spins; i++) {

     // pick an atom from the right class
     int type, ac = siteAlloyClass[i];
     do { type = int( rnd(rng) * alloyDesc[ac].size() ); }
       while( atomsLeft[ac][type] == 0 );

     // assign to site
     occ[i] = type; atomsLeft[ac][type]--;
  }
}
 
template<class RNG>
void WL1dEvecGenerator<RNG>::initializeOccupancies(int inst, int *occ) {

  // Warning: do nothing if walker occupancies read from file
  if( inst < n_initialized_from_file_occ )
    return;

  generateUnsampledOcc(inst, occ);
}


#endif

