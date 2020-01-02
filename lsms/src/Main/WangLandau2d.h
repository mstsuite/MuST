// -*- mode: c++ -*-
#ifndef LSMS_WANG_LANDAU_2d_H
#define LSMS_WANG_LANDAU_2d_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
// we use BOOST for the random number generator
// #include <boost/random.hpp>
#include <random>
#include "../../mjson/json.h"
#include "EvecGenerator.h"
#include "Graph2d.hpp"

class StatesWriter2d
{
public:
  StatesWriter2d(const char *filename=NULL)
  {
    if(filename==NULL) writeFlag=false;
    else {writeFlag=true; of.open(filename);
      of.setf(std::ios::scientific,std::ios::floatfield);
      of.precision(8);}
  }
  ~StatesWriter2d() {if(writeFlag) of.close();}
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

// template<class RNG = boost::mt19937>
template<class RNG = std::mt19937>
class WL2dEvecGenerator : public EvecGenerator
{
 public:
  WL2dEvecGenerator(int num_spins, int num_instances, double **ev_p,
                    const char *init_file_name = NULL,
                    const char *out_file_name = NULL,
                    const char *states_filename = NULL);

  bool determineAcceptance(int instance, double energy);
  bool determineAcceptance(int instance, double energy, double magnetization);

  bool updateHistogram(int instance, double *evecs, bool accepted);
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool accepted);
/*
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization, bool *accepted);
  bool updateHistogram(int instance, double *evecs, double energy, double magnetization)
  {
    bool h;
    return updateHistogram(instance, evecs, energy, magnetization, &h);
  }
  bool updateHistogram(int instance, double *evecs)
  {
    std::cerr << "Need energy for WL2dEvecGenerator\n"; 
    exit(1);
  }
*/

  void generateEvec(int instance, double *evecs, bool accepted);
  //void generateEvec(int instance, double *evecs, double energy);

  void initializeEvec(int instance, double *evecs);

  void generateUnsampledEvec(int instance, double *evecs, double energy)
  {
    initializeEvec(instance, evecs);
    //return false;
  }

  void startSampling(void)
  { sw.writeHeader(gamma, n_walkers, n_spins, evecs_pointer); }

  void writeState(const char *name);

 private:
  int n_walkers;
  int n_spins;
  double ** evecs_pointer;
  int n_initialized_from_file;

  std::string dos_out_name;

  int stepsSinceLastHistogramUpdate;
  int numberOfUpdatesSinceLastBoost;
  int cycleCount;
  int modificationFactorChanges;

  // Random number generator and distribution:
  RNG rng;
  // boost::uniform_real<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> rnd; //, rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);

  /*
  // Histogramm and dos:
  double xMin, xMax, interval;
  int nX;
  double *dos; // std::vector<double> dos;
  int *histo; // std::vector<int> histo;
  int hMinimum;
  */

  double hMinimum;
  Graph2d<double,double> dos, histo;
  Kernel2d<double,double> dosKernel, histoKernel;
  KernelType kernelType;

  unsigned long accept, reject;
  double flatnessCriterion;

  double gamma, gammaFinal;
  int flipPerUpdate, updateCycle;

  // instance specific:
  std::vector<double> positionX, positionY;
  std::vector<bool> out;
  std::vector<int> lastChange;
  std::vector<int> lastAccepted;
  std::vector<double> lastAcceptedEnergy;
  std::vector<double> lastAcceptedMagnetization;
  std::vector<double> lastAcceptedEvec;
  std::vector<double> oldSpin;  // oldSpin[instance*3 + {x=0, y=1, z=2}]

  std::vector<int> numRetentions;

  char *statesFile;
  StatesWriter2d sw;


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

};

template<class RNG>
WL2dEvecGenerator<RNG>::WL2dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
					  const char *init_file_name, const char *out_file_name, const char *states_filename)
  : sw(states_filename)
{
  long nX=-1;
  double intervalEnergy=0.1;
  double intervalMagnetization=0.1;
  double kernelWidthEnergy=1.0;
  double kernelWidthMagnetization=1.0
;
  double xMin=-std::numeric_limits<double>::max();
  double xMax=1.0;
  n_spins=num_spins;
  n_walkers = num_instances;
  n_initialized_from_file = 0;
  evecs_pointer = ev_p;
  positionX.resize(n_walkers);
  positionY.resize(n_walkers);
  out.resize(n_walkers);
  lastChange.resize(n_walkers);
  lastAccepted.resize(n_walkers);
  lastAcceptedEnergy.resize(n_walkers);
  lastAcceptedMagnetization.resize(n_walkers);
  lastAcceptedEvec.resize(3*n_walkers);
  for(int i=0; i<n_walkers; i++)
  {
    lastAccepted[i]=-2;
  }
  for(int i=0; i<3*n_walkers; i++)  lastAcceptedEvec[i]=0.0;

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
  dos_out_name="dos2d.out";
  stepsSinceLastHistogramUpdate=0;
  numberOfUpdatesSinceLastBoost=0;
  modificationFactorChanges=0;
  cycleCount=0;
  hMinimum= 1;     //   10
  accept=reject=0;
  gammaFinal=1.e-6;
  flipPerUpdate=4; //  100
  updateCycle= 1000; // 5*num_instances;  // 1000
  gamma = 1.0;
  flatnessCriterion = 0.75;

  kernelType=Epanechnikov;

  // dos = NULL;   // dos.resize(nX);
  // histo = NULL; // histo.resize(nX);
  dos.setDeltaAndClear(intervalEnergy,intervalMagnetization);
  histo.setDeltaAndClear(intervalEnergy,intervalMagnetization);
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
      std::cerr << "In WL2dEvecGenerator(" << init_file_name << ") parsing failed (bad format)\n";
      exit(1);
    }
  
    for(json_t *it = json_root->child; it != NULL; it=it->next)
    {
      std::string label = it->text;
      if(label=="xMin") xMin=atof(it->child->text);
      else if(label=="xMax") xMax=atof(it->child->text);
      else if(label=="intervalEnergy") intervalEnergy=atof(it->child->text);
      else if(label=="intervalMagnetization") intervalMagnetization=atof(it->child->text);
      else if(label=="kernelWidthEnergy") kernelWidthEnergy=atof(it->child->text);
      else if(label=="kernelWidthMagnetization") kernelWidthMagnetization=atof(it->child->text);
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
	dos.setXRangeAndClear(xMin,xMax,nX);
	histo.setXRangeAndClear(xMin,xMax,nX);
      }
      else if(label=="flipPerUpdate") flipPerUpdate=atoi(it->child->text);
      else if(label=="updateCycle") updateCycle=atoi(it->child->text);
      else if(label=="cycleCount") cycleCount=atoi(it->child->text);
      else if(label=="flatnessCriterion") flatnessCriterion=atof(it->child->text);
      else if(label=="seed") rng.seed(atoi(it->child->text));
      else if(label=="accept") accept=atol(it->child->text);
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
        int ix=0;
        int iy;
        long Ny;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          Ny=atoi(i->text); i=i->next;
          double minY=atof(i->text); i=i->next;
          dos.setYminAndClear(ix,minY,Ny);
          iy=0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
            dos[ix][iy++]=atof(j->text);
          ix++;
        }
        if(ix!=dos.getNx()) {std::cout<<"ERROR #(dos) "<<ix<<" != nX "<<dos.getNx()<<std::endl; exit(1);}
      }
      else if(label=="histo")
      {
        json_t *a = it->child;
        int ix=0;
        int iy;
        long Ny;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          Ny=atoi(i->text); i=i->next;
          double minY=atof(i->text); i=i->next;
          histo.setYminAndClear(ix,minY,Ny);
          iy=0;
          for(json_t *j=i->child; j!=NULL; j=j->next)
            histo[ix][iy++]=atof(j->text);
          ix++;
        }
        if(ix!=histo.getNx()) {std::cout<<"ERROR #(histo) "<<ix<<" != nX "<<histo.getNx()<<std::endl; exit(1);}
      }
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
            positionX[j]=atof(i->text);
            i=i->next;
            positionY[j]=atof(i->text);
            n_initialized++;
          }
        }
      }
      else if(label=="stepsSinceLastHistogramUpdate") stepsSinceLastHistogramUpdate=atoi(it->child->text);
      else if(label=="numberOfUpdatesSinceLastBoost") numberOfUpdatesSinceLastBoost=atoi(it->child->text);
      else if(label=="modificationFactorChanges") modificationFactorChanges=atoi(it->child->text);
      else std::cout<<"WARNING: unknown label: "<<label<<std::endl;
    }

    // set xMax and xMin or interval depending on nX:
    if(nX>0)
    {
      intervalEnergy=(xMax-xMin)/double(nX);
      dos.setXRange(xMin,xMax); histo.setXRange(xMin,xMax);
    } else {
      dos.setDeltaAndClear(intervalEnergy,intervalMagnetization);
      histo.setDeltaAndClear(intervalEnergy,intervalMagnetization);
    }

    json_free_value(&json_root);
  }
  dosKernel.setWidthAndClear(intervalEnergy, intervalMagnetization,
                             kernelWidthEnergy, kernelWidthMagnetization);
  histoKernel.setWidthAndClear(intervalEnergy, intervalMagnetization,
                             kernelWidthEnergy, kernelWidthMagnetization);
  initKernel(kernelType,dosKernel);
  dosKernel.scale(gamma/dosKernel(0.0,0.0));
  initKernel(kernelType,histoKernel);
  histoKernel.scale(1.0/histoKernel(0.0,0.0));
//  histoKernel.scale(interval);

  if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
  std::cout<<"Wang-Landau output will be written to: "<<dos_out_name<<std::endl;
}

template<class RNG>
void WL2dEvecGenerator<RNG>::initializeEvec(int inst, double *evecs)
{
  if(inst>=n_initialized_from_file)
    {
      if(inst==0)
	for(size_t j=0; j<3*n_spins; j+=3)
	  {evecs[j]=0.0; evecs[j+1]=0.0; evecs[j+2]=1.0;}
      else if(inst==1)
        for(size_t j=0; j<3*n_spins; j+=3)
          {evecs[j+2]= (j/3)%2 ? 1.0 : -1.0 ; evecs[j+1]=0.0; evecs[j]=0.0;}
      else
	for(size_t j=0; j<3*n_spins; j+=3)
	  random_evec_1(&evecs[j]);
    }

  out[inst]=false;
}

template<class RNG>
void WL2dEvecGenerator<RNG>::writeState(const char* name)
{
  std::ofstream ofile(name);
  if(ofile)
  {
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(8);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << dos.getMinX() << ",\n";
    ofile<<"\"xMax\" : " << dos.getMaxX() << ",\n";
    ofile<<"\"nX\" : " << dos.getNx() << ",\n";
    ofile<<"\"intervalMagnetization\" : " << dos.getDeltaY() << ",\n";
    ofile<<"\"kernelWidthEnergy\" : " << dosKernel.getWidthX() << ",\n";
    ofile<<"\"kernelWidthMagnetization\" : " << dosKernel.getWidthY() << ",\n";
    std::string kernelName;
    getKernelName(kernelType,kernelName);
    ofile<<"\"kernelType\" : \"" << kernelName << "\",\n";
    ofile<<"\"gamma\" : " << gamma << ",\n";
    ofile<<"\"accept\" : " << accept <<",\n";
    ofile<<"\"reject\" : " << reject <<",\n";
    ofile<<"\"gammaFinal\" : " << gammaFinal << ",\n";
    ofile<<"\"flatnessCriterion\" : "<< flatnessCriterion <<",\n";
    ofile<<"\"flipPerUpdate\" : " << flipPerUpdate << ",\n";
    ofile<<"\"updateCycle\" : " << updateCycle << ",\n";
    ofile<<"\"dos\" : ["<<std::endl;
    for(int ix=0; ix<dos.getNx(); ix++)
    {
      ofile<<dos.getNy(ix)<<", "<<dos.getMinY(ix)<<", ["<<std::endl;
      for(int iy=0; iy<dos.getNy(ix); iy++)
        ofile<<dos[ix][iy]<<((iy==dos.getNy(ix)-1)?"\n":",\n");
      ofile<<((ix==dos.getNx()-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
    ofile<<"\"histo\" : ["<<std::endl;
    for(int ix=0; ix<dos.getNx(); ix++) 
    { 
      ofile<<histo.getNy(ix)<<", "<<dos.getMinY(ix)<<", ["<<std::endl;
      for(int iy=0; iy<histo.getNy(ix); iy++)
        ofile<<histo[ix][iy]<<((iy==histo.getNy(ix)-1)?"\n":",\n");
      ofile<<((ix==histo.getNx()-1)?"]\n":"],\n");
    }
    ofile<<"],\n";
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
      ofile<< positionX[i]<<", "<<positionY[i]<<((i==n_walkers-1)?"\n":",\n");
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
bool WL2dEvecGenerator<RNG>::determineAcceptance(int instance, double energy)
{
  return determineAcceptance(instance, energy, 0.0);
}


template<class RNG>
bool WL2dEvecGenerator<RNG>::determineAcceptance(int instance, double energy, double magnetization)
{

  bool accept_step;
 
  // energy between xMin and xMax ?
  long gridX = dos.idxX(energy);
  // we have to postpone geting the y position until we have checked the validity of gridX
  //  int gridY = dos.idxY(gridX, magnetization);

  stepsSinceLastHistogramUpdate++; // counter[0]++

  // record initial position for step out
  if (lastAccepted[instance] == -2)
  {
    positionX[instance] = lastAcceptedEnergy[instance] = energy;
    positionY[instance] = lastAcceptedMagnetization[instance] = magnetization;
    lastAccepted[instance] = -1;
  }

  // out of energy range
  if (gridX < 0 || gridX > dos.getNx()-1)
  {
    dos.extendToX(energy - dosKernel.getWidthX()); histo.extendToX(energy - histoKernel.getWidthX());
    dos.extendToX(energy + dosKernel.getWidthX()); histo.extendToX(energy + histoKernel.getWidthX());
    gridX = dos.idxX(energy);  // gridY=dos.idxY(gridX, magnetization);
  }
  long gridY = dos.idxY(gridX, magnetization);
  // out of magnetization range
  if (gridY < 0 || gridY > dos.getNy(gridX)-1)
  {
    dos.extendToY(gridX, magnetization - dosKernel.getWidthY());
    dos.extendToY(gridX, magnetization + dosKernel.getWidthY());
    histo.extendToY(gridX, magnetization - histoKernel.getWidthY());
    histo.extendToY(gridX, magnetization + histoKernel.getWidthY());
    gridY = dos.idxY(gridX, magnetization);
  }
 
  numRetentions[instance]++;
  out[instance] = false;
  // ref1X[instance] = gridX; ref1Y[instance] = gridY;
  long refX = dos.idxX(positionX[instance]);
  long refY = dos.idxY(refX,positionY[instance]);

  if (lastAccepted[instance] < 0)
  {
    accept_step = true;
    positionX[instance] = energy;
    positionY[instance] = magnetization;
  }
  else
  {
    if (dos[gridX][gridY] <= dos[refX][refY])
    {
      accept_step = true;
      positionX[instance] = energy;
      positionY[instance] = magnetization;
    } 
    else 
    {
      if(rnd(rng) < exp(dos[refX][refY] - dos[gridX][gridY]))
      {
        accept_step = true;
        positionX[instance] = energy;
	positionY[instance] = magnetization;
      } 
      else 
      {
        accept_step = false;
      }
    } 

  }

  if (verbosity > 0)
    std::cout << "WangLandau 2d EvecGenerator step "
              << modificationFactorChanges << ":" << numberOfUpdatesSinceLastBoost << ":"
              << stepsSinceLastHistogramUpdate << " nX=" << dos.getNx()
              << " [" << dos.getMinX() << ", " << dos.getMaxX() << "] "
              << (accept_step ? "accepted" : "rejected")
              << " Energy = " << energy << ", Magnetization = " << magnetization
	      << ", Instance " << instance << std::endl;

  return accept_step;

}


template<class RNG>
bool WL2dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, bool accepted)
{
  std::cerr << "Need energy and magnetization for updateHistogram in WL2dEvecGenerator.\n";
  std::cerr << "Either wl_lsms or updateHistogram needs to be amended!\n"; 
  exit(1);
}


template<class RNG>
bool WL2dEvecGenerator<RNG>::updateHistogram(int instance, double *evecs, double energy, double magnetization, bool accepted)
{

// Update histogram and DOS
  if (stepsSinceLastHistogramUpdate >= flipPerUpdate)
  {
    stepsSinceLastHistogramUpdate = 0;       // counter[0]
    numberOfUpdatesSinceLastBoost++;         // counter[1]
    cycleCount++;
    // modify the DoS
    if(!out[instance])
    {
      //addKernel(dos, dosKernel, positionX[instance], positionY[instance]);
      //addKernel(histo, histoKernel, positionX[instance], positionY[instance]);
      addKernel(dos, dosKernel, energy, magnetization);
      addKernel(histo, histoKernel, energy, magnetization);
    } 
    else
    {
      std::cerr << "ATTENTION: We should never reach this place in WL2dEvecGenerator!\n";
      exit(1);
    }
  }

// 1. write configuration
// 2. check histogram flatness
  if (cycleCount >= updateCycle)
  {
    cycleCount = 0;
    // syncronizeGraphs(dos,histo);
    writeState("WL2d.state");
    // calculate minimum nonzero histogram
    /*
    int hMin = histo[0];
    int hMax = histo[0];
    double hMean = 0.0;
    for(int i=1; i<nX; i++)
    {
      hMean += double(histo[i]);
      if(histo[i]>=hMax) hMax=histo[i];
      if(histo[i]<hMin) hMin=histo[i];
    }
    hMean /= double(nX);
    double currentFlatnessCriterion = double(hMax-hMin)/hMean;
    
    std::cout <<"# accepence ratio = "<<double(accept)/double(accept+reject)<<"\n";
    std::cout <<"# current flatness = "<<currentFlatnessCriterion<<"\n";
    */

    //if(hMin >= hMinimum && currentFlatnessCriterion<=flatnessCriterion)
    // we look only at the histogram inside the energy interval that was actually sampled
    double hMin = histo.getMinValWithBorders(dosKernel.getWidthX(), dosKernel.getWidthY());
    std::cout << "# acceptence ratio = " << double(accept) / double(accept+reject) << "\n";
    std::cout << "# current histogram minimum = " << hMin << "\n";
    if (hMin >= hMinimum)
    {
      std::cout << "# level " << modificationFactorChanges << " with gamma = " << gamma << " is finished.\n";
      modificationFactorChanges++; // counter[2]

      // write dos;
      // we are using JSON as our file format
      writeState(dos_out_name.data());
 
      // else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";

      // clear the Histogram
      histo.clear();

      // change gamma
      dosKernel.scale(0.5);
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


  if (accepted)
  {
    sw.writeChange(instance, numRetentions[instance], lastAccepted[instance], &lastAcceptedEvec[3*instance], lastAcceptedEnergy[instance]);
    lastAccepted[instance] = lastChange[instance];
    lastAcceptedEnergy[instance] = energy;
    lastAcceptedMagnetization[instance] = magnetization;
    lastAcceptedEvec[  3*instance] = evecs[  3*lastChange[instance]];
    lastAcceptedEvec[1+3*instance] = evecs[1+3*lastChange[instance]];
    lastAcceptedEvec[2+3*instance] = evecs[2+3*lastChange[instance]];
    accept++;
  }
  else reject++;

  if (gamma<gammaFinal) return true;
  else return false;

}

template<class RNG>
void WL2dEvecGenerator<RNG>::generateEvec(int instance, double *evecs, bool accepted)
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


#endif

