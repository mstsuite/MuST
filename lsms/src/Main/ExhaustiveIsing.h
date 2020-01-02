/*
YingWai  (Dec 5, 13):
1. The original generateEvec is split into "determineAcceptance, updateHistogram 
   and generateEvec" in EvecGenerator.h, but they are NOT implemented in here
   (ExhaustiveIsing.h)
2. Here, generateEvec is renamed to updateHistogram, nothing has been changed 
   there otherwise.
3. Main program (wl_lsms.cpp) was rewritten in a way where ExhaustiveIsing 
   was NOT considered. Use with caution!
*/

// -*- mode: c++ -*-
#ifndef EXHAUSTIVEISING_H
#define EXHAUSTIVEISING_H

#define ISING 1

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
// Try MJSON instead:
#include "../../mjson/json.h"
#include "EvecGenerator.h"

class ExhaustiveIsing1dEvecGenerator : public EvecGenerator
{
 public:
  ExhaustiveIsing1dEvecGenerator(int num_spins, int num_instances, double **ev_p,
                                 const char *init_file_name = NULL,
                                 const char *out_file_name = NULL);
  bool updateHistogram(int instance, double *evecs, double energy);
  bool updateHistogram(int instance, double *evecs) {std::cerr<<"Need energy for ExhaustiveIsing1dEvecGenerator\n"; exit(1);}
  void initializeEvec(int instance, double *evecs);
  void writeState(const char *name);
 private:
  int n_walkers;
  int n_spins;
  double ** evecs_pointer;
  int n_initialized_from_file;

  std::string dos_out_name;

  int config_number;

  int stepsSinceLastHistogramUpdate;
  int numberOfUpdatesSinceLastBoost;
  int cycleCount;
  int modificationFactorChanges;

  // Histogramm and dos:
  double xMin, xMax, interval;
  int nX;
  double *dos; // std::vector<double> dos;
  int *histo; // std::vector<int> histo;
  int hMinimum;
  unsigned long accept, reject;
  double flatnessCriterion;

  double gamma, gammaFinal;
  int flipPerUpdate, updateCycle;

  // instance specific:
  std::vector<double> ref0, ref1, position;
  std::vector<bool> out;
  std::vector<int> lastChange;
  std::vector<double> oldSpin;  // oldSpin[instance*3 + {x=0, y=1, z=2}]


  void inline config_from_number(int n, double *ev)
  {
    int j=n;
    for(int i=0; i<n_spins; i++)
    {
      ev[3*i]=ev[3*i+1]=0.0;
      ev[3*i+2]= (j%2 == 0) ? 1.0 : -1.0;
      j=j>>1;  
    }
  }
};

ExhaustiveIsing1dEvecGenerator::ExhaustiveIsing1dEvecGenerator(int num_spins, int num_instances, double ** ev_p,
                                     const char *init_file_name, const char *out_file_name)
{
  const double Huge = std::numeric_limits<double>::max();
  n_spins=num_spins;
  n_walkers = num_instances;
  n_initialized_from_file = 0;
  evecs_pointer = ev_p;
  ref0.resize(n_walkers);  for(int i=0; i<n_walkers; i++) ref0[i]=Huge;
  ref1.resize(n_walkers);
  position.resize(n_walkers);
  out.resize(n_walkers);
  lastChange.resize(n_walkers);
  oldSpin.resize(3*n_walkers);

  nX = -1;
  xMin = -Huge; xMax= 1.0; interval = 0.01; // (xMax-xMin)/double(nX);
  dos_out_name="dos1d.out";
  stepsSinceLastHistogramUpdate=0;
  numberOfUpdatesSinceLastBoost=0;
  modificationFactorChanges=0;
  gamma = 1.0;
  flatnessCriterion = 0.75;
  config_number = 0;

  dos = NULL;   // dos.resize(nX);
  histo = NULL; // histo.resize(nX);
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
      else if(label=="interval") interval=atof(it->child->text);
      else if(label=="gamma") gamma=atof(it->child->text);
      else if(label=="gammaFinal") gammaFinal=atof(it->child->text);
      else if(label=="configurationNumber") config_number=atoi(it->child->text);
      else if(label=="nX")
      {
        nX=atoi(it->child->text);
        if(dos!=NULL) free(dos);
        if(histo!=NULL) free(histo);
        dos=(double *)calloc(nX,sizeof(double));
        histo=(int *)calloc(nX,sizeof(int));
      }
      else if(label=="flipPerUpdate") flipPerUpdate=atoi(it->child->text);
      else if(label=="updateCycle") updateCycle=atoi(it->child->text);
      else if(label=="flatnessCriterion") flatnessCriterion=atof(it->child->text);
//*/
      else if(label=="dos")
      {
        json_t *a = it->child;
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          dos[j++]=atof(i->text);
        }
        if(j!=nX) {std::cout<<"ERROR #(dos) "<<j<<" != nX "<<nX<<std::endl; exit(1);}
      }
      else if(label=="histo")
      {
        json_t *a = it->child;
        int j=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          histo[j++]=atoi(i->text);
        }
        if(j!=nX) {std::cout<<"ERROR #(histo) "<<j<<" != nX "<<nX<<std::endl; exit(1);}
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
      else if(label=="ref")
      {
        json_t *a = it->child;
        int j=0;
        int n_initialized=0;
        for(json_t *i=a->child; i!=NULL; i=i->next)
        {
          if(n_initialized<n_walkers)
          {
            ref0[j]=atof(i->text);
            i=i->next;
            ref1[j++]=atof(i->text);
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
      interval=(xMax-xMin)/double(nX);
    } else {
      nX=-1; xMin=-std::numeric_limits<double>::max(); xMax=1.0;
    }

    json_free_value(&json_root);
  }

  if(out_file_name!=NULL && out_file_name[0]!=0) dos_out_name=out_file_name;
  std::cout<<"Wang-Landau output will be written to: "<<dos_out_name<<std::endl;
}

void ExhaustiveIsing1dEvecGenerator::initializeEvec(int inst, double *evecs)
{
  config_from_number(config_number++, evecs);

  out[inst]=false;
}

void ExhaustiveIsing1dEvecGenerator::writeState(const char* name)
{
  std::ofstream ofile(name);
  if(ofile)
  {
    ofile.setf(std::ios::scientific,std::ios::floatfield);
    ofile.precision(8);
    ofile<<"{\n";
    ofile<<"\"xMin\" : " << xMin << ",\n";
    ofile<<"\"xMax\" : " << xMax << ",\n";
    ofile<<"\"nX\" : " << nX << ",\n";
    ofile<<"\"configurationNumber\" : " << config_number << ",\n";
    ofile<<"\"dos\" : ["<<std::endl;
    for(int i=0; i<nX; i++) ofile<<dos[i]<<((i==nX-1)?"\n":",\n");
    ofile<<"],\n";
    ofile<<"\"histo\" : ["<<std::endl;
    for(int i=0; i<nX; i++) ofile<<histo[i]<<((i==nX-1)?"\n":",\n");
    ofile<<"]\n";
    ofile<<"}\n";
    ofile.close();
  } else std::cerr<<"# CAUTION: DoS output file could not be opened!\n";
} 


bool ExhaustiveIsing1dEvecGenerator::updateHistogram(int instance, double *evecs, double energy)
{
  // energy between xMin and xMax ?
  int grid = int(nX*(energy-xMin)/(xMax-xMin));
  // if(energy < xMin || energy > xMax)
  if(grid < 0)
  {
    if(dos==NULL)
    {
      // first energy
      nX=1;
      // lets put the energy in the middle of the interval
      xMin = energy - 0.5*interval;
      xMax = xMin + interval;
      dos = (double *)calloc(1,sizeof(double));
      histo = (int *)calloc(1,sizeof(int));
      grid=0;
    } else {
      // we need -grid more entries below xMin
      double *p1 = (double *)calloc(nX-grid,sizeof(double));
      int *p2 = (int *)calloc(nX-grid,sizeof(int));
      for(int i=0; i<nX; i++) {p1[i-grid]=dos[i]; p2[i-grid]=histo[i];}
      free(dos); free(histo);
      dos=p1; histo=p2;
      nX = nX-grid;
      xMin += double(grid)*interval;
      grid = 0;
    }
  } else if (grid > nX -1) {
    // we need grid - nX + 1 entries above xMax
    int n = grid - nX + 1;
    double *p1 = (double *)calloc(nX+n,sizeof(double));
    int *p2 = (int *)calloc(nX+n,sizeof(int));
    for(int i=0; i<nX; i++) {p1[i]=dos[i]; p2[i]=histo[i];}
    free(dos); free(histo);
    dos=p1; histo=p2;
    nX += n;
    xMax += double(n)*interval;
    // grid already has the correct value
  }

  stepsSinceLastHistogramUpdate++; // counter[0]++

  if(grid < 0 || grid > nX -1)
  {
    std::cerr<< "ATTENTION: The code should never reach this point, where grid is out of bounds.\n";
    exit(1);
    out[instance]=true;
  } else {
    out[instance]=false;
    ref1[instance] = dos[grid];
  }

  std::cout<<"Exhaustive Ising 1d EvecGenerator step "
           <<config_number<<":"<<numberOfUpdatesSinceLastBoost<<":"
           <<stepsSinceLastHistogramUpdate<<" nX="<<nX<<" ["<<xMin<<", "<<xMax<<"] "
           <<std::endl;

    numberOfUpdatesSinceLastBoost++; // counter[1]
    cycleCount++;
    // modify the DoS
    if(!out[instance])
    {
      dos[grid] += 1.0;
      histo[grid]++;
    } else {
      std::cerr<<"ATTENTION: We should never reach this place in ExhaustiveIsing1dEvecGenerator!\n";
      exit(1);
    }

  if(cycleCount >= updateCycle)
  {
    cycleCount=0;
    writeState("WL1d.state");
  }

  if(gamma<gammaFinal) return true;

  if(config_number>= 1<<n_spins) return true;
  config_from_number(config_number,evecs);
  return false;
}




#endif

