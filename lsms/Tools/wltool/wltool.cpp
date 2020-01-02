#include <stdio.h>
#include <fstream>
#include "../../libjson/json.hpp"

void printMenu(const char*in, const char *on)
{
  printf("Input filename :%s\n",in);
  printf("Output filename:%s\n",on);
  printf("Please choose one of the following commands:\n\n");
  printf(" r - read file           w - write data\n");
  printf(" s - print summary\n");
  printf(" x - extract dos and histogram\n");
  printf(" h - clear histogram     g - step modification factor\n");
  printf("\n q - quit\n");
}

void writeData(JSON::Value &d, const char *on)
{
  std::ofstream os(on);
  os.unsetf(std::ios::floatfield);
  os.setf(std::ios::scientific);
  os.precision(12);
  os<<d;
}

void printSummary(JSON::Value &d)
{
  double gamma=(double)d["gamma"];
  printf("gamma = %lf\n",gamma);
  int changeMode=(int)d["changeMode"];
  printf("changeMode = %d  (",changeMode);
  if(changeMode==0) printf("no modification factor change");
  else if(changeMode < 8)
  {
    if(changeMode & 1)
    {
      printf("minimum of average accepted steps per bin");
      if(changeMode>1) printf(", ");
    }
    if(changeMode & 2)
    {
      printf("histogram minimum");
      if(changeMode>3) printf(", ");
    }
    if(changeMode & 4) printf("flatness test");
  } else if(changeMode < 64) {
    printf("1/t method");
  } else printf("UNKNOWN changeMode");
  printf(")\n");
  printf("modificationFactorChanges=%d\n",(int)d["modificationFactorChanges"]);
  printf("\n");
}

void extractDos(JSON::Value &d, const char *on)
{
  size_t nX=d["nX"].as_int();
  double xMin=(double)d["xMin"];
  double xMax=(double)d["xMax"];
  double xDelta=(xMax-xMin)/((double)nX);

  JSON::Array histo=(JSON::Array)d["histo"];
  JSON::Array dos=(JSON::Array)d["dos"];
  FILE *outf=fopen(on,"w");
  for(size_t i=0; i<nX; i++)
  {
    fprintf(outf,"%4d  %12.8lf  %12.8lf  %12.8lf\n",i,
            xMin+xDelta*(double)i, (double)dos[i], (double)histo[i]);
  }
  fclose(outf);
}

void clearHistogram(JSON::Value &d)
{
// get histogram:
  if(d["histo"].type()==JSON::ARRAY)
  {
    size_t s=((JSON::Array)d["histo"]).size();
    JSON::Array h;
    for(size_t i=0; i<s; i++)
      h.push_back(0);

    d["histo"]=h;
  }
}

void stepGamma(JSON::Value &d)
{
  double gamma=(double)d["gamma"];
  gamma=gamma/2.0;
  d["gamma"]=gamma;
  d["modificationFactorChanges"]=1+(int)d["modificationFactorChanges"];
  d["acceptSinceLastChange"]=0;

  clearHistogram(d);
}

void interact(void)
{
  char cmd=' ';

  JSON::Value data;
  char inName[128];
  char outName[128];
  double gamma=-1.0;

  inName[0]=0; outName[0]=0;

  printMenu(inName,outName);

  while(cmd!='q')
  {
    cmd=getchar();
    switch(cmd)
    {
    case 'r': printf("Filename:"); scanf("%s",inName);
      data=parse_file(inName); break;
    case 'w': printf("Filename:"); scanf("%s",outName);
      writeData(data,outName); break;
    case 'x': printf("Filename:"); scanf("%s",outName);
      extractDos(data,outName); break;
    case 's': printSummary(data); break;
    case 'h': clearHistogram(data); break;
    case 'g': stepGamma(data); break;
    default: printMenu(inName,outName);
    }
  }
}

int main(int argc, char **argv)
{
  if(argc==1) // interactive mode
  {
    interact();
  } else {
    printf("please use the interactive mode.\n");
  }
}
