#include <stdlib.h>
#include "Real.hpp"

// mpicxx -I ../../../include/ plglmax_test.cpp plglmax_new.cpp plglmax.o zeroout.o -lgfortran


extern "C"
{
  void plglmax_(int *lend, Real*costheta, Real*plm);
}
void plglmax_new(int lmax, Real x, Real*plm);

int main(int argc, char *argv[])
{
  Real x;
  Real *plm, *plm_new;
  int lmax;

  lmax = atoi(argv[1]);
  x = atof(argv[2]);
  plm=(Real *)malloc((lmax+1)*(lmax+2)*sizeof(Real)/2);
  plm_new=(Real *)malloc((lmax+1)*(lmax+2)*sizeof(Real)/2);

  plglmax_new(lmax,x,plm_new);
  plglmax_(&lmax,&x,plm);

  printf("i plm[i] plm_new[i] |plm-plm_new|\n\n");

  for(int i=0; i<(lmax+1)*(lmax+2)/2; i++)
    printf("%d %lg %lg %lg\n",i,plm[i],plm_new[i],std::abs(plm[i]-plm_new[i]));


  free(plm); free(plm_new);
  return 0;
}
