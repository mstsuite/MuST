#include "associatedLegendreFunction.hpp"
#include "plglmax_device.hpp"

#include <stdio.h>
#include <cuda.h>

#define NDLM 28
#define LMAX 3

__global__ void wrapAssociatedLegendre(double x, int lmax, double *plmOut)
{
  __shared__ double plm[NDLM];
  associatedLegendreFunctionNormalizedDevice(x,lmax,plm);
  for ( int i = threadIdx.x; i <= LMAX; i += blockDim.x)
    plmOut[i] = plm[i];
}
int main()
{
  double costheta = 0.3;
  int lmax = LMAX;
  double plmC[NDLM];
  double *plmG = 0;
  cudaMallocManaged(&plmG,NDLM * sizeof(double));

  associatedLegendreFunctionNormalized<double>(costheta,lmax,plmC);

  wrapAssociatedLegendre<<<1,1>>>(costheta,lmax,plmG);
  cudaDeviceSynchronize();
  
  for(int i = 0; i <= lmax; i++)
  {
    if(plmC[i] != plmG[i])
    {
      fprintf(stderr, "%d: %lf != %lf\n", i, plmC[i], plmG[i]);
      return -1;
    }
  }

  cudaFree(plmG);
  return 0;
}
