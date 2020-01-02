#include <cmath>
#include "RadialGrid.hpp"

void generateRadialGrid(RadialGrid *g, Real x0, Real h, int N, int jmt, int jws)
{
  g->x_mesh.resize(N);
  g->r_mesh.resize(N);
  g->N=N;
  g->h=h;
  g->jmt=jmt;
  g->jws=jws;
  for(int i=0; i<N; i++)
  {
    g->x_mesh[i]=x0+i*h;
    g->r_mesh[i]=std::exp(g->x_mesh[i]);
  }
}


