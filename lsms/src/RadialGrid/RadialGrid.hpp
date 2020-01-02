#ifndef LSMS_RGRID_H
#define LSMS_RGRID_H

#include <vector>
#include "Real.hpp"

class RadialGrid {
public:
  inline RadialGrid() : N(0), jmt(0), jws(0), h(0.0) {}
  int N,jmt,jws;
  Real h;
  std::vector<Real> r_mesh,x_mesh; 
};

void generateRadialGrid(RadialGrid * g, Real x0, Real h, int N, int jmt, int jws);

#endif
