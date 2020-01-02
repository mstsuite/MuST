#ifndef LSMS_RADIALPOT_H
#define LSMS_RADIALPOT_H

#include "RadialGrid.hpp"
#include "Real.hpp"
#include "Matrix.hpp"

typedef class RadialPotential {
public:
  inline RadialPotential(RadialGrid *_g) {g=_g; vr.resize(g->N,2);}
  inline RadialPotential() {g=new RadialGrid;}
  inline void sync() {vr.resize(g->N,2);}
  RadialGrid *g;
  Matrix<Real> vr;
};

#endif
