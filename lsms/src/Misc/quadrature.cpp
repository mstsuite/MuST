
#include "Real.hpp"
#include "calculateGaussLegendrePoints.hpp"

extern "C"
{
  void gauss_legendre_points_(double *x1, double *x2, double *x, double *w, int *n);
}

void gauss_legendre_points_(double *x1, double *x2, double *x, double *w, int *n)
{
  calculateGaussLegendrePoints(x, w, *n, *x1, *x2);
}

