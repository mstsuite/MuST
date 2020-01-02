#ifndef LSMS_LAPACK_H
#define LSMS_LAPACK_H

#include "Complex.hpp"

extern "C"
{
  void zgetrf_(int *m,int *n, Complex *a, int *lda, int*ipvt,int *info);
  void zgetri_(int *n, Complex *a, int *lda,int *ipvt, Complex *w, int *lw,
               int *info);
  double zlange_(char *c,int *m, int *n, Complex *a, int *lda, double*work);
}

#endif
