/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MBRXC_BG  696 /* Modified Becke-Roussel for band gaps - cuspless hole */

GPU_FUNCTION static double
mbrxc_x_Q(double x, void *_rhs)
{
  double rhs;
  rhs = *((double *)_rhs);

  return pow(1.0 + x, 5.0/3.0)*exp(-2.0*x/3.0) - rhs*(x - 3.0);
}

GPU_FUNCTION
double xc_mgga_x_mbrxc_get_x(double Q)
{
  double rhs, tol, x1, x2;

  tol = 5e-12;
  if(fabs(Q) < 5e-12)
    return 3.0;

  /* build right-hand side of the non-linear equation
     Remember we use a different definition of tau */
  rhs = pow(32.0*M_PI, 2.0/3.0)/(6.0*Q);

  /* starting interval */
  if(rhs > 0.0) {
    /* I checked that the solution is always in this interval */
    x1 = 3.0;
    x2 = 2.0/rhs + 3.0;
  }else{
    x2 = 3.0;
    x1 = -1.0;
  }

  return xc_math_brent(mbrxc_x_Q, x1, x2, tol, 500, &rhs);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_mbrxc_bg.c"
#include "work_mgga.c"


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mbrxc_bg = {
  XC_MGGA_X_MBRXC_BG,
  XC_EXCHANGE,
  "Modified Becke-Roussel for band gaps - cuspless hole",
  XC_FAMILY_MGGA,
  {&xc_ref_Patra2019_045147, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS | XC_FLAGS_DEVELOPMENT,
  1.0e-12,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
