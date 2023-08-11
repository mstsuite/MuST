/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019      Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_AK13  56 /* Armiento & Kuemmel 2013 */

typedef struct{
  double B1, B2;
} gga_x_ak13_params;

static void
gga_x_ak13_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ak13_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ak13.c"
#include "work_gga.c"

#define N_PAR 2
static const char  *names[N_PAR]  = {"_B1", "_B2"};
static const char  *desc[N_PAR]   = {"B1", "B2"};
/* B1 = 3*muGE/5 + 8 pi/15   B2 = muGE - B1 */
static const double par_ak13[N_PAR] =
  {1.74959015598863046792081721182, -1.62613336586517367779736042170};

double xc_gga_ak13_get_asymptotic (double homo)
{
  return xc_gga_ak13_pars_get_asymptotic(homo, par_ak13);
}

double xc_gga_ak13_pars_get_asymptotic (double homo, const double *ext_params)
{
  double Qx, aa, aa2, factor;
  double ak13_B1;

  ak13_B1 = ext_params[0];

  Qx = sqrt(2.0)*ak13_B1/(3.0*CBRT(3.0*M_PI*M_PI));

  aa  = X_FACTOR_C*Qx;
  aa2 = aa*aa;

  factor = (homo < 0.0) ? -1.0 : 1.0;

  return (aa2/2.0)*(1.0 + factor*sqrt(1.0 - 4.0*homo/aa2));
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ak13 = {
  XC_GGA_X_AK13,
  XC_EXCHANGE,
  "Armiento & Kuemmel 2013",
  XC_FAMILY_GGA,
  {&xc_ref_Armiento2013_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_ak13, set_ext_params_cpy},
  gga_x_ak13_init, NULL,
  NULL, work_gga, NULL
};
