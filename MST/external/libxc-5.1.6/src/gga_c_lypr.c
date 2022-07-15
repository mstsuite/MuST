/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_C_LYPR       624  /* Short-range LYP of Ai, Fang and Su */

typedef struct{
  double a;
  double b;
  double c;
  double d;
  double m1;
  double m2;
  double omega;
} gga_c_lypr_params;

#define LYPR_N_PAR 7
static const char  *lypr_names[LYPR_N_PAR]  = {
  "_a", "_b", "_c", "_d", "_m1", "_m2", "_omega"
};
static const char  *lypr_desc[LYPR_N_PAR]   = {
  "Parameter a",
  "Parameter b",
  "Parameter c",
  "Parameter d",
  "Parameter m1",
  "Parameter m2",
  "Range-separation parameter"
};
/* The values of m1 and m2 are from the authors' NWChem implementation */
static const double lypr_values[LYPR_N_PAR] =
  {0.04918, 0.132, 0.2533, 0.349, 0.35/2.29, 2.0/2.29, 0.33};

void xc_gga_c_lypr_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_lypr_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_lypr.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_lypr = {
  XC_GGA_C_LYPR,
  XC_CORRELATION,
  "Short-range LYP by Ai, Fang, and Su",
  XC_FAMILY_GGA,
  {&xc_ref_Ai2021_1207, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {LYPR_N_PAR, lypr_names, lypr_desc, lypr_values, set_ext_params_cpy},
  xc_gga_c_lypr_init, NULL,
  NULL, work_gga, NULL
};
