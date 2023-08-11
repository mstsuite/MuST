/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/* Note: Do not forget to add a correlation (LDA) functional to the
   LB94.

   Note 2: The 160 value is hardcoded in xc.h and libxc_master.F90 to
   define XC_GGA_XC_LB to keep backwards compatibility.

*/
#define XC_GGA_X_LB  160 /* van Leeuwen & Baerends */
#define XC_GGA_X_LBM 182 /* van Leeuwen & Baerends modified*/

typedef struct{
  double alpha;
  double beta;
  double gamma;
} gga_x_lb_params;

/************************************************************************
  Calculates van Leeuwen Baerends functional
************************************************************************/

static void
gga_lb_init(xc_func_type *p)
{
  gga_x_lb_params *params;

  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_lb_params));
  params = (gga_x_lb_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_LB:
    params->alpha = 1.0;
    params->beta  = 0.05;
    params->gamma = 1.0;
    break;
  case XC_GGA_X_LBM:
    params->alpha = 1.19;
    params->beta  = 0.01;
    params->gamma = 1.0;
    break;
  }
}

#define XC_NO_EXC
#include "decl_gga.h"
#include "maple2c/gga_vxc/gga_x_lb.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lb = {
  XC_GGA_X_LB,
  XC_EXCHANGE,
  "van Leeuwen & Baerends",
  XC_FAMILY_GGA,
  {&xc_ref_vanLeeuwen1994_2421, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_lb_init, NULL,
  NULL, work_gga, NULL
};


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lbm = {
  XC_GGA_X_LBM,
  XC_EXCHANGE,
  "van Leeuwen & Baerends modified",
  XC_FAMILY_GGA,
  {&xc_ref_Schipper2000_1344, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_lb_init, NULL,
  NULL, work_gga, NULL
};

