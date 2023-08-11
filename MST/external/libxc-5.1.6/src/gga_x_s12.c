/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_S12G         495 /* Swart 2012 GGA exchange                                    */
#define XC_HYB_GGA_X_S12H     496 /* Swart 2012 GGA hybrid exchange                             */

typedef struct {
  double A, B, C, D, E;
  double bx;
} gga_x_s12_params;

static void
gga_x_s12_init(xc_func_type *p)
{
  gga_x_s12_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_s12_params));
  params = (gga_x_s12_params *) (p->params);

  params->bx  = 1.0; /* we initialize it here */

  if(p->info->number == XC_HYB_GGA_X_S12H)
    xc_hyb_init_hybrid(p, 0.0);
}

#define S12G_N_PAR 5
static const char  *s12g_names[S12G_N_PAR]  = {"_A", "_B", "_C", "_D", "_E"};
static const char  *s12g_desc[S12G_N_PAR]   = {
  "A parameter",
  "B parameter",
  "C parameter",
  "D parameter",
  "E parameter"
};
static const double s12g_values[S12G_N_PAR] = {
  1.03842032, 1.757-1.03842032, 0.00403198, 0.00104596, 0.00594635
};

#define S12H_N_PAR 6
static const char  *s12h_names[S12H_N_PAR]  = {"_A", "_B", "_C", "_D", "_E", "_alpha"};
static const char  *s12h_desc[S12H_N_PAR]   = {
  "A parameter",
  "B parameter",
  "C parameter",
  "D parameter",
  "E parameter",
  "Fraction of exact exchange"
};
static const double s12h_values[S12H_N_PAR] = {
  1.02543951, 1.757-1.02543951, 0.00761554, 0.00211063, 0.00604672, 0.25
};

static void
s12h_set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_s12_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_s12_params *) (p->params);

  params->A    = get_ext_param(p, ext_params, 0);
  params->B    = get_ext_param(p, ext_params, 1);
  params->C    = get_ext_param(p, ext_params, 2);
  params->D    = get_ext_param(p, ext_params, 3);
  params->E    = get_ext_param(p, ext_params, 4);

  p->cam_alpha = get_ext_param(p, ext_params, 5);
  params->bx   = 1.0 - p->cam_alpha;
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_s12.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_s12g = {
  XC_GGA_X_S12G,
  XC_EXCHANGE,
  "Swart 2012 GGA exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {S12G_N_PAR, s12g_names, s12g_desc, s12g_values, set_ext_params_cpy},
  gga_x_s12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_s12h = {
  XC_HYB_GGA_X_S12H,
  XC_EXCHANGE,
  "Swart 2012 hybrid exchange",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {S12H_N_PAR, s12h_names, s12h_desc, s12h_values, s12h_set_ext_params},
  gga_x_s12_init, NULL,
  NULL, work_gga, NULL
};
