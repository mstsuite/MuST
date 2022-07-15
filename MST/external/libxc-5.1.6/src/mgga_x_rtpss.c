/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_RTPSS          299 /* Revised TPSS exchange by Garza, Bell and Head-Gordon */

typedef struct {
  double b, c, e, kappa, mu;
} mgga_x_rtpss_params;

static void
mgga_x_rtpss_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_rtpss_params));
}

#define RTPSS_N_PAR 5
static const char  *rtpss_names[RTPSS_N_PAR]  = {"_b", "_c", "_e", "_kappa", "_mu"};
static const char  *rtpss_desc[RTPSS_N_PAR]   = {
  "b", "c", "e",
  "Asymptotic value of the enhancement function",
  "Coefficient of the 2nd order expansion"
};
static const double rtpss_values[RTPSS_N_PAR] = {
  0.40, 1.59096, 1.537, 0.8040, 0.21951
};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_rtpss.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rtpss = {
  XC_MGGA_X_RTPSS,
  XC_EXCHANGE,
  "TPSS for surface adsorption",
  XC_FAMILY_MGGA,
  {&xc_ref_Garza2018_3083, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {RTPSS_N_PAR, rtpss_names, rtpss_desc, rtpss_values, set_ext_params_cpy},
  mgga_x_rtpss_init, NULL,
  NULL, NULL, work_mgga,
};
