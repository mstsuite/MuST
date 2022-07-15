/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_HLTAPW          699 /* Meta-GGAized PW */

typedef struct{
  double ltafrac;
} mgga_c_ltapw_params;

static void
mgga_c_ltapw_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_ltapw_params));
}

#define N_PAR 1
static const char  *names[N_PAR]   = {"_ltafrac"};
static const char  *desc[N_PAR]    = {"Fraction of LTA density"};
static const double ltapw_values[N_PAR]  = {0.5};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_ltapw.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_hltapw = {
  XC_MGGA_C_HLTAPW,
  XC_CORRELATION,
  "Half-and-half meta-LDAized PW correlation by Lehtola and Marques",
  XC_FAMILY_MGGA,
  {&xc_ref_Lehtola2021_943, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, names, desc, ltapw_values, set_ext_params_cpy},
  mgga_c_ltapw_init, NULL,
  NULL, NULL, work_mgga,
};
