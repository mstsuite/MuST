/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019      Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_C_AM05          135 /* Armiento & Mattsson 05 correlation             */

typedef struct{
  double alpha, gamma;
} gga_c_am05_params;

static void
gga_c_am05_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_am05_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_am05.c"
#include "work_gga.c"

#define AM05_N_PAR 2
static const char  *am05_names[AM05_N_PAR]  = {"_alpha", "_gamma"};
static const char  *am05_desc[AM05_N_PAR]   = {"alpha", "gamma"};
static const double am05_values[AM05_N_PAR] =
  {2.804L, 0.8098L};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_am05 = {
  XC_GGA_C_AM05,
  XC_CORRELATION,
  "Armiento & Mattsson 05",
  XC_FAMILY_GGA,
  {&xc_ref_Armiento2005_085108, &xc_ref_Mattsson2008_084714, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {AM05_N_PAR, am05_names, am05_desc, am05_values, set_ext_params_cpy},
  gga_c_am05_init, NULL,
  NULL, work_gga, NULL
};
