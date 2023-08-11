/*
 Copyright (C) 2015 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MVS          257 /* MVS exchange of Sun, Perdew, and Ruzsinszky */

typedef struct {
  double e1, c1, k0, b;
} mgga_x_mvs_params;

static void
mgga_x_mvs_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_mvs_params));
}

#define MVS_N_PAR 4
static const char  *mvs_names[MVS_N_PAR]  = {"_e1", "_c1", "_k0", "_b"};
static const char  *mvs_desc[MVS_N_PAR]   = {
  "e1 parameter",
  "c1 parameter",
  "k0 parameter",
  "b parameter"
};
static const double mvs_values[MVS_N_PAR] = {-1.6665, 0.7438, 0.174, 0.0233};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_mvs.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mvs = {
  XC_MGGA_X_MVS,
  XC_EXCHANGE,
  "MVS exchange of Sun, Perdew, and Ruzsinszky",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2015_685, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {MVS_N_PAR, mvs_names, mvs_desc, mvs_values, set_ext_params_cpy},
  mgga_x_mvs_init, NULL,
  NULL, NULL, work_mgga,
};
