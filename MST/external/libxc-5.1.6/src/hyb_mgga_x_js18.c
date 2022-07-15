/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_MGGA_X_JS18       705 /* a screened version of TM */

#define JS18_N_PAR 2
static const char  *js18_names[JS18_N_PAR]  = {"_a", "_omega"};
static const char  *js18_desc[JS18_N_PAR]   = {
  "Fraction of short-range Hartree-Fock exchange",
  "Range separation parameter"
};
static const double par_js18[JS18_N_PAR] = {0.1, 0.33};

static void
hyb_mgga_x_js18_init(xc_func_type *p)
{
  xc_hyb_init_sr(p, 0.0, 0.0);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/hyb_mgga_x_js18.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_js18 = {
  XC_HYB_MGGA_X_JS18,
  XC_EXCHANGE,
  "Jana and Samal 2018, screened range-separated TM exchange",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Jana2018_8999, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-14,
  {JS18_N_PAR, js18_names, js18_desc, par_js18, set_ext_params_cpy_cam_sr},
  hyb_mgga_x_js18_init, NULL,
  NULL, NULL, work_mgga
};
