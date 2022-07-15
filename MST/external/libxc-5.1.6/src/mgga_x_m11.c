/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_X_M11       297 /* M11 hybrid exchange functional from Minnesota        */
#define XC_HYB_MGGA_X_REVM11    304 /* revM11 hybrid exchange functional from Minnesota     */

typedef struct{
  const double a[12], b[12];
} mgga_x_m11_params;

#define N_PAR 27
static const char *names[N_PAR] = {"_a0", "_a1", "_a2", "_a3", "_a4", "_a5", "_a6", "_a7", "_a8", "_a9", "_a10", "_a11", "_b0", "_b1", "_b2", "_b3", "_b4", "_b5", "_b6", "_b7", "_b8", "_b9", "_b10", "_b11", "_alpha", "_beta", "_omega"};
static const char *desc[N_PAR] = {"a0 parameter", "a1 parameter", "a2 parameter", "a3 parameter", "a4 parameter", "a5 parameter", "a6 parameter", "a7 parameter", "a8 parameter", "a9 parameter", "a10 parameter", "a11 parameter", "b0 parameter", "b1 parameter", "b2 parameter", "b3 parameter", "b4 parameter", "b5 parameter", "b6 parameter", "b7 parameter", "b8 parameter", "b9 parameter", "b10 parameter", "b11 parameter", "exact exchange", "short-range exchange", "range-separation"};

static const double par_m11[N_PAR] = {
    -0.18399900e+00, -1.39046703e+01,  1.18206837e+01,  3.10098465e+01, -5.19625696e+01,  1.55750312e+01,
    -6.94775730e+00, -1.58465014e+02, -1.48447565e+00,  5.51042124e+01, -1.34714184e+01,  0.00000000e+00,
     0.75599900e+00,  1.37137944e+01, -1.27998304e+01, -2.93428814e+01,  5.91075674e+01, -2.27604866e+01,
    -1.02769340e+01,  1.64752731e+02,  1.85349258e+01, -5.56825639e+01,  7.47980859e+00,  0.00000000e+00,
    1.0, -(1.0 - 0.428), 0.25
};

static const double par_revm11[N_PAR] = {
   -0.3288860885e+00, -8.3888150476e+00,  0.7123891057e+00,  3.6196212952e+00,  4.3941708207e+00,  5.0453345584e+00,
    7.8667061191e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,
    1.1038860885e+00,  8.0476369587e+00, -0.7353624773e+00, -2.4735275550e+00, -4.7319060355e+00, -5.8502502096e+00,
   -7.5059975327e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,
    1.0, -(1.0 - 0.225), 0.40
};

static void
mgga_x_m11_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_m11_params));

  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_m11.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_m11 = {
  XC_HYB_MGGA_X_M11,
  XC_EXCHANGE,
  "Minnesota M11 hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Peverati2011_2810, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-11,
  {N_PAR, names, desc, par_m11, set_ext_params_cpy_cam},
  mgga_x_m11_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_revm11 = {
  XC_HYB_MGGA_X_REVM11,
  XC_EXCHANGE,
  "Revised Minnesota M11 hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Verma2019_2966, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-11,
  {N_PAR, names, desc, par_revm11, set_ext_params_cpy_cam},
  mgga_x_m11_init, NULL,
  NULL, NULL, work_mgga,
};
