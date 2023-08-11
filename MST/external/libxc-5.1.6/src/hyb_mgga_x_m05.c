/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_X_M05      438 /* M05 hybrid exchange functional from Minnesota     */
#define XC_HYB_MGGA_X_M05_2X   439 /* M05-2X hybrid exchange functional from Minnesota  */
#define XC_HYB_MGGA_X_M06_2X   450 /* M06-2X hybrid exchange functional from Minnesota  */

typedef struct{
  const double a[12];
  double csi_HF;
  double cx;
} mgga_x_m05_params;

#define N_PAR 14
static const char  *names[N_PAR]  = {
  "_a0",
  "_a1",
  "_a2",
  "_a3",
  "_a4",
  "_a5",
  "_a6",
  "_a7",
  "_a8",
  "_a9",
  "_a10",
  "_a11",
  "_csi_HF",
  "_cx"
};
static const char  *desc[N_PAR]   = {
  "a0 parameter",
  "a1 parameter",
  "a2 parameter",
  "a3 parameter",
  "a4 parameter",
  "a5 parameter",
  "a6 parameter",
  "a7 parameter",
  "a8 parameter",
  "a9 parameter",
  "a10 parameter",
  "a11 parameter",
  "overall scaling for DFT part",
  "fraction of exact exchange"
};

static const double par_m05[N_PAR] = {
  1.0, 0.08151, -0.43956, -3.22422, 2.01819, 8.79431, -0.00295,
  9.82029, -4.82351, -48.17574, 3.64802, 34.02248,
  1.0 - 0.28,
  0.28
};

static const double par_m05_2x[N_PAR] = {
  1.0, -0.56833, -1.30057, 5.50070, 9.06402, -32.21075, -23.73298,
  70.22996, 29.88614, -60.25778, -13.22205, 15.23694,
  1.0 - 0.56,
  0.56
};

static const double par_m06_2x[N_PAR] = {
  4.600000e-01, -2.206052e-01, -9.431788e-02,  2.164494e+00, -2.556466e+00, -1.422133e+01,
  1.555044e+01,  3.598078e+01, -2.722754e+01, -3.924093e+01,  1.522808e+01,  1.522227e+01,
  1.0, /* the mixing is already included in the params->a */
  0.54
};

static void
mgga_x_m05_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_m05_params));
  xc_hyb_init_hybrid(p, 0.0);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/hyb_mgga_x_m05.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_m05 = {
  XC_HYB_MGGA_X_M05,
  XC_EXCHANGE,
  "Minnesota M05 hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2005_161103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_m05, set_ext_params_cpy_exx},
  mgga_x_m05_init, NULL,
  NULL, NULL, work_mgga,
};


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_m05_2x = {
  XC_HYB_MGGA_X_M05_2X,
  XC_EXCHANGE,
  "Minnesota M05-2X hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2006_364, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_m05_2x, set_ext_params_cpy_exx},
  mgga_x_m05_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_m06_2x = {
  XC_HYB_MGGA_X_M06_2X,
  XC_EXCHANGE,
  "Minnesota M06-2X hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_m06_2x, set_ext_params_cpy_exx},
  mgga_x_m05_init, NULL,
  NULL, NULL, work_mgga,
};
