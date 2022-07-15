/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_SCANL          702 /* SCAN correlation */
#define XC_MGGA_C_SCANL_RVV10    703 /* SCAN correlation + rVV10 correlation */
#define XC_MGGA_C_SCANL_VV10     704 /* SCAN correlation +  VV10 correlation */

#define N_PAR_SCANL 2
static const char *names[N_PAR_SCANL] = {
  "_a", "_b"                 /* parameters of pc07 */
};
static const char *desc[N_PAR_SCANL] = {
  "pc07 a", "pc07 b"
};
static const double par_scanl[N_PAR_SCANL] = {
  1.784720, 0.258304
};

static void
mgga_c_scanl_init(xc_func_type *p)
{
  xc_deorbitalize_init(p, XC_MGGA_C_SCAN, XC_MGGA_K_PC07_OPT);
}

static void
scanl_set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL && p->func_aux != NULL);
  xc_func_set_ext_params(p->func_aux[1], &ext_params[0]);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_scanl = {
  XC_MGGA_C_SCANL,
  XC_CORRELATION,
  "Deorbitalized SCAN (SCAN-L) correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161,&xc_ref_Sun2015_036402, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR_SCANL, names, desc, par_scanl, scanl_set_ext_params},
  mgga_c_scanl_init, NULL,
  NULL, NULL, xc_deorbitalize_func,
};


static void
mgga_c_scan_rvv10_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_C_SCANL};
  static double funcs_coef[1] = {1.0};

  xc_mix_init(p, 1, funcs_id, funcs_coef);

  p->nlc_b = 15.7;
  p->nlc_C = 0.0093;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_scanl_rvv10 = {
  XC_MGGA_C_SCANL_RVV10,
  XC_CORRELATION,
  "SCAN-L + rVV10 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Peng2016_041005, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_VV10 | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_c_scan_rvv10_init, NULL,
  NULL, NULL, NULL
};

static void
mgga_c_scan_vv10_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_C_SCANL};
  static double funcs_coef[1] = {1.0};

  xc_mix_init(p, 1, funcs_id, funcs_coef);

  p->nlc_b = 14.0;
  p->nlc_C = 0.0093;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_scanl_vv10 = {
  XC_MGGA_C_SCANL_VV10,
  XC_CORRELATION,
  "SCAN-L + VV10 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Brandenburg2016_115144, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_VV10 | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_c_scan_vv10_init, NULL,
  NULL, NULL, NULL
};
