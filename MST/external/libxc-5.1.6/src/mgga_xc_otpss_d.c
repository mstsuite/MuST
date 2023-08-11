/*
 Copyright (C) 2006-2013 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_XC_OTPSS_D      64  /* oTPSS_D functional of Goerigk and Grimme   */

static void
mgga_xc_otpss_d_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_TPSS, XC_MGGA_C_TPSS};
  static double funcs_coef[2] = {1.0, 1.0};

  static double par_x_tpss[] = {3.43, 0.75896, 0.165, 0.778, 0.41567, 2.0, 0.0};
  static double par_c_tpss[] = {0.08861, 0.7, 0.59, 0.9269, 0.6225, 2.1540};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_func_set_ext_params(p->func_aux[0], par_x_tpss);
  xc_func_set_ext_params(p->func_aux[1], par_c_tpss);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_xc_otpss_d = {
  XC_MGGA_XC_OTPSS_D,
  XC_EXCHANGE_CORRELATION,
  "oTPSS-D functional of Goerigk and Grimme",
  XC_FAMILY_MGGA,
  {&xc_ref_Goerigk2010_107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_xc_otpss_d_init,
  NULL, NULL, NULL, NULL
};
