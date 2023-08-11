/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_VSXC          232 /* VSxc from Van Voorhis and Scuseria (correlation part) */

typedef struct{
  const double alpha_ss, alpha_ab;
  const double dss[6], dab[6];
} mgga_c_vsxc_params;

#define N_PAR 14
static const char *names[N_PAR] = {"_alpha_ss", "_alpha_os", "_dss0", "_dss1", "_dss2", "_dss3", "_dss4", "_dss5", "_dab0", "_dab1", "_dab2", "_dab3", "_dab4", "_dab5"};
static const char *desc[N_PAR] = {"same-spin alpha", "opposite-spin alpha", "same-spin a parameter", "same-spin b parameter", "same-spin c parameter", "same-spin d parameter", "same-spin e parameter", "same-spin f parameter", "opposite-spin a parameter", "opposite-spin b parameter", "opposite-spin c parameter", "opposite-spin d parameter", "opposite-spin e parameter", "opposite-spin f parameter"};

static const double par_vsxc[N_PAR] = {
  0.00515088, 0.00304966,
  3.270912e-01, -3.228915e-02, -2.942406e-02,  2.134222e-03, -5.451559e-03,  1.577575e-02,
  7.035010e-01,  7.694574e-03,  5.152765e-02,  3.394308e-05, -1.269420e-03,  1.296118e-03
};

static void
mgga_c_vsxc_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_vsxc_params));
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_vsxc.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_vsxc = {
  XC_MGGA_C_VSXC,
  XC_CORRELATION,
  "VSXC (correlation part)",
  XC_FAMILY_MGGA,
  {&xc_ref_VanVoorhis1998_400, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_vsxc, set_ext_params_cpy},
  mgga_c_vsxc_init, NULL,
  NULL, NULL, work_mgga,
};
