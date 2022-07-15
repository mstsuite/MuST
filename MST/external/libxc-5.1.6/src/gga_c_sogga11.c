/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_SOGGA11   152 /* SOGGA11 correlation */
#define XC_GGA_C_SOGGA11_X 159 /* SOGGA11-X correlation  */

typedef struct {
  double sogga11_a[6], sogga11_b[6];
} gga_c_sogga11_params;

#define N_PAR 12
static const char *names[N_PAR] = {"_a0", "_a1", "_a2", "_a3", "_a4", "_a5",
                                   "_b0", "_b1", "_b2", "_b3", "_b4", "_b5"};
static const char *desc[N_PAR] = {"a0", "a1", "a2", "a3", "a4", "a5",
                                  "b0", "b1", "b2", "b3", "b4", "b5"};

static const double par_sogga11[N_PAR] = {
    0.50000, -4.62334, 8.00410, -130.226, 38.2685,  69.5599,
    0.50000, 3.62334,  9.36393, 34.5114,  -18.5684, -0.16519};

static const double par_sogga11_x[N_PAR] = {
    0.50000, 78.2439,  25.7211, -13.8830, -9.87375, -14.1357,
    0.50000, -79.2439, 16.3725, 2.08129,  7.50769,  -10.1861};

static void gga_c_sogga11_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_sogga11_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_sogga11.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_sogga11 = {
  XC_GGA_C_SOGGA11,
  XC_CORRELATION,
  "Second-order generalized gradient approximation 2011",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2011_1991, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_sogga11, set_ext_params_cpy},
  gga_c_sogga11_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_sogga11_x = {
  XC_GGA_C_SOGGA11_X,
  XC_CORRELATION,
  "To be used with HYB_GGA_X_SOGGA11_X",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2011_191102, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_sogga11_x, set_ext_params_cpy},
  gga_c_sogga11_init, NULL,
  NULL, work_gga, NULL
};
