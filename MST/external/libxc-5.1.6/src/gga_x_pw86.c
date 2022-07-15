/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_PW86  108 /* Perdew & Wang 86 */
#define XC_GGA_X_RPW86 144 /* refitted Perdew & Wang 86 */

typedef struct {
  double aa, bb, cc;
} gga_x_pw86_params;

#define N_PAR 3
static const char *names[N_PAR] = {"_aa", "_bb", "_cc"};
static const char *desc[N_PAR] = {"Coefficient of s^2 term",
                                  "Coefficient of s^4 term",
                                  "Coefficient of s^6 term"};

static const double par_pw86[3] = {1.296, 14.0, 0.2};

static const double par_rpw86[3] = {15.0 * 0.1234, 17.33, 0.163};

static void gga_x_pw86_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pw86_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_pw86.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pw86 = {
  XC_GGA_X_PW86,
  XC_EXCHANGE,
  "Perdew & Wang 86",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1986_8800, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_pw86, set_ext_params_cpy},
  gga_x_pw86_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_rpw86 = {
  XC_GGA_X_RPW86,
  XC_EXCHANGE,
  "Refitted Perdew & Wang 86",
  XC_FAMILY_GGA,
  {&xc_ref_Murray2009_2754, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_rpw86, set_ext_params_cpy},
  gga_x_pw86_init, NULL,
  NULL, work_gga, NULL
};
