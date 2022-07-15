/*
 Copyright (C) 2008 Georg Madsen
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_EV93     35 /* Engel and Vosko */
#define XC_GGA_X_ECMV92  215 /* Engel, Chevary, Macdonald, and Vosko */

typedef struct{
  double a1, a2, a3;  /* numerator */
  double b1, b2, b3;  /* denominator */
} gga_x_ev93_params;

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ev93.c"
#include "work_gga.c"

static void
gga_x_ev93_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ev93_params));
}

#define EV93_N_PAR 6
static const char  *ev93_names[EV93_N_PAR]  = {"_a1", "_a2", "_a3", "_b1", "_b2", "_b3"};
static const char  *ev93_desc[EV93_N_PAR]   = {"a1", "a2", "a3", "b1", "b2", "b3"};
static const double ev93_values[EV93_N_PAR] =
  {1.647127, 0.980118, 0.017399, 1.523671, 0.367229, 0.011282};
static const double ecmv92_values[EV93_N_PAR] =
  {27.8428, 11.7683, 0.0, 27.5026, 5.7728, 0.0};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ev93 = {
  XC_GGA_X_EV93,
  XC_EXCHANGE,
  "Engel and Vosko",
  XC_FAMILY_GGA,
  {&xc_ref_Engel1993_13164, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {6, ev93_names, ev93_desc, ev93_values, set_ext_params_cpy},
  gga_x_ev93_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ecmv92 = {
  XC_GGA_X_ECMV92,
  XC_EXCHANGE,
  "Engel, Chevary, Macdonald and Vosko",
  XC_FAMILY_GGA,
  {&xc_ref_Engel1992_7, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {6, ev93_names, ev93_desc, ecmv92_values, set_ext_params_cpy},
  gga_x_ev93_init, NULL,
  NULL, work_gga, NULL
};
