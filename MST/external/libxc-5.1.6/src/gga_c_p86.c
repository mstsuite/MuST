/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_P86          132 /* Perdew 86 */
#define XC_GGA_C_P86_FT       217 /* Perdew 86 with a more accurate value for ftilde */

typedef struct{
  double malpha;
  double mbeta;
  double mgamma;
  double mdelta;
  double aa;
  double bb;
  double ftilde;
} gga_c_p86_params;

#define N_PAR 7
static const char *names[N_PAR]  = {"_malpha", "_mbeta", "_mgamma", "_mdelta", "_aa", "_bb", "_ftilde"};
static const char *desc[N_PAR]   = {"alpha in eq 6", "beta in eq 6", "gamma in eq 6", "delta in eq 6", "linear parameter in eq 6", "constant in the numerator in eq 6", "constant in eq 9"};

/* Original parmeters */
static const double p86_val[N_PAR] = {0.023266, 7.389e-6, 8.723, 0.472, 0.001667, 0.002568, 1.745*0.11};
/* The factor 1.745 is an approximation of (9*Pi)^(1/6) */
static const double p86ft_val[N_PAR] = {0.023266, 7.389e-6, 8.723, 0.472, 0.001667, 0.002568, 1.7454151061251239789*0.11};

static void gga_c_p86_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_p86_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_p86.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_p86 = {
  XC_GGA_C_P86,
  XC_CORRELATION,
  "Perdew 86",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1986_8822, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, p86_val, set_ext_params_cpy},
  gga_c_p86_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_p86_ft = {
  XC_GGA_C_P86_FT,
  XC_CORRELATION,
  "Perdew 86 with more accurate value for ftilde",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1986_8822, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, p86ft_val, set_ext_params_cpy},
  gga_c_p86_init, NULL,
  NULL, work_gga, NULL
};
