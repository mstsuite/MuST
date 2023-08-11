/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_DK          516 /* DePristo and Kress                    */
#define XC_GGA_K_PERDEW      517 /* Perdew                                */
#define XC_GGA_K_VSK         518 /* Vitos, Skriver, and Kollar            */
#define XC_GGA_K_VJKS        519 /* Vitos, Johansson, Kollar, and Skriver */
#define XC_GGA_K_ERNZERHOF   520 /* Ernzerhof */

typedef struct{
  double aa[5], bb[5];
} gga_k_dk_params;

#define N_PAR 10
static const char  *names[N_PAR]  = {
  "_a0", "_a1", "_a2", "_a3", "_a4", "_b0", "_b1", "_b2", "_b3", "_b4"
};
static const char  *desc[N_PAR]   = {
  "constant term in numerator",
  "coefficient for x^2 in numerator",
  "coefficient for x^4 in numerator",
  "coefficient for x^6 in numerator",
  "coefficient for x^8 in numerator",
  "constant term in denominator",
  "coefficient for x^2 in denominator",
  "coefficient for x^4 in denominator",
  "coefficient for x^6 in denominator",
  "coefficient for x^8 in denominator"
};

static void
gga_k_dk_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_dk_params));
}

#define KINS (X2S*X2S) /* conversion to s^2 */
#define KINX (5.0/27.0*KINS) /* conversion to x = (5/27 * s^2) */

/* DK is written in the x variable */
static const double par_dk[N_PAR] = {1.0, 0.95*KINX, 14.281111*KINX*KINX, -19.57962*KINX*KINX*KINX, 26.64765*KINX*KINX*KINX*KINX, 1.0, -0.05*KINX, 9.99802*KINX*KINX, 2.96805*KINX*KINX*KINX, 0.0};
/* Perdew is written in the s variable */
static const double par_perdew[N_PAR] = {1.0, 88.3960*KINS, 16.3683*KINS*KINS, 0.0, 0.0, 1.0, 88.2108*KINS, 0.0, 0.0, 0.0};
/* VSK is written in x */
static const double par_vsk[N_PAR] = {1.0, 0.95*KINX, 0.0, 9*0.396*KINX*KINX*KINX, 0.0, 1.0, -0.05*KINX, 0.396*KINX*KINX, 0.0, 0.0};
/* VJKS is written in s */
static const double par_vjks[N_PAR] = {1.0, 0.8944*KINS, 0.0, -0.0431*KINS*KINS*KINS, 0.0, 1.0, 0.6511*KINS, 0.0431*KINS*KINS, 0.0, 0.0};
/* Ernzerhof is written in s */
static const double par_ernzerhof[N_PAR] = {135.0, 28.0*KINS, 5.0*KINS*KINS, 0.0, 0.0, 135.0, 3.0*KINS, 0.0, 0.0, 0.0};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_dk.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_dk = {
  XC_GGA_K_DK,
  XC_KINETIC,
  "DePristo and Kress",
  XC_FAMILY_GGA,
  {&xc_ref_DePristo1987_438, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_dk, set_ext_params_cpy},
  gga_k_dk_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_perdew = {
  XC_GGA_K_PERDEW,
  XC_KINETIC,
  "Perdew",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1992_79, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_perdew, set_ext_params_cpy},
  gga_k_dk_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_vsk = {
  XC_GGA_K_VSK,
  XC_KINETIC,
  "Vitos, Skriver, and Kollar",
  XC_FAMILY_GGA,
  {&xc_ref_Vitos1998_12611, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_vsk, set_ext_params_cpy},
  gga_k_dk_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_vjks = {
  XC_GGA_K_VJKS,
  XC_KINETIC,
  "Vitos, Johansson, Kollar, and Skriver",
  XC_FAMILY_GGA,
  {&xc_ref_Vitos2000_052511, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_vjks, set_ext_params_cpy},
  gga_k_dk_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_ernzerhof = {
  XC_GGA_K_ERNZERHOF,
  XC_KINETIC,
  "Ernzerhof",
  XC_FAMILY_GGA,
  {&xc_ref_Ernzerhof2000_59, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR, names, desc, par_ernzerhof, set_ext_params_cpy},
  gga_k_dk_init, NULL,
  NULL, work_gga, NULL
};
