/*
 Copyright (C) 2006-2014 L. Talirz, M.A.L. Marques
                    2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_B86          103 /* Becke 86 Xalpha,beta,gamma                      */
#define XC_GGA_X_B86_MGC      105 /* Becke 86 Xalpha,beta,gamma (with mod. grad. correction) */
#define XC_GGA_X_B86_R         41 /* Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction) */
#define XC_GGA_X_OPTB86B_VDW  171 /* Becke 86 reoptimized for use with vdW functional of Dion et al */

typedef struct{
  double beta, gamma, omega;
} gga_x_b86_params;


static void
gga_x_b86_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_b86_params));
}

#define B86_N_PAR 3
static const char  *b86_names[B86_N_PAR]  = {"_beta", "_gamma", "_omega"};
static const char  *b86_desc[B86_N_PAR]   = {
  "Small x limit",
  "Parameter in the denominator",
  "Exponent of denominator"};
static const double b86_values[B86_N_PAR] =
  {0.0036/X_FACTOR_C, 0.004, 1.0};
static const double b86_mgc_values[B86_N_PAR] =
  {0.00375/X_FACTOR_C, 0.007, 4.0/5.0};
static const double b86_r_values[B86_N_PAR] =
  {MU_GE*X2S*X2S, MU_GE*X2S*X2S/0.7114, 4.0/5.0};
static const double b86_optb86b_values[B86_N_PAR] =
  {MU_GE*X2S*X2S, MU_GE*X2S*X2S, 4.0/5.0};


#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_b86.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_b86 = {
  XC_GGA_X_B86,
  XC_EXCHANGE,
  "Becke 86",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1986_4524, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {B86_N_PAR, b86_names, b86_desc, b86_values, set_ext_params_cpy},
  gga_x_b86_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_b86_mgc = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1986_4524, &xc_ref_Becke1986_7184, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {B86_N_PAR, b86_names, b86_desc, b86_mgc_values, set_ext_params_cpy},
  gga_x_b86_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_b86_r = {
  XC_GGA_X_B86_R,
  XC_EXCHANGE,
  "Revised Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  {&xc_ref_Hamada2014_121103, &xc_ref_Becke1986_4524, &xc_ref_Becke1986_7184, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {B86_N_PAR, b86_names, b86_desc, b86_r_values, set_ext_params_cpy},
  gga_x_b86_init, NULL,
  NULL, work_gga, NULL

};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_optb86b_vdw = {
  XC_GGA_X_OPTB86B_VDW,
  XC_EXCHANGE,
  "Becke 86 reoptimized for use with vdW functional of Dion et al",
  XC_FAMILY_GGA,
  {&xc_ref_Klimes2011_195131, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {B86_N_PAR, b86_names, b86_desc, b86_optb86b_values, set_ext_params_cpy},
  gga_x_b86_init, NULL,
  NULL, work_gga, NULL
};
