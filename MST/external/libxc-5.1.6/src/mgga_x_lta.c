/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_LTA          201 /* Local tau approximation of Ernzerhof & Scuseria */
#define XC_MGGA_X_TLDA         685 /* LDA-type exchange with tau-dependent potential */
#define XC_MGGA_X_HLTA         698 /* Half-and-half by Lehtola and Marques */

typedef struct{
  double ltafrac;
} mgga_x_lta_params;

static void
mgga_x_lta_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_lta_params));
}

#define LTA_N_PAR 1
static const char  *lta_names[LTA_N_PAR]   = {"_ltafrac"};
static const char  *lta_desc[LTA_N_PAR]    = {"Fraction of LTA density"};
static const double lta_values[LTA_N_PAR]  = {1.0};
static const double tlda_values[LTA_N_PAR] = {0.25};
static const double hlta_values[LTA_N_PAR] = {0.5};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_lta.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_lta = {
  XC_MGGA_X_LTA,
  XC_EXCHANGE,
  "Local tau approximation",
  XC_FAMILY_MGGA,
  {&xc_ref_Ernzerhof1999_911, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {LTA_N_PAR, lta_names, lta_desc, lta_values, set_ext_params_cpy},
  mgga_x_lta_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_tlda = {
  XC_MGGA_X_TLDA,
  XC_EXCHANGE,
  "LDA-type exchange with tau-dependent potential",
  XC_FAMILY_MGGA,
  {&xc_ref_Eich2014_224107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {LTA_N_PAR, lta_names, lta_desc, tlda_values, set_ext_params_cpy},
  mgga_x_lta_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_hlta = {
  XC_MGGA_X_HLTA,
  XC_EXCHANGE,
  "Half-and-half meta-LDAized LDA exchange by Lehtola and Marques",
  XC_FAMILY_MGGA,
  {&xc_ref_Lehtola2021_943, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {LTA_N_PAR, lta_names, lta_desc, hlta_values, set_ext_params_cpy},
  mgga_x_lta_init, NULL,
  NULL, NULL, work_mgga,
};
