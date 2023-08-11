/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_PBE2         616 /* Two   parameter PBE-like expansion             */
#define XC_GGA_K_PBE3         595 /* Three parameter PBE-like expansion             */
#define XC_GGA_K_PBE4         596 /* Four  parameter PBE-like expansion             */

typedef struct{
  double a;
  double c1, c2, c3;
} gga_k_mpbe_params;

#define N_PAR 4
static const char  *names[N_PAR]  = {"_a", "_c1", "_c2", "_c3"};
static const char  *desc[N_PAR]   = {"a", "c1", "c2", "c3"};
static const double kpbe2_val[N_PAR] = {0.2942, 2.0309, 0.0, 0.0};
static const double kpbe3_val[N_PAR] = {4.1355, -3.7425, 50.258, 0.0};
static const double kpbe4_val[N_PAR] = {1.7107, -7.2333, 61.645, -93.683};

static void
gga_k_mpbe_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_mpbe_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_mpbe.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_pbe2 = {
  XC_GGA_K_PBE2,
  XC_KINETIC,
  "Three parameter PBE-like expansion",
  XC_FAMILY_GGA,
  {&xc_ref_Karasiev2006_111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, kpbe2_val, set_ext_params_cpy},
  gga_k_mpbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_pbe3 = {
  XC_GGA_K_PBE3,
  XC_KINETIC,
  "Three parameter PBE-like expansion",
  XC_FAMILY_GGA,
  {&xc_ref_Karasiev2006_111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, kpbe3_val, set_ext_params_cpy},
  gga_k_mpbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_pbe4 = {
  XC_GGA_K_PBE4,
  XC_KINETIC,
  "Four parameter PBE-like expansion",
  XC_FAMILY_GGA,
  {&xc_ref_Karasiev2006_111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, kpbe4_val, set_ext_params_cpy},
  gga_k_mpbe_init, NULL,
  NULL, work_gga, NULL
};
