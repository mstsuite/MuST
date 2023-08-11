/*
 Copyright (C) 2019 M.A.L. Marques
                    Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_TASK      707 /* TASK exchange of Aschebrock and Kuemmel */
#define XC_MGGA_X_MTASK     724 /* modified TASK exchange */

typedef struct{
  double task_c, task_d, task_h0x;
  double task_anu[3], task_bnu[5];
} mgga_x_task_params;


static void
mgga_x_task_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_task_params));
}

#define TASK_N_PAR 11
static const char  *task_names[TASK_N_PAR]  = {
  "_c", "_d", "_h0x", "_anu0", "_anu1", "_anu2", "_bnu0", "_bnu1", "_bnu2", "_bnu3", "_bnu4"
};
static const char  *task_desc[TASK_N_PAR]   = {
  "Value of the constant in the exponent of g_x",
  "Value of the exponent of g_x(s^2)^c",
  "Value of h_x^0",
  "Coefficient 0 of the Chebyshev expansion for h_x^1",
  "Coefficient 1 of the Chebyshev expansion for h_x^1",
  "Coefficient 2 of the Chebyshev expansion for h_x^1",
  "Coefficient 0 of the Chebyshev expansion for fx(a)",
  "Coefficient 1 of the Chebyshev expansion for fx(a)",
  "Coefficient 2 of the Chebyshev expansion for fx(a)",
  "Coefficient 3 of the Chebyshev expansion for fx(a)",
  "Coefficient 4 of the Chebyshev expansion for fx(a)"
};

static const double task_values[TASK_N_PAR] = {
  4.9479, 10.0, 1.174, 0.938719, -0.076371, -0.0150899, -0.628591, -2.10315, -0.5, 0.103153, 0.128591
};
static const double mtask_values[TASK_N_PAR] = {
  4.9479, 10.0, 1.29, 0.924374, -0.09276847, -0.017143, -0.639572, -2.087488, -0.625, -0.162512, 0.014572
};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_task.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_task = {
  XC_MGGA_X_TASK,
  XC_EXCHANGE,
  "TASK exchange of Aschebrock and Kuemmel",
  XC_FAMILY_MGGA,
  {&xc_ref_Aschebrock2019_033082, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TASK_N_PAR, task_names, task_desc, task_values, set_ext_params_cpy},
  mgga_x_task_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_mtask = {
  XC_MGGA_X_MTASK,
  XC_EXCHANGE,
  "modified TASK exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Neupane2021_, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TASK_N_PAR, task_names, task_desc, mtask_values, set_ext_params_cpy},
  mgga_x_task_init, NULL,
  NULL, NULL, work_mgga
};
