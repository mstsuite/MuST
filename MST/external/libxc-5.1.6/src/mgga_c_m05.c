/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_M05           237 /* M05 correlation functional from Minnesota     */
#define XC_MGGA_C_M05_2X        238 /* M05-2X correlation functional from Minnesota  */
#define XC_MGGA_C_DLDF           37 /* Dispersionless Density Functional             */

typedef struct{
  double gamma_ss, gamma_ab;
  const double css[5], cab[5];
  double Fermi_D_cnst; /* correction term similar to 10.1063/1.2800011 */
} mgga_c_m05_params;

static void
mgga_c_vsxc_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_m05_params));
}

#define M05_N_PAR 13
static const char  *m05_names[M05_N_PAR]  = {
  "_gamma_ss", "_gamma_ab",
  "_css0", "_css1", "_css2", "_css3", "_css4",
  "_cab0", "_cab1", "_cab2", "_cab3", "_cab4",
  "_Fermi_D_cnst"
};
static const char  *m05_desc[M05_N_PAR]   = {
  "gamma_ss", "gamma_ab",
  "css0", "css1", "css2", "css3", "css4",
  "cab0", "cab1", "cab2", "cab3", "cab4",
  "Constant for the correction term similar to 10.1063/1.2800011"
};
static const double m05_values[M05_N_PAR] = {
  0.06, 0.0031,
  1.00000e0,  3.77344e0, -26.04463e0, 30.69913e0, -9.22695e0,
  1.00000e0,  3.78569e0, -14.15261e0, -7.46589e0, 17.94491e0 ,
  1e-10
};
static const double m05_2x_values[M05_N_PAR] = {
  0.06, 0.0031,
  1.00000e0, -3.05430e0,  7.61854e0,  1.47665e0, -11.92365e0,
  1.00000e0,  1.09297e0, -3.79171e0,  2.82810e0, -10.58909e0,
  1e-10
};
static const double m05_dldf_values[M05_N_PAR] = {
  0.06, 0.0031,
  1.00000e0, -2.5960897,   2.2233793, 0.0, 0.0,
  1.00000e0,  5.9515308, -11.1602877, 0.0, 0.0,
  1e-10
};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_m05.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_m05 = {
  XC_MGGA_C_M05,
  XC_CORRELATION,
  "Minnesota M05 correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2005_161103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-15,
  {M05_N_PAR, m05_names, m05_desc, m05_values, set_ext_params_cpy},
  mgga_c_vsxc_init, NULL,
  NULL, NULL, work_mgga
};


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_m05_2x = {
  XC_MGGA_C_M05_2X,
  XC_CORRELATION,
  "Minnesota M05-2X correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2006_364, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-15,
  {M05_N_PAR, m05_names, m05_desc, m05_2x_values, set_ext_params_cpy},
  mgga_c_vsxc_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_dldf = {
  XC_MGGA_C_DLDF,
  XC_CORRELATION,
  "Dispersionless Density Functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Pernal2009_263201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  5.0e-15,
  {M05_N_PAR, m05_names, m05_desc, m05_dldf_values, set_ext_params_cpy},
  mgga_c_vsxc_init, NULL,
  NULL, NULL, work_mgga
};
