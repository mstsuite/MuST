/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_M06_L         233 /* M06-L correlation functional from Minnesota          */
#define XC_MGGA_C_M06_HF        234 /* M06-HF correlation functional from Minnesota         */
#define XC_MGGA_C_M06           235 /* M06 correlation functional from Minnesota            */
#define XC_MGGA_C_M06_2X        236 /* M06-2X correlation functional from Minnesota         */
#define XC_MGGA_C_REVM06_L      294 /* Revised M06-L correlation functional from Minnesota  */
#define XC_MGGA_C_REVM06        306 /* Revised M06 correlation functional from Minnesota    */
#define XC_MGGA_C_M06_SX        311 /* M06-SX correlation functional from Minnesota    */

typedef struct{
  double gamma_ss, gamma_ab, alpha_ss, alpha_ab;
  const double css[5], cab[5], dss[6], dab[6];
  double Fermi_D_cnst; /* correction term similar to 10.1063/1.2800011 */
} mgga_c_m06l_params;

static void
mgga_c_m06l_init(xc_func_type *p)
{
  assert(p != NULL);

  p->n_func_aux  = 1;
  p->func_aux    = (xc_func_type **) libxc_malloc(1*sizeof(xc_func_type *));
  p->func_aux[0] = (xc_func_type *)  libxc_malloc(  sizeof(xc_func_type));

  xc_func_init(p->func_aux[0], XC_LDA_C_PW_MOD, XC_POLARIZED);

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_m06l_params));
}

#define M06L_N_PAR 27
static const char  *m06l_names[M06L_N_PAR]  = {
  "_gamma_ss", "_gamma_ab", "_alpha_ss", "_alpha_ab",
  "_css0", "_css1", "_css2", "_css3", "_css4",
  "_cab0", "_cab1", "_cab2", "_cab3", "_cab4",
  "_dss0", "_dss1", "_dss2", "_dss3", "_dss4", "_dss5",
  "_dab0", "_dab1", "_dab2", "_dab3", "_dab4", "_dab5",
  "_Fermi_D_cnst"
};
static const char  *m06l_desc[M06L_N_PAR]   = {
  "gamma_ss", "gamma_ab", "alpha_ss", "alpha_ab",
  "css0", "css1", "css2", "css3", "css4",
  "cab0", "cab1", "cab2", "cab3", "cab4",
  "dss0", "dss1", "dss2", "dss3", "dss4", "dss5",
  "dab0", "dab1", "dab2", "dab3", "dab4", "dab5",
  "Constant for the correction term similar to 10.1063/1.2800011"
};
static const double m06l_values[M06L_N_PAR] = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  5.349466e-01,  5.396620e-01, -3.161217e+01,  5.149592e+01, -2.919613e+01,
  6.042374e-01,  1.776783e+02, -2.513252e+02,  7.635173e+01, -1.255699e+01,
  4.650534e-01,  1.617589e-01,  1.833657e-01,  4.692100e-04, -4.990573e-03,  0.000000e+00,
  3.957626e-01, -5.614546e-01,  1.403963e-02,  9.831442e-04, -3.577176e-03,  0.000000e+00,
  1e-10
};
static const double m06hf_values[M06L_N_PAR] = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  1.023254e-01, -2.453783e+00,  2.913180e+01, -3.494358e+01,  2.315955e+01,
  1.674634e+00,  5.732017e+01,  5.955416e+01, -2.311007e+02,  1.255199e+02,
  8.976746e-01, -2.345830e-01,  2.368173e-01, -9.913890e-04, -1.146165e-02,  0.000000e+00,
  -6.746338e-01, -1.534002e-01, -9.021521e-02, -1.292037e-03, -2.352983e-04,  0.000000e+00,
  1e-10
};
static const double m06_values[M06L_N_PAR] = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  5.094055e-01, -1.491085e+00,  1.723922e+01, -3.859018e+01,  2.845044e+01,
  3.741539e+00,  2.187098e+02, -4.531252e+02,  2.936479e+02, -6.287470e+01,
  4.905945e-01, -1.437348e-01,  2.357824e-01,  1.871015e-03, -3.788963e-03,  0.000000e+00,
  -2.741539e+00, -6.720113e-01, -7.932688e-02,  1.918681e-03, -2.032902e-03,  0.000000e+00,
  1e-10
};
static const double m062x_values[M06L_N_PAR] = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  3.097855e-01, -5.528642e+00,  1.347420e+01, -3.213623e+01,  2.846742e+01,
  8.833596e-01,  3.357972e+01, -7.043548e+01,  4.978271e+01, -1.852891e+01,
  6.902145e-01,  9.847204e-02,  2.214797e-01, -1.968264e-03, -6.775479e-03,  0.000000e+00,
  1.166404e-01, -9.120847e-02, -6.726189e-02,  6.720580e-05,  8.448011e-04,  0.000000e+00,
  1e-10
};
static const double revm06l_values[M06L_N_PAR] = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  1.227659748,  0.855201283, -3.113346677, -2.239678026,  0.354638962,
  0.344360696, -0.557080242, -2.009821162, -1.857641887, -1.076639864,
  -0.538821292, -0.028296030, 0.023889696,   0.0, 0.0,   -0.002437902,
  0.400714600,  0.015796569, -0.032680984,   0.0, 0.0,    0.001260132,
  1e-10
};
static const double revm06_values[M06L_N_PAR] = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  0.9017224575,  0.2079991827, -1.823747562, -1.384430429, -0.4423253381,
  1.222401598,  0.6613907336, -1.884581043, -2.780360568, -3.068579344,
  -0.1467095900, -0.0001832187007, 0.08484372430, 0.0, 0.0,  0.0002280677172,
  -0.3390666720,  0.003790156384, -0.02762485975, 0.0, 0.0,  0.0004076285162,
  1e-10
};
static const double m06sx_values[M06L_N_PAR] = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  1.17575011057022E+00,  6.58083496678423E-01, -2.78913774852905E+00, -1.18597601856255E+00,  1.16439928209688E+00,
  1.63738167314691E-01, -4.36481171027951E-01, -1.90232628449712E+00, -1.42432902881841E+00, -9.05909137360893E-01,
  8.17322574473352E-02, -2.88531085759385E-02,  9.05917734868130E-02, 0.0, 0.0, -4.86297499082106E-04,
  7.40594619832397E-01,  1.23306511345974E-02, -1.88253421850249E-02, 0.0, 0.0,  4.87276242162303E-04,
  1e-10
};

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_m06l.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_m06_l = {
  XC_MGGA_C_M06_L,
  XC_CORRELATION,
  "Minnesota M06-L correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2006_194101, &xc_ref_Zhao2008_215, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-12,
  {M06L_N_PAR, m06l_names, m06l_desc, m06l_values, set_ext_params_cpy},
  mgga_c_m06l_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_m06_hf = {
  XC_MGGA_C_M06_HF,
  XC_CORRELATION,
  "Minnesota M06-HF correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2006_13126, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-12,
  {M06L_N_PAR, m06l_names, m06l_desc, m06hf_values, set_ext_params_cpy},
  mgga_c_m06l_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_m06 = {
  XC_MGGA_C_M06,
  XC_CORRELATION,
  "Minnesota M06 correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-12,
  {M06L_N_PAR, m06l_names, m06l_desc, m06_values, set_ext_params_cpy},
  mgga_c_m06l_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_m06_2x = {
  XC_MGGA_C_M06_2X,
  XC_CORRELATION,
  "Minnesota M06-2X correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-12,
  {M06L_N_PAR, m06l_names, m06l_desc, m062x_values, set_ext_params_cpy},
  mgga_c_m06l_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revm06_l = {
  XC_MGGA_C_REVM06_L,
  XC_CORRELATION,
  "Minnesota revM06-L correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Wang2017_8487, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-12,
  {M06L_N_PAR, m06l_names, m06l_desc, revm06l_values, set_ext_params_cpy},
  mgga_c_m06l_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revm06 = {
  XC_MGGA_C_REVM06,
  XC_CORRELATION,
  "Revised Minnesota M06 correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Wang2018_10257, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-12,
  {M06L_N_PAR, m06l_names, m06l_desc, revm06_values, set_ext_params_cpy},
  mgga_c_m06l_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_m06_sx = {
  XC_MGGA_C_M06_SX,
  XC_CORRELATION,
  "Minnesota M06-SX correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Wang2020_2294, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-14,
  {M06L_N_PAR, m06l_names, m06l_desc, m06sx_values, set_ext_params_cpy},
  mgga_c_m06l_init, NULL,
  NULL, NULL, work_mgga,
};
