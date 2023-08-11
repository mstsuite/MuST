/*
 Copyright (C) 2008 Lara Ferrigni, Georg Madsen, M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_M06_L         203 /* M06-L exchange functional from Minnesota          */
#define XC_HYB_MGGA_X_M06_HF    444 /* M06-HF hybrid exchange functional from Minnesota  */
#define XC_HYB_MGGA_X_M06       449 /* M06 hybrid exchange functional from Minnesota     */
#define XC_MGGA_X_REVM06_L      293 /* revised M06-L exchange functional from Minnesota  */
#define XC_HYB_MGGA_X_REVM06    305 /* revised M06 hybrid exchange functional from Minnesota  */
#define XC_HYB_MGGA_X_M06_SX    310 /* M06-SX exchange functional from Minnesota  */

typedef struct {
  double a[12], d[6];
} mgga_x_m06l_params;

#define N_PAR_PURE 18
static const char  *pure_names[N_PAR_PURE]  = {
  "_a0",  "_a1",  "_a2",  "_a3",  "_a4",
  "_a5", "_a6", "_a7", "_a8", "_a9",
  "_a10", "_a11", "_d0", "_d1", "_d2",
  "_d3", "_d4", "_d5"
};
static const char  *pure_desc[N_PAR_PURE]   = {
  "_a0 parameter",
  "_a1 parameter",
  "_a2 parameter",
  "_a3 parameter",
  "_a4 parameter",
  "_a5 parameter",
  "_a6 parameter",
  "_a7 parameter",
  "_a8 parameter",
  "_a9 parameter",
  "_a10 parameter",
  "_a11 parameter",
  "_d0 parameter",
  "_d1 parameter",
  "_d2 parameter",
  "_d3 parameter",
  "_d4 parameter",
  "_d5 parameter"
};

#define N_PAR_HYB 19
static const char  *hyb_names[N_PAR_HYB]  = {
  "_a0",  "_a1",  "_a2",  "_a3",  "_a4",
  "_a5", "_a6", "_a7", "_a8", "_a9",
  "_a10", "_a11", "_d0", "_d1", "_d2",
  "_d3", "_d4", "_d5", "_X"
};
static const char  *hyb_desc[N_PAR_HYB]   = {
  "_a0 parameter",
  "_a1 parameter",
  "_a2 parameter",
  "_a3 parameter",
  "_a4 parameter",
  "_a5 parameter",
  "_a6 parameter",
  "_a7 parameter",
  "_a8 parameter",
  "_a9 parameter",
  "_a10 parameter",
  "_a11 parameter",
  "_d0 parameter",
  "_d1 parameter",
  "_d2 parameter",
  "_d3 parameter",
  "_d4 parameter",
  "_d5 parameter",
  "Fraction of exact exchange"
};

#define N_PAR_SR 20
static const char  *sr_names[N_PAR_SR]  = {
  "_a0",  "_a1",  "_a2",  "_a3",  "_a4",
  "_a5", "_a6", "_a7", "_a8", "_a9",
  "_a10", "_a11", "_b0", "_b1", "_b2",
  "b3", "_b4", "_b5", "_beta", "_omega"
};
static const char  *sr_desc[N_PAR_SR]   = {
  "_a0 parameter",
  "_a1 parameter",
  "_a2 parameter",
  "_a3 parameter",
  "_a4 parameter",
  "_a5 parameter",
  "_a6 parameter",
  "_a7 parameter",
  "_a8 parameter",
  "_a9 parameter",
  "_a10 parameter",
  "_a11 parameter",
  "_d0 parameter",
  "_d1 parameter",
  "_d2 parameter",
  "_d3 parameter",
  "_d4 parameter",
  "_d5 parameter",
  "Fraction of short-range exchange",
  "Range separation parameter"
};

static const double par_m06l[N_PAR_PURE] = {
  0.3987756, 0.2548219, 0.3923994, -2.103655, -6.302147, 10.97615,
   30.97273,  -23.18489, -56.73480, 21.60364, 34.21814, -9.049762,
  0.6012244, 0.004748822, -0.008635108, -0.000009308062, 0.00004482811, 0.0};

static const double par_m06hf[N_PAR_HYB] = {
   1.179732e-01, -1.066708e+00, -1.462405e-01,  7.481848e+00,  3.776679e+00, -4.436118e+01,
  -1.830962e+01,  1.003903e+02,  3.864360e+01, -9.806018e+01, -2.557716e+01,  3.590404e+01,
  -1.179732e-01, -2.500000e-03, -1.180065e-02, 0.0, 0.0, 0.0, 1.0};

static const double par_m06[N_PAR_HYB] = {
   5.877943e-01, -1.371776e-01,  2.682367e-01, -2.515898e+00, -2.978892e+00,  8.710679e+00,
   1.688195e+01, -4.489724e+00, -3.299983e+01, -1.449050e+01,  2.043747e+01,  1.256504e+01,
   1.422057e-01, 7.370319e-04, -1.601373e-02, 0.0, 0.0, 0.0, 0.27};

static const double par_revm06l[N_PAR_PURE] = {
  1.423227252,  0.471820438, -0.167555701, -0.250154262,  0.062487588,  0.733501240,
 -2.359736776, -1.436594372,  0.444643793,  1.529925054,  2.053941717, -0.036536031,
 -0.423227252, 0.0, 0.003724234, 0.0, 0.0, 0.0};

static const double par_revm06[N_PAR_HYB] = {
   0.6511394014, -0.1214497763, -0.1367041135,  0.3987218551,  0.6056741356, -2.379738662,
  -1.492098351,   3.031473420,   0.5149637108,  2.633751911,  0.9886749252, -4.243714128,
  -0.05523940140, 0.0, -0.003782631233, 0.0, 0.0, 0.0, 0.4041};

static const double par_m06_sx[N_PAR_SR] = {
  9.96501680264007E-01, 3.01264933631367E-02, -1.03366758333673E-01,
 -1.55653062500239E-01, 7.95768051149902E-03,  8.71986277454856E-02,
 -8.16152625764469E-01, 6.72773006612420E-01,  5.21127186174968E-01,
  3.99466945122217E-01, 5.19400018999204E-01, -9.65261552636835E-01,
 -3.47792307472902E-01, 0.0, -2.70366787478266E-03, 0.0, 0.0, 0.0, 0.335, 0.10};

static void
mgga_x_m06l_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_m06l_params));

  switch(p->info->number) {
  case(XC_HYB_MGGA_X_M06):
  case(XC_HYB_MGGA_X_M06_HF):
  case(XC_HYB_MGGA_X_REVM06):
    xc_hyb_init_hybrid(p, 0.0);
    break;

  case(XC_HYB_MGGA_X_M06_SX):
    xc_hyb_init_sr(p, 0.0, 0.0);
    break;
  }
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_m06l.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_m06_l = {
  XC_MGGA_X_M06_L,
  XC_EXCHANGE,
  "Minnesota M06-L exchange functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2006_194101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_PURE, pure_names, pure_desc, par_m06l, set_ext_params_cpy},
  mgga_x_m06l_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_m06_hf = {
  XC_HYB_MGGA_X_M06_HF,
  XC_EXCHANGE,
  "Minnesota M06-HF hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2006_13126, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_HYB, hyb_names, hyb_desc, par_m06hf, set_ext_params_cpy_exx},
  mgga_x_m06l_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_m06 = {
  XC_HYB_MGGA_X_M06,
  XC_EXCHANGE,
  "Minnesota M06 hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_HYB, hyb_names, hyb_desc, par_m06, set_ext_params_cpy_exx},
  mgga_x_m06l_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_revm06_l = {
  XC_MGGA_X_REVM06_L,
  XC_EXCHANGE,
  "Minnesota revM06-L exchange functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Wang2017_8487, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_PURE, pure_names, pure_desc, par_revm06l, set_ext_params_cpy},
  mgga_x_m06l_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_revm06 = {
  XC_HYB_MGGA_X_REVM06,
  XC_EXCHANGE,
  "Revised Minnesota M06 hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Wang2018_10257, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_HYB, hyb_names, hyb_desc, par_revm06, set_ext_params_cpy_exx},
  mgga_x_m06l_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_m06_sx = {
  XC_HYB_MGGA_X_M06_SX,
  XC_EXCHANGE,
  "Minnesota M06-SX short-range hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Wang2020_2294, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_SR, sr_names, sr_desc, par_m06_sx, set_ext_params_cpy_cam_sr},
  mgga_x_m06l_init, NULL,
  NULL, NULL, work_mgga
};
