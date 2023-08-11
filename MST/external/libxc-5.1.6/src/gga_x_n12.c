/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_N12          82 /* N12 functional from Minnesota    */
#define XC_HYB_GGA_X_N12_SX   81 /* N12-SX functional from Minnesota */
#define XC_GGA_X_GAM          32 /* GAM functional from Minnesota */

typedef struct{
  double CC[4][4];
} gga_x_n12_params;

#define N_PAR_PURE 16
static const char  *pure_names[N_PAR_PURE] = {"_CC00", "_CC01", "_CC02",
  "_CC03", "_CC10", "_CC11", "_CC12", "_CC13", "_CC20", "_CC21", "_CC22",
  "_CC23", "_CC30", "_CC31", "_CC32", "_CC33"};
static const char  *pure_desc[N_PAR_PURE]  = {"_CC00", "_CC01", "_CC02",
  "_CC03", "_CC10", "_CC11", "_CC12", "_CC13", "_CC20", "_CC21", "_CC22",
  "_CC23", "_CC30", "_CC31", "_CC32", "_CC33"};

#define N_PAR_SX 18
static const char  *sx_names[N_PAR_SX] = {"_CC00", "_CC01", "_CC02",
  "_CC03", "_CC10", "_CC11", "_CC12", "_CC13", "_CC20", "_CC21", "_CC22",
  "_CC23", "_CC30", "_CC31", "_CC32", "_CC33", "_beta", "_omega"};
static const char  *sx_desc[N_PAR_SX]  = {"_CC00", "_CC01", "_CC02",
  "_CC03", "_CC10", "_CC11", "_CC12", "_CC13", "_CC20", "_CC21", "_CC22",
  "_CC23", "_CC30", "_CC31", "_CC32", "_CC33",
  "Fraction of short-range exact exchange", "Range separation parameter"};

static const double par_n12[N_PAR_PURE] = {
   1.00000e+00,  5.07880e-01,  1.68233e-01,  1.28887e-01,
   8.60211e-02, -1.71008e+01,  6.50814e+01, -7.01726e+01,
  -3.90755e-01,  5.13392e+01, -1.66220e+02,  1.42738e+02,
   4.03611e-01, -3.44631e+01,  7.61661e+01, -2.41834e+00
};

static const double par_n12_sx[N_PAR_SX] = {
  /* Indices are wrong in the original paper; the first two indices
     need to be flipped */
   6.81116e-01,  1.88858e+00,  1.78590e+00,  8.79456e-01,
  -8.12270e-02, -1.08723e+00, -4.18682e+00, -3.00000e+01,
   5.36236e-01, -5.45678e+00,  3.00000e+01,  5.51105e+01,
  -7.09913e-01,  1.30001e+01, -7.24877e+01,  2.98363e+01,
   0.25, 0.11
};

static const double par_gam[N_PAR_PURE] = {
   1.32730,    0.886102, -5.73833,   8.60197,
  -0.786018,  -4.78787,   3.90989,  -2.11611,
   0.802575,  14.4363,    8.42735,  -6.21552,
  -0.142331, -13.4598,    1.52355, -10.0530
};

static void
gga_x_n12_init(xc_func_type *p)
{
  assert(p != NULL);
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_n12_params));

  if(p->info->number == XC_HYB_GGA_X_N12_SX)
    xc_hyb_init_sr(p, 0.0, 0.0);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_n12.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_n12 = {
  XC_GGA_X_N12,
  XC_EXCHANGE,
  "Minnesota N12 exchange functional",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2012_2310, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR_PURE, pure_names, pure_desc, par_n12, set_ext_params_cpy},
  gga_x_n12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_n12_sx = {
  XC_HYB_GGA_X_N12_SX,
  XC_EXCHANGE,
  "Minnesota N12-SX exchange functional",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Peverati2012_16187, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR_SX, sx_names, sx_desc, par_n12_sx, set_ext_params_cpy_cam_sr},
  gga_x_n12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_gam = {
  XC_GGA_X_GAM,
  XC_EXCHANGE,
  "Minnesota GAM exhange functional",
  XC_FAMILY_GGA,
  {&xc_ref_Yu2015_12146, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {N_PAR_PURE, pure_names, pure_desc, par_gam, set_ext_params_cpy},
  gga_x_n12_init, NULL,
  NULL, work_gga, NULL
};
