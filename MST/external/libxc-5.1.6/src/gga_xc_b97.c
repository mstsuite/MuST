/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_B97     407 /* Becke 97                                 */
#define XC_HYB_GGA_XC_B97_1   408 /* Becke 97-1                               */
#define XC_HYB_GGA_XC_B97_2   410 /* Becke 97-2                               */
#define XC_GGA_XC_B97_D       170 /* Grimme functional to be used with C6 vdW term */
#define XC_HYB_GGA_XC_B97_K   413 /* Boese-Martin for Kinetics                */
#define XC_HYB_GGA_XC_B97_3   414 /* Becke 97-3                               */
#define XC_GGA_XC_HCTH_93     161 /* HCTH functional fitted to  93 molecules  */
#define XC_GGA_XC_HCTH_120    162 /* HCTH functional fitted to 120 molecules  */
#define XC_GGA_XC_HCTH_147    163 /* HCTH functional fitted to 147 molecules  */
#define XC_GGA_XC_HCTH_407    164 /* HCTH functional fitted to 407 molecules  */
#define XC_HYB_GGA_XC_SB98_1A 420 /* Schmider-Becke 98 parameterization 1a    */
#define XC_HYB_GGA_XC_SB98_1B 421 /* Schmider-Becke 98 parameterization 1b    */
#define XC_HYB_GGA_XC_SB98_1C 422 /* Schmider-Becke 98 parameterization 1c    */
#define XC_HYB_GGA_XC_SB98_2A 423 /* Schmider-Becke 98 parameterization 2a    */
#define XC_HYB_GGA_XC_SB98_2B 424 /* Schmider-Becke 98 parameterization 2b    */
#define XC_HYB_GGA_XC_SB98_2C 425 /* Schmider-Becke 98 parameterization 2c    */
#define XC_GGA_XC_B97_GGA1     96 /* Becke 97 GGA-1                           */
#define XC_GGA_XC_HCTH_P14     95 /* HCTH p=1/4                               */
#define XC_GGA_XC_HCTH_P76     94 /* HCTH p=7/6                               */
#define XC_GGA_XC_HCTH_407P    93 /* HCTH/407+                                */
#define XC_HYB_GGA_XC_B97_1P  266 /* version of B97 by Cohen and Handy        */
#define XC_GGA_XC_HLE16       545 /* high local exchange 2016                 */

typedef struct {
  double c_x[5], c_ss[5], c_ab[5];
} gga_xc_b97_params;

#define B97_N_PAR 16
static const char  *b97_names[B97_N_PAR]  = {
  "_cx0",  "_cx1",  "_cx2",  "_cx3",  "_cx4",
  "_css0", "_css1", "_css2", "_css3", "_css4",
  "_cos0", "_cos1", "_cos2", "_cos3", "_cos4",
  "_cxx"};
static const char  *b97_desc[B97_N_PAR]   = {
  "u^0 coefficient for exchange",
  "u^1 coefficient for exchange",
  "u^2 coefficient for exchange",
  "u^3 coefficient for exchange",
  "u^4 coefficient for exchange",
  "u^0 coefficient for same-spin correlation",
  "u^1 coefficient for same-spin correlation",
  "u^2 coefficient for same-spin correlation",
  "u^3 coefficient for same-spin correlation",
  "u^4 coefficient for same-spin correlation",
  "u^0 coefficient for opposite-spin correlation",
  "u^1 coefficient for opposite-spin correlation",
  "u^2 coefficient for opposite-spin correlation",
  "u^3 coefficient for opposite-spin correlation",
  "u^4 coefficient for opposite-spin correlation",
  "coefficient for exact exchange"
};
static const double b97_values[B97_N_PAR] =
  {0.8094, 0.5073, 0.7481, 0.0, 0.0,
   0.1737, 2.3487, -2.4868, 0.0, 0.0,
   0.9454, 0.7471, -4.5961, 0.0, 0.0,
   0.1943};
static const double b97_1_values[B97_N_PAR] =
  {0.789518, 0.573805, 0.660975, 0.0, 0.0,
   0.0820011, 2.71681, -2.87103, 0.0, 0.0,
   0.955689, 0.788552, -5.47869, 0.0, 0.0,
   0.21};
static const double b97_2_values[B97_N_PAR] =
  {0.827642, 0.04784, 1.76125, 0.0, 0.0,
   0.585808, -0.691682, 0.394796, 0.0, 0.0,
   0.999849, 1.40626, -7.4406, 0.0, 0.0,
   0.21};
static const double b97_d_values[B97_N_PAR] =
  {1.08662, -0.52127, 3.25429, 0.0, 0.0,
   0.2234, -1.56208, 1.94293, 0.0, 0.0,
   0.69041, 6.3027, -14.9712, 0.0, 0.0,
   0.0};
static const double b97_k_values[B97_N_PAR] =
  {0.507863, 1.46873, -1.51301, 0.0, 0.0,
   0.12355, 2.65399, -3.20694, 0.0, 0.0,
   1.58613, -6.20977, 6.46106, 0.0, 0.0,
   0.42};
static const double b97_3_values[B97_N_PAR] =
  {0.7334648, 0.292527, 3.338789, -10.51158, 10.60907,
   0.5623649, -1.32298, 6.359191, -7.464002, 1.827082,
   1.13383, -2.811967, 7.431302, -1.969342, -11.74423,
   2.692880E-01};
static const double b97_hcth_93_values[B97_N_PAR] =
  {1.0932, -0.744056, 5.5992, -6.78549, 4.49357,
   0.222601, -0.0338622, -0.012517, -0.802496, 1.55396,
   0.729974, 3.35287, -11.543, 8.08564, -4.47857,
   0.0};
static const double b97_hcth_120_values[B97_N_PAR] =
  {1.09163, -0.747215, 5.07833, -4.10746, 1.17173,
   0.489508, -0.260699, 0.432917, -1.99247, 2.48531,
   0.51473, 6.92982, -24.7073, 23.1098, -11.3234,
   0.0};
static const double b97_hcth_147_values[B97_N_PAR] =
  {1.09025, -0.799194, 5.57212, -5.8676, 3.04544,
   0.562576, 0.0171436, -1.30636, 1.05747, 0.885429,
   0.542352, 7.01464, -28.3822, 35.0329, -20.4284,
   0.0};
static const double b97_hcth_407_values[B97_N_PAR] =
  {1.08184, -0.518339, 3.42562, -2.62901, 2.28855,
   1.18777, -2.40292, 5.61741, -9.17923, 6.24798,
   0.589076, 4.42374, -19.2218, 42.5721, -42.0052,
   0.0};
static const double b97_sb98_1a_values[B97_N_PAR] =
  {0.845975, 0.228183, 0.749949, 0.0, 0.0,
   -0.817637, -0.054676, 0.592163, 0.0, 0.0,
   0.975483, 0.398379, -3.7354, 0.0, 0.0,
   0.229015};
static const double b97_sb98_1b_values[B97_N_PAR] =
  {0.800103, -0.084192, 1.47742, 0.0, 0.0,
   1.44946, -2.37073, 2.13564, 0.0, 0.0,
   0.977621, 0.931199, -4.76973, 0.0, 0.0,
   0.199352};
static const double b97_sb98_1c_values[B97_N_PAR] =
  {0.810936, 0.49609, 0.772385, 0.0, 0.0,
   0.262077, 2.12576, -2.30465, 0.0, 0.0,
   0.939269, 0.898121, -4.91276, 0.0, 0.0,
   0.192416};
static const double b97_sb98_2a_values[B97_N_PAR] =
  {0.7492, 0.402322, 0.620779, 0.0, 0.0,
   1.26686, 1.67146, -1.22565, 0.0, 0.0,
   0.964641, 0.050527, -3.01966, 0.0, 0.0,
   0.232055};
static const double b97_sb98_2b_values[B97_N_PAR] =
  {0.770587, 0.180767, 0.955246, 0.0, 0.0,
   0.170473, 1.24051, -0.862711, 0.0, 0.0,
   0.965362, 0.8633, -4.61778, 0.0, 0.0,
   0.237978};
static const double b97_sb98_2c_values[B97_N_PAR] =
  {0.790194, 0.400271, 0.832857, 0.0, 0.0,
   -0.120163, 2.82332, -2.59412, 0.0, 0.0,
   0.934715, 1.14105, -5.33398, 0.0, 0.0,
   0.219847};
static const double b97_gga1_values[B97_N_PAR] =
  {1.1068, -0.8765, 4.2639, 0.0, 0.0,
   0.4883, -2.117, 2.3235, 0.0, 0.0,
   0.7961, 5.706, -14.982, 0.0, 0.0,
   0.0};
static const double b97_hcth_p14_values[B97_N_PAR] =
  {1.03161, -0.360781, 3.51994, -4.95944, 2.41165,
   2.82414, 0.0318843, -1.78512, 2.39795, -0.876909,
   0.0821827, 4.56466, -13.5529, 13.382, -3.17493,
   0.0};
static const double b97_hcth_p76_values[B97_N_PAR] =
  {1.16525, -0.583033, 2.51769, 3.81278, -5.45906,
   -3.92143, -1.10098, -0.091405, -0.859723, 2.07184,
   0.192949, -5.73335, 50.8757, -135.475, 101.268,
   0.0};
static const double b97_hcth_407p_values[B97_N_PAR] =
  {1.08018, -0.4117, 2.4368, 1.389, -1.3529,
   0.80302, -1.0479, 4.9807, -12.89, 9.6446,
   0.73604, 3.027, -10.075, 20.611, -29.418,
   0.0};
static const double b97_1p_values[B97_N_PAR] =
  {0.8773, 0.2149, 1.5204, 0.0, 0.0,
   0.2228, 1.3678, -1.5068, 0.0, 0.0,
   0.9253, 2.027, -7.3431, 0.0, 0.0,
   0.15};
static const double b97_hle16_values[B97_N_PAR] =
  {1.3523, -0.64792375, 4.282025, -3.2862625, 2.8606875,
   0.593885, -1.20146, 2.808705, -4.589615, 3.12399,
   0.294538, 2.21187, -9.6109, 21.28605, -21.0026,
   0.0};


static void
gga_xc_b97_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_xc_b97_params));

  if(p->info->number == XC_HYB_GGA_XC_B97   ||
     p->info->number == XC_HYB_GGA_XC_B97_1 ||
     p->info->number == XC_HYB_GGA_XC_B97_2 ||
     p->info->number == XC_HYB_GGA_XC_B97_K ||
     p->info->number == XC_HYB_GGA_XC_B97_3 ||
     p->info->number == XC_HYB_GGA_XC_SB98_1A ||
     p->info->number == XC_HYB_GGA_XC_SB98_1B ||
     p->info->number == XC_HYB_GGA_XC_SB98_1C ||
     p->info->number == XC_HYB_GGA_XC_SB98_2A ||
     p->info->number == XC_HYB_GGA_XC_SB98_2B ||
     p->info->number == XC_HYB_GGA_XC_SB98_2C ||
     p->info->number == XC_HYB_GGA_XC_B97_1P){
    xc_hyb_init_hybrid(p, 0.0);
  }

}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_xc_b97.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97 = {
  XC_HYB_GGA_XC_B97,
  XC_EXCHANGE_CORRELATION,
  "Becke 97",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Becke1997_8554, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_1 = {
  XC_HYB_GGA_XC_B97_1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-1",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_1_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_2 = {
  XC_HYB_GGA_XC_B97_2,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-2",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Wilson2001_9233, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_2_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_b97_d = {
  XC_GGA_XC_B97_D,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-D",
  XC_FAMILY_GGA,
  {&xc_ref_Grimme2006_1787, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR-1, b97_names, b97_desc, b97_d_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_k = {
  XC_HYB_GGA_XC_B97_K,
  XC_EXCHANGE_CORRELATION,
  "Boese-Martin for Kinetics",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Boese2004_3405, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_k_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_3 = {
  XC_HYB_GGA_XC_B97_3,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-3",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Keal2005_121103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_3_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_93 = {
  XC_GGA_XC_HCTH_93,
  XC_EXCHANGE_CORRELATION,
  "HCTH/93",
  XC_FAMILY_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hcth_93_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_120 = {
  XC_GGA_XC_HCTH_120,
  XC_EXCHANGE_CORRELATION,
  "HCTH/120",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2000_1670, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hcth_120_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_147 = {
  XC_GGA_XC_HCTH_147,
  XC_EXCHANGE_CORRELATION,
  "HCTH/147",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2000_1670, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hcth_147_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_407 = {
  XC_GGA_XC_HCTH_407,
  XC_EXCHANGE_CORRELATION,
  "HCTH/407",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2001_5497, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hcth_407_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_1a = {
  XC_HYB_GGA_XC_SB98_1A,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1a)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_sb98_1a_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_1b = {
  XC_HYB_GGA_XC_SB98_1B,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1b)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_sb98_1b_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_1c = {
  XC_HYB_GGA_XC_SB98_1C,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1c)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_sb98_1c_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_2a = {
  XC_HYB_GGA_XC_SB98_2A,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2a)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_sb98_2a_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_2b = {
  XC_HYB_GGA_XC_SB98_2B,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2b)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_sb98_2b_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_2c = {
  XC_HYB_GGA_XC_SB98_2C,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2c)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_sb98_2c_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_b97_gga1 = {
  XC_GGA_XC_B97_GGA1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97 GGA-1",
  XC_FAMILY_GGA,
  {&xc_ref_Cohen2000_160, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_gga1_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_p14 = {
  XC_GGA_XC_HCTH_P14,
  XC_EXCHANGE_CORRELATION,
  "HCTH p=1/4",
  XC_FAMILY_GGA,
  {&xc_ref_Menconi2001_3958, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hcth_p14_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_p76 = {
  XC_GGA_XC_HCTH_P76,
  XC_EXCHANGE_CORRELATION,
  "HCTH p=7/6",
  XC_FAMILY_GGA,
  {&xc_ref_Menconi2001_3958, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hcth_p76_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_407p = {
  XC_GGA_XC_HCTH_407P,
  XC_EXCHANGE_CORRELATION,
  "HCTH/407+",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2003_5965, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hcth_407p_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_1p = {
  XC_HYB_GGA_XC_B97_1P,
  XC_EXCHANGE_CORRELATION,
  "version of B97 by Cohen and Handy",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Cohen2000_160, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR, b97_names, b97_desc, b97_1p_values, set_ext_params_cpy_exx},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hle16 = {
  XC_GGA_XC_HLE16,
  XC_EXCHANGE_CORRELATION,
  "high local exchange 2016",
  XC_FAMILY_GGA,
  {&xc_ref_Verma2017_380, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {B97_N_PAR - 1, b97_names, b97_desc, b97_hle16_values, set_ext_params_cpy},
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};
