/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_CAM_B3LYP        433 /* CAM version of B3LYP */
#define XC_HYB_GGA_XC_CAMH_B3LYP       614 /* CAM version of B3LYP tuned for tddft */
#define XC_HYB_GGA_XC_TUNED_CAM_B3LYP  434 /* CAM version of B3LYP tuned for excitations */
#define XC_HYB_GGA_XC_RCAM_B3LYP       610 /* Similar to CAM-B3LYP, but trying to reduce the many-electron self-interaction */
#define XC_HYB_GGA_XC_CAM_PBEH         681 /* CAM version of PBEH */
#define XC_HYB_GGA_XC_CAM_QTP_00       490 /* CAM-QTP-00 */
#define XC_HYB_GGA_XC_CAM_QTP_01       482 /* CAM-QTP-01 */
#define XC_HYB_GGA_XC_CAM_QTP_02       491 /* CAM-QTP-02 */
#define XC_HYB_GGA_XC_LC_QTP           492 /* LC-QTP     */
#define XC_HYB_GGA_XC_MCAM_B3LYP       640 /* Modified CAM-B3LYP */

#define CAMB3_N_PAR 4
static const char  *camb3_names[CAMB3_N_PAR]  = {"_alpha", "_beta", "_omega", "_ac"};
static const char  *camb3_desc[CAMB3_N_PAR]   = {
  "Fraction of Hartree-Fock exchange",
  "Fraction of short-range exact exchange",
  "Range separation parameter",
  "Fraction of LYP correlation"
};

static const double par_cam_b3lyp[CAMB3_N_PAR]       = {0.65, -0.46,   0.33,  0.81};
static const double par_camh_b3lyp[CAMB3_N_PAR]      = {0.50, -0.31,   0.33,  0.81};
static const double par_tuned_cam_b3lyp[CAMB3_N_PAR] = {1.00, -0.9201, 0.15,  0.81};
static const double par_cam_qtp_00[CAMB3_N_PAR]      = {0.91, -0.37,   0.29,  0.8};
static const double par_cam_qtp_01[CAMB3_N_PAR]      = {1.00, -0.77,   0.31,  0.8};
static const double par_cam_qtp_02[CAMB3_N_PAR]      = {1.00, -0.72,   0.335, 1.0};
static const double par_lc_qtp[CAMB3_N_PAR]          = {1.00, -1.00,   0.475, 1.0};
static const double par_mcam_b3lyp[CAMB3_N_PAR]      = {0.38, -0.19,   0.33,  0.81};

static void
set_cam_params(xc_func_type *p, const double *ext_params) {
  p->cam_alpha = get_ext_param(p, ext_params, 0);
  p->cam_beta  = get_ext_param(p, ext_params, 1);
  p->cam_omega = get_ext_param(p, ext_params, 2);
}

static void
camb3_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double alpha, beta, omega, ac;

  assert(p != NULL);
  alpha  = get_ext_param(p, ext_params, 0);
  beta   = get_ext_param(p, ext_params, 1);
  omega  = get_ext_param(p, ext_params, 2);
  ac     = get_ext_param(p, ext_params, 3);

  p->mix_coef[0] = 1.0 - alpha;
  p->mix_coef[1] = -beta;
  p->mix_coef[2] = 1.0 - ac;
  p->mix_coef[3] = ac;

  xc_func_set_ext_params_name(p->func_aux[1], "_omega", omega);
  set_cam_params(p, ext_params);
}

void
xc_hyb_gga_xc_cam_b3lyp_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_GGA_X_B88, XC_GGA_X_ITYH, XC_LDA_C_VWN, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0};
  xc_mix_init(p, 4, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_b3lyp = {
  XC_HYB_GGA_XC_CAM_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAM version of B3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Yanai2004_51, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_cam_b3lyp, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_camh_b3lyp = {
  XC_HYB_GGA_XC_CAMH_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAM version of B3LYP, tuned for TDDFT",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Shao2020_587, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_camh_b3lyp, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_tuned_cam_b3lyp = {
  XC_HYB_GGA_XC_TUNED_CAM_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAM version of B3LYP, tuned for excitations and properties",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Okuno2012_29, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_tuned_cam_b3lyp, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_qtp_00 = {
  XC_HYB_GGA_XC_CAM_QTP_00,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Verma2014_18A534, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_cam_qtp_00, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_qtp_01 = {
  XC_HYB_GGA_XC_CAM_QTP_01,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Jin2016_034107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_cam_qtp_01, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_qtp_02 = {
  XC_HYB_GGA_XC_CAM_QTP_02,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Haiduke2018_184106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_cam_qtp_02, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_qtp = {
  XC_HYB_GGA_XC_LC_QTP,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Haiduke2018_184106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_lc_qtp, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_mcam_b3lyp = {
  XC_HYB_GGA_XC_MCAM_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "Modified CAM-B3LYP by Day, Nguyen and Pachter",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Day2006_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {CAMB3_N_PAR, camb3_names, camb3_desc, par_mcam_b3lyp, camb3_set_ext_params},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};


#define RCAM_N_PAR 4
static const char  *rcam_names[RCAM_N_PAR]  = {"_alpha", "_beta", "_omega", "_ab88"};
static const char  *rcam_desc[RCAM_N_PAR]   = {
  "Fraction of Hartree-Fock exchange",
  "Fraction of short-range exact exchange",
  "Range separation parameter",
  "Fraction of B88 exchange"
};

static const double par_rcam_b3lyp[RCAM_N_PAR] = {0.18352+0.94979, -0.94979, 0.33, 0.95238};

void
xc_hyb_gga_xc_rcam_b3lyp_init(xc_func_type *p)
{
  static int funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_ITYH, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0};
  xc_mix_init(p, 4, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

static void
rcam_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double alpha_libxc, beta_libxc;
  double alpha, beta, omega, cb88;

  assert(p != NULL);
  /* connection is

     libxc_alpha = alpha + beta
     libxc_beta = -beta
  */
  alpha_libxc = get_ext_param(p, ext_params, 0);
  beta_libxc = get_ext_param(p, ext_params, 1);
  alpha  =  alpha_libxc + beta_libxc;
  beta   = -beta_libxc;
  omega  =  get_ext_param(p, ext_params, 2);
  cb88   =  get_ext_param(p, ext_params, 3);

  p->mix_coef[0] = 1.0 - alpha - cb88;
  p->mix_coef[1] = cb88 - beta;
  p->mix_coef[2] = beta;
  p->mix_coef[3] = 1.0;

  xc_func_set_ext_params_name(p->func_aux[2], "_omega", omega);
  set_cam_params(p, ext_params);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_rcam_b3lyp = {
  XC_HYB_GGA_XC_RCAM_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "Similar to CAM-B3LYP, but trying to reduce the many-electron self-interaction",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Cohen2007_191109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {RCAM_N_PAR, rcam_names, rcam_desc, par_rcam_b3lyp, rcam_set_ext_params},
  xc_hyb_gga_xc_rcam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#define CAM_N_PAR 3
static const char  *cam_names[CAM_N_PAR]  = {"_alpha", "_beta", "_omega"};
static const char  *cam_desc[CAM_N_PAR]   = {
  "Fraction of Hartree-Fock exchange",
  "Fraction of short-range exact exchange",
  "Range separation parameter"
};

static const double par_cam_pbeh[CAM_N_PAR] = {0.2, 0.8, 0.7};

static void
cam_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double alpha, beta, omega;

  assert(p != NULL);
  alpha  = get_ext_param(p, ext_params, 0);
  beta   = get_ext_param(p, ext_params, 1);
  omega  = get_ext_param(p, ext_params, 2);

  p->mix_coef[0] = 1.0 - alpha;
  p->mix_coef[1] = -beta;

  xc_func_set_ext_params_name(p->func_aux[1], "_omega", omega);
  set_cam_params(p, ext_params);
}

static void
hyb_gga_xc_cam_pbeh_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_GGA_X_PBE, XC_GGA_X_HJS_PBE, XC_GGA_C_PBE};
  static double funcs_coef[3] = {0.0, 0.0, 1.0}; /* the first two get set by ext_params */

  xc_mix_init(p, 3, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_pbeh = {
  XC_HYB_GGA_XC_CAM_PBEH,
  XC_EXCHANGE_CORRELATION,
  "CAM hybrid screened exchange PBE version",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Chen2018_073803, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {CAM_N_PAR, cam_names, cam_desc, par_cam_pbeh, cam_set_ext_params},
  hyb_gga_xc_cam_pbeh_init,
  NULL, NULL, NULL, NULL
};
