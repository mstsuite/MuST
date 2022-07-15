/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_PBE          101 /* Perdew, Burke & Ernzerhof exchange             */
#define XC_GGA_X_PBE_R        102 /* Perdew, Burke & Ernzerhof exchange (revised)   */
#define XC_GGA_X_PBE_SOL      116 /* Perdew, Burke & Ernzerhof exchange (solids)    */
#define XC_GGA_X_XPBE         123 /* xPBE reparametrization by Xu & Goddard         */
#define XC_GGA_X_PBE_JSJR     126 /* JSJR reparametrization by Pedroza, Silva & Capelle */
#define XC_GGA_X_PBEK1_VDW    140 /* PBE reparametrization for vdW                  */
#define XC_GGA_X_APBE         184 /* mu fixed from the semiclassical neutral atom   */
#define XC_GGA_X_PBE_TCA       59 /* PBE revised by Tognetti et al                  */
#define XC_GGA_X_PBE_MOL       49 /* Del Campo, Gazquez, Trickey and Vela (PBE-like) */
#define XC_GGA_X_LAMBDA_LO_N   45 /* lambda_LO(N) version of PBE                    */
#define XC_GGA_X_LAMBDA_CH_N   44 /* lambda_CH(N) version of PBE                    */
#define XC_GGA_X_LAMBDA_OC2_N  40 /* lambda_OC2(N) version of PBE                   */
#define XC_GGA_X_BCGP          38 /* Burke, Cancio, Gould, and Pittalis             */
#define XC_GGA_X_PBEFE        265 /* PBE for formation energies                     */

typedef struct{
  double kappa, mu;
  double lambda;   /* parameter used in the Odashima & Capelle versions */
} gga_x_pbe_params;


static void
gga_x_pbe_init(xc_func_type *p)
{
  gga_x_pbe_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pbe_params));
  params = (gga_x_pbe_params *) (p->params);

  /* This has to be explicitly initialized here */
  params->lambda = 0.0;
}

#define PBE_N_PAR 2
static const char  *pbe_names[PBE_N_PAR]  = {"_kappa", "_mu"};
static const char  *pbe_desc[PBE_N_PAR]   = {
  "Asymptotic value of the enhancement function",
  "Coefficient of the 2nd order expansion"};

static const double pbe_values[PBE_N_PAR] =
  {0.8040, MU_PBE};
static const double pbe_r_values[PBE_N_PAR] =
  {1.245, MU_PBE};
static const double pbe_sol_values[PBE_N_PAR] =
  {0.804, MU_GE};
static const double pbe_xpbe_values[PBE_N_PAR] =
  {0.91954, 0.23214};
static const double pbe_jsjr_values[PBE_N_PAR] =
  {0.8040, 0.046*M_PI*M_PI/3.0};
static const double pbe_k1_vdw_values[PBE_N_PAR] =
  {1.0, MU_PBE};
static const double pbe_apbe_values[PBE_N_PAR] =
  {0.8040, 0.260};
static const double pbe_tca_values[PBE_N_PAR] =
  {1.227, MU_PBE};
static const double pbe_mol_values[PBE_N_PAR] =
  {0.8040, 0.27583};
static const double pbe_bcgp_values[PBE_N_PAR] =
  {0.8040, 0.249};
static const double pbe_fe_values[PBE_N_PAR] =
  {0.437, 0.346};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_pbe.c"
#include "work_gga.c"


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1996_3865, &xc_ref_Perdew1996_3865_err, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_r = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Revised PBE from Zhang & Yang",
  XC_FAMILY_GGA,
  {&xc_ref_Zhang1998_890, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_r_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_sol = {
  XC_GGA_X_PBE_SOL,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew2008_136406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_sol_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_xpbe = {
  XC_GGA_X_XPBE,
  XC_EXCHANGE,
  "Extended PBE by Xu & Goddard III",
  XC_FAMILY_GGA,
  {&xc_ref_Xu2004_4068, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_xpbe_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_jsjr = {
  XC_GGA_X_PBE_JSJR,
  XC_EXCHANGE,
  "Reparametrized PBE by Pedroza, Silva & Capelle",
  XC_FAMILY_GGA,
  {&xc_ref_Pedroza2009_201106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_jsjr_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbek1_vdw = {
  XC_GGA_X_PBEK1_VDW,
  XC_EXCHANGE,
  "Reparametrized PBE for vdW",
  XC_FAMILY_GGA,
  {&xc_ref_Klimes2010_022201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_k1_vdw_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_apbe = {
  XC_GGA_X_APBE,
  XC_EXCHANGE,
  "mu fixed from the semiclassical neutral atom",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_apbe_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_tca = {
  XC_GGA_X_PBE_TCA,
  XC_EXCHANGE,
  "PBE revised by Tognetti et al",
  XC_FAMILY_GGA,
  {&xc_ref_Tognetti2008_536, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_tca_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_mol = {
  XC_GGA_X_PBE_MOL,
  XC_EXCHANGE,
  "Reparametrized PBE by del Campo, Gazquez, Trickey & Vela",
  XC_FAMILY_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_mol_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_bcgp = {
  XC_GGA_X_BCGP,
  XC_EXCHANGE,
  "Burke, Cancio, Gould, and Pittalis",
  XC_FAMILY_GGA,
  {&xc_ref_Burke2014_4834, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_bcgp_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbefe = {
  XC_GGA_X_PBEFE,
  XC_EXCHANGE,
  "PBE for formation energies",
  XC_FAMILY_GGA,
  {&xc_ref_Perez2015_3844, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_fe_values, set_ext_params_cpy},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#define PBEL_N_PAR 3
static const char  *pbe_lambda_names[PBEL_N_PAR]  = {"_N", "_kappa", "_mu"};
static const char  *pbe_lambda_desc[PBEL_N_PAR]   = {
  "Number of electrons",
  "Asymptotic value of the enhancement function",
  "Coefficient of the 2nd order expansion"};
static const double pbe_lambda_lo_n_values[PBEL_N_PAR] =
  {1e23, MU_PBE, 2.273};
static const double pbe_lambda_ch_n_values[PBEL_N_PAR] =
  {1e23, MU_PBE, 2.215};
static const double pbe_lambda_oc2_n_values[PBEL_N_PAR] =
  {1e23, MU_PBE, 2.00};

static void
pbe_lambda_set_ext_params(xc_func_type *p, const double *ext_params)
{
  const double lambda_1 = 1.48;

  gga_x_pbe_params *params;
  double lambda, N;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_pbe_params *) (p->params);

  N              = get_ext_param(p, ext_params, 0);
  params->mu     = get_ext_param(p, ext_params, 1);
  params->lambda = get_ext_param(p, ext_params, 2);

  lambda = (1.0 - 1.0/N)*params->lambda + lambda_1/N;
  params->kappa = lambda/M_CBRT2 - 1.0;
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lambda_lo_n = {
  XC_GGA_X_LAMBDA_LO_N,
  XC_EXCHANGE,
  "lambda_LO(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBEL_N_PAR, pbe_lambda_names, pbe_lambda_desc, pbe_lambda_lo_n_values, pbe_lambda_set_ext_params},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lambda_ch_n = {
  XC_GGA_X_LAMBDA_CH_N,
  XC_EXCHANGE,
  "lambda_CH(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBEL_N_PAR, pbe_lambda_names, pbe_lambda_desc, pbe_lambda_ch_n_values, pbe_lambda_set_ext_params},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lambda_oc2_n = {
  XC_GGA_X_LAMBDA_OC2_N,
  XC_EXCHANGE,
  "lambda_OC2(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PBEL_N_PAR, pbe_lambda_names, pbe_lambda_desc, pbe_lambda_oc2_n_values, pbe_lambda_set_ext_params},
  gga_x_pbe_init, NULL,
  NULL, work_gga, NULL
};


