/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* for a review on the values of lambda and gamma, please see EV
Ludena and VV Karasiev, in "Reviews of Modern Quantum Chemistry: a
Celebration of the Contributions of Robert G. Parr, edited by KD Sen
(World Scientific, Singapore, 2002), p. 612.
 */

#define XC_GGA_K_TFVW          52  /* Thomas-Fermi plus von Weiszaecker correction */
#define XC_GGA_K_VW            500 /* von Weiszaecker functional */
#define XC_GGA_K_GE2           501 /* Second-order gradient expansion (l = 1/9) */
#define XC_GGA_K_GOLDEN        502 /* TF-lambda-vW form by Golden (l = 13/45) */
#define XC_GGA_K_YT65          503 /* TF-lambda-vW form by Yonei and Tomishima (l = 1/5) */
#define XC_GGA_K_BALTIN        504 /* TF-lambda-vW form by Baltin (l = 5/9) */
#define XC_GGA_K_LIEB          505 /* TF-lambda-vW form by Lieb (l = 0.185909191) */
#define XC_GGA_K_ABSP1         506 /* gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)] */
#define XC_GGA_K_ABSP2         507 /* gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)] */
#define XC_GGA_K_ABSP3         277 /* gamma-TFvW form by Acharya et al [g = 1 - 1.513/N^0.35] */
#define XC_GGA_K_ABSP4         278 /* gamma-TFvW form by Acharya et al [g = l = 1/(1 + 1.332/N^(1/3))] */
#define XC_GGA_K_GR            508 /* gamma-TFvW form by Gazquez and Robles */
#define XC_GGA_K_LUDENA        509 /* gamma-TFvW form by Ludena */
#define XC_GGA_K_GP85          510 /* gamma-TFvW form by Ghosh and Parr */
#define XC_GGA_K_TFVW_OPT      635 /* empirically optimized gamma-TFvW form */

typedef struct{
  double lambda, gamma;
} gga_k_tflw_params;


static void
gga_k_tflw_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_k_tflw_params));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_k_tflw.c"
#include "work_gga.c"

#define TFVW_N_PAR 2
static const char  *tfvw_names[TFVW_N_PAR]  = {"_lambda", "_gamma"};
static const char  *tfvw_desc[TFVW_N_PAR]   = {
  "Lambda", "Gamma"
};
static const double tfvw_values[TFVW_N_PAR]        = {1.0, 1.0};
static const double tfvw_vw_values[TFVW_N_PAR]     = {1.0, 0.0};
static const double tfvw_ge2_values[TFVW_N_PAR]    = {1.0/9.0, 1.0};
static const double tfvw_golden_values[TFVW_N_PAR] = {13.0/45.0, 1.0};
static const double tfvw_yt65_values[TFVW_N_PAR]   = {1.0/5.0, 1.0};
static const double tfvw_baltin_values[TFVW_N_PAR] = {5.0/9.0, 1.0};
static const double tfvw_lieb_values[TFVW_N_PAR]   = {0.185909191, 1.0}; /* 1/5.37897... */
static const double tfvw_opt_values[TFVW_N_PAR]  = {0.599, 0.697};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_tfvw = {
  XC_GGA_K_TFVW,
  XC_KINETIC,
  "Thomas-Fermi plus von Weiszaecker correction",
  XC_FAMILY_GGA,
  {&xc_ref_Weizsacker1935_431, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_vw = {
  XC_GGA_K_VW,
  XC_KINETIC,
  "von Weiszaecker correction to Thomas-Fermi",
  XC_FAMILY_GGA,
  {&xc_ref_Weizsacker1935_431, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_vw_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_ge2 = {
  XC_GGA_K_GE2,
  XC_KINETIC,
  "Second-order gradient expansion of the kinetic energy density",
  XC_FAMILY_GGA,
  {&xc_ref_Kompaneets1956_427, &xc_ref_Kirznits1957_115, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_ge2_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_golden = {
  XC_GGA_K_GOLDEN,
  XC_KINETIC,
  "TF-lambda-vW form by Golden (l = 13/45)",
  XC_FAMILY_GGA,
  {&xc_ref_Golden1957_604, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_golden_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_yt65 = {
  XC_GGA_K_YT65,
  XC_KINETIC,
  "TF-lambda-vW form by Yonei and Tomishima (l = 1/5)",
  XC_FAMILY_GGA,
  {&xc_ref_Yonei1965_1051, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_yt65_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_baltin = {
  XC_GGA_K_BALTIN,
  XC_KINETIC,
  "TF-lambda-vW form by Baltin (l = 5/9)",
  XC_FAMILY_GGA,
  {&xc_ref_Baltin1972_1176, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_baltin_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_lieb = {
  XC_GGA_K_LIEB,
  XC_KINETIC,
  "TF-lambda-vW form by Lieb (l = 0.185909191)",
  XC_FAMILY_GGA,
  {&xc_ref_Lieb1981_603, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_lieb_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_tfvw_opt = {
  XC_GGA_K_TFVW_OPT,
  XC_KINETIC,
  "empirically optimized gamma-TFvW form",
  XC_FAMILY_GGA,
  {&xc_ref_EspinosaLeal2015_31463, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {TFVW_N_PAR, tfvw_names, tfvw_desc, tfvw_opt_values, set_ext_params_cpy},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};


static const char  *N_names[]  = {"N"};
static const char  *N_desc[]   = {"Number of electrons"};
static const double N_values[] = {1.0};

static void
N_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double C0 = CBRT(M_PI/3.0);
  double C1 = CBRT(M_PI*M_PI/36.0)/6.0 - CBRT(M_PI*M_PI/9.0)/4.0;
  double N;
  gga_k_tflw_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_k_tflw_params *) (p->params);

  N = get_ext_param(p, ext_params, 0);

  params->gamma = 1.0;
  params->lambda = 1.0;

  switch(p->info->number){
  case XC_GGA_K_ABSP1:
    params->gamma = 1.0 - 1.412/CBRT(N);
    break;
  case XC_GGA_K_ABSP2:
    params->gamma = 1.0 - 1.332/CBRT(N);
    break;
  case XC_GGA_K_ABSP3:
    params->gamma = 1.0 - 1.513/pow(N, 0.35);
    break;
  case XC_GGA_K_ABSP4:
    params->gamma = 1.0/(1.0 + 1.332/CBRT(N));
    params->lambda = params->gamma;
    break;
  case XC_GGA_K_GR:
    params->gamma = (1.0 - 2.0/N)*(1.0 - C0/CBRT(N) + C1*CBRT(N*N));
    break;
  case XC_GGA_K_LUDENA:
    params->gamma = CBRT(6.0*M_PI)*M_PI*M_PI*(1.0 - 1.0/(N*N));
    break;
  case XC_GGA_K_GP85:
    params->gamma = CBRT(6.0*M_PI*M_PI)*M_PI*M_PI/4.0*
      (1.0 - 1.0/N)*(1.0 + 1.0/N + 6.0/(N*N));
    break;
  }
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_absp1 = {
  XC_GGA_K_ABSP1,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_absp2 = {
  XC_GGA_K_ABSP2,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_absp3 = {
  XC_GGA_K_ABSP3,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = 1 - 1.513/N^0.35]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_absp4 = {
  XC_GGA_K_ABSP4,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = l = 1/(1 + 1.332/N^(1/3))]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_gr = {
  XC_GGA_K_GR,
  XC_KINETIC,
  "gamma-TFvW form by Gazquez and Robles",
  XC_FAMILY_GGA,
  {&xc_ref_Gazquez1982_1467, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_ludena = {
  XC_GGA_K_LUDENA,
  XC_KINETIC,
  "gamma-TFvW form by Ludena",
  XC_FAMILY_GGA,
  {&xc_ref_Ludena1986, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_gp85 = {
  XC_GGA_K_GP85,
  XC_KINETIC,
  "gamma-TFvW form by Ghosh and Parr",
  XC_FAMILY_GGA,
  {&xc_ref_Ghosh1985_3307, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  gga_k_tflw_init, NULL,
  NULL, work_gga, NULL
};
