/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_WIGNER    2   /* Wigner parametrization       */
#define XC_LDA_XC_LP_A   547   /* Lee-Parr reparametrization A */
#define XC_LDA_XC_LP_B   548   /* Lee-Parr reparametrization B */
#define XC_LDA_C_MCWEENY 551   /* McWeeny 76 */
#define XC_LDA_C_BR78    552   /* Brual & Rothstein 78 */
#define XC_LDA_C_OW_LYP  573   /* Wigner with corresponding LYP parameters */
#define XC_LDA_C_OW      574   /* Optimized Wigner */

typedef struct {
  double a, b;
} lda_c_wigner_params;

#define N_PAR 2
static const char  *names[N_PAR]  = {"_a", "_b"};
static const char  *desc[N_PAR]   = {"a parameter", "b parameter"};

static const double val_wigner[N_PAR] = {-0.44, 7.8};
static const double val_lp_a[N_PAR] = {-0.8626*RS_FACTOR, 0.0};
static const double val_lp_b[N_PAR] = {-0.906*RS_FACTOR, 2.1987e-2*RS_FACTOR};
static const double val_mcweeny[N_PAR] = {-RS_FACTOR/2.946, RS_FACTOR*9.652/2.946};
static const double val_br78[N_PAR] = {-RS_FACTOR/21.437, RS_FACTOR*9.810/21.437};
static const double val_ow_lyp[N_PAR] = {-0.04918*RS_FACTOR/0.349, RS_FACTOR/0.349};
static const double val_ow[N_PAR] = {-0.0526*RS_FACTOR/0.349, RS_FACTOR/0.349};

static void
lda_c_wigner_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_wigner_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_wigner.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_wigner = {
  XC_LDA_C_WIGNER,
  XC_CORRELATION,
  "Wigner",
  XC_FAMILY_LDA,
  {&xc_ref_Wigner1938_678, &xc_ref_Stewart1995_4337, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {N_PAR, names, desc, val_wigner, set_ext_params_cpy},
  lda_c_wigner_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_lp_a = {
  XC_LDA_XC_LP_A,
  XC_EXCHANGE_CORRELATION,
  "Lee-Parr reparametrization A",
  XC_FAMILY_LDA,
  {&xc_ref_Lee1990_193, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {N_PAR, names, desc, val_lp_a, set_ext_params_cpy},
  lda_c_wigner_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_lp_b = {
  XC_LDA_XC_LP_B,
  XC_EXCHANGE_CORRELATION,
  "Lee-Parr reparametrization B",
  XC_FAMILY_LDA,
  {&xc_ref_Lee1990_193, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {N_PAR, names, desc, val_lp_b, set_ext_params_cpy},
  lda_c_wigner_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_mcweeny = {
  XC_LDA_C_MCWEENY,
  XC_CORRELATION,
  "McWeeny 76",
  XC_FAMILY_LDA,
  {&xc_ref_McWeeny1976_3, &xc_ref_Brual1978_1177, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {N_PAR, names, desc, val_mcweeny, set_ext_params_cpy},
  lda_c_wigner_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_br78 = {
  XC_LDA_C_BR78,
  XC_CORRELATION,
  "Brual & Rothstein 78",
  XC_FAMILY_LDA,
  {&xc_ref_Brual1978_1177, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {N_PAR, names, desc, val_br78, set_ext_params_cpy},
  lda_c_wigner_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_ow_lyp = {
  XC_LDA_C_OW_LYP,
  XC_CORRELATION,
  "Wigner with corresponding LYP parameters",
  XC_FAMILY_LDA,
  {&xc_ref_Stewart1995_4337, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {N_PAR, names, desc, val_ow_lyp, set_ext_params_cpy},
  lda_c_wigner_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_ow = {
  XC_LDA_C_OW,
  XC_CORRELATION,
  "Optimized Wigner",
  XC_FAMILY_LDA,
  {&xc_ref_Stewart1995_4337, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {N_PAR, names, desc, val_ow, set_ext_params_cpy},
  lda_c_wigner_init, NULL,
  work_lda, NULL, NULL
};

