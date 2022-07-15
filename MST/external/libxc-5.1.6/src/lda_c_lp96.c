/*
 Copyright (C) 2017 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_LP96      289   /* Liu-Parr correlation */
#define XC_LDA_K_LP96      580   /* Liu-Parr kinetic */

typedef struct {
  double C1, C2, C3;
} lda_c_lp96_params;

#define N_PAR 3
static const char *names[N_PAR]  = {
  "_C1", "_C2", "_C3"
};
static const char *desc[N_PAR]  = {
  "C1 parameter", "C2 parameter", "C3 parameter"
};

static const double par_c_lp96[N_PAR] = {-0.0603,   0.0175, -0.00053};
static const double par_k_lp96[N_PAR] = { 0.03777, -0.01002, 0.00039};

static void
lda_c_lp96_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_lp96_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_lp96.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_lp96 = {
  XC_LDA_C_LP96,
  XC_CORRELATION,
  "Liu-Parr correlation",
  XC_FAMILY_LDA,
  {&xc_ref_Liu1996_2211, &xc_ref_Liu2000_29, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-16,
  {N_PAR, names, desc, par_c_lp96, set_ext_params_cpy},
  lda_c_lp96_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_k_lp96 = {
  XC_LDA_K_LP96,
  XC_KINETIC,
  "Liu-Parr kinetic",
  XC_FAMILY_LDA,
  {&xc_ref_Liu1996_2211, &xc_ref_Liu2000_29, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-16,
  {N_PAR, names, desc, par_k_lp96, set_ext_params_cpy},
  lda_c_lp96_init, NULL,
  work_lda, NULL, NULL
};
