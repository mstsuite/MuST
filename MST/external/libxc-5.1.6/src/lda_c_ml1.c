/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_C_ML1    22   /* Modified LSD (version 1) of Proynov and Salahub */
#define XC_LDA_C_ML2    23   /* Modified LSD (version 2) of Proynov and Salahub */

typedef struct {
  double fc, q;
} lda_c_ml1_params;

#define N_PAR 2
static const char  *names[N_PAR]  = {"_fc", "_q"};
/* These should be documented better */
static const char  *desc[N_PAR]   = {"fc", "q"};

static const double par_ml1[N_PAR] = {0.2026, 0.084};
static const double par_ml2[N_PAR] = {0.266, 0.5};

static void
lda_c_ml1_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_ml1_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_ml1.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_ml1 = {
  XC_LDA_C_ML1,
  XC_CORRELATION,
  "Modified LSD (version 1) of Proynov and Salahub",
  XC_FAMILY_LDA,
  {&xc_ref_Proynov1994_7874, &xc_ref_Proynov1998_12616, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_ml1, set_ext_params_cpy},
  lda_c_ml1_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_ml2 = {
  XC_LDA_C_ML2,
  XC_CORRELATION,
  "Modified LSD (version 2) of Proynov and Salahub",
  XC_FAMILY_LDA,
  {&xc_ref_Proynov1994_7874, &xc_ref_Proynov1998_12616, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_ml2, set_ext_params_cpy},
  lda_c_ml1_init, NULL,
  work_lda, NULL, NULL
};
