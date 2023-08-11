/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_K_GDS08_WORKER 100001   /* Combined analytical theory with Monte Carlo sampling */

typedef struct {
  double A, B, C;
} lda_k_gds08_params;

static void
lda_k_gds08_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_k_gds08_params));
}

#define GDS08_N_PAR 3
static const char  *gds08_names[]  = {"_A", "_B", "_C"};
static const char  *gds08_desc[]   = {"linear term", "term proportional to the logarithm of the density", "term proportional to the square of the logarithm"};
static const double gds08_values[] = {0.860, 0.224, 0.0};

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_k_gds08_worker.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_k_gds08_worker = {
  XC_LDA_K_GDS08_WORKER,
  XC_KINETIC,
  "Combined analytical theory with Monte Carlo sampling",
  XC_FAMILY_LDA,
  {&xc_ref_Ghiringhelli2008_073104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {GDS08_N_PAR, gds08_names, gds08_desc, gds08_values, set_ext_params_cpy},
  lda_k_gds08_init, NULL,
  work_lda, NULL, NULL
};
