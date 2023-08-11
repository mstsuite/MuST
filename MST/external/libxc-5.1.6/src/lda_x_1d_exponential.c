/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_X_1D_EXPONENTIAL  600 /* Exchange in 1D for an exponentially screened interaction */

typedef struct{
  double beta;         /* screening parameter beta */
} lda_x_1d_exponential_params;

static void
lda_x_1d_exponential_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_x_1d_exponential_params));
}

GPU_FUNCTION
static inline double FT_inter(double x)
{
  return xc_E1_scaled(x*x);
}

GPU_FUNCTION
static void func1(double *x, int n, void *dummy)
{
  int ii;

  for(ii=0; ii<n; ii++)
    x[ii] = FT_inter(x[ii]);
}

GPU_FUNCTION
static void func2(double *x, int n, void *dummy)
{
  int ii;

  for(ii=0; ii<n; ii++)
    x[ii] = x[ii]*FT_inter(x[ii]);
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_x_1d_exponential.c"
#include "work_lda.c"

static const char  *exp1d_names[]  = {"beta"};
static const char  *exp1d_desc[]   = {"Parameter of the exponential"};
static const double exp1d_values[] = {1.0};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x_1d_exponential = {
  XC_LDA_X_1D_EXPONENTIAL,
  XC_EXCHANGE,
  "Exchange in 1D for an exponentially screened interaction",
  XC_FAMILY_LDA,
  {&xc_ref_Helbig2011_032503, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D | MAPLE2C_FLAGS,
  1e-26,
  {1, exp1d_names, exp1d_desc, exp1d_values, set_ext_params_cpy},
  lda_x_1d_exponential_init, NULL,
  work_lda, NULL, NULL
};
