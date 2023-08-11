/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
Correlation functional by Pittalis, Rasanen & Marques for the 2D electron gas
************************************************************************/

/* TODO: convert this to an (rs, zeta) expression */

#define XC_LDA_C_2D_PRM  16   /* Pittalis, Rasanen & Marques correlation in 2D */

typedef struct{
  double N;
  double c;
} lda_c_2d_prm_params;

/* Initialization */
static void
lda_c_2d_prm_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);

  p->params = libxc_malloc(sizeof(lda_c_2d_prm_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_2d_prm.c"
#include "work_lda.c"

static const char  *N_names[]  = {"N"};
static const char  *N_desc[]   = {"Number of electrons"};
static const double N_values[] = {2.0};

static void
N_set_ext_params(xc_func_type *p, const double *ext_params)
{
  static double prm_q = 3.9274; /* 2.258 */
  lda_c_2d_prm_params *params;

  assert(p != NULL && p->params != NULL);
  params = (lda_c_2d_prm_params *) (p->params);

  params->N = get_ext_param(p, ext_params, 0);

  if(params->N <= 1.0){
    fprintf(stderr, "PRM functional cannot be used for N_electrons <= 1\n");
    exit(1);
  }

  params->c = M_PI/(2.0*(params->N - 1.0)*prm_q*prm_q); /* Eq. (13) */
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_2d_prm = {
  XC_LDA_C_2D_PRM,
  XC_CORRELATION,
  "PRM (for 2D systems)",
  XC_FAMILY_LDA,
  {&xc_ref_Pittalis2008_195322, NULL, NULL, NULL, NULL},
  XC_FLAGS_2D | MAPLE2C_FLAGS,
  1e-32,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  lda_c_2d_prm_init, NULL,
  work_lda, NULL, NULL
};
