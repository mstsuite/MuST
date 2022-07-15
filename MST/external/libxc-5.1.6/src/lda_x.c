
/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019      Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_X             1   /* Exchange                            */
#define XC_LDA_C_XALPHA      6   /* Slater Xalpha                       */
#define XC_LDA_X_RAE       549   /* Rae self-energy corrected exchange  */
#define XC_HYB_LDA_XC_LDA0 177   /* LDA0: hybrid LDA exchange           */

/*
    Slater's Xalpha functional (Exc = alpha Ex)

    Note: this is to be added to the exchange

    This correlation functional, added to the exchange functional, produces
    a total exchange-correlation functional, Exc, equal to 3/2 * alpha * Ex
    Setting alpha equal to one gives the *usual* Slater Xalpha functional,
    whereas alpha equal to 2/3 just leaves the exchange functional unchanged.
*/

/* Range separation
    J. Toulouse, A. Savin, H.-J. Flad, Int. J. of Quant. Chem. 100, 1047-1056 (2004).
*/

typedef struct{
  double alpha;       /* parameter for Xalpha functional */
} lda_x_params;

static void
lda_x_init(xc_func_type *p)
{
  lda_x_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_x_params));
  params = (lda_x_params *) (p->params);

  params->alpha = 1.0;
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_x.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x = {
  XC_LDA_X,
  XC_EXCHANGE,
  "Slater exchange",
  XC_FAMILY_LDA,
  {&xc_ref_Dirac1930_376, &xc_ref_Bloch1929_545, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  lda_x_init, NULL,
  work_lda, NULL, NULL
};

static const char  *xalpha_names[]  = {"alpha"};
static const char  *xalpha_desc[]   = {"X-alpha multiplicative parameter"};
static const double xalpha_values[] = {1.0};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  lda_x_params *params;

  assert(p != NULL && p->params != NULL);
  params = (lda_x_params *)(p->params);

  params->alpha = 1.5*get_ext_param(p, ext_params, 0) - 1.0;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_xalpha = {
  XC_LDA_C_XALPHA,
  XC_CORRELATION,
  "Slater's Xalpha",
  XC_FAMILY_LDA,
  {&xc_ref_Slater1951_385, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, xalpha_names, xalpha_desc, xalpha_values, set_ext_params},
  lda_x_init, NULL,
  work_lda, NULL, NULL
};

static const char  *N_names[]  = {"N"};
static const char  *N_desc[]   = {"Number of electrons"};
static const double N_values[] = {1.0};

static void
N_set_ext_params(xc_func_type *p, const double *ext_params)
{
  lda_x_params *params;
  double N, dx, dx2;

  assert(p != NULL && p->params != NULL);
  params = (lda_x_params *)(p->params);

  N = get_ext_param(p, ext_params, 0);

  dx  = 1.0/CBRT(4.0*N);
  dx2 = dx*dx;
  params->alpha = 1.0 - 8.0/3.0*dx + 2.0*dx2 - dx2*dx2/3.0;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x_rae = {
  XC_LDA_X_RAE,
  XC_EXCHANGE,
  "Rae self-energy corrected exchange",
  XC_FAMILY_LDA,
  {&xc_ref_Rae1973_574, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {1, N_names, N_desc, N_values, N_set_ext_params},
  lda_x_init, NULL,
  work_lda, NULL, NULL
};

/* Patrick Rinke confirmed that this functional only contains
   75% of LDA correlation. */
static void
hyb_lda_xc_lda0_init(xc_func_type *p)
{
  static int    funcs_id[2] = {XC_LDA_X, XC_LDA_C_PW_MOD};
  static double funcs_coef[2] = {1.0 - 0.25, 1.0 - 0.25};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_hybrid(p, 0.25);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_lda_xc_lda0 = {
  XC_HYB_LDA_XC_LDA0,
  XC_EXCHANGE_CORRELATION,
  "LDA hybrid exchange (LDA0)",
  XC_FAMILY_HYB_LDA,
  {&xc_ref_Rinke2012_126404, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  hyb_lda_xc_lda0_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};
