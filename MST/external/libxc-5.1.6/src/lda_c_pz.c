/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
 Correlation energy per particle and potential of a HEG as parametrized
 by
   Perdew & Zunger
   Ortiz & Ballone
************************************************************************/

#define XC_LDA_C_PZ       9   /* Perdew & Zunger              */
#define XC_LDA_C_PZ_MOD  10   /* Perdew & Zunger (Modified)   */
#define XC_LDA_C_OB_PZ   11   /* Ortiz & Ballone (PZ)         */

typedef struct {
  double gamma[2];
  double beta1[2];
  double beta2[2];
  double a[2], b[2], c[2], d[2];
} lda_c_pz_params;

#define N_PAR 14
static const char *names[N_PAR]  = {
  "_gamma0", "_gamma1",
  "_beta10", "_beta11",
  "_beta20", "_beta21",
  "_a0", "_a1",
  "_b0", "_b1",
  "_c0", "_c1",
  "_d0", "_d1"
};

static const char *desc[N_PAR]  = {
  "gamma0 parameter", "gamma1 parameter",
  "beta10 parameter", "beta11 parameter",
  "beta20 parameter", "beta21 parameter",
  "a0 parameter", "a1 parameter",
  "b0 parameter", "b1 parameter",
  "c0 parameter", "c1 parameter",
  "d0 parameter", "d1 parameter"
};

static const double par_pz[N_PAR] = {
  -0.1423, -0.0843,  /* gamma */
   1.0529,  1.3981,  /* beta1 */
   0.3334,  0.2611,  /* beta2 */
   0.0311,  0.01555, /*  a    */
  -0.048,  -0.0269,  /*  b    */
   0.0020,  0.0007,  /*  c    */
  -0.0116, -0.0048   /*  d    */
};

static const double par_pz_mod[N_PAR] = {
  -0.1423, -0.0843,
   1.0529,  1.3981,
   0.3334,  0.2611,
   0.0311,  0.01555,
  -0.048,  -0.0269,
   0.0020191519406228,  0.00069255121311694,
  -0.0116320663789130, -0.00480126353790614
};

static const double par_pz_ob[N_PAR] = {
  -0.103756, -0.065951,
   0.56371,   1.11846,
   0.27358,   0.18797,
   0.031091,  0.015545,
  -0.046644, -0.025599,
   0.00419,   0.00329,  /* the sign of c[0] and c[1] is different from [2], but is consistent
                             with the continuity requirement. There is nothing in [3] about this. */
  -0.00983,  -0.00300
};

static void
lda_c_pz_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_pz_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_pz.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_pz = {
  XC_LDA_C_PZ,
  XC_CORRELATION,
  "Perdew & Zunger",
  XC_FAMILY_LDA,
  {&xc_ref_Perdew1981_5048, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_pz, set_ext_params_cpy},
  lda_c_pz_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_pz_mod = {
  XC_LDA_C_PZ_MOD,
  XC_CORRELATION,
  "Perdew & Zunger (Modified)",
  XC_FAMILY_LDA,
  {&xc_ref_Perdew1981_5048_mod, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_pz_mod, set_ext_params_cpy},
  lda_c_pz_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_ob_pz = {
  XC_LDA_C_OB_PZ,
  XC_CORRELATION,
  "Ortiz & Ballone (PZ parametrization)",
  XC_FAMILY_LDA,
  {&xc_ref_Ortiz1994_1391, &xc_ref_Ortiz1994_1391_err, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_pz_ob, set_ext_params_cpy},
  lda_c_pz_init, NULL,
  work_lda, NULL, NULL
};
