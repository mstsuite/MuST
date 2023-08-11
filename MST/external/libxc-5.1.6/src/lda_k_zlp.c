/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
 LDA parametrization of Vosko, Wilk & Nusair
************************************************************************/

#define XC_LDA_K_ZLP     550   /* kinetic energy version of ZLP */

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_k_zlp.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_k_zlp = {
  XC_LDA_K_ZLP,
  XC_KINETIC,
  "Wigner including kinetic energy contribution",
  XC_FAMILY_LDA,
  {&xc_ref_Fuentealba1995_31, &xc_ref_Zhao1993_918, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  work_lda, NULL, NULL
};
