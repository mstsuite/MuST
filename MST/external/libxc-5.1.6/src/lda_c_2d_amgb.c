/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
 Correlation energy per particle and potentials for a homogeneous electron
 gas in 2D, as parametrized by Attaccalite et al.
************************************************************************/

#define XC_LDA_C_2D_AMGB  15   /* Attaccalite et al             */

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_2d_amgb.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_2d_amgb = {
  XC_LDA_C_2D_AMGB,
  XC_CORRELATION,
  "AMGB (for 2D systems)",
  XC_FAMILY_LDA,
  {&xc_ref_Attaccalite2002_256601, NULL, NULL, NULL, NULL},
  XC_FLAGS_2D | MAPLE2C_FLAGS,
  1e-9,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  work_lda, NULL, NULL
};
