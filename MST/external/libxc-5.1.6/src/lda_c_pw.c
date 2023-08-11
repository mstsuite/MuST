/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
 Correlation energy per-particle and potential of a HEG as parameterized
 by
   J.P. Perdew & Y. Wang
   Ortiz & Ballone

Note that the PW modified corresponds to the version of PW used in the
original PBE routine. This amounts to adding some more digits in some of
the constants of PW.
************************************************************************/

#define XC_LDA_C_PW      12   /* Perdew & Wang                             */
#define XC_LDA_C_PW_MOD  13   /* Perdew & Wang (Modified)                  */
#define XC_LDA_C_OB_PW   14   /* Ortiz & Ballone (PW)                      */
#define XC_LDA_C_PW_RPA  25   /* Perdew & Wang fit of the RPA              */
#define XC_LDA_C_UPW92  683   /* Ruggeri, Rios, and Alavi unrestricted fit */
#define XC_LDA_C_RPW92  684   /* Ruggeri, Rios, and Alavi restricted fit   */

typedef struct {
  double pp[3], a[3], alpha1[3];
  double beta1[3], beta2[3], beta3[3], beta4[3];
  double fz20;
} lda_c_pw_params;

#define PW_N_PAR 22
static const char  *pw_names[PW_N_PAR]  = {"_pp[0]", "_pp[1]", "_pp[2]","_a[0]", "_a[1]", "_a[2]", "_alpha1[0]", "_alpha1[1]", "_alpha1[2]", "_beta1[0]", "_beta1[1]", "_beta1[2]", "_beta2[0]", "_beta2[1]", "_beta2[2]", "_beta3[0]", "_beta3[1]", "_beta3[2]", "_beta4[0]", "_beta4[1]", "_beta4[2]", "_fz20"};
/* These should be documented better */
static const char  *pw_desc[PW_N_PAR]   = {"pp0", "pp1", "pp2","a0", "a1", "a2", "alpha10", "alpha11", "alpha12", "beta10", "beta11", "beta12", "beta20", "beta21", "beta22", "beta30", "beta31", "beta32", "beta40", "beta41", "beta42", "fz20"};

static const double par_pw[PW_N_PAR] = {
  1.0,  1.0,  1.0,
  0.031091,  0.015545,   0.016887,
  0.21370,  0.20548,  0.11125,
  7.5957, 14.1189, 10.357,
  3.5876, 6.1977, 3.6231,
  1.6382, 3.3662, 0.88026,
  0.49294, 0.62517, 0.49671,
  1.709921
};

static const double par_pw_mod[PW_N_PAR] = {
  1.0,  1.0,  1.0,
  0.0310907, 0.01554535, 0.0168869,
  0.21370,  0.20548,  0.11125,
  7.5957, 14.1189, 10.357,
  3.5876, 6.1977, 3.6231,
  1.6382, 3.3662,  0.88026,
  0.49294, 0.62517, 0.49671,
  1.709920934161365617563962776245
};

static const double par_ob[PW_N_PAR] = {
  1.0,  1.0,  1.0,
  0.031091,  0.015545, 0.016887,
  0.026481, 0.022465, 0.11125,
  7.5957, 14.1189, 10.357,
  3.5876, 6.1977, 3.6231,
  -0.46647, -0.56043, 0.88026,
  0.13354, 0.11313, 0.49671,
  1.709921
};

/* The parameters fixed by the low and high-density limits were
   taken from the pw_mod functional, i.e. they contain more
   significant digits than the vanilla pw */
static const double par_pw_rpa[PW_N_PAR] = {
  0.75, 0.75, 1.0,
  0.031091,  0.015545,   0.016887,
  0.082477, 0.035374, 0.028829,
   5.1486, 6.4869, 10.357,
  1.6483, 1.3083, 3.6231,
  0.23647, 0.15180, 0.47990,
  0.20614, 0.082349, 0.12279,
  1.709921
};

static const double par_upw92[PW_N_PAR] = {
  1.0,  1.0,  1.0,
  0.0310907, 0.01554535, 0.0168869,
  0.227012,  0.264193,  0.11125,
  7.5957, 14.1189, 10.357,
  3.5876, 6.1977, 3.6231,
  1.76522, 4.78287,  0.88026,
  0.523918, 0.750424, 0.49671,
  1.709920934161365617563962776245
};

static const double par_rpw92[PW_N_PAR] = {
  1.0,  1.0,  1.0,
  0.0310907, 0.01554535, 0.0168869,
  0.21370,  0.266529,  0.11125,
  7.5957, 14.1189, 10.357,
  3.5876, 6.1977, 3.6231,
  1.6382, 4.86059,  0.88026,
  0.49294, 0.750188, 0.49671,
  1.709920934161365617563962776245
};

static void
lda_c_pw_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_pw_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_pw.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_pw = {
  XC_LDA_C_PW,
  XC_CORRELATION,
  "Perdew & Wang",
  XC_FAMILY_LDA,
  {&xc_ref_Perdew1992_13244, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PW_N_PAR, pw_names, pw_desc, par_pw, set_ext_params_cpy},
  lda_c_pw_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_pw_mod = {
  XC_LDA_C_PW_MOD,
  XC_CORRELATION,
  "Perdew & Wang (modified)",
  XC_FAMILY_LDA,
  {&xc_ref_Perdew1992_13244_mod, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PW_N_PAR, pw_names, pw_desc, par_pw_mod, set_ext_params_cpy},
  lda_c_pw_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_ob_pw = {
  XC_LDA_C_OB_PW,
  XC_CORRELATION,
  "Ortiz & Ballone (PW parametrization)",
  XC_FAMILY_LDA,
  {&xc_ref_Ortiz1994_1391, &xc_ref_Ortiz1994_1391_err, &xc_ref_Perdew1992_13244_mod, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PW_N_PAR, pw_names, pw_desc, par_ob, set_ext_params_cpy},
  lda_c_pw_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_pw_rpa = {
  XC_LDA_C_PW_RPA,
  XC_CORRELATION,
  "Perdew & Wang (fit to the RPA energy)",
  XC_FAMILY_LDA,
  {&xc_ref_Perdew1992_13244, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PW_N_PAR, pw_names, pw_desc, par_pw_rpa, set_ext_params_cpy},
  lda_c_pw_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_upw92 = {
  XC_LDA_C_UPW92,
  XC_CORRELATION,
  "Ruggeri, Rios, and Alavi unrestricted fit",
  XC_FAMILY_LDA,
  {&xc_ref_Ruggeri2018_161105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PW_N_PAR, pw_names, pw_desc, par_upw92, set_ext_params_cpy},
  lda_c_pw_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_rpw92 = {
  XC_LDA_C_RPW92,
  XC_CORRELATION,
  "Ruggeri, Rios, and Alavi restricted fit",
  XC_FAMILY_LDA,
  {&xc_ref_Ruggeri2018_161105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {PW_N_PAR, pw_names, pw_desc, par_rpw92, set_ext_params_cpy},
  lda_c_pw_init, NULL,
  work_lda, NULL, NULL
};
