/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_HL  4  /* Hedin & Lundqvist            */
#define XC_LDA_C_GL  5  /* Gunnarson & Lundqvist        */
#define XC_LDA_C_VBH 17 /* von Barth & Hedin            */

typedef struct {
  double hl_r[2], hl_c[2];
} lda_c_hl_params;

#define N_PAR 4
static const char *names[N_PAR] = {"_r0", "_r1", "_c0", "_c1"};
static const char *desc[N_PAR] = {"r0 parameter", "r1 parameter",
                                  "c0 parameter", "c1 parameter"};

static const double par_hl[N_PAR] = {/* HL unpolarized only*/
                                     21.0, 21.0, 0.0225, 0.0225};

static const double par_gl[N_PAR] = {11.4, 15.9, 0.0333, 0.0203};

static const double par_vbh[N_PAR] = {30.0, 75.0, 0.0252, 0.0127};

static void lda_c_hl_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_hl_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_hl.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_hl = {
  XC_LDA_C_HL,
  XC_CORRELATION,
  "Hedin & Lundqvist",
  XC_FAMILY_LDA,
  {&xc_ref_Hedin1971_2064, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_hl, set_ext_params_cpy},
  lda_c_hl_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_gl = {
  XC_LDA_C_GL,
  XC_CORRELATION,
  "Gunnarson & Lundqvist",
  XC_FAMILY_LDA,
  {&xc_ref_Gunnarsson1976_4274, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_gl, set_ext_params_cpy},
  lda_c_hl_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_vbh = {
  XC_LDA_C_VBH,
  XC_CORRELATION,
  "von Barth & Hedin",
  XC_FAMILY_LDA,
  {&xc_ref_vonBarth1972_1629, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_vbh, set_ext_params_cpy},
  lda_c_hl_init, NULL,
  work_lda, NULL, NULL
};
