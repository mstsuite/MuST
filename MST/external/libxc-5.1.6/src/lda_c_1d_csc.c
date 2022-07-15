/*
 Copyright (C) 2006-2009 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_1D_CSC          18 /* Casula, Sorella, and Senatore 1D correlation     */

typedef struct{
  double para[10], ferro[10];

  int interaction;  /* 0: exponentially screened; 1: soft-Coulomb */
  double bb;         /* screening parameter */

} lda_c_1d_csc_params;

static const double par_para[][10] = { /* paramagnetic */
  /* 0:A    1:B   2:C    3:D    4:E  5:n1   6:n2  7:alpha  8:beta  9:m */
  {  4.66,  0.0,  2.092, 3.735, 0.0, 1.379, 2.0, 23.63,  109.9,    1.837}, /* exponentially screened interaction */
  {  9.5,   0.0,  1.85,  5.64,  0.0, 0.882, 2.0,  5.346,   6.69,   3.110},
  { 16.40,  0.0,  2.90,  6.235, 0.0, 0.908, 2.0,  3.323,   2.23,   3.368},
  { 22.53,  0.0,  2.09,  7.363, 0.0, 0.906, 2.0,  2.029,   0.394,  4.070},
  { 32.1,   0.0,  3.77,  7.576, 0.0, 0.941, 2.0,  1.63,    0.198,  4.086},
  {110.5,   0.0,  7.90,  8.37,  0.0, 1.287, 2.0,  1.399,   0.0481, 4.260},
  {413.0,   0.0, 10.8,   7.99,  0.0, 1.549, 2.0,  1.308,   0.0120, 4.165},
  { 7.40, 1.120, 1.890, 0.0964,  0.0250,   2.0, 3.0, 2.431, 0.0142, 2.922}, /* soft-Coulomb interaction */
  {18.40, 0.0,   7.501, 0.10185, 0.012827, 2.0, 3.0, 1.511, 0.258,  4.424}
};

static const double par_ferro[][10] = { /* ferromagnetic */
  { 5.24, 0.0,   1.568, 0.12856, 0.003201, 2.0, 3.0, 0.0538, 1.56e-5, 2.958}
};

static void
lda_c_1d_csc_init(xc_func_type *p)
{
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_1d_csc_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_1d_csc.c"
#include "work_lda.c"

static const char  *csc_names[]  = {"interaction", "beta"};
static const char  *csc_desc[]   = {"0 (exponentially screened) | 1 (soft-Coulomb)", "Screening parameter"};
static const double csc_values[] = {1, 1.0};

static void
csc_set_ext_params(xc_func_type *p, const double *ext_params)
{
  lda_c_1d_csc_params *params;
  int ii;

  assert(p != NULL && p->params != NULL);
  params = (lda_c_1d_csc_params *)(p->params);

  params->interaction = (int)round(get_ext_param(p, ext_params, 0));
  params->bb = get_ext_param(p, ext_params, 1);

  const double * ppara = NULL;
  const double * pferro = NULL;

  if(params->interaction == 0){
    if      (params->bb == 0.1){
      ppara  = par_para[0];
      pferro = par_para[0];
    }else if(params->bb == 0.3){
      ppara  = par_para[1];
      pferro = par_para[1];
    }else if(params->bb == 0.5){
      ppara  = par_para[2];
      pferro = par_para[2];
    }else if(params->bb == 0.75){
      ppara  = par_para[3];
      pferro = par_para[3];
    }else if(params->bb == 1.0){
      ppara  = par_para[4];
      pferro = par_para[4];
    }else if(params->bb == 2.0){
      ppara  = par_para[5];
      pferro = par_para[5];
    }else if(params->bb == 4.0){
      ppara  = par_para[6];
      pferro = par_para[6];
    }
  }else if(params->interaction == 1){
    if     (params->bb == 0.5){
      ppara  = par_para[7];
      pferro = par_para[7];
    }else if(params->bb == 1.0){
      ppara  = par_para[8];
      pferro = par_ferro[0];
    }
  }

  if(ppara == NULL){
    fprintf(stderr, "Invalid value of parameters (inter,b) = (%d,%f) in lda_c_1d_csc_set_params",
	    params->interaction, params->bb);
    exit(1);
  }

  //we must copy the values (instead of pointing to them) so that they are available on the GPU
  for(ii = 0; ii < 10; ii++){
    params->para[ii] = ppara[ii];
    params->ferro[ii] = pferro[ii];
  }

}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_1d_csc = {
  XC_LDA_C_1D_CSC,
  XC_CORRELATION,
  "Casula, Sorella & Senatore",
  XC_FAMILY_LDA,
  {&xc_ref_Casula2006_245427, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D | MAPLE2C_FLAGS,
  1e-25,
  {2, csc_names, csc_desc, csc_values, csc_set_ext_params},
  lda_c_1d_csc_init, NULL,
  work_lda, NULL, NULL
};
