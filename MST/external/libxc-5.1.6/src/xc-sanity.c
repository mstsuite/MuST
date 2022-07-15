/*
 Copyright (C) 2017 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xc.h"

static int is_cam(xc_func_type *p) {
  return (p->info->flags & XC_FLAGS_HYB_CAM)
    || (p->info->flags & XC_FLAGS_HYB_LC);
}

static int is_yukawa(xc_func_type *p) {
  return (p->info->flags & XC_FLAGS_HYB_CAMY)
    || (p->info->flags & XC_FLAGS_HYB_LCY);
}

static int is_rangesep(xc_func_type *p) {
  return is_cam(p) || is_yukawa(p);
}

int is_hybrid(xc_func_type *p) {
  return (p->info->family == XC_FAMILY_HYB_LDA) ||
    (p->info->family == XC_FAMILY_HYB_GGA) ||
    (p->info->family == XC_FAMILY_HYB_MGGA);
}

int main(void) {
  int i, N, error;
  xc_func_type func;

  /* List of functionals */
  int *flist;
  char *p;

  /* Get list of available functionals */
  N = xc_number_of_functionals();
  flist = (int *) malloc(N*sizeof(int));
  xc_available_functional_numbers(flist);

  /* Loop over functionals */
  for(i=0;i<N;i++) {
    /* Functional id and name */
    int func_id;
    char *fname;
    /* Kind and family of functional */
    char kind[5], family[10];

    printf("Checking functional %i -> %i\n", i, flist[i]);

    func_id = flist[i];
    fname = xc_functional_get_name(func_id);

    /* Initialize functional */
    error = xc_func_init(&func, func_id, XC_UNPOLARIZED);
    if(error) {
      printf("Error initializing functional %i.\n", func_id);
      fflush(stdout);
      return 1;
    }

    /* Check kind is consistent with name */
    switch(func.info->kind) {
    case(XC_EXCHANGE):
      strcpy(kind, "_x_");
      break;

    case(XC_CORRELATION):
      strcpy(kind, "_c_");
      break;

    case(XC_EXCHANGE_CORRELATION):
      strcpy(kind, "_xc_");
      break;

    case(XC_KINETIC):
      strcpy(kind, "_k_");
      break;

    default:
      fprintf(stderr,"Kind %i not handled.\n",func.info->kind);
      return 2;
    }
    p = strstr(fname, kind);
    if(p == NULL)
      printf("Functional %i '%s' name may be inconsistent with its kind '%s'.\n",func_id, fname, kind);

    /* check if hybrid is initialized */
    if(is_rangesep(&func)) {
      if(func.cam_omega == 0.0)
        printf("Range-separated hybrid does not seem to be initialized: omega is zero\n");
      if(func.cam_alpha == 0.0 && func.cam_beta == 0.0)
        printf("Range-separated hybrid does not seem to be initialized: alpha and beta are zero\n");
    } else if(strncmp(fname, "hyb_", 4) == 0) {
      if(func.cam_alpha == 0.0)
        printf("Hybrid does not seem to be initialized: alpha is zero\n");
      if(func.cam_beta != 0.0)
        printf("Hybrid has non-zero beta\n");
      if(func.cam_omega != 0.0)
        printf("Hybrid has non-zero omega\n");
    } else {
      if(func.cam_alpha != 0.0)
        printf("Non-hybrid functional has non-zero alpha\n");
      if(func.cam_beta != 0.0)
        printf("Non-hybrid functional has non-zero beta\n");
      if(func.cam_omega != 0.0)
        printf("Non-hybrid functional '%s' has non-zero omega\n", fname);
    }
    if(is_rangesep(&func) && !is_hybrid(&func))
      printf("Fuctional is range-separated but is not marked hybrid\n");

    /* Check family is consistent with name */
    family[0] = '\0';
    switch(func.info->family) {
    case(XC_FAMILY_LDA):
      strcat(family, "lda_");
      break;

    case(XC_FAMILY_GGA):
      strcat(family, "gga_");
      break;

    case(XC_FAMILY_MGGA):
      strcat(family, "mgga_");
      break;

    case(XC_FAMILY_HYB_LDA):
      strcat(family, "hyb_lda_");
      break;

    case(XC_FAMILY_HYB_GGA):
      strcat(family, "hyb_gga_");
      break;

    case(XC_FAMILY_HYB_MGGA):
      strcat(family, "hyb_mgga_");
      break;

    default:
      fprintf(stderr, "Family %i not handled.\n", func.info->family);
      return 2;
    }

    p = strstr(fname, family);
    if(p != fname)
      printf("Functional %i '%s' name may be inconsistent with its family '%s'.\n",func_id, fname, family);

    /* Check non-local correlation parameters */
    {
      double b, C;
      xc_nlc_coef(&func, &b, &C);

      if(func.info->flags & XC_FLAGS_VV10) {
        if(b == 0.0)
          printf("Functional %i '%s' is supposed to have VV10 but has zero b.\n",func_id, fname);
        if(C == 0.0)
          printf("Functional %i '%s' is supposed to have VV10 but has zero b.\n",func_id, fname);
      } else {
        if(b != 0.0)
          printf("Functional %i '%s' isn't supposed to long-range correlation but has non-zero b.\n",func_id, fname);
        if(C != 0.0)
          printf("Functional %i '%s' isn't supposed to long-range correlation but has non-zero C.\n",func_id, fname);
      }
    }

    /* Free memory */
    free(fname);
    xc_func_end(&func);
  }

  /* Free memory */
  free(flist);

  return 0;
}

