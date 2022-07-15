/*
 Copyright (C) 2014-2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include <ctype.h>
#include "util.h"

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || (id) == XC_FAMILY_HYB_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA || (id) == XC_FAMILY_HYB_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA || (id) == XC_FAMILY_HYB_LDA ||  is_gga(id))

int main(int argc, char **argv)
{
  int i, func_id, error, npar;
  xc_func_type func;
  char *fname;

  if(argc!=2) {
    printf("Usage: %s [ func_id | func_name ]\n",argv[0]);
    return 1;
  }

  /* Print libxc information */
  printf("libxc version %s\n",xc_version_string());
  printf("%s\n", xc_reference());
  printf("doi: %s\n", xc_reference_doi());
  printf("\n");

  /* Is functional defined by a string constant? */
  if(isalpha(argv[1][0]))
    func_id = xc_functional_get_number(argv[1]);
  else
    func_id = atoi(argv[1]);

  /* Initialize functional */
  error = xc_func_init(&func, func_id, XC_UNPOLARIZED);
  if(error) {
    printf("Functional '%s' not found.\n", argv[1]);
    return 1;
  }

  /* Get functional name */
  fname = xc_functional_get_name(func_id);

  /* Print out info */
  printf("%10s: %-20i\n%10s: %-25s\n", "func_id", func_id, "name", fname);
  printf("%10s: %-20s\n%10s: %-25s\n", "family", get_family(&func), "kind", get_kind(&func));
  printf("%10s: %s\n", "comment", func.info->name);

  /* Print out hybrid exchange info */
  if(func.info->family==XC_FAMILY_HYB_LDA || func.info->family==XC_FAMILY_HYB_GGA || func.info->family==XC_FAMILY_HYB_MGGA) {
    /* Range separation? */
    int rangesep=0;
    if(func.info->flags & XC_FLAGS_HYB_CAM)
      rangesep++;
    if(func.info->flags & XC_FLAGS_HYB_CAMY)
      rangesep++;
    if(func.info->flags & XC_FLAGS_HYB_LC)
      rangesep++;
    if(func.info->flags & XC_FLAGS_HYB_LCY)
      rangesep++;

    if(rangesep) {
      double alpha, beta, omega;
      xc_hyb_cam_coef(&func,&omega,&alpha,&beta);
      printf("\nThis is a range-separated hybrid functional with range-separation constant % .3f,\n",omega);
      printf("and %4.1f%% short-range and %4.1f%% long-range exact exchange,\n",(alpha+beta)*100,(alpha)*100);

      if(func.info->flags & XC_FLAGS_HYB_CAM || func.info->flags & XC_FLAGS_HYB_LC)
        printf("using the error function kernel.\n");
      else if(func.info->flags & XC_FLAGS_HYB_CAMY || func.info->flags & XC_FLAGS_HYB_LCY)
        printf("using the Yukawa kernel.\n");
    } else {
      double alpha=xc_hyb_exx_coef(&func);
      printf("\nThis is a global hybrid functional with %4.1f%% of exact exchange.\n",alpha*100);
    }
  } else {
    if(func.info->kind == XC_EXCHANGE || func.info->kind == XC_EXCHANGE_CORRELATION)
      printf("\nThis is a pure functional with no exact exchange.\n");
  }

  printf("\nReference(s):\n");
  for(i = 0; i < 5; i++){
    if(func.info->refs[i] == NULL) break;
    printf("  *) %s\n", func.info->refs[i]->ref);
    if(strlen(func.info->refs[i]->doi) > 0){
       printf("     doi: %s\n", func.info->refs[i]->doi);
    }
  }

  printf("\nImplementation has support for:\n");
  if(func.info->flags & XC_FLAGS_HAVE_EXC)
    printf("  *) energy\n");
  if(func.info->flags & XC_FLAGS_HAVE_VXC)
    printf("  *) first derivative\n");
  if(func.info->flags & XC_FLAGS_HAVE_FXC)
    printf("  *) second derivative\n");
  if(func.info->flags & XC_FLAGS_HAVE_KXC)
    printf("  *) third derivative\n");
  if(func.info->flags & XC_FLAGS_HAVE_KXC)
    printf("  *) fourth derivative\n");

  printf("\nDefault thresholds:\n");
  printf("density: %e\n",func.dens_threshold);
  printf("   zeta: %e\n",func.zeta_threshold);
  if(is_gga(func.info->family))
    printf("  sigma: %e\n",func.sigma_threshold);
  if(is_mgga(func.info->family))
    printf("    tau: %e\n",func.tau_threshold);
  
  /* Query parameters */
  npar = xc_func_info_get_n_ext_params(func.info);
  if(npar > 0) {
    printf("\nFunctional has %i external parameters:\n",npar);
    printf("%3s %13s %8s %s\n","idx","value","name","description");
    for(i = 0; i < npar; i++)
      printf("%3i % e %8s %s\n", i,
             xc_func_info_get_ext_params_default_value(func.info, i),
             xc_func_info_get_ext_params_name(func.info, i),
             xc_func_info_get_ext_params_description(func.info, i));
  } else {
    printf("\nFunctional has no external parameters.\n");
  }

  /* Free memory */
  xc_func_end(&func);
  libxc_free(fname);

  return 0;
}
