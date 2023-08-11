#include <stdio.h>

#include <xc.h>

int main()
{
  xc_func_type func;
  double rho[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
  double sigma[5] = {0.2, 0.3, 0.4, 0.5, 0.6};
  double exc[5];
  int i, vmajor, vminor, vmicro, func_id = 1;

  xc_version(&vmajor, &vminor, &vmicro);
  printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);

  if(xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", func_id);
    return 1;
  }

  switch(func.info->family)
  {
  case XC_FAMILY_LDA:
    xc_lda_exc(&func, 5, rho, exc);
    break;
  case XC_FAMILY_GGA:
  case XC_FAMILY_HYB_GGA:
    xc_gga_exc(&func, 5, rho, sigma, exc);
    break;
  }

  for(i=0; i<5; i+=1){
    printf("%lf %lf\n", rho[i], exc[i]);
  }

  xc_func_end(&func);
}
