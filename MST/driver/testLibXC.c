#include <stdio.h>

#include <xc.h>

int main()
{
  xc_func_type x_func;
  xc_func_type c_func;
  xc_func_type xc_func;
  double rho[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
  double sigma[5] = {0.2, 0.3, 0.4, 0.5, 0.6};
  double ex[5];
  double vx[5];
  double vx2[5];
  double ec[5];
  double vc[5];
  double vc2[5];
  double exc[5];
  double vxc[5];
  double vxc2[5];
  int i, vmajor, vminor, vmicro;
  int x_func_id = 1;
  int c_func_id = 17;
  int xc_func_id = 20;

  xc_version(&vmajor, &vminor, &vmicro);
  printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);

  if(xc_func_init(&x_func, x_func_id, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", x_func_id);
    return 1;
  }
  if(xc_func_init(&c_func, c_func_id, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", c_func_id);
    return 1;
  }

  printf("\nThe functional '%s' is ", x_func.info->name);
  switch (x_func.info->kind) {
  case (XC_EXCHANGE):
    printf("an exchange functional");
    break;
  case (XC_CORRELATION):
    printf("a correlation functional");
    break;
  case (XC_EXCHANGE_CORRELATION):
    printf("an exchange-correlation functional");
    break;
  case (XC_KINETIC):
    printf("a kinetic energy functional");
    break;
  default:
    printf("of unknown kind");
    break;
  }
  printf(", it belongs to the '%s'", x_func.info->name);
  switch (x_func.info->family) {
  case (XC_FAMILY_LDA):
    printf("LDA");
    break;
  case (XC_FAMILY_GGA):
    printf("GGA");
    break;
  case (XC_FAMILY_HYB_GGA):
    printf("Hybrid GGA");
    break;
  case (XC_FAMILY_MGGA):
    printf("MGGA");
    break;
  case (XC_FAMILY_HYB_MGGA):
    printf("Hybrid MGGA");
    break;
  default:
    printf("unknown");
    break;
  }
  printf("' family and is defined in the reference(s):\n");

  for(i = 0; x_func.info->refs[i] != NULL; i++){
    printf("[%d] %s\n", i+1, x_func.info->refs[i]->ref);
  }

  switch(x_func.info->family)
    {
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&x_func, 5, rho, ex, vx);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc_vxc(&x_func, 5, rho, sigma, ex, vx, vx2);
      break;
    }

  for(i=0; i<5; i+=1){
    printf("Exchange functional: %lf %lf %lf\n", rho[i], 2.0e0*vx[i], 2.0e0*ex[i]);
  }

  //=====================================================
  //

  printf("\nThe functional '%s' is ", c_func.info->name);
  switch (c_func.info->kind) {
  case (XC_EXCHANGE):
    printf("an exchange functional");
    break;
  case (XC_CORRELATION):
    printf("a correlation functional");
    break;
  case (XC_EXCHANGE_CORRELATION):
    printf("an exchange-correlation functional");
    break;
  case (XC_KINETIC):
    printf("a kinetic energy functional");
    break;
  default:
    printf("of unknown kind");
    break;
  }
  printf(", it belongs to the '%s'", c_func.info->name);
  switch (c_func.info->family) {
  case (XC_FAMILY_LDA):
    printf("LDA");
    break;
  case (XC_FAMILY_GGA):
    printf("GGA");
    break;
  case (XC_FAMILY_HYB_GGA):
    printf("Hybrid GGA");
    break;
  case (XC_FAMILY_MGGA):
    printf("MGGA");
    break;
  case (XC_FAMILY_HYB_MGGA):
    printf("Hybrid MGGA");
    break;
  default:
    printf("unknown");
    break;
  }
  printf("' family and is defined in the reference(s):\n");

  for(i = 0; c_func.info->refs[i] != NULL; i++){
    printf("[%d] %s\n", i+1, c_func.info->refs[i]->ref);
  }

  switch(c_func.info->family)
    {
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&c_func, 5, rho, ec, vc);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc_vxc(&c_func, 5, rho, sigma, ec, vc, vc2);
      break;
    }

  for(i=0; i<5; i+=1){
    printf("Correlation functional: %lf %lf %lf\n", rho[i], 2.0e0*vc[i], 2.0e0*ec[i]);
  }


  xc_func_end(&x_func);
  xc_func_end(&c_func);

  printf("\n");
  for(i=0; i<5; i+=1){
    printf("Exchange-correlation: %lf %lf %lf\n", rho[i], 2.0e0*(vx[i]+vc[i]), 2.0e0*(ex[i]+ec[i]));
  }
}
