/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <xc.h>
#include <util.h>

void printVar(int nspin, int nn, const char *cvar, double *dvar)
{
  const char *format = "%3d: %20s[%2d] = %#19.12E\n";
  static int n =0;
  int ii;

  printf(format, n++, cvar, 0, dvar[0]);
  if(nspin == XC_UNPOLARIZED)
    return;

  for(ii=1; ii<nn; ii++)
    printf(format, n++, cvar, ii, dvar[ii]);
}

/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  xc_func_type func;
  const xc_func_info_type *info;

  int functional;
  int nspin;

  /* Input */
  double rho[2];         /* rhoa, rhob */
  double sigma[3];       /* sigmaaa, sigmaab, sigmabb */
  double lapl[2];        /* lapla, laplb */
  double tau[2];         /* taua, taub */

  /* Output */
  double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, );

  if(argc != 12){
    printf("Usage:\n%s funct pol rhoa rhob sigmaaa sigmaab sigmabb lapla laplb taua taub\n", argv[0]);
    return 1;
  }

  /* Is functional defined by a string constant? */
  if(isalpha(argv[1][0]))
    functional = xc_functional_get_number(argv[1]);
  else
    functional = atoi(argv[1]);

  nspin      = atoi(argv[2]);
  if(nspin == XC_UNPOLARIZED)
    printf("Unpolarized calculation\n");
  else if(nspin == XC_POLARIZED)
    printf("Polarized calculation\n");
  else {
    printf("Invalid value for pol input.\n");
    exit(1);
  }


  rho[0]     = atof(argv[3]);
  rho[1]     = atof(argv[4]);
  sigma[0]   = atof(argv[5]);
  sigma[1]   = atof(argv[6]);
  sigma[2]   = atof(argv[7]);
  lapl[0]    = atof(argv[8]);
  lapl[1]    = atof(argv[9]);
  tau[0]     = atof(argv[10]);
  tau[1]     = atof(argv[11]);

  if(nspin == 1){
    rho[0]   += rho[1];
    sigma[0] += 2.0*sigma[1] + sigma[2];
    lapl[0]  += lapl[1];
    tau[0]   += tau[1];
  }

  if(xc_func_init(&func, functional, nspin) != 0){
    fprintf(stderr, "Functional '%d' not found\n", functional);
    exit(1);
  }
  info = func.info;

  /* allocate buffers */
  zk MGGA_OUT_PARAMS_NO_EXC(=, ) = NULL;
  xc_mgga_vars_allocate_all(info->family, 1, &(func.dim),
                            1, 1, 1, 1, 1,
                            &zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));

  xc_mgga_evaluate_functional(&func, 1, rho, sigma, lapl, tau,
                              zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));

  /* transform to energy per volume */
  if(nspin == XC_UNPOLARIZED){
    zk[0] *= rho[0];
  }else{
    zk[0] *= rho[0] + rho[1];
  }

  printf(" rhoa= %#0.2E", rho[0]);
  if(nspin == XC_POLARIZED) printf(" rhob= %#0.2E", rho[1]);
  if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
    printf(" sigmaaa= %#0.2E", sigma[0]);
    if(nspin == XC_POLARIZED) printf(" sigmaab= %#0.2E", sigma[1]);
    if(nspin == XC_POLARIZED) printf(" sigmabb= %#0.2E", sigma[2]);
  }
  if(info->family == XC_FAMILY_MGGA){
    printf(" lapla= %#0.2E", lapl[0]);
    if(nspin == XC_POLARIZED) printf(" laplb= %#0.2E",  lapl[1]);
    printf(" taua= %#0.2E", tau[0]);
    if(nspin == XC_POLARIZED) printf(" taub= %#0.2E", tau[1]);
  }
  printf("\n\n");

  if(info->flags & XC_FLAGS_HAVE_EXC){
    printVar(nspin, func.dim.zk, "zk", zk);
    printf("\n");
  }

  if(info->flags & XC_FLAGS_HAVE_VXC){
    printVar(nspin, func.dim.vrho, "vrho", vrho);

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.vsigma, "vsigma", vsigma);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.vlapl, "vlapl", vlapl);
      printVar(nspin, func.dim.vtau, "vtau", vtau);
    }
    printf("\n");
  }

  if(info->flags & XC_FLAGS_HAVE_FXC){
    printVar(nspin, func.dim.v2rho2, "v2rho2", v2rho2);

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v2rhosigma, "v2rhosigma", v2rhosigma);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v2rholapl, "v2rholapl", v2rholapl);
      printVar(nspin, func.dim.v2rhotau, "v2rhotau", v2rhotau);
    }

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v2sigma2, "v2sigma2", v2sigma2);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v2sigmalapl, "v2sigmalapl", v2sigmalapl);
      printVar(nspin, func.dim.v2sigmatau, "v2sigmatau", v2sigmatau);
      printVar(nspin, func.dim.v2lapl2, "v2lapl2", v2lapl2);
      printVar(nspin, func.dim.v2lapltau, "v2lapltau", v2lapltau);
      printVar(nspin, func.dim.v2tau2, "v2tau2", v2tau2);
     }
    printf("\n");
  }

  if(info->flags & XC_FLAGS_HAVE_KXC){
    printVar(nspin, func.dim.v3rho3, "v3rho3", v3rho3);

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v3rho2sigma, "v3rho2sigma", v3rho2sigma);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v3rho2lapl, "v3rho2lapl", v3rho2lapl);
      printVar(nspin, func.dim.v3rho2tau, "v3rho2tau", v3rho2tau);
    }

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v3rhosigma2, "v3rhosigma2", v3rhosigma2);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v3rhosigmalapl, "v3rhosigmalapl", v3rhosigmalapl);
      printVar(nspin, func.dim.v3rhosigmatau, "v3rhosigmatau", v3rhosigmatau);
      printVar(nspin, func.dim.v3rholapl2, "v3rholapl2", v3rholapl2);
      printVar(nspin, func.dim.v3rholapltau, "v3rholapltau", v3rholapltau);
      printVar(nspin, func.dim.v3rhotau2, "v3rhotau2", v3rhotau2);
    }

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v3sigma3, "v3sigma3", v3sigma3);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v3sigma2lapl, "v3sigma2lapl", v3sigma2lapl);
      printVar(nspin, func.dim.v3sigma2tau, "v3sigma2tau", v3sigma2tau);
      printVar(nspin, func.dim.v3sigmalapl2, "v3sigmalapl2", v3sigmalapl2);
      printVar(nspin, func.dim.v3sigmalapltau, "v3sigmalapltau", v3sigmalapltau);
      printVar(nspin, func.dim.v3sigmatau2, "v3sigmatau2", v3sigmatau2);
      printVar(nspin, func.dim.v3lapl3, "v3lapl3", v3lapl3);
      printVar(nspin, func.dim.v3lapl2tau, "v3lapl2tau", v3lapl2tau);
      printVar(nspin, func.dim.v3lapltau2, "v3lapltau2", v3lapltau2);
      printVar(nspin, func.dim.v3tau3, "v3tau3", v3tau3);
    }
    printf("\n");
  }

  if(info->flags & XC_FLAGS_HAVE_LXC){
    printVar(nspin, func.dim.v4rho4, "v4rho4", v4rho4);

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4rho3sigma, "v4rho3sigma", v4rho3sigma);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4rho3lapl, "v4rho3lapl", v4rho3lapl);
      printVar(nspin, func.dim.v4rho3tau, "v4rho3tau", v4rho3tau);
    }

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4rho2sigma2, "v4rho2sigma2", v4rho2sigma2);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4rho2sigmalapl, "v4rho2sigmalapl", v4rho2sigmalapl);
      printVar(nspin, func.dim.v4rho2sigmatau, "v4rho2sigmatau", v4rho2sigmatau);
      printVar(nspin, func.dim.v4rho2lapl2, "v4rho2lapl2", v4rho2lapl2);
      printVar(nspin, func.dim.v4rho2lapltau, "v4rho2lapltau", v4rho2lapltau);
      printVar(nspin, func.dim.v4rho2tau2, "v4rho2tau2", v4rho2tau2);
    }

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4rhosigma3, "v4rhosigma3", v4rhosigma3);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4rhosigma2lapl, "v4rhosigma2lapl", v4rhosigma2lapl);
      printVar(nspin, func.dim.v4rhosigma2tau, "v4rhosigma2tau", v4rhosigma2tau);
      printVar(nspin, func.dim.v4rhosigmalapl2, "v4rhosigmalapl2", v4rhosigmalapl2);
      printVar(nspin, func.dim.v4rhosigmalapltau, "v4rhosigmalapltau", v4rhosigmalapltau);
      printVar(nspin, func.dim.v4rhosigmatau2, "v4rhosigmatau2", v4rhosigmatau2);
      printVar(nspin, func.dim.v4rholapl3, "v4rholapl3", v4rholapl3);
      printVar(nspin, func.dim.v4rholapl2tau, "v4rholapl2tau", v4rholapl2tau);
      printVar(nspin, func.dim.v4rholapltau2, "v4rholapltau2", v4rholapltau2);
      printVar(nspin, func.dim.v4rhotau3, "v4rhotau3", v4rhotau3);
    }

    if(info->family == XC_FAMILY_GGA || info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4sigma4, "v4sigma4", v4sigma4);
    }

    if(info->family == XC_FAMILY_MGGA){
      printVar(nspin, func.dim.v4sigma3lapl, "v4sigma3lapl", v4sigma3lapl);
      printVar(nspin, func.dim.v4sigma3tau, "v4sigma3tau", v4sigma3tau);
      printVar(nspin, func.dim.v4sigma2lapl2, "v4sigma2lapl2", v4sigma2lapl2);
      printVar(nspin, func.dim.v4sigma2lapltau, "v4sigma2lapltau", v4sigma2lapltau);
      printVar(nspin, func.dim.v4sigma2tau2, "v4sigma2tau2", v4sigma2tau2);
      printVar(nspin, func.dim.v4sigmalapl3, "v4sigmalapl3", v4sigmalapl3);
      printVar(nspin, func.dim.v4sigmalapl2tau, "v4sigmalapl2tau", v4sigmalapl2tau);
      printVar(nspin, func.dim.v4sigmalapltau2, "v4sigmalapltau2", v4sigmalapltau2);
      printVar(nspin, func.dim.v4sigmatau3, "v4sigmatau3", v4sigmatau3);
      printVar(nspin, func.dim.v4lapl4, "v4lapl4", v4lapl4);
      printVar(nspin, func.dim.v4lapl3tau, "v4lapl3tau", v4lapl3tau);
      printVar(nspin, func.dim.v4lapl2tau2, "v4lapl2tau2", v4lapl2tau2);
      printVar(nspin, func.dim.v4lapltau3, "v4lapltau3", v4lapltau3);
      printVar(nspin, func.dim.v4tau4, "v4tau4", v4tau4);
    }
  }

  xc_mgga_vars_free_all(zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
  xc_func_end(&func);

  return 0;
}
