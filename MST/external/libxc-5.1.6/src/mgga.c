/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_mgga.c"
#include "funcs_hyb_mgga.c"

void
xc_mgga(const xc_func_type *func, size_t np,
        const double *rho, const double *sigma, const double *lapl, const double *tau,
        double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  assert(func != NULL);
  const xc_dimensions *dim = &(func->dim);

  /* sanity check */
  if(zk != NULL && !(func->info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc\n",
            func->info->name);
    exit(1);
  }

  if(vrho != NULL && !(func->info->flags & XC_FLAGS_HAVE_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc\n",
            func->info->name);
    exit(1);
  }

  if(v2rho2 != NULL && !(func->info->flags & XC_FLAGS_HAVE_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc\n",
            func->info->name);
    exit(1);
  }

  if(v3rho3 != NULL && !(func->info->flags & XC_FLAGS_HAVE_KXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of kxc\n",
            func->info->name);
    exit(1);
  }

  if(v4rho4 != NULL && !(func->info->flags & XC_FLAGS_HAVE_LXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of lxc\n",
            func->info->name);
    exit(1);
  }

  /* initialize output to zero */
  if(zk != NULL)
    libxc_memset(zk, 0, dim->zk*np*sizeof(double));

  if(vrho != NULL){
    assert(vsigma != NULL);
    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
      assert(vlapl != NULL);
    assert(vtau   != NULL);

    libxc_memset(vrho,   0, dim->vrho  *np*sizeof(double));
    libxc_memset(vsigma, 0, dim->vsigma*np*sizeof(double));
    if(vlapl != NULL) /* This is important for deorbitalization */
      libxc_memset(vlapl,  0, dim->vlapl *np*sizeof(double));
    libxc_memset(vtau,   0, dim->vtau  *np*sizeof(double));
  }

  if(v2rho2 != NULL){
    assert(v2rhosigma != NULL);
    assert(v2rhotau   != NULL);
    assert(v2sigma2   != NULL);
    assert(v2sigmatau != NULL);
    assert(v2tau2     != NULL);

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      assert(v2rholapl   != NULL);
      assert(v2sigmalapl != NULL);
      assert(v2lapl2     != NULL);
      assert(v2lapltau   != NULL);
    }

    libxc_memset(v2rho2,     0, dim->v2rho2     *np*sizeof(double));
    libxc_memset(v2rhosigma, 0, dim->v2rhosigma *np*sizeof(double));
    libxc_memset(v2rhotau,   0, dim->v2rhotau   *np*sizeof(double));
    libxc_memset(v2sigma2,   0, dim->v2sigma2   *np*sizeof(double));
    libxc_memset(v2sigmatau, 0, dim->v2sigmatau *np*sizeof(double));
    libxc_memset(v2tau2,     0, dim->v2tau2     *np*sizeof(double));

    if(v2lapl2 != NULL) {  /* This is important for deorbitalization */
      libxc_memset(v2rholapl,   0, dim->v2rholapl  *np*sizeof(double));
      libxc_memset(v2sigmalapl, 0, dim->v2sigmalapl*np*sizeof(double));
      libxc_memset(v2lapl2,     0, dim->v2lapl2    *np*sizeof(double));
      libxc_memset(v2lapltau,   0, dim->v2lapltau  *np*sizeof(double));
    }
  }

  if(v3rho3 != NULL){
    assert(v3rho2sigma   != NULL);
    assert(v3rho2tau     != NULL);
    assert(v3rhosigma2   != NULL);
    assert(v3rhosigmatau != NULL);
    assert(v3rhotau2     != NULL);
    assert(v3sigma3      != NULL);
    assert(v3sigma2tau   != NULL);
    assert(v3sigmatau2   != NULL);
    assert(v3tau3        != NULL);

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      assert(v3rho2lapl     != NULL);
      assert(v3rhosigmalapl != NULL);
      assert(v3rholapl2     != NULL);
      assert(v3rholapltau   != NULL);
      assert(v3sigma2lapl   != NULL);
      assert(v3sigmalapl2   != NULL);
      assert(v3sigmalapltau != NULL);
      assert(v3lapl3        != NULL);
      assert(v3lapl2tau     != NULL);
    }

    libxc_memset(v3rho3,        0, dim->v3rho3       *np*sizeof(double));
    libxc_memset(v3rho2sigma,   0, dim->v3rho2sigma  *np*sizeof(double));
    libxc_memset(v3rho2tau,     0, dim->v3rho2tau    *np*sizeof(double));
    libxc_memset(v3rhosigma2,   0, dim->v3rhosigma2  *np*sizeof(double));
    libxc_memset(v3rhosigmatau, 0, dim->v3rhosigmatau*np*sizeof(double));
    libxc_memset(v3rhotau2,     0, dim->v3rhotau2    *np*sizeof(double));
    libxc_memset(v3sigma3,      0, dim->v3sigma3     *np*sizeof(double));
    libxc_memset(v3sigma2tau,   0, dim->v3sigma2tau  *np*sizeof(double));
    libxc_memset(v3sigmatau2,   0, dim->v3sigmatau2  *np*sizeof(double));
    libxc_memset(v3tau3,        0, dim->v3tau3       *np*sizeof(double));

    if(v3lapl3 != NULL) {  /* This is important for deorbitalization */
      libxc_memset(v3rho2lapl,     0, dim->v3rho2lapl    *np*sizeof(double));
      libxc_memset(v3rhosigmalapl, 0, dim->v3rhosigmalapl*np*sizeof(double));
      libxc_memset(v3rholapl2,     0, dim->v3rholapl2    *np*sizeof(double));
      libxc_memset(v3rholapltau,   0, dim->v3rholapltau  *np*sizeof(double));
      libxc_memset(v3sigma2lapl,   0, dim->v3sigma2lapl  *np*sizeof(double));
      libxc_memset(v3sigmalapl2,   0, dim->v3sigmalapl2  *np*sizeof(double));
      libxc_memset(v3sigmalapltau, 0, dim->v3sigmalapltau*np*sizeof(double));
      libxc_memset(v3lapl3,        0, dim->v3lapl3       *np*sizeof(double));
      libxc_memset(v3lapl2tau,     0, dim->v3lapl2tau    *np*sizeof(double));
    }
  }

  if(v4rho4 != NULL){
    assert(v4rho4         != NULL);
    assert(v4rho3sigma    != NULL);
    assert(v4rho3tau      != NULL);
    assert(v4rho2sigma2   != NULL);
    assert(v4rho2sigmatau != NULL);
    assert(v4rho2tau2     != NULL);
    assert(v4rhosigma3    != NULL);
    assert(v4rhosigma2tau != NULL);
    assert(v4rhosigmatau2 != NULL);
    assert(v4rhotau3      != NULL);
    assert(v4sigma4       != NULL);
    assert(v4sigma3tau    != NULL);
    assert(v4sigma2tau2   != NULL);
    assert(v4sigmatau3    != NULL);
    assert(v4tau4         != NULL);

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      assert(v4rho3lapl        != NULL);
      assert(v4rho2sigmalapl   != NULL);
      assert(v4rho2lapl2       != NULL);
      assert(v4rho2lapltau     != NULL);
      assert(v4rhosigma2lapl   != NULL);
      assert(v4rhosigmalapl2   != NULL);
      assert(v4rhosigmalapltau != NULL);
      assert(v4rholapl3        != NULL);
      assert(v4rholapl2tau     != NULL);
      assert(v4rholapltau2     != NULL);
      assert(v4sigma3lapl      != NULL);
      assert(v4sigma2lapl2     != NULL);
      assert(v4sigma2lapltau   != NULL);
      assert(v4sigmalapl3      != NULL);
      assert(v4sigmalapl2tau   != NULL);
      assert(v4sigmalapltau2   != NULL);
      assert(v4lapl4           != NULL);
      assert(v4lapl3tau        != NULL);
      assert(v4lapl2tau2       != NULL);
      assert(v4lapltau3        != NULL);
    }

    libxc_memset(v4rho4,         0, dim->v4rho4        *np*sizeof(double));
    libxc_memset(v4rho3sigma,    0, dim->v4rho3sigma   *np*sizeof(double));
    libxc_memset(v4rho3tau,      0, dim->v4rho3tau     *np*sizeof(double));
    libxc_memset(v4rho2sigma2,   0, dim->v4rho2sigma2  *np*sizeof(double));
    libxc_memset(v4rho2sigmatau, 0, dim->v4rho2sigmatau*np*sizeof(double));
    libxc_memset(v4rho2tau2,     0, dim->v4rho2tau2    *np*sizeof(double));
    libxc_memset(v4rhosigma3,    0, dim->v4rhosigma3   *np*sizeof(double));
    libxc_memset(v4rhosigma2tau, 0, dim->v4rhosigma2tau*np*sizeof(double));
    libxc_memset(v4rhosigmatau2, 0, dim->v4rhosigmatau2*np*sizeof(double));
    libxc_memset(v4rhotau3,      0, dim->v4rhotau3     *np*sizeof(double));
    libxc_memset(v4sigma4,       0, dim->v4sigma4      *np*sizeof(double));
    libxc_memset(v4sigma3tau,    0, dim->v4sigma3tau   *np*sizeof(double));
    libxc_memset(v4sigma2tau2,   0, dim->v4sigma2tau2  *np*sizeof(double));
    libxc_memset(v4sigmatau3,    0, dim->v4sigmatau3   *np*sizeof(double));
    libxc_memset(v4tau4,         0, dim->v4tau4        *np*sizeof(double));

    if(v4lapl4 != NULL) {  /* This is important for deorbitalization */
      libxc_memset(v4rho3lapl,        0, dim->v4rho3lapl       *np*sizeof(double));
      libxc_memset(v4rho2sigmalapl,   0, dim->v4rho2sigmalapl  *np*sizeof(double));
      libxc_memset(v4rho2lapl2,       0, dim->v4rho2lapl2      *np*sizeof(double));
      libxc_memset(v4rho2lapltau,     0, dim->v4rho2lapltau    *np*sizeof(double));
      libxc_memset(v4rhosigma2lapl,   0, dim->v4rhosigma2lapl  *np*sizeof(double));
      libxc_memset(v4rhosigmalapl2,   0, dim->v4rhosigmalapl2  *np*sizeof(double));
      libxc_memset(v4rhosigmalapltau, 0, dim->v4rhosigmalapltau*np*sizeof(double));
      libxc_memset(v4rholapl3,        0, dim->v4rholapl3       *np*sizeof(double));
      libxc_memset(v4rholapl2tau,     0, dim->v4rholapl2tau    *np*sizeof(double));
      libxc_memset(v4rholapltau2,     0, dim->v4rholapltau2    *np*sizeof(double));
      libxc_memset(v4sigma3lapl,      0, dim->v4sigma3lapl     *np*sizeof(double));
      libxc_memset(v4sigma2lapl2,     0, dim->v4sigma2lapl2    *np*sizeof(double));
      libxc_memset(v4sigma2lapltau,   0, dim->v4sigma2lapltau  *np*sizeof(double));
      libxc_memset(v4sigmalapl3,      0, dim->v4sigmalapl3     *np*sizeof(double));
      libxc_memset(v4sigmalapl2tau,   0, dim->v4sigmalapl2tau  *np*sizeof(double));
      libxc_memset(v4sigmalapltau2,   0, dim->v4sigmalapltau2  *np*sizeof(double));
      libxc_memset(v4lapl4,           0, dim->v4lapl4          *np*sizeof(double));
      libxc_memset(v4lapl3tau,        0, dim->v4lapl3tau       *np*sizeof(double));
      libxc_memset(v4lapl2tau2,       0, dim->v4lapl2tau2      *np*sizeof(double));
      libxc_memset(v4lapltau3,        0, dim->v4lapltau3       *np*sizeof(double));
    }
  }

  /* call functional */
  if(func->info->mgga != NULL)
    func->info->mgga(func, np, rho, sigma, lapl, tau,
                     zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));

  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, sigma, lapl, tau,
                zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
}

/* specializations */
void
xc_mgga_exc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, zk, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void
xc_mgga_exc_vxc(const xc_func_type *p, size_t np,
                const double *rho, const double *sigma, const double *lapl, const double *tau,
                double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void xc_mgga_exc_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2) {
  xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2,
          v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void xc_mgga_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2) {
  xc_mgga(p, np, rho, sigma, lapl, tau, NULL, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2,
          v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void xc_mgga_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                             const double *rho, const double *sigma, const double *lapl, const double *tau,
                             double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                             double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                             double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                             double *v2lapltau, double *v2tau2,
                             double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                             double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                             double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                             double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                             double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                             double *v3tau3) {
  xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2,
          v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,
          v3lapltau2, v3tau3,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void xc_mgga_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                             const double *rho, const double *sigma, const double *lapl, const double *tau,
                             double *vrho, double *vsigma, double *vlapl, double *vtau,
                             double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                             double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                             double *v2lapltau, double *v2tau2,
                             double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                             double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                             double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                             double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                             double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                             double *v3tau3) {
  xc_mgga(p, np, rho, sigma, lapl, tau, NULL, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2,
          v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,
          v3lapltau2, v3tau3,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}


void
xc_mgga_vxc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, NULL, vrho, vsigma, vlapl, vtau,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void
xc_mgga_fxc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
            double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
            double *v2lapltau, double *v2tau2)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, NULL, NULL, NULL, NULL, NULL,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2,
          v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void xc_mgga_kxc(const xc_func_type *p, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau,
                 double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                 double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                 double *v3rholapl2, double *v3rholapltau,  double *v3rhotau2,
                 double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
                 double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmatau2,
                 double *v3lapl3, double *v3lapl2tau, double *v3lapltau2, double *v3tau3)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,
          v3lapltau2, v3tau3,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void xc_mgga_lxc(const xc_func_type *p, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau,
                 double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
                 double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
                 double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
                 double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
                 double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
                 double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
                 double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
                 double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
                 double *v4lapl2tau2, double *v4lapltau3, double *v4tau4)
{
    xc_mgga(p, np, rho, sigma, lapl, tau, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2,
            v4rho2sigmalapl, v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau,
            v4rho2tau2, v4rhosigma3, v4rhosigma2lapl, v4rhosigma2tau,
            v4rhosigmalapl2, v4rhosigmalapltau, v4rhosigmatau2,
            v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3,
            v4sigma4, v4sigma3lapl, v4sigma3tau, v4sigma2lapl2,
            v4sigma2lapltau, v4sigma2tau2, v4sigmalapl3, v4sigmalapl2tau,
            v4sigmalapltau2, v4sigmatau3, v4lapl4, v4lapl3tau,
            v4lapl2tau2, v4lapltau3, v4tau4
            );
}
