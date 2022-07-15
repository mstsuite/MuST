/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola
                2019 X. Andrade

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || (id) == XC_FAMILY_HYB_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA || (id) == XC_FAMILY_HYB_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA || (id) == XC_FAMILY_HYB_LDA ||  is_gga(id))
#define safe_free(pt) if(pt != NULL) libxc_free(pt)

void
xc_mgga_vars_allocate_all(int family, size_t np, const xc_dimensions *dim,
                     int do_zk, int do_vrho, int do_v2rho2, int do_v3rho3, int do_v4rho4,
                     double **zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double **, ))
{
  /* allocate buffers */
  if(do_zk)
    *zk = (double *) libxc_malloc(sizeof(double)*np*dim->zk);

#ifndef XC_DONT_COMPILE_VXC
  if(do_vrho){
    *vrho = (double *) libxc_malloc(sizeof(double)*np*dim->vrho);
    if(is_gga(family)){
      *vsigma = (double *) libxc_malloc(sizeof(double)*np*dim->vsigma);
    }
    if(is_mgga(family)){
      *vlapl = (double *) libxc_malloc(sizeof(double)*np*dim->vlapl);
      *vtau  = (double *) libxc_malloc(sizeof(double)*np*dim->vtau);
    }
  }

#ifndef XC_DONT_COMPILE_FXC
  if(do_v2rho2){
    *v2rho2 = (double *) libxc_malloc(sizeof(double)*np*dim->v2rho2);
    if(is_gga(family)){
      *v2rhosigma  = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhosigma);
      *v2sigma2    = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigma2);
    }
    if(is_mgga(family)){
      *v2rholapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v2rholapl);
      *v2rhotau    = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhotau);
      *v2sigmalapl = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmalapl);
      *v2sigmatau  = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmatau);
      *v2lapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapl2);
      *v2lapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapltau);
      *v2tau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v2tau2);
    }
  }

#ifndef XC_DONT_COMPILE_KXC
  if(do_v3rho3){
    *v3rho3      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho3);
    if(is_gga(family)){
      *v3rho2sigma = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2sigma);
      *v3rhosigma2 = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigma2);
      *v3sigma3    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma3);
    }
    if(is_mgga(family)){
      *v3rho2lapl     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2lapl);
      *v3rho2tau      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2tau);
      *v3rhosigmalapl = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmalapl);
      *v3rhosigmatau  = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmatau);
      *v3rholapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapl2);
      *v3rholapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapltau);
      *v3rhotau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhotau2);
      *v3sigma2lapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2lapl);
      *v3sigma2tau    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2tau);
      *v3sigmalapl2   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapl2);
      *v3sigmalapltau = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapltau);
      *v3sigmatau2    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmatau2);
      *v3lapl3        = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl3);
      *v3lapl2tau     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl2tau);
      *v3lapltau2     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapltau2);
      *v3tau3         = (double *) libxc_malloc(sizeof(double)*np*dim->v3tau3);
    }
  }

#ifndef XC_DONT_COMPILE_LXC
  if(do_v4rho4){
    *v4rho4            = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho4);
    if(is_gga(family)){
      *v4rho3sigma       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3sigma);
      *v4rho2sigma2      = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigma2);
      *v4rhosigma3       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma3);
      *v4sigma4          = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma4);
    }
    if(is_mgga(family)){
      *v4rho3lapl        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3lapl);
      *v4rho3tau         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3tau);
      *v4rho2sigmalapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmalapl);
      *v4rho2sigmatau    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmatau);
      *v4rho2lapl2       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapl2);
      *v4rho2lapltau     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapltau);
      *v4rho2tau2        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2tau2);
      *v4rhosigma2lapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2lapl);
      *v4rhosigma2tau    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2tau);
      *v4rhosigmalapl2   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapl2);
      *v4rhosigmalapltau = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapltau);
      *v4rhosigmatau2    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmatau2);
      *v4rholapl3        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl3);
      *v4rholapl2tau     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl2tau);
      *v4rholapltau2     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapltau2);
      *v4rhotau3         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhotau3);
      *v4sigma3lapl      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3lapl);
      *v4sigma3tau       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3tau);
      *v4sigma2lapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapl2);
      *v4sigma2lapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapltau);
      *v4sigma2tau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2tau2);
      *v4sigmalapl3      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl3);
      *v4sigmalapl2tau   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl2tau);
      *v4sigmalapltau2   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapltau2);
      *v4sigmatau3       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmatau3);
      *v4lapl4           = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl4);
      *v4lapl3tau        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl3tau);
      *v4lapl2tau2       = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl2tau2);
      *v4lapltau3        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapltau3);
      *v4tau4            = (double *) libxc_malloc(sizeof(double)*np*dim->v4tau4);
    }
  }
#endif
#endif
#endif
#endif
}

void
xc_mgga_vars_free_all(double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  /* deallocate internal buffers */
  safe_free(zk);
#ifndef XC_DONT_COMPILE_VXC
  safe_free(vrho);
  safe_free(vsigma);
  safe_free(vlapl);
  safe_free(vtau);

#ifndef XC_DONT_COMPILE_FXC
  safe_free(v2rho2); safe_free(v2rhosigma); safe_free(v2rholapl); safe_free(v2rhotau);
  safe_free(v2sigma2); safe_free(v2sigmalapl); safe_free(v2sigmatau);
  safe_free(v2lapl2); safe_free(v2lapltau); safe_free(v2tau2);

#ifndef XC_DONT_COMPILE_KXC
  safe_free(v3rho3); safe_free(v3rho2sigma); safe_free(v3rho2lapl); safe_free(v3rho2tau);
  safe_free(v3rhosigma2); safe_free(v3rhosigmalapl); safe_free(v3rhosigmatau);
  safe_free(v3rholapl2); safe_free(v3rholapltau); safe_free(v3rhotau2);
  safe_free(v3sigma3); safe_free(v3sigma2lapl); safe_free(v3sigma2tau);
  safe_free(v3sigmalapl2); safe_free(v3sigmalapltau); safe_free(v3sigmatau2);
  safe_free(v3lapl3); safe_free(v3lapl2tau); safe_free(v3lapltau2); safe_free(v3tau3);

#ifndef XC_DONT_COMPILE_LXC
  safe_free(v4rho4); safe_free(v4rho3sigma); safe_free(v4rho3lapl); safe_free(v4rho3tau);
  safe_free(v4rho2sigma2); safe_free(v4rho2sigmalapl); safe_free(v4rho2sigmatau);
  safe_free(v4rho2lapl2); safe_free(v4rho2lapltau); safe_free(v4rho2tau2);
  safe_free(v4rhosigma3); safe_free(v4rhosigma2lapl); safe_free(v4rhosigma2tau);
  safe_free(v4rhosigmalapl2); safe_free(v4rhosigmalapltau); safe_free(v4rhosigmatau2);
  safe_free(v4rholapl3); safe_free(v4rholapl2tau); safe_free(v4rholapltau2); safe_free(v4rhotau3);
  safe_free(v4sigma4); safe_free(v4sigma3lapl); safe_free(v4sigma3tau); safe_free(v4sigma2lapl2);
  safe_free(v4sigma2lapltau); safe_free(v4sigma2tau2); safe_free(v4sigmalapl3); safe_free(v4sigmalapl2tau);
  safe_free(v4sigmalapltau2); safe_free(v4sigmatau3); safe_free(v4lapl4); safe_free(v4lapl3tau);
  safe_free(v4lapl2tau2); safe_free(v4lapltau3); safe_free(v4tau4);
#endif
#endif
#endif
#endif
}

void
xc_mgga_evaluate_functional(const xc_func_type *func, size_t np,
                            const double *rho, const double *sigma, const double *lapl, const double *tau,
                            double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  double *mzk = NULL;
  if(func->info->flags & XC_FLAGS_HAVE_EXC)
    mzk = zk;

  /* Evaluate the functional */
  switch(func->info->family){
  case XC_FAMILY_LDA:
    xc_lda (func, np, rho,
            mzk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  case XC_FAMILY_GGA:
    xc_gga (func, np, rho, sigma,
            mzk GGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  case XC_FAMILY_MGGA:
    xc_mgga(func, np, rho, sigma, lapl, tau,
            mzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  }
}

/* initializes the mixing */
void
xc_deorbitalize_init(xc_func_type *p, int mgga_id, int ked_id)
{
  assert(p != NULL && p->func_aux == NULL);

  /* allocate structures needed for */
  p->n_func_aux = 2;
  p->func_aux   = (xc_func_type **) libxc_malloc(2*sizeof(xc_func_type *));

  p->func_aux[0] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));
  p->func_aux[1] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));

  xc_func_init (p->func_aux[0], mgga_id, p->nspin);
  xc_func_init (p->func_aux[1], ked_id,  p->nspin);
}

void
xc_deorbitalize_func(const xc_func_type *func, size_t np,
                     const double *rho, const double *sigma, const double *lapl, const double *tau,
                     double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  double *mrho, *msigma, *mlapl, *mtau;
  const double *null = NULL;
  double *mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, mgga_);
  double *ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, ked1_);
  double *ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, ked2_);
  size_t ii;
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  /* prepare buffers that will hold the results from the individual functionals */
  mgga_zk MGGA_OUT_PARAMS_NO_EXC(=, mgga_ ) = NULL;
  ked1_zk MGGA_OUT_PARAMS_NO_EXC(=, ked1_ ) = NULL;
  ked2_zk MGGA_OUT_PARAMS_NO_EXC(=, ked2_ ) = NULL;

  /* allocate buffers */
  xc_mgga_vars_allocate_all(func->func_aux[0]->info->family, np, &(func->func_aux[0]->dim),
                       order >= 0, order >= 1, order >= 2, order >= 3, order >= 4,
                       &mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, mgga_));
  xc_mgga_vars_allocate_all(func->func_aux[1]->info->family, np, &(func->func_aux[1]->dim),
                       order >= 0, order >= 1, order >= 2, order >= 3, order >= 4,
                       &ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked1_));
  if(func->nspin == XC_UNPOLARIZED){
    mtau   = (double *) libxc_malloc(sizeof(double)*np);
  }else{
    mrho   = (double *) libxc_malloc(2*sizeof(double)*np);
    msigma = (double *) libxc_malloc(3*sizeof(double)*np);
    mlapl  = (double *) libxc_malloc(2*sizeof(double)*np);
    mtau   = (double *) libxc_malloc(2*sizeof(double)*np);

    xc_mgga_vars_allocate_all(func->func_aux[1]->info->family, np, &(func->func_aux[1]->dim),
                         order >= 0, order >= 1, order >= 2, order >= 3, order >= 4,
                         &ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked2_));
  }

  /* evaluate the kinetic energy functional */
  if(func->nspin == XC_UNPOLARIZED){
    xc_mgga_evaluate_functional(func->func_aux[1], np, rho, sigma, lapl, tau,
                           ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked1_));
  }else{
    for(ii=0; ii<np; ii++){
      mrho  [2*ii] = rho  [2*ii]; mrho  [2*ii+1] = 0.0;
      msigma[3*ii] = sigma[3*ii]; msigma[3*ii+1] = 0.0; msigma[3*ii+2] = 0.0;
      mlapl [2*ii] = lapl [2*ii]; mlapl [2*ii+1] = 0.0;
      mtau  [2*ii] = tau  [2*ii]; mtau  [2*ii+1] = 0.0;
    }
    xc_mgga_evaluate_functional(func->func_aux[1], np, mrho, msigma, mlapl, mtau,
                           ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked1_));

    for(ii=0; ii<np; ii++){
      mrho  [2*ii] = rho  [2*ii + 1];
      msigma[3*ii] = sigma[3*ii + 2];
      mlapl [2*ii] = lapl [2*ii + 1];
      mtau  [2*ii] = tau  [2*ii + 1];
    }
    xc_mgga_evaluate_functional(func->func_aux[1], np, mrho, msigma, mlapl, mtau,
                           ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked2_));
  }

  /* now evaluate the mgga functional */
  if(func->nspin == XC_UNPOLARIZED){
    for(ii=0; ii<np; ii++){
      mtau[ii] = rho[ii]*ked1_zk[ii];
    }
  }else{
    for(ii=0; ii<np; ii++){
      mtau[2*ii    ] = rho[2*ii    ]*ked1_zk[ii];
      mtau[2*ii + 1] = rho[2*ii + 1]*ked2_zk[ii];
    }
  }
  xc_mgga_evaluate_functional(func->func_aux[0], np, rho, sigma, lapl, mtau,
                         mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, mgga_));

  /* now we have to combine the results */
  for(ii=0; ii<np; ii++){
    if(zk != NULL){
      *zk = *mgga_zk;
    }

#ifndef XC_DONT_COMPILE_VXC
    if(vrho != NULL){
#include "maple2c/deorbitalize_1.c"
    }
#ifndef XC_DONT_COMPILE_FXC
    if(v2rho2 != NULL){
#include "maple2c/deorbitalize_2.c"
    }
#ifndef XC_DONT_COMPILE_KXC
    if(v3rho3 != NULL){
#include "maple2c/deorbitalize_3.c"
    }
#ifndef XC_DONT_COMPILE_LXC
    if(v4rho4 != NULL){
#include "maple2c/deorbitalize_4.c"
    }
#endif
#endif
#endif
#endif

    internal_counters_mgga_next(&(func->dim), 0, &null, &null, &null, &null,
                                &zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));
    internal_counters_mgga_next(&(func->func_aux[0]->dim), 0, &null, &null, &null, &null,
                                &mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, mgga_));
    internal_counters_mgga_next(&(func->func_aux[1]->dim), 0, &null, &null, &null, &null,
                                &ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked1_));
    if(func->nspin == XC_POLARIZED){
      internal_counters_mgga_next(&(func->func_aux[1]->dim), 0, &null, &null, &null, &null,
                                  &ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked2_));
    }
  }

  /* move the counters back to zero and deallocate the memory */
  internal_counters_mgga_random(&(func->func_aux[0]->dim), -np, 0, &null, &null, &null, &null,
                                &mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, mgga_));
  xc_mgga_vars_free_all(mgga_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, mgga_));

  internal_counters_mgga_random(&(func->func_aux[1]->dim), -np, 0, &null, &null, &null, &null,
                                &ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked1_));
  xc_mgga_vars_free_all(ked1_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked1_));
  if(func->nspin == XC_POLARIZED){
    internal_counters_mgga_random(&(func->func_aux[1]->dim), -np, 0, &null, &null, &null, &null,
                                  &ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ked2_));
    xc_mgga_vars_free_all(ked2_zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ked2_));

    free(mrho); free(msigma); free(mlapl);
  }
  free(mtau);
}

