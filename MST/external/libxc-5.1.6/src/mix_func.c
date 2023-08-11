/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola
                2019 X. Andrade

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/* initializes the mixing */
void
xc_mix_init(xc_func_type *p, int n_funcs, const int *funcs_id, const double *mix_coef)
{
  int ii;

  assert(p != NULL);
  assert(p->func_aux == NULL && p->mix_coef == NULL);

  /* allocate structures needed for mixed functional */
  p->n_func_aux = n_funcs;
  p->mix_coef   = (double *) libxc_malloc(n_funcs*sizeof(double));
  p->func_aux   = (xc_func_type **) libxc_malloc(n_funcs*sizeof(xc_func_type *));

  for(ii=0; ii<n_funcs; ii++){
    p->mix_coef[ii] = mix_coef[ii];
    p->func_aux[ii] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));
    xc_func_init (p->func_aux[ii], funcs_id[ii], p->nspin);
  }

  /* initialize variables */
  p->cam_alpha = 0.0;
  p->cam_beta  = 0.0;
  p->cam_omega = 0.0;

  p->nlc_b     = 0.0;
  p->nlc_C     = 0.0;
}

#ifdef HAVE_CUDA
__global__ static void add_to_mix_gpu(size_t np, double * dst, double coeff, const double *src){
  size_t ip = blockIdx.x * blockDim.x + threadIdx.x;
  if(ip < np) dst[ip] += coeff*src[ip];
}
#endif

static void add_to_mix(size_t np, double * dst, double coeff, const double *src){
#ifndef HAVE_CUDA
  size_t ip;
  for(ip = 0; ip < np; ip++) dst[ip] += coeff*src[ip];
#else
  size_t nblocks = np/CUDA_BLOCK_SIZE;
  if(np != nblocks*CUDA_BLOCK_SIZE) nblocks++;
  add_to_mix_gpu<<<nblocks, CUDA_BLOCK_SIZE>>>(np, dst, coeff, src);
#endif
}

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || (id) == XC_FAMILY_HYB_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA  || (id) == XC_FAMILY_HYB_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA  || (id) == XC_FAMILY_HYB_LDA ||  is_gga(id))
#define safe_free(pt) if(pt != NULL) libxc_free(pt)
#define sum_var(VAR) add_to_mix(np*dim->VAR, VAR, func->mix_coef[ii], x ## VAR);

void
xc_mix_func(const xc_func_type *func, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  const xc_func_type *aux;
  double *xzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA *, x);
  int ii;

  const xc_dimensions *dim = &(func->dim);

  /* Sanity check: have we claimed the highest possible derivatives?
     First, check for the lowest common derivative (also need to make
     sure the derivatives have been compiled in!)
  */
  int have_vxc = XC_FLAGS_I_HAVE_VXC;
  int have_fxc = XC_FLAGS_I_HAVE_FXC;
  int have_kxc = XC_FLAGS_I_HAVE_KXC;
  int have_lxc = XC_FLAGS_I_HAVE_LXC;
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    if(! (aux->info->flags & XC_FLAGS_HAVE_VXC))
      have_vxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_FXC))
      have_fxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_KXC))
      have_kxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_LXC))
      have_lxc = 0;
  }
  /* Then, for the actual checks */
  assert(have_lxc == (func->info->flags & XC_FLAGS_I_HAVE_LXC));
  assert(have_kxc == (func->info->flags & XC_FLAGS_I_HAVE_KXC));
  assert(have_fxc == (func->info->flags & XC_FLAGS_I_HAVE_FXC));
  assert(have_vxc == (func->info->flags & XC_FLAGS_I_HAVE_VXC));

  /* Sanity check: if component needs the Laplacian, then the mix
     must require it too */
  int need_laplacian = 0;
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
      need_laplacian = XC_FLAGS_NEEDS_LAPLACIAN;
  }
  assert((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) == need_laplacian);

  /* Check compatibility of the individual components */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    /* Sanity check: if component is GGA or meta-GGA, mix functional
       must also be GGA or meta-GGA */
    if(is_gga(aux->info->family))
      assert(is_gga(func->info->family));
    if(is_mgga(aux->info->family) && !is_mgga(func->info->family))
      assert(is_mgga(func->info->family));
    /* Sanity checks: if mix functional has higher derivatives, these
       must also exist in the individual components */
    if(func->info->flags & XC_FLAGS_HAVE_VXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_VXC);
    if(func->info->flags & XC_FLAGS_HAVE_FXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_FXC);
    if(func->info->flags & XC_FLAGS_HAVE_KXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_KXC);
    if(func->info->flags & XC_FLAGS_HAVE_LXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_LXC);
  }

  /* prepare buffers that will hold the results from the individual functionals */
  xzk MGGA_OUT_PARAMS_NO_EXC(=, x) = NULL;

  /* allocate buffers */
  xc_mgga_vars_allocate_all(func->info->family, np, dim,
                            zk != NULL, vrho != NULL, v2rho2 != NULL, v3rho3 != NULL, v4rho4 != NULL,
                            &xzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, x));

  /* Proceed by computing the mix */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];

    /* Evaluate the functional */
    switch(aux->info->family){
    case XC_FAMILY_LDA:
    case XC_FAMILY_HYB_LDA:
      xc_lda(aux, np, rho,
             xzk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, x));
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga(aux, np, rho, sigma,
             xzk GGA_OUT_PARAMS_NO_EXC(XC_COMMA, x));
      break;
    case XC_FAMILY_MGGA:
    case XC_FAMILY_HYB_MGGA:
      xc_mgga(aux, np, rho, sigma, lapl, tau,
              xzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, x));
      break;
    }

    /* Do the mixing */
    if(zk != NULL) {
      sum_var(zk);
    }

 #ifndef XC_DONT_COMPILE_VXC
    if(vrho != NULL) {
      sum_var(vrho);

      if(is_gga(aux->info->family)) {
        sum_var(vsigma);
      }

      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(vlapl);
        }
        sum_var(vtau);
      }
    }

#ifndef XC_DONT_COMPILE_FXC
    if(v2rho2 != NULL){
      sum_var(v2rho2);

      if(is_gga(aux->info->family)) {
        sum_var(v2rhosigma);
        sum_var(v2sigma2);
      }

      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(v2rholapl);
          sum_var(v2sigmalapl);
          sum_var(v2lapl2);
          sum_var(v2lapltau);
        }
        sum_var(v2rhotau);
        sum_var(v2sigmatau);
        sum_var(v2tau2);
      }
    }

#ifndef XC_DONT_COMPILE_KXC
    if(v3rho3 != NULL){
      sum_var(v3rho3);

      if(is_gga(aux->info->family)) {
        sum_var(v3rho2sigma);
        sum_var(v3rhosigma2);
        sum_var(v3sigma3);
      }

      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(v3rho2lapl);
          sum_var(v3rhosigmalapl);
          sum_var(v3rholapl2);
          sum_var(v3rholapltau);
          sum_var(v3sigma2lapl);
          sum_var(v3sigmalapl2);
          sum_var(v3sigmalapltau);
          sum_var(v3lapl3);
          sum_var(v3lapl2tau);
          sum_var(v3lapltau2);
        }
        sum_var(v3rho2tau);
        sum_var(v3rhosigmatau);
        sum_var(v3rhotau2);
        sum_var(v3sigma2tau);
        sum_var(v3sigmatau2);
        sum_var(v3tau3);
      }
    }

#ifndef XC_DONT_COMPILE_LXC
    if(v4rho4 != NULL){
      sum_var(v4rho4);

      if(is_gga(aux->info->family)) {
        sum_var(v4rho3sigma);
        sum_var(v4rho2sigma2);
        sum_var(v4rhosigma3);
        sum_var(v4sigma4);
      }
      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(v4rho3lapl);
          sum_var(v4rho2sigmalapl);
          sum_var(v4rho2lapl2);
          sum_var(v4rho2lapltau);
          sum_var(v4rhosigma2lapl);
          sum_var(v4rhosigmalapl2);
          sum_var(v4rhosigmalapltau);
          sum_var(v4rholapl3);
          sum_var(v4rholapl2tau);
          sum_var(v4rholapltau2);
          sum_var(v4sigma3lapl);
          sum_var(v4sigma2lapl2);
          sum_var(v4sigma2lapltau);
          sum_var(v4sigmalapl3);
          sum_var(v4sigmalapl2tau);
          sum_var(v4sigmalapltau2);
          sum_var(v4lapl4);
          sum_var(v4lapl3tau);
          sum_var(v4lapl2tau2);
          sum_var(v4lapltau3);
        }
        sum_var(v4rho3tau);
        sum_var(v4rho2sigmatau);
        sum_var(v4rho2tau2);
        sum_var(v4rhosigma2tau);
        sum_var(v4rhosigmatau2);
        sum_var(v4rhotau3);
        sum_var(v4sigma3tau);
        sum_var(v4sigma2tau2);
        sum_var(v4sigmatau3);
        sum_var(v4tau4);
      }
    }
#endif
#endif
#endif
#endif
  } /* end functional loop */

  /* deallocate internal buffers */
  xc_mgga_vars_free_all(xzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, x));
}

