/*
 Copyright (C) 2006-2018 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_lda.c
 * @brief This file is to be included in LDA functionals.
 */

#ifdef XC_DEBUG
#define __USE_GNU
#include <fenv.h>
#endif

/* hack to avoid compiler warnings */
#define NOARG

#ifdef XC_NO_EXC
#define OUT_PARAMS LDA_OUT_PARAMS_NO_EXC(XC_COMMA, )
#else
#define OUT_PARAMS XC_COMMA zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, )
#endif

#ifdef HAVE_CUDA
__global__ static void
work_lda_gpu(const XC(func_type) *p, int order, size_t np, const double *rho,
             double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));
#endif

/**
 * @param[in,out] func_type: pointer to functional structure
 */

static void
work_lda(const XC(func_type) *p, size_t np, const double *rho,
         double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

#ifdef XC_DEBUG
  /* This throws an exception when floating point errors are encountered */
  /*feenableexcept(FE_DIVBYZERO | FE_INVALID);*/
#endif

#ifdef HAVE_CUDA

  //make a copy of 'p' since it might be in host-only memory
  XC(func_type) * pcuda = (XC(func_type) *) libxc_malloc(sizeof(XC(func_type)));

  *pcuda = *p;

  size_t nblocks = np/CUDA_BLOCK_SIZE;
  if(np != nblocks*CUDA_BLOCK_SIZE) nblocks++;

  work_lda_gpu<<<nblocks, CUDA_BLOCK_SIZE>>>(pcuda, order, np, rho,
                                             zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, ));

  libxc_free(pcuda);

#else
  size_t ip;
  double dens;
  double my_rho[2] = {0.0, 0.0};

  for(ip = 0; ip < np; ip++){
    /* Screen low density */
    dens = (p->nspin == XC_POLARIZED) ? rho[0]+rho[1] : rho[0];
    if(dens >= p->dens_threshold) {
      /* sanity check of input parameters */
      my_rho[0] = m_max(p->dens_threshold, rho[0]);
      if(p->nspin == XC_POLARIZED){
        my_rho[1] = m_max(p->dens_threshold, rho[1]);
      }
      if(p->nspin == XC_UNPOLARIZED)
        func_unpol(p, order, my_rho OUT_PARAMS);
      else if(p->nspin == XC_POLARIZED)
        func_pol  (p, order, my_rho OUT_PARAMS);
    }

    /* check for NaNs */
#ifdef XC_DEBUG
    {
      size_t ii;
      const xc_dimensions *dim = &(p->dim);
      int is_OK = 1;

      if(zk != NULL)
        is_OK = is_OK & isfinite(*zk);

      if(vrho != NULL){
        for(ii=0; ii < dim->vrho; ii++)
          is_OK = is_OK && isfinite(vrho[ii]);
      }

      if(!is_OK){
        printf("Problem in the evaluation of the functional\n");
        if(p->nspin == XC_UNPOLARIZED){
          printf("./xc-get_data %d 1 %le 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
                 p->info->number, *rho);
        }else{
          printf("./xc-get_data %d 2 %le %le 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
                 p->info->number, rho[0], rho[1]);
        }
      }
    }
#endif

    internal_counters_lda_next(&(p->dim), 0, &rho, &zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));
  }   /* for(ip) */

#endif

}

#ifdef HAVE_CUDA

__global__ static void
work_lda_gpu(const XC(func_type) *p, int order, size_t np, const double *rho,
             double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  double my_rho[2] = {0.0, 0.0};
  double dens;
  size_t ip = blockIdx.x*blockDim.x + threadIdx.x;

  if(ip >= np) return;

  internal_counters_lda_random(&(p->dim), ip, 0, &rho, &zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));

  dens = (p->nspin == XC_POLARIZED) ? rho[0]+rho[1] : rho[0];
  if(dens >= p->dens_threshold) {
    my_rho[0] = m_max(p->dens_threshold, rho[0]);
    if(p->nspin == XC_POLARIZED){
      my_rho[1] = m_max(p->dens_threshold, rho[1]);
    }
    if(p->nspin == XC_UNPOLARIZED)
      func_unpol(p, order, my_rho OUT_PARAMS);
    else if(p->nspin == XC_POLARIZED)
      func_pol  (p, order, my_rho OUT_PARAMS);
  }
}

#endif
