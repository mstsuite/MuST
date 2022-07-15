/*
 Copyright (C) 2006-2018 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_mgga.c
 * @brief This file is to be included in MGGA functionals.
 */

#ifdef XC_DEBUG
#define __USE_GNU
#include <fenv.h>
#endif

/* hack to avoid compiler warnings */
#define NOARG

#ifdef XC_NO_EXC
#define OUT_PARAMS MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, )
#else
#define OUT_PARAMS XC_COMMA zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, )
#endif

#ifdef HAVE_CUDA
__global__ static void
work_mgga_gpu(const XC(func_type) *p, int order, size_t np, const double *rho, const double *sigma, const double *lapl, const double *tau,
              double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));
#endif

/**
 * @param[in,out] func_type: pointer to functional structure
 */
static void
work_mgga(const XC(func_type) *p, size_t np,
          const double *rho, const double *sigma, const double *lapl, const double *tau,
          double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
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
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

#ifdef HAVE_CUDA

  //make a copy of 'p' since it might be in host-only memory
  XC(func_type) * pcuda = (XC(func_type) *) libxc_malloc(sizeof(XC(func_type)));

  *pcuda = *p;

  size_t nblocks = np/CUDA_BLOCK_SIZE;
  if(np != nblocks*CUDA_BLOCK_SIZE) nblocks++;

  work_mgga_gpu<<<nblocks, CUDA_BLOCK_SIZE>>>(pcuda, order, np, rho, sigma, lapl, tau,
                                              zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));

  libxc_free(pcuda);

#else
  size_t ip;
  double dens;
  double my_rho[2]={0.0, 0.0};
  double my_sigma[3]={0.0, 0.0, 0.0};
  double my_tau[2]={0.0, 0.0};

  for(ip = 0; ip < np; ip++){
    //fprintf(stderr,
    //        "%14.10le %14.10le %14.10le %14.10le %14.10le %14.10le %14.10le %14.10le %14.10le\n",
    //        rho[0], rho[1], sigma[0], sigma[1], sigma[2], lapl[0], lapl[1], tau[0], tau[1]);
    //fflush(stderr);
  
    /* Screen low densities */
    dens = (p->nspin == XC_POLARIZED) ? rho[0]+rho[1] : rho[0];
    if(dens >= p->dens_threshold) {
      /* sanity check of input parameters */
      my_rho[0] = m_max(p->dens_threshold, rho[0]);
      /* Many functionals shamelessly divide by tau, so we set a reasonable threshold */
      /* skip all checks on tau for the kinetic functionals */
      if(p->info->family != XC_KINETIC)
        my_tau[0] = m_max(p->tau_threshold, tau[0]);
      my_sigma[0] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[0]);
      /* The Fermi hole curvature 1 - xs^2/(8*ts) must be positive */
      if(p->info->family != XC_KINETIC)
        my_sigma[0] = m_min(my_sigma[0], 8.0*my_rho[0]*my_tau[0]);
      /* lapl can have any values */
      if(p->nspin == XC_POLARIZED){
        double s_ave;

        my_rho[1] = m_max(p->dens_threshold, rho[1]);
        if(p->info->family != XC_KINETIC)
          my_tau[1] = m_max(p->tau_threshold, tau[1]);
        my_sigma[2] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[2]);
        if(p->info->family != XC_KINETIC)
          my_sigma[2] = m_min(my_sigma[2], 8.0*my_rho[1]*my_tau[1]);

        my_sigma[1] = sigma[1];
        s_ave = 0.5*(my_sigma[0] + my_sigma[2]);
        /* | grad n |^2 = |grad n_up + grad n_down|^2 > 0 */
        my_sigma[1] = (my_sigma[1] >= -s_ave ? my_sigma[1] : -s_ave);
        /* Since |grad n_up - grad n_down|^2 > 0 we also have */
        my_sigma[1] = (my_sigma[1] <= +s_ave ? my_sigma[1] : +s_ave);
      }

      if(p->nspin == XC_UNPOLARIZED)
        func_unpol(p, order, my_rho, my_sigma, lapl, my_tau OUT_PARAMS);
      else if(p->nspin == XC_POLARIZED)
        func_pol  (p, order, my_rho, my_sigma, lapl, my_tau OUT_PARAMS);
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
        for(ii=0; ii < dim->vsigma; ii++)
          is_OK = is_OK && isfinite(vsigma[ii]);
        if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
          for(ii=0; ii < dim->vlapl; ii++)
            is_OK = is_OK && isfinite(vlapl[ii]);
        for(ii=0; ii < dim->vtau; ii++)
          is_OK = is_OK && isfinite(vtau[ii]);
      }

      if(!is_OK){
        printf("Problem in the evaluation of the functional\n");
        if(p->nspin == XC_UNPOLARIZED){
          printf("./xc-get_data %d 1 ", p->info->number);
          if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
            printf("%le 0.0 %le 0.0 0.0 %le 0.0 %le 0.0\n",
                   *rho, *sigma, *lapl, *tau);
          else
            printf("%le 0.0 %le 0.0 0.0 0.0 0.0 %le 0.0\n",
                  *rho, *sigma, *tau);
        }else{
          printf("./xc-get_data %d 2 ", p->info->number);
          if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
            printf("%le %le %le %le %le %le %le %le %le\n",
                   rho[0], rho[1], sigma[0], sigma[1], sigma[2], lapl[0], lapl[1], tau[0], tau[1]);
          else
            printf("%le %le %le %le %le 0.0 0.0 %le %le\n",
                   rho[0], rho[1], sigma[0], sigma[1], sigma[2], tau[0], tau[1]);
        }
      }
    }
#endif

    internal_counters_mgga_next(&(p->dim), 0, &rho, &sigma, &lapl, &tau,
                                &zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));
  }   /* for(ip) */

#endif

}

#ifdef HAVE_CUDA
__global__ static void
work_mgga_gpu(const XC(func_type) *p, int order, size_t np,
              const double *rho, const double *sigma, const double *lapl, const double *tau,
              double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{

  size_t ip = blockIdx.x * blockDim.x + threadIdx.x;
  double my_rho[2] = {0.0, 0.0};
  double my_sigma[3] = {0.0, 0.0, 0.0};
  double my_tau[2] = {0.0, 0.0};
  double dens;

  if(ip >= np) return;

  internal_counters_mgga_random(&(p->dim), ip, 0, &rho, &sigma, &lapl, &tau,
                                &zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));

  /* Screen low densities */
  dens = (p->nspin == XC_POLARIZED) ? rho[0]+rho[1] : rho[0];
  if(dens >= p->dens_threshold) {
    /* sanity check of input parameters */
    my_rho[0] = m_max(p->dens_threshold, rho[0]);
    /* Many functionals shamelessly divide by tau, so we set a reasonable threshold */
    if(p->info->family != XC_KINETIC)
      my_tau[0] = m_max(p->tau_threshold, tau[0]);
    /* The Fermi hole curvature 1 - xs^2/(8*ts) must be positive */
    my_sigma[0] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[0]);
    if(p->info->family != XC_KINETIC)
      my_sigma[0] = m_min(my_sigma[0], 8.0*my_rho[0]*my_tau[0]);
    /* lapl can have any values */
    if(p->nspin == XC_POLARIZED){
      double s_ave;

      my_rho[1]   = m_max(p->dens_threshold, rho[1]);
      if(p->info->family != XC_KINETIC)
        my_tau[1]   = m_max(p->tau_threshold, tau[1]);
      my_sigma[2] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[2]);
      if(p->info->family != XC_KINETIC)
        my_sigma[2] = m_min(my_sigma[2], 8.0*my_rho[1]*my_tau[1]);

      my_sigma[1] = sigma[1];
      s_ave = 0.5*(my_sigma[0] + my_sigma[2]);
      /* | grad n |^2 = |grad n_up + grad n_down|^2 > 0 */
      my_sigma[1] = (my_sigma[1] >= -s_ave ? my_sigma[1] : -s_ave);
      /* Since |grad n_up - grad n_down|^2 > 0 we also have */
      my_sigma[1] = (my_sigma[1] <= +s_ave ? my_sigma[1] : +s_ave);
    }

    if(p->nspin == XC_UNPOLARIZED)
      func_unpol(p, order, my_rho, my_sigma, lapl, my_tau OUT_PARAMS);
    else if(p->nspin == XC_POLARIZED)
      func_pol  (p, order, my_rho, my_sigma, lapl, my_tau OUT_PARAMS);
  }
}
#endif
