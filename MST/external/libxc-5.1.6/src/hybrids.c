/*
 Copyright (C) 2020 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* some specializations */
void
xc_hyb_init_hybrid(xc_func_type *p, double alpha)
{
  p->cam_alpha = alpha;
  p->cam_beta = 0.0;
  p->cam_omega = 0.0;
}

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

void
xc_hyb_init_sr(xc_func_type *p, double beta, double omega)
{
  assert(is_rangesep(p));
  p->cam_alpha = 0.0;
  p->cam_beta = beta;
  p->cam_omega = omega;
}

void
xc_hyb_init_cam(xc_func_type *p, double alpha, double beta, double omega)
{
  assert(is_cam(p));
  p->cam_alpha = alpha;
  p->cam_beta = beta;
  p->cam_omega = omega;
}

void
xc_hyb_init_camy(xc_func_type *p, double alpha, double beta, double omega)
{
  assert(is_yukawa(p));
  p->cam_alpha = alpha;
  p->cam_beta = beta;
  p->cam_omega = omega;
}

/*------------------------------------------------------*/
/* returns the mixing coefficient for the hybrid functions */
double
xc_hyb_exx_coef(const xc_func_type *p)
{
  assert(p!=NULL);
  return p->cam_alpha;
}

/* returns the CAM parameters for screened hybrids */
void
xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta)
{
  assert(p!=NULL);
  *omega = p->cam_omega;
  *alpha = p->cam_alpha;
  *beta = p->cam_beta;
}
