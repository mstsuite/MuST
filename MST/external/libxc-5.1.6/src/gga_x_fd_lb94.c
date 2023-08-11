/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_FD_LB94     604 /* Functional derivative recovered from the stray LB94 potential */
#define XC_GGA_X_FD_REVLB94  605 /* Revised FD_LB94 */

typedef struct{
  double beta;         /* screening parameter beta */
} gga_x_fd_lb94_params;

#define N_PAR 1
static const char *names[N_PAR] = {"_beta"};
static const char *desc[N_PAR] = {"beta parameter"};

static const double lb94_par[N_PAR] = {0.05};
static const double revlb94_par[N_PAR] = {0.004};

static void
gga_x_fd_lb94_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_fd_lb94_params));
}

GPU_FUNCTION
static inline double FT_inter(int n, double x)
{
  static double fd_beta = 0.05, fd_csi = M_CBRT2;

  double mlog;

  mlog = (n == 0) ? 1 : log(x);

  return -3/4 * fd_beta*fd_csi*mlog /
    (1 + 3*fd_beta*fd_csi*x*log(fd_csi*x + sqrt(fd_csi*fd_csi*x*x + 1)));
}

GPU_FUNCTION
static void func0(double *x, int n, void *dummy)
{
  int ii;

  for(ii=0; ii<n; ii++)
    x[ii] = FT_inter(0, x[ii]);
}

GPU_FUNCTION
static void func1(double *x, int n, void *dummy)
{
  int ii;

  for(ii=0; ii<n; ii++)
    x[ii] = FT_inter(1, x[ii]);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_fd_lb94.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_fd_lb94 = {
  XC_GGA_X_FD_LB94,
  XC_EXCHANGE,
  "Functional derivative recovered from the stray LB94 potential",
  XC_FAMILY_GGA,
  {&xc_ref_Gaiduk2011_012509, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, lb94_par, set_ext_params_cpy},
  gga_x_fd_lb94_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_fd_revlb94 = {
  XC_GGA_X_FD_REVLB94,
  XC_EXCHANGE,
  "Revised FD_LB94",
  XC_FAMILY_GGA,
  {&xc_ref_Gaiduk2011_012509, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, revlb94_par, set_ext_params_cpy},
  gga_x_fd_lb94_init, NULL,
  NULL, work_gga, NULL
};
