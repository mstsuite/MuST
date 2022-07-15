/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_lda.c"
#include "funcs_hyb_lda.c"

/* get the lda functional */
void
xc_lda(const xc_func_type *func, size_t np, const double *rho,
       double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
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

  /* initialize output */
  if(zk != NULL)
    libxc_memset(zk,     0, np*sizeof(double)*dim->zk);

  if(vrho != NULL)
    libxc_memset(vrho,   0, np*sizeof(double)*dim->vrho);

  if(v2rho2 != NULL)
    libxc_memset(v2rho2, 0, np*sizeof(double)*dim->v2rho2);

  if(v3rho3 != NULL)
    libxc_memset(v3rho3, 0, np*sizeof(double)*dim->v3rho3);

  if(v4rho4 != NULL)
    libxc_memset(v4rho4, 0, np*sizeof(double)*dim->v4rho4);

  /* call the LDA routines */
  if(func->info->lda != NULL)
    func->info->lda(func, np, rho, zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, ));

  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, NULL, NULL, NULL, zk, vrho, NULL, NULL, NULL,
                v2rho2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                v3rho3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                v4rho4, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL);
}


/* specializations */
void
xc_lda_exc(const xc_func_type *p, size_t np, const double *rho, double *zk)
{
  xc_lda(p, np, rho, zk, NULL, NULL, NULL, NULL);
}

void
xc_lda_exc_vxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho)
{
  xc_lda(p, np, rho, zk, vrho, NULL, NULL, NULL);
}

void
xc_lda_exc_vxc_fxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2)
{
  xc_lda(p, np, rho, zk, vrho, v2rho2, NULL, NULL);
}

void
xc_lda_vxc_fxc(const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2)
{
  xc_lda(p, np, rho, NULL, vrho, v2rho2, NULL, NULL);
}

void
xc_lda_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_lda(p, np, rho, zk, vrho, v2rho2, v3rho3, NULL);
}

void
xc_lda_vxc_fxc_kxc(const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_lda(p, np, rho, NULL, vrho, v2rho2, v3rho3, NULL);
}

void
xc_lda_vxc(const xc_func_type *p, size_t np, const double *rho, double *vrho)
{
  xc_lda(p, np, rho, NULL, vrho, NULL, NULL, NULL);
}

void
xc_lda_fxc(const xc_func_type *p, size_t np, const double *rho, double *v2rho2)
{
  xc_lda(p, np, rho, NULL, NULL, v2rho2, NULL, NULL);
}

void
xc_lda_kxc(const xc_func_type *p, size_t np, const double *rho, double *v3rho3)
{
  xc_lda(p, np, rho, NULL, NULL, NULL, v3rho3, NULL);
}

void
xc_lda_lxc(const xc_func_type *p, size_t np, const double *rho, double *v4rho4)
{
  xc_lda(p, np, rho, NULL, NULL, NULL, NULL, v4rho4);
}
