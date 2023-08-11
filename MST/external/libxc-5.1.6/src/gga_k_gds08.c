/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_GDS08     591 /* Combined analytical theory with Monte Carlo sampling */
#define XC_GGA_K_GHDS10    592 /* As GDS08 but for an electron gas with spin */
#define XC_GGA_K_GHDS10R   593 /* Reparametrized GHDS10 */
#define XC_GGA_K_TKVLN     594 /* Trickey, Karasiev, and Vela */


static void
gga_k_gds08_init(xc_func_type *p)
{
  static int    funcs_id  [2] = {XC_GGA_K_VW, XC_LDA_K_GDS08_WORKER};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_gds08 = {
  XC_GGA_K_GDS08,
  XC_KINETIC,
  "Combined analytical theory with Monte Carlo sampling",
  XC_FAMILY_GGA,
  {&xc_ref_Ghiringhelli2008_073104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_gds08_init, NULL,
  NULL, NULL, NULL
};

static void
gga_k_ghds10_init(xc_func_type *p)
{
  static int    funcs_id  [2] = {XC_GGA_K_TFVW, XC_LDA_K_GDS08_WORKER};
  static double funcs_coef[2] = {1.0, 1.0};

  static double par_k_gds08[] = {1.02, 0.163, 0.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[1], par_k_gds08);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_ghds10 = {
  XC_GGA_K_GHDS10,
  XC_KINETIC,
  "As GDS08 but for an electron gas with spin",
  XC_FAMILY_GGA,
  {&xc_ref_Ghiringhelli2010_014106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_ghds10_init, NULL,
  NULL, NULL, NULL
};

static void
gga_k_ghds10r_init(xc_func_type *p)
{
  static int    funcs_id  [2] = {XC_GGA_K_TFVW, XC_LDA_K_GDS08_WORKER};
  static double funcs_coef[2] = {1.0, 1.0};

  static double par_k_gds08[] = {0.61434e-1, 0.61317e-2, 0.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[1], par_k_gds08);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_ghds10r = {
  XC_GGA_K_GHDS10R,
  XC_KINETIC,
  "Reparametrized GHDS10",
  XC_FAMILY_GGA,
  {&xc_ref_Trickey2011_075146, &xc_ref_Ghiringhelli2010_014106, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_ghds10r_init, NULL,
  NULL, NULL, NULL
};

static void
gga_k_tkvln_init(xc_func_type *p)
{
  static int    funcs_id  [2] = {XC_GGA_K_TFVW, XC_LDA_K_GDS08_WORKER};
  static double funcs_coef[2] = {1.0, 1.0};

  static double par_k_gds08[] = {0.45960e-1, 0.65545e-2, 0.23131e-3};

  xc_mix_init(p, 2, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[1], par_k_gds08);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_k_tkvln = {
  XC_GGA_K_TKVLN,
  XC_KINETIC,
  "Trickey, Karasiev, and Vela",
  XC_FAMILY_GGA,
  {&xc_ref_Trickey2011_075146, &xc_ref_Ghiringhelli2010_014106, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  gga_k_tkvln_init, NULL,
  NULL, NULL, NULL
};
