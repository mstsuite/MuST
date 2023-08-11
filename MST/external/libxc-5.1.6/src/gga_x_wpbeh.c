/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_WPBEH 524 /* short-range version of the PBE */

/* This implementation follows the one from espresso, that, in turn,
   follows the one of the thesis of Jochen Heyd. Analytic derivatives
   are only implemented in espresso though. These implementations can
   be found in:

   vasp: xclib_grad.F, MODULE wpbe, and in particular SUBROUTINE EXCHWPBE_R
   espresso: flib/functionals.f90, SUBROUTINE wpbe_analy_erfc_approx_grad

   very important details can be found in references:

   *) J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
      Erratum: J. Chem. Phys. 124, 219906 (2006).
   *) M Ernzerhof and JP Perdew, J. Chem. Phys. 109, 3313 (1998)
   *) J Heyd and GE Scuseria, J. Chem. Phys. 120, 7274 (2004)

   Also the whole mess with the rescaling of s is explained in

   *) TM Henderson, AF Izmaylov, G Scalmani, and GE Scuseria,
      J. Chem. Phys. 131, 044108 (2009)
*/

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_wpbeh.c"
#include "work_gga.c"

static void
xc_gga_x_wpbeh_init(xc_func_type *p)
{
  xc_hyb_init_hybrid(p, 0.0);
}

static const char  *omega_names[]  = {"_omega"};
static const char  *omega_desc[]   = {"screening parameter"};
/* The default value is actually PBEh */
static const double omega_values[] = {0.0};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_wpbeh = {
  XC_GGA_X_WPBEH,
  XC_EXCHANGE,
  "short-range part of the PBE (default w=0 gives PBEh)",
  XC_FAMILY_GGA,
  {&xc_ref_Heyd2003_8207, &xc_ref_Heyd2003_8207_err, &xc_ref_Ernzerhof1998_3313, &xc_ref_Heyd2004_7274, &xc_ref_Henderson2009_044108},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_omega},
  xc_gga_x_wpbeh_init, NULL,
  NULL, work_gga, NULL
};
