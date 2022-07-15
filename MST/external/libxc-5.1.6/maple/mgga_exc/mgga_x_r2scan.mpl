(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_r2scan_params *params;

  assert(p->params != NULL);
  params = (mgga_x_r2scan_params * )(p->params);
*)

$include "mgga_x_rscan.mpl"
$include "mgga_x_scan.mpl"

(* eqn S6 *)
r2scan_alpha := (x, t) -> (t - x^2/8)/(K_FACTOR_C + params_a_eta*x^2/8):

(* f(alpha) replaced with a polynomial for alpha in [0, 2.5], eqn S7 *)
r2scan_f_alpha_neg := a -> exp(-params_a_c1*a/(1 - a)):
r2scan_f_alpha := (a, ff) -> my_piecewise5(a <= 0, r2scan_f_alpha_neg(m_min(a, 0)), a <= 2.5, rscan_f_alpha_small(m_min(a, 2.5), ff), rscan_f_alpha_large(m_max(a, 2.5))):

(* eqn S11 *)
Cn := 20/27 + params_a_eta*5/3:
(* eqn S12 *)
C2 := ff -> -add(i*ff[9-i], i=1..8) * (1-scan_h0x):

(* eqn S10; this is analogous to scan_y *)
r2scan_x := (p, ff) -> (Cn*C2(ff)*exp(-p^2/params_a_dp2^4)+MU_GE)*p:

r2scan_f := (x, u, t) -> (scan_h1x(r2scan_x(scan_p(x), rscan_fx)) + r2scan_f_alpha(r2scan_alpha(x, t), rscan_fx) * (scan_h0x - scan_h1x(r2scan_x(scan_p(x), rscan_fx))))*scan_gx(x):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(r2scan_f, rs, z, xs0, xs1, u0, u1, t0, t1):
