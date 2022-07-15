(*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_rscan_params *params;

  assert(p->params != NULL);
  params = (mgga_x_rscan_params * )(p->params);
*)

$include "mgga_x_scan.mpl"

(* Coefficients of the rSCAN switching function from SI, in reversed(!) order: 7, 6, ..., 0 *)
rscan_fx := [-0.023185843322, 0.234528941479, -0.887998041597, 1.451297044490, -0.663086601049, -0.4445555, -0.667, 1]:

np53 := (rs, z) -> n_spin(rs,z)^(5/3):

(* First regularization: tau^u -> tau^u + tau_r.
   tau_r gets scaled by 2^(2/3) due to the spin scaling. *)
rscan_alpha0 := (rs, z, x, t) -> (np53(rs,z)*m_max(t - x^2/8, 0))/(np53(rs,z)*K_FACTOR_C + params_a_taur/2):

(* Second regularization: alpha -> alpha^3/(alpha^2 + alpha_r) *)
rscan_alpha := (rs, z, x, t) -> rscan_alpha0(rs, z, x, t)^3/(rscan_alpha0(rs, z, x ,t)^2 + params_a_alphar):

(* f(alpha) replaced with a polynomial for alpha in [0, 2.5] *)
rscan_f_alpha_small := (a, ff) -> add(ff[8-i]*a^i, i=0..7):
rscan_f_alpha_large := a -> -params_a_d*exp(params_a_c2/(1 - a)):
rscan_f_alpha := (a, ff) -> my_piecewise3(a <= 2.5, rscan_f_alpha_small(m_min(a, 2.5), ff), rscan_f_alpha_large(m_max(a, 2.5))):
rscan_f   := (rs, z, x, u, t) -> (scan_h1x(scan_y(x, rscan_alpha(rs, z, x, t)))*(1 - rscan_f_alpha(rscan_alpha(rs, z, x, t), rscan_fx))
  + scan_h0x*rscan_f_alpha(rscan_alpha(rs, z, x, t), rscan_fx))*scan_gx(x):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange_nsp(rscan_f, rs, z, xs0, xs1, u0, u1, t0, t1):
