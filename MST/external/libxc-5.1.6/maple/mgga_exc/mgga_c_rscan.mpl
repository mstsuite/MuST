(*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "mgga_c_scan.mpl"

(* Coefficients of the rSCAN switching function from SI, in reversed(!) order: 7, 6, ..., 0 *)
rscan_fc := [-0.051848879792, 0.516884468372, -1.915710236206, 3.061560252175, -1.535685604549, -0.4352, -0.64, 1]:

np53 := rs -> n_total(rs)^(5/3):
(* First regularization: tau^u -> tau^u + tau^r *)

rscan_alpha0 := (rs, z, xt, ts0, ts1) ->
  (np53(rs)*m_max(t_total(z, ts0, ts1) - xt^2/8, 0))/((K_FACTOR_C*np53(rs) + 2^(2/3)*params_a_taur)*t_total(z, 1, 1)):

(* Second regularization: alpha -> alpha^3/(alpha^2 + alpha_r) *)
rscan_alpha := (rs, z, xt, ts0, ts1) -> rscan_alpha0(rs, z, xt, ts0, ts1)^3/(rscan_alpha0(rs, z, xt, ts0, ts1)^2 + params_a_alphar):

(* f(alpha) replaced with a polynomial for alpha in [0, 2.5] *)
rscan_f_alpha_small := (a,ff) -> add(ff[8-i]*a^i, i=0..7):
rscan_f_alpha_large := a -> -params_a_d*exp(params_a_c2/(1 - a)):
rscan_f_alpha := (a, ff) -> my_piecewise3( a <= 2.5, rscan_f_alpha_small(m_min(a, 2.5),ff), rscan_f_alpha_large(m_max(a, 2.5)) ):

(* set parameters of f_alpha *)
params_a_alphar := 1e-3:
params_a_taur := 1e-4:

rscan_f := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  f_pbe(rs, z, xt, xs0, xs1) + rscan_f_alpha(rscan_alpha(rs, z, xt, ts0, ts1), rscan_fc)*(
    + scan_e0(rs, z, X2S*2^(1/3)*xt)
    - f_pbe(rs, z, xt, xs0, xs1)
  ):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  rscan_f(rs, z, xt, xs0, xs1, ts0, ts1):
