(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

params_a_mu := MU_GE:
params_a_b  := 0.4:
params_a_c  := 2.14951:
params_a_e  := 1.987:

vt84_gamma := 0.000023:
tpss_ff    := z -> 3:
tpss_kappa := (x, t) ->
  1/(vt84_gamma/params_a_mu^2 + vt84_gamma/params_a_mu + 1):

$include "tpss_x.mpl"

(* Equation (8) *)

vt84_f   := (x, u, t) -> 1
    + tpss_fx(x, t)*exp(-vt84_gamma*tpss_fx(x, t)/params_a_mu)/(1 + tpss_fx(x, t))
    + (1 - exp(-vt84_gamma*tpss_fx(x, t)^2/params_a_mu^2))*(params_a_mu/tpss_fx(x, t) - 1):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(vt84_f, rs, z, xs0, xs1, u0, u1, t0, t1):
