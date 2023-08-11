(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_tpss_params *params;

  assert(p->params != NULL);
  params = (mgga_x_tpss_params * ) (p->params);
*)

$ifdef mgga_x_revtpss_params
params_a_b      := 0.40:
params_a_c      := 2.35203946:
params_a_e      := 2.16769874:
params_a_kappa  := 0.804:
params_a_mu     := 0.14:
params_a_BLOC_a := 3.0:
params_a_BLOC_b := 0.0:
$endif

tpss_ff     := z -> params_a_BLOC_a + params_a_BLOC_b*z:
tpss_kappa := (x, t) -> params_a_kappa:

$include "tpss_x.mpl"

(* Equation (5) *)

tpss_a1  := (x, t) -> tpss_kappa(x, t)/(tpss_kappa(x, t) + tpss_fx(x, t)):
tpss_f   := (x, u, t) -> 1 + tpss_kappa(x, t)*(1 - tpss_a1(x, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(tpss_f, rs, z, xs0, xs1, u0, u1, t0, t1):
