(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

sa_a := 2.413:
sa_b := 0.348:

params_a_b      := 0.40:
params_a_c      := 1.59096:
params_a_e      := 1.537:
params_a_mu     := 0.21951:

sa_alpha           := (x, t) -> (t - x^2/8)/K_FACTOR_C:

(* Equation (8) *)

tpss_ff    := z -> 2:
tpss_kappa := (x, t) -> 2*Pi/(3*sqrt(5)) * \
          sqrt(sa_alpha(x, t) + 1)/sqrt(sa_a + log(sa_alpha(x, t) + sa_b)):

$include "tpss_x.mpl"

sa_a1  := (x, t) -> tpss_kappa(x, t)/(tpss_kappa(x, t) + tpss_fx(x, t)):
sa_f   := (x, u, t) -> 1 + tpss_kappa(x, t)*(1 - sa_a1(x, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(sa_f, rs, z, xs0, xs1, u0, u1, t0, t1):
