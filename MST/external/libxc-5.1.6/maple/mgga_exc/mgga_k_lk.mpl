(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_k_lk_params *params;

  assert(p->params != NULL);
  params = (mgga_k_lk_params * ) (p->params);
*)

lk_p := x -> X2S^2*x^2:
lk_q := u -> X2S^2*u:

(* Equation (10) *)
lk_delta := (p, q) -> 8/81*q^2 - 1/9*p*q + 8/243*p^2:
(* Equation (15) *)
lk_f0 := (x1, x2) -> 1 + params_a_kappa*(2 - 1/(1+x1/params_a_kappa) - 1/(1+x2/params_a_kappa)):
(* Equation (16) *)
lk_x1 := (p, q) -> 5/27*p + lk_delta(p,q) + (5/27*p)^2/params_a_kappa:
(* Equation (17) *)
lk_x2 := (p, q) -> 2*(5/27*p)*lk_delta(p,q)/params_a_kappa + (5/27*p)^3/params_a_kappa^2:

(* Full functional *)
lk_f := (x, u) ->
  lk_f0(lk_x1(lk_p(x),lk_q(u)), lk_x2(lk_p(x),lk_q(u))):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(lk_f, rs, z, xs0, xs1, u0, u1):
