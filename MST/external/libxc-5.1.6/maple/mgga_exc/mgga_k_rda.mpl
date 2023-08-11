(*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_k_rda_params *params;

  assert(p->params != NULL);
  params = (mgga_k_rda_params * ) (p->params);
*)

rda_s := x -> X2S*x:
rda_p := u -> X2S^2*u:

(* Equation (61) *)
rda_k4 := (s, p, b) -> sqrt(s^4 + b*p^2):
(* Equation (63) *)
rda_k2 := (s, p, b) -> s^2 + b*p:

(* Equation (71); first term is von Weiszäcker according to equation (13) *)
rda_f0 := (s, p) ->
       5/3*s^2 + params_a_A0
       + params_a_A1 * (rda_k4(s,p,params_a_a) / (1 + params_a_beta1*rda_k4(s,p,params_a_a)))^2
       + params_a_A2 * (rda_k4(s,p,params_a_b) / (1 + params_a_beta2*rda_k4(s,p,params_a_b)))^4
       + params_a_A3 * (rda_k2(s,p,params_a_c) / (1 + params_a_beta3*rda_k2(s,p,params_a_c))):

(* Complete functional *)
rda_f := (xs, us) -> rda_f0(rda_s(xs), rda_p(us)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(rda_f, rs, z, xs0, xs1, u0, u1):
