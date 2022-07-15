(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_k_pc07_params *params;

  assert(p->params != NULL);
  params = (mgga_k_pc07_params * ) (p->params);
*)

pc07_p := x -> X2S^2*x^2:
pc07_q := u -> X2S^2*u:

(* Equation (15) redefined with decaying exponentials to avoid inf/inf situations *)
pc07_fab0 := z -> exp(-params_a_a*params_a_b/z) * (1+exp(-params_a_a/(params_a_a-z)))^params_a_b/(exp(-params_a_a/z) + exp(-params_a_a/(params_a_a-z)))^params_a_b:
(* The function is ill-behaving at z=0 and z=a.
   However, it also goes to 0 and 1 very quickly near these points.
   pc07_fab0(params_a_a/40) is ~ 1.6e-17, which is smaller than machine epsilon.
*)
pc07_thr := 1/40:
pc07_zlo := pc07_thr*params_a_a:
pc07_zhi := (1-pc07_thr)*params_a_a:
pc07_fab := z -> my_piecewise5(z<=pc07_zlo, 0, z>=pc07_zhi, 1, pc07_fab0( m_min(pc07_zhi, m_max(pc07_zlo, z)) ) ):

(* Equation (7) *)
pc07_Delta := (x, u) ->
  8*pc07_q(u)^2/81 - pc07_p(x)*pc07_q(u)/9 + 8*pc07_p(x)^2/243:

pc07_f_W    := x -> 5*pc07_p(x)/3:

(* Equation (8) *)
pc07_GE4  := (x, u) ->
  1 + 5*pc07_p(x)/27 + 20*pc07_q(u)/9 + pc07_Delta(x, u):

(* Equation (11) *)
pc07_GE4_M := (x, u) ->
  pc07_GE4(x, u)/sqrt(1 + pc07_Delta(x, u)^2/(1 + pc07_f_W(x))^2):

(* Equation (17) *)
pc07_f := (x, u) ->
  pc07_f_W(x) + (pc07_GE4_M(x, u) - pc07_f_W(x))*pc07_fab(pc07_GE4_M(x, u) - pc07_f_W(x)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_kinetic(pc07_f, rs, z, xs0, xs1, u0, u1):
