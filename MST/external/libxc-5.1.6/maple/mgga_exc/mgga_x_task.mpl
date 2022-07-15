(*
 Copyright (C) 2019 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_task_params *params;

  assert(p->params != NULL);
  params = (mgga_x_task_params * )(p->params);
*)

task_alpha := (x, t) -> (t/K_FACTOR_C) * m_max(1 - x^2/(8*t), 1e-10):

task_gx0 := x -> 1 - exp(-params_a_task_c*x^(-1/4)):
task_gx := x -> my_piecewise3(x > 0, task_gx0(m_max(x, 0)), 0):

task_hx1 := r -> simplify(add(params_a_task_anu[i+1]*ChebyshevT(i, (r - 1)/(r + 1)), i=0..2)):

task_fx  := r -> simplify(add(params_a_task_bnu[i+1]*ChebyshevT(i, (r - 1)/(r + 1)), i=0..4)):

task_f0 := (s, a) -> params_a_task_h0x*task_gx(s^2) +
  (1.0 - task_fx(a))*(task_hx1(s^2) - params_a_task_h0x)*task_gx(s^2)^params_a_task_d:

task_f := (x, u, t) -> task_f0(X2S*x, task_alpha(x, t)):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(task_f, rs, z, xs0, xs1, u0, u1, t0, t1):
