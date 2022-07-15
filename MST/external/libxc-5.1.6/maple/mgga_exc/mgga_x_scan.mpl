(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_scan_params *params;

  assert(p->params != NULL);
  params = (mgga_x_scan_params * )(p->params);
*)

scan_p     := x -> X2S^2*x^2:
scan_alpha := (x, t) -> (t - x^2/8)/K_FACTOR_C:

(* The interpolating functions are nasty for a -> 1, so we need to
   truncate them. The natural choice is to cut off the functions to
   zero when the exponential term reaches machine epsilon.

   The left cutoff is |log epsilon|/(|log epsilon| + c1) < 1
   and the right one is (|log epsilon| + c2)/|log epsilon| > 1,
   so we don't even really need the step function.
*)
scan_f_alpha_left0 := a -> exp(-params_a_c1*a/(1 - a)):
scan_f_alpha_left_cutoff := -log(DBL_EPSILON)/(-log(DBL_EPSILON) + params_a_c1):
scan_f_alpha_left := a -> my_piecewise3(a > scan_f_alpha_left_cutoff, 0, scan_f_alpha_left0(m_min(scan_f_alpha_left_cutoff, a))):

scan_f_alpha_right0 := a -> -params_a_d*exp(params_a_c2/(1 - a)):
scan_f_alpha_right_cutoff := (-log(DBL_EPSILON/abs(params_a_d)) + params_a_c2)/(-log(DBL_EPSILON/abs(params_a_d))):
scan_f_alpha_right := a -> my_piecewise3(a < scan_f_alpha_right_cutoff, 0, scan_f_alpha_right0(m_max(scan_f_alpha_right_cutoff, a))):
scan_f_alpha := a -> my_piecewise3(
  a <= 1, scan_f_alpha_left(a), scan_f_alpha_right(a)
  ):

scan_h1x := x -> 1 + params_a_k1*(1 - params_a_k1/(params_a_k1 + x)):

scan_b2 := sqrt(5913/405000):
scan_b1 := (511/13500)/(2*scan_b2):
scan_b3 := 1/2:
scan_b4 := MU_GE^2/params_a_k1 - 1606/18225 - scan_b1^2:
scan_y  := (x, a) -> MU_GE*scan_p(x) + scan_b4*scan_p(x)^2*exp(-scan_b4*scan_p(x)/MU_GE)
  + (scan_b1*scan_p(x) + scan_b2*(1 - a)*exp(-scan_b3*(1 - a)^2))^2:

scan_a1 := 4.9479:
scan_gx := x -> 1 - exp(-scan_a1/sqrt(X2S*x)):

scan_h0x := 1.174:
scan_f   := (x, u, t) -> (scan_h1x(scan_y(x, scan_alpha(x, t)))*(1 - scan_f_alpha(scan_alpha(x, t)))
  + scan_h0x*scan_f_alpha(scan_alpha(x, t)))*scan_gx(x):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) -> mgga_exchange(scan_f, rs, z, xs0, xs1, u0, u1, t0, t1):
