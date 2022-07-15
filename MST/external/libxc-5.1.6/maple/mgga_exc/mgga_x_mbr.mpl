(*
 Copyright (C) 2020 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* prefix:
  mgga_x_mbr_params *params;

  assert(p->params != NULL);
  params = (mgga_x_mbr_params * ) (p->params);
*)

(* replace: "br89_x\(" -> "xc_mgga_x_br89_get_x(" *)

$include "mgga_x_br89.mpl"

params_a_at := 0:

(*Equation 15. The three equations below are from mgga_x_tm.mpl; note
that the numerical factors in the denominator in eqn 15 just
correspond to going from 'x' to 's'.*)
tm_p  := x -> (X2S*x)^2:
tm_y  := x -> (2*params_a_lambda - 1)^2 * tm_p(x):
tm_f0 := x -> (1 + 10*(70/27)*tm_y(x) + params_a_beta*tm_y(x)^2)^(1/10):

(* definition below equation 16 *)
mbr_D := (ts, xs) -> 2*ts - 1/4 * (2*params_a_lambda - 1)^2 * xs^2:

(* k_\sigma is not defined in the paper, but Subrata Jana
said on GitLab that k_\sigma = (6\pi^2\rho_\sigma)^{1/3} *)
k_sigma := (6*Pi^2)^(1/3):

(* Equation 18. Note that there's a typo in the paper: since Becke's
tau is two times the kinetic energy density, there should be a factor
of two in front of tau uniform as well. *)
br89_Q := (x, u, t) ->
  1/6*(
  + 6*(params_a_lambda^2 - params_a_lambda + 1/2)*(2*t - 2*K_FACTOR_C - 1/36*x^2)
  + 6/5*k_sigma^2*(tm_f0(x)^2 - 1)
  - 2*params_a_gamma*mbr_D(t, x)
  ):
