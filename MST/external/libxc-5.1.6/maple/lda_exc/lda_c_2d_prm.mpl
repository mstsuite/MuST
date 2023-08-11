(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_2d_prm_params *params;

  assert(p->params != NULL);
  params = (lda_c_2d_prm_params * )(p->params);

  assert(params->N > 1);
*)

$define xc_dimensions_2d

prm_q := 3.9274: (* 2.258 *)

(* Equation (4) *)
beta := rs -> prm_q/(sqrt(Pi)*rs):

phi  := rs -> beta(rs)/(beta(rs) + sqrt(Pi)/2):

f0 := rs ->
  + sqrt(Pi)*beta(rs)*(phi(rs) - 1)/(2*sqrt(2 + params_a_c)) (* original version has (phi-1)^2 *)
  + phi(rs)*(phi(rs) - 1)/(2 + params_a_c)
  + sqrt(Pi)*phi(rs)*phi(rs)/(4*beta(rs)*(2 + params_a_c)^1.5)
  + sqrt(Pi)*beta(rs)*(phi(rs) - 1)/sqrt(1 + params_a_c)
  + phi(rs)/(1 + params_a_c):

f  := (rs, z) -> f0(rs)*Pi/(2*prm_q*prm_q):
