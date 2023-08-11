(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_k_ol2_params *params;

  assert(p->params != NULL);
  params = (gga_k_ol2_params * )(p->params);
*)

ol2_f := x ->
  + params_a_aa
  + params_a_bb*x^2/72.0
  + params_a_cc*x/(2^(1/3) + 4*x):

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(ol2_f, rs, zeta, xs0, xs1):
