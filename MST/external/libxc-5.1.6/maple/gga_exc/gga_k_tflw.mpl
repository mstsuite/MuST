(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_k_tflw_params *params;

  assert(p->params != NULL);
  params = (gga_k_tflw_params * )(p->params);
*)

tflw_f := x -> params_a_gamma + (params_a_lambda/8)*x^2/K_FACTOR_C:

f := (rs, zeta, xt, xs0, xs1) -> gga_kinetic(tflw_f, rs, zeta, xs0, xs1):
