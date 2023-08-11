(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_dk87_params *params;

  assert(p->params != NULL);
  params = (gga_x_dk87_params * )(p->params);
*)

dk87_betag := 7/(432*Pi*(6*Pi^2)^(1/3))/X_FACTOR_C:

dk87_f := x -> 1 + dk87_betag*x^2*(1 + params_a_a1*x^params_a_alpha)/(1 + params_a_b1*x^2):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(dk87_f, rs, z, xs0, xs1):
