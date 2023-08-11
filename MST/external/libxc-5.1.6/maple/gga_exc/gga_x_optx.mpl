(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_optx_params *params;

  assert(p->params != NULL);
  params = (gga_x_optx_params * )(p->params);
*)

optx_f := x-> params_a_a + params_a_b*(params_a_gamma*x^2/(1 + params_a_gamma*x^2))^2:

f := (rs, z, xt, xs0, xs1) -> gga_exchange(optx_f, rs, z, xs0, xs1):
