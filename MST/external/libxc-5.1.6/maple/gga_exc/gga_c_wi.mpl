(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_wi_params *params;

  assert(p->params != NULL);
  params = (gga_c_wi_params * )(p->params);
*)

f_num := xt -> params_a_a + params_a_b*xt^2*exp(-params_a_k*xt^2):
f_den := (rs, xt) -> params_a_c + rs*(1 + params_a_d*(4*Pi/3)^(1/3)*xt^(7/2)):

f := (rs, zeta, xt, xs0, xs1) -> f_num(xt)/f_den(rs, xt):
