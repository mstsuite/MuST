(*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_s12_params *params;

  assert(p->params != NULL);
  params = (gga_x_s12_params * )(p->params);
*)

s12g_f := x -> params_a_bx*(params_a_A + params_a_B*(1 - 1/(1 + params_a_C*x^2 + params_a_D*x^4))*(1 - 1/(1 + params_a_E*x^2))):

f := (rs, z, xt, xs0, xs1) -> gga_exchange(s12g_f, rs, z, xs0, xs1):

