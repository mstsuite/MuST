(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_x_kt_params *params;

  assert(p->params != NULL);
  params = (gga_x_kt_params * )(p->params);
*)

kt_fx := (rs, z, xs) ->
   1 - params_a_gamma/X_FACTOR_C*n_spin(rs, z)^(4/3)*xs^2/(n_spin(rs, z)^(4/3) + params_a_delta):

(* we want energy per particle *)
f := (rs, zeta, xt, xs0, xs1) -> gga_exchange_nsp(kt_fx, rs, zeta, xs0, xs1):
