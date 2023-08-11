(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_optc_params *params;

  assert(p->params != NULL);
  params = (gga_c_optc_params * )(p->params);
*)

$include "gga_c_pw91.mpl"

optc_f2 := (rs, z, xt, xs0, xs1) ->
  + f_pw91(rs*(2/(1 + z))^(1/3),  1, xs0, xs0, 0)*opz_pow_n( z,1)/2
  + f_pw91(rs*(2/(1 - z))^(1/3), -1, xs1, 0, xs1)*opz_pow_n(-z,1)/2:

f  := (rs, z, xt, xs0, xs1) ->
  + params_a_c1*f_pw91(rs, z, xt, xs0, xs1) + (params_a_c2 - params_a_c1)*optc_f2(rs, z, xt, xs0, xs1):
