(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_zpbeint_params *params;

  assert(p->params != NULL);
  params = (gga_c_zpbeint_params * )(p->params);
*)

params_a_gamma := (1 - log(2))/Pi^2:
params_a_BB    := 1:
$include "gga_c_pbe.mpl"

ff := (z, t) -> mphi(z)^(params_a_alpha*t^3):

f  := (rs, z, xt, xs0, xs1) ->
  f_pw(rs, z) + ff(z, tp(rs, z, xt))*fH(rs, z, tp(rs, z, xt)):
