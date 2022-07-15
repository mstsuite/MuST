(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_bmk_params *params;

  assert(p->params != NULL);
  params = (gga_c_bmk_params * )(p->params);
*)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

$include "b97.mpl"

f := (rs, z, xt, xs0, xs1) ->
  b97_f(f_pw, 0.2, params_a_c_ss, 0.006, params_a_c_ab,
        rs, z, xs0, xs1):
