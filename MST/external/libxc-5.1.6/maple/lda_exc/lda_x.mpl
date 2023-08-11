(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_x_params *params;

  assert(p->params != NULL);
  params = (lda_x_params * )(p->params);
*)

$ifdef lda_x_params
params_a_alpha := 1:
$endif

f_lda_x := (rs, z) ->
  + params_a_alpha*my_piecewise3(screen_dens(rs,  z), 0, lda_x_spin(rs,  z))
  + params_a_alpha*my_piecewise3(screen_dens(rs, -z), 0, lda_x_spin(rs, -z)):
f       := (rs, z) -> f_lda_x(rs, z):
