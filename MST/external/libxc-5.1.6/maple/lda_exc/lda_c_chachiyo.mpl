(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_chachiyo_params *params;

  assert(p->params != NULL);
  params = (lda_c_chachiyo_params * )(p->params);
*)

e0 := rs -> params_a_ap*log(1 + params_a_bp/rs + params_a_bp/rs^2):
e1 := rs -> params_a_af*log(1 + params_a_bf/rs + params_a_bf/rs^2):

f_chachiyo := (rs, zeta) -> e0(rs) + (e1(rs) - e0(rs))*f_zeta(zeta):
f := (rs, zeta) -> f_chachiyo(rs, zeta):
