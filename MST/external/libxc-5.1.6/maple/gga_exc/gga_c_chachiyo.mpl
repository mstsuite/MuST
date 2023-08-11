(*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_chachiyo_params *params;

  assert(p->params != NULL);
  params = (gga_c_chachiyo_params * )(p->params);
*)

(* Functional is based on Chachiyo correlation with modified spin scaling *)
$include "lda_c_chachiyo_mod.mpl"

(* Reduced gradient parameter *)
cha_t := (rs, xt) -> (Pi/3)^(1/6) / 4 * n_total(rs)^(1/6) * xt:

(* The full functional that agrees with the given reference values is *)
f_chachiyo_gga := (rs, z, xt, xs0, xs1) -> f_chachiyo(rs, z) * (1 + cha_t(rs, xt)^2)^(params_a_h / f_chachiyo(rs, z)):

f  := (rs, z, xt, xs0, xs1) -> f_chachiyo_gga(rs, z, xt, xs0, xs1):
