(*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_chachiyo_mod_params *params;

  assert(p->params != NULL);
  params = (lda_c_chachiyo_mod_params * )(p->params);
*)

(* Functional is based on Chachiyo correlation *)
$include "lda_c_chachiyo.mpl"
(* .. but with a different scaling function *)
g := z -> (opz_pow_n(z,2/3) + opz_pow_n(-z,2/3))/2:
g_zeta := zeta -> 2*(1 - g(zeta)^3):

f_chachiyo := (rs, zeta) -> e0(rs) + (e1(rs) - e0(rs))*g_zeta(zeta):
f := (rs, zeta) -> f_chachiyo(rs, zeta):
