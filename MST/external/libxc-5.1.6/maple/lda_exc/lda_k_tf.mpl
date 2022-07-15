(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_k_tf_params *params;

  assert(p->params != NULL);
  params = (lda_k_tf_params * )(p->params);
*)

f_zeta_k := z -> 1/2*(opz_pow_n(z,5/3) + opz_pow_n(-z,5/3)):

f := (rs, zeta) -> params_a_ax*f_zeta_k(zeta)/rs^2:
