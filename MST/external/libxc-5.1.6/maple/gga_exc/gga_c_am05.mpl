(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_am05_params *params;

  assert(p->params != NULL);
  params = (gga_c_am05_params * )(p->params);
*)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

XX := s -> 1/(1 + params_a_alpha*s^2):
ff := s -> XX(s) + params_a_gamma*(1 - XX(s)):

f := (rs, z, xt, xs0, xs1) -> f_pw(rs, z)*(
  + opz_pow_n( z,1)/2 * ff(X2S*xs0)
  + opz_pow_n(-z,1)/2 * ff(X2S*xs1)
):
