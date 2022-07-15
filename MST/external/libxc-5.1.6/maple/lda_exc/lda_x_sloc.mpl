(*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_x_sloc_params *params;

  assert(p->params != NULL);
  params = (lda_x_sloc_params * )(p->params);
*)

# https://onlinelibrary.wiley.com/doi/full/10.1002/qua.25312

f_sloc := (rs, z) ->
  -params_a_a/(2*(params_a_b + 1)) * n_total(rs)^params_a_b *
  (opz_pow_n(z,params_a_b + 1) + opz_pow_n(-z,params_a_b + 1)):

f := (rs, z) -> f_sloc(rs, z):
