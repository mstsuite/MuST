(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);
*)

$define xc_dimensions_1d

f := (rs, zeta) -> \
 (params_a_a1 + params_a_a2*n_total(rs) + params_a_a3*n_total(rs)^2) * n_total(rs)^params_a_alpha:
