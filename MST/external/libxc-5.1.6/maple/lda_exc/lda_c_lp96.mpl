(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)
(* prefix:
  lda_c_lp96_params *params;

  assert(p->params != NULL);
  params = (lda_c_lp96_params * )(p->params);
*)

f := (rs, zeta) -> params_a_C1 + params_a_C2*n_total(rs)^(-1/3) + params_a_C3*n_total(rs)^(-2/3):

