(*
 Copyright (C) 2017 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: lda_exc *)

$include "lda_x_yukawa.mpl"
$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

bn05_A  := 3.4602:
bn05_C0 := 3.2:
bn05_C1 := -0.9:

f := (rs, z) ->
  f_lda_x_yukawa(rs, z) + f_pw(rs, z)*bn05_A/(bn05_C0 + bn05_C1*rs + rs^2):
