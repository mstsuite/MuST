(*
 Copyright (C) 2017 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_c_p86vwn_params *params;

  assert(p->params != NULL);
  params = (gga_c_p86vwn_params * )(p->params);
*)

$include "gga_c_p86.mpl"
$include "lda_c_vwn.mpl"

f_p86 := (rs, z, xt, xs0, xs1) ->
  f_vwn(rs, z) + p86_H(rs, z, xt):

f := (rs, z, xt, xs0, xs1) ->
  f_p86(rs, z, xt, xs0, xs1):
