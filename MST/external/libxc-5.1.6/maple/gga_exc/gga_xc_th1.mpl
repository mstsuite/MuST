(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)
(* prefix:
  gga_xc_th1_params *params;

  assert(p->params != NULL);
  params = (gga_xc_th1_params * )(p->params);
*)

params_a_n := 21:

params_a_a := [
    7/6,  8/6,  9/6, 10/6,  8/6,  9/6, 10/6,
   11/6,  9/6, 10/6, 11/6, 12/6,  9/6, 10/6,
   11/6, 12/6,  7/6,  8/6,  9/6, 10/6, 1
]:

params_a_b := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0]:
params_a_c := [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]:
params_a_d := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0]:

$include "th.mpl"

f := (rs, z, xt, xs0, xs1) -> f_th(rs, z, xt, xs0, xs1):
