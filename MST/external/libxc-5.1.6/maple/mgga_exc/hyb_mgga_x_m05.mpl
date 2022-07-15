(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_x_m05_params *params;

  assert(p->params != NULL);
  params = (mgga_x_m05_params * )(p->params);
*)

$define gga_x_pbe_params
$include "gga_x_pbe.mpl"

m05_f := (x, u, t) ->
  + params_a_csi_HF*pbe_f(x)*mgga_series_w(params_a_a, 12, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(m05_f, rs, z, xs0, xs1, u0, u1, t0, t1):
