(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)
(* prefix:
  mgga_c_m06l_params *params;

  assert(p->params != NULL);
  params = (mgga_c_m06l_params * )(p->params);
*)

$include "mgga_c_vsxc.mpl"
$include "mgga_c_m05.mpl"

m06l_f := (rs, z, xs0, xs1, ts0, ts1) ->
  + m05_f(rs, z, xs0, xs1, ts0, ts1)
  + vsxc_f(rs, z, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  m06l_f(rs, z, xs0, xs1, ts0, ts1):

