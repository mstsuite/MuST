(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

pkzb_c := 0.53:

pkzb_perp := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  + (1 + pkzb_c*(t_total(z, xs0^2, xs1^2)/(8*t_total(z, ts0, ts1)))^2)
  * f_pbe(rs, z, xt, xs0, xs1):

pkzb_par  := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  - (1 + pkzb_c)*(
    + (xs0^2/(8*ts0))^2*gga_stoll_par(f_pbe, rs,  z, xs0,  1)
    + (xs1^2/(8*ts1))^2*gga_stoll_par(f_pbe, rs, -z, xs1, -1)
):

pkzb_f := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  pkzb_perp(rs, z, xt, xs0, xs1, ts0, ts1) + pkzb_par(rs, z, xt, xs0, xs1, ts0, ts1):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  pkzb_f(rs, z, xt, xs0, xs1, ts0, ts1):
