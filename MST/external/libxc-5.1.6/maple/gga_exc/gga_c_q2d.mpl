(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_c_pbe_params
$include "gga_c_pbe.mpl"

$include "lda_c_2d_amgb.mpl"
$undef xc_dimensions_2d

rs2D_factor := 1.704:
q2d_dd := 1e6:

q2d_rs2D := (rs, xt) -> rs2D_factor*rs*sqrt(X2S*xt)/RS_FACTOR:

q2d_fac := t -> t^4*(1 + t^2)/(q2d_dd + t^6):

q2d_f := (rs, z, xt, xs0, xs1) ->
  (1 - q2d_fac(tt(rs, z, xt)))*f_pbe(rs, z, xt, xs0, xs1) + q2d_fac(tt(rs, z, xt))*f_amgb(q2d_rs2D(rs, xt), z):

f := (rs, z, xt, xs0, xs1) -> q2d_f(rs, z, xt, xs0, xs1):
