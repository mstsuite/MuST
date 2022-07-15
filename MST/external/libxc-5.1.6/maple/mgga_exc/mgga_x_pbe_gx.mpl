(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define mgga_x_gx_params
$include "mgga_x_gx.mpl"

pbe_gx_mu := 0.001015549:
pbe_gx := x -> 1/(1 + pbe_gx_mu*x^2):

pbe_gx_f := (x, u, t) ->
  gx_f_a(gx_alpha(x, t)) * pbe_gx(x):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(pbe_gx_f, rs, z, xs0, xs1, u0, u1, t0, t1):
