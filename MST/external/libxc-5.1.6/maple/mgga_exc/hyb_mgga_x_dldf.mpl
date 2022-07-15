(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

params_a_kappa := 4.8827323:
params_a_mu    := 0.3511128:
$include "gga_x_pbe.mpl"

dldf_a := [1, -0.1637571, -0.1880028, -0.4490609, -0.0082359]:
dldf_csi_HF := 1 - 0.6144129:

dldf_f := (x, u, t) ->
  + dldf_csi_HF*pbe_f(x)*mgga_series_w(dldf_a, 5, t):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(dldf_f, rs, z, xs0, xs1, u0, u1, t0, t1):
