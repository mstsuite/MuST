(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define xc_dimensions_2d

_2d_b88_beta := 0.007:
_2d_b88_csi  := 8: (* for harmonic potentials *)

_2d_b88_f := x -> 1 + _2d_b88_beta/X_FACTOR_2D_C*x^2/(1 + _2d_b88_csi*_2d_b88_beta*x*arcsinh(x)):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange(_2d_b88_f, rs, zeta, xs0, xs1):
