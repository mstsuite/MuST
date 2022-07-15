(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

(* replace: "BesselI\(0, " -> "xc_bessel_I0(" *)
(* replace: "BesselI\(1, " -> "xc_bessel_I1(" *)

$define xc_dimensions_2d

prhg07_C := (x, u, t) -> (u - 4*t + x^2/2)/4:

(* This is the solution of solve((y-1)*exp(y) = x/Pi) *)
prhg07_y := x -> LambertW(m_max(x/Pi, -0.9999999999) * exp(-1)) + 1:

prhg07_v := y -> Pi/X_FACTOR_2D_C * BesselI(0, y/2):

prhg07_f := (x, u, t) ->
  prhg07_v(prhg07_y(prhg07_C(x, u, t)))/2:

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  mgga_exchange(prhg07_f, rs, z, xs0, xs1, u0, u1, t0, t1):
