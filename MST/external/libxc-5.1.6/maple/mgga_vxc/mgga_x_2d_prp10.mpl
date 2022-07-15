(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_vxc *)

(* replace: "BesselI\(0, " -> "xc_bessel_I0(" *)
(* replace: "BesselI\(1, " -> "xc_bessel_I1(" *)

$define xc_dimensions_2d
$include "mgga_x_2d_prhg07.mpl"

# we add here a threshold of 1e-10 for either tau or the Fermi curvature
prp10_f := (rs, z, x, u, t) -> -(X_FACTOR_2D_C*prhg07_v(prhg07_y(prhg07_C(x, u, t)))
  - (2*sqrt(2)/(3*Pi)) * sqrt(2*m_max(t - x^2/8, 1e-10)))*n_spin(rs, z)^(1/2):

f := (rs, z, xt, xs0, xs1, u0, u1, t0, t1) ->
  prp10_f(rs, z, xs0, u0, t0):
