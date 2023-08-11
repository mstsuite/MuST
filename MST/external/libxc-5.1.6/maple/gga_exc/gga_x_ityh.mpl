(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: gga_exc *)

$define gga_x_b88_params
$include "gga_x_b88.mpl"
$include "lda_x_erf.mpl"

ityh_enhancement := xs -> b88_f(xs):
ityh_attenuation := a  -> attenuation_erf(a):

ityh_k_GGA := (rs, z, xs) -> sqrt(9*Pi/(2*X_FACTOR_C*ityh_enhancement(xs))) * n_spin(rs, z)^(1/3):
ityh_f_aa  := (rs, z, xs) -> ityh_attenuation(p_a_hyb_omega_0_/(2*ityh_k_GGA(rs, z, xs))):

ityh_f := (rs, z, xs) -> ityh_f_aa(rs, z, xs) * ityh_enhancement(xs):

f := (rs, zeta, xt, xs0, xs1) -> gga_exchange_nsp(ityh_f, rs, zeta, xs0, xs1):
