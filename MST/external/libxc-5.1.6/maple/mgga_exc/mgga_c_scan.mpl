(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$include "gga_c_scan_e0.mpl"
$include "mgga_x_scan.mpl"

scan_b1c := 0.0285764:
scan_b2c := 0.0889:
scan_b3c := 0.125541:
scan_eclda0 := rs -> -scan_b1c/(1 + scan_b2c*sqrt(rs) + scan_b3c*rs):

scan_chi_infty := 0.12802585262625815:
scan_g_infty := s -> 1/(1 + 4*scan_chi_infty*s^2)^(1/4):

(* in the paper it is 2.3631 *)
scan_G_cnst := 2.363:
scan_Gc := z -> (1 - scan_G_cnst*(2^(1/3) - 1)*f_zeta(z))*(1 - z^12):

scan_H0 := (rs, s) ->
  scan_b1c*log(1 + (exp(-scan_eclda0(rs)/scan_b1c) - 1)*(1 - scan_g_infty(s))):
scan_e0 := (rs, z, s) ->
  (scan_eclda0(rs) + scan_H0(rs, s))*scan_Gc(z):

scan_alpha := (z, xt, ts0, ts1) ->
  (t_total(z, ts0, ts1) - xt^2/8)/(K_FACTOR_C*t_total(z, 1, 1)):

(* set parameters of f_alpha *)
params_a_c1 := 0.64:
params_a_c2 := 1.5:
params_a_d  := 0.7:

scan_f := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  f_pbe(rs, z, xt, xs0, xs1) + scan_f_alpha(scan_alpha(z, xt, ts0, ts1))*(
    + scan_e0(rs, z, X2S*2^(1/3)*xt)
    - f_pbe(rs, z, xt, xs0, xs1)
  ):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  scan_f(rs, z, xt, xs0, xs1, ts0, ts1):
