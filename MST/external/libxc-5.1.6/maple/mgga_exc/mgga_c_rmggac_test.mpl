(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

$define lda_c_pw_params
$define lda_c_pw_modified_params
$include "lda_c_pw.mpl"

rmggac_gamma1 := 0.08:
rmggac_gamma2 := 0.3:
rmggac_g := (alpha, s) ->
  (1 + rmggac_gamma1)*alpha/(rmggac_gamma1 + alpha + rmggac_gamma2*s^2):

rmggac_f2 := (alpha, s) ->
  3*rmggac_g(alpha, s)^3/(1 + rmggac_g(alpha, s)^3 + rmggac_g(alpha, s)^6):
rmggac_f1 := (alpha, s) ->
  1 - rmggac_f2(alpha, s):

rmggac_gamma := 0.031091:

(* from gga_c_regtpss *)
beta_a := 0.066724550603149220:
beta_b := 0.1:
beta_c := 0.1778:
mbeta := (rs, t) -> beta_a*(1 + beta_b*rs)/(1 + beta_c*rs):
(* from mgga_c_r2scan *)
rmggac_w1 := (rs, z) -> exp(-f_pw(rs, z)/(rmggac_gamma*mphi(z)^3)) - 1:
(* from gga_c_pbe *)
A := (rs, z, t) -> mbeta(rs, t)/(rmggac_gamma * rmggac_w1(rs, z)):
(* from gga_c_scan_e0 *)
scan_e0_g := (rs, z, t) -> (1 + 4*A(rs, z, t)*t^2)^(-1/4):
(* from mgga_c_r2scan *)
rmggac_H1 := (rs, z, t) -> rmggac_gamma*mphi(z)^3*log(1 + rmggac_w1(rs, z) * (1 - scan_e0_g(rs, z, t))):

(* from mgga_c_scan *)
scan_alpha := (z, xt, ts0, ts1) ->
  (t_total(z, ts0, ts1) - xt^2/8)/(K_FACTOR_C*t_total(z, 1, 1)):
scan_b1c := 0.0285764:
scan_b2c := 0.0889:
scan_b3c := 0.125541:
scan_eclda0 := rs -> -scan_b1c/(1 + scan_b2c*sqrt(rs) + scan_b3c*rs):

scan_chi_infty := 0.12802585262625815:
scan_g_infty := s -> 1/(1 + 4*scan_chi_infty*s^2)^(1/4):
scan_G_cnst := 2.3631:
scan_Gc := z -> (1 - scan_G_cnst*(2^(1/3) - 1)*f_zeta(z))*(1 - z^12):

scan_H0 := (rs, s) ->
  scan_b1c*log(1 + (exp(-scan_eclda0(rs)/scan_b1c) - 1)*(1 - scan_g_infty(s))):
scan_e0 := (rs, z, s) ->
  (scan_eclda0(rs) + scan_H0(rs, s))*scan_Gc(z):

(* define the functional *)
rmggac_eps1 := (rs, z, t) ->
  (f_pw(rs, z) +  rmggac_H1(rs, z, t)):

rmggac_f := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  + scan_e0(rs, z, X2S*2^(1/3)*xt)
      * rmggac_f1(scan_alpha(z, xt, ts0, ts1), X2S*2^(1/3)*xt) 
  + rmggac_eps1(rs, z, tt(rs, z, xt))
      * rmggac_f2(scan_alpha(z, xt, ts0, ts1), X2S*2^(1/3)*xt):

f := (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  rmggac_f(rs, z, xt, xs0, xs1, ts0, ts1):
