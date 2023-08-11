(*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*)

(* type: mgga_exc *)

b97mv_par_n := 6:

b97mv_gamma_x := 0.004:
b97mv_par_x := [
    [  0.85,  0, 0],
    [  1.007, 0, 1],
    [  0.259, 1, 0],
    [  0,   0, 0],
    [  0,   0, 0],
    [  0,   0, 0]
]:

b97mv_gamma_ss := 0.2:
b97mv_par_ss := [
    [  0.443,  0, 0],
    [ -1.437,  0, 4],
    [ -4.535,  1, 0],
    [ -3.39,   2, 0],
    [  4.278,  4, 3],
    [  0,    0, 0]
]:

b97mv_gamma_os := 0.006:
b97mv_par_os := [
    [  1.000,  0, 0],
    [  1.358,  1, 0],
    [  2.924,  2, 0],
    [ -8.812,  2, 1],
    [ -1.39,   6, 0],
    [  9.142,  6, 1]
]:

$include "lda_x_erf.mpl"
$include "b97mv.mpl"

wb97mv_f := (rs, z, xs0, xs1, ts0, ts1) ->
   my_piecewise3(screen_dens_zeta(rs,  z), 0, (1 + z)/2 * lda_x_erf_spin(rs*(2/(1 + z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_par_x,  b97mv_par_n, xs0, ts0, 0))
+  my_piecewise3(screen_dens_zeta(rs, -z), 0, (1 - z)/2 * lda_x_erf_spin(rs*(2/(1 - z))^(1/3),  1)
    * b97mv_g(b97mv_gamma_x,  b97mv_wx_ss, b97mv_par_x,  b97mv_par_n, xs1, ts1, 0)):

f :=  (rs, z, xt, xs0, xs1, us0, us1, ts0, ts1) ->
  wb97mv_f(rs, z, xs0, xs1, ts0, ts1) +
  b97mv_f(rs, z, xs0, xs1, ts0, ts1):
